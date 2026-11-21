#!/usr/bin/env python3
"""4-uMLIP Ensemble Single-Point Physical Evaluator for 198 Semiconductor Molecules.

Evaluates the 198 benchmark & patent-expansion molecules (160 unified benchmark +
38 patent/literature expansion) across four heterogeneous uMLIPs on identical
3D geometries:
  1. ``SevenNet-0`` (precomputed LBFGS-relaxed baseline in DB JSONs)
  2. ``MACE-OMol`` (OMol25 molecular/organometallic equivariant GNN, via ``.venv_mace``)
  3. ``MACE-MP-0`` (MPtrj universal equivariant GNN, via ``.venv_mace``)
  4. ``ORB-v3`` (non-equivariant fast GNN, via ``.venv_orb``)

Aggregates per-model ``BDE_min`` (# units: kcal/mol), ``E_coh`` (# units: eV/atom),
and ``Delta E_dimer`` (# units: kcal/mol) via :func:`compute_ensemble_sapis_metrics`
into ``database/external/umlip_ensemble_198.json``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

EV_TO_KCAL: float = 23.060548  # units: (kcal/mol) / eV

GEOMS_PATH = os.path.join(PROJECT_ROOT, "database", "external", "ensemble_geoms_198.json")
MACE_OUT_PATH = os.path.join(PROJECT_ROOT, "database", "external", "ensemble_mace_198.json")
ORB_OUT_PATH = os.path.join(PROJECT_ROOT, "database", "external", "ensemble_orb_198.json")
FINAL_OUT_PATH = os.path.join(PROJECT_ROOT, "database", "external", "umlip_ensemble_198.json")


def _build_dimer_positions(pos: np.ndarray, gap_ang: float = 3.2) -> np.ndarray:
    """Create a deterministic non-overlapping dimer by shifting a copy along +x (# units: Å)."""
    pos = np.asarray(pos, dtype=float)
    span_x = float(np.max(pos[:, 0]) - np.min(pos[:, 0])) if len(pos) > 1 else 1.5
    shift = np.array([span_x + gap_ang, 0.0, 0.0], dtype=float)  # units: Å
    return np.vstack([pos, pos + shift])


def _extract_ligand_split(mol_h: Any, c_idx: int) -> Optional[Tuple[List[int], List[int]]]:
    """Split a molecule into one detached ligand subtree and the remaining core+ligands."""
    c_atom = mol_h.GetAtomWithIdx(c_idx)
    nbs = [nb.GetIdx() for nb in c_atom.GetNeighbors()]
    if not nbs:
        return None
    n_atoms = mol_h.GetNumAtoms()
    # Pick the smallest non-trivial connected component attached to c_idx (avoids ring breakage)
    best_lig: Optional[List[int]] = None
    for nb_idx in nbs:
        visited = set()
        stack = [nb_idx]
        hits_other_nb = False
        while stack:
            cur = stack.pop()
            if cur == c_idx or cur in visited:
                continue
            visited.add(cur)
            for nxt in mol_h.GetAtomWithIdx(cur).GetNeighbors():
                ni = nxt.GetIdx()
                if ni == c_idx and cur != nb_idx:
                    hits_other_nb = True
                elif ni != c_idx and ni not in visited:
                    stack.append(ni)
        if len(visited) < n_atoms - 1:
            if best_lig is None or (not hits_other_nb and len(visited) < len(best_lig)):
                best_lig = sorted(visited)
    if not best_lig:
        return None
    lig_set = set(best_lig)
    rem = [i for i in range(n_atoms) if i not in lig_set]
    return best_lig, rem


def export_shared_geometries() -> None:
    """Export shared 3D monomer, fragment, and dimer geometries for all 198 molecules."""
    from rdkit import Chem
    from element_physics import embed_etkdg_3d, find_central_atom_intrinsic, get_element_properties
    from tmc_sa_scorer import build_coordination_template_3d

    with open(os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json"), "r", encoding="utf-8") as fh:
        u160 = json.load(fh)["records"]

    with open(os.path.join(PROJECT_ROOT, "database", "external", "patent_lit_expansion_candidates.json"), "r", encoding="utf-8") as fh:
        pat_raw = json.load(fh)
    pat38 = pat_raw.get("records", pat_raw.get("candidates", pat_raw))

    with open(os.path.join(PROJECT_ROOT, "database", "external", "patent_expansion_sevenn_features.json"), "r", encoding="utf-8") as fh:
        pat7_raw = json.load(fh)
    pat7_list = pat7_raw.get("records", pat7_raw.get("results", pat7_raw))
    if isinstance(pat7_list, dict):
        pat7_list = list(pat7_list.values())
    pat7_by_id = {r.get("id"): r for r in pat7_list if isinstance(r, dict)}

    combined: List[Dict[str, Any]] = []
    seen_ids = set()

    for r in u160:
        rid = str(r["id"])
        seen_ids.add(rid)
        combined.append({
            "id": rid,
            "source_set": "unified_160",
            "name": r.get("name", ""),
            "smiles": r.get("smiles", ""),
            "base_sa": float(r.get("ligand_sa") if r.get("ligand_sa") is not None else r.get("sascore", 2.0)),
            "vbur_val": float(r.get("vbur_percent", 50.0)),
            "v_min": float(r.get("v_min_adaptive", 44.0)),
            "is_d_block": bool(get_element_properties(r.get("center_atom", "C"), 6).get("is_d_block", False)) if r.get("center_atom") else False,
            "sevenn_0": {
                "bde_min": r.get("bde_min_kcal_mol"),
                "e_coh": r.get("e_coh_ev_per_atom"),
                "de_dimer": r.get("delta_e_dimer_kcal_mol"),
            },
            "tmc_sa_v2": r.get("tmc_sa_v2"),
            "unified_pass": r.get("unified_pass"),
            "gate_reasons": r.get("gate_reasons") or [],
        })

    for r in pat38:
        rid = str(r["id"])
        if rid in seen_ids:
            continue
        p7 = pat7_by_id.get(rid, {})
        combined.append({
            "id": rid,
            "source_set": "patent_38",
            "name": r.get("name", ""),
            "smiles": r.get("smiles", ""),
            "base_sa": float(r.get("ligand_sa") if r.get("ligand_sa") is not None else r.get("sascore", 2.0)),
            "vbur_val": float(p7.get("vbur_percent") or r.get("vbur_d2_7net") or r.get("vbur_d2_percent") or 50.0),
            "v_min": 44.0,
            "sevenn_0": {
                "bde_min": p7.get("bde_min_kcal_mol", r.get("bde_min_kcal_mol")),
                "e_coh": p7.get("e_coh_ev_per_atom", r.get("e_coh_ev_per_atom")),
                "de_dimer": p7.get("delta_e_dimer_kcal_mol", r.get("delta_e_dimer_kcal_mol")),
            },
            "tmc_sa_v2": r.get("tmc_sa_7net", r.get("tmc_sa_no3d")),
            "unified_pass": r.get("unified_pass_7net", r.get("unified_pass")),
            "gate_reasons": [],
        })

    geoms: List[Dict[str, Any]] = []
    for item in combined:
        smi = item["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        c_idx, c_sym, c_z = find_central_atom_intrinsic(mol)
        props = get_element_properties(c_sym, c_z)
        try:
            emb = embed_etkdg_3d(smi)
        except Exception:
            emb = None
        if emb is None:
            emb = build_coordination_template_3d(mol, c_idx)
        if emb is None:
            continue
        mol_h, pos, syms, c_idx_h, c_sym, c_z = emb
        split = _extract_ligand_split(mol_h, c_idx_h)
        dimer_pos = _build_dimer_positions(pos, gap_ang=3.2)

        # Count hydrogens directly bonded to center atom
        c_atom = mol_h.GetAtomWithIdx(c_idx_h)
        center_h = sum(1 for nb in c_atom.GetNeighbors() if nb.GetSymbol() == "H")

        geoms.append({
            **item,
            "center_sym": c_sym,
            "center_z": int(c_z),
            "is_d_block": bool(props.get("is_d_block", False)),
            "is_coord_center": bool(props.get("is_coord_center", False)),
            "chi_pauling": float(props.get("chi_pauling", 2.0)),
            "center_h": int(center_h),
            "syms": syms,
            "pos": np.round(pos, 4).tolist(),  # units: Å
            "dimer_syms": syms + syms,
            "dimer_pos": np.round(dimer_pos, 4).tolist(),  # units: Å
            "lig_idxs": split[0] if split else [],
            "rem_idxs": split[1] if split else [],
        })

    with open(GEOMS_PATH, "w", encoding="utf-8") as fh:
        json.dump({"count": len(geoms), "records": geoms}, fh, ensure_ascii=False)
    print(f"Exported {len(geoms)} shared 3D geometries to {GEOMS_PATH}")


def _eval_single_point_metrics(
    rec: Dict[str, Any],
    calc: Any,
    atom_ref_cache: Dict[str, float],
    set_spin_charge: bool = False,
) -> Dict[str, Optional[float]]:
    """Compute single-point BDE_min (kcal/mol), E_coh (eV/atom), and Delta E_dimer (kcal/mol)."""
    from ase import Atoms

    def _sp(syms: List[str], pos: List[List[float]]) -> float:
        atoms = Atoms(symbols=syms, positions=np.asarray(pos, dtype=float))
        if set_spin_charge:
            atoms.info["charge"] = 0
            # Odd total electron count -> doublet (spin=2), even -> singlet (spin=1)
            tot_z = int(np.sum(atoms.get_atomic_numbers()))
            atoms.info["spin"] = 2 if (tot_z % 2 == 1) else 1
        atoms.calc = calc
        return float(atoms.get_potential_energy())  # units: eV

    syms = rec["syms"]
    pos = rec["pos"]
    e_mono = _sp(syms, pos)  # units: eV

    # Cohesive energy (atomization energy per atom, eV/atom) — skip gracefully if 1-atom 0-edge unsupported
    e_coh: Optional[float] = None
    try:
        e_atoms_sum = 0.0
        for s in syms:
            if s not in atom_ref_cache:
                atom_ref_cache[s] = _sp([s], [[0.0, 0.0, 0.0]])
            e_atoms_sum += atom_ref_cache[s]
        e_coh = round((e_atoms_sum - e_mono) / max(1, len(syms)), 4)  # units: eV/atom
    except Exception:
        e_coh = None

    # BDE_min via rigid fragment single-point dissociation (kcal/mol)
    bde_kcal: Optional[float] = None
    lig_idxs = rec.get("lig_idxs") or []
    rem_idxs = rec.get("rem_idxs") or []
    if lig_idxs and rem_idxs:
        try:
            e_lig = _sp([syms[i] for i in lig_idxs], [pos[i] for i in lig_idxs])
            e_rem = _sp([syms[i] for i in rem_idxs], [pos[i] for i in rem_idxs])
            bde_kcal = round((e_lig + e_rem - e_mono) * EV_TO_KCAL, 3)  # units: kcal/mol
        except Exception:
            bde_kcal = None

    # Dimerization interaction energy (kcal/mol)
    de_dimer: Optional[float] = None
    try:
        e_dimer = _sp(rec["dimer_syms"], rec["dimer_pos"])
        de_dimer = round((e_dimer - 2.0 * e_mono) * EV_TO_KCAL, 3)  # units: kcal/mol
    except Exception:
        de_dimer = None

    return {
        "bde_min": bde_kcal,
        "e_coh": e_coh,
        "de_dimer": de_dimer,
    }


def run_eval_mace() -> None:
    """Evaluate MACE-OMol and MACE-MP-0 on all 198 shared geometries (run inside .venv_mace)."""
    from mace.calculators import mace_mp, mace_omol

    with open(GEOMS_PATH, "r", encoding="utf-8") as fh:
        geoms = json.load(fh)["records"]

    t0 = time.time()
    calc_omol = mace_omol(model="extra_large", device="cpu", default_dtype="float32")
    calc_mp = mace_mp(model="small", device="cpu", default_dtype="float32")
    print(f"Loaded MACE-OMol + MACE-MP-0 (float32) in {time.time() - t0:.1f}s")

    cache_omol: Dict[str, float] = {}
    cache_mp: Dict[str, float] = {}
    results: Dict[str, Dict[str, Any]] = {}
    if os.path.exists(MACE_OUT_PATH):
        try:
            with open(MACE_OUT_PATH, "r", encoding="utf-8") as fh:
                results = json.load(fh)
        except Exception:
            results = {}

    for idx, rec in enumerate(geoms, start=1):
        rid = rec["id"]
        entry: Dict[str, Any] = dict(results.get(rid, {}))
        if "MACE-OMol" not in entry:
            entry.pop("MACE-OMol_err", None)
            try:
                entry["MACE-OMol"] = _eval_single_point_metrics(rec, calc_omol, cache_omol, set_spin_charge=True)
            except Exception as exc:
                entry["MACE-OMol_err"] = str(exc)
        if "MACE-MP-0" not in entry:
            entry.pop("MACE-MP-0_err", None)
            try:
                entry["MACE-MP-0"] = _eval_single_point_metrics(rec, calc_mp, cache_mp, set_spin_charge=False)
            except Exception as exc:
                entry["MACE-MP-0_err"] = str(exc)
        results[rid] = entry
        if idx % 20 == 0 or idx == len(geoms):
            with open(MACE_OUT_PATH, "w", encoding="utf-8") as fh:
                json.dump(results, fh, indent=2, ensure_ascii=False)
            print(f"[MACE] {idx}/{len(geoms)} done ({time.time() - t0:.1f}s)")


def run_eval_orb() -> None:
    """Evaluate ORB-v3 on all 198 shared geometries (run inside .venv_orb)."""
    from orb_models.forcefield import pretrained
    from orb_models.forcefield.calculator import ORBCalculator

    with open(GEOMS_PATH, "r", encoding="utf-8") as fh:
        geoms = json.load(fh)["records"]

    t0 = time.time()
    orb_name = "ORB-v3-omol"
    if hasattr(pretrained, "orb_v3_conservative_omol"):
        orbff = pretrained.orb_v3_conservative_omol(device="cpu", precision="float32-high", compile=False)
    elif hasattr(pretrained, "orb_v3_direct_omol"):
        orbff = pretrained.orb_v3_direct_omol(device="cpu", precision="float32-high", compile=False)
    elif hasattr(pretrained, "orb_v3_conservative_inf_omat"):
        orb_name = "ORB-v3-omat"
        orbff = pretrained.orb_v3_conservative_inf_omat(device="cpu", precision="float32-high", compile=False)
    else:
        orb_name = "ORB-v2"
        orbff = pretrained.orb_v2(device="cpu")
    calc_orb = ORBCalculator(orbff, device="cpu")
    print(f"Loaded {orb_name} in {time.time() - t0:.1f}s")

    cache_orb: Dict[str, float] = {}
    results: Dict[str, Dict[str, Any]] = {}
    is_omol = "omol" in orb_name.lower()

    for idx, rec in enumerate(geoms, start=1):
        rid = rec["id"]
        try:
            results[rid] = {orb_name: _eval_single_point_metrics(rec, calc_orb, cache_orb, set_spin_charge=is_omol)}
        except Exception as exc:
            results[rid] = {f"{orb_name}_err": str(exc)}
        if idx % 20 == 0 or idx == len(geoms):
            with open(ORB_OUT_PATH, "w", encoding="utf-8") as fh:
                json.dump({"model_name": orb_name, "records": results}, fh, indent=2, ensure_ascii=False)
            print(f"[{orb_name}] {idx}/{len(geoms)} done ({time.time() - t0:.1f}s)")


def run_aggregate() -> None:
    """Aggregate SevenNet-0 + MACE-OMol + MACE-MP-0 + ORB-v3 into database/external/umlip_ensemble_198.json."""
    from tmc_sa_scorer import compute_ensemble_sapis_metrics

    with open(GEOMS_PATH, "r", encoding="utf-8") as fh:
        geoms = json.load(fh)["records"]
    mace_res = json.load(open(MACE_OUT_PATH, "r", encoding="utf-8")) if os.path.exists(MACE_OUT_PATH) else {}
    orb_raw = json.load(open(ORB_OUT_PATH, "r", encoding="utf-8")) if os.path.exists(ORB_OUT_PATH) else {}
    orb_model_name = orb_raw.get("model_name", "ORB-v3-omol")
    orb_res = orb_raw.get("records", {})

    out_records: Dict[str, Dict[str, Any]] = {}
    grade_counts: Dict[str, int] = {"High": 0, "Medium": 0, "Low (DFT Flag)": 0}

    for rec in geoms:
        rid = rec["id"]
        model_metrics: Dict[str, Dict[str, Optional[float]]] = {}
        s7 = rec.get("sevenn_0") or {}
        if any(v is not None for v in s7.values()):
            model_metrics["SevenNet-0"] = s7
        m_entry = mace_res.get(rid, {})
        if "MACE-OMol" in m_entry:
            model_metrics["MACE-OMol"] = m_entry["MACE-OMol"]
        if "MACE-MP-0" in m_entry:
            model_metrics["MACE-MP-0"] = m_entry["MACE-MP-0"]
        o_entry = orb_res.get(rid, {})
        if orb_model_name in o_entry:
            model_metrics[orb_model_name] = o_entry[orb_model_name]

        ens = compute_ensemble_sapis_metrics(
            base_sa=float(rec["base_sa"]),
            vbur_val=float(rec["vbur_val"]),
            v_min=float(rec["v_min"]),
            is_coord_center=bool(rec["is_coord_center"]),
            model_metrics=model_metrics,
            is_d_block=bool(rec["is_d_block"]),
            chi_pauling=float(rec["chi_pauling"]),
            center_h=int(rec["center_h"]),
        )
        grade_counts[ens["confidence_grade"]] = grade_counts.get(ens["confidence_grade"], 0) + 1
        out_records[rid] = {
            "id": rid,
            "name": rec["name"],
            "source_set": rec["source_set"],
            "per_model_metrics": model_metrics,
            **ens,
        }

    payload = {
        "metadata": {
            "models": ["SevenNet-0", "MACE-OMol", "MACE-MP-0", orb_model_name],
            "n_molecules": len(out_records),
            "confidence_summary": grade_counts,
        },
        "records": out_records,
    }
    with open(FINAL_OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print(f"Saved 4-uMLIP ensemble results ({len(out_records)} molecules) to {FINAL_OUT_PATH}")
    print("Confidence distribution:", grade_counts)


def main() -> None:
    parser = argparse.ArgumentParser(description="4-uMLIP Ensemble Runner for 198 Semiconductor Molecules")
    parser.add_argument("--mode", choices=["export-geoms", "eval-mace", "eval-orb", "aggregate"], required=True)
    args = parser.parse_args()
    if args.mode == "export-geoms":
        export_shared_geometries()
    elif args.mode == "eval-mace":
        run_eval_mace()
    elif args.mode == "eval-orb":
        run_eval_orb()
    elif args.mode == "aggregate":
        run_aggregate()


if __name__ == "__main__":
    main()
