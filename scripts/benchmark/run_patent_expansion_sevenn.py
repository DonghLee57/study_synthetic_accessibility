#!/usr/bin/env python3
"""SevenNet-0 Physical Descriptor Enrichment for Patent & Literature Expansion DB.

Reads ``database/external/patent_lit_expansion_candidates.json`` (38 curated
KR/JP/CN/US/TW patent and metal-centered literature candidates), executes
``SevenNet-0`` LBFGS monomer relaxation, non-overlapping dimerization, fragment/stretch
BDE evaluation, Design 2 ``%V_bur^D2`` & open attack void cone on relaxed coordinates,
and standard-state formation enthalpy per coordination bond (``Delta E_f^std / v_M``),
saving outputs to ``database/external/patent_expansion_sevenn_features.json``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List

import torch
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sevenn.calculator import SevenNetCalculator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from element_physics import (  # noqa: E402
    calc_buried_volume_sphere,
    calc_design2_steric_electronic_metrics,
    compute_coordination_state,
    compute_formation_energy_metrics,
    compute_intrinsic_topology,
    compute_refined_pi_hapticity,
    compute_stat_mech_tdelta_s_dimer,
    find_central_atom_intrinsic,
    get_element_properties,
)
from run_phase2_sevenn_features import compute_precursor_physical_features  # noqa: E402

IN_PATH = os.path.join(PROJECT_ROOT, "database", "external", "patent_lit_expansion_candidates.json")
OUT_PATH = os.path.join(PROJECT_ROOT, "database", "external", "patent_expansion_sevenn_features.json")


def _write_payload(results: List[Dict[str, Any]], elapsed_s: float) -> None:
    """Persist computed SevenNet features deterministically."""
    by_el: Dict[str, int] = {}
    for r in results:
        el = str(r.get("center_atom") or "-")
        by_el[el] = by_el.get(el, 0) + 1
    payload = {
        "metadata": {
            "framework": "SevenNet (7net-0 LBFGS Relaxed) + Design 1-2 State-Adapted Sphere",
            "source_database": "database/external/patent_lit_expansion_candidates.json",
            "total_molecules_computed": len(results),
            "n_precursors": sum(1 for r in results if r.get("target_type") == "precursor"),
            "n_inhibitors": sum(1 for r in results if r.get("target_type") == "inhibitor"),
            "counts_by_element": by_el,
            "fallback_count": sum(1 for r in results if r.get("is_fallback")),
            "computation_time_seconds": round(elapsed_s, 1),
        },
        "features": results,
    }
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def run_patent_expansion_sevenn(num_threads: int = 6, resume: bool = True) -> None:
    """Compute SevenNet-0 physical features for all 38 patent/literature expansion records.

    Physics/Algorithm:
        Reuses ``compute_precursor_physical_features(..., return_coords=True)`` for
        monomer LBFGS relaxation (``fmax = 0.04 eV/A``), complementary fragment / constrained
        ring-stretch BDE (kcal/mol), cohesive energy (eV/atom), and dimerization interaction
        (``Delta E_dimer``, kcal/mol), then evaluates state-adapted Design 2 steric crowding
        (``%V_bur^D2``, %) and thermodynamic formation energy per bond directly on the
        SevenNet-0 relaxed coordinates.
    """
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    torch.set_num_threads(num_threads)

    with open(IN_PATH, "r", encoding="utf-8") as fh:
        candidates = json.load(fh)["records"]

    cached_by_id: Dict[str, Dict[str, Any]] = {}
    if resume and os.path.exists(OUT_PATH):
        with open(OUT_PATH, "r", encoding="utf-8") as fh:
            prev = json.load(fh)
        for rec in prev.get("features", []):
            if rec.get("id") and not rec.get("is_fallback"):
                cached_by_id[rec["id"]] = rec

    print(
        f"Loaded {len(candidates)} patent/literature candidates "
        f"(cached={len(cached_by_id)}, to_compute={len(candidates) - len(cached_by_id)}).",
        flush=True,
    )

    calc = SevenNetCalculator(model="7net-0", device="cpu")
    frag_cache: dict = {}
    results_by_id = dict(cached_by_id)
    results: List[Dict[str, Any]] = []
    t_start = time.time()

    for idx, item in enumerate(candidates, 1):
        eid = item["id"]
        if eid in cached_by_id:
            results.append(cached_by_id[eid])
            print(f"[{idx:02d}/{len(candidates):02d}] {eid} | {item['name']:22s} | cached", flush=True)
            continue

        t0 = time.time()
        smi = item["smiles"]
        p_data = compute_precursor_physical_features(smi, calc, frag_cache=frag_cache, return_coords=True)
        mol_no_h = Chem.MolFromSmiles(smi)
        assert mol_no_h is not None, f"Invalid SMILES: {smi}"
        formula = rdMolDescriptors.CalcMolFormula(mol_no_h)

        # Intrinsic element physics and Design 2 geometry directly on SevenNet-relaxed coordinates
        topo = compute_intrinsic_topology(smi)
        c_idx_nh, c_sym, c_z = find_central_atom_intrinsic(mol_no_h)
        props = get_element_properties(c_sym, c_z)
        pi_h = compute_refined_pi_hapticity(mol_no_h, c_idx_nh)
        cstate = compute_coordination_state(mol_no_h, c_idx_nh, pi_hapticity=pi_h)

        mol_rel = p_data["mol_relaxed"]
        pos_rel = p_data["pos_aligned"]  # units: A
        syms_rel = p_data["symbols"]
        c_idx_rel = int(p_data["center_idx"])
        c_deg = max(1, mol_rel.GetAtomWithIdx(c_idx_rel).GetDegree())

        vbur_rprobe_7n = calc_buried_volume_sphere(
            pos_rel, syms_rel, c_idx_rel, radius_ang=float(props["r_probe_ang"]), grid_step_ang=0.25
        )  # units: %
        d2_7n = calc_design2_steric_electronic_metrics(
            mol_rel, pos_rel, syms_rel, c_idx_rel, radius_d2_ang=float(cstate["r_probe_d2_ang"]), grid_step_ang=0.28
        )
        form_metrics = compute_formation_energy_metrics(float(p_data["e_eq"]), syms_rel, c_deg)
        td_s_dimer_373 = compute_stat_mech_tdelta_s_dimer(pos_rel, syms_rel, temperature_k=373.15)  # units: kcal/mol
        dg_dimer_373 = round(float(p_data["delta_e_dimer_kcal_mol"]) + td_s_dimer_373, 2)  # units: kcal/mol

        dt = time.time() - t0
        rec = {
            "id": eid,
            "name": item["name"],
            "smiles": smi,
            "formula": formula,
            "target_type": item["target_type"],
            "center_atom": c_sym,
            "center_atom_p2": p_data["center_atom"],
            "is_d_block": bool(topo["is_d_block"]),
            "is_fallback": bool(p_data["is_fallback"]),
            "mono_relax_steps": int(p_data["mono_relax_steps"]),
            "mono_fmax_ev_ang": float(p_data["mono_fmax_ev_ang"]),  # units: eV/A
            "dimer_init_gap_ang": float(p_data["dimer_init_gap_ang"]),  # units: A
            "dimer_relax_steps": int(p_data["dimer_relax_steps"]),
            "bde_mode": p_data["bde_mode"],
            "e_pot_ev": float(p_data["e_eq"]),  # units: eV
            "e_coh_ev_per_atom": float(p_data["e_coh"]),  # units: eV/atom
            "bde_min_kcal_mol": float(p_data["bde_min"]),  # units: kcal/mol
            "vbur_percent": float(p_data["vbur_percent"]),  # units: % (3.5 A sphere on 7net-0 coords)
            "vbur_rprobe_7net": float(vbur_rprobe_7n),  # units: % (R_probe sphere on 7net-0 coords)
            "vbur_d2_7net": float(d2_7n["vbur_d2_percent"]),  # units: % (Design 2 R_probe^D2 on 7net-0 coords)
            "void_cone_deg_7net": float(d2_7n["void_cone_deg"]),  # units: deg
            "delta_e_dimer_kcal_mol": float(p_data["delta_e_dimer_kcal_mol"]),  # units: kcal/mol
            "td_s_dimer_373_kcal_mol": float(td_s_dimer_373),  # units: kcal/mol
            "delta_g_dimer_373_kcal_mol": dg_dimer_373,  # units: kcal/mol
            "delta_e_complex_kcal_mol": float(p_data["delta_e_complex_kcal_mol"]),  # units: kcal/mol
            "delta_e_form_ev_per_atom": float(form_metrics["delta_e_form_ev_per_atom"]),  # units: eV/atom
            "delta_e_form_kcal_mol": float(form_metrics["delta_e_form_kcal_mol"]),  # units: kcal/mol
            "delta_e_form_per_bond_kcal_mol": float(form_metrics["delta_e_form_per_bond_kcal_mol"]),  # units: kcal/mol
            "radius_of_gyration_ang": float(p_data["rg"]),  # units: A
            "mol_volume_ang3": float(p_data["vol_mol"]),  # units: A^3
        }
        results.append(rec)
        results_by_id[eid] = rec
        print(
            f"[{idx:02d}/{len(candidates):02d}] {eid} | {item['name']:22s} | "
            f"VburD2:{rec['vbur_d2_7net']:5.1f}% | BDE:{rec['bde_min_kcal_mol']:6.1f} | "
            f"ecoh:{rec['e_coh_ev_per_atom']:5.3f} | dEdim:{rec['delta_e_dimer_kcal_mol']:6.1f} | "
            f"dEf/vM:{rec['delta_e_form_per_bond_kcal_mol']:6.1f} | fb:{rec['is_fallback']} | {dt:.1f}s",
            flush=True,
        )

        ordered = [results_by_id[c["id"]] for c in candidates if c["id"] in results_by_id]
        _write_payload(ordered, time.time() - t_start)

    _write_payload(results, time.time() - t_start)
    print(f"Completed SevenNet-0 calculation for {len(results)} candidates -> {OUT_PATH}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--threads", type=int, default=6, help="PyTorch CPU threads")
    parser.add_argument("--no-resume", action="store_true", help="Recompute all entries without cache")
    args = parser.parse_args()
    run_patent_expansion_sevenn(num_threads=args.threads, resume=not args.no_resume)


if __name__ == "__main__":
    main()
