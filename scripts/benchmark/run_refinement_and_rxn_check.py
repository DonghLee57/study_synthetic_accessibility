#!/usr/bin/env python3
"""
Refinement (Items 1 & 4) + Synthesis Reaction Energy Consistency Check (Item 2)
-------------------------------------------------------------------------------
1. Item 1: Element covalent-radius & pi-hapticity adaptive %V_bur lower bound V_min.
2. Item 4: Topological pi_hapticity_index and chelate_ring_count descriptors to
   eliminate false flags on eta-diene/Cp organometallics (Ru(IMBCH), Ru(EBCHD)).
3. Item 2: Evaluate actual stoichiometric synthesis reaction energies (Delta E_syn)
   using SevenNet-0 vs MACE-MP0 across SiCl4 alcoholysis and Ti/Hf salt metathesis.

Outputs:
  - database/benchmark_refinement_and_rxn_check.json
"""

from __future__ import annotations

import json
import math
import os
from typing import Any, Dict, List

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE, LBFGS
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CENTER_ATOMIC_NUMS = {13, 14, 22, 27, 31, 33, 40, 42, 44, 45, 49, 50, 72, 73, 74, 77}
TM_ATOMIC_NUMS = {22, 27, 40, 42, 44, 45, 72, 73, 74, 77}
METAL_CENTER_SET = {"Ti", "Zr", "Hf", "Al", "Ru", "Co", "W", "Mo", "Ta", "Sn", "Cu", "Si"}
EV_TO_KCAL_MOL = 23.06054783


def sigmoid(x: float) -> float:
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    z = math.exp(x)
    return z / (1.0 + z)


def compute_topology_descriptors(smiles: str) -> Dict[str, Any]:
    """Item 4: Compute pi_hapticity_index, chelate_ring_count, and refined coordination flags."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    pi_hapticity_index = 0
    chelate_ring_count = 0

    for ring in atom_rings:
        ring_set = set(ring)
        has_center = any(mol.GetAtomWithIdx(i).GetAtomicNum() in CENTER_ATOMIC_NUMS for i in ring)
        if has_center:
            has_hetero_donor = any(
                mol.GetAtomWithIdx(i).GetAtomicNum() in (7, 8, 15, 16) for i in ring
            )
            if has_hetero_donor:
                chelate_ring_count += 1

    # Count metal-bonded carbons that are in a ring containing sp2/aromatic carbons or unsaturation
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in CENTER_ATOMIC_NUMS:
            continue
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 6 and nb.IsInRing():
                if nb.GetIsAromatic() or nb.GetHybridization() == Chem.HybridizationType.SP2:
                    pi_hapticity_index += 1
                else:
                    # Check if this ring carbon shares a ring with any sp2/aromatic carbon (allylic/eta-diene)
                    nb_idx = nb.GetIdx()
                    for ring in atom_rings:
                        if nb_idx in ring:
                            if any(
                                mol.GetAtomWithIdx(ri).GetAtomicNum() == 6
                                and (
                                    mol.GetAtomWithIdx(ri).GetIsAromatic()
                                    or mol.GetAtomWithIdx(ri).GetHybridization() == Chem.HybridizationType.SP2
                                )
                                for ri in ring
                            ):
                                pi_hapticity_index += 1
                                break

    tm_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() in TM_ATOMIC_NUMS]
    is_homoleptic_tm_sp3 = (
        pi_hapticity_index == 0
        and any(
            len(a.GetNeighbors()) > 0
            and all(
                nb.GetAtomicNum() == 6 and nb.GetHybridization() == Chem.HybridizationType.SP3
                for nb in a.GetNeighbors()
            )
            for a in tm_atoms
        )
    )

    patt_oh = Chem.MolFromSmarts("[#13,#14,#22,#27,#31,#33,#40,#42,#44,#45,#49,#50,#72,#73,#74,#77]-[OX2H1]")
    metal_oh_count = len(mol.GetSubstructMatches(patt_oh)) if patt_oh else 0
    aro_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    coord_flagged = is_homoleptic_tm_sp3 or (metal_oh_count >= 2) or (len(tm_atoms) > 0 and aro_rings >= 5)

    return {
        "pi_hapticity_index": pi_hapticity_index,
        "chelate_ring_count": chelate_ring_count,
        "is_homoleptic_tm_sp3_refined": is_homoleptic_tm_sp3,
        "metal_oh_count": metal_oh_count,
        "aromatic_rings": aro_rings,
        "coord_flagged": coord_flagged,
    }


def evaluate_adaptive_tmc_sa_v2(m: Dict[str, Any], topo: Dict[str, Any]) -> Dict[str, Any]:
    """Item 1: Element covalent-radius & pi-hapticity adaptive V_min for %V_bur."""
    bde = float(m.get("bde_min_kcal_mol", 60.0))
    ecoh = float(m.get("e_coh_ev_per_atom", 4.10))
    vbur = float(m.get("vbur_percent", 70.0))
    d_dim = float(m.get("delta_e_dimer_kcal_mol", 10.0))
    center = str(m.get("center_atom", "C"))
    has_metal_center = center in METAL_CENTER_SET

    # Adaptive V_min: Si (small R_cov = 1.11 A) -> 34.0%; pi-haptic ring complexes -> 34.0%; standard d/p-metal -> 44.0%
    if center == "Si" or topo["pi_hapticity_index"] > 0:
        v_min = 34.0
    else:
        v_min = 44.0

    pen_bde = 3.5 * sigmoid((45.0 - bde) / 2.0)
    pen_coh = 3.5 * sigmoid((3.10 - ecoh) / 0.12)
    pen_labile = max(pen_bde, pen_coh)

    if has_metal_center:
        pen_crowd = 3.5 * sigmoid((vbur - 83.5) / 1.0)
        pen_olig = 3.5 * sigmoid((v_min - vbur) / 2.0) * sigmoid((-1.5 - d_dim) / 1.0)
    else:
        pen_crowd = 0.0
        pen_olig = 0.0

    phys_pen_v2 = round(pen_labile + pen_crowd + pen_olig, 3)

    # Same ligand SA + general SMARTS chem flag from Phase 2, plus refined coordination gate
    base_sa = float(m.get("ligand_sa", m.get("br_sascore", 2.0)))
    if m.get("target_type") == "Inhibitor" and not has_metal_center:
        base_sa = float(m.get("br_sascore", base_sa))
    tmc_sa_v2 = round(min(10.0, max(1.0, base_sa + phys_pen_v2)), 3)

    # Check general SMARTS reactivity flag (-COOH or pure alkane)
    mol = Chem.MolFromSmiles(m["smiles"])
    patt_cooh = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    has_cooh = bool(mol.HasSubstructMatch(patt_cooh)) if patt_cooh else False
    atom_syms = {a.GetSymbol() for a in mol.GetAtoms()}
    is_alkane = (atom_syms <= {"C", "H"}) and (rdMolDescriptors.CalcNumRings(mol) == 0)

    passed_v2 = (tmc_sa_v2 <= 4.0) and (not has_cooh) and (not is_alkane) and (not topo["coord_flagged"])
    return {
        "v_min_used": v_min,
        "phys_pen_v2": phys_pen_v2,
        "tmc_sa_v2": tmc_sa_v2,
        "passed_v2": passed_v2,
    }


def _smiles_to_ase_relaxed(smiles: str, calc: Any) -> float:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=42)
    if AllChem.UFFHasAllMoleculeParams(mol_h):
        try:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=200)
        except Exception:
            pass
    conf = mol_h.GetConformer()
    symbols = [a.GetSymbol() for a in mol_h.GetAtoms()]
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_h.GetNumAtoms())])
    atoms = Atoms(symbols=symbols, positions=pos, cell=[40.0, 40.0, 40.0], pbc=False)
    atoms.center()
    atoms.calc = calc
    FIRE(atoms, logfile=None).run(fmax=0.2, steps=40)
    LBFGS(atoms, maxstep=0.05, logfile=None).run(fmax=0.05, steps=80)
    return float(atoms.get_potential_energy())


def _licl_rocksalt_energy_per_fu(calc: Any) -> float:
    """Compute relaxed rock-salt LiCl(s) energy per LiCl formula unit (eV)."""
    licl_bulk = bulk("LiCl", "rocksalt", a=5.14)
    licl_bulk.calc = calc
    ucf = FrechetCellFilter(licl_bulk)
    LBFGS(ucf, maxstep=0.05, logfile=None).run(fmax=0.05, steps=40)
    return float(licl_bulk.get_potential_energy()) / (len(licl_bulk) / 2.0)


def check_item2_reaction_energies() -> Dict[str, Any]:
    """Item 2: Compare SevenNet-0 vs MACE-MP0 on balanced synthesis reactions."""
    from sevenn.calculator import SevenNetCalculator

    calcs: Dict[str, Any] = {"SevenNet-0": SevenNetCalculator(model="7net-0", device="cpu")}
    mace_status = "unavailable"
    try:
        from mace.calculators import mace_mp
        calcs["MACE-MP0"] = mace_mp(model="small", dispersion=False, default_dtype="float32", device="cpu")
        mace_status = "MACE-MP0 (small) loaded"
    except Exception as e:
        mace_status = f"MACE-MP0 load failed: {e}"

    reactions = [
        {
            "id": "SYN_TEOS",
            "name": "SiCl4 + 4 EtOH -> Si(OEt)4 + 4 HCl",
            "reactants": [("Cl[Si](Cl)(Cl)Cl", 1), ("CCO", 4)],
            "products": [("CCO[Si](OCC)(OCC)OCC", 1), ("Cl", 4)],
        },
        {
            "id": "SYN_TDMAT",
            "name": "TiCl4 + 4 LiNMe2 -> Ti(NMe2)4 + 4 LiCl(gas)",
            "reactants": [("Cl[Ti](Cl)(Cl)Cl", 1), ("[Li]N(C)C", 4)],
            "products": [("CN(C)[Ti](N(C)C)(N(C)C)N(C)C", 1), ("[Li]Cl", 4)],
        },
        {
            "id": "SYN_TEMAH",
            "name": "HfCl4 + 4 LiNMeEt -> Hf(NMeEt)4 + 4 LiCl(gas)",
            "reactants": [("Cl[Hf](Cl)(Cl)Cl", 1), ("[Li]N(C)CC", 4)],
            "products": [("CCN(C)[Hf](N(C)CC)(N(C)CC)N(C)CC", 1), ("[Li]Cl", 4)],
        },
    ]

    unique_smis = sorted({s for rxn in reactions for side in (rxn["reactants"], rxn["products"]) for s, _ in side})
    rxn_results: Dict[str, Any] = {"mace_status": mace_status, "reactions": []}

    energies: Dict[str, Dict[str, float]] = {cname: {} for cname in calcs}
    licl_solid_e: Dict[str, float] = {}
    for cname, calc in calcs.items():
        for smi in unique_smis:
            energies[cname][smi] = _smiles_to_ase_relaxed(smi, calc)
        licl_solid_e[cname] = _licl_rocksalt_energy_per_fu(calc)

    for rxn in reactions:
        row: Dict[str, Any] = {"id": rxn["id"], "name": rxn["name"], "delta_e_kcal_mol": {}}
        has_licl = any(smi == "[Li]Cl" for smi, _ in rxn["products"])
        if has_licl:
            row["delta_e_with_solid_licl_kcal_mol"] = {}
        for cname in calcs:
            e_r = sum(coef * energies[cname][smi] for smi, coef in rxn["reactants"])
            e_p = sum(coef * energies[cname][smi] for smi, coef in rxn["products"])
            row["delta_e_kcal_mol"][cname] = round((e_p - e_r) * EV_TO_KCAL_MOL, 2)
            if has_licl:
                e_p_solid = sum(
                    coef * (licl_solid_e[cname] if smi == "[Li]Cl" else energies[cname][smi])
                    for smi, coef in rxn["products"]
                )
                row["delta_e_with_solid_licl_kcal_mol"][cname] = round((e_p_solid - e_r) * EV_TO_KCAL_MOL, 2)
        rxn_results["reactions"].append(row)

    return rxn_results


def main() -> None:
    p2_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json")
    exp_db_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json")
    with open(p2_path, "r", encoding="utf-8") as f:
        p2_rows = json.load(f)["molecules_evaluated"]
    with open(exp_db_path, "r", encoding="utf-8") as f:
        exp_prec = json.load(f)["precursors_102"]["records"]

    # Items 1 & 4 on 55 audited molecules
    v2_counts = {"Precursor": [0, 0], "Inhibitor": [0, 0], "Negative Control": [0, 0]}
    improved_precursors = []
    for m in p2_rows:
        topo = compute_topology_descriptors(m["smiles"])
        v2 = evaluate_adaptive_tmc_sa_v2(m, topo)
        tt = m["target_type"]
        v2_counts[tt][1] += 1
        if v2["passed_v2"]:
            v2_counts[tt][0] += 1
        if tt == "Precursor" and (not m["tuned_with_chem_pass"]) and v2["passed_v2"]:
            improved_precursors.append({"id": m["id"], "name": m["name"], "old_sa": m["tmc_sa_phys"], "new_sa": v2["tmc_sa_v2"], "v_min": v2["v_min_used"]})

    # Item 4 on 75 expanded Ir/Rh/As/Ga/In/Ru precursors
    baseline_smis = {m["smiles"] for m in p2_rows if m["target_type"] == "Precursor"}
    exp_75 = [r for r in exp_prec if r["smiles"] not in baseline_smis]
    exp_75_flagged = []
    pi_haptic_75 = 0
    chelate_75 = 0
    for r in exp_75:
        topo = compute_topology_descriptors(r["smiles"])
        if topo["pi_hapticity_index"] > 0:
            pi_haptic_75 += 1
        if topo["chelate_ring_count"] > 0:
            chelate_75 += 1
        if topo["coord_flagged"]:
            exp_75_flagged.append(r["name"])

    # Item 2 reaction energy check
    rxn_check = check_item2_reaction_energies()

    out = {
        "items_1_and_4": {
            "cohort_55_v2": {
                "precursor_pass": f"{v2_counts['Precursor'][0]}/{v2_counts['Precursor'][1]}",
                "precursor_tpr_pct": round(v2_counts["Precursor"][0] / v2_counts["Precursor"][1] * 100.0, 1),
                "inhibitor_pass": f"{v2_counts['Inhibitor'][0]}/{v2_counts['Inhibitor'][1]}",
                "inhibitor_tpr_pct": round(v2_counts["Inhibitor"][0] / v2_counts["Inhibitor"][1] * 100.0, 1),
                "neg_control_pass": f"{v2_counts['Negative Control'][0]}/{v2_counts['Negative Control'][1]}",
                "neg_control_tnr_pct": round((v2_counts["Negative Control"][1] - v2_counts["Negative Control"][0]) / v2_counts["Negative Control"][1] * 100.0, 1),
            },
            "newly_rescued_precursors": improved_precursors,
            "expanded_75_topology": {
                "total": len(exp_75),
                "pi_haptic_complexes": pi_haptic_75,
                "chelate_complexes": chelate_75,
                "coord_flagged_count": len(exp_75_flagged),
                "unflagged_count": len(exp_75) - len(exp_75_flagged),
                "coord_flagged_names": exp_75_flagged,
            },
        },
        "item_2_reaction_energy_check": rxn_check,
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_refinement_and_rxn_check.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print("Saved ->", out_path)
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
