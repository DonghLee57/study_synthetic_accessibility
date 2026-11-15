#!/usr/bin/env python3
"""
Reusable Semiconductor Precursor Descriptors, Coordination Gates & TMC-SA Scorer
---------------------------------------------------------------------------------
Shared computational core consolidating:
  - 2D physicochemical + SMARTS coordination/reactivity descriptors (no ID/label leakage)
  - Intrinsic periodic-property + Design 2 multi-state coordination/topological gate
  - Pro-ligand separation & SA scoring (`compute_ligand_sa`)
  - Continuous Fermi-Dirac physics-coupled TMC-SA scoring (`compute_tmc_sa_score`)
  - Design 1 + Design 2 3D geometry & steric-electronic descriptor evaluation
"""

from __future__ import annotations

import math
import os
import sys
from typing import Any, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, RDConfig, rdMolDescriptors

sys.path.insert(0, os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore

from element_physics import (
    calc_buried_volume_sphere,
    calc_design2_steric_electronic_metrics,
    compute_coordination_state,
    compute_refined_pi_hapticity,
    compute_stat_mech_tdelta_s_dimer,
    embed_etkdg_3d,
    find_central_atom_intrinsic,
    get_element_properties,
    is_coordination_center_atom,
    is_d_block_atomic_num,
)

# General physical/chemical SMARTS patterns (non-target-ID specific)
PATT_BETA_H_ALKYL = Chem.MolFromSmarts("[#22,#24,#27,#40,#42,#44,#45,#72,#73,#74,#77]-[CH2,CH1]-[CH3,CH2,CH1]")
PATT_METAL_OH = Chem.MolFromSmarts("[#13,#14,#22,#24,#27,#31,#33,#40,#42,#44,#45,#49,#50,#72,#73,#74,#77]-[OX2H1]")
PATT_CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def sigmoid(x: float) -> float:
    """Numerically stable logistic Fermi-Dirac barrier mapping."""
    if x >= 40.0:
        return 1.0
    if x <= -40.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def compute_chemical_coord_descriptors(smiles: str) -> Dict[str, float]:
    """
    Compute 17 RDKit physicochemical, SMARTS reactivity, and center-coordination
    descriptors without element-ID hardcoding.

    Physics/Algorithm:
        Captures molecular volatility proxies (MW, logP, TPSA, HBD/HBA), metal-to-ligand
        sigma vs. heteroatom bond fractions, homoleptic sp3 alkyl instability, saturated
        alkane inertness, carboxylic protonolysis sites, and strained 3/4-membered
        unsaturated silacycles.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"

    mw = float(Descriptors.MolWt(mol))  # units: g/mol
    logp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))  # units: A^2
    hbd = float(Lipinski.NumHDonors(mol))
    hba = float(Lipinski.NumHAcceptors(mol))
    rotb = float(Lipinski.NumRotatableBonds(mol))
    heavy = float(mol.GetNumHeavyAtoms())
    aro_rings = float(rdMolDescriptors.CalcNumAromaticRings(mol))
    rings = float(rdMolDescriptors.CalcNumRings(mol))

    center_atoms = [a for a in mol.GetAtoms() if is_coordination_center_atom(a.GetSymbol(), a.GetAtomicNum())]
    n_centers = len(center_atoms)
    metal_c_bonds = 0
    metal_hetero_bonds = 0
    total_center_bonds = 0
    for ca in center_atoms:
        for nb in ca.GetNeighbors():
            total_center_bonds += 1
            z = nb.GetAtomicNum()
            if z == 6:
                metal_c_bonds += 1
            elif z in (7, 8, 9, 15, 16, 17, 35, 53):
                metal_hetero_bonds += 1

    metal_c_sigma_frac = float(metal_c_bonds / total_center_bonds) if total_center_bonds > 0 else 0.0
    metal_hetero_frac = float(metal_hetero_bonds / total_center_bonds) if total_center_bonds > 0 else 0.0
    is_homoleptic_mc = 1.0 if (n_centers > 0 and total_center_bonds > 0 and metal_c_bonds == total_center_bonds) else 0.0

    tm_atoms = [a for a in mol.GetAtoms() if is_d_block_atomic_num(a.GetAtomicNum())]
    is_homoleptic_tm_sp3 = 1.0 if any(
        len(a.GetNeighbors()) > 0
        and all(nb.GetAtomicNum() == 6 and nb.GetHybridization() == Chem.HybridizationType.SP3 for nb in a.GetNeighbors())
        for a in tm_atoms
    ) else 0.0

    beta_h_count = float(len(mol.GetSubstructMatches(PATT_BETA_H_ALKYL))) if PATT_BETA_H_ALKYL else 0.0
    metal_oh_count = float(len(mol.GetSubstructMatches(PATT_METAL_OH))) if PATT_METAL_OH else 0.0
    cooh_count = float(len(mol.GetSubstructMatches(PATT_CARBOXYLIC_ACID))) if PATT_CARBOXYLIC_ACID else 0.0

    atom_syms = {a.GetSymbol() for a in mol.GetAtoms()}
    all_single = all(b.GetBondType() == Chem.BondType.SINGLE for b in mol.GetBonds())
    is_inert_alkane = 1.0 if (atom_syms <= {"C", "H"} and all_single and mol.GetNumBonds() > 0) else 0.0

    ri = mol.GetRingInfo()
    strained_si_ring = 0.0
    for ring in ri.AtomRings():
        if len(ring) <= 4 and "Si" in {mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring}:
            if any(
                b.GetBondType() != Chem.BondType.SINGLE
                for b in mol.GetBonds()
                if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring
            ):
                strained_si_ring = 1.0

    return {
        "mw": mw,
        "logp": logp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rot_bonds": rotb,
        "heavy_atoms": heavy,
        "aromatic_rings": aro_rings,
        "ring_count": rings,
        "metal_c_sigma_frac": metal_c_sigma_frac,
        "metal_hetero_frac": metal_hetero_frac,
        "is_homoleptic_mc": is_homoleptic_mc,
        "is_homoleptic_tm_sp3": is_homoleptic_tm_sp3,
        "beta_h_alkyl_count": beta_h_count,
        "metal_oh_count": metal_oh_count,
        "cooh_count": cooh_count,
        "is_inert_alkane": is_inert_alkane,
        "strained_si_ring": strained_si_ring,
    }


def evaluate_coord_chem_gate(
    smiles: str,
    sa_score: float,
    desc: Dict[str, float],
    topo: Dict[str, Any],
    include_design2: bool = False,
) -> Dict[str, Any]:
    """
    Evaluate uniform coordination, reactivity, and topology gating rules alongside
    the synthetic accessibility threshold ``TMC-SA <= 4.0``.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"
    has_tm = bool(topo.get("is_d_block")) or any(is_d_block_atomic_num(a.GetAtomicNum()) for a in mol.GetAtoms())

    reasons: List[str] = []
    if sa_score > 4.0:
        reasons.append("sa_gt_4")
    if desc["cooh_count"] > 0:
        reasons.append("carboxylic_acid")
    if desc["is_inert_alkane"] > 0:
        reasons.append("inert_alkane")
    if topo["is_homoleptic_tm_sp3_refined"] or (desc["is_homoleptic_tm_sp3"] > 0 and topo["pi_hapticity_index"] == 0):
        reasons.append("homoleptic_tm_sp3")
    if desc["metal_oh_count"] >= 2:
        reasons.append("poly_hydroxide")
    if has_tm and desc["aromatic_rings"] >= 3 and desc["metal_c_sigma_frac"] >= 0.9:
        reasons.append("overcrowded_poly_aryl_tm")
    if desc["strained_si_ring"] > 0:
        reasons.append("strained_si_ring")
    if include_design2 and topo.get("catenated_main_group_hydride"):
        reasons.append("catenated_main_group_hydride")
    if include_design2 and topo.get("pi_conjugated_main_group_hydride"):
        reasons.append("pi_conjugated_main_group_hydride")

    passed = len(reasons) == 0
    return {"passed": passed, "rejected": not passed, "reasons": reasons}


def compute_coord_risk_count(desc: Dict[str, float], topo: Dict[str, Any], is_coord_center: bool) -> int:
    """Count independent intrinsic chemical/coordination failure flags (excluding ``sa_gt_4``)."""
    return (
        int(desc["cooh_count"] > 0)
        + int(desc["is_inert_alkane"] > 0)
        + int(topo["is_homoleptic_tm_sp3_refined"] or (desc["is_homoleptic_tm_sp3"] > 0 and topo["pi_hapticity_index"] == 0))
        + int(desc["metal_oh_count"] >= 2)
        + int(desc["aromatic_rings"] >= 3 and desc["metal_c_sigma_frac"] >= 0.9 and is_coord_center)
        + int(desc["strained_si_ring"] > 0)
        + int(bool(topo.get("catenated_main_group_hydride")))
        + int(bool(topo.get("pi_conjugated_main_group_hydride")))
    )


def compute_ligand_sa(
    smi: str,
    sa_val: float,
    br_val: float,
    br_scorer: Any = None,
) -> float:
    """
    Decompose coordination complex into detached organic pro-ligands and score the
    rate-limiting ligand fragment with ``min(BR-SAScore, SAScore)``.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return sa_val
    idxs = sorted(
        [a.GetIdx() for a in mol.GetAtoms() if is_coordination_center_atom(a.GetSymbol(), a.GetAtomicNum())],
        reverse=True,
    )
    if idxs:
        rw = Chem.RWMol(mol)
        for c_idx in idxs:
            rw.RemoveAtom(c_idx)
        frags = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
        lig_scores: List[float] = []
        for frag in frags:
            try:
                for a in frag.GetAtoms():
                    a.SetNumRadicalElectrons(0)
                    a.SetNoImplicit(False)
                    a.SetNumExplicitHs(0)
                Chem.SanitizeMol(frag)
                l_smi = Chem.MolToSmiles(frag)
                if l_smi and l_smi not in ("[H]", "Cl", "O"):
                    cand: List[float] = []
                    if br_scorer is not None:
                        try:
                            bs, _ = br_scorer.calculateScore(l_smi)
                            cand.append(float(bs))
                        except Exception:
                            pass
                    l_mol = Chem.MolFromSmiles(l_smi)
                    if l_mol is not None:
                        cand.append(float(sascorer.calculateScore(l_mol)))
                    if cand:
                        lig_scores.append(min(cand))
            except Exception:
                pass
        return round(max(lig_scores), 3) if lig_scores else 1.00
    return round(min(sa_val, br_val) if br_val < 9.9 else sa_val, 3)


def compute_tmc_sa_score(
    base_sa: float,
    vbur_val: float,
    v_min: float,
    is_coord_center: bool,
    bde_val: Optional[float] = None,  # units: kcal/mol
    ecoh_val: Optional[float] = None,  # units: eV/atom
    ddim_val: Optional[float] = None,  # units: kcal/mol
    is_d_block: bool = False,
    ve_deficit: float = 0.0,
    chi_pauling: float = 2.0,
    center_h: int = 0,
) -> float:
    """
    Compute continuous Transition-Metal Complex Synthetic Accessibility (``TMC-SA``)
    combining detached pro-ligand synthetic complexity with smooth Fermi-Dirac physical
    penalties for bond lability, steric overcrowding, and oligomerization.

    Physics/Algorithm:
        $$\text{TMC-SA} = \operatorname{clip}\!\left(S_{\text{base}} + \text{pen}_{\text{labile}} + \text{pen}_{\text{crowd}} + \text{pen}_{\text{olig}} + \text{pen}_{\text{hydride}},\; 1.0,\; 10.0\right)$$
    """
    pen_crowd = 3.5 * sigmoid((vbur_val - 83.5) / 1.0) if is_coord_center else 0.0
    if bde_val is not None and ecoh_val is not None and ddim_val is not None:
        pen_bde = 3.5 * sigmoid((45.0 - bde_val) / 2.0)
        pen_coh = 3.5 * sigmoid((3.10 - ecoh_val) / 0.12)
        pen_labile = max(pen_bde, pen_coh)
        pen_olig = (3.5 * sigmoid((v_min - vbur_val) / 2.0) * sigmoid((-1.5 - ddim_val) / 1.0)) if is_coord_center else 0.0
        return round(min(10.0, max(1.0, base_sa + pen_labile + pen_crowd + pen_olig)), 3)

    hydride_pen = 0.0
    if (not is_d_block) and is_coord_center and (bde_val is not None) and (0.0 <= bde_val < 30.0) and (ve_deficit > 0.0) and (chi_pauling <= 2.05) and (center_h >= 1):
        hydride_pen = 3.5 * sigmoid((35.0 - bde_val) / 2.0)
    return round(min(10.0, max(1.0, base_sa + pen_crowd + hydride_pen)), 3)


def compute_design12_geometry_metrics(smiles: str) -> Dict[str, Any]:
    """
    Generate deterministic ETKDGv3 3D conformer and compute Design 1 (covalent probe sphere)
    and Design 2 (state-adapted effective radius, polarity/hydride-weighted %V_bur^D2,
    open attack void cone, and 373.15 K dimerization entropy) descriptors.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        props = get_element_properties("C", 6)
        return {
            "vbur_r35": 0.0, "vbur_rprobe": 0.0, "vbur_d2_percent": 0.0, "void_cone_deg": 180.0,
            "td_s_dimer_373": 0.0, "center_sym": "C", "center_z": 6, "center_degree": 1, "embed_ok": False,
            "chi_eff": props["chi_pauling"], "delta_chi_mean": 0.0, "ionicity_mean": 0.0,
            "r_eff_ang": props["r_cov_ang"], "r_probe_d2_ang": props["r_probe_ang"],
            "cn_eff": 1.0, "n_ve": 4.0, "ve_deficit": 0.5, "coord_propensity": 0.0,
            **props,
        }

    emb = embed_etkdg_3d(smiles)
    c_idx, c_sym, c_z = find_central_atom_intrinsic(mol)
    props = get_element_properties(c_sym, c_z)
    pi_h = compute_refined_pi_hapticity(mol, c_idx)
    cstate = compute_coordination_state(mol, c_idx, pi_hapticity=pi_h)
    c_deg = max(1, mol.GetAtomWithIdx(c_idx).GetDegree() + mol.GetAtomWithIdx(c_idx).GetTotalNumHs())

    if emb is None:
        return {
            "vbur_r35": 0.0, "vbur_rprobe": 0.0, "vbur_d2_percent": 0.0, "void_cone_deg": 180.0,
            "td_s_dimer_373": 0.0, "center_sym": c_sym, "center_z": c_z, "center_degree": c_deg, "embed_ok": False,
            **cstate, **props,
        }

    mol_h, pos, syms, c_idx_h, c_sym, c_z = emb
    c_deg = max(1, mol_h.GetAtomWithIdx(c_idx_h).GetDegree())
    v35 = calc_buried_volume_sphere(pos, syms, c_idx_h, radius_ang=3.5, grid_step_ang=0.25)
    vrp = calc_buried_volume_sphere(pos, syms, c_idx_h, radius_ang=props["r_probe_ang"], grid_step_ang=0.25)
    d2 = calc_design2_steric_electronic_metrics(
        mol_h, pos, syms, c_idx_h, radius_d2_ang=float(cstate["r_probe_d2_ang"]), grid_step_ang=0.28,
    )
    td_s = compute_stat_mech_tdelta_s_dimer(pos, syms, temperature_k=373.15)
    return {
        "vbur_r35": float(v35),
        "vbur_rprobe": float(vrp),
        "vbur_d2_percent": float(d2["vbur_d2_percent"]),
        "void_cone_deg": float(d2["void_cone_deg"]),
        "td_s_dimer_373": float(td_s),
        "syms": syms,
        "center_sym": c_sym,
        "center_z": c_z,
        "center_degree": c_deg,
        "embed_ok": True,
        **cstate,
        **props,
    }
