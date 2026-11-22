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
from collections import deque
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
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


# Simple inorganic ligands excluded from organic pro-ligand SA complexity scoring
_SIMPLE_INORGANIC_LIGANDS = frozenset({
    "[H]", "H", "F", "Cl", "Br", "I", "O", "S",
    "[F]", "[Cl]", "[Br]", "[I]", "[O]", "[S]",
    "[C-]#[O+]", "[O+]#[C-]", "C#O",
})

# d8 centers preferring square-planar (D4h) geometry at CN = 4
_SQUARE_PLANAR_D8_Z = frozenset({45, 46, 77, 78, 79})  # Rh, Pd, Ir, Pt, Au


def _ideal_coordination_unit_vectors(cn: int, center_z: int) -> np.ndarray:
    """Return ideal unit direction vectors for coordination number ``cn``.

    Physics/Algorithm:
        Generates canonical VSEPR / ligand-field polyhedral directions:
        linear (CN=2), trigonal planar (CN=3), square planar (CN=4, d8) or
        tetrahedral (CN=4), trigonal bipyramidal (CN=5), octahedral (CN=6),
        and deterministic Fibonacci-sphere vectors (CN>=7).

    Returns:
        Array of shape ``(cn, 3)`` with unit-norm vectors (dimensionless).
    """
    if cn <= 1:
        return np.array([[1.0, 0.0, 0.0]], dtype=float)
    if cn == 2:
        return np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float)
    if cn == 3:
        ang = [0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0]
        return np.array([[math.cos(a), math.sin(a), 0.0] for a in ang], dtype=float)
    if cn == 4:
        if center_z in _SQUARE_PLANAR_D8_Z:
            return np.array(
                [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
                dtype=float,
            )
        inv_s3 = 1.0 / math.sqrt(3.0)
        return np.array(
            [
                [inv_s3, inv_s3, inv_s3],
                [inv_s3, -inv_s3, -inv_s3],
                [-inv_s3, inv_s3, -inv_s3],
                [-inv_s3, -inv_s3, inv_s3],
            ],
            dtype=float,
        )
    if cn == 5:
        ang = [0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0]
        vecs = [[math.cos(a), math.sin(a), 0.0] for a in ang]
        vecs.extend([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
        return np.array(vecs, dtype=float)
    if cn == 6:
        return np.array(
            [
                [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0], [0.0, 0.0, -1.0],
            ],
            dtype=float,
        )
    # Deterministic Fibonacci sphere for CN >= 7
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    vecs = []
    for i in range(cn):
        y = 1.0 - (2.0 * i + 1.0) / float(cn)
        r = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden_angle * i
        vecs.append([math.cos(theta) * r, y, math.sin(theta) * r])
    return np.array(vecs, dtype=float)


def build_coordination_template_3d(
    mol: Chem.Mol,
    c_idx: int,
) -> Optional[Tuple[Chem.Mol, np.ndarray, List[str], int, str, int]]:
    """Build deterministic 3D polyhedral coordinates when RDKit ETKDG/UFF fails.

    Physics/Algorithm:
        Places the central coordination atom at the origin, distributes first-shell
        coordinated neighbors along ideal VSEPR/ligand-field unit vectors at
        covalent sum bond lengths ``d = R_cov(M) + R_cov(L)`` (# units: Å), and
        propagates outer-shell ligand atoms outward via BFS (using 1.13 Å for
        C≡O triple bonds and covalent sums otherwise) to avoid steric overlap.

    Args:
        mol: RDKit molecule (heavy-atom graph).
        c_idx: Atom index of the coordination center.

    Returns:
        Tuple ``(mol_h, pos, syms, c_idx_h, c_sym, c_z)`` matching
        ``embed_etkdg_3d``, or ``None`` if invalid.
    """
    if mol is None or c_idx < 0 or c_idx >= mol.GetNumAtoms():
        return None
    mol_h = Chem.AddHs(mol)
    n_atoms = mol_h.GetNumAtoms()
    c_atom = mol_h.GetAtomWithIdx(c_idx)
    c_sym = c_atom.GetSymbol()
    c_z = int(c_atom.GetAtomicNum())
    c_rcov = float(get_element_properties(c_sym, c_z)["r_cov_ang"])  # units: Å

    pos = np.zeros((n_atoms, 3), dtype=float)  # units: Å
    placed = {c_idx}
    neighbors = [nb.GetIdx() for nb in c_atom.GetNeighbors()]
    u_vecs = _ideal_coordination_unit_vectors(len(neighbors), c_z)

    # Queue entries: (atom_idx, outward_unit_vec)
    queue: deque[Tuple[int, np.ndarray]] = deque()
    for nb_idx, u_dir in zip(neighbors, u_vecs):
        nb_atom = mol_h.GetAtomWithIdx(nb_idx)
        nb_rcov = float(get_element_properties(nb_atom.GetSymbol(), int(nb_atom.GetAtomicNum()))["r_cov_ang"])  # units: Å
        bond_len = max(1.10, c_rcov + nb_rcov)  # units: Å
        pos[nb_idx] = u_dir * bond_len
        placed.add(nb_idx)
        queue.append((nb_idx, u_dir))

    while queue:
        parent_idx, parent_dir = queue.popleft()
        p_atom = mol_h.GetAtomWithIdx(parent_idx)
        p_rcov = float(get_element_properties(p_atom.GetSymbol(), int(p_atom.GetAtomicNum()))["r_cov_ang"])  # units: Å
        children = [nb.GetIdx() for nb in p_atom.GetNeighbors() if nb.GetIdx() not in placed]
        if not children:
            continue
        # Orthonormal basis perpendicular to parent_dir
        ref = np.array([0.0, 0.0, 1.0]) if abs(parent_dir[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        e1 = np.cross(parent_dir, ref)
        e1 /= max(1e-9, float(np.linalg.norm(e1)))
        e2 = np.cross(parent_dir, e1)

        n_ch = len(children)
        for k, ch_idx in enumerate(children):
            ch_atom = mol_h.GetAtomWithIdx(ch_idx)
            bond = mol_h.GetBondBetweenAtoms(parent_idx, ch_idx)
            is_triple = (bond is not None) and (bond.GetBondType() == Chem.BondType.TRIPLE)
            ch_rcov = float(get_element_properties(ch_atom.GetSymbol(), int(ch_atom.GetAtomicNum()))["r_cov_ang"])  # units: Å
            b_len = 1.13 if is_triple else max(0.96, p_rcov + ch_rcov)  # units: Å
            if n_ch == 1:
                ch_dir = parent_dir
            else:
                phi = 2.0 * math.pi * k / float(n_ch)
                # Cone half-angle ~60 deg outward
                ch_dir = 0.5 * parent_dir + 0.8660254 * (math.cos(phi) * e1 + math.sin(phi) * e2)
                ch_dir /= max(1e-9, float(np.linalg.norm(ch_dir)))
            pos[ch_idx] = pos[parent_idx] + ch_dir * b_len
            placed.add(ch_idx)
            queue.append((ch_idx, ch_dir))

    # Attach 3D conformer onto mol_h so downstream RDKit geometry queries work
    conf = Chem.Conformer(n_atoms)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, (float(pos[i, 0]), float(pos[i, 1]), float(pos[i, 2])))
    mol_h.RemoveAllConformers()
    mol_h.AddConformer(conf, assignId=True)
    syms = [a.GetSymbol() for a in mol_h.GetAtoms()]
    return mol_h, pos, syms, c_idx, c_sym, c_z


def relax_3d_structure_with_sevenn(
    mol_h: Chem.Mol,
    pos: np.ndarray,
    syms: List[str],
    fmax_ev_ang: float = 0.04,
    max_steps: int = 80,
    sevenn_calc: Any = None,
) -> Tuple[Chem.Mol, np.ndarray, Dict[str, Any]]:
    """Relax 3D molecular coordinates via ASE LBFGS with SevenNet-0.

    Physics/Algorithm:
        Minimizes the SevenNet-0 potential energy surface ``E(R)`` (# units: eV)
        until ``max_i ||-∇_{r_i} E|| < fmax_ev_ang`` (# units: eV/Å) or
        ``max_steps`` iterations, updating the RDKit conformer in-place.

    Args:
        mol_h: RDKit molecule with explicit hydrogens.
        pos: Initial Cartesian coordinates, shape ``(N, 3)`` (# units: Å).
        syms: Element symbols of length ``N``.
        fmax_ev_ang: Force convergence threshold (# units: eV/Å).
        max_steps: Maximum LBFGS optimization steps.
        sevenn_calc: Optional pre-initialized ``SevenNetCalculator`` instance.

    Returns:
        Tuple ``(mol_h, relaxed_pos, relax_meta)``.
    """
    from ase import Atoms
    from ase.optimize import LBFGS
    from sevenn.calculator import SevenNetCalculator

    calc = sevenn_calc if sevenn_calc is not None else SevenNetCalculator(model="7net-0", device="cpu")
    atoms = Atoms(symbols=syms, positions=np.asarray(pos, dtype=float))
    atoms.calc = calc
    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=float(fmax_ev_ang), steps=int(max_steps))

    relaxed_pos = np.asarray(atoms.get_positions(), dtype=float)  # units: Å
    energy_ev = float(atoms.get_potential_energy())  # units: eV
    forces = np.asarray(atoms.get_forces(), dtype=float)  # units: eV/Å
    final_fmax = float(np.max(np.linalg.norm(forces, axis=1))) if len(forces) else 0.0  # units: eV/Å

    conf = Chem.Conformer(len(syms))
    for i in range(len(syms)):
        conf.SetAtomPosition(i, (float(relaxed_pos[i, 0]), float(relaxed_pos[i, 1]), float(relaxed_pos[i, 2])))
    mol_h.RemoveAllConformers()
    mol_h.AddConformer(conf, assignId=True)
    meta = {
        "relax_converged": bool(final_fmax <= fmax_ev_ang),
        "relax_steps": int(opt.nsteps),
        "relax_fmax_ev_ang": round(final_fmax, 5),  # units: eV/Å
        "mono_energy_ev": round(energy_ev, 5),  # units: eV
    }
    return mol_h, relaxed_pos, meta


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
                if l_smi and l_smi not in _SIMPLE_INORGANIC_LIGANDS:
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


def compute_ensemble_sapis_metrics(
    base_sa: float,
    vbur_val: float,
    v_min: float,
    is_coord_center: bool,
    model_metrics: Dict[str, Dict[str, Optional[float]]],
    is_d_block: bool = False,
    ve_deficit: float = 0.0,
    chi_pauling: float = 2.0,
    center_h: int = 0,
    sigma_high_kcal: float = 5.0,
    sigma_med_kcal: float = 12.0,
) -> Dict[str, Any]:
    """Aggregate multi-uMLIP single-point physical metrics into a robust SAPIS ensemble score and UQ grade.

    Physics/Algorithm:
        Given per-model physical descriptors ``{model_name: {'bde_min', 'e_coh', 'de_dimer'}}``
        evaluated on a shared relaxed 3D geometry (e.g., 4-uMLIP ensemble of
        ``SevenNet-0``, ``MACE-OMol``, ``ORB-v3-omol``, and ``FAIRChem UMA-S``):
          1. Computes outlier-robust median physical estimates
             ``BDE_hat = median_m(BDE_min^(m))``, ``E_coh_hat = median_m(E_coh^(m))``,
             ``dE_dimer_hat = median_m(dE_dimer^(m))`` (# units: kcal/mol, eV/atom, kcal/mol)
             so a single-model crystal-bias underbinding outlier (e.g., SevenNet-0 on Mo(CO)6)
             does not corrupt the representative score.
          2. Evaluates representative ``sapis_median`` and per-model scores to report the
             conservative uncertainty interval ``[sapis_min, sapis_max]``.
          3. Quantifies epistemic uncertainty via ``sigma_bde = std_m(BDE_min^(m))`` and
             ``iqr_bde = Q75 - Q25`` (# units: kcal/mol), assigning a 3-tier confidence grade:
             ``High`` (<= 5.0), ``Medium`` (5.0-12.0), or ``Low (DFT Flag)`` (> 12.0).

    Args:
        base_sa: Track 1 ligand/organic synthetic accessibility score (1.0-10.0).
        vbur_val: Buried volume percentage (%V_bur).
        v_min: Adaptive minimum buried volume threshold (%).
        is_coord_center: Whether the central atom is a coordination center.
        model_metrics: Mapping from uMLIP model name to metric dict
            (``bde_min`` in kcal/mol, ``e_coh`` in eV/atom, ``de_dimer`` in kcal/mol).
        is_d_block: Whether the center is a d-block transition metal.
        ve_deficit: Valence electron deficit (dimensionless).
        chi_pauling: Pauling electronegativity of the central atom.
        center_h: Number of hydrogens directly bonded to the central atom.
        sigma_high_kcal: Upper std threshold for High confidence (# units: kcal/mol).
        sigma_med_kcal: Upper std threshold for Medium confidence (# units: kcal/mol).

    Returns:
        Dictionary with median physical metrics, representative ``sapis_median``,
        conservative interval ``[sapis_min, sapis_max]``, per-model scores,
        ``sigma_bde_kcal_mol``, ``iqr_bde_kcal_mol``, ``confidence_grade``,
        and ``dft_verification_flag``.
    """
    bdes: List[float] = []
    ecohs: List[float] = []
    ddims: List[float] = []
    per_model_sapis: Dict[str, float] = {}

    for m_name, m_vals in model_metrics.items():
        b_v = m_vals.get("bde_min")
        e_v = m_vals.get("e_coh")
        d_v = m_vals.get("de_dimer")
        if b_v is not None:
            bdes.append(float(b_v))
        if e_v is not None:
            ecohs.append(float(e_v))
        if d_v is not None:
            ddims.append(float(d_v))
        per_model_sapis[m_name] = compute_tmc_sa_score(
            base_sa=base_sa,
            vbur_val=vbur_val,
            v_min=v_min,
            is_coord_center=is_coord_center,
            bde_val=float(b_v) if b_v is not None else None,
            ecoh_val=float(e_v) if e_v is not None else None,
            ddim_val=float(d_v) if d_v is not None else None,
            is_d_block=is_d_block,
            ve_deficit=ve_deficit,
            chi_pauling=chi_pauling,
            center_h=center_h,
        )

    bde_med = float(np.median(bdes)) if bdes else None  # units: kcal/mol
    ecoh_med = float(np.median(ecohs)) if ecohs else None  # units: eV/atom
    ddim_med = float(np.median(ddims)) if ddims else None  # units: kcal/mol

    sapis_med = compute_tmc_sa_score(
        base_sa=base_sa,
        vbur_val=vbur_val,
        v_min=v_min,
        is_coord_center=is_coord_center,
        bde_val=bde_med,
        ecoh_val=ecoh_med,
        ddim_val=ddim_med,
        is_d_block=is_d_block,
        ve_deficit=ve_deficit,
        chi_pauling=chi_pauling,
        center_h=center_h,
    )

    scores_list = list(per_model_sapis.values())
    sapis_min = round(min(scores_list), 3) if scores_list else sapis_med
    sapis_max = round(max(scores_list), 3) if scores_list else sapis_med
    sapis_spread = round(sapis_max - sapis_min, 3)
    sigma_sapis = round(float(np.std(scores_list)), 3) if len(scores_list) >= 2 else 0.0

    if len(bdes) >= 2:
        sigma_bde = float(np.std(bdes))  # units: kcal/mol
        q75, q25 = np.percentile(bdes, [75, 25])
        iqr_bde = float(q75 - q25)  # units: kcal/mol
    else:
        sigma_bde = 0.0
        iqr_bde = 0.0

    n_m = len(scores_list)
    n_pass = sum(1 for v in scores_list if v <= 4.0)
    n_agree = max(n_pass, n_m - n_pass)

    # Decision-aware 3-tier confidence grading:
    # - High: unanimous verdict (N_agree == M) AND tight score spread (delta_SAPIS <= 0.30)
    # - Low (DFT Flag): split verdict (N_agree <= M // 2, e.g. 2:2 Pass vs Fail)
    # - Medium: majority consensus (e.g. 3/4 with 1-model outlier) OR unanimous with delta_SAPIS > 0.30
    if n_m >= 2 and n_agree == n_m and sapis_spread <= 0.30:
        conf_grade = "High"
        dft_flag = False
    elif n_m >= 2 and n_agree <= (n_m // 2):
        conf_grade = "Low (DFT Flag)"
        dft_flag = True
    else:
        conf_grade = "Medium"
        dft_flag = False

    return {
        "n_models": len(model_metrics),
        "sapis_median": sapis_med,
        "sapis_min": sapis_min,
        "sapis_max": sapis_max,
        "sapis_spread": sapis_spread,
        "verdict_agreement": f"{n_agree}/{n_m}" if n_m else "0/0",
        "sigma_sapis": sigma_sapis,
        "per_model_sapis": per_model_sapis,
        "bde_min_median_kcal_mol": round(bde_med, 3) if bde_med is not None else None,
        "e_coh_median_ev_atom": round(ecoh_med, 4) if ecoh_med is not None else None,
        "de_dimer_median_kcal_mol": round(ddim_med, 3) if ddim_med is not None else None,
        "sigma_bde_kcal_mol": round(sigma_bde, 3),
        "iqr_bde_kcal_mol": round(iqr_bde, 3),
        "confidence_grade": conf_grade,
        "dft_verification_flag": dft_flag,
    }


def compute_design12_geometry_metrics(
    smiles: str,
    relax_with_sevenn: bool = False,
    sevenn_calc: Any = None,
    fmax_ev_ang: float = 0.04,
    max_steps: int = 80,
) -> Dict[str, Any]:
    """
    Generate deterministic ETKDGv3 (or polyhedral coordination template fallback) 3D conformer,
    optionally relax with SevenNet-0 LBFGS, and compute Design 1 (covalent probe sphere)
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
    embed_source = "etkdg_3d"
    c_idx, c_sym, c_z = find_central_atom_intrinsic(mol)
    props = get_element_properties(c_sym, c_z)
    pi_h = compute_refined_pi_hapticity(mol, c_idx)
    cstate = compute_coordination_state(mol, c_idx, pi_hapticity=pi_h)
    c_deg = max(1, mol.GetAtomWithIdx(c_idx).GetDegree() + mol.GetAtomWithIdx(c_idx).GetTotalNumHs())

    if emb is None:
        emb = build_coordination_template_3d(mol, c_idx)
        embed_source = "coord_template_3d"

    if emb is None:
        return {
            "vbur_r35": 0.0, "vbur_rprobe": 0.0, "vbur_d2_percent": 0.0, "void_cone_deg": 180.0,
            "td_s_dimer_373": 0.0, "center_sym": c_sym, "center_z": c_z, "center_degree": c_deg,
            "embed_ok": False, "embed_source": "failed",
            **cstate, **props,
        }

    mol_h, pos, syms, c_idx_h, c_sym, c_z = emb
    relax_meta: Dict[str, Any] = {}
    if relax_with_sevenn:
        mol_h, pos, relax_meta = relax_3d_structure_with_sevenn(
            mol_h, pos, syms, fmax_ev_ang=fmax_ev_ang, max_steps=max_steps, sevenn_calc=sevenn_calc,
        )
        embed_source = f"{embed_source}+7net_lbfgs"

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
        "embed_source": embed_source,
        **relax_meta,
        **cstate,
        **props,
    }
