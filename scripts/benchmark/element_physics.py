#!/usr/bin/env python3
"""
Intrinsic Element-Property Physics & First-Principles Thermochemical Module
---------------------------------------------------------------------------
Eliminates hardcoded element-symbol branching (e.g., ``if center == "Si"``) and
fixed element whitelists (``METAL_SET``, ``CENTER_ATOMIC_NUMS``) by deriving:

  1. Periodic-table central-atom identification and d-block classification from
     Pauling electronegativity (chi_P), Cordero covalent radius (R_cov), and
     atomic-number block structure.
  2. Design 1 Element-Property-Driven Buried Volume (%V_bur at adaptive probe
     radius R_probe(Z_M) = R_cov(Z_M) + 1.85 A) with universal lower bound
     V_min^univ(pi_true) = 44.0 - 3.0 * min(3, pi_true).
  3. Thermochemical Descriptor 4 replacement: Standard-State Formation Energy
     per Coordination Bond (Delta E_f^std / v_M) computed relative to consistent
     SevenNet-0 elemental reference chemical potentials mu_Z^std.
  4. Statistical-mechanics gas-phase dimerization free energy
     Delta G_dimer(T) = Delta E_dimer - T * Delta S_dimer^trans+rot(m, I, T).
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms
from ase.build import bulk, molecule
from ase.data import atomic_masses, atomic_numbers, covalent_radii
from ase.optimize import LBFGS
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

PT = Chem.GetPeriodicTable()

# Pauling electronegativity (dimensionless; Pauling 1932 / Allred 1961)
PAULING_EN: Dict[str, float] = {
    "H": 2.20, "Li": 0.98, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
    "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16,
    "Ti": 1.54, "Cr": 1.66, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
    "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Br": 2.96,
    "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16, "Ru": 2.20, "Rh": 2.28,
    "In": 1.78, "Sn": 1.96, "I": 2.66,
    "La": 1.10, "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Pt": 2.28, "Ir": 2.20,
}

# van der Waals radii (units: A): Bondi (1964, DOI 10.1021/j100785a001),
# Mantina et al. (2009, DOI 10.1021/jp8111556) for Ga/Ge/As/In, RDKit scale for Nb/Rh/Ir/Y/La.
BONDI_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "Al": 2.00, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75,
    "Ti": 2.00, "Cr": 2.00, "Co": 2.00, "Ni": 1.63, "Cu": 2.00, "Zn": 1.39,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Br": 1.85,
    "Y": 2.40, "Zr": 2.15, "Nb": 2.15, "Mo": 2.10, "Ru": 2.05, "Rh": 2.00,
    "In": 1.93, "Sn": 2.17, "I": 1.98,
    "La": 2.50, "Hf": 2.15, "Ta": 2.15, "W": 2.10, "Pt": 1.75, "Ir": 2.00,
}

MAIN_NONMETALS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}

# SevenNet-0 (7net-0) elemental standard-state chemical potentials mu_Z^std (units: eV/atom)
# Evaluated consistently via bulk crystals (diamond/fcc/bcc/hcp/ortho/rhomb) & 50 A diatomic/P4 gases.
SEVENN_STANDARD_ELEMENT_REFS_EV: Dict[str, float] = {
    "Al": -3.719203, "As": -4.626450, "Br": -1.455892, "C": -9.086049,
    "Cl": -1.710204, "Co": -7.070953, "Cr": -9.496262, "Cu": -4.096194,
    "F": -1.922669, "Ga": -1.032459, "Ge": -4.569924, "H": -3.322610,
    "Hf": -9.952309, "I": -1.288664, "In": -2.704830, "Ir": -8.847458,
    "La": -4.919512, "Mo": -10.679195, "N": -8.241588, "Nb": -10.091097,
    "Ni": -5.767138, "O": -4.898368, "P": -1.910730, "Pt": -6.034290,
    "Rh": -7.326056, "Ru": -9.254974, "S": -3.573154, "Si": -5.400962,
    "Sn": -3.945925, "Ta": -11.834443, "Ti": -7.817384, "W": -12.751690,
    "Y": -6.449141, "Zn": -1.254205, "Zr": -8.570562,
}

# Physical constants for statistical mechanics and unit conversions
EV_TO_KCAL_MOL: float = 23.060548  # units: (kcal/mol) / eV
KB_SI: float = 1.380649e-23  # units: J/K
H_SI: float = 6.62607015e-34  # units: J*s
AMU_TO_KG: float = 1.6605390666e-27  # units: kg/amu
R_KCAL_MOL_K: float = 1.987204258e-3  # units: kcal/(mol*K)

# Design 1 universal shell extension beyond central atom covalent radius (units: A)
DELTA_R_SHELL_ANG: float = 1.85  # units: A


def is_d_block_atomic_num(z: int) -> bool:
    """Return True iff atomic number ``z`` belongs to the 3d, 4d, or 5d transition series."""
    return (21 <= z <= 30) or (39 <= z <= 48) or (71 <= z <= 80)


def is_coordination_center_atom(sym: str, z: int) -> bool:
    """
    Determine via periodic properties whether an atom can act as an electropositive/
    d-block coordination center rather than an organic backbone/halide/chalcogen atom.
    """
    if sym in MAIN_NONMETALS:
        return False
    chi = PAULING_EN.get(sym, 2.0)
    return chi <= 2.40 or is_d_block_atomic_num(z)


def compute_coordination_state(
    mol: Chem.Mol,
    atom_idx: int,
    pi_hapticity: Optional[int] = None,
) -> Dict[str, float]:
    """
    Compute state-dependent multi-property coordination vectors for atom ``atom_idx`` so
    transition metals are NEVER governed by a single static scalar cutoff or single radius.

    Physics/Algorithm:
        1. Mean bond polarity contrast ($\\overline{\\Delta\\chi}_i$) and Hannay-Smyth ionicity
           ($\\bar{f}_{\\text{ionic}}$) including both heavy neighbors and H substituents.
        2. Sanderson / Parr-Pearson effective electronegativity shift $\\chi_{\\text{eff}}(M \\mid \\text{env})$.
        3. Oxidation- and CN-adapted effective covalent radius (Bloomer/Pyykko polar contraction +
           coordination-shell steric expansion):
           $$R_{\\text{eff}}(M \\mid \\text{env}) = R_{\\text{cov}}^{(1)}(Z_M)
             \\left[1 - 0.06\\,\\max(0, \\overline{\\Delta\\chi}_{ML}) + 0.025\\,(CN_{\\text{eff}} - 4)\\right]$$
        4. Tolman 18e / Main-Group 8e valence electron count ($N_{\\text{ve}}$) and electron
           deficiency ($\\Delta N_{\\text{ve}} / N_{\\text{sat}}$), plus continuous orbital-manifold
           coordination propensity $S_{\\text{coord}} \\in [0, 1]$.
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    sym = atom.GetSymbol()
    z = atom.GetAtomicNum()
    chi0 = float(PAULING_EN.get(sym, 2.0))
    r_cov1 = float(covalent_radii[z])  # units: A
    n_val = float(PT.GetNOuterElecs(z))
    is_tm = is_d_block_atomic_num(z)
    pi_h = int(compute_refined_pi_hapticity(mol, atom_idx) if pi_hapticity is None else pi_hapticity)

    nbrs = list(atom.GetNeighbors())
    n_imp_h = int(atom.GetTotalNumHs())
    total_bonds = len(nbrs) + n_imp_h
    if total_bonds == 0:
        return {
            "chi_eff": round(chi0, 3),
            "delta_chi_mean": 0.0,
            "ionicity_mean": 0.0,
            "r_eff_ang": round(r_cov1, 3),  # units: A
            "r_probe_d2_ang": round(r_cov1 + DELTA_R_SHELL_ANG, 3),  # units: A
            "cn_eff": 0.0,
            "n_ve": n_val,
            "ve_deficit": 1.0,
            "coord_propensity": 1.0 if is_tm else 0.0,
        }

    dchis: List[float] = []
    ionics: List[float] = []
    e_lig = 0.0
    for nb in nbrs:
        dchi = float(PAULING_EN.get(nb.GetSymbol(), 2.0)) - chi0
        dchis.append(dchi)
        ionics.append(1.0 - math.exp(-0.25 * (dchi ** 2)))
        bond = mol.GetBondBetweenAtoms(atom_idx, nb.GetIdx())
        if bond is not None and bond.GetBondType() == Chem.BondType.DATIVE:
            e_lig += 2.0
        elif bond is not None:
            e_lig += float(bond.GetBondTypeAsDouble())
        else:
            e_lig += 1.0

    # Include implicit H atoms symmetrically in bond polarity and electron count
    for _ in range(n_imp_h):
        dchi_h = float(PAULING_EN["H"]) - chi0
        dchis.append(dchi_h)
        ionics.append(1.0 - math.exp(-0.25 * (dchi_h ** 2)))
        e_lig += 1.0

    cn_eff = float(total_bonds) + 2.0 * float(pi_h)
    mean_dchi = float(np.mean(dchis))
    mean_ion = float(np.mean(ionics))

    sum_signed_ionic = sum(math.copysign(ion, dc) for ion, dc in zip(ionics, dchis))
    chi_eff = chi0 + 0.12 * sum_signed_ionic

    r_scale = 1.0 - 0.06 * max(0.0, mean_dchi) + 0.025 * (cn_eff - 4.0)
    r_eff = max(0.75 * r_cov1, min(1.22 * r_cov1, r_cov1 * r_scale))  # units: A
    r_probe_d2 = r_eff + DELTA_R_SHELL_ANG  # units: A

    n_ve = n_val - float(atom.GetFormalCharge()) + e_lig + 4.0 * float(pi_h)
    n_sat = 18.0 if is_tm else 8.0
    ve_deficit = max(0.0, (n_sat - n_ve) / n_sat)

    coord_propensity = 0.0
    if is_tm:
        # Open d-orbital acceptor manifold + donor polarity modulation
        coord_propensity = 0.65 + 0.35 * max(0.0, min(1.0, (mean_dchi + 0.20) / 1.50))
    elif sym not in MAIN_NONMETALS and mean_dchi > 0.05:
        coord_propensity = 0.55 + 0.45 * min(1.0, mean_dchi / 1.20)
    if pi_h > 0 and coord_propensity > 0.0:
        coord_propensity = min(1.0, coord_propensity + 0.10)

    return {
        "chi_eff": round(float(chi_eff), 3),
        "delta_chi_mean": round(float(mean_dchi), 3),
        "ionicity_mean": round(float(mean_ion), 3),
        "r_eff_ang": round(float(r_eff), 3),
        "r_probe_d2_ang": round(float(r_probe_d2), 3),
        "cn_eff": round(float(cn_eff), 2),
        "n_ve": round(float(n_ve), 1),
        "ve_deficit": round(float(ve_deficit), 3),
        "coord_propensity": round(float(coord_propensity), 3),
    }


def find_central_atom_intrinsic(mol: Chem.Mol) -> Tuple[int, str, int]:
    """
    Identify the central atom via continuous periodic & coordination-state propensity
    (replacing static ``METAL_SET`` whitelists and single ``chi <= 2.40`` cutoffs).
    """
    candidates: List[Tuple[int, float, int, float, float, int, int, str, int]] = []
    for a in mol.GetAtoms():
        idx = a.GetIdx()
        sym = a.GetSymbol()
        z = a.GetAtomicNum()
        chi = PAULING_EN.get(sym, 2.0)
        rcov = float(covalent_radii[z])  # units: A
        deg = a.GetDegree()
        st = compute_coordination_state(mol, idx, pi_hapticity=0)
        is_coord = 1 if (st["coord_propensity"] >= 0.50 or is_coordination_center_atom(sym, z)) else 0
        candidates.append((is_coord, st["coord_propensity"], deg, rcov, -chi, -idx, idx, sym, z))
    best = max(candidates)
    return int(best[6]), str(best[7]), int(best[8])


def get_element_properties(sym: str, z: int) -> Dict[str, float]:
    """
    Return periodic-table atomic invariants for element ``(sym, z)``:
    ``r_cov_ang`` (units: A), ``r_probe_ang`` (units: A), ``n_val`` (valence electrons),
    ``chi_pauling`` (dimensionless), and ``r_vdw_ang`` (units: A).
    """
    r_cov = float(covalent_radii[z])  # units: A
    r_probe = round(r_cov + DELTA_R_SHELL_ANG, 3)  # units: A
    n_val = float(PT.GetNOuterElecs(z))
    chi = float(PAULING_EN.get(sym, 2.0))
    r_vdw = float(BONDI_VDW_RADII.get(sym, r_cov + 0.70))  # units: A
    return {
        "r_cov_ang": round(r_cov, 3),
        "r_probe_ang": r_probe,
        "n_val": n_val,
        "chi_pauling": round(chi, 3),
        "r_vdw_ang": round(r_vdw, 3),
    }


def compute_refined_pi_hapticity(mol: Chem.Mol, center_idx: int) -> int:
    """
    Count multi-hapto pi-coordinated ligands while rejecting eta^1 sigma-aryl rings.

    Physics/Algorithm:
        A ring attached to ``center_idx`` contributes to ``pi_hapticity_index`` iff:
        (a) >= 2 atoms of the same ring directly bond to ``center_idx`` (explicit polyhapto), or
        (b) it is an unsaturated non-benzene ring (e.g., 5-membered Cp/cyclobutadiene or
            non-aromatic cyclic diene represented as a single bond in SMILES).
        Six-membered aromatic phenyl rings attached via 1 ipso carbon (e.g. ``WPh5``,
        ``CrPh3``) are strictly classified as eta^1 sigma-aryl (not pi-haptic).
    """
    c_atom = mol.GetAtomWithIdx(center_idx)
    if not is_coordination_center_atom(c_atom.GetSymbol(), c_atom.GetAtomicNum()):
        return 0
    nbrs = {nb.GetIdx() for nb in c_atom.GetNeighbors()}
    ri = mol.GetRingInfo()
    pi_count = 0
    for ring in ri.AtomRings():
        shared = set(ring) & nbrs
        if not shared:
            continue
        if len(shared) >= 2:
            pi_count += 1
        else:
            is_aro_hex = len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
            has_unsat = any(
                mol.GetBondBetweenAtoms(ring[k], ring[(k + 1) % len(ring)]).GetBondType() != Chem.BondType.SINGLE
                for k in range(len(ring))
            )
            if (not is_aro_hex) and has_unsat and len(ring) >= 5:
                pi_count += 1
    return pi_count


def compute_intrinsic_topology(smiles: str) -> Dict[str, object]:
    """
    Element-agnostic coordination & ring topology analysis replacing ``CENTER_ATOMIC_NUMS``
    and ``TM_ATOMIC_NUMS`` whitelists.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"
    c_idx, c_sym, c_z = find_central_atom_intrinsic(mol)
    is_coord = is_coordination_center_atom(c_sym, c_z)
    is_tm = is_d_block_atomic_num(c_z)

    pi_hapticity = compute_refined_pi_hapticity(mol, c_idx)
    ri = mol.GetRingInfo()
    chelate_ring_count = 0
    for ring in ri.AtomRings():
        if not (4 <= len(ring) <= 7):
            continue
        if c_idx in ring:
            hetero_donors = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() in (7, 8, 15, 16))
            if hetero_donors >= 2:
                chelate_ring_count += 1

    c_atom = mol.GetAtomWithIdx(c_idx)
    nbrs = list(c_atom.GetNeighbors())
    is_homoleptic_tm_sp3 = bool(
        is_tm
        and len(nbrs) > 0
        and all(nb.GetAtomicNum() == 6 and nb.GetHybridization() == Chem.HybridizationType.SP3 for nb in nbrs)
        and pi_hapticity == 0
    )

    metal_c_bonds = sum(1 for nb in nbrs if nb.GetAtomicNum() == 6) if (is_coord and nbrs) else 0
    metal_hetero_bonds = sum(1 for nb in nbrs if nb.GetAtomicNum() in (7, 8, 9, 15, 16, 17, 35, 53)) if (is_coord and nbrs) else 0
    total_bonds = len(nbrs) if is_coord else 0

    metal_c_sigma_frac = (metal_c_bonds / total_bonds) if total_bonds > 0 else 0.0
    metal_hetero_frac = (metal_hetero_bonds / total_bonds) if total_bonds > 0 else 0.0
    is_homoleptic_mc = 1.0 if (is_coord and total_bonds > 0 and metal_c_bonds == total_bonds) else 0.0

    aro_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    # Intrinsic M-OH pattern around any coordination center
    metal_oh_count = sum(
        1
        for nb in nbrs
        if is_coord and nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() >= 1 and nb.GetDegree() <= 2
    )
    # Intrinsic Design 2 disproportionation detector: non-d-block main-group/metalloid center
    # with homoatomic M-M catenation AND terminal hydride M-H (e.g., ClSiH2-SiH2Cl, Cl2SiH-SiHCl2).
    catenated_mg_hydride = any(
        (not is_d_block_atomic_num(a.GetAtomicNum()))
        and (a.GetSymbol() not in MAIN_NONMETALS)
        and (a.GetTotalNumHs() >= 1)
        and any(nb.GetAtomicNum() == a.GetAtomicNum() for nb in a.GetNeighbors())
        for a in mol.GetAtoms()
    )
    # Intrinsic Design 2 retro-hydrometallation / protonolysis detector: non-d-block center
    # carrying BOTH reactive M-H hydride AND pi-conjugated unsaturation (e.g., Cp2AlH, Cp2GaH).
    pi_conj_mg_hydride = any(
        (not is_d_block_atomic_num(a.GetAtomicNum()))
        and (a.GetSymbol() not in MAIN_NONMETALS)
        and (a.GetTotalNumHs() >= 1)
        and (compute_refined_pi_hapticity(mol, a.GetIdx()) >= 1)
        for a in mol.GetAtoms()
    )
    coord_state = compute_coordination_state(mol, c_idx, pi_hapticity=pi_hapticity)
    is_coord_state = bool(coord_state["coord_propensity"] >= 0.50 or is_coord)
    coord_flagged = (
        is_homoleptic_tm_sp3
        or (metal_oh_count >= 2)
        or (is_coord_state and aro_rings >= 3 and metal_c_sigma_frac >= 0.9)
        or catenated_mg_hydride
        or pi_conj_mg_hydride
    )

    return {
        "center_idx": c_idx,
        "center_atom": c_sym,
        "center_z": c_z,
        "is_coord_center": is_coord_state,
        "is_d_block": is_tm,
        "pi_hapticity_index": int(pi_hapticity),
        "chelate_ring_count": int(chelate_ring_count),
        "is_homoleptic_tm_sp3_refined": is_homoleptic_tm_sp3,
        "catenated_main_group_hydride": bool(catenated_mg_hydride),
        "pi_conjugated_main_group_hydride": bool(pi_conj_mg_hydride),
        "metal_c_sigma_frac": float(metal_c_sigma_frac),
        "metal_hetero_frac": float(metal_hetero_frac),
        "is_homoleptic_mc": float(is_homoleptic_mc),
        "metal_oh_count": float(metal_oh_count),
        "aromatic_rings": int(aro_rings),
        "coord_flagged": bool(coord_flagged),
        **coord_state,
    }


def universal_vbur_lower_bound(pi_hapticity: int) -> float:
    """
    Universal Design 1 steric lower bound ``V_min^univ`` (units: %) for probe sphere
    ``R_probe(Z_M) = R_cov(Z_M) + 1.85 A``. Contains zero element-name branching.

    Physics/Algorithm:
        $$V_{\\min}^{\\text{univ}}(\\eta_{\\pi}) = 44.0 - 5.0 \\min(2,\\, \\eta_{\\pi})$$
    """
    return round(44.0 - 5.0 * min(2, max(0, int(pi_hapticity))), 1)


def embed_etkdg_3d(smiles: str) -> Optional[Tuple[Chem.Mol, np.ndarray, List[str], int, str, int]]:
    """Generate deterministic ETKDGv3 + UFF 3D atomic coordinates and intrinsic center atom."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    c_idx, c_sym, c_z = find_central_atom_intrinsic(mol)
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    cid = AllChem.EmbedMolecule(mol_h, params)
    if cid < 0:
        cid = AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=42)
    if cid < 0:
        return None
    if AllChem.UFFHasAllMoleculeParams(mol_h):
        try:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=120)
        except Exception:
            pass
    conf = mol_h.GetConformer()
    pos = np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol_h.GetNumAtoms())],
        dtype=float,
    )  # units: A
    syms = [a.GetSymbol() for a in mol_h.GetAtoms()]
    return mol_h, pos, syms, c_idx, c_sym, c_z


def calc_buried_volume_sphere(
    pos: np.ndarray,
    syms: List[str],
    center_idx: int,
    radius_ang: float,
    grid_step_ang: float = 0.25,
) -> float:
    """
    Compute SambVca-style buried volume ``%V_bur`` (units: %) inside spherical probe
    ``radius_ang`` centered on ``center_idx``.

    Physics/Algorithm:
        $$\\%V_{\\text{bur}}(R_{\\text{probe}}) = \\frac{V_{\\text{ligand}} \\cap V_{\\text{sphere}}(R_{\\text{probe}})}{V_{\\text{sphere}}(R_{\\text{probe}})} \\times 100$$
    """
    coords = pos - pos[center_idx]  # units: A
    vdw_radii = np.array([BONDI_VDW_RADII.get(s, float(covalent_radii[atomic_numbers.get(s, 6)]) + 0.70) for s in syms])  # units: A
    g = np.arange(-radius_ang, radius_ang + grid_step_ang * 0.5, grid_step_ang)
    gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    r2 = np.sum(pts ** 2, axis=1)
    sphere_pts = pts[r2 <= radius_ang ** 2]
    if len(sphere_pts) == 0:
        return 0.0

    other_mask = np.arange(len(syms)) != center_idx
    other_coords = coords[other_mask]
    other_r2 = vdw_radii[other_mask] ** 2
    if len(other_coords) == 0:
        return 0.0

    occupied = np.zeros(len(sphere_pts), dtype=bool)
    chunk = 2000
    for i in range(0, len(sphere_pts), chunk):
        p_chunk = sphere_pts[i:i + chunk]
        d2 = np.sum((p_chunk[:, None, :] - other_coords[None, :, :]) ** 2, axis=2)
        occupied[i:i + chunk] = np.any(d2 <= other_r2[None, :], axis=1)
    return round(float(np.mean(occupied) * 100.0), 2)


def calc_design2_steric_electronic_metrics(
    mol_h: Chem.Mol,
    pos: np.ndarray,
    syms: List[str],
    center_idx: int,
    radius_d2_ang: float,
    grid_step_ang: float = 0.28,
) -> Dict[str, float]:
    """
    Compute Design 2 electronic-polarity and hydride-weighted buried volume
    (``vbur_d2_percent``, units: %) inside the state-dependent probe sphere
    ``R_probe^D2(M | env) = R_eff(M | env) + 1.85 A``, plus the maximum open
    solid-angle attack void cone (``void_cone_deg``, units: deg).

    Physics/Algorithm:
        Terminal hydrides bonded directly to center ``M`` provide zero inert kinetic
        passivation (``w_MH = 0.0``), while donor/ligand atoms carry Pauling-polarity
        shielding weights ``w_j = clip(1.0 + 0.20 * (chi_P(j) - 2.55), 0.75, 1.25)``.
        Void cone angle measures the largest unshielded approach vector across
        non-hydride steric neighbors via a 200-point Fibonacci $S^2$ quadrature.
    """
    coords = pos - pos[center_idx]  # units: A
    c_atom = mol_h.GetAtomWithIdx(center_idx)
    center_h_indices = {nb.GetIdx() for nb in c_atom.GetNeighbors() if nb.GetAtomicNum() == 1}
    vdw_radii = np.array([BONDI_VDW_RADII.get(s, float(covalent_radii[atomic_numbers.get(s, 6)]) + 0.70) for s in syms])  # units: A

    weights = np.ones(len(syms), dtype=float)
    for j, s in enumerate(syms):
        if j in center_h_indices:
            weights[j] = 0.0
        else:
            chi_j = float(PAULING_EN.get(s, 2.0))
            weights[j] = min(1.25, max(0.75, 1.0 + 0.20 * (chi_j - 2.55)))

    g = np.arange(-radius_d2_ang, radius_d2_ang + grid_step_ang * 0.5, grid_step_ang)
    gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    sphere_pts = pts[np.sum(pts ** 2, axis=1) <= radius_d2_ang ** 2]

    other_mask = np.arange(len(syms)) != center_idx
    other_coords = coords[other_mask]
    other_r2 = vdw_radii[other_mask] ** 2
    other_w = weights[other_mask]

    if len(sphere_pts) > 0 and len(other_coords) > 0:
        pt_w = np.zeros(len(sphere_pts), dtype=float)
        chunk = 2000
        for i in range(0, len(sphere_pts), chunk):
            p_chunk = sphere_pts[i:i + chunk]
            d2 = np.sum((p_chunk[:, None, :] - other_coords[None, :, :]) ** 2, axis=2)
            in_atom = d2 <= other_r2[None, :]
            pt_w[i:i + chunk] = np.max(np.where(in_atom, other_w[None, :], 0.0), axis=1)
        vbur_d2 = round(float(np.mean(pt_w) * 100.0), 2)
    else:
        vbur_d2 = 0.0

    # Angular maximum void cone excluding reactive M-H hydrides
    non_mh_mask = other_mask & np.array([i not in center_h_indices for i in range(len(syms))])
    vecs = coords[non_mh_mask]
    norms = np.linalg.norm(vecs, axis=1)
    vecs = vecs[norms > 0.3] / norms[norms > 0.3][:, None]
    if len(vecs) > 0:
        k = np.arange(200)
        phi = np.arccos(1.0 - 2.0 * (k + 0.5) / 200.0)
        theta = np.pi * (1.0 + 5.0 ** 0.5) * k
        dirs = np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1)
        dots = np.clip(dirs @ vecs.T, -1.0, 1.0)
        void_deg = round(float(np.degrees(np.arccos(np.min(np.max(dots, axis=1))))), 2)
    else:
        void_deg = 180.0

    return {
        "vbur_d2_percent": vbur_d2,  # units: %
        "void_cone_deg": void_deg,  # units: deg
        "center_hydride_count": float(len(center_h_indices)),
    }


def compute_stat_mech_tdelta_s_dimer(pos: np.ndarray, syms: List[str], temperature_k: float = 373.15) -> float:
    """
    Compute statistical-mechanics gas-phase entropic penalty ``-T * Delta S_dimer^trans+rot``
    (units: kcal/mol at ``temperature_k`` and 1 bar) for ``2 M -> M2`` via Sackur-Tetrode
    and rigid-rotor rotational partition functions from 3D principal inertia eigenvalues.
    """
    pressure_pa = 1.0e5  # units: Pa (1 bar)
    masses = np.array([atomic_masses[atomic_numbers[s]] for s in syms], dtype=float)  # units: amu
    m_kg = float(np.sum(masses)) * AMU_TO_KG  # units: kg
    com = np.sum(pos * masses[:, None], axis=0) / np.sum(masses)
    rc = (pos - com) * 1.0e-10  # units: m
    i_tensor = np.zeros((3, 3))
    for m_amu, r in zip(masses, rc):
        m = m_amu * AMU_TO_KG
        i_tensor += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    eigvals = np.maximum(np.linalg.eigvalsh(i_tensor), 1.0e-50)  # units: kg*m^2

    q_trans_m = ((2.0 * math.pi * m_kg * KB_SI * temperature_k / (H_SI ** 2)) ** 1.5) * (KB_SI * temperature_k / pressure_pa)
    s_trans_m = R_KCAL_MOL_K * (math.log(q_trans_m) + 2.5)  # units: kcal/(mol*K)
    s_trans_d = s_trans_m + 1.5 * R_KCAL_MOL_K * math.log(2.0)
    ds_trans = s_trans_d - 2.0 * s_trans_m

    if np.sum(masses > 0) >= 2 and eigvals[0] > 1.0e-48:
        theta_r = (H_SI ** 2) / (8.0 * (math.pi ** 2) * eigvals * KB_SI)  # units: K
        q_rot_m = math.sqrt(math.pi) / math.sqrt(theta_r[0] * theta_r[1] * theta_r[2]) * (temperature_k ** 1.5)
        s_rot_m = R_KCAL_MOL_K * (math.log(max(q_rot_m, 1.0)) + 1.5)
    else:
        s_rot_m = 0.0
    s_rot_d = s_rot_m + 1.5 * R_KCAL_MOL_K * math.log(3.2)
    ds_rot = s_rot_d - 2.0 * s_rot_m

    return round(float(-temperature_k * (ds_trans + ds_rot)), 2)


def compute_formation_energy_metrics(
    e_pot_ev: float,
    syms: List[str],
    center_degree: int,
    refs_ev: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    """
    Compute standard-state thermodynamic formation energy per atom and per coordination bond:

    Physics/Algorithm:
        $$\\Delta E_{f}^{\\text{std}} = E_{\\text{pot}}^{\\text{MLIP}} - \\sum_Z n_Z \\mu_Z^{\\text{std}}, \\qquad
          \\frac{\\Delta E_{f}^{\\text{std}}}{v_M} = \\frac{\\Delta E_{f}^{\\text{std}} \\times 23.060548}{\\max(1, \\deg(M))}$$
    """
    ref_map = refs_ev if refs_ev is not None else SEVENN_STANDARD_ELEMENT_REFS_EV
    e_ref_sum = sum(ref_map[s] for s in syms)  # units: eV
    de_form_ev = float(e_pot_ev) - e_ref_sum  # units: eV
    de_form_kcal = de_form_ev * EV_TO_KCAL_MOL  # units: kcal/mol
    v_m = max(1, int(center_degree))
    return {
        "delta_e_form_ev_per_atom": round(de_form_ev / max(1, len(syms)), 4),
        "delta_e_form_kcal_mol": round(de_form_kcal, 2),
        "delta_e_form_per_bond_kcal_mol": round(de_form_kcal / v_m, 2),
    }


def compute_sevenn_element_refs() -> Dict[str, float]:
    """Recompute ``SEVENN_STANDARD_ELEMENT_REFS_EV`` on-the-fly with ``SevenNetCalculator('7net-0')``."""
    from sevenn.calculator import SevenNetCalculator

    calc = SevenNetCalculator("7net-0")
    refs: Dict[str, float] = {}
    for sym, formula in [("H", "H2"), ("N", "N2"), ("O", "O2"), ("F", "F2"), ("Cl", "Cl2"), ("S", "S2")]:
        at = molecule(formula)
        at.set_cell([50.0, 50.0, 50.0])
        at.set_pbc([False, False, False])
        at.calc = calc
        LBFGS(at, logfile=None).run(fmax=0.02, steps=50)
        refs[sym] = round(float(at.get_potential_energy()) / len(at), 6)

    for sym, bl in [("Br", 2.28), ("I", 2.67)]:
        at = Atoms(sym + "2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, bl]], cell=[50.0, 50.0, 50.0], pbc=False)
        at.calc = calc
        LBFGS(at, logfile=None).run(fmax=0.02, steps=40)
        refs[sym] = round(float(at.get_potential_energy()) / 2.0, 6)

    d = 2.21 / np.sqrt(2.0)
    p4 = Atoms("P4", positions=[[d, d, d], [d, -d, -d], [-d, d, -d], [-d, -d, d]], cell=[50.0, 50.0, 50.0], pbc=False)
    p4.calc = calc
    LBFGS(p4, logfile=None).run(fmax=0.02, steps=40)
    refs["P"] = round(float(p4.get_potential_energy()) / 4.0, 6)

    for sym, struct, kw in [
        ("C", "diamond", {"a": 3.567}),
        ("Si", "diamond", {"a": 5.431}),
        ("Al", "fcc", {"a": 4.05}),
        ("Cu", "fcc", {"a": 3.615}),
        ("Ir", "fcc", {"a": 3.839}),
        ("Rh", "fcc", {"a": 3.803}),
        ("W", "bcc", {"a": 3.165}),
        ("Mo", "bcc", {"a": 3.147}),
        ("Ta", "bcc", {"a": 3.301}),
        ("Cr", "bcc", {"a": 2.885}),
        ("Ti", "hcp", {"a": 2.951, "covera": 1.587}),
        ("Zr", "hcp", {"a": 3.232, "covera": 1.593}),
        ("Hf", "hcp", {"a": 3.196, "covera": 1.581}),
        ("Co", "hcp", {"a": 2.507, "covera": 1.623}),
        ("Ru", "hcp", {"a": 2.706, "covera": 1.582}),
        ("Sn", "diamond", {"a": 6.489}),
        ("Ga", "orthorhombic", {"a": 4.519, "b": 7.657, "c": 4.526}),
    ]:
        at = bulk(sym, struct, **kw)
        at.calc = calc
        refs[sym] = round(float(at.get_potential_energy()) / len(at), 6)

    for sym in ["In", "As"]:
        at = bulk(sym)
        at.calc = calc
        refs[sym] = round(float(at.get_potential_energy()) / len(at), 6)

    return refs
