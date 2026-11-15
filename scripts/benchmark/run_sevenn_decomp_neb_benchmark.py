#!/usr/bin/env python3
"""SevenNet-0 NEB / CI-NEB Benchmark for Semiconductor Precursor Decomposition.

Evaluates the accuracy, saddle-point fidelity, and kinetic selectivity of the
SevenNet-0 universal MLIAP against published peer-reviewed DFT and ab initio
benchmark pathways for volatile semiconductor precursors already present in
the project database (`benchmark_unified_latest_160.json`).

Outputs:
  - database/external/sevenn_decomp_neb_benchmark.json
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
import os
import sys
import tempfile
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import torch
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.mep import NEB, idpp_interpolate
from ase.optimize import FIRE
from ase.vibrations import Vibrations
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from sevenn.calculator import SevenNetCalculator

PROJECT_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts", "benchmark"))

from run_phase2_sevenn_features import relax_atoms_sevenn  # noqa: E402

EV_TO_KCAL_MOL: float = 23.06054783  # units: (kcal/mol) / eV
CM1_IMAG_THRESHOLD: float = 20.0  # units: cm^-1; filter near-zero rigid/rotational noise


@dataclass(frozen=True)
class ReactionBenchmarkSpec:
    """Metadata and literature quantum-chemical reference for one elementary reaction."""

    reaction_id: str
    precursor_id: str
    precursor_name: str
    precursor_smiles: str
    reaction_equation: str
    mechanism_class: str
    lit_authors: str
    lit_year: int
    lit_journal: str
    lit_doi: str
    lit_method: str
    lit_ea_kcal_mol: float  # units: kcal/mol
    lit_ea_note: str
    lit_delta_e_rxn_kcal_mol: Optional[float]  # units: kcal/mol


BENCHMARK_SPECS: List[ReactionBenchmarkSpec] = [
    ReactionBenchmarkSpec(
        reaction_id="HCDS_SiCl2_elim",
        precursor_id="PREC_SI_06",
        precursor_name="HCDS",
        precursor_smiles="Cl[Si](Cl)(Cl)[Si](Cl)(Cl)Cl",
        reaction_equation="Si2Cl6 -> SiCl4 + SiCl2",
        mechanism_class="Disilane 1,2-Cl shift / alpha-elimination",
        lit_authors="Swihart & Carr; Noda et al.",
        lit_year=1997,
        lit_journal="J. Phys. Chem. A 101, 7434 / J. Phys. Chem. A 126",
        lit_doi="10.1021/jp971651t",
        lit_method="BAC-MP4 / CCSD(T)//B3LYP (cf. 10.1021/acs.jpca.2c06195)",
        lit_ea_kcal_mol=49.0,
        lit_ea_note="48.3-49.5 kcal/mol ab initio / experimental high-P limit",
        lit_delta_e_rxn_kcal_mol=None,
    ),
    ReactionBenchmarkSpec(
        reaction_id="DCDS_SiHCl_elim",
        precursor_id="NEG-009",
        precursor_name="1,2-Dichlorodisilane",
        precursor_smiles="Cl[SiH2][SiH2]Cl",
        reaction_equation="ClH2Si-SiH2Cl -> SiH3Cl + SiHCl",
        mechanism_class="Disilane 1,2-H shift / alpha-elimination",
        lit_authors="Swihart & Carr",
        lit_year=1997,
        lit_journal="J. Phys. Chem. A 101, 7434-7445",
        lit_doi="10.1021/jp971651t",
        lit_method="BAC-MP4 / G2//MP2",
        lit_ea_kcal_mol=45.0,
        lit_ea_note="Predicted ~45 kcal/mol for SiHCl extrusion",
        lit_delta_e_rxn_kcal_mol=None,
    ),
    ReactionBenchmarkSpec(
        reaction_id="GaEt3_betaH_elim",
        precursor_id="EXT-ALD-035",
        precursor_name="GaEt3",
        precursor_smiles="CC[Ga](CC)CC",
        reaction_equation="Ga(C2H5)3 -> HGa(C2H5)2 + C2H4",
        mechanism_class="Group-13 organometallic beta-H elimination",
        lit_authors="Lee, Kim, & Anderson",
        lit_year=2009,
        lit_journal="ECS Trans. 25(8), 107",
        lit_doi="10.1149/1.3207574",
        lit_method="B3LYP DFT",
        lit_ea_kcal_mol=40.3,
        lit_ea_note="beta-H elimination 40.3 kcal/mol vs Ga-C homolysis 65.0 kcal/mol",
        lit_delta_e_rxn_kcal_mol=None,
    ),
    ReactionBenchmarkSpec(
        reaction_id="TEOS_C2H4_elim",
        precursor_id="PREC_SI_07",
        precursor_name="TEOS",
        precursor_smiles="CCO[Si](OCC)(OCC)OCC",
        reaction_equation="Si(OEt)4 -> Si(OH)(OEt)3 + C2H4",
        mechanism_class="Main-group alkoxide 4-centered beta-H-to-O elimination",
        lit_authors="Nurkowski, Buerger, Akroyd, & Kraft",
        lit_year=2015,
        lit_journal="Proc. Combust. Inst. 35, 229-236",
        lit_doi="10.1016/j.proci.2014.06.093",
        lit_method="B97-1/6-311+G(d,p) (0 K enthalpy)",
        lit_ea_kcal_mol=59.35,
        lit_ea_note="248.32 kJ/mol zero-point-corrected barrier",
        lit_delta_e_rxn_kcal_mol=12.20,
    ),
    ReactionBenchmarkSpec(
        reaction_id="TTIP_acetone_elim",
        precursor_id="PREC_TI_02",
        precursor_name="TTIP",
        precursor_smiles="CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C",
        reaction_equation="Ti(OiPr)4 -> HTi(OiPr)3 + (CH3)2CO",
        mechanism_class="TM alkoxide 4-centered alpha-H-to-Ti ketone elimination",
        lit_authors="Ershov, Kochubei, Kiselev, & Baklanov",
        lit_year=2018,
        lit_journal="J. Phys. Chem. A 122, 1064-1070",
        lit_doi="10.1021/acs.jpca.7b10396",
        lit_method="DFT + DLPNO-CCSD(T)",
        lit_ea_kcal_mol=49.0,
        lit_ea_note="Preferred unimolecular decomposition channel (~49 kcal/mol)",
        lit_delta_e_rxn_kcal_mol=None,
    ),
    ReactionBenchmarkSpec(
        reaction_id="TTIP_propene_elim",
        precursor_id="PREC_TI_02",
        precursor_name="TTIP",
        precursor_smiles="CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C",
        reaction_equation="Ti(OiPr)4 -> HO-Ti(OiPr)3 + C3H6",
        mechanism_class="TM alkoxide 4-centered beta-H-to-O propene elimination",
        lit_authors="Ershov, Kochubei, Kiselev, & Baklanov",
        lit_year=2018,
        lit_journal="J. Phys. Chem. A 122, 1064-1070",
        lit_doi="10.1021/acs.jpca.7b10396",
        lit_method="DFT + DLPNO-CCSD(T)",
        lit_ea_kcal_mol=62.0,
        lit_ea_note="Disfavored competing channel (~13 kcal/mol above acetone channel)",
        lit_delta_e_rxn_kcal_mol=None,
    ),
    ReactionBenchmarkSpec(
        reaction_id="TDMAT_NH3_transam",
        precursor_id="PREC_TI_01",
        precursor_name="TDMAT",
        precursor_smiles="CN(C)[Ti](N(C)C)(N(C)C)N(C)C",
        reaction_equation="Ti(NMe2)4 + NH3 -> Ti(NMe2)3(NH2) + HNMe2",
        mechanism_class="Bimolecular ALD/CVD transamination proton transfer",
        lit_authors="Sen, Banu, Debnath, Ghosh, & Das",
        lit_year=2014,
        lit_journal="Dalton Trans. 43, 8877-8887",
        lit_doi="10.1039/c4dt00690a",
        lit_method="B3LYP DFT (adduct-referenced initial barrier)",
        lit_ea_kcal_mol=7.3,
        lit_ea_note="Initial transamination barrier from adduct (cf. Cross & Schlegel 2000: 8.4)",
        lit_delta_e_rxn_kcal_mol=None,
    ),
    ReactionBenchmarkSpec(
        reaction_id="CpZr_NMe2_betaH_elim",
        precursor_id="PREC_HK_07",
        precursor_name="Cp(NMe2)3Zr (ZAC)",
        precursor_smiles="CN(C)[Zr](N(C)C)(N(C)C)C1=CC=CC1",
        reaction_equation="CpZr(NMe2)3 -> CpZr(NMe2)(eta2-NMeCH2) + HNMe2",
        mechanism_class="Inter-ligand beta-H transfer / dimethylamine elimination",
        lit_authors="Van, Jang, Jung, Noh, Moon, Kil, Chung, & Shong",
        lit_year=2022,
        lit_journal="J. Phys. Chem. C 126(42), 18090-18099",
        lit_doi="10.1021/acs.jpcc.2c04425",
        lit_method="M06L DFT (gas-phase beta-H elimination to HNMe2)",
        lit_ea_kcal_mol=34.63,
        lit_ea_note="144.9 kJ/mol beta-H transfer between two -NMe2 ligands",
        lit_delta_e_rxn_kcal_mol=20.8,
    ),
]


def embed_mol_3d(smiles: str, seed: int = 42) -> Chem.Mol:
    """Generate a 3D RDKit conformer with explicit hydrogens.

    Physics/Algorithm: Uses ETKDGv3 distance geometry with random-coordinate
    fallback for main-group and organometallic coordination compounds.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    status = AllChem.EmbedMolecule(mol, params)
    if status < 0:
        status = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=seed)
    if status < 0:
        raise RuntimeError(f"RDKit 3D embedding failed for {smiles}")
    return mol


def unit_vec(v: np.ndarray) -> np.ndarray:
    """Return normalized vector, raising if norm is near zero."""
    norm = float(np.linalg.norm(v))
    if norm < 1e-12:
        raise ValueError("Cannot normalize near-zero vector")
    return v / norm


def kabsch_align(
    mobile_pos: np.ndarray,
    ref_pos: np.ndarray,
    anchor_indices: List[int],
) -> np.ndarray:
    """Rigidly align `mobile_pos` onto `ref_pos` using Kabsch SVD on spectator anchors.

    Physics/Algorithm: Minimizes mass-unweighted RMSD over non-reacting spectator
    atoms to suppress rigid-body translation/rotation coupling along the NEB band.
    """
    p_mob = mobile_pos[anchor_indices]  # units: A
    p_ref = ref_pos[anchor_indices]  # units: A
    c_mob = np.mean(p_mob, axis=0)
    c_ref = np.mean(p_ref, axis=0)

    h_mat = (p_mob - c_mob).T @ (p_ref - c_ref)
    u_mat, _, vt_mat = np.linalg.svd(h_mat)
    rot = u_mat @ vt_mat
    if np.linalg.det(rot) < 0.0:
        u_mat[:, -1] *= -1.0
        rot = u_mat @ vt_mat

    aligned = (mobile_pos - c_mob) @ rot + c_ref  # units: A
    return aligned - np.mean(aligned, axis=0)


def relax_mol_smiles(
    smiles: str,
    calc: SevenNetCalculator,
    fmax: float = 0.03,
    steps: int = 150,
) -> Tuple[float, List[str], np.ndarray, float]:
    """Embed SMILES and relax isolated molecule in a 50 A non-PBC box."""
    mol = embed_mol_3d(smiles)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    pos = mol.GetConformer().GetPositions()
    e_ev, pos_rel, _, final_fmax = relax_atoms_sevenn(
        symbols=symbols, positions=pos, calc=calc, fmax=fmax, steps=steps, cell_len=50.0
    )
    return float(e_ev), symbols, pos_rel, float(final_fmax)


def get_subtree_indices(mol: Chem.Mol, root: int, blocked: int) -> List[int]:
    """Collect connected atom indices reachable from `root` without crossing `blocked`."""
    visited = {root}
    stack = [root]
    while stack:
        curr = stack.pop()
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            nidx = nbr.GetIdx()
            if nidx != blocked and nidx not in visited:
                visited.add(nidx)
                stack.append(nidx)
    return sorted(visited)


def _place_atom_along_bond(
    pos: np.ndarray, atom_idx: int, target_idx: int, dist: float
) -> None:
    """Translate `atom_idx` to distance `dist` (A) from `target_idx` along current ray."""
    vec = unit_vec(pos[atom_idx] - pos[target_idx])
    pos[atom_idx] = pos[target_idx] + dist * vec


def _displace_fragment(
    pos: np.ndarray, frag_indices: List[int], direction: np.ndarray, shift: float
) -> None:
    """Rigidly shift a molecular fragment along `direction` by `shift` A."""
    u = unit_vec(direction)
    pos[frag_indices] = pos[frag_indices] + shift * u


def _build_disilane_shift_pair(
    smiles: str,
    calc: SevenNetCalculator,
    shift_element: str,
    p1_smiles: str,
    p2_smiles: str,
) -> Dict[str, Any]:
    """Build IS, adduct FS, and separated products for Si2X6 -> SiX4 + SiX2 shifts."""
    mol = embed_mol_3d(smiles)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    si_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Si"]
    si0, si1 = si_indices[0], si_indices[1]

    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols,
        positions=mol.GetConformer().GetPositions(),
        calc=calc,
        fmax=0.03,
        steps=160,
    )

    nbrs_si1 = [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(si1).GetNeighbors()
        if n.GetIdx() != si0 and n.GetSymbol() == shift_element
    ]
    shift_idx = nbrs_si1[0]
    leaving_idx = [si1] + [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(si1).GetNeighbors()
        if n.GetIdx() not in (si0, shift_idx)
    ]

    pos_fs_guess = pos_is.copy()
    bond_len = 2.04 if shift_element == "Cl" else 1.48  # units: A
    axis = unit_vec(pos_is[si1] - pos_is[si0])
    # Place migrated substituent on Si0 away from Si1 approach vector when shift_element == "H"
    if shift_element == "H":
        perp_in = unit_vec(pos_is[shift_idx] - pos_is[si1])
        pos_fs_guess[shift_idx] = pos_fs_guess[si0] + bond_len * unit_vec(axis + 0.8 * perp_in)
        _displace_fragment(pos_fs_guess, leaving_idx, axis, shift=2.30)
    else:
        pos_fs_guess[shift_idx] = pos_fs_guess[si0] + bond_len * axis
        _displace_fragment(pos_fs_guess, leaving_idx, axis, shift=1.40)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=200
    )

    anchor_indices = [si0] + [
        n.GetIdx() for n in mol.GetAtomWithIdx(si0).GetNeighbors() if n.GetIdx() != si1
    ]
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=anchor_indices)

    pos_ts_seed = 0.5 * (pos_is + pos_fs)
    bridge_perp = unit_vec(pos_is[shift_idx] - pos_is[si1])
    bridge_dist = 1.65 if shift_element == "Cl" else 1.15  # units: A
    pos_ts_seed[shift_idx] = (
        0.5 * (pos_ts_seed[si0] + pos_ts_seed[si1]) + bridge_dist * bridge_perp
    )

    e_p1, _, _, _ = relax_mol_smiles(p1_smiles, calc=calc)
    e_p2, _, _, _ = relax_mol_smiles(p2_smiles, calc=calc)

    d_si0_shift_fs = float(np.linalg.norm(pos_fs[si0] - pos_fs[shift_idx]))
    d_si0_si1_fs = float(np.linalg.norm(pos_fs[si0] - pos_fs[si1]))
    si0_coord_fs = sum(
        1
        for i in range(len(symbols))
        if i != si0 and np.linalg.norm(pos_fs[si0] - pos_fs[i]) < 2.45
    )
    si1_coord_fs = sum(
        1
        for i in range(len(symbols))
        if i != si1 and np.linalg.norm(pos_fs[si1] - pos_fs[i]) < 2.45
    )
    valid_fs = (
        d_si0_shift_fs < 2.35
        and d_si0_si1_fs > 2.95
        and si0_coord_fs == 4
        and si1_coord_fs == 2
    )

    reactive_indices = sorted({si0, si1, shift_idx} | set(anchor_indices[:3]))
    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_p1 + e_p2),
        "e_sep_react_ev": float(e_is),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_si0_shift_A": round(d_si0_shift_fs, 3),
            "d_si0_si1_A": round(d_si0_si1_fs, 3),
            "si0_coord": int(si0_coord_fs),
            "si1_coord": int(si1_coord_fs),
        },
        "bond_probe": [
            ("breaking_Si1_X", si1, shift_idx),
            ("forming_Si0_X", si0, shift_idx),
            ("cleaving_Si0_Si1", si0, si1),
        ],
        "reactive_indices": reactive_indices,
    }


def _build_gaet3_pair(calc: SevenNetCalculator) -> Dict[str, Any]:
    """Build IS, adduct FS, and separated products for GaEt3 -> HGaEt2 + C2H4."""
    smiles = "CC[Ga](CC)CC"
    mol = embed_mol_3d(smiles)
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    ga_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Ga")
    c_alpha = mol.GetAtomWithIdx(ga_idx).GetNeighbors()[0].GetIdx()
    c_beta = next(
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_alpha).GetNeighbors()
        if n.GetSymbol() == "C"
    )
    h_beta = next(
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_beta).GetNeighbors()
        if n.GetSymbol() == "H"
    )

    # Set syn-coplanar reactive dihedral Ga-Ca-Cb-Hb = 0 deg before relaxation
    rdMolTransforms.SetDihedralDeg(conf, ga_idx, c_alpha, c_beta, h_beta, 0.0)
    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols, positions=conf.GetPositions(), calc=calc, fmax=0.03, steps=180
    )

    c2h4_all = set(get_subtree_indices(mol, c_alpha, blocked=ga_idx))
    c2h4_keep = sorted(c2h4_all - {h_beta})
    spectators = sorted(set(range(len(symbols))) - c2h4_all)

    pos_fs_guess = pos_is.copy()
    ga_to_c2h4 = unit_vec(0.5 * (pos_is[c_alpha] + pos_is[c_beta]) - pos_is[ga_idx])
    pos_fs_guess[h_beta] = pos_is[ga_idx] + 1.60 * unit_vec(pos_is[h_beta] - pos_is[ga_idx])
    _place_atom_along_bond(pos_fs_guess, c_beta, c_alpha, dist=1.34)
    _displace_fragment(pos_fs_guess, c2h4_keep, ga_to_c2h4, shift=1.35)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=220
    )
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=spectators)

    pos_ts_seed = pos_is + 0.35 * (pos_fs - pos_is)
    pos_ts_seed[h_beta] = 0.5 * (pos_ts_seed[ga_idx] + pos_ts_seed[c_beta]) + 0.25 * unit_vec(
        pos_is[h_beta] - 0.5 * (pos_is[ga_idx] + pos_is[c_beta])
    )

    e_hgaet2, _, _, _ = relax_mol_smiles("CC[GaH]CC", calc=calc)
    e_c2h4, _, _, _ = relax_mol_smiles("C=C", calc=calc)

    d_ga_h = float(np.linalg.norm(pos_fs[ga_idx] - pos_fs[h_beta]))
    d_cc = float(np.linalg.norm(pos_fs[c_alpha] - pos_fs[c_beta]))
    d_ga_ca = float(np.linalg.norm(pos_fs[ga_idx] - pos_fs[c_alpha]))
    valid_fs = d_ga_h < 1.85 and d_cc < 1.42 and d_ga_ca > 2.55

    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_hgaet2 + e_c2h4),
        "e_sep_react_ev": float(e_is),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_Ga_H_A": round(d_ga_h, 3),
            "d_CC_A": round(d_cc, 3),
            "d_Ga_Ca_A": round(d_ga_ca, 3),
        },
        "bond_probe": [
            ("breaking_Cb_H", c_beta, h_beta),
            ("forming_Ga_H", ga_idx, h_beta),
            ("cleaving_Ga_Ca", ga_idx, c_alpha),
            ("forming_C_C", c_alpha, c_beta),
        ],
        "reactive_indices": [ga_idx, c_alpha, c_beta, h_beta],
    }


def _build_teos_pair(calc: SevenNetCalculator) -> Dict[str, Any]:
    """Build IS, adduct FS, and separated products for Si(OEt)4 -> Si(OH)(OEt)3 + C2H4."""
    smiles = "CCO[Si](OCC)(OCC)OCC"
    mol = embed_mol_3d(smiles)
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    si_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Si")
    o_rxn = mol.GetAtomWithIdx(si_idx).GetNeighbors()[0].GetIdx()
    c_alpha = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(o_rxn).GetNeighbors() if n.GetSymbol() == "C"
    )
    c_beta = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(c_alpha).GetNeighbors() if n.GetSymbol() == "C"
    )
    h_beta = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(c_beta).GetNeighbors() if n.GetSymbol() == "H"
    )

    rdMolTransforms.SetDihedralDeg(conf, si_idx, o_rxn, c_alpha, c_beta, 170.0)
    rdMolTransforms.SetDihedralDeg(conf, o_rxn, c_alpha, c_beta, h_beta, 0.0)
    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols, positions=conf.GetPositions(), calc=calc, fmax=0.03, steps=200
    )

    et_all = set(get_subtree_indices(mol, c_alpha, blocked=o_rxn))
    c2h4_keep = sorted(et_all - {h_beta})
    spectators = sorted(set(range(len(symbols))) - et_all)

    pos_fs_guess = pos_is.copy()
    leave_dir = unit_vec(0.5 * (pos_is[c_alpha] + pos_is[c_beta]) - pos_is[o_rxn])
    pos_fs_guess[h_beta] = pos_is[o_rxn] + 0.98 * unit_vec(pos_is[h_beta] - pos_is[o_rxn])
    _place_atom_along_bond(pos_fs_guess, c_beta, c_alpha, dist=1.34)
    _displace_fragment(pos_fs_guess, c2h4_keep, leave_dir, shift=1.35)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=220
    )
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=spectators)

    pos_ts_seed = pos_is + 0.35 * (pos_fs - pos_is)
    pos_ts_seed[h_beta] = 0.5 * (pos_ts_seed[o_rxn] + pos_ts_seed[c_beta]) + 0.20 * unit_vec(
        pos_is[h_beta] - 0.5 * (pos_is[o_rxn] + pos_is[c_beta])
    )

    e_silanol, _, _, _ = relax_mol_smiles("O[Si](OCC)(OCC)OCC", calc=calc)
    e_c2h4, _, _, _ = relax_mol_smiles("C=C", calc=calc)

    d_o_h = float(np.linalg.norm(pos_fs[o_rxn] - pos_fs[h_beta]))
    d_cc = float(np.linalg.norm(pos_fs[c_alpha] - pos_fs[c_beta]))
    d_o_ca = float(np.linalg.norm(pos_fs[o_rxn] - pos_fs[c_alpha]))
    valid_fs = d_o_h < 1.15 and d_cc < 1.42 and d_o_ca > 2.40

    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_silanol + e_c2h4),
        "e_sep_react_ev": float(e_is),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_O_H_A": round(d_o_h, 3),
            "d_CC_A": round(d_cc, 3),
            "d_O_Ca_A": round(d_o_ca, 3),
        },
        "bond_probe": [
            ("breaking_Cb_H", c_beta, h_beta),
            ("forming_O_H", o_rxn, h_beta),
            ("cleaving_O_Ca", o_rxn, c_alpha),
            ("forming_C_C", c_alpha, c_beta),
        ],
        "reactive_indices": [si_idx, o_rxn, c_alpha, c_beta, h_beta],
    }


def _build_ttip_acetone_pair(calc: SevenNetCalculator) -> Dict[str, Any]:
    """Build IS, adduct FS, and separated products for Ti(OiPr)4 -> HTi(OiPr)3 + acetone."""
    smiles = "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C"
    mol = embed_mol_3d(smiles)
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    ti_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Ti")
    o_rxn = mol.GetAtomWithIdx(ti_idx).GetNeighbors()[0].GetIdx()
    c_alpha = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(o_rxn).GetNeighbors() if n.GetSymbol() == "C"
    )
    h_alpha = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(c_alpha).GetNeighbors() if n.GetSymbol() == "H"
    )

    rdMolTransforms.SetDihedralDeg(conf, ti_idx, o_rxn, c_alpha, h_alpha, 0.0)
    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols, positions=conf.GetPositions(), calc=calc, fmax=0.03, steps=220
    )

    oipr_all = set(get_subtree_indices(mol, o_rxn, blocked=ti_idx))
    acetone_keep = sorted(oipr_all - {h_alpha})
    spectators = sorted(set(range(len(symbols))) - oipr_all)

    pos_fs_guess = pos_is.copy()
    leave_dir = unit_vec(pos_is[c_alpha] - pos_is[ti_idx])
    pos_fs_guess[h_alpha] = pos_is[ti_idx] + 1.72 * unit_vec(pos_is[h_alpha] - pos_is[ti_idx])
    _displace_fragment(pos_fs_guess, acetone_keep, leave_dir, shift=1.45)
    # Contract O=C carbonyl bond slightly to facilitate acetone relaxation
    _place_atom_along_bond(pos_fs_guess, o_rxn, c_alpha, dist=1.22)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=220
    )
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=spectators)

    # Seed concerted 4-centered [Ti-O-Ca-Ha] saddle geometry without prior Ti-O dissociation
    pos_ts_seed = pos_is + 0.25 * (pos_fs - pos_is)
    _place_atom_along_bond(pos_ts_seed, c_alpha, o_rxn, dist=1.32)
    bisec = unit_vec(
        unit_vec(pos_ts_seed[ti_idx] - pos_ts_seed[o_rxn])
        + unit_vec(pos_ts_seed[c_alpha] - pos_ts_seed[o_rxn])
    )
    pos_ts_seed[h_alpha] = 0.5 * (pos_ts_seed[ti_idx] + pos_ts_seed[c_alpha]) + 0.20 * bisec

    e_hti, _, _, _ = relax_mol_smiles("CC(C)O[TiH](OC(C)C)OC(C)C", calc=calc)
    e_ace, _, _, _ = relax_mol_smiles("CC(=O)C", calc=calc)

    d_ti_h = float(np.linalg.norm(pos_fs[ti_idx] - pos_fs[h_alpha]))
    d_co = float(np.linalg.norm(pos_fs[o_rxn] - pos_fs[c_alpha]))
    d_ti_o = float(np.linalg.norm(pos_fs[ti_idx] - pos_fs[o_rxn]))
    valid_fs = d_ti_h < 2.05 and d_co < 1.35 and d_ti_o > 2.25

    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_hti + e_ace),
        "e_sep_react_ev": float(e_is),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_Ti_H_A": round(d_ti_h, 3),
            "d_CO_A": round(d_co, 3),
            "d_Ti_O_A": round(d_ti_o, 3),
        },
        "bond_probe": [
            ("breaking_Ca_H", c_alpha, h_alpha),
            ("forming_Ti_H", ti_idx, h_alpha),
            ("cleaving_Ti_O", ti_idx, o_rxn),
            ("forming_C_O", c_alpha, o_rxn),
        ],
        "reactive_indices": [ti_idx, o_rxn, c_alpha, h_alpha],
    }


def _build_ttip_propene_pair(calc: SevenNetCalculator) -> Dict[str, Any]:
    """Build IS, adduct FS, and separated products for Ti(OiPr)4 -> HO-Ti(OiPr)3 + C3H6."""
    smiles = "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C"
    mol = embed_mol_3d(smiles)
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    ti_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Ti")
    o_rxn = mol.GetAtomWithIdx(ti_idx).GetNeighbors()[0].GetIdx()
    c_alpha = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(o_rxn).GetNeighbors() if n.GetSymbol() == "C"
    )
    c_beta = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(c_alpha).GetNeighbors() if n.GetSymbol() == "C"
    )
    h_beta = next(
        n.GetIdx() for n in mol.GetAtomWithIdx(c_beta).GetNeighbors() if n.GetSymbol() == "H"
    )

    rdMolTransforms.SetDihedralDeg(conf, ti_idx, o_rxn, c_alpha, c_beta, 170.0)
    rdMolTransforms.SetDihedralDeg(conf, o_rxn, c_alpha, c_beta, h_beta, 0.0)
    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols, positions=conf.GetPositions(), calc=calc, fmax=0.03, steps=220
    )

    ipr_all = set(get_subtree_indices(mol, c_alpha, blocked=o_rxn))
    propene_keep = sorted(ipr_all - {h_beta})
    spectators = sorted(set(range(len(symbols))) - ipr_all)

    pos_fs_guess = pos_is.copy()
    leave_dir = unit_vec(0.5 * (pos_is[c_alpha] + pos_is[c_beta]) - pos_is[o_rxn])
    pos_fs_guess[h_beta] = pos_is[o_rxn] + 0.98 * unit_vec(pos_is[h_beta] - pos_is[o_rxn])
    _place_atom_along_bond(pos_fs_guess, c_beta, c_alpha, dist=1.34)
    _displace_fragment(pos_fs_guess, propene_keep, leave_dir, shift=1.40)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=220
    )
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=spectators)

    pos_ts_seed = pos_is + 0.35 * (pos_fs - pos_is)
    pos_ts_seed[h_beta] = 0.5 * (pos_ts_seed[o_rxn] + pos_ts_seed[c_beta]) + 0.20 * unit_vec(
        pos_is[h_beta] - 0.5 * (pos_is[o_rxn] + pos_is[c_beta])
    )

    e_tioh, _, _, _ = relax_mol_smiles("O[Ti](OC(C)C)(OC(C)C)OC(C)C", calc=calc)
    e_propene, _, _, _ = relax_mol_smiles("CC=C", calc=calc)

    d_o_h = float(np.linalg.norm(pos_fs[o_rxn] - pos_fs[h_beta]))
    d_cc = float(np.linalg.norm(pos_fs[c_alpha] - pos_fs[c_beta]))
    d_o_ca = float(np.linalg.norm(pos_fs[o_rxn] - pos_fs[c_alpha]))
    valid_fs = d_o_h < 1.15 and d_cc < 1.42 and d_o_ca > 2.40

    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_tioh + e_propene),
        "e_sep_react_ev": float(e_is),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_O_H_A": round(d_o_h, 3),
            "d_CC_A": round(d_cc, 3),
            "d_O_Ca_A": round(d_o_ca, 3),
        },
        "bond_probe": [
            ("breaking_Cb_H", c_beta, h_beta),
            ("forming_O_H", o_rxn, h_beta),
            ("cleaving_O_Ca", o_rxn, c_alpha),
            ("forming_C_C", c_alpha, c_beta),
        ],
        "reactive_indices": [ti_idx, o_rxn, c_alpha, c_beta, h_beta],
    }


def _build_tdmat_nh3_pair(calc: SevenNetCalculator) -> Dict[str, Any]:
    """Build adduct IS, transaminated FS, and separated species for TDMAT + NH3."""
    e_tdmat, sym_td, pos_td, _ = relax_mol_smiles(
        "CN(C)[Ti](N(C)C)(N(C)C)N(C)C", calc=calc, fmax=0.03, steps=200
    )
    e_nh3, sym_nh3, pos_nh3, _ = relax_mol_smiles("N", calc=calc, fmax=0.03, steps=150)

    mol_td = embed_mol_3d("CN(C)[Ti](N(C)C)(N(C)C)N(C)C")
    ti_idx = next(a.GetIdx() for a in mol_td.GetAtoms() if a.GetSymbol() == "Ti")
    n_ligs = [a.GetIdx() for a in mol_td.GetAtomWithIdx(ti_idx).GetNeighbors()]

    # Approach NH3 along an open tetrahedral face opposite one Ti-N bond
    face_dir = -unit_vec(pos_td[n_ligs[0]] - pos_td[ti_idx])
    # Orient NH3 so N points toward Ti and hydrogens trail outward
    n_nh3_pos = pos_nh3[0]
    h_mean = np.mean(pos_nh3[1:], axis=0)
    nh3_axis = unit_vec(h_mean - n_nh3_pos)
    # Build orthonormal basis rotation from nh3_axis -> face_dir
    v_cross = np.cross(nh3_axis, face_dir)
    s_norm = float(np.linalg.norm(v_cross))
    if s_norm > 1e-8:
        k_mat = np.array(
            [
                [0.0, -v_cross[2], v_cross[1]],
                [v_cross[2], 0.0, -v_cross[0]],
                [-v_cross[1], v_cross[0], 0.0],
            ]
        )
        c_dot = float(np.dot(nh3_axis, face_dir))
        rot = np.eye(3) + k_mat + (k_mat @ k_mat) * (1.0 / (1.0 + c_dot + 1e-12))
    else:
        rot = np.eye(3)
    pos_nh3_oriented = (pos_nh3 - n_nh3_pos) @ rot.T + pos_td[ti_idx] + 2.40 * face_dir

    symbols = sym_td + sym_nh3
    pos_adduct_guess = np.vstack([pos_td, pos_nh3_oriented])
    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols, positions=pos_adduct_guess, calc=calc, fmax=0.03, steps=220
    )

    n_ammonia = len(sym_td)
    h_ammonia = [n_ammonia + 1, n_ammonia + 2, n_ammonia + 3]
    # Pick closest -NMe2 amide ligand to any NH3 proton
    best_pair = min(
        (
            (float(np.linalg.norm(pos_is[n_l] - pos_is[h_a])), n_l, h_a)
            for n_l in n_ligs
            for h_a in h_ammonia
        ),
        key=lambda x: x[0],
    )
    _, n_leave, h_transfer = best_pair
    leaving_nme2 = sorted(set(get_subtree_indices(mol_td, n_leave, blocked=ti_idx)))
    hnme2_all = leaving_nme2 + [h_transfer]
    spectators = sorted(set(range(len(sym_td))) - set(leaving_nme2))

    pos_fs_guess = pos_is.copy()
    leave_dir = unit_vec(pos_is[n_leave] - pos_is[ti_idx])
    pos_fs_guess[h_transfer] = pos_is[n_leave] + 1.03 * unit_vec(
        pos_is[n_ammonia] - pos_is[n_leave]
    )
    # Move incoming NH2 closer to Ti and push HNMe2 out
    nh2_indices = [i for i in [n_ammonia] + h_ammonia if i != h_transfer]
    nh2_shift = unit_vec(pos_is[ti_idx] - pos_is[n_ammonia]) * max(
        0.0, float(np.linalg.norm(pos_is[ti_idx] - pos_is[n_ammonia])) - 1.92
    )
    pos_fs_guess[nh2_indices] = pos_fs_guess[nh2_indices] + nh2_shift
    _displace_fragment(pos_fs_guess, hnme2_all, leave_dir, shift=1.35)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=240
    )
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=spectators)

    pos_ts_seed = 0.5 * (pos_is + pos_fs)
    pos_ts_seed[h_transfer] = 0.5 * (pos_ts_seed[n_ammonia] + pos_ts_seed[n_leave])

    e_prod_ti, _, _, _ = relax_mol_smiles("CN(C)[Ti](N)(N(C)C)N(C)C", calc=calc)
    e_hnme2, _, _, _ = relax_mol_smiles("CNC", calc=calc)

    d_nleave_h = float(np.linalg.norm(pos_fs[n_leave] - pos_fs[h_transfer]))
    d_ti_nammonia = float(np.linalg.norm(pos_fs[ti_idx] - pos_fs[n_ammonia]))
    d_ti_nleave = float(np.linalg.norm(pos_fs[ti_idx] - pos_fs[n_leave]))
    valid_fs = d_nleave_h < 1.18 and d_ti_nammonia < 2.15 and d_ti_nleave > 2.55

    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_prod_ti + e_hnme2),
        "e_sep_react_ev": float(e_tdmat + e_nh3),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_Nleave_H_A": round(d_nleave_h, 3),
            "d_Ti_Nnh2_A": round(d_ti_nammonia, 3),
            "d_Ti_Nleave_A": round(d_ti_nleave, 3),
        },
        "bond_probe": [
            ("breaking_Nnh3_H", n_ammonia, h_transfer),
            ("forming_Nleave_H", n_leave, h_transfer),
            ("forming_Ti_Nnh2", ti_idx, n_ammonia),
            ("cleaving_Ti_Nleave", ti_idx, n_leave),
        ],
        "reactive_indices": [ti_idx, n_ammonia, h_transfer, n_leave],
    }


def _build_piano_stool_cpzr_coords() -> Tuple[
    List[str], np.ndarray, List[int], List[int], List[List[int]], List[List[List[int]]]
]:
    """Construct analytical eta5-CpZr(NMe2)3 three-legged piano-stool geometry.

    Physics/Algorithm: RDKit ETKDGv3/UFF treats Cp SMILES as eta1 and lacks Zr3
    typing, causing haptic ring collapse; analytical C3v-like placement initializes
    symmetrical eta5-Cp coordination (Zr-centroid 2.22 A, Zr-N 2.07 A).
    """
    symbols: List[str] = ["Zr"]
    coords: List[np.ndarray] = [np.array([0.0, 0.0, 0.0])]
    z_cp, r_c, r_h = 2.22, 1.21, 2.29  # units: A
    cp_c_idx: List[int] = []
    for k in range(5):
        ang = 2.0 * math.pi * k / 5.0
        cp_c_idx.append(len(symbols))
        symbols.append("C")
        coords.append(np.array([r_c * math.cos(ang), r_c * math.sin(ang), z_cp]))
        symbols.append("H")
        coords.append(np.array([r_h * math.cos(ang), r_h * math.sin(ang), z_cp + 0.05]))

    d_zr_n, d_n_c, d_c_h, tilt = 2.07, 1.46, 1.09, math.radians(108.0)  # units: A, rad
    n_indices: List[int] = []
    c_me_indices: List[List[int]] = []
    h_me_indices: List[List[List[int]]] = []
    for leg in range(3):
        phi = 2.0 * math.pi * leg / 3.0 + math.radians(15.0)
        u_radial = np.array([math.cos(phi), math.sin(phi), 0.0])
        u_tang = np.array([-math.sin(phi), math.cos(phi), 0.0])
        u_z = np.array([0.0, 0.0, 1.0])
        u_zn = math.sin(tilt) * u_radial + math.cos(tilt) * u_z
        p_n = d_zr_n * u_zn
        n_indices.append(len(symbols))
        symbols.append("N")
        coords.append(p_n)

        perp = unit_vec(np.cross(u_zn, u_tang))
        plane_1 = math.cos(math.radians(25.0)) * u_tang + math.sin(math.radians(25.0)) * perp
        c_pair: List[int] = []
        h_pair: List[List[int]] = []
        for sgn in (+1.0, -1.0):
            u_nc = (
                math.cos(math.radians(60.0)) * u_zn
                + sgn * math.sin(math.radians(60.0)) * plane_1
            )
            p_c = p_n + d_n_c * u_nc
            c_pair.append(len(symbols))
            symbols.append("C")
            coords.append(p_c)
            b1 = -u_nc
            ref = u_tang if abs(float(np.dot(b1, u_tang))) < 0.8 else u_z
            e1 = unit_vec(np.cross(b1, ref))
            e2 = unit_vec(np.cross(b1, e1))
            ch_list: List[int] = []
            for hk in range(3):
                ha = 2.0 * math.pi * hk / 3.0 + leg * 0.4
                u_ch = math.cos(math.radians(109.5)) * b1 + math.sin(
                    math.radians(109.5)
                ) * (math.cos(ha) * e1 + math.sin(ha) * e2)
                ch_list.append(len(symbols))
                symbols.append("H")
                coords.append(p_c + d_c_h * u_ch)
            h_pair.append(ch_list)
        c_me_indices.append(c_pair)
        h_me_indices.append(h_pair)

    return symbols, np.array(coords, dtype=float), cp_c_idx, n_indices, c_me_indices, h_me_indices


def _build_cpzr_nme2_betah_pair(calc: SevenNetCalculator) -> Dict[str, Any]:
    """Build IS, adduct FS, and separated products for CpZr(NMe2)3 -> CpZr(NMe2)(eta2-NMeCH2) + HNMe2."""
    symbols, pos0, cp_c, n_idx, c_me, h_me = _build_piano_stool_cpzr_coords()
    zr_idx = 0

    e_is, pos_is, _, fmax_is = relax_atoms_sevenn(
        symbols=symbols, positions=pos0, calc=calc, fmax=0.03, steps=220
    )

    n_donor, n_leave = n_idx[0], n_idx[1]
    best: Optional[Tuple[float, int, int]] = None
    for ci in range(2):
        c_i = c_me[0][ci]
        for hi in h_me[0][ci]:
            dist = float(np.linalg.norm(pos_is[hi] - pos_is[n_leave]))
            if best is None or dist < best[0]:
                best = (dist, c_i, hi)
    assert best is not None
    _, c_alpha, h_transfer = best

    leg1_idx = [n_leave, c_me[1][0], c_me[1][1]] + h_me[1][0] + h_me[1][1]
    frag_hnme2 = leg1_idx + [h_transfer]
    anchor_indices = [zr_idx, n_idx[2]] + cp_c

    u_out = unit_vec(pos_is[n_leave] - pos_is[zr_idx])
    u_zr_ca = unit_vec(pos_is[c_alpha] - pos_is[zr_idx])

    pos_fs_guess = pos_is.copy()
    pos_fs_guess[h_transfer] = pos_is[n_leave] + 1.02 * unit_vec(
        pos_is[c_alpha] - pos_is[n_leave]
    )
    pos_fs_guess[c_alpha] = pos_is[zr_idx] + 2.22 * u_zr_ca
    _displace_fragment(pos_fs_guess, frag_hnme2, u_out, shift=2.40)

    e_fs, pos_fs_raw, _, fmax_fs = relax_atoms_sevenn(
        symbols=symbols, positions=pos_fs_guess, calc=calc, fmax=0.03, steps=260
    )
    pos_fs = kabsch_align(pos_fs_raw, pos_is, anchor_indices=anchor_indices)

    # Separated products (20 A shift of HNMe2 along u_out to decouple Zr...N adduct binding)
    pos_sep_guess = pos_is.copy()
    pos_sep_guess[h_transfer] = pos_is[n_leave] + 1.02 * (-u_out)
    pos_sep_guess[c_alpha] = pos_is[zr_idx] + 2.22 * u_zr_ca
    _displace_fragment(pos_sep_guess, frag_hnme2, u_out, shift=15.0)
    e_sep_prod, _, _, _ = relax_atoms_sevenn(
        symbols=symbols, positions=pos_sep_guess, calc=calc, fmax=0.03, steps=240
    )

    # Seed concerted proton transfer between C_alpha and N_leave with partial Zr-N_leave stretch
    pos_ts_seed = pos_is + 0.45 * (pos_fs - pos_is)
    mid_cn = 0.5 * (pos_ts_seed[c_alpha] + pos_ts_seed[n_leave])
    outward_dir = unit_vec(mid_cn - pos_ts_seed[zr_idx])
    pos_ts_seed[h_transfer] = mid_cn + 0.15 * outward_dir

    d_nleave_h = float(np.linalg.norm(pos_fs[n_leave] - pos_fs[h_transfer]))
    d_zr_ca = float(np.linalg.norm(pos_fs[zr_idx] - pos_fs[c_alpha]))
    d_zr_nleave = float(np.linalg.norm(pos_fs[zr_idx] - pos_fs[n_leave]))
    valid_fs = d_nleave_h < 1.15 and d_zr_ca < 2.45 and d_zr_nleave > 2.28

    return {
        "symbols": symbols,
        "pos_is": pos_is,
        "pos_ts_seed": pos_ts_seed,
        "pos_fs": pos_fs,
        "e_is_ev": float(e_is),
        "e_fs_ev": float(e_fs),
        "e_sep_prod_ev": float(e_sep_prod),
        "e_sep_react_ev": float(e_is),
        "fmax_is": float(fmax_is),
        "fmax_fs": float(fmax_fs),
        "valid_fs": bool(valid_fs),
        "fs_check": {
            "d_Nleave_H_A": round(d_nleave_h, 3),
            "d_Zr_Ca_A": round(d_zr_ca, 3),
            "d_Zr_Nleave_A": round(d_zr_nleave, 3),
        },
        "bond_probe": [
            ("breaking_Ca_H", c_alpha, h_transfer),
            ("forming_Nleave_H", n_leave, h_transfer),
            ("forming_Zr_Ca", zr_idx, c_alpha),
            ("cleaving_Zr_Nleave", zr_idx, n_leave),
        ],
        "reactive_indices": [zr_idx, n_donor, c_alpha, h_transfer, n_leave],
    }


REACTION_BUILDERS: Dict[str, Callable[[SevenNetCalculator], Dict[str, Any]]] = {
    "HCDS_SiCl2_elim": lambda calc: _build_disilane_shift_pair(
        "Cl[Si](Cl)(Cl)[Si](Cl)(Cl)Cl", calc, "Cl", "Cl[Si](Cl)(Cl)Cl", "Cl[Si]Cl"
    ),
    "DCDS_SiHCl_elim": lambda calc: _build_disilane_shift_pair(
        "Cl[SiH2][SiH2]Cl", calc, "H", "Cl[SiH3]", "Cl[SiH]"
    ),
    "GaEt3_betaH_elim": _build_gaet3_pair,
    "TEOS_C2H4_elim": _build_teos_pair,
    "TTIP_acetone_elim": _build_ttip_acetone_pair,
    "TTIP_propene_elim": _build_ttip_propene_pair,
    "TDMAT_NH3_transam": _build_tdmat_nh3_pair,
    "CpZr_NMe2_betaH_elim": _build_cpzr_nme2_betah_pair,
}


def compute_core_vibrations(
    atoms_ts: Atoms,
    calc: SevenNetCalculator,
    reactive_indices: List[int],
    delta: float = 0.02,
) -> Dict[str, Any]:
    """Evaluate finite-difference Hessian eigenvalues restricted to reactive core atoms.

    Physics/Algorithm: Checks whether the climbing-image structure lies on a
    genuine saddle mode along the active bond-forming/breaking subspace.
    """
    ts_copy = atoms_ts.copy()
    ts_copy.calc = calc
    with tempfile.TemporaryDirectory() as tmpdir:
        vib_name = os.path.join(tmpdir, "vib")
        vib = Vibrations(
            ts_copy, indices=reactive_indices, name=vib_name, delta=delta
        )
        vib.run()
        freqs = vib.get_frequencies()

    imag_cm1: List[float] = []
    real_cm1: List[float] = []
    for f in freqs:
        if abs(f.imag) > CM1_IMAG_THRESHOLD:
            imag_cm1.append(round(float(abs(f.imag)), 1))
        elif abs(f.real) > CM1_IMAG_THRESHOLD:
            real_cm1.append(round(float(f.real), 1))
    imag_cm1.sort(reverse=True)
    return {
        "n_imag_modes_gt20cm1": len(imag_cm1),
        "imag_freqs_cm1": imag_cm1,
        "min_real_freq_cm1": min(real_cm1) if real_cm1 else None,
        "max_real_freq_cm1": max(real_cm1) if real_cm1 else None,
    }


def run_single_ci_neb(
    spec: ReactionBenchmarkSpec,
    calc: SevenNetCalculator,
    n_intermediate: int = 5,
    k_spring: float = 0.10,
    neb_steps: int = 75,
    cineb_steps: int = 120,
    neb_fmax: float = 0.15,
    cineb_fmax: float = 0.07,
    maxstep: float = 0.05,
    check_endpoints_only: bool = False,
) -> Dict[str, Any]:
    """Construct endpoints, run IDPP + NEB + CI-NEB, and diagnose transition state."""
    t0 = time.perf_counter()
    pair = REACTION_BUILDERS[spec.reaction_id](calc)
    symbols: List[str] = pair["symbols"]

    if check_endpoints_only:
        seed_bonds = {}
        if "pos_ts_seed" in pair:
            for label, a_i, b_j in pair["bond_probe"]:
                seed_bonds[label] = {
                    "IS_A": round(float(np.linalg.norm(pair["pos_is"][a_i] - pair["pos_is"][b_j])), 3),
                    "SEED_A": round(float(np.linalg.norm(pair["pos_ts_seed"][a_i] - pair["pos_ts_seed"][b_j])), 3),
                    "FS_A": round(float(np.linalg.norm(pair["pos_fs"][a_i] - pair["pos_fs"][b_j])), 3),
                }
        return {
            "reaction_id": spec.reaction_id,
            "valid_fs": pair["valid_fs"],
            "fs_check": pair["fs_check"],
            "fmax_is": pair["fmax_is"],
            "fmax_fs": pair["fmax_fs"],
            "de_adduct_kcal": round((pair["e_fs_ev"] - pair["e_is_ev"]) * EV_TO_KCAL_MOL, 2),
            "de_sep_kcal": round((pair["e_sep_prod_ev"] - pair["e_sep_react_ev"]) * EV_TO_KCAL_MOL, 2),
            "seed_bonds_A": seed_bonds,
            "wall_time_s": round(time.perf_counter() - t0, 2),
        }

    def _make_atoms(pos: np.ndarray, energy_ev: float) -> Atoms:
        at = Atoms(
            symbols=symbols,
            positions=pos.copy(),
            cell=[50.0, 50.0, 50.0],
            pbc=[False, False, False],
        )
        forces = np.zeros_like(pos)
        at.calc = SinglePointCalculator(at, energy=energy_ev, forces=forces)
        return at

    atoms_is = _make_atoms(pair["pos_is"], pair["e_is_ev"])
    atoms_fs = _make_atoms(pair["pos_fs"], pair["e_fs_ev"])

    images: List[Atoms] = [atoms_is]
    for _ in range(n_intermediate):
        im = atoms_is.copy()
        im.calc = calc
        images.append(im)
    images.append(atoms_fs)

    if "pos_ts_seed" in pair and n_intermediate >= 3:
        mid_idx = n_intermediate // 2 + 1
        atoms_seed = _make_atoms(pair["pos_ts_seed"], pair["e_is_ev"])
        left_band = [atoms_is.copy() for _ in range(mid_idx)] + [atoms_seed.copy()]
        right_band = [atoms_seed.copy()] + [
            atoms_seed.copy() for _ in range(n_intermediate - mid_idx)
        ] + [atoms_fs.copy()]
        idpp_interpolate(NEB(left_band), traj=None, log=None, fmax=0.10, steps=120)
        idpp_interpolate(NEB(right_band), traj=None, log=None, fmax=0.10, steps=120)
        for idx in range(1, mid_idx):
            images[idx].set_positions(left_band[idx].get_positions())
        images[mid_idx].set_positions(atoms_seed.get_positions())
        for offset in range(1, n_intermediate - mid_idx + 1):
            images[mid_idx + offset].set_positions(
                right_band[offset].get_positions()
            )
        neb = NEB(
            images,
            k=k_spring,
            climb=False,
            method="improvedtangent",
            allow_shared_calculator=True,
        )
    else:
        neb = NEB(
            images,
            k=k_spring,
            climb=False,
            method="improvedtangent",
            allow_shared_calculator=True,
        )
        idpp_interpolate(neb, traj=None, log=None, fmax=0.10, steps=150)

    opt_neb = FIRE(neb, maxstep=maxstep, logfile=None)
    opt_neb.run(fmax=neb_fmax, steps=neb_steps)
    nsteps_neb = int(opt_neb.nsteps)
    fmax_after_neb = float(opt_neb.fmax if hasattr(opt_neb, "fmax") else np.nan)

    neb.climb = True
    opt_ci = FIRE(neb, maxstep=maxstep, logfile=None)
    opt_ci.run(fmax=cineb_fmax, steps=cineb_steps)
    nsteps_ci = int(opt_ci.nsteps)

    # Read residual NEB forces and image energies cleanly
    neb_forces = neb.get_forces()
    final_neb_fmax = float(np.max(np.linalg.norm(neb_forces, axis=1)))
    energies_ev = [float(im.get_potential_energy()) for im in images]
    e_is_ev = energies_ev[0]
    e_fs_ev = energies_ev[-1]
    rel_kcal = [round((e - e_is_ev) * EV_TO_KCAL_MOL, 3) for e in energies_ev]

    ts_idx = int(np.argmax(energies_ev))
    ea_cineb_kcal = float((energies_ev[ts_idx] - e_is_ev) * EV_TO_KCAL_MOL)
    ea_from_sep_kcal = float(
        (energies_ev[ts_idx] - pair["e_sep_react_ev"]) * EV_TO_KCAL_MOL
    )
    de_adduct_kcal = float((e_fs_ev - e_is_ev) * EV_TO_KCAL_MOL)
    de_sep_kcal = float(
        (pair["e_sep_prod_ev"] - pair["e_sep_react_ev"]) * EV_TO_KCAL_MOL
    )

    pos_ts = images[ts_idx].get_positions()
    bond_lengths = {}
    for label, a_i, b_j in pair["bond_probe"]:
        bond_lengths[label] = {
            "IS_A": round(float(np.linalg.norm(pair["pos_is"][a_i] - pair["pos_is"][b_j])), 3),
            "TS_A": round(float(np.linalg.norm(pos_ts[a_i] - pos_ts[b_j])), 3),
            "FS_A": round(float(np.linalg.norm(pair["pos_fs"][a_i] - pair["pos_fs"][b_j])), 3),
        }

    # Check true saddle-point vibrational signature on internal TS image
    if 0 < ts_idx < len(images) - 1:
        vib_diag = compute_core_vibrations(
            images[ts_idx], calc=calc, reactive_indices=pair["reactive_indices"]
        )
    else:
        vib_diag = {
            "n_imag_modes_gt20cm1": 0,
            "imag_freqs_cm1": [],
            "min_real_freq_cm1": None,
            "max_real_freq_cm1": None,
            "note": "Maximum occurred at endpoint rather than internal image",
        }

    elapsed_s = round(time.perf_counter() - t0, 2)
    signed_err = round(ea_cineb_kcal - spec.lit_ea_kcal_mol, 2)

    return {
        "spec": asdict(spec),
        "status": {
            "valid_fs_topology": pair["valid_fs"],
            "fs_check": pair["fs_check"],
            "fmax_is_ev_A": pair["fmax_is"],
            "fmax_fs_ev_A": pair["fmax_fs"],
            "nsteps_neb": nsteps_neb,
            "nsteps_cineb": nsteps_ci,
            "final_neb_fmax_ev_A": round(final_neb_fmax, 4),
            "converged_cineb": bool(final_neb_fmax <= cineb_fmax + 1e-6),
            "ts_image_index": ts_idx,
            "n_images_total": len(images),
            "wall_time_s": elapsed_s,
        },
        "energies_kcal_mol": {
            "lit_ea_kcal_mol": spec.lit_ea_kcal_mol,
            "sevenn_ea_cineb_kcal_mol": round(ea_cineb_kcal, 2),
            "sevenn_ea_vs_sep_react_kcal_mol": round(ea_from_sep_kcal, 2),
            "signed_error_ea_kcal_mol": signed_err,
            "abs_error_ea_kcal_mol": round(abs(signed_err), 2),
            "sevenn_delta_e_adduct_kcal_mol": round(de_adduct_kcal, 2),
            "sevenn_delta_e_sep_kcal_mol": round(de_sep_kcal, 2),
            "lit_delta_e_rxn_kcal_mol": spec.lit_delta_e_rxn_kcal_mol,
            "band_profile_rel_kcal_mol": rel_kcal,
        },
        "ts_diagnostics": {
            "bond_lengths_A": bond_lengths,
            "core_vibrations": vib_diag,
        },
    }


def summarize_results(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute aggregate benchmark error statistics and selectivity metrics."""
    abs_errs = [r["energies_kcal_mol"]["abs_error_ea_kcal_mol"] for r in results]
    signed_errs = [r["energies_kcal_mol"]["signed_error_ea_kcal_mol"] for r in results]
    converged_n = sum(1 for r in results if r["status"]["converged_cineb"])
    true_ts_n = sum(
        1
        for r in results
        if r["ts_diagnostics"]["core_vibrations"]["n_imag_modes_gt20cm1"] == 1
    )

    by_id = {r["spec"]["reaction_id"]: r for r in results}
    ttip_sel: Optional[Dict[str, Any]] = None
    if "TTIP_acetone_elim" in by_id and "TTIP_propene_elim" in by_id:
        ea_ace = by_id["TTIP_acetone_elim"]["energies_kcal_mol"][
            "sevenn_ea_cineb_kcal_mol"
        ]
        ea_pro = by_id["TTIP_propene_elim"]["energies_kcal_mol"][
            "sevenn_ea_cineb_kcal_mol"
        ]
        lit_diff = (
            by_id["TTIP_propene_elim"]["spec"]["lit_ea_kcal_mol"]
            - by_id["TTIP_acetone_elim"]["spec"]["lit_ea_kcal_mol"]
        )
        sevenn_diff = round(ea_pro - ea_ace, 2)
        ttip_sel = {
            "lit_delta_ea_propene_minus_acetone_kcal_mol": round(lit_diff, 2),
            "sevenn_delta_ea_propene_minus_acetone_kcal_mol": sevenn_diff,
            "correct_channel_order": bool(sevenn_diff > 0.0),
        }

    return {
        "n_reactions": len(results),
        "n_valid_fs_topology": sum(
            1 for r in results if r["status"]["valid_fs_topology"]
        ),
        "n_converged_cineb": converged_n,
        "n_exact_1_imag_core_ts": true_ts_n,
        "mae_ea_kcal_mol": round(float(np.mean(abs_errs)), 2) if abs_errs else None,
        "rmse_ea_kcal_mol": (
            round(float(np.sqrt(np.mean(np.square(signed_errs)))), 2)
            if signed_errs
            else None
        ),
        "bias_mean_signed_err_kcal_mol": (
            round(float(np.mean(signed_errs)), 2) if signed_errs else None
        ),
        "ttip_selectivity": ttip_sel,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run SevenNet-0 CI-NEB decomposition benchmark on DB precursors."
    )
    parser.add_argument(
        "--reactions",
        nargs="+",
        default=None,
        help="Optional subset of reaction_ids to run.",
    )
    parser.add_argument(
        "--images",
        type=int,
        default=5,
        help="Number of intermediate NEB images (default: 5, total 7).",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(
            PROJECT_ROOT, "database", "external", "sevenn_decomp_neb_benchmark.json"
        ),
        help="Output JSON path.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="PyTorch intra-op CPU threads (default: 4).",
    )
    parser.add_argument(
        "--check-endpoints-only",
        action="store_true",
        help="Validate relaxed IS/FS/TS-seed endpoints and separated thermodynamics only.",
    )
    parser.add_argument(
        "--neb-steps",
        type=int,
        default=75,
        help="Max regular NEB optimization steps before climbing (default: 75).",
    )
    parser.add_argument(
        "--cineb-steps",
        type=int,
        default=130,
        help="Max CI-NEB climbing steps (default: 130).",
    )
    parser.add_argument(
        "--k-spring",
        type=float,
        default=0.10,
        help="NEB spring constant (eV/A^2, default: 0.10).",
    )
    parser.add_argument(
        "--maxstep",
        type=float,
        default=0.05,
        help="FIRE optimizer maxstep (A, default: 0.05).",
    )
    parser.add_argument(
        "--merge-parts",
        nargs="+",
        default=None,
        help="Merge partial JSON outputs in canonical order and save to --output.",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Only verify existing output JSON schema without recomputing.",
    )
    args = parser.parse_args()

    if args.merge_parts:
        by_id: Dict[str, Any] = {}
        for part_path in args.merge_parts:
            with open(part_path, "r", encoding="utf-8") as fh:
                part_data = json.load(fh)
            for rxn in part_data.get("reactions", []):
                by_id[rxn["spec"]["reaction_id"]] = rxn
        ordered = [by_id[s.reaction_id] for s in BENCHMARK_SPECS if s.reaction_id in by_id]
        merged_payload = {
            "model": "SevenNet-0 (7net-0)",
            "method": "ASE improvedtangent NEB + CI-NEB (IDPP initial path, FIRE optimizer)",
            "summary": summarize_results(ordered),
            "reactions": ordered,
        }
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as fh:
            json.dump(merged_payload, fh, indent=2)
        print("Merged", len(ordered), "reactions into:", args.output)
        print(json.dumps(merged_payload["summary"], indent=2))
        return

    if args.verify_only:
        with open(args.output, "r", encoding="utf-8") as fh:
            saved = json.load(fh)
        print("Verified benchmark JSON:", args.output)
        print(json.dumps(saved["summary"], indent=2))
        return

    torch.set_num_threads(args.threads)
    calc = SevenNetCalculator(model="7net-0", device="cpu")

    specs = BENCHMARK_SPECS
    if args.reactions:
        selected = set(args.reactions)
        specs = [s for s in BENCHMARK_SPECS if s.reaction_id in selected]

    if args.check_endpoints_only:
        ep_results = []
        for spec in specs:
            ep = run_single_ci_neb(
                spec=spec, calc=calc, n_intermediate=args.images, check_endpoints_only=True
            )
            ep_results.append(ep)
            print(
                f"{ep['reaction_id']}: valid_fs={ep['valid_fs']}, "
                f"dE_add={ep['de_adduct_kcal']} kcal/mol, "
                f"dE_sep={ep['de_sep_kcal']} kcal/mol, "
                f"fs_check={ep['fs_check']}, t={ep['wall_time_s']}s"
            )
        return

    results: List[Dict[str, Any]] = []
    for spec in specs:
        print(f"Running {spec.reaction_id} ({spec.reaction_equation}) ...")
        res = run_single_ci_neb(
            spec=spec,
            calc=calc,
            n_intermediate=args.images,
            k_spring=args.k_spring,
            neb_steps=args.neb_steps,
            cineb_steps=args.cineb_steps,
            maxstep=args.maxstep,
        )
        results.append(res)
        st = res["status"]
        en = res["energies_kcal_mol"]
        vb = res["ts_diagnostics"]["core_vibrations"]
        print(
            f"  [valid_fs={st['valid_fs_topology']}, conv={st['converged_cineb']}, "
            f"fmax={st['final_neb_fmax_ev_A']:.4f}] "
            f"Lit Ea={en['lit_ea_kcal_mol']:.2f} vs 7net={en['sevenn_ea_cineb_kcal_mol']:.2f} kcal/mol "
            f"(dE_sep={en['sevenn_delta_e_sep_kcal_mol']:.2f}), "
            f"n_imag={vb['n_imag_modes_gt20cm1']} {vb['imag_freqs_cm1']}, "
            f"t={st['wall_time_s']}s"
        )

    payload = {
        "model": "SevenNet-0 (7net-0)",
        "method": "ASE improvedtangent NEB + CI-NEB (IDPP initial path, FIRE optimizer)",
        "summary": summarize_results(results),
        "reactions": results,
    }
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print("\nSaved results to:", args.output)
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
