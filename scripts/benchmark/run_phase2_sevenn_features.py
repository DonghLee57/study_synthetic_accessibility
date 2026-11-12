#!/usr/bin/env python3
"""
Phase 2 (Refined): SevenNet Atomistic Physical Features for Precursor SA Evaluation
-----------------------------------------------------------------------------------
Focused strictly on Precursor Physical SA Pillars (excluding surface selectivity):
 1. Cohesion energy per atom (e_coh, eV/atom)
 2. Minimum bond dissociation/stretch energy (Delta E_stretch^min, kcal/mol)
 3. 3D Steric Buried Volume (%V_bur, %) via Cavallo SambVca sphere (R = 3.5 A)
 4. Monomer vs Dimerization interaction energy (Delta E_dimer, kcal/mol)
 5. Thermodynamic complexation feasibility proxy (Delta E_complex, kcal/mol)
 6. Molecular radius of gyration (Rg, A) & van der Waals volume (V_mol, A^3)

Saves calculated physical features to database/external/sevenn_physical_features.json.
"""

import os
import sys
import json
import time
from collections import deque
import numpy as np
import torch

# Force unbuffered output so progress prints immediately
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)

from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import LBFGS
from sevenn.calculator import SevenNetCalculator
from rdkit import Chem
from rdkit.Chem import AllChem

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Isolated single-atom reference energies in SevenNet-0 (eV)
ISOLATED_ATOM_ENERGIES = {
    'H': -0.9187,
    'C': -2.0688,
    'N': -2.8482,
    'O': -1.9414,
    'F': -1.6172,
    'Si': -1.3102,
    'P': -1.9107,
    'S': -0.8525,
    'Cl': -0.6802,
    'Ti': -2.7495,
    'Al': -0.8508,
    'Cu': -0.6695,
    'Zr': -2.7250,
    'Ru': -2.6209,
    'Hf': -3.6573,
    'Ta': -3.9477,
    'W': -4.2950,
    'Co': -1.6691,
    'Mo': -4.5894,
    'Sn': -0.9625,
    'Br': -0.6409,
    'I': -0.6888
}

BONDI_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'Si': 2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Ti': 2.00,
    'Al': 2.00, 'Zr': 2.15, 'Hf': 2.15, 'Ru': 2.05, 'Co': 2.00,
    'W': 2.10, 'Mo': 2.10, 'Ta': 2.15, 'Sn': 2.17, 'Cu': 2.00,
    'Br': 1.85, 'I': 1.98
}

METAL_SET = {'Ti', 'Zr', 'Hf', 'Al', 'Ru', 'Co', 'W', 'Mo', 'Ta', 'Sn', 'Cu', 'Si'}

def find_central_atom_idx(mol):
    """Identify the central metal or metalloid atom in precursor."""
    atoms = mol.GetAtoms()
    for a in atoms:
        sym = a.GetSymbol()
        if sym in METAL_SET and sym != 'Si':
            return a.GetIdx(), sym
    for a in atoms:
        if a.GetSymbol() == 'Si':
            return a.GetIdx(), 'Si'
    max_deg = -1
    best_idx = 0
    best_sym = 'C'
    for a in atoms:
        deg = a.GetDegree()
        if deg > max_deg:
            max_deg = deg
            best_idx = a.GetIdx()
            best_sym = a.GetSymbol()
    return best_idx, best_sym

def calc_buried_volume(mol, center_idx, radius=3.5, grid_step=0.25):
    """
    Compute 3D Steric Buried Volume (%V_bur) around the central atom
    using Cavallo SambVca spherical integration.
    """
    conf = mol.GetConformer()
    pos = conf.GetPositions()
    c_pos = pos[center_idx]
    
    x = np.arange(-radius, radius + grid_step, grid_step)
    y = np.arange(-radius, radius + grid_step, grid_step)
    z = np.arange(-radius, radius + grid_step, grid_step)
    gx, gy, gz = np.meshgrid(x, y, z, indexing='ij')
    points = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=-1)
    
    dists = np.linalg.norm(points, axis=1)
    sphere_mask = dists <= radius
    sphere_points = points[sphere_mask] + c_pos
    total_pts = len(sphere_points)
    
    occupied = np.zeros(total_pts, dtype=bool)
    for i, a in enumerate(mol.GetAtoms()):
        if i == center_idx:
            continue
        sym = a.GetSymbol()
        rvdw = BONDI_RADII.get(sym, 1.70)
        a_pos = pos[i]
        d_pts = np.linalg.norm(sphere_points - a_pos, axis=1)
        occupied |= (d_pts <= rvdw)
        
    vbur = (np.sum(occupied) / total_pts) * 100.0
    return round(float(vbur), 2)

def align_principal_axes(positions: np.ndarray) -> np.ndarray:
    """Center and rotate coordinates so principal axis z has minimal spatial extent."""
    pos_c = positions - np.mean(positions, axis=0)
    if len(pos_c) <= 1:
        return pos_c
    cov = np.cov(pos_c.T)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    rot = evecs[:, order]
    if np.linalg.det(rot) < 0:
        rot[:, 2] *= -1
    return pos_c @ rot


def relax_atoms_sevenn(
    symbols: list,
    positions: np.ndarray,
    calc: SevenNetCalculator,
    fmax: float = 0.04,
    steps: int = 100,
    cell_len: float = 50.0,
    constraint=None,
):
    """Run SevenNet-0 LBFGS geometry relaxation and return (energy_ev, centered_pos, nsteps, final_fmax)."""
    pos_c = positions - np.mean(positions, axis=0)
    atoms = Atoms(
        symbols=symbols,
        positions=pos_c,
        cell=[cell_len, cell_len, cell_len],
        pbc=[False, False, False],
    )
    if constraint is not None:
        atoms.set_constraint(constraint)
    atoms.calc = calc
    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    energy_ev = float(atoms.get_potential_energy())
    forces = atoms.get_forces()
    final_fmax = float(np.max(np.linalg.norm(forces, axis=1)))
    rel_pos = atoms.get_positions() - np.mean(atoms.get_positions(), axis=0)
    return energy_ev, rel_pos, int(opt.nsteps), round(final_fmax, 4)


def split_acyclic_bond(mol: Chem.Mol, a1: int, a2: int):
    """Partition atom indices into two disjoint connected subtrees across acyclic bond (a1, a2)."""
    adj = {i: [] for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (u == a1 and v == a2) or (u == a2 and v == a1):
            continue
        adj[u].append(v)
        adj[v].append(u)
    visited = {a1}
    queue = deque([a1])
    while queue:
        curr = queue.popleft()
        for nbr in adj[curr]:
            if nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)
    frag1 = sorted(visited)
    frag2 = sorted(set(range(mol.GetNumAtoms())) - visited)
    return frag1, frag2


def canonical_fragment_key(mol: Chem.Mol, indices: list, symbols: list):
    """Build deterministic topological cache key for an isolated fragment subgraph."""
    try:
        idx_set = set(indices)
        idx_map = {old_i: new_i for new_i, old_i in enumerate(indices)}
        rw = Chem.RWMol()
        for old_i in indices:
            rw.AddAtom(Chem.Atom(symbols[old_i]))
        for bond in mol.GetBonds():
            u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if u in idx_set and v in idx_set:
                rw.AddBond(idx_map[u], idx_map[v], bond.GetBondType())
        smi = Chem.MolToSmiles(rw.GetMol(), canonical=True)
        return (smi, tuple(sorted(symbols[i] for i in indices)))
    except Exception:
        return None


def compute_dimerization_energy(symbols: list, positions_aligned: np.ndarray, e_mono: float, calc: SevenNetCalculator):
    """
    Compute relaxed dimerization interaction energy (Delta E_dimer, kcal/mol).
    Guarantees non-overlapping separation along z via Delta_z = z_span + d_vdW (d_vdW >= 2.8 A).
    """
    vdw_gap = 2.8  # units: A
    pos1 = positions_aligned.copy()
    pos2 = positions_aligned.copy()
    z_span = float(np.max(pos1[:, 2]) - np.min(pos1[:, 2]))
    z_shift = z_span + vdw_gap
    pos2[:, 0] *= -1
    pos2[:, 1] *= -1
    pos2[:, 2] += z_shift

    init_gap = float(np.min(np.linalg.norm(pos1[:, None, :] - pos2[None, :, :], axis=-1)))
    e_dimer, _, dim_steps, _ = relax_atoms_sevenn(
        symbols=symbols + symbols,
        positions=np.vstack([pos1, pos2]),
        calc=calc,
        fmax=0.04,
        steps=35,
        cell_len=60.0,
    )
    dE_kcal = (e_dimer - 2.0 * e_mono) * 23.0605  # units: kcal/mol
    return round(float(dE_kcal), 2), round(init_gap, 3), dim_steps


def compute_precursor_physical_features(smi: str, calc: SevenNetCalculator, frag_cache: dict = None):
    """Evaluate SevenNet-0 BFGS/LBFGS relaxed conformer, fragment BDE, %V_bur, and dimer energy."""
    if frag_cache is None:
        frag_cache = {}
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return None
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    embed_res = AllChem.EmbedMolecule(mol, params)
    if embed_res < 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    positions = conf.GetPositions()

    # 1. Relax monomer geometry with SevenNet-0 LBFGS
    try:
        e_eq, pos_rel, mono_steps, mono_fmax = relax_atoms_sevenn(
            symbols=symbols,
            positions=positions,
            calc=calc,
            fmax=0.04,
            steps=100,
            cell_len=50.0,
        )
        pos_aligned = align_principal_axes(pos_rel)
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, pos_aligned[i].tolist())
    except Exception:
        e_eq = np.nan
        pos_aligned = align_principal_axes(positions)
        mono_steps = 0
        mono_fmax = 999.0

    # 2. Cohesion energy per atom (e_coh, eV/atom)
    if not np.isnan(e_eq):
        sum_isolated = sum(ISOLATED_ATOM_ENERGIES.get(sym, -1.0) for sym in symbols)
        binding_e = sum_isolated - e_eq
        e_coh = float(binding_e / len(symbols))
    else:
        e_coh = 4.10

    # 3. Identify central metal/metalloid atom
    center_idx, center_sym = find_central_atom_idx(mol)

    # 4. 3D Steric Buried Volume (%V_bur) on SevenNet-relaxed coordinates
    vbur = calc_buried_volume(mol, center_idx, radius=3.5, grid_step=0.25)

    # 5. Bond dissociation energy (BDE_min) via symmetry-deduplicated core bond cleavage:
    #    - Acyclic bonds: true relaxed complementary fragment dissociation E(ML_{n-1}) + E(L) - E(ML_n)
    #    - Cyclic bonds: constrained FixAtoms ring stretch after exocyclic subtree displacement
    bdes = []
    bde_modes = []
    if not np.isnan(e_eq):
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        seen_sigs = set()
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            sym1, sym2 = symbols[a1], symbols[a2]
            if sym1 == 'H' or sym2 == 'H':
                continue
            is_core_bond = (a1 == center_idx or a2 == center_idx or sym1 in METAL_SET or sym2 in METAL_SET)
            if not is_core_bond and not (sym1 in {'N', 'O', 'S', 'Cl'} and sym2 in {'N', 'O', 'S', 'Cl'}):
                continue
            sig = tuple(sorted([ranks[a1], ranks[a2]]))
            if sig in seen_sigs:
                continue
            seen_sigs.add(sig)

            if not b.IsInRing():
                f1, f2 = split_acyclic_bond(mol, a1, a2)
                e_frags = []
                ok = True
                for f_idx in (f1, f2):
                    fkey = canonical_fragment_key(mol, f_idx, symbols)
                    if fkey is not None and fkey in frag_cache:
                        e_frags.append(frag_cache[fkey])
                    else:
                        try:
                            ef, _, _, _ = relax_atoms_sevenn(
                                symbols=[symbols[i] for i in f_idx],
                                positions=pos_aligned[f_idx],
                                calc=calc,
                                fmax=0.04,
                                steps=45,
                                cell_len=45.0,
                            )
                            if fkey is not None:
                                frag_cache[fkey] = ef
                            e_frags.append(ef)
                        except Exception:
                            ok = False
                            break
                if ok and len(e_frags) == 2:
                    dE_kcal = (e_frags[0] + e_frags[1] - e_eq) * 23.0605  # units: kcal/mol
                    bdes.append(float(dE_kcal))
                    bde_modes.append("fragment_dissociation")
            else:
                new_pos = pos_aligned.copy()
                vec = new_pos[a2] - new_pos[a1]
                dist0 = float(np.linalg.norm(vec))
                if dist0 < 0.1:
                    continue
                unit_v = vec / dist0
                # Shift a2 along with its non-ring exocyclic branch to prevent local steric collision
                move_set = {a2}
                q_exo = deque()
                for nbr in mol.GetAtomWithIdx(a2).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni != a1 and not nbr.IsInRing():
                        move_set.add(ni)
                        q_exo.append(ni)
                while q_exo:
                    curr = q_exo.popleft()
                    for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                        ni = nbr.GetIdx()
                        if ni not in move_set and ni != a1 and not nbr.IsInRing():
                            move_set.add(ni)
                            q_exo.append(ni)
                for m_idx in move_set:
                    new_pos[m_idx] += unit_v * 1.2
                try:
                    e_str, _, _, _ = relax_atoms_sevenn(
                        symbols=symbols,
                        positions=new_pos,
                        calc=calc,
                        fmax=0.05,
                        steps=25,
                        cell_len=50.0,
                        constraint=FixAtoms(indices=[a1, a2]),
                    )
                    dE_kcal = (e_str - e_eq) * 23.0605  # units: kcal/mol
                    bdes.append(float(dE_kcal))
                    bde_modes.append("constrained_ring_stretch")
                except Exception:
                    pass
    bde_min = float(min(bdes)) if bdes else 60.0
    if "fragment_dissociation" in bde_modes and "constrained_ring_stretch" in bde_modes:
        bde_mode_label = "hybrid_fragment_and_ring"
    elif bde_modes:
        bde_mode_label = bde_modes[0]
    else:
        bde_mode_label = "default_unreactive"

    # 6. Monomer vs Dimerization interaction energy (Delta E_dimer) via SevenNet-0
    is_fallback = bool(np.isnan(e_eq))
    if not is_fallback:
        try:
            dE_dimer, dimer_init_gap, dimer_steps = compute_dimerization_energy(symbols, pos_aligned, e_eq, calc)
        except Exception:
            dE_dimer = -5.0
            dimer_init_gap = 0.0
            dimer_steps = 0
            is_fallback = True
    else:
        dE_dimer = 15.0
        dimer_init_gap = 0.0
        dimer_steps = 0

    # 7. Molecular volume & Radius of gyration on relaxed coordinates
    masses = [mol.GetAtomWithIdx(i).GetMass() for i in range(mol.GetNumAtoms())]
    total_mass = sum(masses)
    rg = float(np.sqrt(sum(m * np.sum(r**2) for m, r in zip(masses, pos_aligned)) / total_mass))
    vol_mol = float(4.0 / 3.0 * np.pi * (rg + 1.8)**3)

    # 8. Empirical Algebraic Complexation Proxy (Delta E_complex_proxy, kcal/mol)
    # NOTE: Derived deterministically from %V_bur, dE_dimer, and bde_min; not an explicit salt-metathesis DFT/MLIP calculation.
    steric_penalty = max(0.0, (vbur - 82.0) * 3.5) if vbur > 82.0 else 0.0
    dimer_penalty = max(0.0, -dE_dimer * 2.0) if dE_dimer < 0.0 else 0.0
    e_complex_proxy = round(float(-25.0 + steric_penalty + dimer_penalty + max(0.0, (50.0 - bde_min) * 1.5)), 2)

    return {
        "symbols": symbols,
        "center_atom": center_sym,
        "is_fallback": is_fallback,
        "mono_relax_steps": mono_steps,
        "mono_fmax_ev_ang": mono_fmax,
        "dimer_init_gap_ang": dimer_init_gap,
        "dimer_relax_steps": dimer_steps,
        "bde_mode": bde_mode_label,
        "e_eq": round(e_eq, 3) if not np.isnan(e_eq) else -100.0,
        "e_coh": round(e_coh, 3),
        "bde_min": round(bde_min, 2),
        "vbur_percent": vbur,
        "delta_e_dimer_kcal_mol": dE_dimer,
        "delta_e_complex_kcal_mol": e_complex_proxy,
        "rg": round(rg, 2),
        "vol_mol": round(vol_mol, 1)
    }


def run_recalculation():
    print("=" * 80)
    print("PHASE 2 RECALCULATION: RELAXED SEVENNET-0 + NON-OVERLAPPING DIMER + FRAGMENT BDE")
    print("=" * 80)

    torch.set_num_threads(4)
    phase1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
    with open(phase1_path, "r", encoding="utf-8") as f:
        phase1_data = json.load(f)
    molecules = phase1_data["results"]
    print(f"Loaded {len(molecules)} molecules from Phase 1 baseline (torch threads={torch.get_num_threads()}).")

    print("Initializing SevenNet-0 Calculator on CPU...")
    calc = SevenNetCalculator(model='7net-0', device='cpu')

    features_list = []
    frag_cache = {}
    t_start = time.time()

    for idx, mol_item in enumerate(molecules, 1):
        m_id = mol_item["id"]
        name = mol_item["name"]
        smi = mol_item["smiles"]
        t_type = mol_item["target_type"]
        is_valid = mol_item["is_fab_valid"]

        t0 = time.time()
        p_data = compute_precursor_physical_features(smi, calc, frag_cache=frag_cache)
        t_el = time.time() - t0

        is_monomer_volatile = (60.0 <= p_data["vbur_percent"] <= 82.0) and (p_data["delta_e_dimer_kcal_mol"] >= -5.0)
        is_thermally_labile = bool(p_data["bde_min"] < 40.0 or p_data["e_coh"] < 3.0)

        record = {
            "id": m_id,
            "name": name,
            "smiles": smi,
            "formula": mol_item["formula"],
            "domain": mol_item["domain"],
            "target_type": t_type,
            "is_fab_valid": is_valid,
            "center_atom": p_data["center_atom"],
            "is_fallback": p_data["is_fallback"],
            "mono_relax_steps": p_data["mono_relax_steps"],
            "mono_fmax_ev_ang": p_data["mono_fmax_ev_ang"],
            "dimer_init_gap_ang": p_data["dimer_init_gap_ang"],
            "dimer_relax_steps": p_data["dimer_relax_steps"],
            "bde_mode": p_data["bde_mode"],
            "e_pot_ev": p_data["e_eq"],
            "e_coh_ev_per_atom": p_data["e_coh"],
            "bde_min_kcal_mol": p_data["bde_min"],
            "vbur_percent": p_data["vbur_percent"],
            "delta_e_dimer_kcal_mol": p_data["delta_e_dimer_kcal_mol"],
            "delta_e_complex_kcal_mol": p_data["delta_e_complex_kcal_mol"],
            "radius_of_gyration_ang": p_data["rg"],
            "mol_volume_ang3": p_data["vol_mol"],
            "is_monomer_volatile": is_monomer_volatile,
            "is_thermally_labile": is_thermally_labile
        }
        features_list.append(record)
        print(
            f"[{idx:02d}/55] {m_id:10s} | {name:26s} | %Vbur:{p_data['vbur_percent']:5.1f}% | "
            f"gap:{p_data['dimer_init_gap_ang']:.2f}A dE_dim:{p_data['delta_e_dimer_kcal_mol']:6.1f} | "
            f"BDE:{p_data['bde_min']:5.1f} | st:{p_data['mono_relax_steps']:3d} | fb:{p_data['is_fallback']} | {t_el:.1f}s",
            flush=True
        )

    t_total = time.time() - t_start
    fallback_cnt = sum(1 for r in features_list if r["is_fallback"])
    print(f"\nSuccessfully recalculated physical features for 55 molecules in {t_total:.1f}s (fallbacks={fallback_cnt}, cached_frags={len(frag_cache)}).")

    out_dir = os.path.join(PROJECT_ROOT, "database", "external")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "sevenn_physical_features.json")

    out_payload = {
        "metadata": {
            "framework": "SevenNet (7net-0 LBFGS Relaxed) + Cavallo SambVca Sphere",
            "relaxation_method": "ASE LBFGS (Monomer fmax=0.04 eV/A, <=100 steps; Dimer fmax=0.04 eV/A, <=35 steps)",
            "dimer_shift_rule": "z_shift = z_span + 2.8 A along thinnest principal axis (non-overlapping gap >= 2.8 A)",
            "bde_method": "Relaxed separated fragment dissociation E(ML_{n-1})+E(L)-E(ML_n) for acyclic bonds; FixAtoms exocyclic-shifted stretch for cyclic bonds",
            "pillars": [
                "1. Bond Dissociation Energy (BDE_min, kcal/mol)",
                "2. Cohesion Energy (e_coh, eV/atom)",
                "3. 3D Steric Buried Volume (%V_bur, %)",
                "4. Monomer vs Dimerization Interaction (Delta E_dimer, kcal/mol)",
                "5. Complexation Feasibility Proxy (Delta E_complex, kcal/mol)"
            ],
            "total_molecules_computed": len(features_list),
            "fallback_count": fallback_cnt,
            "computation_time_seconds": round(t_total, 1)
        },
        "features": features_list
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_payload, f, indent=2, ensure_ascii=False)

    print(f"Saved refined physical features to: {out_path}")

if __name__ == "__main__":
    run_recalculation()
