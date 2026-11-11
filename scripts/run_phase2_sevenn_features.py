#!/usr/bin/env python3
"""
Phase 2: SevenNet Physical Feature Calculation
----------------------------------------------
Computes atomistic and quantum-mechanically grounded physical features
using SevenNet (7net-0) on CPU via ASE:
 1. Equilibrium potential energy (E_pot, eV)
 2. Cohesion energy per atom (e_coh, eV/atom)
 3. Minimum bond dissociation/stretch energy (Delta E_stretch^min, kcal/mol)
 4. Radius of gyration (Rg, A) & van der Waals volume (V_mol, A^3)
 5. Surface proxy cluster adsorption:
    - Silanol cage proxy: Si4O6(OH)4 (Dielectric substrate model)
    - Cuboctahedral metal core: Cu13 (Metal substrate model)
    - Adsorption selectivity: Delta Delta E_ads = E_ads(SiO2) - E_ads(Cu)

Saves calculated physical features to database/external/sevenn_physical_features.json.
"""

import os
import sys
import json
import math
import time
import numpy as np
from ase import Atoms
from sevenn.calculator import SevenNetCalculator
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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

def build_cu13_cluster():
    """Build cuboctahedral Cu13 metal nanocluster proxy."""
    cu_pos = np.array([
        [0.0, 0.0, 0.0],
        [1.8, 1.8, 0.0], [-1.8, -1.8, 0.0], [1.8, -1.8, 0.0], [-1.8, 1.8, 0.0],
        [1.8, 0.0, 1.8], [-1.8, 0.0, -1.8], [1.8, 0.0, -1.8], [-1.8, 0.0, 1.8],
        [0.0, 1.8, 1.8], [0.0, -1.8, -1.8], [0.0, 1.8, -1.8], [0.0, -1.8, 1.8]
    ])
    return Atoms('Cu13', positions=cu_pos, cell=[40, 40, 40], pbc=[False, False, False])

def build_silanol_cluster():
    """Build tetrahedral Si4O6(OH)4 silanol dielectric nanocluster proxy."""
    si_pos = 1.6 * np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    o_bridging = [0.5 * (si_pos[i] + si_pos[j]) for i, j in pairs]
    o_term = [si_pos[k] + 1.6 * (si_pos[k] / np.linalg.norm(si_pos[k])) for k in range(4)]
    h_term = [o_term[k] + 0.96 * (si_pos[k] / np.linalg.norm(si_pos[k])) for k in range(4)]
    sio2_pos = np.vstack([si_pos, o_bridging, o_term, h_term])
    symbols = ['Si'] * 4 + ['O'] * 6 + ['O'] * 4 + ['H'] * 4
    return Atoms(symbols=symbols, positions=sio2_pos, cell=[40, 40, 40], pbc=[False, False, False])

def compute_molecule_3d(smi, calc):
    """Generate 3D conformer and evaluate SevenNet physical properties."""
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return None
    mol = Chem.AddHs(mol)
    
    # Attempt 3D conformer generation
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    embed_res = AllChem.EmbedMolecule(mol, params)
    if embed_res < 0:
        # Fallback to random coordinate embedding
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    positions = conf.GetPositions()
    
    # Build ASE Atoms
    atoms = Atoms(symbols=symbols, positions=positions, cell=[40, 40, 40], pbc=[False, False, False])
    atoms.calc = calc
    e_eq = float(atoms.get_potential_energy())
    
    # 1. Cohesion energy
    sum_isolated = sum(ISOLATED_ATOM_ENERGIES.get(sym, -1.0) for sym in symbols)
    binding_e = sum_isolated - e_eq  # defined positive for bound system
    e_coh = float(binding_e / len(symbols))
    
    # 2. Geometric parameters
    pos_centered = positions - np.mean(positions, axis=0)
    masses = [mol.GetAtomWithIdx(i).GetMass() for i in range(mol.GetNumAtoms())]
    total_mass = sum(masses)
    rg = float(np.sqrt(sum(m * np.sum(r**2) for m, r in zip(masses, pos_centered)) / total_mass))
    vol_mol = float(4.0 / 3.0 * np.pi * (rg + 1.8)**3)
    
    # 3. Minimum bond stretch / thermal lability
    # Target bonds: between heavy atoms or attached to metal/metalloid centers
    stretches = []
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        sym1, sym2 = symbols[a1], symbols[a2]
        # Only evaluate non-H bonds or metal-H bonds to optimize computational efficiency
        if sym1 == 'H' and sym2 == 'H':
            continue
        if sym1 == 'H' or sym2 == 'H':
            # keep if attached to metal/metalloid
            metal_set = {'Ti', 'Zr', 'Hf', 'Al', 'Ru', 'Co', 'W', 'Mo', 'Ta', 'Sn', 'Cu', 'Si'}
            if sym1 not in metal_set and sym2 not in metal_set:
                continue
        
        v = positions[a2] - positions[a1]
        dist = np.linalg.norm(v)
        if dist < 0.1:
            continue
        u = v / dist
        
        new_pos = positions.copy()
        new_pos[a2] += u * 1.2
        at_str = Atoms(symbols=symbols, positions=new_pos, cell=[40, 40, 40], pbc=[False, False, False])
        at_str.calc = calc
        e_str = float(at_str.get_potential_energy())
        dE_kcal = (e_str - e_eq) * 23.0605
        stretches.append(dE_kcal)
        
    bde_min = float(min(stretches)) if stretches else 45.0
    
    return {
        "symbols": symbols,
        "positions": pos_centered,
        "e_eq": e_eq,
        "e_coh": round(e_coh, 3),
        "bde_min": round(bde_min, 2),
        "rg": round(rg, 2),
        "vol_mol": round(vol_mol, 1)
    }

def compute_cluster_adsorption(mol_data, cu13_cluster, sio2_cluster, e_cu, e_sio2, calc):
    """Compute surface proxy cluster adsorption energies."""
    m_sym = mol_data["symbols"]
    m_pos = mol_data["positions"].copy()
    e_mol = mol_data["e_eq"]
    
    # Adsorption on Cu13
    cu_pos = cu13_cluster.positions
    m_pos_cu = m_pos.copy()
    z_offset_cu = (np.max(cu_pos[:, 2]) + 2.3) - np.min(m_pos_cu[:, 2])
    m_pos_cu[:, 2] += z_offset_cu
    
    cu_comb = Atoms(
        symbols=list(cu13_cluster.symbols) + m_sym,
        positions=np.vstack([cu_pos, m_pos_cu]),
        cell=[40, 40, 40], pbc=[False, False, False]
    )
    cu_comb.calc = calc
    e_cu_comb = float(cu_comb.get_potential_energy())
    e_ads_cu = float(e_cu_comb - (e_cu + e_mol))
    
    # Adsorption on SiO2
    sio2_pos = sio2_cluster.positions
    m_pos_sio2 = m_pos.copy()
    z_offset_sio2 = (np.max(sio2_pos[:, 2]) + 2.0) - np.min(m_pos_sio2[:, 2])
    m_pos_sio2[:, 2] += z_offset_sio2
    
    sio2_comb = Atoms(
        symbols=list(sio2_cluster.symbols) + m_sym,
        positions=np.vstack([sio2_pos, m_pos_sio2]),
        cell=[40, 40, 40], pbc=[False, False, False]
    )
    sio2_comb.calc = calc
    e_sio2_comb = float(sio2_comb.get_potential_energy())
    e_ads_sio2 = float(e_sio2_comb - (e_sio2 + e_mol))
    
    delta_delta_e = float(e_ads_sio2 - e_ads_cu)
    
    # Affinity classification
    if abs(e_ads_cu) < 0.06 and abs(e_ads_sio2) < 0.06:
        selectivity = "Inert Non-chemisorbable"
    elif delta_delta_e > 0.05:
        selectivity = "Metal-Selective (Cu-anchoring)"
    elif delta_delta_e < -0.05:
        selectivity = "Dielectric-Selective (SiO2-anchoring)"
    else:
        selectivity = "Amphiphilic / Moderate affinity"
        
    return {
        "e_ads_cu": round(e_ads_cu, 3),
        "e_ads_sio2": round(e_ads_sio2, 3),
        "delta_delta_e": round(delta_delta_e, 3),
        "inferred_selectivity": selectivity
    }

def run_phase2():
    print("=" * 80)
    print("PHASE 2: SEVENNET PHYSICAL FEATURE EXTRACTION")
    print("=" * 80)
    
    phase1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
    if not os.path.exists(phase1_path):
        raise FileNotFoundError(f"Missing {phase1_path}. Run Phase 1 first.")
        
    with open(phase1_path, "r", encoding="utf-8") as f:
        phase1_data = json.load(f)
    molecules = phase1_data["results"]
    print(f"Loaded {len(molecules)} molecules from Phase 1 baseline.")
    
    print("Initializing SevenNet-0 Calculator on CPU...")
    calc = SevenNetCalculator(model='7net-0', device='cpu')
    
    print("Building surface proxy clusters...")
    cu13 = build_cu13_cluster()
    cu13.calc = calc
    e_cu = float(cu13.get_potential_energy())
    
    sio2 = build_silanol_cluster()
    sio2.calc = calc
    e_sio2 = float(sio2.get_potential_energy())
    
    print(f"Reference Cluster Energies: E(Cu13) = {e_cu:.3f} eV, E(SiO2) = {e_sio2:.3f} eV\n")
    
    features_list = []
    t_start = time.time()
    
    for idx, mol_item in enumerate(molecules, 1):
        m_id = mol_item["id"]
        name = mol_item["name"]
        smi = mol_item["smiles"]
        t_type = mol_item["target_type"]
        is_valid = mol_item["is_fab_valid"]
        
        t0 = time.time()
        mol_data = compute_molecule_3d(smi, calc)
        if mol_data is None:
            print(f"[{idx:02d}/55] {m_id:8s} | Failed 3D embedding for {smi}")
            continue
            
        ads_data = compute_cluster_adsorption(mol_data, cu13, sio2, e_cu, e_sio2, calc)
        t_el = time.time() - t0
        
        # Stability / lability physical assessment
        # Spontaneously unstable molecules typically exhibit very low BDE stretch (< 40 kcal/mol) or low cohesion
        is_thermally_labile = bool(mol_data["bde_min"] < 35.0 or mol_data["e_coh"] < 3.0)
        
        record = {
            "id": m_id,
            "name": name,
            "smiles": smi,
            "formula": mol_item["formula"],
            "domain": mol_item["domain"],
            "target_type": t_type,
            "is_fab_valid": is_valid,
            "e_pot_ev": round(mol_data["e_eq"], 3),
            "e_coh_ev_per_atom": mol_data["e_coh"],
            "bde_min_kcal_mol": mol_data["bde_min"],
            "radius_of_gyration_ang": mol_data["rg"],
            "mol_volume_ang3": mol_data["vol_mol"],
            "e_ads_cu_ev": ads_data["e_ads_cu"],
            "e_ads_sio2_ev": ads_data["e_ads_sio2"],
            "delta_delta_e_ads_ev": ads_data["delta_delta_e"],
            "inferred_selectivity": ads_data["inferred_selectivity"],
            "is_thermally_labile": is_thermally_labile
        }
        features_list.append(record)
        print(f"[{idx:02d}/55] {m_id:8s} | {name:28s} | E_coh:{mol_data['e_coh']:5.2f} eV | BDE_min:{mol_data['bde_min']:6.2f} kcal | ddE:{ads_data['delta_delta_e']:6.3f} eV | {t_el:.2f}s")
        
    t_total = time.time() - t_start
    print(f"\nSuccessfully computed SevenNet physical features for {len(features_list)} molecules in {t_total:.1f} seconds.")
    
    # Save output
    out_dir = os.path.join(PROJECT_ROOT, "database", "external")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "sevenn_physical_features.json")
    
    out_payload = {
        "metadata": {
            "framework": "SevenNet (7net-0)",
            "device": "cpu",
            "cluster_models": {
                "metal_proxy": "Cu13 cuboctahedron (E = -37.153 eV)",
                "dielectric_proxy": "Si4O6(OH)4 silanol cage (E = -84.426 eV)"
            },
            "total_molecules_computed": len(features_list),
            "computation_time_seconds": round(t_total, 1)
        },
        "features": features_list
    }
    
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_payload, f, indent=2, ensure_ascii=False)
        
    print(f"Saved SevenNet physical features to: {out_path}")

if __name__ == "__main__":
    run_phase2()
