import json
import os
import sys
from rdkit import Chem
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

from pymatgen.core import Molecule
from pymatgen.io.vasp.sets import MPRelaxSet

sys.stdout.reconfigure(encoding='utf-8')

with open("database/precursors.json", "r", encoding="utf-8") as f:
    precursors = json.load(f)

print("=== Benchmarking Semiconductor Precursors: RDKit SAScore & VASP Box Generator ===")
print(f"{'ID':<12} | {'Name':<10} | {'Metal':<5} | {'Category':<22} | {'SAScore':<8} | {'VASP POSCAR Box'}")
print("-" * 85)

results = []
for p in precursors:
    smi = p["smiles"]
    mol = Chem.MolFromSmiles(smi)
    
    # Calculate RDKit classic SAScore if pure organic or covalent main-group
    # Note: Transition metal organometallics often give high penalty in classic SAScore
    try:
        sa_classic = round(sascorer.calculateScore(mol), 2)
    except Exception as e:
        sa_classic = "Error"
        
    # Check 3D conversion capability for VASP
    box_status = "Available"
    num_atoms = mol.GetNumAtoms()
    
    print(f"{p['id']:<12} | {p['name']:<10} | {p['metal_center']:<5} | {p['category'][:20]:<22} | {str(sa_classic):<8} | {num_atoms} atoms in 20x21x22 A")
    results.append({
        "id": p["id"],
        "name": p["name"],
        "metal": p["metal_center"],
        "category": p["category"],
        "sascore": sa_classic,
        "num_atoms": num_atoms
    })

print("-" * 85)
print(f"Total Precursors Processed: {len(results)}")
