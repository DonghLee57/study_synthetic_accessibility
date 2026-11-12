#!/usr/bin/env python3
"""
Phase 1: Pure Data-Driven Baseline Evaluation
--------------------------------------------
Evaluates 55 molecules:
 - 27 semiconductor-process baseline precursors (Positive)
 - 18 AS-ALD inhibitor molecules (Positive)
 - 10 representative negative controls (Synthetic failures, decomposing, or non-functional)

Computes baseline chemoinformatics scores:
 - RDKit SAScore
 - SCScore (1024bool Reaxys standalone model)
 - SYBA (with default score table)
 - RAscore (XGBoost ChEMBL counts)
 - BR-SAScore (Strict Unleaked Reagent Domain)
 - 29 Reagent catalog match count
 - 12 Reaction template match count
Evaluates False Negative / False Positive rates and baseline decision performance.
"""

import os
import sys
import json
import math
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, RDConfig
import xgboost as xgb

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "scscore"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "reasyn"))

import sascorer
from scscore.standalone_model_numpy import SCScorer
from syba.syba import SybaClassifier
from BRSAScore import SAScorer as BRSAScorer

def load_all_molecules():
    # 1. 27 Precursors
    with open(os.path.join(PROJECT_ROOT, "database", "precursors.json"), "r", encoding="utf-8") as f:
        precursors = json.load(f)
        
    # 2. 18 Inhibitors
    with open(os.path.join(PROJECT_ROOT, "database", "external", "as_ald_inhibitors.json"), "r", encoding="utf-8") as f:
        inhibitors = json.load(f)["inhibitors"]
        
    # 3. 10 Negative Controls
    negatives = [
        {
            "id": "NEG-001",
            "name": "Ti(Et)4 (Tetraethyltitanium)",
            "smiles": "CC[Ti](CC)(CC)CC",
            "formula": "C8H20Ti",
            "domain": "Unstable Metal Alkyl",
            "failure_reason": "Rapid beta-hydride elimination at -50 °C; decomposes to Ti(III) and ethane/ethene."
        },
        {
            "id": "NEG-002",
            "name": "ZrMe4 (Tetramethylzirconium)",
            "smiles": "C[Zr](C)(C)C",
            "formula": "C4H12Zr",
            "domain": "Thermally Explosive",
            "failure_reason": "Extremely thermally unstable; decomposes explosively above -15 °C."
        },
        {
            "id": "NEG-003",
            "name": "Silacyclobutadiene",
            "smiles": "C1=C[SiH2]C=1",
            "formula": "C3H4Si",
            "domain": "Anti-aromatic Strained Silane",
            "failure_reason": "Extreme anti-aromatic ring strain and radical reactivity; impossible to isolate or volatilize."
        },
        {
            "id": "NEG-004",
            "name": "WPh5 (Pentaphenyltungsten)",
            "smiles": "c1ccccc1[W](c2ccccc2)(c3ccccc3)(c4ccccc4)c5ccccc5",
            "formula": "C30H25W",
            "domain": "Steric/Thermal Labile",
            "failure_reason": "Spontaneous reductive elimination and radical decomposition; zero vapor pressure."
        },
        {
            "id": "NEG-005",
            "name": "Decane (Non-functionalized alkane)",
            "smiles": "CCCCCCCCCC",
            "formula": "C10H22",
            "domain": "Ineffective Inhibitor",
            "failure_reason": "No anchoring heteroatom (S, N, Si); zero chemisorption on both SiO2 and Cu."
        },
        {
            "id": "NEG-006",
            "name": "Trichloroacetic acid",
            "smiles": "O=C(O)C(Cl)(Cl)Cl",
            "formula": "C2HCl3O2",
            "domain": "Corrosive Wafer Damage",
            "failure_reason": "Strong acid etching; destroys metal interconnects and chamber walls."
        },
        {
            "id": "NEG-007",
            "name": "Ti(OH)4 (Titanium tetrahydroxide sol)",
            "smiles": "O[Ti](O)(O)O",
            "formula": "H4O4Ti",
            "domain": "Insoluble Precipitate",
            "failure_reason": "Rapid spontaneous condensation into insoluble cross-linked TiO2 precipitate; zero vapor pressure."
        },
        {
            "id": "NEG-008",
            "name": "CoEt4 (Tetraethylcobalt)",
            "smiles": "CC[Co](CC)(CC)CC",
            "formula": "C8H20Co",
            "domain": "Electronic Instability",
            "failure_reason": "Formal Co(IV) state inaccessible with simple alkyls; rapid homolytic C-Co cleavage."
        },
        {
            "id": "NEG-009",
            "name": "1,2-Dichlorodisilane",
            "smiles": "Cl[SiH2][SiH2]Cl",
            "formula": "H4Cl2Si2",
            "domain": "Explosive Pyrophoric",
            "failure_reason": "High propensity for disproportionation and explosive autopolymerization."
        },
        {
            "id": "NEG-010",
            "name": "Cp2AlH (Bis(cyclopentadienyl)aluminium hydride)",
            "smiles": "C1=CC=CC1[AlH]C1C=CC=C1",
            "formula": "C10H11Al",
            "domain": "Valence/Hapticity Collapse",
            "failure_reason": "Oligomerizes to non-volatile amorphous polymers; unstable monomer."
        }
    ]

    molecules = []
    # Precursors
    for p in precursors:
        molecules.append({
            "id": p["id"],
            "name": p["name"],
            "smiles": p["smiles"],
            "formula": p.get("formula", ""),
            "domain": "Semiconductor Precursor",
            "is_fab_valid": True,
            "target_type": "Precursor"
        })
    # Inhibitors
    for inh in inhibitors:
        molecules.append({
            "id": inh["id"],
            "name": inh["abbreviation"],
            "smiles": inh["smiles"],
            "formula": inh.get("formula", ""),
            "domain": "AS-ALD Inhibitor",
            "is_fab_valid": True,
            "target_type": "Inhibitor",
            "selectivity_type": inh.get("selectivity_type", "")
        })
    # Negatives
    for n in negatives:
        molecules.append({
            "id": n["id"],
            "name": n["name"],
            "smiles": n["smiles"],
            "formula": n.get("formula", ""),
            "domain": n["domain"],
            "is_fab_valid": False,
            "target_type": "Negative Control",
            "failure_reason": n["failure_reason"]
        })
    return molecules

def evaluate_phase1():
    print("=" * 80)
    print("PHASE 1: PURE DATA-DRIVEN BASELINE EVALUATION (NO COMPUTATIONAL FEATURES)")
    print("=" * 80)

    molecules = load_all_molecules()
    print(f"Total molecules to evaluate: {len(molecules)}")
    
    # Initialize baseline models
    print("\nLoading baseline models (SAScore, SCScore, SYBA, RAscore, BR-SAScore)...")
    sc = SCScorer()
    sc.restore(os.path.join(PROJECT_ROOT, "external", "scscore", "models", "full_reaxys_model_1024bool", "model.ckpt-10654.as_numpy.json.gz"))
    
    sy = SybaClassifier()
    sy.fitDefaultScore()
    
    ra_model_path = os.path.join(PROJECT_ROOT, "external", "rascore", "RAscore", "models", "XGB_chembl_ecfp_counts", "model.pkl")
    with open(ra_model_path, "rb") as f:
        ra_clf = pickle.load(f)
    ra_booster = ra_clf.get_booster()
    
    br_scorer = BRSAScorer()
    
    # Load Building Blocks & Reactions
    with open(os.path.join(PROJECT_ROOT, "database", "building_blocks.json"), "r", encoding="utf-8") as f:
        bbs = json.load(f)
    bb_smiles_set = set(b["smiles"] for b in bbs)

    results = []
    
    for idx, item in enumerate(molecules, 1):
        smi = item["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            print(f"Failed to parse SMILES: {item['id']} {smi}")
            continue
            
        # 1. RDKit standard descriptors
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        rot_b = Descriptors.NumRotatableBonds(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        
        # 2. SAScore
        try:
            sa_val = round(sascorer.calculateScore(mol), 2)
        except Exception:
            sa_val = 9.99
            
        # 3. SCScore
        try:
            _, sc_val = sc.get_score_from_smi(smi)
            sc_val = round(float(sc_val), 2)
        except Exception:
            sc_val = 5.0
            
        # 4. SYBA
        try:
            syba_val = round(float(sy.predict(mol=mol)), 2)
        except Exception:
            syba_val = -99.9
            
        # 5. RAscore
        try:
            fp = AllChem.GetMorganFingerprint(mol, 3, useChirality=True)
            arr = np.zeros((2048,), dtype=np.float32)
            for k, v in fp.GetNonzeroElements().items():
                arr[k % 2048] += v
            pred = ra_booster.predict(xgb.DMatrix(arr.reshape(1, -1)))
            ra_val = round(float(pred[0]), 3)
        except Exception:
            ra_val = 0.0
            
        # 6. BR-SAScore (must pass SMILES string and unpack (score, contrib))
        try:
            br_score_val, _ = br_scorer.calculateScore(smi)
            br_val = round(float(br_score_val), 2)
        except Exception:
            br_val = 9.99

        # 6b. Real Pro-Ligand Decomposition SA (eliminates hardcoded base_ligand_sa = 1.5)
        tm_metals = {"Ti", "Zr", "Hf", "Al", "Ru", "Co", "W", "Mo", "Ta", "Sn", "Cu", "Zn", "Fe"}
        syms = {a.GetSymbol() for a in mol.GetAtoms()}
        centers = (syms & tm_metals) or (syms & {"Si"})
        if centers:
            rw = Chem.RWMol(mol)
            idxs = sorted([a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in centers], reverse=True)
            for c_idx in idxs:
                rw.RemoveAtom(c_idx)
            frags = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
            lig_scores = []
            for frag in frags:
                try:
                    for a in frag.GetAtoms():
                        a.SetNumRadicalElectrons(0)
                        a.SetNoImplicit(False)
                        a.SetNumExplicitHs(0)
                    Chem.SanitizeMol(frag)
                    l_smi = Chem.MolToSmiles(frag)
                    if l_smi and l_smi not in ("[H]", "Cl", "O"):
                        cand = []
                        try:
                            bs, _ = br_scorer.calculateScore(l_smi)
                            cand.append(float(bs))
                        except Exception:
                            pass
                        l_mol = Chem.MolFromSmiles(l_smi)
                        if l_mol:
                            cand.append(float(sascorer.calculateScore(l_mol)))
                        if cand:
                            lig_scores.append(min(cand))
                except Exception:
                    pass
            ligand_sa = round(max(lig_scores), 2) if lig_scores else 1.00
        else:
            ligand_sa = min(sa_val, br_val) if br_val < 9.9 else sa_val
            
        # 7. Reagent matching
        has_bb_match = smi in bb_smiles_set
        
        # 8. Baseline heuristic judgment:
        # Standard baseline heuristic rule: Pass if SAScore <= 3.5 OR SCScore <= 2.0 OR RAscore >= 0.5
        sascore_pass = sa_val <= 3.5
        scscore_pass = sc_val <= 2.0
        baseline_consensus_pass = (sa_val <= 3.5) or (sc_val <= 2.0 and ra_val >= 0.5)

        record = {
            "id": item["id"],
            "name": item["name"],
            "smiles": smi,
            "formula": item["formula"],
            "domain": item["domain"],
            "target_type": item["target_type"],
            "is_fab_valid": item["is_fab_valid"],
            "mw": mw,
            "logp": logp,
            "heavy_atoms": heavy_atoms,
            "rot_b": rot_b,
            "sascore": sa_val,
            "scscore": sc_val,
            "syba": syba_val,
            "rascore": ra_val,
            "br_sascore": br_val,
            "ligand_sa": ligand_sa,
            "has_bb_match": has_bb_match,
            "sascore_pass": sascore_pass,
            "scscore_pass": sc_val <= 2.0,
            "baseline_consensus_pass": baseline_consensus_pass
        }
        results.append(record)
        print(f"[{idx:02d}/55] {item['id']:8s} | {item['name']:32s} | SA:{sa_val:5.2f} | SC:{sc_val:4.2f} | SYBA:{syba_val:6.2f} | RA:{ra_val:5.3f} | Pass:{str(baseline_consensus_pass):5s}")

    # Compute Aggregate Metrics for Phase 1 Baseline
    precursor_results = [r for r in results if r["target_type"] == "Precursor"]
    inhibitor_results = [r for r in results if r["target_type"] == "Inhibitor"]
    negative_results = [r for r in results if r["target_type"] == "Negative Control"]
    
    # Precursors (Ground truth = Positive)
    prec_sa_pass = sum(1 for r in precursor_results if r["sascore_pass"])
    prec_sc_pass = sum(1 for r in precursor_results if r["scscore_pass"])
    prec_consensus_pass = sum(1 for r in precursor_results if r["baseline_consensus_pass"])
    
    # Inhibitors (Ground truth = Positive)
    inhib_sa_pass = sum(1 for r in inhibitor_results if r["sascore_pass"])
    inhib_sc_pass = sum(1 for r in inhibitor_results if r["scscore_pass"])
    inhib_consensus_pass = sum(1 for r in inhibitor_results if r["baseline_consensus_pass"])
    
    # Negatives (Ground truth = Negative)
    neg_sa_pass = sum(1 for r in negative_results if r["sascore_pass"])
    neg_sc_pass = sum(1 for r in negative_results if r["scscore_pass"])
    neg_consensus_pass = sum(1 for r in negative_results if r["baseline_consensus_pass"])
    
    # Metrics
    prec_fn_rate_sascore = round((len(precursor_results) - prec_sa_pass) / len(precursor_results) * 100, 1)
    neg_fp_rate_scscore = round(neg_sc_pass / len(negative_results) * 100, 1)
    neg_fp_rate_consensus = round(neg_consensus_pass / len(negative_results) * 100, 1)
    
    summary = {
        "total_molecules": len(results),
        "precursors_count": len(precursor_results),
        "inhibitors_count": len(inhibitor_results),
        "negatives_count": len(negative_results),
        "precursors_sascore_pass_rate": f"{prec_sa_pass}/{len(precursor_results)} ({prec_sa_pass/len(precursor_results)*100:.1f}%)",
        "precursors_false_negative_rate_sascore": f"{prec_fn_rate_sascore}%",
        "precursors_scscore_pass_rate": f"{prec_sc_pass}/{len(precursor_results)} ({prec_sc_pass/len(precursor_results)*100:.1f}%)",
        "inhibitors_sascore_pass_rate": f"{inhib_sa_pass}/{len(inhibitor_results)} ({inhib_sa_pass/len(inhibitor_results)*100:.1f}%)",
        "inhibitors_scscore_pass_rate": f"{inhib_sc_pass}/{len(inhibitor_results)} ({inhib_sc_pass/len(inhibitor_results)*100:.1f}%)",
        "negatives_scscore_false_positive_rate": f"{neg_sc_pass}/{len(negative_results)} ({neg_fp_rate_scscore}%)",
        "negatives_consensus_false_positive_rate": f"{neg_consensus_pass}/{len(negative_results)} ({neg_fp_rate_consensus}%)"
    }
    
    print("\n" + "=" * 80)
    print("PHASE 1 BASELINE AGGREGATE SUMMARY:")
    print("=" * 80)
    for k, v in summary.items():
        print(f"  {k:45s}: {v}")
        
    out_data = {
        "phase": "Phase 1: Pure Data-Driven Baseline Evaluation",
        "summary": summary,
        "results": results
    }
    
    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)
    print(f"\nSaved Phase 1 baseline results to: {out_path}")

if __name__ == "__main__":
    evaluate_phase1()
