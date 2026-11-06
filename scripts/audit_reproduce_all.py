#!/usr/bin/env python3
"""
Reproducible Benchmark Verification & Audit Script
--------------------------------------------------
Recalculates all model scores and verifies empirical benchmark metrics
across all 27 semiconductor precursors.
"""

import copy
import csv
import json
import math
import os
import pickle
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig, rdMolDescriptors
import xgboost as xgb

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "scscore"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "reasyn"))

import sascorer
from scscore.standalone_model_numpy import SCScorer
from syba.syba import SybaClassifier
from BRSAScore import SAScorer as BRSAScorer
from reasyn.chem.featurize import _smiles_vocab, _smiles_token_pattern, _smiles_token_to_id


def run_full_audit():
    print("=" * 80)
    print("AUDIT VERIFICATION SUITE: EMPIRICAL RECALCULATION ACROSS ALL 27 PRECURSORS")
    print("=" * 80)

    # Load Database
    with open(os.path.join(PROJECT_ROOT, "database", "precursors.json"), "r", encoding="utf-8") as f:
        precursors = json.load(f)
    assert len(precursors) == 27, f"Expected 27 precursors, found {len(precursors)}"

    with open(os.path.join(PROJECT_ROOT, "database", "precursor_real_benchmark_results.json"), "r", encoding="utf-8") as f:
        stored_results = {r["id"]: r for r in json.load(f)}

    # 1. Recalculate SAScore (RDKit)
    print("\n[CHECK 1/8] Recalculating RDKit SAScore...")
    sa_pass = 0
    for p in precursors:
        mol = Chem.MolFromSmiles(p["smiles"])
        score = round(sascorer.calculateScore(mol), 2)
        expected = stored_results[p["id"]]["sascore"]
        assert abs(score - expected) < 1e-2, f"SAScore mismatch for {p['name']}: {score} vs {expected}"
        if score <= 3.5:
            sa_pass += 1
    print(f"  -> Recalculated SAScore for all 27 precursors. Pass (<=3.5): {sa_pass}/27 (Expected: 3/27, 11.1%)")
    assert sa_pass == 3, f"Expected 3 passes, got {sa_pass}"

    # 2. Recalculate SCScore (NumPy Standalone)
    print("\n[CHECK 2/8] Recalculating SCScore (Reaxys NumPy Standalone)...")
    sc = SCScorer()
    sc_weight_path = os.path.join(PROJECT_ROOT, "external", "scscore", "models", "full_reaxys_model_1024bool", "model.ckpt-10654.as_numpy.json.gz")
    sc.restore(sc_weight_path)
    sc_pass = 0
    for p in precursors:
        _, score = sc.get_score_from_smi(p["smiles"])
        score = round(float(score), 2)
        expected = stored_results[p["id"]]["scscore"]
        assert abs(score - expected) < 1e-2, f"SCScore mismatch for {p['name']}: {score} vs {expected}"
        if score <= 2.0:
            sc_pass += 1
    print(f"  -> Recalculated SCScore for all 27 precursors. Pass (<=2.0): {sc_pass}/27 (Expected: 23/27, 85.2%)")
    assert sc_pass == 23, f"Expected 23 passes, got {sc_pass}"

    # 3. Recalculate SYBA (Naive Bayesian)
    print("\n[CHECK 3/8] Recalculating SYBA (with syba.csv.gz)...")
    sy = SybaClassifier()
    sy.fitDefaultScore()
    sy_pass = 0
    for p in precursors:
        mol = Chem.MolFromSmiles(p["smiles"])
        score = round(float(sy.predict(mol=mol)), 2)
        expected = stored_results[p["id"]]["syba"]
        assert abs(score - expected) < 1e-2, f"SYBA mismatch for {p['name']}: {score} vs {expected}"
        if score > 0.0:
            sy_pass += 1
    print(f"  -> Recalculated SYBA for all 27 precursors. Pass (>0.0): {sy_pass}/27 (Expected: 23/27, 85.2%)")
    assert sy_pass == 23, f"Expected 23 passes, got {sy_pass}"

    # 4. Recalculate RAscore (XGBoost Booster)
    print("\n[CHECK 4/8] Recalculating RAscore (XGBoost Booster)...")
    ra_model_path = os.path.join(PROJECT_ROOT, "external", "rascore", "RAscore", "models", "XGB_chembl_ecfp_counts", "model.pkl")
    with open(ra_model_path, "rb") as f:
        ra_clf = pickle.load(f)
    ra_booster = ra_clf.get_booster()
    ra_pass = 0
    for p in precursors:
        mol = Chem.MolFromSmiles(p["smiles"])
        fp = AllChem.GetMorganFingerprint(mol, 3, useChirality=True)
        arr = np.zeros((2048,), dtype=np.float32)
        for k, v in fp.GetNonzeroElements().items():
            arr[k % 2048] += v
        pred = ra_booster.predict(xgb.DMatrix(arr.reshape(1, -1)))
        score = round(float(pred[0]), 3)
        expected = stored_results[p["id"]]["rascore"]
        assert abs(score - expected) < 1e-2, f"RAscore mismatch for {p['name']}: {score} vs {expected}"
        if score >= 0.5:
            ra_pass += 1
    print(f"  -> Recalculated RAscore for all 27 precursors. Pass (>=0.5): {ra_pass}/27 (Expected: 27/27, 100.0%)")
    assert ra_pass == 27, f"Expected 27 passes, got {ra_pass}"

    # 5. Recalculate BR-SAScore (Default vs Custom Semiconductor Domain)
    print("\n[CHECK 5/8] Recalculating BR-SAScore (Default & Custom Domain)...")
    br_default = BRSAScorer()
    br_custom = BRSAScorer()
    custom_fscores = copy.deepcopy(br_default._fscores)

    def add_bits(csv_path, is_rxn=False):
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if not row or not row[0].strip():
                    continue
                line = row[0].strip()
                if is_rxn:
                    for side in line.split(">>"):
                        for smi in side.split("."):
                            m = Chem.MolFromSmiles(smi.strip())
                            if m:
                                bi = {}
                                rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                                for b in bi:
                                    custom_fscores[b] = 0.0
                else:
                    m = Chem.MolFromSmiles(line)
                    if m:
                        bi = {}
                        rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                        for b in bi:
                            custom_fscores[b] = 0.0

    add_bits(os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_bb.csv"), is_rxn=False)
    add_bits(os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_rxn.csv"), is_rxn=True)
    br_custom._fscores = custom_fscores

    br_def_pass = 0
    br_cust_pass = 0
    for p in precursors:
        smi = p["smiles"]
        s_def, _ = br_default.calculateScore(smi)
        s_def = round(float(s_def), 2)
        s_cust, _ = br_custom.calculateScore(smi)
        s_cust = round(float(s_cust), 2)

        exp_def = stored_results[p["id"]]["br_sascore_default"]
        exp_cust = stored_results[p["id"]]["br_sascore_custom"]
        assert abs(s_def - exp_def) < 1e-2, f"BR-SAScore Default mismatch for {p['name']}: {s_def} vs {exp_def}"
        assert abs(s_cust - exp_cust) < 1e-2, f"BR-SAScore Custom mismatch for {p['name']}: {s_cust} vs {exp_cust}"

        if s_def <= 3.5:
            br_def_pass += 1
        if s_cust <= 3.5:
            br_cust_pass += 1

    print(f"  -> BR-SAScore Default Pass: {br_def_pass}/27 (Expected: 3/27, 11.1%)")
    print(f"  -> BR-SAScore Custom Pass:  {br_cust_pass}/27 (Expected: 24/27, 88.9%)")
    assert br_def_pass == 3, f"Expected 3 default passes, got {br_def_pass}"
    assert br_cust_pass == 24, f"Expected 24 custom passes, got {br_cust_pass}"

    # 6. Verify BioNeMo ReaSyn Tokenizer Roundtrip & Decoding Offset Bug Analysis
    print("\n[CHECK 6/8] Auditing BioNeMo ReaSyn Tokenizer & Decoding Offset...")
    reasyn_1_to_tok = {i: t for i, t in enumerate(_smiles_vocab, start=1)}
    reasyn_exact_count = 0
    for p in precursors:
        smi = p["smiles"]
        tokens = _smiles_token_pattern.findall(smi)
        token_ids = [_smiles_token_to_id[t] for t in tokens]
        reconstructed = "".join(reasyn_1_to_tok[tid] for tid in token_ids)
        assert reconstructed == smi, f"ReaSyn roundtrip mismatch for {p['name']}: {smi} vs {reconstructed}"
        reasyn_exact_count += 1
    print(f"  -> BioNeMo ReaSyn Tokenizer: 100% Exact Roundtrip ({reasyn_exact_count}/27) with 0 unknown tokens.")
    print("  -> Audit Note: Prior 'false' flag was caused by decoding offset in ReaSyn's decode_smiles()")
    print("     (_smiles_id_to_token starts at MOL_END+1=4 while _smiles_token_to_id starts at 1).")

    # 7. GASA Element Set & Silent Mis-encoding Verification
    print("\n[CHECK 7/8] Auditing GASA Element Support & Silent Mis-encoding...")
    gasa_allowed = {
        'C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na',
        'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb',
        'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu',
        'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb'
    }
    unsupported_count = 0
    unsupported_elements_found = set()
    for p in precursors:
        mol = Chem.MolFromSmiles(p["smiles"])
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        unsupported = [a for a in atoms if a not in gasa_allowed]
        if unsupported:
            unsupported_count += 1
            unsupported_elements_found.update(unsupported)
    print(f"  -> GASA Unsupported Precursors: {unsupported_count}/27 (Expected: 7/27 on {sorted(list(unsupported_elements_found))})")
    assert unsupported_count == 7, f"Expected 7 unsupported precursors, got {unsupported_count}"
    assert unsupported_elements_found == {'Hf', 'Ta', 'Ru', 'W', 'Mo'}, f"Unexpected element set: {unsupported_elements_found}"
    print("  -> Audit Note: dgllife.utils.atom_type_one_hot(..., encode_unknown=False) does not raise an unhandled")
    print("     exception; it emits an all-zero vector (Silent Mis-encoding / Feature Loss).")

    # 8. DFRscore & ChemProjector / SynNet Reaction Manifold Audit
    print("\n[CHECK 8/8] Auditing DFRscore Node Encoding & ChemProjector/SynNet Manifolds...")
    dfr_allowed = {"C", "N", "O", "F", "S", "Cl", "Br", "I", "B", "P"}
    dfr_collapse_count = sum(1 for p in precursors if any(a.GetSymbol() not in dfr_allowed for a in Chem.MolFromSmiles(p["smiles"]).GetAtoms()))
    assert dfr_collapse_count == 27, f"Expected 27 collapsed precursors in DFRscore, got {dfr_collapse_count}"
    print(f"  -> DFRscore Metal Center Collapse to 'ELSE': {dfr_collapse_count}/27 (100.0%)")

    with open(os.path.join(PROJECT_ROOT, "external", "chemprojector", "data", "reaction_templates_hb.txt"), "r", encoding="utf-8") as f:
        cp_rules = f.read()
    metals = ['Si', 'Ti', 'Zr', 'Hf', 'Ru', 'W', 'Mo', 'Al', 'Co', 'Ta']
    cp_metals = [m for m in metals if m in cp_rules]
    assert len(cp_metals) == 0, f"Unexpected metals found in ChemProjector rules: {cp_metals}"
    print("  -> ChemProjector HB Templates: 0 metal coordination rules found (100% organic drug templates).")

    print("\n" + "=" * 80)
    print(">>> ALL 8 EMPIRICAL RECALCULATION & DIAGNOSTIC CHECKS PASSED WITH 100% INTEGRITY! <<<")
    print("=" * 80)


if __name__ == "__main__":
    run_full_audit()
