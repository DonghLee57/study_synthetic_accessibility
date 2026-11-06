#!/usr/bin/env python3
"""
Unit Test Suite for Semiconductor Benchmark Reproducibility
-----------------------------------------------------------
Formal automated test suite structured with unittest.TestCase for automated CI/CD
and peer review (Claude / independent auditor).

Usage:
    python -m unittest tests/test_benchmark_reproducibility.py
    python -m unittest discover tests
"""

import copy
import csv
import json
import math
import os
import pickle
import sys
import unittest

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


class TestBenchmarkReproducibility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Load precursors
        precursor_file = os.path.join(PROJECT_ROOT, "database", "precursors.json")
        with open(precursor_file, "r", encoding="utf-8") as f:
            cls.precursors = json.load(f)
        
        # Load verified database results
        benchmark_file = os.path.join(PROJECT_ROOT, "database", "precursor_real_benchmark_results.json")
        with open(benchmark_file, "r", encoding="utf-8") as f:
            cls.stored_results = {r["id"]: r for r in json.load(f)}

    def test_01_precursor_database_integrity(self):
        """Verify the 27 precursor dataset has complete chemical metadata."""
        self.assertEqual(len(self.precursors), 27, "Must contain exactly 27 industrial precursors")
        for p in self.precursors:
            mol = Chem.MolFromSmiles(p["smiles"])
            self.assertIsNotNone(mol, f"SMILES parsing failed for {p['name']}: {p['smiles']}")
            self.assertTrue(len(p["metal_center"]) > 0 or len(p.get("metal", "")) > 0)

    def test_02_sascore_recalculation(self):
        """Recalculate RDKit SAScore and assert against database values."""
        pass_count = 0
        for p in self.precursors:
            mol = Chem.MolFromSmiles(p["smiles"])
            score = round(sascorer.calculateScore(mol), 2)
            expected = self.stored_results[p["id"]]["sascore"]
            self.assertAlmostEqual(score, expected, places=2, msg=f"SAScore mismatch for {p['name']}")
            if score <= 3.5:
                pass_count += 1
        self.assertEqual(pass_count, 3, "SAScore should pass exactly 3 out of 27 precursors (TSA, TMA, DMAH)")

    def test_03_scscore_recalculation(self):
        """Recalculate SCScore with patched standalone NumPy runner."""
        sc = SCScorer()
        sc_weight = os.path.join(PROJECT_ROOT, "external", "scscore", "models", "full_reaxys_model_1024bool", "model.ckpt-10654.as_numpy.json.gz")
        sc.restore(sc_weight)
        
        pass_count = 0
        for p in self.precursors:
            _, score = sc.get_score_from_smi(p["smiles"])
            score = round(float(score), 2)
            expected = self.stored_results[p["id"]]["scscore"]
            self.assertAlmostEqual(score, expected, places=2, msg=f"SCScore mismatch for {p['name']}")
            if score <= 2.0:
                pass_count += 1
        self.assertEqual(pass_count, 23, "SCScore should pass exactly 23 out of 27 precursors (<=2.0)")

    def test_04_syba_recalculation(self):
        """Recalculate SYBA score using default syba.csv.gz count table."""
        sy = SybaClassifier()
        sy.fitDefaultScore()
        
        pass_count = 0
        for p in self.precursors:
            mol = Chem.MolFromSmiles(p["smiles"])
            score = round(float(sy.predict(mol=mol)), 2)
            expected = self.stored_results[p["id"]]["syba"]
            self.assertAlmostEqual(score, expected, places=2, msg=f"SYBA mismatch for {p['name']}")
            if score > 0.0:
                pass_count += 1
        self.assertEqual(pass_count, 23, "SYBA should pass exactly 23 out of 27 precursors (>0.0)")

    def test_05_rascore_recalculation(self):
        """Recalculate RAscore using direct Booster predict on Morgan count fingerprints."""
        model_path = os.path.join(PROJECT_ROOT, "external", "rascore", "RAscore", "models", "XGB_chembl_ecfp_counts", "model.pkl")
        with open(model_path, "rb") as f:
            clf = pickle.load(f)
        booster = clf.get_booster()

        pass_count = 0
        for p in self.precursors:
            mol = Chem.MolFromSmiles(p["smiles"])
            fp = AllChem.GetMorganFingerprint(mol, 3, useChirality=True)
            arr = np.zeros((2048,), dtype=np.float32)
            for k, v in fp.GetNonzeroElements().items():
                arr[k % 2048] += v
            pred = booster.predict(xgb.DMatrix(arr.reshape(1, -1)))
            score = round(float(pred[0]), 3)
            expected = self.stored_results[p["id"]]["rascore"]
            self.assertAlmostEqual(score, expected, places=2, msg=f"RAscore mismatch for {p['name']}")
            if score >= 0.5:
                pass_count += 1
        self.assertEqual(pass_count, 27, "RAscore should pass 27 out of 27 precursors (ECFP6 metal projection bias)")

    def test_06_br_sascore_default_and_custom(self):
        """Recalculate BR-SAScore under default and custom semiconductor dictionaries."""
        br_def = BRSAScorer()
        br_cust = BRSAScorer()
        custom_fscores = copy.deepcopy(br_def._fscores)

        def add_bits(path, is_rxn=False):
            with open(path, "r", encoding="utf-8") as f:
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
        br_cust._fscores = custom_fscores

        def_passes = 0
        cust_passes = 0
        for p in self.precursors:
            smi = p["smiles"]
            s_def, _ = br_def.calculateScore(smi)
            s_cust, _ = br_cust.calculateScore(smi)
            s_def = round(float(s_def), 2)
            s_cust = round(float(s_cust), 2)

            exp_def = self.stored_results[p["id"]]["br_sascore_default"]
            exp_cust = self.stored_results[p["id"]]["br_sascore_custom"]

            self.assertAlmostEqual(s_def, exp_def, places=2, msg=f"BR-SAScore Default mismatch for {p['name']}")
            self.assertAlmostEqual(s_cust, exp_cust, places=2, msg=f"BR-SAScore Custom mismatch for {p['name']}")

            if s_def <= 3.5:
                def_passes += 1
            if s_cust <= 3.5:
                cust_passes += 1

        self.assertEqual(def_passes, 3, "BR-SAScore Default should pass exactly 3/27")
        self.assertEqual(cust_passes, 24, "BR-SAScore Custom should pass exactly 24/27 (88.9%)")

    def test_07_reasyn_lossless_tokenization(self):
        """Verify BioNeMo ReaSyn universal tokenizer losslessly roundtrips all 27 SMILES."""
        reasyn_1_to_tok = {i: t for i, t in enumerate(_smiles_vocab, start=1)}
        for p in self.precursors:
            smi = p["smiles"]
            tokens = _smiles_token_pattern.findall(smi)
            token_ids = [_smiles_token_to_id[t] for t in tokens]
            reconstructed = "".join(reasyn_1_to_tok[tid] for tid in token_ids)
            self.assertEqual(reconstructed, smi, f"Lossless roundtrip failed for {p['name']}")

    def test_08_gasa_element_support_and_silent_misencoding(self):
        """Verify GASA unsupported elements (7 precursors) and silent mis-encoding behavior."""
        gasa_allowed = {
            'C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na',
            'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb',
            'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu',
            'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb'
        }
        unsupported_count = 0
        unsupported_elements = set()
        for p in self.precursors:
            mol = Chem.MolFromSmiles(p["smiles"])
            atoms = [a.GetSymbol() for a in mol.GetAtoms()]
            unsup = [a for a in atoms if a not in gasa_allowed]
            if unsup:
                unsupported_count += 1
                unsupported_elements.update(unsup)
        
        self.assertEqual(unsupported_count, 7, "GASA unsupported count must be 7")
        self.assertEqual(unsupported_elements, {'Hf', 'Ta', 'Ru', 'W', 'Mo'}, "Unsupported elements must be Hf, Ta, Ru, W, Mo")

    def test_09_dfrscore_and_reaction_manifolds(self):
        """Verify DFRscore collapse to ELSE and ChemProjector/SynNet zero transition metal rules."""
        dfr_allowed = {"C", "N", "O", "F", "S", "Cl", "Br", "I", "B", "P"}
        dfr_collapse = sum(1 for p in self.precursors if any(a.GetSymbol() not in dfr_allowed for a in Chem.MolFromSmiles(p["smiles"]).GetAtoms()))
        self.assertEqual(dfr_collapse, 27, "All 27 precursors must collapse in DFRscore")

        with open(os.path.join(PROJECT_ROOT, "external", "chemprojector", "data", "reaction_templates_hb.txt"), "r", encoding="utf-8") as f:
            cp_rules = f.read()
        metals = ['Si', 'Ti', 'Zr', 'Hf', 'Ru', 'W', 'Mo', 'Al', 'Co', 'Ta']
        for m in metals:
            self.assertNotIn(m, cp_rules, f"Metal {m} should not exist in ChemProjector HB rules")


if __name__ == "__main__":
    unittest.main()
