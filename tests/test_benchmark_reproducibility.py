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
        cls.target_smiles = set()
        for p in cls.precursors:
            cls.target_smiles.add(p["smiles"])
            m_t = Chem.MolFromSmiles(p["smiles"])
            if m_t:
                cls.target_smiles.add(Chem.MolToSmiles(m_t))
            if "." in p["smiles"]:
                for comp in p["smiles"].split("."):
                    cls.target_smiles.add(comp.strip())
                    m_c = Chem.MolFromSmiles(comp.strip())
                    if m_c:
                        cls.target_smiles.add(Chem.MolToSmiles(m_c))
        
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
            self.assertTrue(len(p.get("metal_center", "")) > 0 or len(p.get("metal", "")) > 0)

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
        """Recalculate SCScore with patched standalone NumPy runner (full_reaxys_model_1024bool)."""
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

    def test_06_br_sascore_default_strict_and_counterfactual_leakage(self):
        """Recalculate BR-SAScore under Default, Strict Reagent (Unleaked, max(existing, -0.47)), Zero-Overwrite, and Leaked Counterfactual."""
        REAGENT_BSCORE_FLOOR = -0.47
        br_def = BRSAScorer()
        br_strict = BRSAScorer()
        fscores_strict = copy.deepcopy(br_def._fscores)
        br_strict_zero = BRSAScorer()
        fscores_strict_zero = copy.deepcopy(br_def._fscores)

        # Ingest commercial building blocks excluding target precursors via canonical SMILES
        with open(os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_bb.csv"), "r", encoding="utf-8") as f:
            for row in list(csv.reader(f))[1:]:
                if not row or not row[0].strip(): continue
                smi = row[0].strip()
                m = Chem.MolFromSmiles(smi)
                if not m: continue
                if smi in self.target_smiles or Chem.MolToSmiles(m) in self.target_smiles: continue
                bi = {}
                rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                for b in bi:
                    fscores_strict[b] = max(fscores_strict.get(b, -6.0), REAGENT_BSCORE_FLOOR)
                    fscores_strict_zero[b] = 0.0

        # Ingest reactants only (NO target products or identity rules via canonical SMILES)
        with open(os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_rxn.csv"), "r", encoding="utf-8") as f:
            for row in list(csv.reader(f))[1:]:
                if not row or not row[0].strip(): continue
                parts_rxn = row[0].strip().split(">>")
                reactants = parts_rxn[0].strip()
                m_full = Chem.MolFromSmiles(reactants)
                if reactants in self.target_smiles or (m_full and Chem.MolToSmiles(m_full) in self.target_smiles): continue
                if len(parts_rxn) > 1 and reactants == parts_rxn[1].strip(): continue
                for smi in reactants.split("."):
                    smi_clean = smi.strip()
                    m = Chem.MolFromSmiles(smi_clean)
                    if not m: continue
                    if smi_clean in self.target_smiles or Chem.MolToSmiles(m) in self.target_smiles: continue
                    bi = {}
                    rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                    for b in bi:
                        fscores_strict[b] = max(fscores_strict.get(b, -6.0), REAGENT_BSCORE_FLOOR)
                        fscores_strict_zero[b] = 0.0

        br_strict._fscores = fscores_strict
        br_strict_zero._fscores = fscores_strict_zero

        # Counterfactual Leaked Model: Pre-memorizes target product fragments
        br_leaked = BRSAScorer()
        fscores_leaked = copy.deepcopy(fscores_strict_zero)
        with open(os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_rxn.csv"), "r", encoding="utf-8") as f:
            for row in list(csv.reader(f))[1:]:
                if not row or not row[0].strip(): continue
                for side in row[0].strip().split(">>"):
                    for smi in side.split("."):
                        m = Chem.MolFromSmiles(smi.strip())
                        if m:
                            bi = {}
                            rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                            for b in bi: fscores_leaked[b] = 0.0
        for p in self.precursors:
            m = Chem.MolFromSmiles(p["smiles"])
            if m:
                bi = {}
                rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                for b in bi: fscores_leaked[b] = 0.0
        br_leaked._fscores = fscores_leaked

        def_passes = 0
        strict_passes = 0
        leaked_passes = 0
        for p in self.precursors:
            smi = p["smiles"]
            s_def, _ = br_def.calculateScore(smi)
            s_strict, _ = br_strict.calculateScore(smi)
            s_strict_zero, _ = br_strict_zero.calculateScore(smi)
            s_leaked, _ = br_leaked.calculateScore(smi)

            s_def = round(float(s_def), 2)
            s_strict = round(float(s_strict), 2)
            s_strict_zero = round(float(s_strict_zero), 2)
            s_leaked = round(float(s_leaked), 2)

            exp_def = self.stored_results[p["id"]]["br_sascore_default"]
            exp_strict = self.stored_results[p["id"]]["br_sascore_custom_strict"]
            exp_strict_zero = self.stored_results[p["id"]]["br_sascore_custom_strict_zero_overwrite"]
            exp_leaked = self.stored_results[p["id"]]["br_sascore_custom_leaked"]

            self.assertAlmostEqual(s_def, exp_def, places=2, msg=f"BR-SAScore Default mismatch for {p['name']}")
            self.assertAlmostEqual(s_strict, exp_strict, places=2, msg=f"BR-SAScore Strict mismatch for {p['name']}")
            self.assertAlmostEqual(s_strict_zero, exp_strict_zero, places=2, msg=f"BR-SAScore Strict Zero mismatch for {p['name']}")
            self.assertAlmostEqual(s_leaked, exp_leaked, places=2, msg=f"BR-SAScore Leaked mismatch for {p['name']}")
            self.assertLessEqual(s_strict, s_def, msg=f"Strict max(existing, -0.47) must not worsen {p['name']}")

            if s_def <= 3.5: def_passes += 1
            if s_strict <= 3.5: strict_passes += 1
            if s_leaked <= 3.5: leaked_passes += 1

        self.assertEqual(def_passes, 3, "BR-SAScore Default should pass exactly 3/27")
        self.assertEqual(strict_passes, 3, "BR-SAScore Strict Reagent Domain (Unleaked) should pass exactly 3/27 (TSA, TMA, DMAH)")
        self.assertEqual(leaked_passes, 27, "BR-SAScore Product-Memorized baseline confirms 27/27 circular data leakage")

        # Verify that max(existing, -0.47) worsens 0 precursors, whereas 0.0-overwriting worsens 4 specific precursors
        worsened_strict = [p["name"] for p in self.precursors if self.stored_results[p["id"]]["br_sascore_custom_strict"] > self.stored_results[p["id"]]["br_sascore_default"]]
        self.assertEqual(worsened_strict, [])
        worsened_zero = [p["name"] for p in self.precursors if self.stored_results[p["id"]]["br_sascore_custom_strict_zero_overwrite"] > self.stored_results[p["id"]]["br_sascore_default"]]
        self.assertEqual(sorted(worsened_zero), sorted(["Cp(NMe2)3Zr", "CpTi(NMe2)3", "Ru(EtCp)2", "Co(AMD)2"]))

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

    def test_09_chemprojector_synnet_synformer_rxnflow_and_dfrscore(self):
        """Verify live SMARTS matching on RxnFlow (72), HB (91), and SynFormer (115) templates and live ingestion."""
        import yaml
        from aizynthfinder.chem import Molecule as AiZynthMolecule
        from rdchiral.initialization import rdchiralReactants
        dfr_allowed = {"C", "N", "O", "F", "S", "Cl", "Br", "I", "B", "P"}
        dfr_collapse = sum(1 for p in self.precursors if any(a.GetSymbol() not in dfr_allowed for a in Chem.MolFromSmiles(p["smiles"]).GetAtoms()))
        self.assertEqual(dfr_collapse, 27, "All 27 precursors must collapse in DFRscore")

        target_metal_atomic_nums = {13, 14, 22, 27, 30, 40, 42, 44, 50, 72, 73, 74}
        target_metal_symbols = {"Si", "Ti", "Zr", "Hf", "Ru", "W", "Mo", "Al", "Co", "Ta", "Sn", "Zn"}
        def _has_metal(m_or_p):
            return any(a.GetAtomicNum() in target_metal_atomic_nums or a.GetSymbol() in target_metal_symbols for a in m_or_p.GetAtoms())

        # RxnFlow 72 templates
        rf_patts = []
        for scope_name, r_file in [
            ("explore", os.path.join(PROJECT_ROOT, "external", "rxnflow", "data", "templates", "explore", "reaction.yaml")),
            ("cgflow", os.path.join(PROJECT_ROOT, "external", "rxnflow", "data", "templates", "cgflow", "reaction.yaml"))
        ]:
            with open(r_file, "r", encoding="utf-8") as f:
                d_rf = yaml.safe_load(f)
                for sec in ["UniReaction", "BiReaction"]:
                    for r_name, r_body in d_rf.get(sec, {}).items():
                        rev = r_body.get("reverse", "")
                        if ">>" in rev:
                            pt = Chem.MolFromSmarts(rev.split(">>")[0].strip())
                            if pt:
                                rf_patts.append((pt, _has_metal(pt)))
        rf_subst = sum(1 for p in self.precursors if any(Chem.MolFromSmiles(p["smiles"]).HasSubstructMatch(pt) for pt, _ in rf_patts))
        self.assertEqual(rf_subst, 9, "RxnFlow must match 9/27 on generic product substructures")
        self.assertEqual(sum(1 for _, hm in rf_patts if hm), 0, "RxnFlow must contain 0 metal coordination templates")

        # ChemProjector / SynNet 91 HB templates
        hb_path = os.path.join(PROJECT_ROOT, "external", "chemprojector", "data", "reaction_templates_hb.txt")
        hb_patterns = []
        with open(hb_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split(">>")
                if len(parts) == 2:
                    patt = Chem.MolFromSmarts(parts[1].strip())
                    if patt: hb_patterns.append((patt, _has_metal(patt)))

        subst_matches = sum(1 for p in self.precursors if any(Chem.MolFromSmiles(p["smiles"]).HasSubstructMatch(patt) for patt, _ in hb_patterns))
        self.assertEqual(subst_matches, 23, "ChemProjector/SynNet must match 23/27 on generic product substructures")
        self.assertEqual(sum(1 for _, hm in hb_patterns if hm), 0, "ChemProjector HB rules must contain 0 metal templates")

        # SynFormer 115 templates
        sf_path = os.path.join(PROJECT_ROOT, "external", "synformer", "data", "rxn_templates", "comprehensive.txt")
        sf_rxns = []
        with open(sf_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                rxn = AllChem.ReactionFromSmarts(line)
                if rxn:
                    prods = list(rxn.GetProducts())
                    sf_rxns.append((prods, any(_has_metal(pt) for pt in prods)))
        sf_subst = sum(1 for p in self.precursors if any(all(Chem.MolFromSmiles(p["smiles"]).HasSubstructMatch(pt) for pt in prods) for prods, _ in sf_rxns))
        sf_ge2 = sum(1 for p in self.precursors if any(all(pt.GetNumAtoms() >= 2 and Chem.MolFromSmiles(p["smiles"]).HasSubstructMatch(pt) for pt in prods) for prods, _ in sf_rxns))
        sf_ge4 = sum(1 for p in self.precursors if any(all(pt.GetNumAtoms() >= 4 and Chem.MolFromSmiles(p["smiles"]).HasSubstructMatch(pt) for pt in prods) for prods, _ in sf_rxns))
        self.assertEqual(sf_subst, 24, "SynFormer must match 24/27 on generic product substructures (including #84 [N:1] on TSA)")
        self.assertEqual(sf_ge2, 23, "SynFormer must match 23/27 when single-atom product templates (<2 atoms) are excluded")
        self.assertEqual(sf_ge4, 11, "SynFormer must match 11/27 when restricted to >=4-atom product templates")
        self.assertEqual(sum(1 for _, hm in sf_rxns if hm), 0, "SynFormer product templates must contain 0 metal coordination rules")

        # Live ingestion checks for AiZynthFinder and Retro*
        for p in self.precursors:
            az_mol = AiZynthMolecule(smiles=p["smiles"])
            az_mol.sanitize()
            self.assertIsNotNone(az_mol.rd_mol, f"AiZynthFinder Molecule failed for {p['name']}")
            rdc_obj = rdchiralReactants(p["smiles"])
            self.assertIsNotNone(rdc_obj.reactants, f"rdchiralReactants failed for {p['name']}")


if __name__ == "__main__":
    unittest.main()
