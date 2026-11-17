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
import importlib.util
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
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
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
            self.assertEqual(
                rdMolDescriptors.CalcMolFormula(mol),
                p["formula"],
                f"Formula mismatch for {p['name']}: JSON={p['formula']} vs RDKit={rdMolDescriptors.CalcMolFormula(mol)}",
            )

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
        self.assertEqual(pass_count, 2, "SAScore should pass 2 out of 27 precursors (TEOS, 4DMAS) after connected metallocene SMILES fix")

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
        self.assertEqual(pass_count, 21, "SCScore should pass 21 out of 27 precursors (<=2.0) after Ru(EtCp)2 & CCTBA SMILES fix")

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

        # Verify that max(existing, -0.47) worsens 0 precursors, whereas 0.0-overwriting worsens 3 Cp/amidinate precursors
        worsened_strict = [p["name"] for p in self.precursors if self.stored_results[p["id"]]["br_sascore_custom_strict"] > self.stored_results[p["id"]]["br_sascore_default"]]
        self.assertEqual(worsened_strict, [])
        worsened_zero = [p["name"] for p in self.precursors if self.stored_results[p["id"]]["br_sascore_custom_strict_zero_overwrite"] > self.stored_results[p["id"]]["br_sascore_default"]]
        self.assertEqual(sorted(worsened_zero), sorted(["Cp(NMe2)3Zr", "CpTi(NMe2)3", "Co(AMD)2"]))

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

    @staticmethod
    def _load_script_module(rel_path: str, module_name: str):
        full_path = os.path.join(PROJECT_ROOT, rel_path)
        spec = importlib.util.spec_from_file_location(module_name, full_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load module from {full_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_10_phase1_55_molecule_baseline_reproducibility(self):
        """Verify 55-molecule Phase 1 baseline integrity, formula parity, and deterministic SAScore/BR-SAScore metrics."""
        p1_mod = self._load_script_module(
            os.path.join("scripts", "benchmark", "run_phase1_baseline.py"),
            "repo_run_phase1_baseline",
        )
        expected_mols = p1_mod.load_all_molecules()
        self.assertEqual(len(expected_mols), 55)

        p1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
        with open(p1_path, "r", encoding="utf-8") as f:
            p1_json = json.load(f)
        rows = p1_json["results"]
        self.assertEqual(len(rows), 55)

        by_id = {r["id"]: r for r in rows}
        self.assertEqual([m["id"] for m in expected_mols], [r["id"] for r in rows])

        br = BRSAScorer()
        for m in expected_mols:
            r = by_id[m["id"]]
            self.assertEqual(r["smiles"], m["smiles"])
            mol = Chem.MolFromSmiles(m["smiles"])
            self.assertIsNotNone(mol, f"Invalid SMILES in Phase 1: {m['id']}")
            self.assertEqual(
                rdMolDescriptors.CalcMolFormula(mol),
                m["formula"],
                f"Formula mismatch in 55-molecule dataset for {m['id']} ({m['name']})",
            )
            sa_calc = round(sascorer.calculateScore(mol), 2)
            self.assertAlmostEqual(sa_calc, r["sascore"], places=2)
            br_calc, _ = br.calculateScore(m["smiles"])
            self.assertAlmostEqual(round(float(br_calc), 2), r["br_sascore"], places=2)

    def test_11_sevenn_physical_features_integrity(self):
        """Verify sevenn_physical_features.json completeness, 0 fallbacks, relaxed invariants, and dimer gap >= 2.8 A."""
        p1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
        with open(p1_path, "r", encoding="utf-8") as f:
            p1_rows = json.load(f)["results"]

        p2_path = os.path.join(PROJECT_ROOT, "database", "external", "sevenn_physical_features.json")
        with open(p2_path, "r", encoding="utf-8") as f:
            p2_json = json.load(f)

        meta = p2_json["metadata"]
        feats = p2_json["features"]
        self.assertEqual(meta["total_molecules_computed"], 55)
        self.assertEqual(len(feats), 55)
        self.assertEqual(meta.get("fallback_count", -1), 0)
        self.assertEqual([r["id"] for r in p1_rows], [x["id"] for x in feats])

        for p1_item, feat in zip(p1_rows, feats):
            self.assertEqual(feat["smiles"], p1_item["smiles"])
            self.assertFalse(feat.get("is_fallback", True), f"Fallback triggered for {feat['id']}")
            self.assertNotEqual(feat["e_pot_ev"], -100.0)
            self.assertGreaterEqual(feat.get("mono_relax_steps", 0), 1)
            self.assertGreaterEqual(feat.get("dimer_init_gap_ang", 0.0), 2.8)
            self.assertTrue(0.0 <= feat["vbur_percent"] <= 100.0)
            self.assertGreater(feat["e_coh_ev_per_atom"], 0.0)
            for key in ("e_pot_ev", "e_coh_ev_per_atom", "bde_min_kcal_mol", "vbur_percent", "delta_e_dimer_kcal_mol", "radius_of_gyration_ang", "mol_volume_ang3"):
                self.assertTrue(math.isfinite(float(feat[key])), f"Non-finite {key} for {feat['id']}")

    def test_12_model_tuning_cv_loocv_reproducibility(self):
        """Verify deterministic reproduction of run_model_tuning_comparison.py analytical & CV/LOOCV outputs."""
        tuning_mod = self._load_script_module(
            os.path.join("scripts", "benchmark", "run_model_tuning_comparison.py"),
            "repo_run_model_tuning_comparison",
        )
        evaluate_machine_learning_models = tuning_mod.evaluate_machine_learning_models
        evaluate_unbiased_analytical_models = tuning_mod.evaluate_unbiased_analytical_models
        load_merged_data = tuning_mod.load_merged_data

        tuned_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json")
        with open(tuned_path, "r", encoding="utf-8") as f:
            stored = json.load(f)

        merged, _ = load_merged_data()
        ana_live = evaluate_unbiased_analytical_models(merged)
        ml_live = evaluate_machine_learning_models(merged)

        self.assertEqual(ana_live["overall_uniform_physics_only"], stored["analytical_comparison"]["overall_uniform_physics_only"])
        self.assertEqual(ana_live["overall_physics_plus_smarts_chem"], stored["analytical_comparison"]["overall_physics_plus_smarts_chem"])
        self.assertEqual(ml_live["full_dataset_n55"], stored["machine_learning_comparison"]["full_dataset_n55"])
        self.assertEqual(ml_live["precursor_subset_n35"], stored["machine_learning_comparison"]["precursor_subset_n35"])
        self.assertEqual(ml_live["all_sample_phys_importance_sum_pct"], stored["machine_learning_comparison"]["all_sample_phys_importance_sum_pct"])
        self.assertNotIn("delta_e_complex_kcal_mol", ml_live["all_sample_rf_feature_importances"])

    def test_13_expanded_target_elements_coverage_and_sevenn_readiness(self):
        """Verify Ir/Rh/As/Ga/In/Ru identifier curation, expanded master DB, SevenNet atom constants, and computed features."""
        focus_elements = {"Ir", "Rh", "As", "Ga", "In", "Ru"}
        sevenn_mod = self._load_script_module(
            os.path.join("scripts", "benchmark", "run_phase2_sevenn_features.py"),
            "repo_run_phase2_sevenn_features",
        )
        for el in focus_elements:
            self.assertIn(el, sevenn_mod.ISOLATED_ATOM_ENERGIES)
            self.assertIn(el, sevenn_mod.BONDI_RADII)
            self.assertIn(el, sevenn_mod.METAL_SET)

        curated_path = os.path.join(PROJECT_ROOT, "database", "external", "curated_target_element_identifiers.json")
        with open(curated_path, "r", encoding="utf-8") as f:
            curated = json.load(f)["identifiers"]
        self.assertGreaterEqual(len(curated), 75)
        curated_els = set()
        for name, item in curated.items():
            curated_els.add(item["center_element"])
            mol = Chem.MolFromSmiles(item["smiles"])
            self.assertIsNotNone(mol, f"Invalid curated SMILES for {name}")
            self.assertEqual(rdMolDescriptors.CalcMolFormula(mol), item["formula"], f"Formula mismatch for {name}")
        self.assertEqual(curated_els, focus_elements)

        master_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_semiconductor_precursors_master.json")
        with open(master_path, "r", encoding="utf-8") as f:
            master = json.load(f)
        self.assertGreaterEqual(len(master), 205)
        master_focus_with_smi = {r["center_element"] for r in master if r.get("center_element") in focus_elements and r.get("smiles")}
        self.assertEqual(master_focus_with_smi, focus_elements)

        exp_sevenn_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_elements_sevenn_features.json")
        if os.path.exists(exp_sevenn_path):
            with open(exp_sevenn_path, "r", encoding="utf-8") as f:
                exp_sevenn = json.load(f)
            feats = exp_sevenn["features"]
            self.assertEqual(len(feats), 75)
            self.assertEqual(exp_sevenn["metadata"].get("fallback_count", -1), 0)
            for feat in feats:
                self.assertFalse(feat.get("is_fallback", True), f"Fallback triggered for {feat['id']}")
                self.assertGreaterEqual(feat.get("dimer_init_gap_ang", 0.0), 2.8)
                self.assertGreater(feat["e_coh_ev_per_atom"], 0.0)
                for key in ("e_pot_ev", "e_coh_ev_per_atom", "bde_min_kcal_mol", "vbur_percent", "delta_e_dimer_kcal_mol"):
                    self.assertTrue(math.isfinite(float(feat[key])), f"Non-finite {key} for {feat['id']}")

        exp_neg_sevenn_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_neg_sevenn_features.json")
        if os.path.exists(exp_neg_sevenn_path):
            with open(exp_neg_sevenn_path, "r", encoding="utf-8") as f:
                exp_neg_sevenn = json.load(f)
            neg_feats = exp_neg_sevenn["features"]
            self.assertEqual(len(neg_feats), 30)
            self.assertEqual(exp_neg_sevenn["metadata"].get("fallback_count", -1), 0)
            for feat in neg_feats:
                self.assertFalse(feat.get("is_fallback", True), f"Fallback triggered for {feat['id']}")
                self.assertGreater(feat["e_coh_ev_per_atom"], 0.0)

    def test_14_expanded_18_model_semiconductor_db_evaluation(self):
        """Verify 18-model evaluation across 102 SMILES precursors, 18 inhibitors, and 29 reagents."""
        eval_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json")
        self.assertTrue(os.path.exists(eval_path), "Missing database/benchmark_expanded_models_db.json")
        with open(eval_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        meta = data["metadata"]
        self.assertEqual(meta["total_models"], 18)
        self.assertEqual(meta["precursors_smiles_verified"], 102)
        self.assertEqual(meta["inhibitors_smiles_verified"], 18)
        self.assertEqual(meta["reagents_smiles_verified"], 29)
        self.assertEqual(meta["total_smiles_evaluated"], 149)

        p_counts = data["precursors_102"]["counts"]
        self.assertEqual(p_counts["total"], 102)
        self.assertEqual(p_counts["rdkit_valid"], 102)
        self.assertEqual(p_counts["sascore_pass"], 2)
        self.assertEqual(p_counts["gasa_supported"] + p_counts["gasa_silent_misencoded"], 102)
        self.assertEqual(p_counts["dfrscore_else_collapsed"], 102)
        self.assertEqual(p_counts["synnet_hb_core"], 0)
        self.assertEqual(p_counts["synformer_core"], 0)
        self.assertEqual(p_counts["rxnflow_core"], 0)
        self.assertEqual(len(data["precursors_102"]["records"]), 102)

    def test_15_neg_control_descriptor_and_model_search(self):
        """Verify expanded descriptors, OneClassSVM_RBF, and physics coordination gate results."""
        path = os.path.join(PROJECT_ROOT, "database", "benchmark_neg_control_model_search.json")
        self.assertTrue(os.path.exists(path), "Missing database/benchmark_neg_control_model_search.json")
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        self.assertIn("is_homoleptic_tm_sp3", data["new_descriptors_added"])
        gate = data["physics_coord_gate_eval"]
        self.assertEqual(gate["cohort_55"]["precursor_27"]["with_coord_gate_pass"], 22)
        self.assertEqual(gate["cohort_55"]["inhibitor_18"]["with_coord_gate_pass"], 16)
        self.assertEqual(gate["cohort_55"]["negative_control_10"]["with_coord_gate_pass"], 0)
        self.assertEqual(gate["negative_control_tnr_specificity"], 1.0)
        self.assertEqual(data["cohort_130_chem_coord"]["best_model_by_cv5_bal_acc"], "OneClassSVM_RBF")

    def test_16_refinement_and_rxn_check(self):
        """Verify Item 1 adaptive V_min, Item 4 pi-hapticity/chelate gate, and Item 2 SevenNet-0 reaction energies."""
        path = os.path.join(PROJECT_ROOT, "database", "benchmark_refinement_and_rxn_check.json")
        self.assertTrue(os.path.exists(path), "Missing database/benchmark_refinement_and_rxn_check.json")
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        v2 = data["items_1_and_4"]["cohort_55_v2"]
        self.assertEqual(v2["precursor_pass"], "25/27")
        self.assertEqual(v2["neg_control_pass"], "0/10")
        self.assertEqual(v2["neg_control_tnr_pct"], 100.0)

        exp75 = data["items_1_and_4"]["expanded_75_topology"]
        self.assertEqual(exp75["total"], 75)
        self.assertEqual(exp75["coord_flagged_count"], 0)
        self.assertEqual(exp75["unflagged_count"], 75)

        rxns = {r["id"]: r for r in data["item_2_reaction_energy_check"]["reactions"]}
        self.assertLess(rxns["SYN_TEOS"]["delta_e_kcal_mol"]["SevenNet-0"], 0.0)
        self.assertLess(rxns["SYN_TDMAT"]["delta_e_with_solid_licl_kcal_mol"]["SevenNet-0"], -100.0)
        self.assertLess(rxns["SYN_TEMAH"]["delta_e_with_solid_licl_kcal_mol"]["SevenNet-0"], -100.0)

    def test_17_expanded_neg_control_validation(self):
        """Verify Proposal 5 failure-mechanism-stratified 40 Negative Control & 160-cohort validation."""
        db_path = os.path.join(PROJECT_ROOT, "database", "expanded_negative_controls.json")
        res_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_neg_control_validation.json")
        self.assertTrue(os.path.exists(db_path), "Missing database/expanded_negative_controls.json")
        self.assertTrue(os.path.exists(res_path), "Missing database/benchmark_expanded_neg_control_validation.json")
        with open(res_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        cs = data["cohort_summary"]
        self.assertEqual(cs["total_molecules"], 160)
        self.assertEqual(cs["positive_total"], 120)
        self.assertEqual(cs["negative_controls_total"], 40)

        base = data["negative_control_40_baselines"]
        self.assertEqual(base["sascore_rejected"], "18/40")
        self.assertEqual(base["syba_rejected"], "4/40")
        self.assertEqual(base["coord_chem_gate_rejected"], "39/40")
        self.assertEqual(base["coord_chem_gate_tnr_pct"], 97.5)

        self.assertEqual(len(data["by_failure_mechanism"]), 6)
        ml = data["ml_cv5_cohort_160"]
        self.assertGreaterEqual(ml["SVC_RBF"]["cv5_balanced_acc"], 0.90)
        self.assertGreaterEqual(ml["LogisticRegression"]["cv5_tnr_specificity"], 0.90)

        lofmo = data["lofmo_extrapolation"]
        self.assertLess(lofmo["LogisticRegression"]["mean_unseen_tnr"], 0.10)
        self.assertLess(lofmo["SVC_RBF"]["mean_unseen_tnr"], 0.25)
        self.assertGreaterEqual(lofmo["OneClassSVM_RBF"]["mean_unseen_tnr"], 0.80)
        self.assertGreaterEqual(lofmo["Coord_Chem_Gate"]["mean_unseen_tnr"], 0.95)

    def test_18_unified_latest_160_analysis(self):
        """Verify unified latest 160-molecule evaluation JSON and clean HTML terminology/plots."""
        res_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
        self.assertTrue(os.path.exists(res_path), "Missing database/benchmark_unified_latest_160.json")
        with open(res_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        ds = data["dataset_summary"]
        self.assertEqual(ds["total_molecules"], 160)
        self.assertEqual(ds["precursors_total"], 102)
        self.assertEqual(ds["inhibitors_total"], 18)
        self.assertEqual(ds["negative_controls_total"], 40)
        self.assertEqual(ds["sevenn_relaxed_count"], 160)
        self.assertEqual(ds["etkdg_3d_count"], 0)
        self.assertEqual(ds["sevenn_relaxed_count"] + ds["etkdg_3d_count"], 160)
        self.assertEqual(ds["feature_count"], 38)
        for req_feat in (
            "vbur_percent",
            "vbur_d2_percent",
            "void_cone_deg",
            "r_cov_ang",
            "r_probe_ang",
            "r_eff_ang",
            "r_probe_d2_ang",
            "n_val",
            "chi_pauling",
            "chi_eff",
            "cn_eff",
            "n_ve",
            "ve_deficit",
            "coord_propensity",
        ):
            self.assertIn(req_feat, ds["feature_names"])
        self.assertNotIn("v_min_adaptive", ds["feature_names"])
        self.assertEqual(len(data["records"]), 160)

        v2 = data["unified_v2_summary"]
        self.assertEqual(v2["precursors_all_tp"], "95/102")
        self.assertAlmostEqual(v2["precursors_all_tpr_pct"], 93.1, places=1)
        self.assertEqual(v2["audited_precursors_tp"], "25/27")
        self.assertEqual(v2["expanded_precursors_tp"], "70/75")
        self.assertEqual(v2["negative_all_tn"], "40/40")
        self.assertAlmostEqual(v2["negative_all_tnr_pct"], 100.0, places=1)
        self.assertEqual(v2["negative_fp_ids"], [])
        self.assertNotIn("NEG-019", v2["negative_fp_ids"])
        self.assertNotIn("NEG-039", v2["negative_fp_ids"])
        self.assertNotIn("NEG-040", v2["negative_fp_ids"])

        lofmo = data["lofmo_extrapolation_160"]
        self.assertLess(lofmo["LogisticRegression"]["mean_unseen_tnr"], 0.10)
        self.assertGreaterEqual(lofmo["OneClassSVM_RBF"]["mean_unseen_tnr"], 0.85)
        self.assertEqual(lofmo["Unified_TMC_SA_v2_Gate"]["mean_unseen_tnr"], 1.0)

        html_path = os.path.join(PROJECT_ROOT, "index.html")
        with open(html_path, "r", encoding="utf-8") as f:
            html_txt = f.read()
        self.assertNotIn("전구체", html_txt)
        self.assertNotIn("코호트", html_txt)
        self.assertNotIn("cohort", html_txt.lower())
        self.assertNotIn("24 Features", html_txt)
        self.assertNotIn("27 Features", html_txt)
        self.assertIn("38 Features", html_txt)
        self.assertIn("40/40", html_txt)
        self.assertIn("[Descriptor 4 Physical Validation]", html_txt)
        self.assertIn("전이금속 다중상태 배위 인식 원리", html_txt)
        self.assertIn("[Unified Plot 1]", html_txt)
        self.assertIn("[Unified Plot 2]", html_txt)
        for idx, toc_step in enumerate((
            "반도체 프리커서 SA 모델 평가 기준 (Confusion Matrix)",
            "4세대 AI Model 체계와 반도체 프리커서 적합성 판정",
            "기존 모델 오분류 분석 및 SevenNet-0 물리 Descriptor 체계",
            "통합 160종 물리-배위 모델(TMC-SA) 평가 결과",
        ), start=1):
            self.assertEqual(html_txt.count(f">{idx}. {toc_step}<"), 1, f"Storyline TOC missing: {toc_step}")
            self.assertEqual(html_txt.count(f'<h3 class="storyline-step-title">{toc_step}</h3>'), 1, f"Storyline body heading mismatch: {toc_step}")

    def test_19_design1_element_physics_and_descriptor4_replacement(self):
        """Verify Design 1+2 multi-state TM coordination physics and Descriptor 4 orthogonality."""
        bench_dir = os.path.join(PROJECT_ROOT, "scripts", "benchmark")
        if bench_dir not in sys.path:
            sys.path.insert(0, bench_dir)
        from element_physics import (
            DELTA_R_SHELL_ANG,
            compute_intrinsic_topology,
            get_element_properties,
            universal_vbur_lower_bound,
        )

        self.assertEqual(universal_vbur_lower_bound(0), 44.0)
        self.assertEqual(universal_vbur_lower_bound(1), 39.0)
        self.assertEqual(universal_vbur_lower_bound(2), 34.0)
        self.assertEqual(universal_vbur_lower_bound(5), 34.0)

        # Element-property center detection on CrPh3 without METAL_SET whitelist
        topo_cr = compute_intrinsic_topology("[Cr](c1ccccc1)(c1ccccc1)c1ccccc1")
        self.assertEqual(topo_cr["center_atom"], "Cr")
        self.assertTrue(topo_cr["is_coord_center"])
        self.assertTrue(topo_cr["is_d_block"])
        self.assertTrue(topo_cr["coord_flagged"])

        # Design 2 intrinsic hydride disproportionation & pi-conjugated hydride gates
        topo_si2cl4h2 = compute_intrinsic_topology("Cl[SiH](Cl)[SiH](Cl)Cl")
        self.assertTrue(topo_si2cl4h2["catenated_main_group_hydride"])
        topo_cp2gah = compute_intrinsic_topology("[GaH](C1C=CC=C1)C1C=CC=C1")
        self.assertTrue(topo_cp2gah["pi_conjugated_main_group_hydride"])

        # Continuous covalent-radius scaling: Si < Ti < Zr
        p_si = get_element_properties("Si", 14)
        p_ti = get_element_properties("Ti", 22)
        p_zr = get_element_properties("Zr", 40)
        self.assertLess(p_si["r_cov_ang"], p_ti["r_cov_ang"])
        self.assertLess(p_ti["r_cov_ang"], p_zr["r_cov_ang"])
        self.assertAlmostEqual(DELTA_R_SHELL_ANG, 1.85, places=2)

        # Verify Design 1+2 specs & Descriptor 4 replacement in unified 160 JSON
        res_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
        with open(res_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        spec = data["design1_model_spec"]
        self.assertIn("No Element-Symbol Branching", spec["mode"])
        self.assertIn("design2_model_spec", data)
        self.assertIn("tm_coordination_vector", data["design2_model_spec"])

        d4_85 = data["descriptor4_audit_85"]["full_85"]
        d4_89 = data["descriptor4_audit_85"]["full_89"]
        loocv_160 = data["descriptor4_audit_85"]["loocv_89_marginal_contribution"]
        self.assertEqual(d4_89["n_total"], 160)
        self.assertGreater(d4_85["delta_e_form_per_bond"]["auroc"], 0.70)
        self.assertGreater(d4_85["delta_e_form_per_bond"]["auroc"], d4_85["legacy_proxy"]["auroc"] + 0.20)
        self.assertGreater(d4_89["delta_e_form_per_bond"]["auroc"], d4_89["legacy_proxy"]["auroc"] + 0.25)
        self.assertLess(abs(d4_89["delta_e_form_per_bond"]["corr_bde_vbur_ddim"][0]), 0.25)
        self.assertGreater(loocv_160["D1_D3_plus_dEform_plus_dGdimer"]["auroc"], loocv_160["D1_D3_only"]["auroc"])

        by_id = {r["id"]: r for r in data["records"]}
        self.assertTrue(by_id["PREC_SI_01"]["unified_pass"])
        self.assertEqual(by_id["PREC_SI_01"]["cn_eff"], 4.0)
        self.assertEqual(by_id["PREC_SI_01"]["n_ve"], 8.0)
        self.assertTrue(by_id["PREC_MET_03"]["unified_pass"])
        self.assertTrue(by_id["EXT-ALD-115"]["unified_pass"])

        # State-dependent TM adaptation: high-ox RuO4 contracts R_eff and raises chi_eff vs Ru(EtCp)2
        self.assertLess(by_id["EXT-ALD-115"]["r_eff_ang"], by_id["PREC_MET_03"]["r_eff_ang"])
        self.assertGreater(by_id["EXT-ALD-115"]["chi_eff"], by_id["PREC_MET_03"]["chi_eff"])

        # All 40 Negative Controls blocked including CrPh3, Cl2SiH-SiHCl2, and Cp2GaH
        self.assertFalse(by_id["NEG-019"]["unified_pass"])
        self.assertFalse(by_id["NEG-039"]["unified_pass"])
        self.assertFalse(by_id["NEG-040"]["unified_pass"])

    def test_20_ligand_coord_html_section_and_report_sync(self):
        """Verify Section 18 multi-coordination & ligand compositional stats module and HTML sync."""
        bench_dir = os.path.join(PROJECT_ROOT, "scripts", "benchmark")
        if bench_dir not in sys.path:
            sys.path.insert(0, bench_dir)
        from ligand_coord_html_section import (
            analyze_precursor_coordination_and_ligands,
            build_ligand_coord_storyline_html,
            build_ligand_coord_technical_html,
        )

        unified_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
        stats = analyze_precursor_coordination_and_ligands(unified_path)
        self.assertEqual(stats["n_precursors"], 102)
        self.assertEqual(stats["n_unique_metals"], 16)
        self.assertEqual(stats["n_unique_ligands"], 51)
        self.assertEqual(stats["homoleptic_count"] + stats["heteroleptic_count"], 102)
        self.assertEqual(stats["n_shared_ligands"] + stats["n_singleton_ligands"], 51)

        ab = {r["key"]: r for r in stats["ablation_rows"]}
        self.assertGreater(ab["vbur_percent"]["coup_cv5_r2"], ab["vbur_percent"]["disc_cv5_r2"] + 0.5)
        self.assertGreater(ab["delta_e_form_per_bond_kcal_mol"]["cont_lomo_r2"], 0.6)

        story_html = build_ligand_coord_storyline_html(unified_path, stats=stats)
        tech_html = build_ligand_coord_technical_html(unified_path, stats=stats)
        for chunk in (story_html, tech_html):
            for ctrl in ("\x07", "\x08", "\x0c", "{{", "\t"):
                self.assertNotIn(ctrl, chunk)
            self.assertNotIn("전구체", chunk)
            self.assertNotIn("코호트", chunk)
            self.assertNotIn("cohort", chunk.lower())

        for fname in ("index.html", "synthetic_accessibility_report.html"):
            fpath = os.path.join(PROJECT_ROOT, fname)
            with open(fpath, "r", encoding="utf-8") as fh:
                html_txt = fh.read()
            self.assertIn("4-2. 전이금속 다중 배위수(", html_txt)
            self.assertIn('id="tm-coord-ligand-stats"', html_txt)
            self.assertIn("Sec 1/22", html_txt)
            self.assertIn("18. 전이금속 다중 배위수 결정 및 리간드 조합 통계", html_txt)
            self.assertIn("CCTBA", html_txt)
            self.assertIn("Co(AMD)2", html_txt)

    def test_21_inverse_design_model_section_and_benchmark(self):
        """Verify Section 19 CoOrd-Inv inverse-design model, out-of-fold recall, and HTML sync."""
        bench_dir = os.path.join(PROJECT_ROOT, "scripts", "benchmark")
        if bench_dir not in sys.path:
            sys.path.insert(0, bench_dir)
        from inverse_design_model import (
            build_inverse_design_storyline_html,
            build_inverse_design_technical_html,
            compute_inverse_design_analysis,
        )

        unified_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
        inv = compute_inverse_design_analysis(unified_path)
        self.assertEqual(inv["n_metals"], 16)
        self.assertGreater(inv["n_candidates"], 30000)
        self.assertGreater(inv["n_hard_pass"], 3000)
        self.assertGreaterEqual(inv["retro_cv5"]["coupled"]["n"], 50)
        # CoOrd-Inv coupled surrogate must strictly outperform discrete Bag(L) baseline
        self.assertLess(inv["retro_cv5"]["coupled"]["median_pct"], inv["retro_cv5"]["discrete"]["median_pct"])
        self.assertLess(inv["retro_lomo5"]["coupled"]["median_pct"], inv["retro_lomo5"]["discrete"]["median_pct"])
        self.assertGreater(inv["retro_cv5"]["coupled"]["top5_pct"], inv["retro_cv5"]["discrete"]["top5_pct"])
        self.assertEqual(len(inv["novel_top"]), 10)

        story_html = build_inverse_design_storyline_html(unified_path, inv_stats=inv)
        tech_html = build_inverse_design_technical_html(unified_path, inv_stats=inv)
        for chunk in (story_html, tech_html):
            for ctrl in ("\x07", "\x08", "\x0c", "{{", "\t"):
                self.assertNotIn(ctrl, chunk)
            self.assertNotIn("전구체", chunk)
            self.assertNotIn("코호트", chunk)
            self.assertNotIn("cohort", chunk.lower())

        for fname in ("index.html", "synthetic_accessibility_report.html"):
            fpath = os.path.join(PROJECT_ROOT, fname)
            with open(fpath, "r", encoding="utf-8") as fh:
                html_txt = fh.read()
            self.assertIn("Step 4-3", html_txt)
            self.assertIn('id="precursor-inverse-model-synthesis"', html_txt)
            self.assertIn("19. 반도체 프리커서 역설계 모델 통합 (CoOrd-Inv)", html_txt)
            self.assertNotIn("<!-- INVERSE_DESIGN_STORYLINE_PLACEHOLDER -->", html_txt)
            self.assertNotIn("<!-- INVERSE_DESIGN_TECHNICAL_PLACEHOLDER -->", html_txt)

    def test_22_mechanistic_ml_model_and_generalization(self):
        """Verify Section 20 mechanistic M+L models, 4-regime generalization, and HTML sync."""
        bench_dir = os.path.join(PROJECT_ROOT, "scripts", "benchmark")
        if bench_dir not in sys.path:
            sys.path.insert(0, bench_dir)
        from mechanistic_ml_model import (
            build_mechanistic_storyline_html,
            build_mechanistic_technical_html,
            evaluate_mechanistic_ml_models,
        )

        unified_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
        neb_path = os.path.join(PROJECT_ROOT, "database", "external", "sevenn_decomp_neb_benchmark.json")
        stats = evaluate_mechanistic_ml_models(unified_path, neb_path)

        self.assertEqual(stats["n_precs"], 102)
        self.assertEqual(stats["n_homo"], 54)
        self.assertEqual(stats["n_hetero"], 48)
        self.assertEqual(stats["n_viab"], 142)
        self.assertGreater(stats["ef_full_r2"], 0.70)
        self.assertGreater(stats["vb_full_r2"], 0.75)

        reg_map = {r["key"]: r for r in stats["reg_rows"]}
        # Physics laws must strictly beat RF in Homo->Hetero extrapolation for ΔEf/vM and %V_bur
        self.assertGreater(
            reg_map["delta_e_form_per_bond_kcal_mol"]["homo_het_phys_r2"],
            reg_map["delta_e_form_per_bond_kcal_mol"]["homo_het_rf_r2"],
        )
        self.assertGreater(
            reg_map["vbur_percent"]["homo_het_phys_r2"],
            reg_map["vbur_percent"]["homo_het_rf_r2"],
        )
        # BDE_min must remain negative on Homo->Hetero (honest failure boundary)
        self.assertLess(reg_map["bde_min_kcal_mol"]["homo_het_rf_r2"], 0.0)

        cls_2d = stats["cls_results"]["only_2d_ml"]
        cls_full = stats["cls_results"]["full_ml_3d"]
        self.assertGreater(cls_2d["cv5"]["auroc"], 0.98)
        self.assertGreater(cls_2d["lomo5"]["auroc"], 0.95)
        self.assertGreater(cls_2d["lolo5"]["auroc"], 0.90)
        self.assertGreaterEqual(cls_2d["lolo5"]["auroc"], cls_full["lolo5"]["auroc"])
        self.assertEqual(len(stats["motif_rows"]), 5)
        for m in stats["motif_rows"]:
            self.assertNotEqual(m["neb_anchor"].strip(), "-")
            self.assertNotIn("-;", m["neb_anchor"])
            self.assertNotIn("; -", m["neb_anchor"])
            self.assertRegex(m["neb_anchor"], r"\d+\.\d{2}")

        story_html = build_mechanistic_storyline_html(unified_path, neb_path, mech_stats=stats)
        tech_html = build_mechanistic_technical_html(unified_path, neb_path, mech_stats=stats)
        for chunk in (story_html, tech_html):
            for ctrl in ("\x07", "\x08", "\x0c", "{{", "\t"):
                self.assertNotIn(ctrl, chunk)
            self.assertNotIn("전구체", chunk)
            self.assertNotIn("코호트", chunk)
            self.assertNotIn("cohort", chunk.lower())

        for fname in ("index.html", "synthetic_accessibility_report.html"):
            fpath = os.path.join(PROJECT_ROOT, fname)
            with open(fpath, "r", encoding="utf-8") as fh:
                html_txt = fh.read()
            self.assertIn("Step 4-4", html_txt)
            self.assertIn('id="ml-mechanistic-generalization"', html_txt)
            self.assertIn("20. 중심원자&ndash;리간드 물리&middot;메커니즘 모델 및 일반화 검증", html_txt)
            self.assertNotIn("<!-- MECHANISTIC_ML_STORYLINE_PLACEHOLDER -->", html_txt)
            self.assertNotIn("<!-- MECHANISTIC_ML_TECHNICAL_PLACEHOLDER -->", html_txt)

    def test_23_analysis_history_and_patent_expansion(self):
        """Verify Section 21 analysis history and Section 22 / Step 4-5 KR/JP/CN/US/TW expansion."""
        bench_dir = os.path.join(PROJECT_ROOT, "scripts", "benchmark")
        if bench_dir not in sys.path:
            sys.path.insert(0, bench_dir)
        from history_and_expansion_html_section import (
            build_analysis_history_technical_html,
            build_patent_expansion_storyline_html,
            build_patent_expansion_technical_html,
        )

        unified_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
        neb_path = os.path.join(PROJECT_ROOT, "database", "external", "sevenn_decomp_neb_benchmark.json")
        exp_path = os.path.join(PROJECT_ROOT, "database", "external", "patent_lit_expansion_candidates.json")

        exp_7n_path = os.path.join(PROJECT_ROOT, "database", "external", "patent_expansion_sevenn_features.json")
        with open(exp_7n_path, "r", encoding="utf-8") as fh:
            exp_7n_data = json.load(fh)
        self.assertEqual(exp_7n_data["metadata"]["total_molecules_computed"], 38)
        self.assertEqual(exp_7n_data["metadata"]["fallback_count"], 0)
        self.assertEqual(len(exp_7n_data["features"]), 38)

        with open(exp_path, "r", encoding="utf-8") as fh:
            exp_data = json.load(fh)
        s = exp_data["summary"]
        self.assertEqual(s["n_total"], 38)
        self.assertEqual(s["n_precursors"], 28)
        self.assertEqual(s["n_inhibitors"], 10)
        self.assertEqual(s["n_new_elements_added"], 8)
        self.assertEqual(s["n_merged_elements"], 24)
        self.assertEqual(s["n_inhibitors_pass"], 10)
        self.assertEqual(s["n_precursors_pass"], 24)
        self.assertEqual(s["n_sevenn_computed"], 38)
        self.assertEqual(s["n_sevenn_fallback"], 0)
        self.assertEqual(s["n_precursors_pass_7net"], 15)
        self.assertEqual(s["n_inhibitors_pass_7net"], 8)

        # Ensure zero canonical SMILES duplication with existing 160-molecule DB & verify 7net rescue
        from rdkit import Chem
        with open(unified_path, "r", encoding="utf-8") as fh:
            ex_cans = {
                Chem.MolToSmiles(Chem.MolFromSmiles(r["smiles"]))
                for r in json.load(fh)["records"]
                if Chem.MolFromSmiles(r["smiles"]) is not None
            }
        rec_by_name = {rec["name"]: rec for rec in exp_data["records"]}
        for rec in exp_data["records"]:
            self.assertTrue(rec["embed_ok"])
            self.assertNotIn(rec["canonical_smiles"], ex_cans)
            self.assertIsNotNone(rec.get("tmc_sa_7net"))
            self.assertIsNotNone(rec.get("bde_min_kcal_mol"))
            self.assertIsNotNone(rec.get("delta_e_form_per_bond_kcal_mol"))
            self.assertFalse(rec.get("sevenn_is_fallback"))
        self.assertTrue(rec_by_name["PDMAT"]["pass_7net"])
        self.assertTrue(rec_by_name["Ge(NMe2)4"]["pass_7net"])

        hist_html = build_analysis_history_technical_html(unified_path, neb_path, exp_path)
        exp_story = build_patent_expansion_storyline_html(exp_path, unified_path)
        exp_tech = build_patent_expansion_technical_html(exp_path, unified_path)
        for chunk in (hist_html, exp_story, exp_tech):
            for ctrl in ("\x07", "\x08", "\x0c", "{{", "\t"):
                self.assertNotIn(ctrl, chunk)
            self.assertNotIn("전구체", chunk)
            self.assertNotIn("코호트", chunk)
            self.assertNotIn("cohort", chunk.lower())

        for fname in ("index.html", "synthetic_accessibility_report.html"):
            fpath = os.path.join(PROJECT_ROOT, fname)
            with open(fpath, "r", encoding="utf-8") as fh:
                html_txt = fh.read()
            self.assertIn("4-5. 한&middot;중&middot;일&middot;미&middot;대만 특허", html_txt)
            self.assertIn('id="analysis-history-roadmap"', html_txt)
            self.assertIn('id="patent-literature-expansion"', html_txt)
            self.assertIn("21. 연구&middot;분석 이력 및 단계별 검증 로드맵", html_txt)
            self.assertIn("22. 한&middot;중&middot;일&middot;미&middot;대만 특허 및 중심금속별 문헌 DB 확장", html_txt)
            self.assertNotIn("<!-- ANALYSIS_HISTORY_TECHNICAL_PLACEHOLDER -->", html_txt)
            self.assertNotIn("<!-- PATENT_EXPANSION_STORYLINE_PLACEHOLDER -->", html_txt)
            self.assertNotIn("<!-- PATENT_EXPANSION_TECHNICAL_PLACEHOLDER -->", html_txt)
            self.assertIn("38 Patent &amp; Lit Expansion", html_txt) if "38 Patent &amp; Lit Expansion" in html_txt else self.assertIn("38 Patent & Lit Expansion", html_txt)
            self.assertIn('data-dbcat="patent_expansion"', html_txt)

        # Verify Canonical Deduplicated DB Explorer payload inside index.html
        import re as _re
        m_db = _re.search(r'<script[^>]*id="db-explorer-payload"[^>]*>([\s\S]*?)</script>', html_txt)
        self.assertIsNotNone(m_db)
        db_payload = json.loads(m_db.group(1).replace("<\\/", "</"))
        db_recs = db_payload["records"]
        db_meta = db_payload["meta"]
        self.assertEqual(db_meta["n_patent_total"], 38)
        self.assertEqual(db_meta["n_backfilled_atomiclimits"], 14)
        self.assertEqual(db_meta["n_net_new_patent"], 24)
        self.assertEqual(db_meta["n_merged_dual_role_reagents"], 10)
        self.assertEqual(db_meta["n_merged_notation_aliases"], 2)
        self.assertEqual(len(db_recs), 286)

        # Assert zero canonical SMILES duplicates across all non-reaction molecular DB records
        non_rxn_cs = []
        for rec in db_recs:
            if rec.get("db_type") == "reaction":
                continue
            smi = rec.get("smiles") or ""
            if smi:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    non_rxn_cs.append(Chem.MolToSmiles(mol, isomericSmiles=False))
        self.assertEqual(len(non_rxn_cs), len(set(non_rxn_cs)))

        # Verify MoCl5 appears exactly once as unified canonical entity linking EXT-ALD-058, BB_11, and EXP_02
        mocl5_hits = [rec for rec in db_recs if rec.get("formula") == "Cl5Mo" or "mocl5" in str(rec.get("name", "")).lower()]
        self.assertEqual(len(mocl5_hits), 1)
        self.assertEqual(mocl5_hits[0]["id"], "EXT-ALD-058")
        self.assertIn("BB_11", mocl5_hits[0].get("secondary_ids", []))
        self.assertIn("EXP_02", mocl5_hits[0].get("secondary_ids", []))
        self.assertTrue(mocl5_hits[0]["has_2d"])

    def test_24_br_sascore_theory_tab(self):
        """Verify Background Theory tab builder (Section 1 2009 SAScore issues + Section 2 BR-SAScore) and HTML wiring."""
        br_mod = self._load_script_module(
            os.path.join("scripts", "benchmark", "br_sascore_theory_section.py"),
            "repo_br_sascore_theory_section",
        )
        vanilla = br_mod.compute_vanilla_sascore_issues()
        c_map = {(r["complex"], r["rep"]): r["sa_2009"] for r in vanilla["complex_rows"]}
        self.assertAlmostEqual(c_map[("Hf(OMe)3(acac)", "6원 킬레이트 배위형 (Chelate ->)")], 5.2957, places=3)
        self.assertAlmostEqual(c_map[("Hf(acac)4", "8배위 4중 킬레이트 (Chelate ->)")], 5.4854, places=3)
        self.assertAlmostEqual(c_map[("HfCp(OMe)(acac)2", "2중 킬레이트 + eta5 다중배위 Cp")], 6.4305, places=3)
        l_map = {r["pro_name"]: r["sa_2009"] for r in vanilla["ligand_rows"]}
        self.assertAlmostEqual(l_map["MeOH (메탄올)"], 2.7010, places=3)
        self.assertAlmostEqual(l_map["Hacac (keto형)"], 2.1725, places=3)
        self.assertAlmostEqual(l_map["HCp (시클로펜타디엔)"], 3.3069, places=3)

        stats = br_mod.compute_br_pickle_statistics()
        self.assertEqual(stats["B"]["count"], 244572)
        self.assertEqual(stats["R"]["count"], 175186)
        self.assertEqual(stats["BR"]["count"], 331332)
        self.assertEqual(
            stats["overlap"]["b_only"] + stats["overlap"]["r_only"] + stats["overlap"]["both"],
            stats["overlap"]["union"],
        )

        examples = br_mod.compute_br_worked_examples()
        by_name = {ex["name"]: ex for ex in examples}
        self.assertAlmostEqual(by_name["Aspirin"]["br_sascore"], 2.09, places=2)
        self.assertGreater(by_name["TDMAH (Hf(NMe2)4)"]["br_sascore"], 8.0)
        self.assertGreater(by_name["Hf(acac)4"]["br_sascore"], 8.5)
        self.assertAlmostEqual(by_name["MeOH (Methanol)"]["br_sascore"], 1.01, places=2)

        theory_html = br_mod.build_br_sascore_theory_html()
        for ctrl in ("\x07", "\x08", "\x0c", "{{", "\t"):
            self.assertNotIn(ctrl, theory_html)
        self.assertNotIn("전구체", theory_html)
        self.assertNotIn("코호트", theory_html)
        self.assertNotIn("cohort", theory_html.lower())
        self.assertIn('id="brt-sec1"', theory_html)
        self.assertIn('id="brt-sec2"', theory_html)
        self.assertNotIn('id="brt-sec3"', theory_html)
        self.assertIn("5.179 &lt; 5.491 &lt; 6.74", theory_html)
        self.assertIn("Morgan Fingerprint(r&le;2) 출현 빈도", theory_html)
        self.assertIn("재구성 방식별 명시적 계산 비교 (RDKit 실측 항 대입):", theory_html)
        self.assertNotIn("ECFP", theory_html)
        self.assertNotIn("홉", theory_html)
        self.assertNotIn("착물", theory_html)
        self.assertNotIn("금속 복합체", theory_html)
        self.assertNotIn("명시적 수소", theory_html)
        self.assertIn(r"\tag{2}", theory_html)
        self.assertIn(r"\tag{3}", theory_html)
        self.assertIn(r"\tag{4}", theory_html)
        self.assertIn(r"\tag{5}", theory_html)
        self.assertNotIn("코퍼스", theory_html)
        self.assertNotIn("말뭉치", theory_html)
        self.assertIn("반도체용 프리커서에 그대로 적용하면", theory_html)
        self.assertIn("2.3 학습&middot;구축 데이터셋 해부", theory_html)
        self.assertIn("2.6 베이스 모델로 BR-SAScore (2024) 선정 근거", theory_html)
        self.assertIn("2. BR-SAScore (2024) 모델의 원리 및 베이스 모델 선정 근거", theory_html)
        self.assertIn('id="card-brt-pipeline"', theory_html)
        self.assertIn("s_3=0.000", theory_html)
        self.assertIn("합성이 용이(Easy)", theory_html)
        self.assertIn("CpTi(NMe2)3", theory_html)
        self.assertNotIn("CN(C)[Ti](N(C)C)(N(C)C)C1=CC=CC1", theory_html)
        self.assertIn("금속 인접 미등재", theory_html)
        self.assertIn("Cp 시약 확장", theory_html)
        self.assertIn("7.2029", theory_html)
        self.assertIn("6.9458", theory_html)
        self.assertIn("8.9772", theory_html)
        self.assertIn("BR-SAScore 수치 해석 기준", theory_html)
        self.assertIn("*Pro-ligand 분리 (Pro-ligand Decoupling)", theory_html)
        self.assertIn("331,332", theory_html)
        self.assertNotIn('id="brTheorySimWidget"', theory_html)
        self.assertIn("BR-SAScore 베이스 선정 5대 결정적 근거", theory_html)
        self.assertIn("Strict 28/29 (96.6%)", theory_html)
        self.assertIn("0/102 M-L 규칙", theory_html)

        for fname in ("index.html", "synthetic_accessibility_report.html"):
            fpath = os.path.join(PROJECT_ROOT, fname)
            with open(fpath, "r", encoding="utf-8") as fh:
                html_txt = fh.read()
            self.assertIn('id="btn-tab-brtheory"', html_txt)
            self.assertIn("배경이론", html_txt)
            self.assertIn('id="toc-brtheory"', html_txt)
            self.assertIn('id="tab-content-brtheory"', html_txt)
            self.assertNotIn("Background Theory &bull; Base SA Model Foundations", html_txt)
            self.assertNotIn("반도체 SA 배경이론: 2009 SAScore 재구성 이슈부터 BR-SAScore까지", html_txt)
            self.assertIn("반도체 SA 모델을 이해하기 위한 기초 이론 가이드입니다.", html_txt)
            self.assertIn("1. SAScore (2009) 모델의 원리 및 SMILES 재구성 이슈", html_txt)
            self.assertIn("*Morgan Fingerprint", html_txt)
            self.assertNotIn("Morgan Fingerprint(Radius 0~2)란?", html_txt)
            self.assertIn("1.1 SAScore (2009) 모델의 수식과 정규화 규격", html_txt)
            self.assertIn("1.2 SMILES 재구성 이슈", html_txt)
            self.assertIn("단독(미배위) 프로리간드", html_txt)
            self.assertNotIn("무료 프로리간드", html_txt)
            self.assertIn("1.4 Summary", html_txt)
            self.assertIn("renderMermaidInBrtheoryTab", html_txt)
            self.assertNotIn("updateBrSimulator", html_txt)
            self.assertNotIn("<!-- BR_SASCORE_THEORY_PLACEHOLDER -->", html_txt)
            self.assertIn('id="brt-sec1"', html_txt)
            self.assertIn('id="brt-sec2"', html_txt)
            self.assertNotIn('id="brt-sec3"', html_txt)

    def test_25_html_modularization_and_compression(self):
        """Verify Option A modular HTML architecture, shared html_common utilities, and payload compression."""
        hc_mod = self._load_script_module(
            os.path.join("scripts", "benchmark", "html_common.py"),
            "repo_html_common",
        )
        svg_tdmah = hc_mod.render_mol2d_svg("CN(C)[Hf](N(C)C)(N(C)C)N(C)C")
        self.assertIn("<svg class='rdk'", svg_tdmah)
        self.assertNotIn("xmlns:rdkit", svg_tdmah)
        self.assertNotIn("class='bond-", svg_tdmah)

        tpl_dir = os.path.join(PROJECT_ROOT, "scripts", "benchmark", "templates")
        for asset in (
            "report_styles.css",
            "tab_storyline.html",
            "tab_technical.html",
            "tab_database.html",
            "report_main.js",
            "report_neb.js",
        ):
            self.assertTrue(os.path.exists(os.path.join(tpl_dir, asset)), f"Missing template asset: {asset}")

        with open(os.path.join(PROJECT_ROOT, "template.html"), "r", encoding="utf-8") as fh:
            tpl_lines = fh.read().splitlines()
        self.assertLessEqual(len(tpl_lines), 350)

        with open(os.path.join(PROJECT_ROOT, "generate_html.py"), "r", encoding="utf-8") as fh:
            gen_lines = fh.read().splitlines()
        self.assertLessEqual(len(gen_lines), 200)

        for fname in ("index.html", "synthetic_accessibility_report.html"):
            fpath = os.path.join(PROJECT_ROOT, fname)
            with open(fpath, "r", encoding="utf-8") as fh:
                html_txt = fh.read()
            self.assertLess(len(html_txt.encode("utf-8")), 2800000)
            self.assertIn("svg.rdk path:not([fill])", html_txt)
            self.assertIn(".u01{", html_txt)
            for placeholder in (
                "<!-- REPORT_STYLES_PLACEHOLDER -->",
                "<!-- TAB_STORYLINE_SHELL_PLACEHOLDER -->",
                "<!-- TAB_TECHNICAL_SHELL_PLACEHOLDER -->",
                "<!-- TAB_DATABASE_SHELL_PLACEHOLDER -->",
                "<!-- REPORT_MAIN_JS_PLACEHOLDER -->",
                "<!-- REPORT_NEB_JS_PLACEHOLDER -->",
            ):
                self.assertNotIn(placeholder, html_txt)


if __name__ == "__main__":
    unittest.main()
