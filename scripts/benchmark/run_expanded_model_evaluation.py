#!/usr/bin/env python3
"""
Expanded 18-Model Evaluation Engine on Expanded Semiconductor Database
-----------------------------------------------------------------------
Evaluates all 18 SA / Retrosynthesis / Pathway-Generative / Domain ML models
across:
  1. 102 SMILES-verified semiconductor precursors (expanded_semiconductor_precursors_master.json)
  2. 18 AS-ALD surface inhibitors (as_ald_inhibitors.json)
  3. 29 commercial reagents / building blocks (building_blocks.json)

Outputs:
  - database/benchmark_expanded_models_db.json
"""

from __future__ import annotations

import copy
import csv
import json
import os
import pickle
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import xgboost as xgb
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, RDConfig, rdMolDescriptors

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "scscore"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "rascore"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "reasyn"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "external", "synformer"))

import sascorer
from aizynthfinder.chem import Molecule as AiZynthMolecule
from BRSAScore import SAScorer as BRSAScorer
from rdchiral.initialization import rdchiralReactants
from reasyn.chem.featurize import _smiles_token_pattern, _smiles_token_to_id, _smiles_vocab
from scscore.standalone_model_numpy import SCScorer
from syba.syba import SybaClassifier
from synformer.chem.featurize import (
    _smiles_token_max as synformer_token_max,
    tokenize_smiles as synformer_tokenize_smiles,
)

# Expanded metal & metalloid set (includes baseline 12 + curated Ir, Rh, As, Ga, In)
EXPANDED_CENTER_SYMBOLS: Set[str] = {
    "Al", "As", "Co", "Ga", "Hf", "In", "Ir", "Mo",
    "Rh", "Ru", "Si", "Sn", "Ta", "Ti", "W", "Zr",
}
EXPANDED_CENTER_ATOMIC_NUMS: Set[int] = {
    13, 14, 22, 27, 31, 33, 40, 42, 44, 45, 49, 50, 72, 73, 74, 77,
}

# GASA allowable_set from external/gasa/model/gasa_utils.py (43 elements, encode_unknown=False)
GASA_ALLOWED_ELEMENTS: Set[str] = {
    "C", "N", "O", "S", "F", "Si", "P", "Cl", "Br", "Mg", "Na",
    "Ca", "Fe", "As", "Al", "I", "B", "V", "K", "Tl", "Yb", "Sb",
    "Sn", "Ag", "Pd", "Co", "Se", "Ti", "Zn", "H", "Li", "Ge", "Cu",
    "Au", "Ni", "Cd", "In", "Mn", "Zr", "Cr", "Pt", "Hg", "Pb",
}

# DFRscore organic-only vocabulary (10 elements; all others collapse to 'ELSE')
DFR_ALLOWED_ELEMENTS: Set[str] = {"C", "N", "O", "F", "S", "Cl", "Br", "I", "B", "P"}

# MPScore (Jelfs Group, JCIM 2021) organic porous-cage precursor domain elements
MPSCORE_ORGANIC_ELEMENTS: Set[str] = {"C", "H", "N", "O", "F", "S", "B", "P", "Cl", "Br", "I"}

REAGENT_BSCORE_FLOOR = -0.47


def _mol_contains_target_metal(mol_or_patt: Optional[Chem.Mol]) -> bool:
    if mol_or_patt is None:
        return False
    for atom in mol_or_patt.GetAtoms():
        if atom.GetAtomicNum() in EXPANDED_CENTER_ATOMIC_NUMS or atom.GetSymbol() in EXPANDED_CENTER_SYMBOLS:
            return True
    return False


def _canonical_smiles_set(smiles_list: List[str]) -> Set[str]:
    out: Set[str] = set()
    for smi in smiles_list:
        if not smi:
            continue
        out.add(smi)
        m = Chem.MolFromSmiles(smi)
        if m:
            out.add(Chem.MolToSmiles(m))
        if "." in smi:
            for comp in smi.split("."):
                c = comp.strip()
                out.add(c)
                mc = Chem.MolFromSmiles(c)
                if mc:
                    out.add(Chem.MolToSmiles(mc))
    return out


def load_cohorts() -> Dict[str, List[Dict[str, Any]]]:
    master_path = os.path.join(
        PROJECT_ROOT, "database", "external", "expanded_semiconductor_precursors_master.json"
    )
    with open(master_path, "r", encoding="utf-8") as f:
        master_all = json.load(f)
    precursors_102 = [p for p in master_all if p.get("smiles")]

    inh_path = os.path.join(PROJECT_ROOT, "database", "external", "as_ald_inhibitors.json")
    with open(inh_path, "r", encoding="utf-8") as f:
        inh_raw = json.load(f).get("inhibitors", [])
    inhibitors_18 = [x for x in inh_raw if x.get("smiles")]

    bb_path = os.path.join(PROJECT_ROOT, "database", "building_blocks.json")
    with open(bb_path, "r", encoding="utf-8") as f:
        bb_raw = json.load(f)
    reagents_29 = [x for x in bb_raw if x.get("smiles")]

    return {
        "master_total_count": len(master_all),
        "precursors": precursors_102,
        "inhibitors": inhibitors_18,
        "reagents": reagents_29,
    }


def init_models(target_smiles: Set[str]) -> Dict[str, Any]:
    # 1. SCScore
    sc = SCScorer()
    sc.restore(
        os.path.join(
            PROJECT_ROOT,
            "external",
            "scscore",
            "models",
            "full_reaxys_model_1024bool",
            "model.ckpt-10654.as_numpy.json.gz",
        )
    )

    # 2. SYBA
    sy = SybaClassifier()
    sy.fitDefaultScore()

    # 3. RAscore
    ra_model_path = os.path.join(
        PROJECT_ROOT, "external", "rascore", "RAscore", "models", "XGB_chembl_ecfp_counts", "model.pkl"
    )
    with open(ra_model_path, "rb") as f:
        ra_clf = pickle.load(f)
    ra_booster = ra_clf.get_booster()

    # 4. BR-SAScore (Default & Strict Reagent Domain Unleaked)
    br_default = BRSAScorer()
    br_custom_strict = BRSAScorer()
    fscores_strict = copy.deepcopy(br_default._fscores)

    semi_bb_path = os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_bb.csv")
    semi_rxn_path = os.path.join(PROJECT_ROOT, "database", "br_sascore_inputs", "semi_rxn.csv")

    with open(semi_bb_path, "r", encoding="utf-8") as f:
        for row in list(csv.reader(f))[1:]:
            if not row or not row[0].strip():
                continue
            smi = row[0].strip()
            m = Chem.MolFromSmiles(smi)
            if not m:
                continue
            if smi in target_smiles or Chem.MolToSmiles(m) in target_smiles:
                continue
            bi: Dict[int, Any] = {}
            rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
            for b in bi:
                fscores_strict[b] = max(fscores_strict.get(b, -6.0), REAGENT_BSCORE_FLOOR)

    with open(semi_rxn_path, "r", encoding="utf-8") as f:
        for row in list(csv.reader(f))[1:]:
            if not row or not row[0].strip():
                continue
            parts_rxn = row[0].strip().split(">>")
            reactants = parts_rxn[0].strip()
            m_full = Chem.MolFromSmiles(reactants)
            if reactants in target_smiles or (m_full and Chem.MolToSmiles(m_full) in target_smiles):
                continue
            if len(parts_rxn) > 1 and reactants == parts_rxn[1].strip():
                continue
            for comp in reactants.split("."):
                smi_clean = comp.strip()
                m = Chem.MolFromSmiles(smi_clean)
                if not m:
                    continue
                if smi_clean in target_smiles or Chem.MolToSmiles(m) in target_smiles:
                    continue
                bi = {}
                rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
                for b in bi:
                    fscores_strict[b] = max(fscores_strict.get(b, -6.0), REAGENT_BSCORE_FLOOR)

    br_custom_strict._fscores = fscores_strict

    # 5. RxnFlow 72 reverse templates
    rxnflow_patterns: List[Tuple[str, Chem.Mol, bool]] = []
    for scope_name, rel_path in [
        ("explore", os.path.join("external", "rxnflow", "data", "templates", "explore", "reaction.yaml")),
        ("cgflow", os.path.join("external", "rxnflow", "data", "templates", "cgflow", "reaction.yaml")),
    ]:
        t_path = os.path.join(PROJECT_ROOT, rel_path)
        if os.path.exists(t_path):
            with open(t_path, "r", encoding="utf-8") as f:
                d = yaml.safe_load(f)
            for sec in ("UniReaction", "BiReaction"):
                for rule_name, rule_body in d.get(sec, {}).items():
                    rev_smarts = rule_body.get("reverse", "")
                    if ">>" in rev_smarts:
                        prod_smarts = rev_smarts.split(">>")[0].strip()
                        patt = Chem.MolFromSmarts(prod_smarts)
                        if patt:
                            rxnflow_patterns.append(
                                (f"{scope_name}:{rule_name}", patt, _mol_contains_target_metal(patt))
                            )

    # 6. HB-91 templates (ChemProjector / SynNet)
    hb_path = os.path.join(PROJECT_ROOT, "external", "chemprojector", "data", "reaction_templates_hb.txt")
    hb_patterns: List[Tuple[int, Chem.Mol, bool]] = []
    if os.path.exists(hb_path):
        with open(hb_path, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(">>")
                if len(parts) == 2:
                    patt = Chem.MolFromSmarts(parts[1].strip())
                    if patt:
                        hb_patterns.append((i, patt, _mol_contains_target_metal(patt)))

    # 7. SynFormer 115 templates
    synformer_path = os.path.join(
        PROJECT_ROOT, "external", "synformer", "data", "rxn_templates", "comprehensive.txt"
    )
    synformer_reactions: List[Tuple[int, List[Chem.Mol], bool]] = []
    if os.path.exists(synformer_path):
        with open(synformer_path, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rxn = AllChem.ReactionFromSmarts(line)
                if rxn:
                    prods = list(rxn.GetProducts())
                    has_metal = any(_mol_contains_target_metal(pt) for pt in prods) or any(
                        _mol_contains_target_metal(rt) for rt in rxn.GetReactants()
                    )
                    synformer_reactions.append((i, prods, has_metal))

    reasyn_vocab_1_to_tok = {i: t for i, t in enumerate(_smiles_vocab, start=1)}

    return {
        "sc": sc,
        "sy": sy,
        "ra_booster": ra_booster,
        "br_default": br_default,
        "br_custom_strict": br_custom_strict,
        "rxnflow_patterns": rxnflow_patterns,
        "hb_patterns": hb_patterns,
        "synformer_reactions": synformer_reactions,
        "reasyn_vocab_1_to_tok": reasyn_vocab_1_to_tok,
    }


def evaluate_cohort(items: List[Dict[str, Any]], models: Dict[str, Any]) -> Dict[str, Any]:
    n = len(items)
    counts: Dict[str, int] = {
        "total": n,
        "rdkit_valid": 0,
        "sascore_pass": 0,
        "syba_pass": 0,
        "scscore_pass": 0,
        "rascore_pass": 0,
        "retrognn_ecfp_active": 0,
        "mpscore_ecfp4_active": 0,
        "mpscore_organic_in_domain": 0,
        "br_sascore_default_pass": 0,
        "br_sascore_strict_pass": 0,
        "gasa_supported": 0,
        "gasa_silent_misencoded": 0,
        "deepsa_ro5_pass": 0,
        "deepsa_bracket_fragmented": 0,
        "dfrscore_organic_only": 0,
        "dfrscore_else_collapsed": 0,
        "synnet_morgan_active": 0,
        "synnet_hb_substruct": 0,
        "synnet_hb_core": 0,
        "chemprojector_z100_pass": 0,
        "chemprojector_hb_substruct": 0,
        "chemprojector_hb_core": 0,
        "synformer_vocab_pass": 0,
        "synformer_substruct": 0,
        "synformer_ge2_substruct": 0,
        "synformer_core": 0,
        "rxnflow_substruct": 0,
        "rxnflow_core": 0,
        "reasyn_vocab_pass": 0,
        "reasyn_lossless_roundtrip": 0,
        "aizynth_ingest_pass": 0,
        "retro_star_rdchiral_ingest_pass": 0,
        "synflownet_rdkit_ingest_pass": 0,
    }
    records: List[Dict[str, Any]] = []
    gasa_oov_elements: Set[str] = set()
    rascore_vals: List[float] = []
    sascore_vals: List[float] = []
    scscore_vals: List[float] = []
    syba_vals: List[float] = []

    for item in items:
        smi = item["smiles"]
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES in cohort: {smi}"
        counts["rdkit_valid"] += 1

        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        atomic_nums = [a.GetAtomicNum() for a in mol.GetAtoms()]
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)

        # 1. SAScore
        sa_val = round(float(sascorer.calculateScore(mol)), 2)
        sascore_vals.append(sa_val)
        if sa_val <= 3.5:
            counts["sascore_pass"] += 1

        # 2. SYBA
        sy_val = round(float(models["sy"].predict(mol=mol)), 2)
        syba_vals.append(sy_val)
        if sy_val > 0.0:
            counts["syba_pass"] += 1

        # 3. SCScore
        _, sc_raw = models["sc"].get_score_from_smi(smi)
        sc_val = round(float(sc_raw), 2)
        scscore_vals.append(sc_val)
        if sc_val <= 2.0:
            counts["scscore_pass"] += 1

        # 4. RAscore
        fp3 = AllChem.GetMorganFingerprint(mol, 3, useChirality=True)
        arr = np.zeros((2048,), dtype=np.float32)
        for k, v in fp3.GetNonzeroElements().items():
            arr[k % 2048] += v
        ra_val = round(float(models["ra_booster"].predict(xgb.DMatrix(arr.reshape(1, -1)))[0]), 3)
        rascore_vals.append(ra_val)
        if ra_val >= 0.5:
            counts["rascore_pass"] += 1

        # 5. RetroGNN & 6. MPScore (Featurizer / Domain Audit; weights absent)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        if fp2.GetNumOnBits() > 0:
            counts["retrognn_ecfp_active"] += 1
            counts["mpscore_ecfp4_active"] += 1
            counts["synnet_morgan_active"] += 1
        if all(a in MPSCORE_ORGANIC_ELEMENTS for a in atoms):
            counts["mpscore_organic_in_domain"] += 1

        # 7. BR-SAScore
        br_def_val = round(float(models["br_default"].calculateScore(smi)[0]), 2)
        if br_def_val <= 3.5:
            counts["br_sascore_default_pass"] += 1
        br_str_val = round(float(models["br_custom_strict"].calculateScore(smi)[0]), 2)
        if br_str_val <= 3.5:
            counts["br_sascore_strict_pass"] += 1

        # 8. GASA element check
        unsupported_gasa = [a for a in atoms if a not in GASA_ALLOWED_ELEMENTS]
        if not unsupported_gasa:
            counts["gasa_supported"] += 1
        else:
            counts["gasa_silent_misencoded"] += 1
            gasa_oov_elements.update(unsupported_gasa)

        # 9. DeepSA Ro5 & bracket tokenization check
        ro5_viol = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
        if ro5_viol == 0:
            counts["deepsa_ro5_pass"] += 1
        if "[" in smi:
            counts["deepsa_bracket_fragmented"] += 1

        # 10. DFRscore element check
        if any(a not in DFR_ALLOWED_ELEMENTS for a in atoms):
            counts["dfrscore_else_collapsed"] += 1
        else:
            counts["dfrscore_organic_only"] += 1

        # 11. ChemProjector Z<=100 & HB-91
        if all(z <= 100 for z in atomic_nums):
            counts["chemprojector_z100_pass"] += 1
        hb_sub = False
        hb_core = False
        for _, patt, has_metal in models["hb_patterns"]:
            if mol.HasSubstructMatch(patt):
                hb_sub = True
                if has_metal:
                    hb_core = True
        if hb_sub:
            counts["synnet_hb_substruct"] += 1
            counts["chemprojector_hb_substruct"] += 1
        if hb_core:
            counts["synnet_hb_core"] += 1
            counts["chemprojector_hb_core"] += 1

        # 12. SynFormer vocab + 115 comprehensive templates
        sf_tokens = synformer_tokenize_smiles(smi)
        if not any(tid > synformer_token_max for tid in sf_tokens):
            counts["synformer_vocab_pass"] += 1
        sf_sub = False
        sf_ge2 = False
        sf_core = False
        for _, prods, has_metal in models["synformer_reactions"]:
            if all(mol.HasSubstructMatch(pt) for pt in prods):
                sf_sub = True
                if all(pt.GetNumAtoms() >= 2 for pt in prods):
                    sf_ge2 = True
                if has_metal:
                    sf_core = True
        if sf_sub:
            counts["synformer_substruct"] += 1
        if sf_ge2:
            counts["synformer_ge2_substruct"] += 1
        if sf_core:
            counts["synformer_core"] += 1

        # 13. RxnFlow 72 templates
        rf_sub = False
        rf_core = False
        for _, patt, has_metal in models["rxnflow_patterns"]:
            if mol.HasSubstructMatch(patt):
                rf_sub = True
                if has_metal:
                    rf_core = True
        if rf_sub:
            counts["rxnflow_substruct"] += 1
        if rf_core:
            counts["rxnflow_core"] += 1

        # 14. ReaSyn 118-element tokenizer
        tokens = _smiles_token_pattern.findall(smi)
        token_ids = [_smiles_token_to_id.get(t, len(_smiles_vocab) + 1) for t in tokens]
        if not any(tid > len(_smiles_vocab) for tid in token_ids):
            counts["reasyn_vocab_pass"] += 1
        reconstructed = "".join(models["reasyn_vocab_1_to_tok"].get(tid, "") for tid in token_ids)
        if reconstructed == smi:
            counts["reasyn_lossless_roundtrip"] += 1

        # 15. AiZynthFinder ingestion
        try:
            az = AiZynthMolecule(smiles=smi)
            az.sanitize()
            if az.rd_mol is not None:
                counts["aizynth_ingest_pass"] += 1
        except Exception:
            pass

        # 16. Retro* rdchiral ingestion
        try:
            rdc = rdchiralReactants(smi)
            if rdc is not None and rdc.reactants is not None:
                counts["retro_star_rdchiral_ingest_pass"] += 1
        except Exception:
            pass

        # 17. SynFlowNet RDKit Mol ingestion
        counts["synflownet_rdkit_ingest_pass"] += 1

        records.append({
            "id": item.get("id", item.get("name", "")),
            "name": item.get("name", ""),
            "center_element": item.get("center_element", ""),
            "smiles": smi,
            "sascore": sa_val,
            "syba": sy_val,
            "scscore": sc_val,
            "rascore": ra_val,
            "br_sascore_default": br_def_val,
            "br_sascore_strict": br_str_val,
            "gasa_supported": len(unsupported_gasa) == 0,
            "gasa_unsupported_elements": sorted(set(unsupported_gasa)),
            "dfrscore_collapsed": any(a not in DFR_ALLOWED_ELEMENTS for a in atoms),
            "mpscore_organic_in_domain": all(a in MPSCORE_ORGANIC_ELEMENTS for a in atoms),
            "ro5_violations": ro5_viol,
            "reasyn_lossless": reconstructed == smi,
        })

    return {
        "counts": counts,
        "gasa_oov_elements": sorted(gasa_oov_elements),
        "distributions": {
            "sascore_mean": round(float(np.mean(sascore_vals)), 3) if sascore_vals else None,
            "syba_mean": round(float(np.mean(syba_vals)), 3) if syba_vals else None,
            "scscore_mean": round(float(np.mean(scscore_vals)), 3) if scscore_vals else None,
            "rascore_mean": round(float(np.mean(rascore_vals)), 3) if rascore_vals else None,
            "rascore_min": round(float(np.min(rascore_vals)), 3) if rascore_vals else None,
            "rascore_max": round(float(np.max(rascore_vals)), 3) if rascore_vals else None,
        },
        "records": records,
    }


def main() -> None:
    cohorts = load_cohorts()
    all_target_smiles = _canonical_smiles_set([p["smiles"] for p in cohorts["precursors"]])
    models = init_models(all_target_smiles)

    precursor_eval = evaluate_cohort(cohorts["precursors"], models)
    inhibitor_eval = evaluate_cohort(cohorts["inhibitors"], models)
    reagent_eval = evaluate_cohort(cohorts["reagents"], models)

    out = {
        "metadata": {
            "total_models": 18,
            "master_precursor_total": cohorts["master_total_count"],
            "precursors_smiles_verified": len(cohorts["precursors"]),
            "inhibitors_smiles_verified": len(cohorts["inhibitors"]),
            "reagents_smiles_verified": len(cohorts["reagents"]),
            "total_smiles_evaluated": (
                len(cohorts["precursors"]) + len(cohorts["inhibitors"]) + len(cohorts["reagents"])
            ),
        },
        "precursors_102": precursor_eval,
        "inhibitors_18": inhibitor_eval,
        "reagents_29": reagent_eval,
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print("=" * 72)
    print(f"Expanded 18-Model Evaluation Saved -> {out_path}")
    print("Precursors (102 SMILES):", json.dumps(precursor_eval["counts"], indent=2))
    print("GASA OOV elements (102):", precursor_eval["gasa_oov_elements"])
    print("Distributions (102):", precursor_eval["distributions"])
    print("Inhibitors (18 SMILES):", json.dumps(inhibitor_eval["counts"], indent=2))
    print("Reagents (29 SMILES):", json.dumps(reagent_eval["counts"], indent=2))


if __name__ == "__main__":
    main()
