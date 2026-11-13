#!/usr/bin/env python3
"""
Negative Control Separation: Expanded Computational Descriptors & Model Search
------------------------------------------------------------------------------
Adds physically grounded coordination, volatility, and thermal-instability
descriptors (without target-ID hardcoding or label leakage) and benchmarks
multiple linear, kernel, neural, and tree-ensemble classifiers across:
  1. 55-Molecule Audited Cohort (Base SA + SevenNet-0 Physics + New Chem/Coord Descriptors)
  2. 85-Molecule SevenNet-0 Cohort (55 Audited + 30 Expanded Ir/Rh/As/Ga/In/Ru)
  3. 130-Molecule Expanded Cohort (55 Audited + 75 Expanded Ir/Rh/As/Ga/In/Ru, no 7net needed)

Outputs:
  - database/benchmark_neg_control_model_search.json
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import numpy as np
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, IsolationForest, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, recall_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import OneClassSVM, SVC

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CENTER_ATOMIC_NUMS = {13, 14, 22, 27, 31, 33, 40, 42, 44, 45, 49, 50, 72, 73, 74, 77}
TM_ATOMIC_NUMS = {22, 27, 40, 42, 44, 45, 72, 73, 74, 77}

# Physically grounded SMARTS patterns (non-target-specific general chemistry)
PATT_BETA_H_ALKYL = Chem.MolFromSmarts("[#22,#27,#40,#42,#44,#45,#72,#73,#74,#77]-[CH2,CH1]-[CH3,CH2,CH1]")
PATT_METAL_OH = Chem.MolFromSmarts("[#13,#14,#22,#27,#31,#33,#40,#42,#44,#45,#49,#50,#72,#73,#74,#77]-[OX2H1]")
PATT_CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def compute_new_chem_coord_descriptors(smiles: str) -> Dict[str, float]:
    """Compute physically grounded coordination, volatility, and instability descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"

    mw = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd = float(Lipinski.NumHDonors(mol))
    hba = float(Lipinski.NumHAcceptors(mol))
    rotb = float(Lipinski.NumRotatableBonds(mol))
    heavy = float(mol.GetNumHeavyAtoms())
    aro_rings = float(rdMolDescriptors.CalcNumAromaticRings(mol))
    rings = float(rdMolDescriptors.CalcNumRings(mol))

    # Center-atom coordination environment
    center_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() in CENTER_ATOMIC_NUMS]
    n_centers = len(center_atoms)
    metal_c_bonds = 0
    metal_hetero_bonds = 0  # N, O, F, Cl, Br, I, S, P
    total_center_bonds = 0
    for ca in center_atoms:
        for nb in ca.GetNeighbors():
            total_center_bonds += 1
            z = nb.GetAtomicNum()
            if z == 6:
                metal_c_bonds += 1
            elif z in (7, 8, 9, 15, 16, 17, 35, 53):
                metal_hetero_bonds += 1

    metal_c_sigma_frac = float(metal_c_bonds / total_center_bonds) if total_center_bonds > 0 else 0.0
    metal_hetero_frac = float(metal_hetero_bonds / total_center_bonds) if total_center_bonds > 0 else 0.0
    # Homoleptic all-carbon metal/metalloid complex indicator (all center bonds are to C)
    is_homoleptic_mc = 1.0 if (n_centers > 0 and total_center_bonds > 0 and metal_c_bonds == total_center_bonds) else 0.0
    # Homoleptic d-block sp3 sigma-alkyl indicator (distinguishes M(R)_n from pi-Cp/CO/diene organometallics)
    tm_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() in TM_ATOMIC_NUMS]
    is_homoleptic_tm_sp3 = 1.0 if any(
        len(a.GetNeighbors()) > 0
        and all(nb.GetAtomicNum() == 6 and nb.GetHybridization() == Chem.HybridizationType.SP3 for nb in a.GetNeighbors())
        for a in tm_atoms
    ) else 0.0

    beta_h_count = float(len(mol.GetSubstructMatches(PATT_BETA_H_ALKYL))) if PATT_BETA_H_ALKYL else 0.0
    metal_oh_count = float(len(mol.GetSubstructMatches(PATT_METAL_OH))) if PATT_METAL_OH else 0.0
    cooh_count = float(len(mol.GetSubstructMatches(PATT_CARBOXYLIC_ACID))) if PATT_CARBOXYLIC_ACID else 0.0

    # Pure alkane check (only C and H, no rings, no centers)
    atom_syms = {a.GetSymbol() for a in mol.GetAtoms()}
    is_inert_alkane = 1.0 if (atom_syms <= {"C", "H"} and rings == 0 and n_centers == 0) else 0.0

    return {
        "mw": mw,
        "logp": logp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rot_bonds": rotb,
        "heavy_atoms": heavy,
        "aromatic_rings": aro_rings,
        "ring_count": rings,
        "metal_c_sigma_frac": metal_c_sigma_frac,
        "metal_hetero_frac": metal_hetero_frac,
        "is_homoleptic_mc": is_homoleptic_mc,
        "is_homoleptic_tm_sp3": is_homoleptic_tm_sp3,
        "beta_h_alkyl_count": beta_h_count,
        "metal_oh_count": metal_oh_count,
        "cooh_count": cooh_count,
        "is_inert_alkane": is_inert_alkane,
    }


def _build_models(scale_pos_weight: float) -> Dict[str, Any]:
    return {
        "LogisticRegression": make_pipeline(
            StandardScaler(), LogisticRegression(class_weight="balanced", C=1.0, max_iter=2000, random_state=42)
        ),
        "SVC_RBF": make_pipeline(
            StandardScaler(), SVC(kernel="rbf", class_weight="balanced", C=2.0, probability=True, random_state=42)
        ),
        "MLP_NeuralNet": make_pipeline(
            StandardScaler(), MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=1500, alpha=0.5, random_state=42)
        ),
        "RandomForest": RandomForestClassifier(
            n_estimators=200, max_depth=6, class_weight="balanced", random_state=42
        ),
        "ExtraTrees": ExtraTreesClassifier(
            n_estimators=200, max_depth=6, class_weight="balanced", random_state=42
        ),
        "GradientBoosting": GradientBoostingClassifier(
            n_estimators=120, max_depth=3, learning_rate=0.05, random_state=42
        ),
        "XGBoost_Weighted": xgb.XGBClassifier(
            n_estimators=120, max_depth=3, learning_rate=0.05,
            scale_pos_weight=scale_pos_weight, eval_metric="logloss", random_state=42
        ),
    }


def _eval_one_class_cv(X: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
    """Evaluate One-Class Anomaly Detectors (trained ONLY on Positive samples) via 5-Fold CV."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    out: Dict[str, Any] = {}
    for name, base_clf in [
        ("IsolationForest_OneClass", IsolationForest(contamination=0.12, random_state=42)),
        ("OneClassSVM_RBF", OneClassSVM(nu=0.15, kernel="rbf", gamma="scale")),
    ]:
        scores = np.zeros(len(y), dtype=float)
        preds = np.zeros(len(y), dtype=int)
        for tr, te in skf.split(X, y):
            scaler = StandardScaler()
            X_tr_pos = scaler.fit_transform(X[tr][y[tr] == 1])
            X_te = scaler.transform(X[te])
            clf = base_clf.__class__(**base_clf.get_params())
            clf.fit(X_tr_pos)
            scores[te] = clf.decision_function(X_te)
            preds[te] = (clf.predict(X_te) == 1).astype(int)
        out[name] = {
            "cv5_roc_auc": round(float(roc_auc_score(y, scores)), 3),
            "cv5_balanced_acc": round(float(balanced_accuracy_score(y, preds)), 3),
            "cv5_tpr_sensitivity": round(float(recall_score(y, preds, pos_label=1)), 3),
            "cv5_tnr_specificity": round(float(recall_score(y, preds, pos_label=0)), 3),
        }
    return out


def _eval_cohort(X: np.ndarray, y: np.ndarray, feat_names: List[str], run_loocv: bool = True) -> Dict[str, Any]:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    loo = LeaveOneOut()
    spw = float(np.sum(y == 0) / max(1, np.sum(y == 1)))
    models = _build_models(spw)
    results: Dict[str, Any] = {}

    for mname, clf in models.items():
        # Stratified 5-Fold CV
        probs_5 = cross_val_predict(clf, X, y, cv=skf, method="predict_proba")[:, 1]
        preds_5 = (probs_5 >= 0.5).astype(int)
        entry: Dict[str, Any] = {
            "cv5_roc_auc": round(float(roc_auc_score(y, probs_5)), 3),
            "cv5_balanced_acc": round(float(balanced_accuracy_score(y, preds_5)), 3),
            "cv5_tpr_sensitivity": round(float(recall_score(y, preds_5, pos_label=1)), 3),
            "cv5_tnr_specificity": round(float(recall_score(y, preds_5, pos_label=0)), 3),
        }

        if run_loocv:
            probs_loo = cross_val_predict(clf, X, y, cv=loo, method="predict_proba")[:, 1]
            preds_loo = (probs_loo >= 0.5).astype(int)
            entry["loocv_roc_auc"] = round(float(roc_auc_score(y, probs_loo)), 3)
            entry["loocv_balanced_acc"] = round(float(balanced_accuracy_score(y, preds_loo)), 3)
            entry["loocv_tpr_sensitivity"] = round(float(recall_score(y, preds_loo, pos_label=1)), 3)
            entry["loocv_tnr_specificity"] = round(float(recall_score(y, preds_loo, pos_label=0)), 3)

        results[mname] = entry

    results.update(_eval_one_class_cv(X, y))

    # Feature importances from ExtraTrees & RandomForest fitted on full cohort (for interpretability)
    rf_full = RandomForestClassifier(n_estimators=200, max_depth=6, class_weight="balanced", random_state=42).fit(X, y)
    top_feats = sorted(
        [(feat_names[i], round(float(rf_full.feature_importances_[i]), 4)) for i in range(len(feat_names))],
        key=lambda x: x[1],
        reverse=True,
    )[:10]

    best_model = max(results, key=lambda k: (results[k]["cv5_balanced_acc"], results[k]["cv5_roc_auc"]))
    return {
        "n_samples": int(len(y)),
        "n_pos": int(np.sum(y == 1)),
        "n_neg": int(np.sum(y == 0)),
        "n_features": int(len(feat_names)),
        "best_model_by_cv5_bal_acc": best_model,
        "models": results,
        "rf_top10_features": top_feats,
    }


def main() -> None:
    p1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
    p2_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json")
    exp_db_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json")
    exp_7n_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_elements_sevenn_features.json")

    with open(p1_path, "r", encoding="utf-8") as f:
        p1_by_id = {r["id"]: r for r in json.load(f)["results"]}
    with open(p2_path, "r", encoding="utf-8") as f:
        p2_rows = json.load(f)["molecules_evaluated"]
    with open(exp_db_path, "r", encoding="utf-8") as f:
        exp_prec = json.load(f)["precursors_102"]["records"]
    with open(exp_7n_path, "r", encoding="utf-8") as f:
        exp_7n = json.load(f)["features"]

    exp_prec_by_id = {r["id"]: r for r in exp_prec}
    exp_prec_by_smi = {r["smiles"]: r for r in exp_prec}

    # Build 55-molecule audited + expanded 85-molecule physics cohort
    rows_55: List[Dict[str, Any]] = []
    for m in p2_rows:
        r1 = p1_by_id[m["id"]]
        new_desc = compute_new_chem_coord_descriptors(m["smiles"])
        vbur = float(m["vbur_percent"])
        bde = float(m["bde_min_kcal_mol"])
        dimer = float(m["delta_e_dimer_kcal_mol"])
        ecoh = float(m["e_coh_ev_per_atom"])
        rg = float(m["radius_of_gyration_ang"])
        vol = float(m["mol_volume_ang3"])
        feats = {
            "sascore": float(r1["sascore"]),
            "syba": float(r1["syba"]),
            "scscore": float(r1["scscore"]),
            "rascore": float(r1["rascore"]),
            "br_sascore": float(r1["br_sascore"]),
            "e_coh_ev_per_atom": ecoh,
            "bde_min_kcal_mol": bde,
            "vbur_percent": vbur,
            "delta_e_dimer_kcal_mol": dimer,
            "radius_of_gyration_ang": rg,
            "mol_volume_ang3": vol,
            "vbur_x_bde": vbur * bde,
            **new_desc,
        }
        rows_55.append({
            "id": m["id"],
            "name": m["name"],
            "label": 1 if m["target_type"] in ("Precursor", "Inhibitor") else 0,
            "features": feats,
        })

    rows_85 = list(rows_55)
    for ef in exp_7n:
        rec = exp_prec_by_id.get(ef["id"]) or exp_prec_by_smi.get(ef["smiles"])
        assert rec is not None
        new_desc = compute_new_chem_coord_descriptors(ef["smiles"])
        vbur = float(ef["vbur_percent"])
        bde = float(ef["bde_min_kcal_mol"])
        feats = {
            "sascore": float(rec["sascore"]),
            "syba": float(rec["syba"]),
            "scscore": float(rec["scscore"]),
            "rascore": float(rec["rascore"]),
            "br_sascore": float(rec["br_sascore_default"]),
            "e_coh_ev_per_atom": float(ef["e_coh_ev_per_atom"]),
            "bde_min_kcal_mol": bde,
            "vbur_percent": vbur,
            "delta_e_dimer_kcal_mol": float(ef["delta_e_dimer_kcal_mol"]),
            "radius_of_gyration_ang": float(ef["radius_of_gyration_ang"]),
            "mol_volume_ang3": float(ef["mol_volume_ang3"]),
            "vbur_x_bde": vbur * bde,
            **new_desc,
        }
        rows_85.append({
            "id": ef["id"],
            "name": ef["name"],
            "label": 1,
            "features": feats,
        })

    # Build 130-molecule cohort (Base SA + New Chem/Coord Descriptors, no 7net needed)
    baseline_smis = {m["smiles"] for m in p2_rows if m["target_type"] == "Precursor"}
    rows_130: List[Dict[str, Any]] = []
    for r in rows_55:
        f_no_7n = {
            k: v for k, v in r["features"].items()
            if k not in ("e_coh_ev_per_atom", "bde_min_kcal_mol", "vbur_percent", "delta_e_dimer_kcal_mol", "radius_of_gyration_ang", "mol_volume_ang3", "vbur_x_bde")
        }
        rows_130.append({"id": r["id"], "name": r["name"], "label": r["label"], "features": f_no_7n})
    for rec in exp_prec:
        if rec["smiles"] in baseline_smis:
            continue
        new_desc = compute_new_chem_coord_descriptors(rec["smiles"])
        rows_130.append({
            "id": rec["id"],
            "name": rec["name"],
            "label": 1,
            "features": {
                "sascore": float(rec["sascore"]),
                "syba": float(rec["syba"]),
                "scscore": float(rec["scscore"]),
                "rascore": float(rec["rascore"]),
                "br_sascore": float(rec["br_sascore_default"]),
                **new_desc,
            },
        })

    def _to_xy(rows: List[Dict[str, Any]]):
        names = list(rows[0]["features"].keys())
        X = np.array([[r["features"][k] for k in names] for r in rows], dtype=float)
        y = np.array([r["label"] for r in rows], dtype=int)
        return X, y, names

    X55, y55, fn55 = _to_xy(rows_55)
    X85, y85, fn85 = _to_xy(rows_85)
    X130, y130, fn130 = _to_xy(rows_130)

    # Evaluate Physics-Constrained Coordination Gate (TMC-SA_phys + SMARTS + Coordination/Volatility Gate)
    def _passes_coord_gate(smi: str) -> bool:
        d = compute_new_chem_coord_descriptors(smi)
        mol = Chem.MolFromSmiles(smi)
        has_tm = any(a.GetAtomicNum() in TM_ATOMIC_NUMS for a in mol.GetAtoms())
        flagged = (d["is_homoleptic_tm_sp3"] == 1.0) or (d["metal_oh_count"] >= 2.0) or (has_tm and d["aromatic_rings"] >= 5.0)
        return not flagged

    gate_55 = {}
    for tt, key in [("Precursor", "precursor_27"), ("Inhibitor", "inhibitor_18"), ("Negative Control", "negative_control_10")]:
        sub = [m for m in p2_rows if m["target_type"] == tt]
        old_p = sum(1 for m in sub if m["tuned_with_chem_pass"])
        new_p = sum(1 for m in sub if m["tuned_with_chem_pass"] and _passes_coord_gate(m["smiles"]))
        gate_55[key] = {"total": len(sub), "prior_tmc_sa_chem_pass": old_p, "with_coord_gate_pass": new_p}

    exp_75 = [r for r in exp_prec if r["smiles"] not in baseline_smis]
    exp_75_flagged = [r["name"] for r in exp_75 if not _passes_coord_gate(r["smiles"])]

    out = {
        "new_descriptors_added": [
            "metal_c_sigma_frac", "metal_hetero_frac", "is_homoleptic_mc", "is_homoleptic_tm_sp3",
            "beta_h_alkyl_count", "metal_oh_count", "cooh_count", "is_inert_alkane",
            "hbd", "hba", "tpsa", "aromatic_rings", "vbur_x_bde",
        ],
        "cohort_55_full": _eval_cohort(X55, y55, fn55, run_loocv=True),
        "cohort_85_full": _eval_cohort(X85, y85, fn85, run_loocv=True),
        "cohort_130_chem_coord": _eval_cohort(X130, y130, fn130, run_loocv=False),
        "physics_coord_gate_eval": {
            "cohort_55": gate_55,
            "negative_control_tnr_specificity": round((10 - gate_55["negative_control_10"]["with_coord_gate_pass"]) / 10.0, 3),
            "expanded_75_unflagged": 75 - len(exp_75_flagged),
            "expanded_75_flagged_names": exp_75_flagged,
        },
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_neg_control_model_search.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print("Saved ->", out_path)
    for cname in ("cohort_55_full", "cohort_85_full", "cohort_130_chem_coord"):
        c = out[cname]
        print(f"\n=== {cname} (N={c['n_samples']}, pos={c['n_pos']}, neg={c['n_neg']}, feats={c['n_features']}) ===")
        print("Best model:", c["best_model_by_cv5_bal_acc"])
        for mname, metrics in c["models"].items():
            print(f"  {mname:20s}: {metrics}")
        print("Top 5 features:", c["rf_top10_features"][:5])


if __name__ == "__main__":
    main()
