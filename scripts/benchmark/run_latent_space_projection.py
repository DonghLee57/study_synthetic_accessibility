#!/usr/bin/env python3
"""
Latent Space Redefinition & Separability Analysis (PCA & t-SNE)
----------------------------------------------------------------
Evaluates whether redefining the descriptor space via PCA and t-SNE
improves separation between synthetically accessible semiconductor molecules
(Precursors & AS-ALD Inhibitors) and Negative Controls compared to any
single descriptor alone.

Cohorts:
  1. 130-Molecule Expanded Cohort (55 Audited + 75 Expanded Ir/Rh/As/Ga/In/Ru)
     using 12 Base SA + RDKit Physicochemical descriptors.
  2. 85-Molecule SevenNet-0 Physics Cohort (55 Audited + 30 Expanded Ir/Rh/As/Ga/In/Ru)
     using 11 Base SA + SevenNet-0 3D Physical descriptors.

Outputs:
  - database/benchmark_latent_space_projections.json
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, silhouette_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _rdkit_phys_feats(smiles: str) -> Dict[str, float]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES: {smiles}"
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "heavy_atoms": float(mol.GetNumHeavyAtoms()),
        "ring_count": float(rdMolDescriptors.CalcNumRings(mol)),
        "rot_bonds": float(Lipinski.NumRotatableBonds(mol)),
    }


def _evaluate_separability(
    X: np.ndarray, y: np.ndarray, feature_names: List[str], pca_coords: np.ndarray, tsne_coords: np.ndarray
) -> Dict[str, Any]:
    """Compute Stratified 5-Fold CV ROC-AUC & Balanced Accuracy for single features vs PCA/t-SNE."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    single_results: Dict[str, Dict[str, float]] = {}
    for idx, fname in enumerate(feature_names):
        col = X[:, [idx]]
        clf = LogisticRegression(class_weight="balanced", random_state=42, max_iter=1000)
        probs = cross_val_predict(clf, col, y, cv=skf, method="predict_proba")[:, 1]
        preds = (probs >= 0.5).astype(int)
        single_results[fname] = {
            "cv5_roc_auc": round(float(roc_auc_score(y, probs)), 3),
            "cv5_balanced_acc": round(float(balanced_accuracy_score(y, preds)), 3),
        }

    # Best single descriptor by ROC-AUC
    best_single_name = max(single_results, key=lambda k: single_results[k]["cv5_roc_auc"])
    best_single = single_results[best_single_name]

    # PCA 2D & Full-Space Logistic + k-NN
    def _eval_2d(coords: np.ndarray) -> Dict[str, float]:
        lr = LogisticRegression(class_weight="balanced", random_state=42, max_iter=1000)
        probs_lr = cross_val_predict(lr, coords, y, cv=skf, method="predict_proba")[:, 1]
        preds_lr = (probs_lr >= 0.5).astype(int)

        knn = KNeighborsClassifier(n_neighbors=5)
        probs_knn = cross_val_predict(knn, coords, y, cv=skf, method="predict_proba")[:, 1]
        preds_knn = (probs_knn >= 0.5).astype(int)

        return {
            "lr_cv5_roc_auc": round(float(roc_auc_score(y, probs_lr)), 3),
            "lr_cv5_balanced_acc": round(float(balanced_accuracy_score(y, preds_lr)), 3),
            "knn5_cv5_roc_auc": round(float(roc_auc_score(y, probs_knn)), 3),
            "knn5_cv5_balanced_acc": round(float(balanced_accuracy_score(y, preds_knn)), 3),
            "silhouette_score": round(float(silhouette_score(coords, y)), 3),
        }

    lr_full = LogisticRegression(class_weight="balanced", random_state=42, max_iter=1000)
    X_scaled = StandardScaler().fit_transform(X)
    probs_full = cross_val_predict(lr_full, X_scaled, y, cv=skf, method="predict_proba")[:, 1]
    preds_full = (probs_full >= 0.5).astype(int)

    return {
        "single_descriptors": single_results,
        "best_single_descriptor": {"name": best_single_name, **best_single},
        "pca_2d": _eval_2d(pca_coords),
        "tsne_2d": _eval_2d(tsne_coords),
        "full_scaled_space": {
            "lr_cv5_roc_auc": round(float(roc_auc_score(y, probs_full)), 3),
            "lr_cv5_balanced_acc": round(float(balanced_accuracy_score(y, preds_full)), 3),
        },
    }


def run_projections() -> None:
    p1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
    p2_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json")
    exp_db_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json")
    exp_7n_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_elements_sevenn_features.json")

    with open(p1_path, "r", encoding="utf-8") as f:
        p1_rows = json.load(f)["results"]
    with open(p2_path, "r", encoding="utf-8") as f:
        p2_rows = json.load(f)["molecules_evaluated"]
    with open(exp_db_path, "r", encoding="utf-8") as f:
        exp_db = json.load(f)
    with open(exp_7n_path, "r", encoding="utf-8") as f:
        exp_7n_feats = json.load(f)["features"]

    p1_by_id = {r["id"]: r for r in p1_rows}
    p2_by_id = {r["id"]: r for r in p2_rows}
    exp_prec_by_id = {r["id"]: r for r in exp_db["precursors_102"]["records"]}
    exp_prec_by_smi = {r["smiles"]: r for r in exp_db["precursors_102"]["records"]}

    # -------------------------------------------------------------------------
    # Cohort 1: 130-Molecule Expanded Cohort (55 Audited + 75 New Precursors)
    # -------------------------------------------------------------------------
    cohort_130: List[Dict[str, Any]] = []
    seen_ids = set()

    # Add 55 audited molecules first
    for m in p2_rows:
        mid = m["id"]
        seen_ids.add(mid)
        r1 = p1_by_id[mid]
        phys = _rdkit_phys_feats(m["smiles"])
        tt = m["target_type"]
        group = (
            "Baseline Precursor (27)"
            if tt == "Precursor"
            else ("AS-ALD Inhibitor (18)" if tt == "Inhibitor" else "Negative Control (10)")
        )
        cohort_130.append({
            "id": mid,
            "name": m["name"],
            "group": group,
            "element": m.get("metal", ""),
            "label": 1 if tt in ("Precursor", "Inhibitor") else 0,
            "features": {
                "sascore": float(r1["sascore"]),
                "syba": float(r1["syba"]),
                "scscore": float(r1["scscore"]),
                "rascore": float(r1["rascore"]),
                "br_sascore": float(r1["br_sascore"]),
                **phys,
            },
        })

    # Add 75 expanded Ir/Rh/As/Ga/In/Ru precursors (non-overlapping with 27 baseline)
    baseline_smis = {m["smiles"] for m in p2_rows if m["target_type"] == "Precursor"}
    for rec in exp_db["precursors_102"]["records"]:
        if rec["smiles"] in baseline_smis:
            continue
        phys = _rdkit_phys_feats(rec["smiles"])
        cohort_130.append({
            "id": rec["id"],
            "name": rec["name"],
            "group": "Expanded Precursor (75)",
            "element": rec.get("center_element", ""),
            "label": 1,
            "features": {
                "sascore": float(rec["sascore"]),
                "syba": float(rec["syba"]),
                "scscore": float(rec["scscore"]),
                "rascore": float(rec["rascore"]),
                "br_sascore": float(rec["br_sascore_default"]),
                **phys,
            },
        })

    feat_names_130 = list(cohort_130[0]["features"].keys())
    X_130 = np.array([[row["features"][k] for k in feat_names_130] for row in cohort_130], dtype=float)
    y_130 = np.array([row["label"] for row in cohort_130], dtype=int)

    X_130_scaled = StandardScaler().fit_transform(X_130)
    pca_130 = PCA(n_components=2, random_state=42)
    pca_130_coords = pca_130.fit_transform(X_130_scaled)
    tsne_130 = TSNE(n_components=2, perplexity=20, init="pca", learning_rate="auto", random_state=42)
    tsne_130_coords = tsne_130.fit_transform(X_130_scaled)

    sep_130 = _evaluate_separability(X_130, y_130, feat_names_130, pca_130_coords, tsne_130_coords)
    points_130 = []
    for i, row in enumerate(cohort_130):
        points_130.append({
            "id": row["id"],
            "name": row["name"],
            "group": row["group"],
            "element": row["element"],
            "label": int(row["label"]),
            "pc1": round(float(pca_130_coords[i, 0]), 4),
            "pc2": round(float(pca_130_coords[i, 1]), 4),
            "tsne1": round(float(tsne_130_coords[i, 0]), 4),
            "tsne2": round(float(tsne_130_coords[i, 1]), 4),
        })

    # -------------------------------------------------------------------------
    # Cohort 2: 85-Molecule SevenNet-0 Physics Cohort (55 Audited + 30 Expanded)
    # -------------------------------------------------------------------------
    cohort_85: List[Dict[str, Any]] = []
    for m in p2_rows:
        mid = m["id"]
        r1 = p1_by_id[mid]
        tt = m["target_type"]
        group = (
            "Baseline Precursor (27)"
            if tt == "Precursor"
            else ("AS-ALD Inhibitor (18)" if tt == "Inhibitor" else "Negative Control (10)")
        )
        cohort_85.append({
            "id": mid,
            "name": m["name"],
            "group": group,
            "element": m.get("metal", ""),
            "label": 1 if tt in ("Precursor", "Inhibitor") else 0,
            "features": {
                "sascore": float(r1["sascore"]),
                "syba": float(r1["syba"]),
                "scscore": float(r1["scscore"]),
                "rascore": float(r1["rascore"]),
                "br_sascore": float(r1["br_sascore"]),
                "e_coh_ev_per_atom": float(m["e_coh_ev_per_atom"]),
                "bde_min_kcal_mol": float(m["bde_min_kcal_mol"]),
                "vbur_percent": float(m["vbur_percent"]),
                "delta_e_dimer_kcal_mol": float(m["delta_e_dimer_kcal_mol"]),
                "radius_of_gyration_ang": float(m["radius_of_gyration_ang"]),
                "mol_volume_ang3": float(m["mol_volume_ang3"]),
            },
        })

    for ef in exp_7n_feats:
        rec = exp_prec_by_id.get(ef["id"]) or exp_prec_by_smi.get(ef["smiles"])
        assert rec is not None, f"Missing SA record for {ef['id']}"
        cohort_85.append({
            "id": ef["id"],
            "name": ef["name"],
            "group": "Expanded Precursor (30)",
            "element": ef.get("center_element_expected", ef.get("center_atom", "")),
            "label": 1,
            "features": {
                "sascore": float(rec["sascore"]),
                "syba": float(rec["syba"]),
                "scscore": float(rec["scscore"]),
                "rascore": float(rec["rascore"]),
                "br_sascore": float(rec["br_sascore_default"]),
                "e_coh_ev_per_atom": float(ef["e_coh_ev_per_atom"]),
                "bde_min_kcal_mol": float(ef["bde_min_kcal_mol"]),
                "vbur_percent": float(ef["vbur_percent"]),
                "delta_e_dimer_kcal_mol": float(ef["delta_e_dimer_kcal_mol"]),
                "radius_of_gyration_ang": float(ef["radius_of_gyration_ang"]),
                "mol_volume_ang3": float(ef["mol_volume_ang3"]),
            },
        })

    feat_names_85 = list(cohort_85[0]["features"].keys())
    X_85 = np.array([[row["features"][k] for k in feat_names_85] for row in cohort_85], dtype=float)
    y_85 = np.array([row["label"] for row in cohort_85], dtype=int)

    X_85_scaled = StandardScaler().fit_transform(X_85)
    pca_85 = PCA(n_components=2, random_state=42)
    pca_85_coords = pca_85.fit_transform(X_85_scaled)
    tsne_85 = TSNE(n_components=2, perplexity=18, init="pca", learning_rate="auto", random_state=42)
    tsne_85_coords = tsne_85.fit_transform(X_85_scaled)

    sep_85 = _evaluate_separability(X_85, y_85, feat_names_85, pca_85_coords, tsne_85_coords)
    points_85 = []
    for i, row in enumerate(cohort_85):
        points_85.append({
            "id": row["id"],
            "name": row["name"],
            "group": row["group"],
            "element": row["element"],
            "label": int(row["label"]),
            "pc1": round(float(pca_85_coords[i, 0]), 4),
            "pc2": round(float(pca_85_coords[i, 1]), 4),
            "tsne1": round(float(tsne_85_coords[i, 0]), 4),
            "tsne2": round(float(tsne_85_coords[i, 1]), 4),
        })

    def _top_loadings(pca_obj: PCA, names: List[str]) -> Dict[str, List[Tuple[str, float]]]:
        out = {}
        for comp_idx, comp_name in enumerate(("PC1", "PC2")):
            pairs = sorted(
                [(names[j], round(float(pca_obj.components_[comp_idx, j]), 3)) for j in range(len(names))],
                key=lambda x: abs(x[1]),
                reverse=True,
            )
            out[comp_name] = pairs[:4]
        return out

    output = {
        "cohort_130_sa_rdkit": {
            "total": len(cohort_130),
            "features": feat_names_130,
            "pca_explained_variance_ratio": [round(float(x), 4) for x in pca_130.explained_variance_ratio_],
            "pca_top_loadings": _top_loadings(pca_130, feat_names_130),
            "separability": sep_130,
            "points": points_130,
        },
        "cohort_85_sa_sevenn_physics": {
            "total": len(cohort_85),
            "features": feat_names_85,
            "pca_explained_variance_ratio": [round(float(x), 4) for x in pca_85.explained_variance_ratio_],
            "pca_top_loadings": _top_loadings(pca_85, feat_names_85),
            "separability": sep_85,
            "points": points_85,
        },
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_latent_space_projections.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print("Saved ->", out_path)
    print("Cohort 130 PCA var:", output["cohort_130_sa_rdkit"]["pca_explained_variance_ratio"])
    print("Cohort 130 Separability:", json.dumps({
        k: sep_130[k] for k in ("best_single_descriptor", "pca_2d", "tsne_2d", "full_scaled_space")
    }, indent=2))
    print("Cohort 85 PCA var:", output["cohort_85_sa_sevenn_physics"]["pca_explained_variance_ratio"])
    print("Cohort 85 Separability:", json.dumps({
        k: sep_85[k] for k in ("best_single_descriptor", "pca_2d", "tsne_2d", "full_scaled_space")
    }, indent=2))


if __name__ == "__main__":
    run_projections()
