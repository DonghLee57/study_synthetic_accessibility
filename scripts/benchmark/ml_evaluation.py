#!/usr/bin/env python3
"""
Reusable ML Training, Cross-Validation, LOFMO Extrapolation & Descriptor Audits
-------------------------------------------------------------------------------
Shared evaluation utilities for JSON I/O, class-weighted supervised & one-class
classifiers, Leave-One-Failure-Mechanism-Out (LOFMO) extrapolation, Descriptor 4
thermodynamic orthogonality audits, and PCA / t-SNE latent space projections.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Tuple

import numpy as np
import xgboost as xgb
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import balanced_accuracy_score, recall_score, roc_auc_score, silhouette_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import OneClassSVM, SVC


def load_json(path: str) -> Any:
    """Load UTF-8 JSON payload from disk."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_json(path: str, payload: Any) -> None:
    """Write formatted UTF-8 JSON payload to disk, creating parent directories."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def rate_pct(num: int, den: int) -> float:
    """Compute percentage ``100 * num / den`` rounded to 1 decimal place."""
    return round(100.0 * num / den, 1) if den else 0.0


def build_ml_classifiers(scale_pos_weight: float, random_state: int = 42) -> Dict[str, Any]:
    """Instantiate class-balanced supervised classifiers."""
    return {
        "LogisticRegression": make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight="balanced", random_state=random_state),
        ),
        "SVC_RBF": make_pipeline(
            StandardScaler(),
            SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=random_state),
        ),
        "RandomForest": RandomForestClassifier(
            n_estimators=200, class_weight="balanced", random_state=random_state
        ),
        "XGBoost_Weighted": xgb.XGBClassifier(
            n_estimators=150,
            max_depth=3,
            learning_rate=0.05,
            scale_pos_weight=scale_pos_weight,
            eval_metric="logloss",
            random_state=random_state,
        ),
    }


def evaluate_binary_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_score: np.ndarray) -> Dict[str, float]:
    """Compute ROC AUC, balanced accuracy, TPR (Recall), and TNR (Specificity)."""
    return {
        "cv5_roc_auc": round(float(roc_auc_score(y_true, y_score)), 4),
        "cv5_balanced_acc": round(float(balanced_accuracy_score(y_true, y_pred)), 4),
        "cv5_tpr_sensitivity": round(float(recall_score(y_true, y_pred, pos_label=1)), 4),
        "cv5_tnr_specificity": round(float(recall_score(y_true, y_pred, pos_label=0)), 4),
    }


def evaluate_cv5_and_lofmo(
    X: np.ndarray,
    y: np.ndarray,
    row_mechs: List[str],
    random_state: int = 42,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Benchmark supervised models + OneClassSVM_RBF under Stratified 5-Fold CV and
    Leave-One-Failure-Mechanism-Out (LOFMO) extrapolation.
    """
    cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    scale_pos = (n_neg / n_pos) if n_pos else 1.0

    ml_cv5: Dict[str, Any] = {}
    for mname, clf in build_ml_classifiers(scale_pos, random_state=random_state).items():
        prob = cross_val_predict(clf, X, y, cv=cv5, method="predict_proba")[:, 1]
        pred = (prob >= 0.5).astype(int)
        ml_cv5[mname] = evaluate_binary_metrics(y, pred, prob)

    oc_scores = np.zeros(len(y), dtype=float)
    oc_preds = np.zeros(len(y), dtype=int)
    for tr_idx, te_idx in cv5.split(X, y):
        pipe = make_pipeline(StandardScaler(), OneClassSVM(kernel="rbf", nu=0.08, gamma="scale"))
        pipe.fit(X[tr_idx][y[tr_idx] == 1])
        oc_scores[te_idx] = pipe.decision_function(X[te_idx])
        oc_preds[te_idx] = (pipe.predict(X[te_idx]) == 1).astype(int)
    ml_cv5["OneClassSVM_RBF"] = evaluate_binary_metrics(y, oc_preds, oc_scores)

    unique_mechs = sorted({m for m in row_mechs if m != "pos"})
    lofmo: Dict[str, Any] = {}
    for mname, fn in [
        ("LogisticRegression", lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced", random_state=random_state))),
        ("SVC_RBF", lambda: make_pipeline(StandardScaler(), SVC(kernel="rbf", class_weight="balanced", random_state=random_state))),
        ("RandomForest", lambda: RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=random_state)),
    ]:
        per_mech = {}
        for hold in unique_mechs:
            tr = np.array([m != hold for m in row_mechs])
            te = np.array([m == hold for m in row_mechs])
            clf = fn()
            clf.fit(X[tr], y[tr])
            per_mech[hold] = round(float(np.mean(clf.predict(X[te]) == 0)), 4)
        lofmo[mname] = {
            "mean_unseen_tnr": round(float(np.mean(list(per_mech.values()))), 4),
            "per_unseen_mechanism_tnr": per_mech,
        }

    oc_full = make_pipeline(StandardScaler(), OneClassSVM(kernel="rbf", nu=0.08, gamma="scale")).fit(X[y == 1])
    oc_mech = {
        hold: round(float(np.mean(oc_full.predict(X[np.array([m == hold for m in row_mechs])]) == -1)), 4)
        for hold in unique_mechs
    }
    lofmo["OneClassSVM_RBF"] = {
        "mean_unseen_tnr": round(float(np.mean(list(oc_mech.values()))), 4),
        "per_unseen_mechanism_tnr": oc_mech,
    }
    return ml_cv5, lofmo


def compute_latent_projections(X: np.ndarray, y: np.ndarray, random_state: int = 42) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Run standardized 2D PCA & t-SNE projections and measure linear/kernel separability."""
    X_std = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=random_state)
    X_pca = pca.fit_transform(X_std)

    tsne = TSNE(n_components=2, perplexity=25, random_state=random_state, init="pca", learning_rate="auto", max_iter=1000)
    X_tsne = tsne.fit_transform(X_std)

    sil_pca = float(silhouette_score(X_pca, y))
    sil_tsne = float(silhouette_score(X_tsne, y))
    linear_sep_pca = float(
        make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", random_state=random_state))
        .fit(X_pca, y)
        .score(X_pca, y)
    )
    pca_lr_bal = float(
        balanced_accuracy_score(
            y, LogisticRegression(class_weight="balanced", random_state=random_state).fit(X_pca, y).predict(X_pca)
        )
    )
    full_svc_bal = float(
        balanced_accuracy_score(
            y, SVC(kernel="rbf", class_weight="balanced", random_state=random_state).fit(X_std, y).predict(X_std)
        )
    )
    stats = {
        "pca_explained_variance_ratio": [round(float(v), 4) for v in pca.explained_variance_ratio_],
        "pca_silhouette": round(sil_pca, 4),
        "tsne_silhouette": round(sil_tsne, 4),
        "pca_2d_balanced_acc": round(pca_lr_bal, 4),
        "full_space_svc_balanced_acc": round(full_svc_bal, 4),
        "pca_linear_acc": round(linear_sep_pca, 4),
    }
    return X_pca, X_tsne, stats


def evaluate_descriptor4_audit(records: List[Dict[str, Any]], random_state: int = 42) -> Dict[str, Any]:
    """
    Audit orthogonality and discriminative power (AUROC & marginal LOOCV contribution)
    of Descriptor 4 (``delta_e_form_per_bond_kcal_mol``) vs. the legacy algebraic proxy
    across MLIP-relaxed subsets (89 total, 85 core+expanded precursors, 55 audited core).
    """
    s89_recs = [r for r in records if r["delta_e_form_per_bond_kcal_mol"] is not None]
    s85_recs = [r for r in s89_recs if not (r["group"] == "Negative Control (40)" and r["origin"] != "audited_baseline_10")]
    s55_recs = [r for r in s85_recs if r["group"] in ("Audited Precursor (27)", "AS-ALD Inhibitor (18)") or r["origin"] == "audited_baseline_10"]

    def _audit_sub(sub: List[Dict[str, Any]]) -> Dict[str, Any]:
        y_sub = np.array([r["label"] for r in sub], dtype=int)
        stats: Dict[str, Any] = {"n_total": len(sub), "n_pos": int(np.sum(y_sub == 1)), "n_neg": int(np.sum(y_sub == 0))}
        for key, field, higher_is_pos in [
            ("legacy_proxy", "delta_e_complex_legacy_kcal_mol", False),
            ("delta_e_form_per_bond", "delta_e_form_per_bond_kcal_mol", False),
            ("delta_g_dimer_373", "delta_g_dimer_373_kcal_mol", True),
            ("bde_min", "bde_min_kcal_mol", True),
            ("e_coh", "e_coh_ev_per_atom", False),
        ]:
            vals = np.array([float(r[field]) for r in sub], dtype=float)
            auc = roc_auc_score(y_sub, vals if higher_is_pos else -vals)
            rb = float(np.corrcoef(vals, [float(r["bde_min_kcal_mol"]) for r in sub])[0, 1])
            rv = float(np.corrcoef(vals, [float(r["vbur_r35_percent"]) for r in sub])[0, 1])
            rd = float(np.corrcoef(vals, [float(r["delta_e_dimer_kcal_mol"]) for r in sub])[0, 1])
            stats[key] = {
                "auroc": round(float(auc), 3),
                "pos_mean": round(float(vals[y_sub == 1].mean()), 3),
                "pos_std": round(float(vals[y_sub == 1].std()), 3),
                "neg_mean": round(float(vals[y_sub == 0].mean()), 3),
                "neg_std": round(float(vals[y_sub == 0].std()), 3),
                "corr_bde_vbur_ddim": [round(rb, 2), round(rv, 2), round(rd, 2)],
            }
        return stats

    def _loocv_sub(sub_recs: List[Dict[str, Any]]) -> Dict[str, Any]:
        y_s = np.array([r["label"] for r in sub_recs], dtype=int)
        loo = LeaveOneOut()
        res: Dict[str, Any] = {}
        for label, cols in [
            ("D1_D3_only", ["bde_min_kcal_mol", "e_coh_ev_per_atom", "vbur_percent", "delta_e_dimer_kcal_mol"]),
            ("D1_D3_plus_legacy_proxy", ["bde_min_kcal_mol", "e_coh_ev_per_atom", "vbur_percent", "delta_e_dimer_kcal_mol", "delta_e_complex_legacy_kcal_mol"]),
            ("D1_D3_plus_delta_e_form_per_bond", ["bde_min_kcal_mol", "e_coh_ev_per_atom", "vbur_percent", "delta_e_dimer_kcal_mol", "delta_e_form_per_bond_kcal_mol"]),
            ("D1_D3_plus_dEform_plus_dGdimer", ["bde_min_kcal_mol", "e_coh_ev_per_atom", "vbur_percent", "delta_e_dimer_kcal_mol", "delta_e_form_per_bond_kcal_mol", "delta_g_dimer_373_kcal_mol"]),
        ]:
            X_s = np.array([[float(r[c]) for c in cols] for r in sub_recs], dtype=float)
            pipe = make_pipeline(
                StandardScaler(),
                LogisticRegression(class_weight="balanced", max_iter=1000, random_state=random_state),
            )
            preds = cross_val_predict(pipe, X_s, y_s, cv=loo)
            probs = cross_val_predict(pipe, X_s, y_s, cv=loo, method="predict_proba")[:, 1]
            res[label] = {
                "balanced_acc_pct": round(100.0 * float(balanced_accuracy_score(y_s, preds)), 1),
                "tpr_pct": round(100.0 * float(recall_score(y_s, preds, pos_label=1)), 1),
                "tnr_pct": round(100.0 * float(recall_score(y_s, preds, pos_label=0)), 1),
                "auroc": round(float(roc_auc_score(y_s, probs)), 3),
            }
        return res

    return {
        "audited_55": _audit_sub(s55_recs),
        "full_85": _audit_sub(s85_recs),
        "full_89": _audit_sub(s89_recs),
        "loocv_55_marginal_contribution": _loocv_sub(s55_recs),
        "loocv_89_marginal_contribution": _loocv_sub(s89_recs),
    }
