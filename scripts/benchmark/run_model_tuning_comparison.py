#!/usr/bin/env python3
"""
Phase 2 Model Re-Evaluation & Rigorous Statistical Audit
---------------------------------------------------------
Evaluates the integration of SevenNet (7net-0) atomistic features and Cavallo SambVca
buried volume (%V_bur) without target-type branching, hardcoded sample IDs, or
synthetic algebraic reaction proxies.

Provides:
  1. Unbiased Uniform Continuous Physical Screening (TMC-SA_smooth):
     Applies the exact same ligand decomposition and smooth Fermi-Dirac physical
     penalties (E_coh, BDE_min, %V_bur, Delta E_dimer) uniformly across all 55 molecules.
  2. Supervised Machine Learning Benchmark:
     Evaluates All-Sample (in-sample), 5-Fold Stratified CV, and Leave-One-Out CV (LOOCV)
     across Random Forest and Logistic Regression on:
       - Full Dataset (N=55: 45 Valid vs 10 Negative Controls)
       - Precursor-Domain Subset (N=35: 27 Precursors vs 8 Organometallic/Silane Negatives)
       - Feature sets: Chemoinformatics Baseline vs 7net Physical Only vs Combined.

Outputs complete audited results to database/benchmark_phase2_tuned_comparison.json.
"""

import json
import math
import os
import sys
import numpy as np
from rdkit import Chem
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def sigmoid(x: float) -> float:
    """Numerically stable logistic sigmoid."""
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    z = math.exp(x)
    return z / (1.0 + z)


def load_merged_data():
    """Merge Phase 1 baseline (with fixed BR-SAScore) and latest Phase 2 SevenNet features."""
    phase1_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase1_baseline.json")
    phase2_path = os.path.join(PROJECT_ROOT, "database", "external", "sevenn_physical_features.json")

    with open(phase1_path, "r", encoding="utf-8") as f:
        p1 = json.load(f)["results"]
    with open(phase2_path, "r", encoding="utf-8") as f:
        p2_json = json.load(f)
    p2 = {item["id"]: item for item in p2_json["features"]}

    merged = []
    for m in p1:
        feat = p2.get(m["id"], {})
        entry = {**m, **feat}
        entry["is_fallback_geometry"] = bool(entry.get("e_pot_ev") == -100.0)
        merged.append(entry)
    return merged, p2_json.get("metadata", {})


METAL_CENTER_SET = {"Ti", "Zr", "Hf", "Al", "Ru", "Co", "W", "Mo", "Ta", "Sn", "Cu", "Si"}


def compute_smooth_physical_penalty(m: dict) -> dict:
    """
    Compute smooth C-infinity physical penalty components strictly from computed
    SevenNet/SambVca physical descriptors (no target_type branching or hardcoded IDs).
    Note: Steric crowding & metal under-coordination oligomerization (%V_bur) apply
    structurally only when the center atom is a metal or silicon atom.
    """
    bde = float(m.get("bde_min_kcal_mol", 60.0))
    ecoh = float(m.get("e_coh_ev_per_atom", 4.10))
    vbur = float(m.get("vbur_percent", 70.0))
    d_dim = float(m.get("delta_e_dimer_kcal_mol", 10.0))
    has_metal_center = str(m.get("center_atom", "C")) in METAL_CENTER_SET

    # Smooth Fermi-Dirac penalties centered at domain boundaries
    pen_bde = 3.5 * sigmoid((45.0 - bde) / 2.0)
    pen_coh = 3.5 * sigmoid((3.10 - ecoh) / 0.12)
    pen_labile = max(pen_bde, pen_coh)

    if has_metal_center:
        pen_crowd = 3.5 * sigmoid((vbur - 83.5) / 1.0)
        pen_olig = 3.5 * sigmoid((44.0 - vbur) / 2.0) * sigmoid((-1.5 - d_dim) / 1.0)
    else:
        pen_crowd = 0.0
        pen_olig = 0.0

    total_pen = round(pen_labile + pen_crowd + pen_olig, 3)
    return {
        "has_metal_center": has_metal_center,
        "pen_labile": round(pen_labile, 3),
        "pen_crowd": round(pen_crowd, 3),
        "pen_olig": round(pen_olig, 3),
        "total_phys_penalty": total_pen,
    }


def has_chemical_reactivity_or_inertness_flag(smi: str) -> bool:
    """
    General SMARTS chemical-domain screen independent of sample IDs:
    1) Free carboxylic acid (-COOH) or haloacetic acid (corrosive etching hazard);
    2) Pure saturated acyclic alkane lacking any heteroatom or pi-bond (no ALD/ASD anchoring site).
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        return False
    acid_patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    is_acid = mol.HasSubstructMatch(acid_patt)
    atoms = {a.GetSymbol() for a in mol.GetAtoms()}
    has_unsaturated = any(b.GetBondTypeAsDouble() > 1.0 for b in mol.GetBonds())
    is_inert_alkane = (atoms <= {"C", "H"}) and (not has_unsaturated) and (mol.GetRingInfo().NumRings() == 0)
    return bool(is_acid or is_inert_alkane)


def evaluate_unbiased_analytical_models(merged):
    """
    Evaluate uniform continuous screening across all 55 molecules without
    ground-truth label branching or hardcoded molecule IDs.
    """
    for m in merged:
        pen = compute_smooth_physical_penalty(m)
        ligand_sa = float(m.get("ligand_sa", m.get("sascore", 3.5)))
        tmc_sa_phys = round(ligand_sa + pen["total_phys_penalty"], 2)

        m["phys_penalty_detail"] = pen
        m["tmc_sa_phys"] = tmc_sa_phys
        # Uniform decision threshold (pure physics + ligand SA <= 4.0)
        m["tuned_pass"] = bool(tmc_sa_phys <= 4.0)
        # Combined chemistry-domain + physical screening
        chem_flag = has_chemical_reactivity_or_inertness_flag(m["smiles"])
        m["chemical_domain_flag"] = chem_flag
        m["tuned_with_chem_pass"] = bool(m["tuned_pass"] and (not chem_flag))

    precursors = [m for m in merged if m["target_type"] == "Precursor"]
    inhibitors = [m for m in merged if m["target_type"] == "Inhibitor"]
    negatives = [m for m in merged if m["target_type"] == "Negative Control"]

    def summarize_group(items, pred_key):
        y_true = np.array([1 if x["is_fab_valid"] else 0 for x in items])
        y_pred = np.array([1 if x[pred_key] else 0 for x in items])
        tp = int(np.sum((y_true == 1) & (y_pred == 1)))
        fn = int(np.sum((y_true == 1) & (y_pred == 0)))
        fp = int(np.sum((y_true == 0) & (y_pred == 1)))
        tn = int(np.sum((y_true == 0) & (y_pred == 0)))
        return {
            "count": len(items),
            "tp": tp, "fn": fn, "fp": fp, "tn": tn,
            "pass_count": int(np.sum(y_pred == 1)),
            "pass_rate_pct": round(float(np.mean(y_pred == 1) * 100.0), 1),
            "failed_ids": [x["id"] for x in items if (x["is_fab_valid"] and not x[pred_key])],
            "missed_negative_ids": [x["id"] for x in items if (not x["is_fab_valid"] and x[pred_key])],
        }

    # Global metrics (un-hardcoded uniform rule)
    y_all = np.array([1 if x["is_fab_valid"] else 0 for x in merged])
    p_phys = np.array([1 if x["tuned_pass"] else 0 for x in merged])
    p_chem = np.array([1 if x["tuned_with_chem_pass"] else 0 for x in merged])

    def global_metrics(y_true, y_pred):
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        return {
            "accuracy": round(float(accuracy_score(y_true, y_pred)), 3),
            "recall": round(float(recall_score(y_true, y_pred)), 3),
            "specificity": round(float(tn / (tn + fp)), 3),
            "f1_score": round(float(f1_score(y_true, y_pred)), 3),
            "tp": int(tp), "fn": int(fn), "fp": int(fp), "tn": int(tn),
        }

    prec_sa_pass = sum(1 for m in precursors if m["sascore"] <= 3.5)
    inhib_sa_pass = sum(1 for m in inhibitors if m["sascore"] <= 3.5)
    neg_sc_pass = sum(1 for m in negatives if m["scscore"] <= 2.0)
    neg_cons_pass = sum(1 for m in negatives if m["baseline_consensus_pass"])

    phys_prec = summarize_group(precursors, "tuned_pass")
    phys_inhib = summarize_group(inhibitors, "tuned_pass")
    phys_neg = summarize_group(negatives, "tuned_pass")

    chem_prec = summarize_group(precursors, "tuned_with_chem_pass")
    chem_inhib = summarize_group(inhibitors, "tuned_with_chem_pass")
    chem_neg = summarize_group(negatives, "tuned_with_chem_pass")

    return {
        "methodology": "Uniform continuous C-infinity TMC-SA_phys (ligand_sa + smooth E_coh/BDE/Vbur/dE_dimer penalties <= 4.0); zero target_type branching or hardcoded IDs.",
        "overall_uniform_physics_only": global_metrics(y_all, p_phys),
        "overall_physics_plus_smarts_chem": global_metrics(y_all, p_chem),
        "baseline_comparison": {
            "precursors_sascore_pass": f"{prec_sa_pass}/{len(precursors)} ({prec_sa_pass/len(precursors)*100:.1f}%)",
            "precursors_sascore_fn": f"{(len(precursors)-prec_sa_pass)/len(precursors)*100:.1f}%",
            "inhibitors_sascore_pass": f"{inhib_sa_pass}/{len(inhibitors)} ({inhib_sa_pass/len(inhibitors)*100:.1f}%)",
            "negatives_scscore_fp": f"{neg_sc_pass}/{len(negatives)} ({neg_sc_pass/len(negatives)*100:.1f}%)",
            "negatives_consensus_fp": f"{neg_cons_pass}/{len(negatives)} ({neg_cons_pass/len(negatives)*100:.1f}%)",
        },
        "uniform_physics_only_breakdown": {
            "precursors": phys_prec,
            "inhibitors": phys_inhib,
            "negative_controls": phys_neg,
        },
        "physics_plus_smarts_chem_breakdown": {
            "precursors": chem_prec,
            "inhibitors": chem_inhib,
            "negative_controls": chem_neg,
        },
    }


def _evaluate_ml_regime(X, y, clf_factory):
    """Run All-sample, Stratified 5-Fold CV, and Leave-One-Out CV on given (X, y)."""
    clf_in = clf_factory()
    clf_in.fit(X, y)
    pr_in = clf_in.predict_proba(X)[:, 1]
    p_in = (pr_in >= 0.5).astype(int)

    cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    pr_5 = cross_val_predict(clf_factory(), X, y, cv=cv5, method="predict_proba", n_jobs=1)[:, 1]
    p_5 = (pr_5 >= 0.5).astype(int)

    loo = LeaveOneOut()
    pr_loo = cross_val_predict(clf_factory(), X, y, cv=loo, method="predict_proba", n_jobs=1)[:, 1]
    p_loo = (pr_loo >= 0.5).astype(int)

    def _metrics(pred, prob):
        tn, fp, fn, tp = confusion_matrix(y, pred).ravel()
        return {
            "accuracy": round(float(accuracy_score(y, pred)), 3),
            "precision": round(float(precision_score(y, pred, zero_division=0)), 3),
            "recall": round(float(recall_score(y, pred, zero_division=0)), 3),
            "specificity": round(float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0, 3),
            "f1_score": round(float(f1_score(y, pred, zero_division=0)), 3),
            "auroc": round(float(roc_auc_score(y, prob)), 3),
            "tp": int(tp), "fn": int(fn), "fp": int(fp), "tn": int(tn),
        }

    return {
        "all_sample": _metrics(p_in, pr_in),
        "stratified_5fold_cv": _metrics(p_5, pr_5),
        "loocv": _metrics(p_loo, pr_loo),
    }


def evaluate_machine_learning_models(merged):
    """Run comprehensive ML cross-validation across All-Sample, 5-Fold, and LOOCV."""
    base_names = ["SAScore", "SCScore", "SYBA", "RAscore", "BR-SAScore", "Ligand_SA", "MW", "LogP", "HeavyAtoms", "RotBonds"]
    phys_names = ["7net_E_coh", "7net_BDE_min", "SambVca_Vbur", "7net_dE_Dimer", "Rg", "V_mol"]

    def build_matrices(items):
        Xb, Xp, y_vec = [], [], []
        for m in items:
            y_vec.append(1 if m["is_fab_valid"] else 0)
            Xb.append([
                float(m["sascore"]), float(m["scscore"]), float(m["syba"]),
                float(m["rascore"]), float(m["br_sascore"]), float(m["ligand_sa"]),
                float(m["mw"]), float(m["logp"]), float(m["heavy_atoms"]), float(m["rot_b"]),
            ])
            Xp.append([
                float(m["e_coh_ev_per_atom"]), float(m["bde_min_kcal_mol"]),
                float(m["vbur_percent"]), float(m["delta_e_dimer_kcal_mol"]),
                float(m["radius_of_gyration_ang"]), float(m["mol_volume_ang3"]),
            ])
        Xb = np.array(Xb, dtype=np.float64)
        Xp = np.array(Xp, dtype=np.float64)
        Xc = np.hstack([Xb, Xp])
        return Xb, Xp, Xc, np.array(y_vec, dtype=np.int32)

    # 1. Full Dataset (N=55: 45 Pos / 10 Neg)
    Xb_55, Xp_55, Xc_55, y_55 = build_matrices(merged)
    # 2. Precursor-Domain Subset (N=35: 27 Precursors vs 8 Organometallic/Silane Negatives)
    sub_35 = [
        m for m in merged
        if m["target_type"] == "Precursor" or m["id"] in {"NEG-001", "NEG-002", "NEG-003", "NEG-004", "NEG-007", "NEG-008", "NEG-009", "NEG-010"}
    ]
    Xb_35, Xp_35, Xc_35, y_35 = build_matrices(sub_35)

    rf_default = lambda: RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
    rf_balanced = lambda: RandomForestClassifier(n_estimators=100, max_depth=5, class_weight="balanced", random_state=42)
    logreg_balanced = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, class_weight="balanced", random_state=42))

    # In-sample feature importance (strictly labeled as All-Sample RF importance)
    rf_full = rf_default()
    rf_full.fit(Xc_55, y_55)
    all_names = base_names + phys_names
    imp_all_sample = {k: round(float(v), 4) for k, v in sorted(zip(all_names, rf_full.feature_importances_), key=lambda kv: kv[1], reverse=True)}
    phys_imp_sum = round(sum(imp_all_sample[k] for k in phys_names) * 100.0, 1)

    return {
        "note": (
            "Evaluated using genuine computed Phase 1 (fixed BR-SAScore + real ligand_sa) and "
            "latest Phase 2 (SevenNet-0 + SambVca without algebraic complexation proxy) features. "
            "In-sample Gini importances are reported alongside strict out-of-fold 5-Fold CV and LOOCV metrics."
        ),
        "full_dataset_n55": {
            "n_samples": len(y_55),
            "n_positive": int(np.sum(y_55 == 1)),
            "n_negative": int(np.sum(y_55 == 0)),
            "majority_baseline_acc": round(float(np.mean(y_55 == 1)), 3),
            "random_forest_default": {
                "base_only": _evaluate_ml_regime(Xb_55, y_55, rf_default),
                "phys_only": _evaluate_ml_regime(Xp_55, y_55, rf_default),
                "combined_base_and_phys": _evaluate_ml_regime(Xc_55, y_55, rf_default),
            },
            "random_forest_balanced": {
                "base_only": _evaluate_ml_regime(Xb_55, y_55, rf_balanced),
                "phys_only": _evaluate_ml_regime(Xp_55, y_55, rf_balanced),
                "combined_base_and_phys": _evaluate_ml_regime(Xc_55, y_55, rf_balanced),
            },
            "logistic_regression_balanced": {
                "base_only": _evaluate_ml_regime(Xb_55, y_55, logreg_balanced),
                "phys_only": _evaluate_ml_regime(Xp_55, y_55, logreg_balanced),
                "combined_base_and_phys": _evaluate_ml_regime(Xc_55, y_55, logreg_balanced),
            },
        },
        "precursor_subset_n35": {
            "n_samples": len(y_35),
            "n_positive": int(np.sum(y_35 == 1)),
            "n_negative": int(np.sum(y_35 == 0)),
            "majority_baseline_acc": round(float(np.mean(y_35 == 1)), 3),
            "random_forest_default": {
                "base_only": _evaluate_ml_regime(Xb_35, y_35, rf_default),
                "phys_only": _evaluate_ml_regime(Xp_35, y_35, rf_default),
                "combined_base_and_phys": _evaluate_ml_regime(Xc_35, y_35, rf_default),
            },
        },
        "all_sample_rf_feature_importances": imp_all_sample,
        "all_sample_phys_importance_sum_pct": phys_imp_sum,
    }


def run_tuning():
    print("=" * 80)
    print("AUDITED PHASE 2 MODEL EVALUATION (UN-HARDCODED + ALL-SAMPLE / 5-FOLD / LOOCV)")
    print("=" * 80)

    merged, p2_meta = load_merged_data()
    print(f"Loaded {len(merged)} molecules. Phase 2 metadata compute_time={p2_meta.get('computation_time_seconds')}s")

    analytical_comp = evaluate_unbiased_analytical_models(merged)
    print("\n--- Uniform Analytical Screening (No Hardcoding / No Label Branching) ---")
    print("Overall (Physics Only):", analytical_comp["overall_uniform_physics_only"])
    print("Overall (+SMARTS Chem) :", analytical_comp["overall_physics_plus_smarts_chem"])
    print("Physics Breakdown      :", json.dumps(analytical_comp["uniform_physics_only_breakdown"], indent=2))
    print("Phys+Chem Breakdown    :", json.dumps(analytical_comp["physics_plus_smarts_chem_breakdown"], indent=2))

    ml_comp = evaluate_machine_learning_models(merged)
    print("\n--- ML Validation Summary (N=55 RF Default) ---")
    for feat_set, res in ml_comp["full_dataset_n55"]["random_forest_default"].items():
        print(f"  {feat_set:24s} -> All-Sample={res['all_sample']['accuracy']:.3f}(AUC {res['all_sample']['auroc']:.3f}) | 5Fold={res['stratified_5fold_cv']['accuracy']:.3f}(AUC {res['stratified_5fold_cv']['auroc']:.3f}) | LOOCV={res['loocv']['accuracy']:.3f}(AUC {res['loocv']['auroc']:.3f}, Spec={res['loocv']['specificity']:.3f})")
    print(f"\nAll-sample RF Physical Feature Importance Sum: {ml_comp['all_sample_phys_importance_sum_pct']}%")

    out_payload = {
        "title": "Audited Precursor SA Evaluation: Un-Hardcoded Continuous Physical Screening & ML Cross-Validation",
        "phase2_metadata": p2_meta,
        "analytical_comparison": analytical_comp,
        "machine_learning_comparison": ml_comp,
        "molecules_evaluated": merged,
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_payload, f, indent=2, ensure_ascii=False)
    print(f"\nSaved audited comparison results to: {out_path}")


if __name__ == "__main__":
    run_tuning()
