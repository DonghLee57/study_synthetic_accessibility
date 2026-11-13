#!/usr/bin/env python3
"""
Unified Latest Model Evaluation Across the Full 160-Molecule Evaluation Set
---------------------------------------------------------------------------
Consolidates all latest architectural improvements into a single canonical analysis
on the full updated evaluation dataset (160 molecules: 27 Audited Precursors +
75 Expanded Precursors + 18 AS-ALD Inhibitors + 40 Failure-Mechanism-Stratified
Negative Controls across 6 failure classes):

  1. Latest TMC-SA v2 + Coordination/Topology Gate + Adaptive %V_bur (V_min = 34% / 44%)
  2. Unified 3D %V_bur coverage across all 160 molecules:
       - Exact SevenNet-0 LBFGS relaxed 3D %V_bur, BDE_min, DeltaE_dimer, E_coh (85 mols)
       - Deterministic RDKit ETKDGv3 3D %V_bur for remaining molecules (75 mols)
  3. 160-molecule 2D Latent Space Projections (PCA & t-SNE) + Separation Metrics
  4. Unified Supervised & One-Class Classifier Comparison (CV-5 & LOFMO Extrapolation)
  5. Six Failure-Mechanism-Stratified Negative Control TNR Breakdown (40 Negative Controls)

Outputs:
  - database/benchmark_unified_latest_160.json
"""

from __future__ import annotations

import json
import os
import sys
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import xgboost as xgb
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, RDConfig, rdMolDescriptors
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import balanced_accuracy_score, recall_score, roc_auc_score, silhouette_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import OneClassSVM, SVC

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(RDConfig.RDContribDir, "SA_Score"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts", "benchmark"))

import sascorer  # type: ignore
from syba.syba import SybaClassifier  # type: ignore
from run_phase2_sevenn_features import calc_buried_volume, find_central_atom_idx
from run_refinement_and_rxn_check import compute_topology_descriptors, evaluate_adaptive_tmc_sa_v2
from run_expanded_neg_control_validation import compute_expanded_descriptors, evaluate_coord_chem_gate

try:
    from BRSAScore import SAScorer as BRSAScorer  # type: ignore
except Exception:
    BRSAScorer = None


def compute_ligand_sa(
    smi: str,
    sa_val: float,
    br_val: float,
    br_scorer: Any = None,
) -> float:
    """Pro-ligand decomposition SA (exact run_phase1_baseline logic extended to Ga/In/As/Rh/Ir)."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return sa_val
    tm_metals = {"Ti", "Zr", "Hf", "Al", "Ru", "Co", "W", "Mo", "Ta", "Sn", "Cu", "Zn", "Fe", "Ir", "Rh", "Ga", "In", "As"}
    syms = {a.GetSymbol() for a in mol.GetAtoms()}
    centers = (syms & tm_metals) or (syms & {"Si"})
    if centers:
        rw = Chem.RWMol(mol)
        idxs = sorted([a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in centers], reverse=True)
        for c_idx in idxs:
            rw.RemoveAtom(c_idx)
        frags = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
        lig_scores: List[float] = []
        for frag in frags:
            try:
                for a in frag.GetAtoms():
                    a.SetNumRadicalElectrons(0)
                    a.SetNoImplicit(False)
                    a.SetNumExplicitHs(0)
                Chem.SanitizeMol(frag)
                l_smi = Chem.MolToSmiles(frag)
                if l_smi and l_smi not in ("[H]", "Cl", "O"):
                    cand: List[float] = []
                    if br_scorer is not None:
                        try:
                            bs, _ = br_scorer.calculateScore(l_smi)
                            cand.append(float(bs))
                        except Exception:
                            pass
                    l_mol = Chem.MolFromSmiles(l_smi)
                    if l_mol is not None:
                        cand.append(float(sascorer.calculateScore(l_mol)))
                    if cand:
                        lig_scores.append(min(cand))
            except Exception:
                pass
        return round(max(lig_scores), 3) if lig_scores else 1.00
    return round(min(sa_val, br_val) if br_val < 9.9 else sa_val, 3)


def compute_etkdg_vbur(smiles: str) -> Tuple[float, str]:
    """Generate deterministic RDKit ETKDGv3 3D conformer and compute %V_bur (R=3.5 A)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0, "C"
    center_idx, center_sym = find_central_atom_idx(mol)
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    cid = AllChem.EmbedMolecule(mol_h, params)
    if cid < 0:
        cid = AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=42)
    if cid >= 0:
        try:
            AllChem.UFFOptimizeMolecule(mol_h, maxIters=120)
        except Exception:
            pass
        vbur = calc_buried_volume(mol_h, center_idx, radius=3.5, grid_step=0.25)
        return float(vbur), center_sym
    # Fallback minimal estimate if 3D embedding fails for unusual organometallic connectivity
    return 0.0, center_sym


def main() -> None:
    p2_path = os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json")
    exp_db_path = os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json")
    neg_db_path = os.path.join(PROJECT_ROOT, "database", "expanded_negative_controls.json")
    exp_7n_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_elements_sevenn_features.json")
    refine_path = os.path.join(PROJECT_ROOT, "database", "benchmark_refinement_and_rxn_check.json")

    with open(p2_path, "r", encoding="utf-8") as f:
        p2_rows = json.load(f)["molecules_evaluated"]
    with open(exp_db_path, "r", encoding="utf-8") as f:
        exp_db = json.load(f)
    with open(neg_db_path, "r", encoding="utf-8") as f:
        neg_db = json.load(f)
    with open(exp_7n_path, "r", encoding="utf-8") as f:
        exp_7n_rows = json.load(f).get("features", [])
    with open(refine_path, "r", encoding="utf-8") as f:
        refine_res = json.load(f)

    syba_clf = SybaClassifier()
    syba_clf.fitDefaultScore()
    br_scorer = BRSAScorer() if BRSAScorer is not None else None

    # Index SevenNet-0 computed features by canonical/raw SMILES & ID
    sevenn_by_smi: Dict[str, Dict[str, Any]] = {}
    for r in p2_rows:
        sevenn_by_smi[r["smiles"]] = r
    for r in exp_7n_rows:
        smi = r.get("smiles", "")
        if smi and smi not in sevenn_by_smi:
            sevenn_by_smi[smi] = r

    # 1. Assemble full 160-molecule Evaluation Set
    raw_entries: List[Dict[str, Any]] = []
    audited_pos_smis = set()

    for r in p2_rows:
        if r["target_type"] in ("Precursor", "Inhibitor"):
            audited_pos_smis.add(r["smiles"])
            raw_entries.append({
                "id": r["id"],
                "name": r["name"],
                "smiles": r["smiles"],
                "group": "Audited Precursor (27)" if r["target_type"] == "Precursor" else "AS-ALD Inhibitor (18)",
                "target_type": r["target_type"],
                "label": 1,
                "failure_mechanism": "pos",
                "origin": r.get("domain", r["target_type"]),
            })

    for r in exp_db["precursors_102"]["records"]:
        if r["smiles"] not in audited_pos_smis:
            raw_entries.append({
                "id": r["id"],
                "name": r["name"],
                "smiles": r["smiles"],
                "group": "Expanded Precursor (75)",
                "target_type": "Precursor",
                "label": 1,
                "failure_mechanism": "pos",
                "origin": r.get("element", "Expanded"),
            })

    for r in neg_db["records"]:
        raw_entries.append({
            "id": r["id"],
            "name": r["name"],
            "smiles": r["smiles"],
            "group": "Negative Control (40)",
            "target_type": "Negative Control",
            "label": 0,
            "failure_mechanism": r["failure_mechanism"],
            "origin": r["origin"],
        })

    assert len(raw_entries) == 160, f"Expected 160 entries, got {len(raw_entries)}"

    # 2. Compute unified descriptors, 3D %V_bur, and latest TMC-SA v2 + Coord/Topology Gate
    records: List[Dict[str, Any]] = []
    mech_stats: Dict[str, Dict[str, int]] = defaultdict(lambda: {
        "total": 0,
        "sascore_rejected": 0,
        "syba_rejected": 0,
        "unified_rejected": 0,
    })

    for entry in raw_entries:
        smi = entry["smiles"]
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES: {smi}"
        formula = rdMolDescriptors.CalcMolFormula(mol)
        sa_val = float(sascorer.calculateScore(mol))
        syba_val = float(syba_clf.predict(smi))
        br_val = float(br_scorer.calculateScore(smi)[0]) if br_scorer is not None else sa_val

        desc = compute_expanded_descriptors(smi)
        topo = compute_topology_descriptors(smi)
        phys_rec = sevenn_by_smi.get(smi)

        # Compute pro-ligand SA cleanly (use stored audited ligand_sa for 55-set exact parity)
        if phys_rec is not None and "ligand_sa" in phys_rec:
            lig_sa = float(phys_rec["ligand_sa"])
        else:
            lig_sa = compute_ligand_sa(smi, sa_val, br_val, br_scorer=br_scorer)

        # Coordination/topology gate checks on structural chemistry + ligand_sa
        gate = evaluate_coord_chem_gate(smi, lig_sa, desc, topo)

        if phys_rec is not None:
            vbur_val = float(phys_rec["vbur_percent"])
            bde_val: Optional[float] = float(phys_rec["bde_min_kcal_mol"])
            ddim_val: Optional[float] = float(phys_rec["delta_e_dimer_kcal_mol"])
            ecoh_val: Optional[float] = float(phys_rec["e_coh_ev_per_atom"])
            center_sym = str(phys_rec.get("center_atom", phys_rec.get("center_element_expected", "")))
            if not center_sym:
                _, center_sym = find_central_atom_idx(mol)
            vbur_source = "SevenNet-0-LBFGS"
            if "ligand_sa" in phys_rec:
                v2_eval = evaluate_adaptive_tmc_sa_v2(phys_rec, topo)
                tmc_sa_v2 = float(v2_eval["tmc_sa_v2"])
                v_min = float(v2_eval["v_min_used"])
                tmc_pass_v2 = (tmc_sa_v2 <= 4.0) and gate["passed"]
            else:
                v_min = 34.0 if (center_sym == "Si" or topo["pi_hapticity_index"] > 0) else 44.0
                # For expanded elements, include smooth steric crowding penalty >83.5% + ligand SA
                steric_pen = 3.5 / (1.0 + np.exp(-max(-20.0, min(20.0, vbur_val - 83.5)))) if center_sym in ("Ir", "Rh", "As", "Ga", "In", "Ru", "Ti", "Zr", "Hf", "Al", "Co", "W", "Mo", "Ta", "Sn", "Cu", "Si") else 0.0
                tmc_sa_v2 = round(min(10.0, lig_sa + float(steric_pen)), 3)
                tmc_pass_v2 = (tmc_sa_v2 <= 4.0) and gate["passed"]
        else:
            vbur_val, center_sym = compute_etkdg_vbur(smi)
            bde_val = None
            ddim_val = None
            ecoh_val = None
            vbur_source = "RDKit-ETKDGv3-3D"
            v_min = 34.0 if (center_sym == "Si" or topo["pi_hapticity_index"] > 0) else 44.0
            has_center = center_sym in {"Ti", "Zr", "Hf", "Al", "Ru", "Co", "W", "Mo", "Ta", "Sn", "Cu", "Si", "Ir", "Rh", "Ga", "In", "As"}
            steric_pen = 3.5 / (1.0 + np.exp(-max(-20.0, min(20.0, vbur_val - 83.5)))) if has_center else 0.0
            tmc_sa_v2 = round(min(10.0, lig_sa + float(steric_pen)), 3)
            tmc_pass_v2 = (tmc_sa_v2 <= 4.0) and gate["passed"]

        coord_risk_count = (
            int(desc["cooh_count"] > 0)
            + int(desc["is_inert_alkane"] > 0)
            + int(topo["is_homoleptic_tm_sp3_refined"] or (desc["is_homoleptic_tm_sp3"] > 0 and topo["pi_hapticity_index"] == 0))
            + int(desc["metal_oh_count"] >= 2)
            + int(desc["aromatic_rings"] >= 3 and desc["metal_c_sigma_frac"] >= 0.9 and any(a.GetSymbol() not in ("C", "H", "O", "N", "F", "Cl", "Br", "I", "S", "P") for a in mol.GetAtoms()))
            + int(desc["strained_si_ring"] > 0)
        )

        # Continuous Risk Index combines SA deviation + steric crowding + coordination/reactivity red flags
        risk_index = round(
            max(0.0, (sa_val - 4.0) * 0.5)
            + 2.5 * coord_risk_count
            + (1.2 if (vbur_val > 83.5 and center_sym not in ("C", "N", "O", "F")) else 0.0),
            3,
        )

        mech = entry["failure_mechanism"]
        if entry["label"] == 0:
            mech_stats[mech]["total"] += 1
            if sa_val > 4.0:
                mech_stats[mech]["sascore_rejected"] += 1
            if syba_val < 0.0:
                mech_stats[mech]["syba_rejected"] += 1
            if not tmc_pass_v2:
                mech_stats[mech]["unified_rejected"] += 1

        records.append({
            **entry,
            "formula": formula,
            "center_atom": center_sym,
            "sascore": round(sa_val, 3),
            "syba": round(syba_val, 3),
            "br_sascore": round(br_val, 3),
            "ligand_sa": round(lig_sa, 3),
            "vbur_percent": round(vbur_val, 2),
            "v_min_adaptive": round(v_min, 1),
            "vbur_margin": round(vbur_val - v_min, 2),
            "vbur_source": vbur_source,
            "bde_min_kcal_mol": round(bde_val, 2) if bde_val is not None else None,
            "delta_e_dimer_kcal_mol": round(ddim_val, 2) if ddim_val is not None else None,
            "e_coh_ev_per_atom": round(ecoh_val, 3) if ecoh_val is not None else None,
            "pi_hapticity_index": int(topo["pi_hapticity_index"]),
            "chelate_ring_count": int(topo["chelate_ring_count"]),
            "is_homoleptic_tm_sp3_refined": bool(topo["is_homoleptic_tm_sp3_refined"]),
            "tmc_sa_v2": round(tmc_sa_v2, 3),
            "unified_pass": bool(tmc_pass_v2),
            "gate_reasons": gate["reasons"],
            "coord_risk_count": int(coord_risk_count),
            "risk_index": risk_index,
            "descriptors": desc,
        })

    # 3. Group & Overall Summary Metrics
    precursors_all = [r for r in records if r["target_type"] == "Precursor"]
    prec_aud = [r for r in records if r["group"] == "Audited Precursor (27)"]
    prec_exp = [r for r in records if r["group"] == "Expanded Precursor (75)"]
    inhib_all = [r for r in records if r["target_type"] == "Inhibitor"]
    pos_all = [r for r in records if r["label"] == 1]
    neg_all = [r for r in records if r["label"] == 0]

    def _rate(num: int, den: int) -> float:
        return round(100.0 * num / den, 1) if den else 0.0

    prec_tp = sum(1 for r in precursors_all if r["unified_pass"])
    prec_aud_tp = sum(1 for r in prec_aud if r["unified_pass"])
    prec_exp_tp = sum(1 for r in prec_exp if r["unified_pass"])
    inhib_tp = sum(1 for r in inhib_all if r["unified_pass"])
    pos_tp = sum(1 for r in pos_all if r["unified_pass"])
    neg_tn = sum(1 for r in neg_all if not r["unified_pass"])
    neg_fp = [r["id"] for r in neg_all if r["unified_pass"]]

    # Baseline SA / SYBA on Negative Controls & Positives
    neg_sas_tn = sum(1 for r in neg_all if r["sascore"] > 4.0)
    neg_syba_tn = sum(1 for r in neg_all if r["syba"] < 0.0)
    prec_sas_tp = sum(1 for r in precursors_all if r["sascore"] <= 3.5)

    # 4. Unified 160-Molecule Descriptor Feature Matrix -> PCA & t-SNE + CV-5 & LOFMO
    desc_keys = sorted(records[0]["descriptors"].keys())
    feature_names = ["sascore", "syba", "vbur_percent", "v_min_adaptive", "pi_hapticity_index", "chelate_ring_count"] + desc_keys

    X = np.array(
        [
            [
                r["sascore"],
                r["syba"],
                r["vbur_percent"],
                r["v_min_adaptive"],
                float(r["pi_hapticity_index"]),
                float(r["chelate_ring_count"]),
            ]
            + [float(r["descriptors"][k]) for k in desc_keys]
            for r in records
        ],
        dtype=float,
    )
    y = np.array([r["label"] for r in records], dtype=int)

    X_std = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_std)

    tsne = TSNE(n_components=2, perplexity=25, random_state=42, init="pca", learning_rate="auto", max_iter=1000)
    X_tsne = tsne.fit_transform(X_std)

    sil_pca = float(silhouette_score(X_pca, y))
    sil_tsne = float(silhouette_score(X_tsne, y))
    linear_sep_pca = float(
        make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", random_state=42))
        .fit(X_pca, y)
        .score(X_pca, y)
    )
    # Balanced accuracy in 2D PCA & Full Space
    pca_lr_bal = float(balanced_accuracy_score(
        y,
        LogisticRegression(class_weight="balanced", random_state=42).fit(X_pca, y).predict(X_pca),
    ))
    full_svc_bal = float(balanced_accuracy_score(
        y,
        SVC(kernel="rbf", class_weight="balanced", random_state=42).fit(X_std, y).predict(X_std),
    ))

    for i, r in enumerate(records):
        r["pc1"] = round(float(X_pca[i, 0]), 4)
        r["pc2"] = round(float(X_pca[i, 1]), 4)
        r["tsne1"] = round(float(X_tsne[i, 0]), 4)
        r["tsne2"] = round(float(X_tsne[i, 1]), 4)

    # 5. CV-5 and LOFMO Benchmarks across 160 molecules using expanded + adaptive %V_bur/topology features
    cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    scale_pos = n_neg / n_pos

    models = {
        "LogisticRegression": make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42)),
        "SVC_RBF": make_pipeline(StandardScaler(), SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=42)),
        "RandomForest": RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42),
        "XGBoost_Weighted": xgb.XGBClassifier(
            n_estimators=150, max_depth=3, learning_rate=0.05,
            scale_pos_weight=scale_pos, eval_metric="logloss", random_state=42,
        ),
    }

    ml_cv5: Dict[str, Any] = {}
    for mname, clf in models.items():
        prob = cross_val_predict(clf, X, y, cv=cv5, method="predict_proba")[:, 1]
        pred = (prob >= 0.5).astype(int)
        ml_cv5[mname] = {
            "cv5_roc_auc": round(float(roc_auc_score(y, prob)), 4),
            "cv5_balanced_acc": round(float(balanced_accuracy_score(y, pred)), 4),
            "cv5_tpr_sensitivity": round(float(recall_score(y, pred, pos_label=1)), 4),
            "cv5_tnr_specificity": round(float(recall_score(y, pred, pos_label=0)), 4),
        }

    # OneClassSVM_RBF trained strictly on Positive folds
    oc_scores = np.zeros(len(y), dtype=float)
    oc_preds = np.zeros(len(y), dtype=int)
    for tr_idx, te_idx in cv5.split(X, y):
        pipe = make_pipeline(StandardScaler(), OneClassSVM(kernel="rbf", nu=0.08, gamma="scale"))
        pipe.fit(X[tr_idx][y[tr_idx] == 1])
        oc_scores[te_idx] = pipe.decision_function(X[te_idx])
        oc_preds[te_idx] = (pipe.predict(X[te_idx]) == 1).astype(int)

    ml_cv5["OneClassSVM_RBF"] = {
        "cv5_roc_auc": round(float(roc_auc_score(y, oc_scores)), 4),
        "cv5_balanced_acc": round(float(balanced_accuracy_score(y, oc_preds)), 4),
        "cv5_tpr_sensitivity": round(float(recall_score(y, oc_preds, pos_label=1)), 4),
        "cv5_tnr_specificity": round(float(recall_score(y, oc_preds, pos_label=0)), 4),
    }

    # LOFMO Extrapolation
    row_mechs = [r["failure_mechanism"] for r in records]
    unique_mechs = sorted({m for m in row_mechs if m != "pos"})
    lofmo: Dict[str, Any] = {}
    for mname, fn in [
        ("LogisticRegression", lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42))),
        ("SVC_RBF", lambda: make_pipeline(StandardScaler(), SVC(kernel="rbf", class_weight="balanced", random_state=42))),
        ("RandomForest", lambda: RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42)),
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
    oc_mech = {}
    for hold in unique_mechs:
        te = np.array([m == hold for m in row_mechs])
        oc_mech[hold] = round(float(np.mean(oc_full.predict(X[te]) == -1)), 4)
    lofmo["OneClassSVM_RBF"] = {
        "mean_unseen_tnr": round(float(np.mean(list(oc_mech.values()))), 4),
        "per_unseen_mechanism_tnr": oc_mech,
    }

    mech_summary: Dict[str, Any] = {}
    for mech, st in mech_stats.items():
        tot = st["total"]
        mech_summary[mech] = {
            **st,
            "sascore_tnr_pct": _rate(st["sascore_rejected"], tot),
            "syba_tnr_pct": _rate(st["syba_rejected"], tot),
            "unified_tnr_pct": _rate(st["unified_rejected"], tot),
        }

    lofmo["Unified_TMC_SA_v2_Gate"] = {
        "mean_unseen_tnr": round(float(np.mean([mech_summary[m]["unified_tnr_pct"] / 100.0 for m in unique_mechs])), 4),
        "per_unseen_mechanism_tnr": {m: round(mech_summary[m]["unified_tnr_pct"] / 100.0, 4) for m in unique_mechs},
    }

    out = {
        "dataset_summary": {
            "total_molecules": len(records),
            "positive_total": len(pos_all),
            "precursors_total": len(precursors_all),
            "audited_precursors": len(prec_aud),
            "expanded_precursors": len(prec_exp),
            "inhibitors_total": len(inhib_all),
            "negative_controls_total": len(neg_all),
            "sevenn_relaxed_count": sum(1 for r in records if r["vbur_source"] == "SevenNet-0-LBFGS"),
            "etkdg_3d_count": sum(1 for r in records if r["vbur_source"] == "RDKit-ETKDGv3-3D"),
            "feature_count": len(feature_names),
            "feature_names": feature_names,
        },
        "unified_v2_summary": {
            "precursors_all_tp": f"{prec_tp}/{len(precursors_all)}",
            "precursors_all_tpr_pct": _rate(prec_tp, len(precursors_all)),
            "precursors_all_fnr_pct": round(100.0 - _rate(prec_tp, len(precursors_all)), 1),
            "audited_precursors_tp": f"{prec_aud_tp}/{len(prec_aud)}",
            "audited_precursors_tpr_pct": _rate(prec_aud_tp, len(prec_aud)),
            "expanded_precursors_tp": f"{prec_exp_tp}/{len(prec_exp)}",
            "expanded_precursors_tpr_pct": _rate(prec_exp_tp, len(prec_exp)),
            "inhibitors_tp": f"{inhib_tp}/{len(inhib_all)}",
            "inhibitors_tpr_pct": _rate(inhib_tp, len(inhib_all)),
            "positive_all_tp": f"{pos_tp}/{len(pos_all)}",
            "positive_all_tpr_pct": _rate(pos_tp, len(pos_all)),
            "negative_all_tn": f"{neg_tn}/{len(neg_all)}",
            "negative_all_tnr_pct": _rate(neg_tn, len(neg_all)),
            "negative_all_fpr_pct": round(100.0 - _rate(neg_tn, len(neg_all)), 1),
            "negative_fp_ids": neg_fp,
            "baseline_sascore_precursors_tp": f"{prec_sas_tp}/{len(precursors_all)}",
            "baseline_sascore_precursors_tpr_pct": _rate(prec_sas_tp, len(precursors_all)),
            "baseline_sascore_neg_tn": f"{neg_sas_tn}/{len(neg_all)}",
            "baseline_sascore_neg_tnr_pct": _rate(neg_sas_tn, len(neg_all)),
            "baseline_syba_neg_tn": f"{neg_syba_tn}/{len(neg_all)}",
            "baseline_syba_neg_tnr_pct": _rate(neg_syba_tn, len(neg_all)),
        },
        "latent_space_160": {
            "pca_explained_variance_ratio": [round(float(v), 4) for v in pca.explained_variance_ratio_],
            "pca_silhouette": round(sil_pca, 4),
            "tsne_silhouette": round(sil_tsne, 4),
            "pca_2d_balanced_acc": round(pca_lr_bal, 4),
            "full_space_svc_balanced_acc": round(full_svc_bal, 4),
            "pca_linear_acc": round(linear_sep_pca, 4),
        },
        "by_failure_mechanism": mech_summary,
        "ml_cv5_160": ml_cv5,
        "lofmo_extrapolation_160": lofmo,
        "item2_reaction_energy_summary": refine_res["item_2_reaction_energy_check"]["reactions"],
        "records": records,
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(f"Saved unified 160-molecule analysis -> {out_path}")
    print(json.dumps({k: v for k, v in out.items() if k not in ("records", "item2_reaction_energy_summary")}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
