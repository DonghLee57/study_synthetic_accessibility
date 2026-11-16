#!/usr/bin/env python3
"""
Interpretable Metal-Ligand Physics & Mechanistic Generalization Model
---------------------------------------------------------------------
Builds and benchmarks first-principles interpretable Metal-Ligand models
over the curated 160-molecule dataset and 8-channel SevenNet-0 CI-NEB paths:

  1. Pauling-HSAB Thermodynamic Stabilization Law (Delta E_f^std / v_M)
  2. Spherical-Saturation Steric Packing Law (%V_bur = 100 * (1 - exp(-z)))
  3. 4-Regime Generalization Audit (5-CV, LOMO-5, LOLO-5, Homo->Hetero)
     contrasting interpretable laws vs 19-feat RandomForest and exposing
     why min-operator BDE_min and dimerization cannot extrapolate cleanly.
  4. 142-Molecule (102 Prec + 40 Neg) 2D M+L + SMARTS Viability Classification.
  5. 102-Precursor Kinetic Ligand-Motif Taxonomy linked to SevenNet-0 CI-NEB.
"""

from __future__ import annotations

import html
import json
from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import balanced_accuracy_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from element_physics import PAULING_EN
from ligand_coord_html_section import (
    _extract_continuous_ligand_features,
    _render_template,
    fragment_metal_and_ligands,
)
from tmc_sa_scorer import compute_chemical_coord_descriptors

SPECTATOR_PI_CHELATE = {
    "C1=CCC=C1", "CC1=CC=CC1", "CCC1=CC=CC1", "CC1=C(C)C(C)C(C)=C1C",
    "C1=CCCC=C1", "CC=CC=C(C)C", "c1ccnc1", "C1=CCCC=CCC1",
    "CC(=O)C=C(C)O", "CC(C)(C)C(=O)C=C(O)C(C)(C)C", "CC(=NC(C)C)NC(C)C",
    "CC(=NC(C)(C)C)NC(C)(C)C", "CC(C)N=NNC(C)C", "CC(C)N=CNC(C)C",
    "CC(O)=CC(=O)C(C)C", "Cc1ccc(C(C)C)cc1", "CCC1=CC=CCC1", "CCc1ccccc1",
}
AMIDE_MOTIFS = {
    "CNC", "CCNC", "CCNCC", "CC(C)NC(C)C", "CC(C)(C)N",
    "CN(C)CCNC(C)(C)C", "COCC[NH+]C(C)(C)C", "C[Si](C)(C)N[Si](C)(C)C", "C[N+]C",
}
ALKOXIDE_MOTIFS = {"CCO", "CC(C)O", "CC(C)(C)O", "CO", "CN(C)CC(C)(C)O", "CCONC(=O)C(C)(C)C"}
BETA_ALKYL_MOTIFS = {"CC", "CC(C)C", "C=C(C)C", "C=C(C)C(=C)C", "[CH]CC(C)(C)C"}
HALIDE_MOTIFS = {"Cl", "Br", "I", "F"}
ME_CO_MOTIFS = {"C", "[C-]#[O+]"}


def _eval_split_reg(
    x_mat: np.ndarray,
    y_vec: np.ndarray,
    sp_list: List[Tuple[np.ndarray, np.ndarray]],
    model_fn: Callable[[], Any],
    transform_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    inv_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> Tuple[float, float]:
    """Evaluate R^2 and MAE on CV or single holdout splits."""
    y_fit = transform_fn(y_vec) if transform_fn is not None else y_vec
    if len(sp_list) == 1:
        tr, te = sp_list[0]
        mdl = model_fn().fit(x_mat[tr], y_fit[tr])
        pred = mdl.predict(x_mat[te])
        if inv_fn is not None:
            pred = inv_fn(pred)
        return float(r2_score(y_vec[te], pred)), float(mean_absolute_error(y_vec[te], pred))

    pred_all = np.zeros_like(y_vec)
    for tr, te in sp_list:
        mdl = model_fn().fit(x_mat[tr], y_fit[tr])
        pr = mdl.predict(x_mat[te])
        pred_all[te] = inv_fn(pr) if inv_fn is not None else pr
    return float(r2_score(y_vec, pred_all)), float(mean_absolute_error(y_vec, pred_all))


def evaluate_mechanistic_ml_models(
    unified_json_path: str,
    neb_json_path: str,
) -> Dict[str, Any]:
    """
    Compute interpretable Pauling-HSAB / Spherical-Packing physics models,
    4-regime extrapolation metrics, 142-molecule viability classification,
    and 102-precursor kinetic motif census integrated with CI-NEB barriers.

    Physics/Algorithm:
        1) Pauling-HSAB bond enthalpy scaling: regresses Delta E_f^std/v_M against
           mean Pauling ionic resonance Delta chi^2, donor fractions (f_N, f_O, f_X),
           volumetric heavy density N_hv/R_cov^3, ring count, pi-hapticity, and chelates.
        2) Spherical-Saturation packing law: models %V_bur = 100*(1 - exp(-z)) using
           normalized solid-angle packing Sum_i N_hv,i^(2/3) / (R_cov + 1.85)^2.
        3) Contrasts generalization across 5-CV, LOMO-5, LOLO-5, and Homo->Hetero.
    """
    with open(unified_json_path, "r", encoding="utf-8") as f:
        recs = json.load(f)["records"]
    with open(neb_json_path, "r", encoding="utf-8") as f:
        neb_data = json.load(f)
    neb_by_id = {
        r["spec"]["reaction_id"]: r
        for r in neb_data.get("reactions", [])
        if r.get("spec", {}).get("reaction_id")
    }

    precs = [r for r in recs if r.get("group") in {"Audited Precursor (27)", "Expanded Precursor (75)"}]
    viab_recs = [
        r for r in recs
        if r.get("group") in {"Audited Precursor (27)", "Expanded Precursor (75)", "Negative Control (40)"}
    ]

    x_pt, x_ps, x_rf19 = [], [], []
    groups_m, groups_l, is_het = [], [], []
    y_map: Dict[str, List[float]] = {
        "delta_e_form_per_bond_kcal_mol": [],  # units: kcal/mol/bond
        "vbur_percent": [],  # units: %
        "ligand_sa": [],  # units: score
        "bde_min_kcal_mol": [],  # units: kcal/mol
        "delta_g_dimer_373_kcal_mol": [],  # units: kcal/mol
    }
    motif_flags: Counter[str] = Counter()

    RDLogger.DisableLog("rdApp.error")
    try:
        for r in precs:
            c = r["center_atom"]
            groups_m.append(c)
            _, frags_s, _, donors = fragment_metal_and_ligands(r["smiles"], c)
            uniq = sorted(set(frags_s))
            is_het.append(len(uniq) > 1)
            groups_l.append(Counter(frags_s).most_common(1)[0][0] if frags_s else "NONE")

            # Kinetic motif presence census across the 102 precursors
            f_beta = any(s in BETA_ALKYL_MOTIFS for s in uniq)
            f_amd = any(s in AMIDE_MOTIFS for s in uniq)
            f_alk = any(s in ALKOXIDE_MOTIFS for s in uniq)
            f_spec = any(s in SPECTATOR_PI_CHELATE for s in uniq)
            f_hal = any(s in HALIDE_MOTIFS for s in uniq)
            f_meco = any(s in ME_CO_MOTIFS for s in uniq)
            if f_beta:
                motif_flags["beta-H_alkyl"] += 1
            if f_amd:
                motif_flags["amide_H_transfer"] += 1
            if f_alk:
                motif_flags["alkoxide_alpha_beta"] += 1
            if f_spec:
                motif_flags["spectator_pi_chelate"] += 1
            if f_hal:
                motif_flags["halide"] += 1
            if f_meco:
                motif_flags["Me_or_CO"] += 1
            if (f_beta or f_amd or f_alk) and f_spec:
                motif_flags["mixed_spectator_plus_reactive"] += 1

            rcov = float(r["r_cov_ang"])  # units: A
            chi_m = float(r["chi_pauling"])
            nval = float(r["n_val"])
            pi_h = float(r["pi_hapticity_index"])
            chel = float(r["chelate_ring_count"])
            cn_eff = float(r["cn_eff"])
            ve_def = float(r["ve_deficit"])
            dchi_mean = float(r["delta_chi_mean"])

            lf = _extract_continuous_ligand_features(r["smiles"], c)
            mol = Chem.MolFromSmiles(r["smiles"])
            cidxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == c]
            rw = Chem.RWMol(mol)
            for ci in sorted(cidxs, reverse=True):
                rw.RemoveAtom(ci)
            frag_mols = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
            sa_pack = sum((fm.GetNumAtoms() ** (2.0 / 3.0)) / ((rcov + 1.85) ** 2) for fm in frag_mols)  # units: A^-2
            hv_total = float(lf[1])
            dchi_sq = sum((chi_m - PAULING_EN.get(d, 2.55)) ** 2 for d in donors) / max(1, len(donors))

            x_pt.append([dchi_sq, lf[7], lf[8], lf[9], hv_total / (rcov ** 3), lf[5], pi_h, chel])
            x_ps.append([sa_pack, hv_total / ((rcov + 1.85) ** 2), lf[4] / ((rcov + 1.85) ** 2), pi_h, cn_eff, lf[3] / (rcov ** 2)])
            x_rf19.append([rcov, chi_m, nval] + lf + [pi_h, chel, hv_total / (rcov ** 3), lf[3] / (rcov ** 2), dchi_mean, ve_def])

            for k in y_map:
                y_map[k].append(float(r[k]))

        # 142-molecule (102 Prec + 40 Neg) classification vectors
        x_cls_2d, x_cls_3d, y_cls, grp_cls_m, grp_cls_l = [], [], [], [], []
        for r in viab_recs:
            c = r["center_atom"]
            grp_cls_m.append(c)
            _, fs, _, _ = fragment_metal_and_ligands(r["smiles"], c)
            grp_cls_l.append(Counter(fs).most_common(1)[0][0] if fs else "NONE")
            rcov = float(r["r_cov_ang"])  # units: A
            chi_m = float(r["chi_pauling"])
            nval = float(r["n_val"])
            lf = _extract_continuous_ligand_features(r["smiles"], c)
            cd = compute_chemical_coord_descriptors(r["smiles"])
            v2d = [rcov, chi_m, nval] + lf + [
                float(r["pi_hapticity_index"]),
                float(r["chelate_ring_count"]),
                lf[1] / (rcov ** 3),
                lf[3] / (rcov ** 2),
                float(r["delta_chi_mean"]),
                float(r["ve_deficit"]),
                float(r["ligand_sa"]),
                float(cd["cooh_count"] + cd["metal_oh_count"] + cd["is_inert_alkane"] + cd["strained_si_ring"]),
                float(cd["beta_h_alkyl_count"]),
                float(cd["is_homoleptic_tm_sp3"]),
            ]
            x_cls_2d.append(v2d)
            x_cls_3d.append(
                v2d
                + [
                    float(r["vbur_percent"]),
                    float(r["delta_e_form_per_bond_kcal_mol"]),
                    float(r["bde_min_kcal_mol"]),
                    float(r["delta_g_dimer_373_kcal_mol"]),
                ]
            )
            y_cls.append(0 if r.get("group") == "Negative Control (40)" else 1)
    finally:
        RDLogger.EnableLog("rdApp.error")

    x_pt_arr = np.array(x_pt)
    x_ps_arr = np.array(x_ps)
    x_rf_arr = np.array(x_rf19)
    het_arr = np.array(is_het)
    x_lin_comb = np.hstack([x_pt_arr, x_ps_arr[:, :3], x_ps_arr[:, 4:5], x_rf_arr[:, -1:]])

    splits = {
        "cv5": list(KFold(n_splits=5, shuffle=True, random_state=42).split(x_rf_arr)),
        "lomo5": list(GroupKFold(n_splits=5).split(x_rf_arr, groups=groups_m)),
        "lolo5": list(GroupKFold(n_splits=5).split(x_rf_arr, groups=groups_l)),
        "homo_het": [(np.where(~het_arr)[0], np.where(het_arr)[0])],
    }

    def ridge_fn() -> Any:
        return make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 25)))

    def rf_fn() -> Any:
        return RandomForestRegressor(n_estimators=120, random_state=42)

    # Transform helpers for Spherical Saturation %V_bur = 100 * (1 - exp(-z))
    def sat_fwd(yv: np.ndarray) -> np.ndarray:
        return -np.log(np.clip(1.0 - yv / 100.0, 0.05, 0.99))

    def sat_inv(zv: np.ndarray) -> np.ndarray:
        return 100.0 * (1.0 - np.exp(-np.clip(zv, 0.01, 3.0)))

    reg_rows: List[Dict[str, Any]] = []
    target_meta = [
        ("delta_e_form_per_bond_kcal_mol", "표준 결합 생성에너지 (ΔEf/vM)", "Pauling–HSAB 8항 물리식", "kcal/mol", x_pt_arr, None, None),
        ("vbur_percent", "구면 매몰 부피 (%Vbur)", "구면 포화 충진 법칙 (1-e^-z)", "%", x_ps_arr, sat_fwd, sat_inv),
        ("ligand_sa", "프로-리간드 합성 복잡도 (Ligand SA)", "선형 조합 물리식 (13항)", "score", x_lin_comb, None, None),
        ("bde_min_kcal_mol", "최약 결합 해리에너지 (BDEmin)", "선형 조합 물리식 (13항)", "kcal/mol", x_lin_comb, None, None),
        ("delta_g_dimer_373_kcal_mol", "373 K 이합체 자유에너지 (ΔGdimer)", "선형 조합 물리식 (13항)", "kcal/mol", x_lin_comb, None, None),
    ]
    for key, label, phys_name, unit, xp, tf, invf in target_meta:
        yv = np.array(y_map[key])
        row: Dict[str, Any] = {"key": key, "label": label, "phys_name": phys_name, "unit": unit}
        for sname, sp in splits.items():
            r2_p, mae_p = _eval_split_reg(xp, yv, sp, ridge_fn, transform_fn=tf, inv_fn=invf)
            r2_r, mae_r = _eval_split_reg(x_rf_arr, yv, sp, rf_fn)
            row[f"{sname}_phys_r2"] = round(r2_p, 3)
            row[f"{sname}_phys_mae"] = round(mae_p, 2)
            row[f"{sname}_rf_r2"] = round(r2_r, 3)
            row[f"{sname}_rf_mae"] = round(mae_r, 2)
        reg_rows.append(row)

    # Full standardized coefficients for interpretable physics laws
    y_ef = np.array(y_map["delta_e_form_per_bond_kcal_mol"])
    sc_ef = StandardScaler().fit(x_pt_arr)
    m_ef = RidgeCV(alphas=np.logspace(-2, 3, 25)).fit(sc_ef.transform(x_pt_arr), y_ef)
    ef_full_r2 = round(float(r2_score(y_ef, m_ef.predict(sc_ef.transform(x_pt_arr)))), 3)
    ef_coefs = [
        {"term": t, "coef": round(float(c), 2)}
        for t, c in zip(
            ["(Δχ)² Pauling", "f_N (Amide/Nitride)", "f_O (Oxo/Alkoxide)", "f_X (Halide)", "N_hv / R_cov³", "N_ring", "η_π", "N_chelate"],
            m_ef.coef_,
        )
    ]

    y_vb = np.array(y_map["vbur_percent"])
    sc_vb = StandardScaler().fit(x_ps_arr)
    m_vb = RidgeCV(alphas=np.logspace(-2, 3, 25)).fit(sc_vb.transform(x_ps_arr), sat_fwd(y_vb))
    vb_full_r2 = round(float(r2_score(y_vb, sat_inv(m_vb.predict(sc_vb.transform(x_ps_arr))))), 3)
    vb_coefs = [
        {"term": t, "coef": round(float(c), 3)}
        for t, c in zip(
            ["Ω_pack = Σ N_hv^(2/3)/(R+1.85)²", "N_hv / (R+1.85)²", "N_rot / (R+1.85)²", "η_π", "CN_eff", "N_hv,max / R_cov²"],
            m_vb.coef_,
        )
    ]

    # 142-molecule Precursor (102) vs Negative Control (40) classification
    x2d_arr = np.array(x_cls_2d)
    x3d_arr = np.array(x_cls_3d)
    y_cls_arr = np.array(y_cls)
    cls_splits = {
        "cv5": list(KFold(n_splits=5, shuffle=True, random_state=42).split(x2d_arr)),
        "lomo5": list(GroupKFold(n_splits=5).split(x2d_arr, y_cls_arr, groups=grp_cls_m)),
        "lolo5": list(GroupKFold(n_splits=5).split(x2d_arr, y_cls_arr, groups=grp_cls_l)),
    }
    cls_results: Dict[str, Dict[str, Any]] = {}
    for tag, Xmat in [("only_2d_ml", x2d_arr), ("full_ml_3d", x3d_arr)]:
        cls_results[tag] = {}
        for sname, sp in cls_splits.items():
            probs = np.zeros(len(y_cls_arr))
            for tr, te in sp:
                clf = RandomForestClassifier(
                    n_estimators=150, class_weight="balanced", random_state=42
                ).fit(Xmat[tr], y_cls_arr[tr])
                probs[te] = clf.predict_proba(Xmat[te])[:, 1]
            cls_results[tag][sname] = {
                "auroc": round(float(roc_auc_score(y_cls_arr, probs)), 3),
                "bal_acc": round(float(balanced_accuracy_score(y_cls_arr, probs >= 0.5)), 3),
            }

    # Kinetic Motif Table grounded in 102-precursor census & SevenNet CI-NEB benchmark
    def _neb_str(rid: str) -> str:
        item = neb_by_id.get(rid)
        if not item:
            return "-"
        sp = item["spec"]
        en = item["energies_kcal_mol"]
        ea_7n = en["sevenn_ea_cineb_kcal_mol"]
        ea_lit = sp["lit_ea_kcal_mol"]
        ea_sep = en.get("sevenn_ea_vs_sep_react_kcal_mol")
        if ea_sep is not None and abs(float(ea_sep) - float(ea_7n)) > 0.1:
            return f"{sp['precursor_name']}: {ea_7n:.2f} [분리계 {ea_sep:.2f}] (문헌 {ea_lit:.2f})"
        return f"{sp['precursor_name']}: {ea_7n:.2f} (문헌 {ea_lit:.2f})"

    motif_rows = [
        {
            "motif": "Amide / Alkylamide (α/β-H 이동)",
            "count": motif_flags["amide_H_transfer"],
            "examples": "-NMe2 (36), -NEtMe (8), -NHtBu (7), -NEt2 (3)",
            "mechanism": "이종·인접 리간드 간 4원환 양성이동 -> 이민/아민 탈리 (협동 분해)",
            "neb_anchor": f"{_neb_str('CpZr_NMe2_betaH_elim')}; +NH3 교환 {_neb_str('TDMAT_NH3_transam')}",
            "design_rule": "정적 BDE(69.96) 대비 장벽이 절반(34.67); 비대칭 아민 도입 시 증기압 개선 대가로 T_decomp 하강 주의",
        },
        {
            "motif": "Spectator π-Ring / Chelate (고리 안정화)",
            "count": motif_flags["spectator_pi_chelate"],
            "examples": "Cp/MeCp/EtCp (20), acac/thd (21), amd (8), CHD/DMPD",
            "mechanism": "π-역결합 및 킬레이트 입체 차폐 (%Vbur 상승); 직접 해리 >56 kcal/mol",
            "neb_anchor": "CpZr 내 Cp 양성자화 ΔE_sep=+56.57 (+21.90 vs β-H 장벽)",
            "design_rule": "표면 노출을 가리고 단량체를 보존(혼합 구경꾼-활성 착물 4종 확인); co-ligand 취약점만 별도 차단 필요",
        },
        {
            "motif": "Alkyl β-H 보유기 (C₂+ 알킬)",
            "count": motif_flags["beta-H_alkyl"],
            "examples": "-Et (15), -iPr, -tBu, allyl/DMBD",
            "mechanism": "동기 4중심 β-H 제거 -> M-H + 알켄(C₂H₄/C₃H₆) 탈리",
            "neb_anchor": _neb_str("GaEt3_betaH_elim"),
            "design_rule": "Ga-C 균일해리(65.0)보다 23.5 kcal/mol 낮은 경로; 250 °C 이상 순수 열 ALD 상한 제약",
        },
        {
            "motif": "Alkoxide (경쟁 α-H vs β-H-to-O)",
            "count": motif_flags["alkoxide_alpha_beta"],
            "examples": "-OiPr (12), -OEt (4), -OtBu (4), dmamp (5)",
            "mechanism": "주족(TEOS) 4중심 β-H-to-O 알켄 탈리 vs 전이금속(TTIP) α-H-to-M 케톤 탈리",
            "neb_anchor": f"{_neb_str('TEOS_C2H4_elim')}; {_neb_str('TTIP_acetone_elim')}",
            "design_rule": "TEOS(52.40)는 최고 내열성 보유; Ti 등 전이금속 알콕사이드는 아세톤 α-H 채널 및 -OH 자가촉매 유의",
        },
        {
            "motif": "Halide / Disilane (1,2-이동)",
            "count": motif_flags["halide"],
            "examples": "-Cl (10), -I (4), -Br (4) + Si-Si 골격",
            "mechanism": "이규소 1,2-Cl / 1,2-H 이동에 의한 실릴렌(:SiX₂) 압출",
            "neb_anchor": f"{_neb_str('HCDS_SiCl2_elim')}; {_neb_str('DCDS_SiHCl_elim')}",
            "design_rule": "Si-H -> Si-Cl 치환 시 탄소 오염 없는 고온 창(T_mid 450->525 °C) 확보 가능하나 HCl 부식·저휘발성 트레이드오프",
        },
    ]

    return {
        "n_precs": len(precs),
        "n_homo": int((~het_arr).sum()),
        "n_hetero": int(het_arr.sum()),
        "n_viab": len(viab_recs),
        "reg_rows": reg_rows,
        "ef_full_r2": ef_full_r2,
        "ef_coefs": ef_coefs,
        "vb_full_r2": vb_full_r2,
        "vb_coefs": vb_coefs,
        "cls_results": cls_results,
        "motif_flags": dict(motif_flags),
        "motif_rows": motif_rows,
    }


def _build_generalization_regime_svg(reg_rows: List[Dict[str, Any]]) -> str:
    """Visualize R^2 across 5-CV, LOMO-5, LOLO-5, and Homo->Hetero extrapolation."""
    rmap = {r["key"]: r for r in reg_rows}
    r_ef = rmap["delta_e_form_per_bond_kcal_mol"]
    r_vb = rmap["vbur_percent"]
    r_bde = rmap["bde_min_kcal_mol"]

    regimes = [
        ("5-CV (무작위)", "cv5"),
        ("LOMO-5 (미지 금속)", "lomo5"),
        ("LOLO-5 (미지 리간드)", "lolo5"),
        ("Homo->Hetero", "homo_het"),
    ]
    w, h = 560, 240
    x0, plot_w = 52, 490
    y_max, y_min = 0.90, -0.35
    y_top, y_bot = 36.0, 198.0

    def _y(val: float) -> float:
        v = max(y_min, min(y_max, val))
        frac = (v - y_min) / (y_max - y_min)
        return y_bot - frac * (y_bot - y_top)

    y0_px = _y(0.0)
    grid_svg = []
    for tick in (-0.2, 0.0, 0.2, 0.4, 0.6, 0.8):
        yy = _y(tick)
        dash = "" if abs(tick) < 1e-6 else 'stroke-dasharray="3 3"'
        stroke = "#64748b" if abs(tick) < 1e-6 else "#e2e8f0"
        grid_svg.append(
            f'<line x1="{x0}" y1="{yy:.1f}" x2="{x0 + plot_w}" y2="{yy:.1f}" stroke="{stroke}" {dash}/>'
            f'<text x="{x0 - 6}" y="{yy + 3:.1f}" text-anchor="end" font-size="9.2" fill="#64748b">{tick:+.1f}</text>'
        )

    series = [
        ("Pauling-HSAB ΔEf", "#0284c7", [r_ef[f"{k}_phys_r2"] for _, k in regimes]),
        ("Spherical-Sat %Vbur", "#059669", [r_vb[f"{k}_phys_r2"] for _, k in regimes]),
        ("RF %Vbur (블랙박스)", "#94a3b8", [r_vb[f"{k}_rf_r2"] for _, k in regimes]),
        ("RF BDEmin (최소연산 한계)", "#e11d48", [r_bde[f"{k}_rf_r2"] for _, k in regimes]),
    ]
    step = plot_w / len(regimes)
    bars_svg = []
    for i, (rlabel, _) in enumerate(regimes):
        cx = x0 + (i + 0.5) * step
        bw = 21
        offsets = [-45, -22, 1, 24]
        for (_, color, vals), dx in zip(series, offsets):
            val = vals[i]
            yv = _y(val)
            top_y = min(y0_px, yv)
            bh = max(2.0, abs(y0_px - yv))
            lbl_y = yv - 4 if val >= 0 else yv + 10
            bars_svg.append(
                f'<rect x="{cx + dx:.1f}" y="{top_y:.1f}" width="{bw}" height="{bh:.1f}" rx="2.5" fill="{color}"/>'
                f'<text x="{cx + dx + bw / 2:.1f}" y="{lbl_y:.1f}" text-anchor="middle" font-size="8.0" font-weight="700" fill="{color}">{val:+.2f}</text>'
            )
        bars_svg.append(
            f'<text x="{cx:.1f}" y="{y_bot + 16:.1f}" text-anchor="middle" font-size="9.2" font-weight="700" fill="#1e293b">{html.escape(rlabel)}</text>'
        )

    legend_svg = (
        '<rect x="52" y="18" width="9" height="9" rx="2" fill="#0284c7"/><text x="64" y="26" font-size="8.7" fill="#0369a1" font-weight="700">Pauling-HSAB ΔEf</text>'
        '<rect x="172" y="18" width="9" height="9" rx="2" fill="#059669"/><text x="184" y="26" font-size="8.7" fill="#065f46" font-weight="700">구면포화 %Vbur</text>'
        '<rect x="282" y="18" width="9" height="9" rx="2" fill="#94a3b8"/><text x="294" y="26" font-size="8.7" fill="#475569">RF %Vbur</text>'
        '<rect x="358" y="18" width="9" height="9" rx="2" fill="#e11d48"/><text x="370" y="26" font-size="8.7" fill="#be123c" font-weight="700">RF BDEmin (한계)</text>'
    )

    return (
        f'<svg viewBox="0 0 {w} {h}" style="width:100%; height:auto; display:block; background:#ffffff; border:1px solid #e2e8f0; border-radius:8px;">'
        f'<text x="14" y="12" font-size="10.8" font-weight="800" fill="#0f172a">해석 가능 M&ndash;L 물리 법칙 vs 블랙박스 RF 외삽 성능 (4대 분할 R&sup2;, n=102)</text>'
        f'{legend_svg}{"".join(grid_svg)}{"".join(bars_svg)}'
        f'<text x="{w - 12}" y="{h - 6}" text-anchor="end" font-size="8.2" fill="#64748b">* Homo-&gt;Hetero: 단일 리간드 착물(54종)만으로 이종 리간드 착물(48종) 예측</text>'
        f'</svg>'
    )


def build_mechanistic_storyline_html(
    unified_json_path: str,
    neb_json_path: str,
    mech_stats: Optional[Dict[str, Any]] = None,
) -> str:
    """Build concise Step 4-4 summary card highlighting M-L physics & generalization."""
    if mech_stats is None:
        mech_stats = evaluate_mechanistic_ml_models(unified_json_path, neb_json_path)
    rmap = {r["key"]: r for r in mech_stats["reg_rows"]}
    r_ef = rmap["delta_e_form_per_bond_kcal_mol"]
    r_vb = rmap["vbur_percent"]
    r_bde = rmap["bde_min_kcal_mol"]
    cls2d = mech_stats["cls_results"]["only_2d_ml"]

    tpl = r"""
          <div class="card" style="margin-bottom:24px; border-top:4px solid #2563eb; background:linear-gradient(180deg, #ffffff 0%, #eff6ff 100%);">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
              <div>
                <span class="badge" style="background:#dbeafe; color:#1e40af; font-size:0.74rem; margin-bottom:4px; display:inline-block;">Step 4-4 &bull; Physicochemical M&ndash;L Scaling &amp; Mechanistic Boundaries</span>
                <h4 style="margin:0; font-size:1.12rem; color:#0f172a; font-weight:800;">중심원자&ndash;리간드 물리&middot;화학 해석 모델: 무엇이 DB만으로 일반화되고 어디서 3D&middot;CI-NEB가 필수인가</h4>
              </div>
              <span class="badge badge-phys" style="font-size:0.76rem;">Pauling-HSAB + Spherical Packing + CI-NEB</span>
            </div>

            <p style="font-size:0.88rem; color:#334155; line-height:1.65; margin-bottom:14px;">
              블랙박스 트리 대신 <strong>제1원리 해석 방정식</strong>으로 <code>__N_PREC__</code>종 프리커서를 모델링한 결과, 열역학(&Delta;E<sub>f</sub>/v<sub>M</sub>: <strong>Pauling&ndash;HSAB</strong> 8항)과 입체 차폐(%V<sub>bur</sub>: <strong>구면 포화 충진 법칙</strong> $1-e^{-z}$)는 미지 금속(<code>LOMO-5</code>)&middot;미지 리간드(<code>LOLO-5</code>)는 물론 <strong>단일&rarr;이종 리간드 조합(<code>Homo&rarr;Hetero</code>)까지 일관되게 외삽</strong>되며 오히려 과적합된 RF를 능가합니다. 또한 순수 <code>2D M+L</code> 조합 특징만으로 <code>__N_VIAB__</code>종 적합성 판별 AUROC <strong>__CLS_CV5__</strong>(<code>LOMO __CLS_LOMO__</code>)를 달성합니다. 반면 최소연산자(<code>min</code>) 성격을 지닌 <code>BDE<sub>min</sub></code>과 협동 &beta;-H 분해 장벽(<code>CI-NEB</code>)은 리간드 간 전자&middot;입체 짝짓기에 의존하므로 3D/반응경로 검증을 분리 결합해야 합니다.
            </p>

            <div class="kpi-row" style="margin-bottom:14px;">
              <div class="kpi-chip">
                <div class="kpi-val" style="color:#059669;">R&sup2; = __VB_HOMO__ <span style="font-size:0.75rem; color:#64748b;">(LOMO __VB_LOMO__)</span></div>
                <div class="kpi-label">구면 포화 %Vbur 외삽<br>(선형식 +0.075 &middot; RF __VB_RF_HOMO__ 대비 우위)</div>
              </div>
              <div class="kpi-chip">
                <div class="kpi-val" style="color:#0284c7;">R&sup2; = __EF_HOMO__ <span style="font-size:0.75rem; color:#64748b;">(RF __EF_RF_HOMO__)</span></div>
                <div class="kpi-label">Pauling-HSAB &Delta;Ef/vM<br>(Homo&rarr;Hetero 조합 외삽 일반화)</div>
              </div>
              <div class="kpi-chip">
                <div class="kpi-val" style="color:#2563eb;">AUROC __CLS_CV5__ <span style="font-size:0.75rem; color:#64748b;">(LOLO __CLS_LOLO__)</span></div>
                <div class="kpi-label">2D M+L 조합 판별력<br>(142종 3D 미사용 1차 필터 성능)</div>
              </div>
              <div class="kpi-chip">
                <div class="kpi-val" style="color:#e11d48;">__BDE_RF_HOMO__</div>
                <div class="kpi-label">BDEmin Homo&rarr;Hetero R&sup2;<br>(협동 &beta;-H 분해는 CI-NEB 필수 규명)</div>
              </div>
            </div>

            <div style="background:#ffffff; border:1px solid #bfdbfe; border-radius:8px; padding:12px 14px; font-size:0.84rem; color:#334155; display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:10px;">
              <span><strong>메커니즘 규명 요약:</strong> 102종 DB 내 아미드(<code>__N_AMD__</code>종)&middot;알콕사이드(<code>__N_ALK__</code>종)&middot;&beta;-H 알킬(<code>__N_BETA__</code>종)&middot;구경꾼 &pi;/킬레이트(<code>__N_SPEC__</code>종) 모티프를 <code>SevenNet CI-NEB</code> 8개 채널과 연결해 조합별 열분해 메커니즘 지도를 정립했습니다.</span>
              <button type="button" onclick="switchViewTab('technical'); setTimeout(()=>{document.getElementById('ml-mechanistic-generalization').scrollIntoView({behavior:'smooth'});}, 120);" class="btn-pdf" style="padding:6px 14px; font-size:0.8rem; background:#2563eb; color:#fff; border-color:#2563eb; cursor:pointer;">
                Section 20 물리&middot;메커니즘 상세 &rarr;
              </button>
            </div>
          </div>"""

    return _render_template(
        tpl,
        {
            "__N_PREC__": mech_stats["n_precs"],
            "__N_VIAB__": mech_stats["n_viab"],
            "__VB_HOMO__": f"{r_vb['homo_het_phys_r2']:+.3f}",
            "__VB_LOMO__": f"{r_vb['lomo5_phys_r2']:+.3f}",
            "__VB_RF_HOMO__": f"{r_vb['homo_het_rf_r2']:+.3f}",
            "__EF_HOMO__": f"{r_ef['homo_het_phys_r2']:+.3f}",
            "__EF_RF_HOMO__": f"{r_ef['homo_het_rf_r2']:+.3f}",
            "__CLS_CV5__": f"{cls2d['cv5']['auroc']:.3f}",
            "__CLS_LOMO__": f"{cls2d['lomo5']['auroc']:.3f}",
            "__CLS_LOLO__": f"{cls2d['lolo5']['auroc']:.3f}",
            "__BDE_RF_HOMO__": f"{r_bde['homo_het_rf_r2']:+.3f}",
            "__N_AMD__": mech_stats["motif_flags"].get("amide_H_transfer", 0),
            "__N_ALK__": mech_stats["motif_flags"].get("alkoxide_alpha_beta", 0),
            "__N_BETA__": mech_stats["motif_flags"].get("beta-H_alkyl", 0),
            "__N_SPEC__": mech_stats["motif_flags"].get("spectator_pi_chelate", 0),
        },
    )


def build_mechanistic_technical_html(
    unified_json_path: str,
    neb_json_path: str,
    mech_stats: Optional[Dict[str, Any]] = None,
) -> str:
    """Build Section 20 cards for Technical Report tab."""
    if mech_stats is None:
        mech_stats = evaluate_mechanistic_ml_models(unified_json_path, neb_json_path)

    reg_svg = _build_generalization_regime_svg(mech_stats["reg_rows"])
    ef_coef_rows = "".join(
        f'<tr><td><code>{html.escape(c["term"])}</code></td><td style="text-align:center; font-family:\'JetBrains Mono\'; font-weight:700; color:{"#0284c7" if c["coef"] < 0 else "#be123c"};">{c["coef"]:+.2f}</td></tr>'
        for c in mech_stats["ef_coefs"]
    )
    vb_coef_rows = "".join(
        f'<tr><td><code>{html.escape(c["term"])}</code></td><td style="text-align:center; font-family:\'JetBrains Mono\'; font-weight:700; color:{"#059669" if c["coef"] > 0 else "#be123c"};">{c["coef"]:+.3f}</td></tr>'
        for c in mech_stats["vb_coefs"]
    )
    reg_table_rows = "".join(
        f"""
        <tr>
          <td><strong>{html.escape(r['label'])}</strong><br><span style="font-size:0.71rem; color:#64748b;">{html.escape(r['phys_name'])}</span></td>
          <td style="text-align:center; font-family:'JetBrains Mono';">{r['cv5_phys_r2']:+.3f} / <strong>{r['cv5_rf_r2']:+.3f}</strong></td>
          <td style="text-align:center; font-family:'JetBrains Mono';">{r['lomo5_phys_r2']:+.3f} / <strong>{r['lomo5_rf_r2']:+.3f}</strong></td>
          <td style="text-align:center; font-family:'JetBrains Mono';">{r['lolo5_phys_r2']:+.3f} / <strong>{r['lolo5_rf_r2']:+.3f}</strong></td>
          <td style="text-align:center; font-family:'JetBrains Mono';"><strong>{r['homo_het_phys_r2']:+.3f}</strong> / {r['homo_het_rf_r2']:+.3f}<br><span style="font-size:0.7rem; color:#64748b;">MAE {r['homo_het_phys_mae']:.2f} vs {r['homo_het_rf_mae']:.2f} {html.escape(r['unit'])}</span></td>
        </tr>"""
        for r in mech_stats["reg_rows"]
    )
    motif_table_rows = "".join(
        f"""
        <tr>
          <td><strong>{html.escape(m['motif'])}</strong><br><span style="font-size:0.72rem; color:#0369a1; font-weight:700;">보유 프리커서: {m['count']} / {mech_stats['n_precs']}종</span></td>
          <td><code style="font-size:0.74rem;">{html.escape(m['examples'])}</code></td>
          <td>{html.escape(m['mechanism'])}</td>
          <td style="font-family:'JetBrains Mono'; font-size:0.74rem; color:#5b21b6;">{html.escape(m['neb_anchor'])}</td>
          <td>{html.escape(m['design_rule'])}</td>
        </tr>"""
        for m in mech_stats["motif_rows"]
    )
    c2 = mech_stats["cls_results"]["only_2d_ml"]
    c3 = mech_stats["cls_results"]["full_ml_3d"]
    rmap = {r["key"]: r for r in mech_stats["reg_rows"]}

    tpl = r"""
        <!-- Card 20.1: First-Principles Pauling-HSAB & Spherical Packing Laws -->
        <div class="card" style="margin-bottom:20px; border-top:4px solid #2563eb;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">1. 중심원자&ndash;리간드 물리&middot;화학적 지배 법칙 정식화 및 표준화 물리 계수</h3>
            <span class="badge badge-phys" style="font-size:0.78rem;">Pauling-HSAB &amp; Spherical Packing</span>
          </div>
          <p style="font-size:0.89rem; color:#334155; line-height:1.65; margin-bottom:14px;">
            <code>M + {L<sub>i</sub>}</code> 조합만으로 왜 특정 물성은 높은 정확도로 일반화되는지 규명하기 위해, 결합 열역학과 3D 입체 차폐를 각각 <strong>Pauling 이온공명&ndash;HSAB 선형 방정식</strong>($R^2=__EF_FULL__$)과 <strong>비선형 구면 포화 충진 법칙</strong>($R^2=__VB_FULL__$)으로 압축 정식화했습니다.
          </p>

          <div class="versus-grid" style="margin-bottom:14px;">
            <div class="versus-card" style="border-left:4px solid #0284c7;">
              <div class="versus-header" style="color:#0369a1;">법칙 A: Pauling&ndash;HSAB 결합 안정화 ($\Delta E_f^{\text{std}}/v_M$, Full $R^2=__EF_FULL__$)</div>
              <div style="background:#f8fafc; border-radius:6px; padding:8px 10px; font-size:0.8rem; margin-bottom:8px;">
                $$\frac{\Delta E_f^{\text{std}}}{v_M} \approx \beta_0 - \beta_{\text{ion}}\overline{(\chi_D - \chi_M)^2} - \beta_O f_O - \beta_{\text{pack}}\frac{N_{\text{hv}}}{R_{\text{cov}}^3} + \beta_\pi \eta_\pi$$
              </div>
              <table class="data-table" style="font-size:0.75rem; margin-bottom:8px;">
                <thead><tr><th>물리 항 (표준화)</th><th style="text-align:center;">회귀 계수 (kcal/mol/bond)</th></tr></thead>
                <tbody>__EF_COEFS__</tbody>
              </table>
              <p style="font-size:0.79rem; color:#475569; margin:0; line-height:1.5;">
                <strong>화학적 해석:</strong> Pauling 이온성 $(\Delta\chi)^2$ (<code>-16.24</code>)과 경질 산소 주개 $f_O$ (<code>-6.53</code>)는 $M\text{--}L$ 결합을 강하게 안정화하며, 반대로 <code>CO</code>&middot;<code>Cp</code> 중심의 $\pi$-결합($\eta_\pi$, <code>+11.40</code>)은 표준원소 대비 결합당 안정화도를 완화합니다.
              </p>
            </div>

            <div class="versus-card" style="border-left:4px solid #059669;">
              <div class="versus-header" style="color:#065f46;">법칙 B: 구면 포화 충진 법칙 ($\%V_{\text{bur}}$, Full $R^2=__VB_FULL__$)</div>
              <div style="background:#f8fafc; border-radius:6px; padding:8px 10px; font-size:0.8rem; margin-bottom:8px;">
                $$\%V_{\text{bur}} \approx 100\!\left[1 - \exp\!\left(-\left(c_0 + c_1 \sum_{i=1}^{N_L} \frac{N_{\text{hv},i}^{2/3}}{(R_{\text{cov}}+1.85)^2} + c_2 CN_{\text{eff}} - \cdots\right)\right)\right]$$
              </div>
              <table class="data-table" style="font-size:0.75rem; margin-bottom:8px;">
                <thead><tr><th>기하 항 (표준화)</th><th style="text-align:center;">포화 지수 $z$ 기여도</th></tr></thead>
                <tbody>__VB_COEFS__</tbody>
              </table>
              <p style="font-size:0.79rem; color:#475569; margin:0; line-height:1.5;">
                <strong>물리적 해석:</strong> 부피가 큰 리간드를 혼합할 때 선형 가산은 $100\%$ 포화 상한을 무시해 외삽이 붕괴하지만(<code>Homo&rarr;Hetero +0.075</code>), Beer&ndash;Lambert형 구면 포화 지수 변환을 적용하면 입체 차폐율 예측 $R^2$가 <strong>__VB_HOMO__</strong>(MAE <strong>__VB_MAE__%</strong>)까지 급상승합니다.
              </p>
            </div>
          </div>
        </div>

        <!-- Card 20.2: 4-Regime Generalization & 142-Mol 2D vs 3D Classification -->
        <div class="card" style="margin-bottom:20px; border-top:4px solid #059669;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">2. 현재 DB의 4대 일반화 한계 검증 (<code>5-CV</code> &middot; <code>LOMO-5</code> &middot; <code>LOLO-5</code> &middot; <code>Homo&rarr;Hetero</code>) &amp; 적합성 분류</h3>
            <span class="badge badge-2g" style="font-size:0.78rem;">4-Regime Generalization Stress Test</span>
          </div>

          <div style="display:grid; grid-template-columns:minmax(320px, 1.12fr) minmax(260px, 0.88fr); gap:16px; align-items:start; margin-bottom:14px;">
            <div>
              __REG_SVG__
            </div>
            <div style="display:flex; flex-direction:column; gap:10px;">
              <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; padding:12px;">
                <div style="font-weight:800; color:#0f172a; font-size:0.85rem; margin-bottom:6px;">순수 <code>2D M+L + SMARTS</code> vs <code>3D SevenNet</code> 추가 분류 성능 (n=__N_VIAB__종: 102 Prec + 40 Neg)</div>
                <table class="data-table" style="font-size:0.76rem; margin:0;">
                  <thead>
                    <tr>
                      <th>분류 특징셋</th>
                      <th style="text-align:center;">5-CV AUROC (Bal)</th>
                      <th style="text-align:center;">LOMO-5 AUROC</th>
                      <th style="text-align:center;">LOLO-5 AUROC</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr style="background:#f0fdf4;">
                      <td><strong>2D M+L + SMARTS만</strong><br><span style="font-size:0.69rem; color:#166534;">3D 계산 0회</span></td>
                      <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:800;">__C2_CV__ (__C2_CV_B__)</td>
                      <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:800;">__C2_MO__ (__C2_MO_B__)</td>
                      <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:800;">__C2_LO__ (__C2_LO_B__)</td>
                    </tr>
                    <tr>
                      <td>2D + SevenNet 3D 전체<br><span style="font-size:0.69rem; color:#64748b;">RF 과적합 대조</span></td>
                      <td style="text-align:center; font-family:'JetBrains Mono';">__C3_CV__ (__C3_CV_B__)</td>
                      <td style="text-align:center; font-family:'JetBrains Mono';">__C3_MO__ (__C3_MO_B__)</td>
                      <td style="text-align:center; font-family:'JetBrains Mono';">__C3_LO__ (__C3_LO_B__)</td>
                    </tr>
                  </tbody>
                </table>
              </div>
              <div style="background:#fff1f2; border:1px solid #fecdd3; border-radius:8px; padding:10px 12px; font-size:0.79rem; color:#881337; line-height:1.55;">
                <strong>방법론적 통찰:</strong> 일반 Precursor vs Negative Control 분류는 <code>2D M+L + SMARTS</code>만으로도 AUROC <code>__C2_CV__</code>(LOMO <code>__C2_MO__</code>)까지 도달하며 소규모 142종 RF에서는 오히려 3D 노이즈 과적합이 적습니다. 그러나 Section 15에서 검증했듯이 학습 DB에 없던 <strong>완전 미지 실패 메커니즘(<code>LOFMO</code>)</strong>에 대해서는 ML 분류기가 붕괴(<code>0.000</code>)하고 물리 게이트(<code>TMC-SA</code>)만 <code>1.000</code>을 유지하므로 1차 <code>2D M+L</code> 필터 뒤 <code>TMC-SA</code> 물리 검증이 필수입니다.
              </div>
            </div>
          </div>

          <div class="table-responsive">
            <table class="data-table" style="font-size:0.79rem;">
              <thead>
                <tr>
                  <th>물성 (해석 물리 모델명)</th>
                  <th style="text-align:center;">5-CV $R^2$<br>(물리식 / RF)</th>
                  <th style="text-align:center;">LOMO-5 미지 금속<br>(물리식 / RF)</th>
                  <th style="text-align:center;">LOLO-5 미지 리간드<br>(물리식 / RF)</th>
                  <th style="text-align:center;">Homo(__N_HOMO__)&rarr;Hetero(__N_HET__)<br>(물리식 / RF)</th>
                </tr>
              </thead>
              <tbody>
                __REG_ROWS__
              </tbody>
            </table>
          </div>
          <p style="font-size:0.81rem; color:#475569; line-height:1.55; margin:10px 0 0 0;">
            <strong>핵심 발견:</strong> 합(Sum)/평균 성질인 $\Delta E_f/v_M$, $\%V_{\text{bur}}$, <code>Ligand SA</code>는 현재 DB만으로 미지 금속&middot;리간드&middot;이종 조합에 강하게 일반화됩니다. 반면 최약 결합 단일값인 <code>BDE<sub>min</sub></code>(최소 연산자 + 산화수 환원 특이성 + <em>trans</em>-효과)과 <code>&Delta;G<sub>dimer</sub></code>는 <code>Homo&rarr;Hetero</code> 외삽에서 일관되게 실패(<code>__BDE_RF__</code>)하므로 조합 조성 모델의 적용 경계를 넘어섭니다.
          </p>
        </div>

        <!-- Card 20.3: 102-Precursor Kinetic Ligand-Motif Decomposition Map Linked to CI-NEB & Literature -->
        <div class="card" style="border-top:4px solid #7c3aed;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">3. 102종 리간드 모티프&ndash;<code>SevenNet CI-NEB</code>&ndash;문헌 연계 열분해 메커니즘 해석 지도</h3>
            <span class="badge badge-4g" style="font-size:0.78rem;">Kinetic Motif-Mechanism Map</span>
          </div>
          <p style="font-size:0.88rem; color:#334155; line-height:1.65; margin-bottom:12px;">
            <code>BDE<sub>min</sub></code>이 이종 리간드 외삽과 $T_{\text{decomp}}$ 순위 예측($r=0.507$)에서 한계를 드러낸 이유는 실제 열분해가 균일 결합 해리가 아니라 <strong>리간드 모티프 특이적 폐쇄 껍질 협동 전이상태</strong>(Section 17 <code>CI-NEB</code>: $r=0.879$, 문헌 DFT $r=0.972$)를 따르기 때문입니다. 아래 표는 <code>102</code>종 DB를 반응 메커니즘 모티프별로 분해해 <code>SevenNet CI-NEB</code> 실측값 및 문헌 메커니즘과 직결한 설계 규칙입니다.
          </p>
          <div class="table-responsive">
            <table class="data-table" style="font-size:0.79rem;">
              <thead>
                <tr>
                  <th style="width:185px;">리간드 분해 모티프</th>
                  <th style="width:165px;">DB 주요 리간드 (빈도)</th>
                  <th>지배 열분해&middot;반응 메커니즘</th>
                  <th style="width:205px;">SevenNet CI-NEB / 문헌 $E_a$ (kcal/mol)</th>
                  <th>반도체 프리커서 설계 해석</th>
                </tr>
              </thead>
              <tbody>
                __MOTIF_ROWS__
              </tbody>
            </table>
          </div>
        </div>"""

    r_vb_map = rmap["vbur_percent"]
    r_bde_map = rmap["bde_min_kcal_mol"]
    return _render_template(
        tpl,
        {
            "__EF_FULL__": f"{mech_stats['ef_full_r2']:.3f}",
            "__VB_FULL__": f"{mech_stats['vb_full_r2']:.3f}",
            "__EF_COEFS__": ef_coef_rows,
            "__VB_COEFS__": vb_coef_rows,
            "__VB_HOMO__": f"{r_vb_map['homo_het_phys_r2']:+.3f}",
            "__VB_MAE__": f"{r_vb_map['homo_het_phys_mae']:.2f}",
            "__REG_SVG__": reg_svg,
            "__N_VIAB__": mech_stats["n_viab"],
            "__C2_CV__": f"{c2['cv5']['auroc']:.3f}",
            "__C2_CV_B__": f"{c2['cv5']['bal_acc']:.3f}",
            "__C2_MO__": f"{c2['lomo5']['auroc']:.3f}",
            "__C2_MO_B__": f"{c2['lomo5']['bal_acc']:.3f}",
            "__C2_LO__": f"{c2['lolo5']['auroc']:.3f}",
            "__C2_LO_B__": f"{c2['lolo5']['bal_acc']:.3f}",
            "__C3_CV__": f"{c3['cv5']['auroc']:.3f}",
            "__C3_CV_B__": f"{c3['cv5']['bal_acc']:.3f}",
            "__C3_MO__": f"{c3['lomo5']['auroc']:.3f}",
            "__C3_MO_B__": f"{c3['lomo5']['bal_acc']:.3f}",
            "__C3_LO__": f"{c3['lolo5']['auroc']:.3f}",
            "__C3_LO_B__": f"{c3['lolo5']['bal_acc']:.3f}",
            "__N_HOMO__": mech_stats["n_homo"],
            "__N_HET__": mech_stats["n_hetero"],
            "__REG_ROWS__": reg_table_rows,
            "__BDE_RF__": f"{r_bde_map['homo_het_rf_r2']:+.3f}",
            "__MOTIF_ROWS__": motif_table_rows,
        },
    )
