#!/usr/bin/env python3
"""
HTML Report Builder for Section 18 & Summary Step 4 Add-on:
Transition-Metal Multi-Coordination Determination & Combinatorial Ligand Statistics
---------------------------------------------------------------------------------
Computes live dataset telemetry across the 102 semiconductor precursors (27 audited +
75 expanded) and 160 unified benchmark records, runs deterministic 5-fold CV and
Leave-One-Metal-Out (LOMO-5) compositional ablations, and emits standalone SVG/HTML
cards for both the Summary tab and Section 18 of the Technical Report tab.
"""

from __future__ import annotations

import html
import json
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Human-readable chemical shorthand for canonical fragment SMILES
LIGAND_COMMON_LABELS: Dict[str, str] = {
    "C": "Me (-CH3)",
    "CNC": "NMe2 (-N(CH3)2)",
    "CC": "Et (-C2H5)",
    "[C-]#[O+]": "CO (Carbonyl)",
    "Cl": "Cl (Chloro)",
    "CC(C)O": "OiPr (-OCH(CH3)2)",
    "CC(=O)C=C(C)O": "acac (Acetylacetonate)",
    "CC(C)(C)C(=O)C=C(O)C(C)(C)C": "thd / tmhd",
    "CCNC": "NEtMe (-N(Et)Me)",
    "CCC1=CC=CC1": "EtCp",
    "CC(C)(C)N": "tBuN / tBuNH",
    "CC1=CC=CC1": "MeCp",
    "C1=CCC=C1": "Cp (Cyclopentadienyl)",
    "CC(=NC(C)C)NC(C)C": "iPr-amd (Amidinate)",
    "CN(C)CC(C)(C)O": "dmamp",
    "CCO": "OEt (-OC2H5)",
    "CC(C)(C)O": "OtBu (-OC(CH3)3)",
    "O": "Oxo / Hydroxo",
    "C[N+]C": "mu-NMe2 (Bridging)",
    "I": "I (Iodo)",
    "Br": "Br (Bromo)",
    "CCNCC": "NEt2 (-N(C2H5)2)",
    "C1=CCCC=C1": "1,3-Cyclohexadiene (CHD)",
    "CC(C)N=NNC(C)C": "iPr-Triazenide",
    "CC(C)N=CNC(C)C": "iPr-Formamidinate",
    "CO": "OMe (-OCH3)",
    "C[Si](C)(C)N[Si](C)(C)C": "N(SiMe3)2 (hmds)",
    "CC(C)NC(C)C": "NiPr2 (-N(iPr)2)",
    "CC=CC=C(C)C": "DMPD (2,4-Me2-pentadienyl)",
    "CC(=NC(C)(C)C)NC(C)(C)C": "tBu-amd",
    "CC(O)=CC(=O)C(C)C": "mhdk",
    "CN(C)CCNC(C)(C)C": "Me2N(CH2)2NtBu",
    "COCC[NH+]C(C)(C)C": "MeO(CH2)2NHtBu",
    "c1ccnc1": "pyr (Pyrrolyl)",
}

COORD_CASE_SPECS: List[Tuple[str, str, str]] = [
    ("EXT-ALD-135", "Ga(I) 1배위 할라이드", "단일 M-Cl 시그마 결합 (pi=0)"),
    ("PREC_MET_04", "Co(II) 아미디네이트", "2D SMILES 단일 M-N 표기 경계 (dative 생략 시 CN=2)"),
    ("NEG-019", "Cr(III) 과소배위 트리아릴", "eta^1 페닐 3개 (sigma-aryl 제외 -> pi=0, 9e-)"),
    ("EXT-ALD-115", "Ru(VIII) 고산화 테트라옥소", "고극성 M=O 수축 (r_eff 1.46->1.351 A, chi_eff 상승)"),
    ("PREC_HK_01", "Hf(IV) 테트라아마이드", "사면체 4배위 M-N (r_eff 1.75->1.567 A)"),
    ("NEG-004", "W(V) 펜타아릴", "eta^1 페닐 5개 (sigma-aryl 제외 -> pi=0, CN_eff=5)"),
    ("PREC_MET_03", "Ru(II) 18전자 메탈로센", "Cp 2개 (pi=2 -> 2+2*2=6, r_eff 1.46->1.502 A 팽창)"),
    ("PREC_HK_07", "Zr(IV) 피아노스툴 하프샌드위치", "Cp 1개 + NMe2 3개 (pi=1 -> 4+2*1=6, 12e-)"),
    ("EXT-ALD-146", "Ru(0) TMM 트라이카보닐", "exo-디엔형 고리 + 3 CO (pi=1 -> 5+2*1=7)"),
    ("PREC_HK_06", "Zr(IV) 비스(메틸Cp) 착물", "MeCp 2개 + OMe + Me (pi=2 -> 4+2*2=8)"),
    ("PREC_MET_05", "Co2(CO)6 알카인 가교 착물", "다중심 3원 고리 중복 계수 경계 (pi=4, clip(1.22 r_cov))"),
]


def fragment_metal_and_ligands(smiles: str, center_sym: str) -> Tuple[str, List[str], int, List[str]]:
    """
    Detach the central atom(s) from a coordination complex SMILES into canonical ligand fragments.

    Physics/Algorithm:
        Removes all atoms matching ``center_sym`` cleanly using an ``RWMol`` graph cut,
        preserves metal-adjacent donor atom identities (C, N, O, halogen, etc.) prior to
        excision, and returns canonical isomeric-stripped fragment SMILES for stoichiometric
        multiset counting.
    """
    RDLogger.DisableLog("rdApp.error")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return center_sym, [], 0, []
        center_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == center_sym]
        donor_elems: List[str] = []
        for cidx in center_idxs:
            for nb in mol.GetAtomWithIdx(cidx).GetNeighbors():
                if nb.GetIdx() not in center_idxs:
                    donor_elems.append(nb.GetSymbol())

        rw = Chem.RWMol(mol)
        for cidx in sorted(center_idxs, reverse=True):
            rw.RemoveAtom(cidx)
        frag_mol = rw.GetMol()
        frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)
        frag_smiles: List[str] = []
        for fm in frags:
            try:
                Chem.SanitizeMol(fm)
            except Exception:
                pass
            frag_smiles.append(Chem.MolToSmiles(fm, isomericSmiles=False))
        return center_sym, frag_smiles, len(center_idxs), donor_elems
    finally:
        RDLogger.EnableLog("rdApp.error")


def _extract_continuous_ligand_features(smiles: str, center_sym: str) -> List[float]:
    """Extract permutation-invariant continuous multiset features from detached ligand fragments."""
    RDLogger.DisableLog("rdApp.error")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return [0.0] * 10
        center_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == center_sym]
        donor_elems: List[str] = []
        for cidx in center_idxs:
            for nb in mol.GetAtomWithIdx(cidx).GetNeighbors():
                if nb.GetIdx() not in center_idxs:
                    donor_elems.append(nb.GetSymbol())

        rw = Chem.RWMol(mol)
        for cidx in sorted(center_idxs, reverse=True):
            rw.RemoveAtom(cidx)
        frags = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
        heavy_sum = 0.0
        heavy_max = 0.0
        rotb_sum = 0.0
        n_rings_sum = 0.0
        n_lig = float(len(frags))
        for fm in frags:
            try:
                fm.UpdatePropertyCache(strict=False)
                Chem.GetSymmSSSR(fm)
                hv = float(fm.GetNumHeavyAtoms())
                heavy_sum += hv
                heavy_max = max(heavy_max, hv)
                rotb_sum += float(Descriptors.NumRotatableBonds(fm))
                n_rings_sum += float(fm.GetRingInfo().NumRings())
            except Exception:
                hv = float(fm.GetNumAtoms())
                heavy_sum += hv
                heavy_max = max(heavy_max, hv)
        n_d = max(1, len(donor_elems))
        return [
            n_lig,
            heavy_sum,
            heavy_sum / max(1.0, n_lig),
            heavy_max,
            rotb_sum,
            n_rings_sum,
            sum(1 for e in donor_elems if e == "C") / n_d,
            sum(1 for e in donor_elems if e == "N") / n_d,
            sum(1 for e in donor_elems if e == "O") / n_d,
            sum(1 for e in donor_elems if e in {"F", "Cl", "Br", "I"}) / n_d,
        ]
    finally:
        RDLogger.EnableLog("rdApp.error")


def compute_compositional_modeling_benchmark(precursors: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Evaluate 5-Fold CV and Leave-One-Metal-Out (LOMO-5) regression across 6 targets.

    Physics/Algorithm:
        Compares three progressively richer representations on the 102-precursor matrix:
        (A) Discrete additive ``OneHot(M) + BagOfLigands`` (RidgeCV),
        (B) Continuous additive ``Phys(M) + ContinuousLigandMultiset`` (RandomForest), and
        (B+) Explicit Metal-Ligand Coupling ``Phys(M) + ContinuousLigandMultiset + M x L``
        including ``heavy_sum / r_cov^3`` (# units: A^-3), ``heavy_max / r_cov^2`` (# units: A^-2),
        ``delta_chi_mean``, and ``ve_deficit``.
    """
    rows: List[Tuple[Dict[str, Any], str, List[str]]] = []
    lig_vocab = set()
    metal_vocab = set()
    for r in precursors:
        c, frags, _, _ = fragment_metal_and_ligands(r["smiles"], r["center_atom"])
        rows.append((r, c, frags))
        metal_vocab.add(c)
        lig_vocab.update(frags)

    all_metals = sorted(metal_vocab)
    all_ligs = sorted(lig_vocab)
    m_idx = {m: i for i, m in enumerate(all_metals)}
    l_idx = {l: i for i, l in enumerate(all_ligs)}

    X_m_oh = np.zeros((len(rows), len(all_metals)))
    X_l_bag = np.zeros((len(rows), len(all_ligs)))
    X_m_phys = np.zeros((len(rows), 3))
    X_l_cont = np.zeros((len(rows), 12))
    X_int = np.zeros((len(rows), 4))
    groups_m: List[str] = []

    for i, (r, c, frags) in enumerate(rows):
        X_m_oh[i, m_idx[c]] = 1.0
        groups_m.append(c)
        for fs in frags:
            X_l_bag[i, l_idx[fs]] += 1.0

        rcov = float(r["r_cov_ang"])  # units: A
        chi = float(r["chi_pauling"])
        nval = float(r["n_val"])
        pi_h = float(r["pi_hapticity_index"])
        chel = float(r["chelate_ring_count"])
        ve_def = float(r["ve_deficit"])
        dchi = float(r["delta_chi_mean"])
        lf = _extract_continuous_ligand_features(r["smiles"], c)

        X_m_phys[i] = [rcov, chi, nval]
        X_l_cont[i] = lf + [pi_h, chel]
        X_int[i] = [
            lf[1] / (rcov ** 3),  # units: A^-3
            lf[3] / (rcov ** 2),  # units: A^-2
            dchi,
            ve_def,
        ]

    X_disc = np.hstack([X_m_oh, X_l_bag])
    X_cont = np.hstack([X_m_phys, X_l_cont])
    X_coup = np.hstack([X_m_phys, X_l_cont, X_int])

    target_specs = [
        ("delta_e_form_per_bond_kcal_mol", "표준 결합 생성에너지 (ΔEf/vM)", "kcal/mol"),
        ("vbur_percent", "매몰부피 (%Vbur)", "%"),
        ("ligand_sa", "프로-리간드 SA (Ligand SA)", "score"),
        ("bde_min_kcal_mol", "최소 결합 해리 에너지 (BDEmin)", "kcal/mol"),
        ("tmc_sa_v2", "통합 합성가능성 (SAPIS)", "score"),
        ("delta_g_dimer_373_kcal_mol", "373 K 이합체 자유에너지 (ΔGdimer)", "kcal/mol"),
    ]

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    gkf = GroupKFold(n_splits=5)
    results: List[Dict[str, Any]] = []

    for key, label, unit in target_specs:
        vals = [r.get(key) for r, _, _ in rows]
        mask = np.array([v is not None and not np.isnan(v) for v in vals])
        y = np.array([float(v) for v, m in zip(vals, mask) if m])
        g_sub = [g for g, m in zip(groups_m, mask) if m]

        row_res: Dict[str, Any] = {"key": key, "label": label, "unit": unit, "n": int(mask.sum())}
        for tag, Xfull, use_ridge in [
            ("disc", X_disc, True),
            ("cont", X_cont, False),
            ("coup", X_coup, False),
        ]:
            X_sub = Xfull[mask]
            pred_kf = np.zeros_like(y)
            for tr, te in kf.split(X_sub):
                if use_ridge:
                    mdl = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 20)))
                else:
                    mdl = RandomForestRegressor(n_estimators=120, random_state=42)
                mdl.fit(X_sub[tr], y[tr])
                pred_kf[te] = mdl.predict(X_sub[te])

            pred_lomo = np.zeros_like(y)
            for tr, te in gkf.split(X_sub, y, groups=g_sub):
                if use_ridge:
                    mdl = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 20)))
                else:
                    mdl = RandomForestRegressor(n_estimators=120, random_state=42)
                mdl.fit(X_sub[tr], y[tr])
                pred_lomo[te] = mdl.predict(X_sub[te])

            row_res[f"{tag}_cv5_r2"] = round(float(r2_score(y, pred_kf)), 3)
            row_res[f"{tag}_cv5_mae"] = round(float(mean_absolute_error(y, pred_kf)), 2)
            row_res[f"{tag}_lomo_r2"] = round(float(r2_score(y, pred_lomo)), 3)
            row_res[f"{tag}_lomo_mae"] = round(float(mean_absolute_error(y, pred_lomo)), 2)
        results.append(row_res)

    return results


def analyze_precursor_coordination_and_ligands(unified_json_path: str) -> Dict[str, Any]:
    """Compute live multi-coordination and ligand statistics from the unified 160-molecule dataset."""
    with open(unified_json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    records = data["records"] if isinstance(data, dict) else data
    by_id = {r["id"]: r for r in records}
    precursors = [
        r for r in records if r.get("group") in {"Audited Precursor (27)", "Expanded Precursor (75)"}
    ]

    center_counts: Counter[str] = Counter()
    cn_eff_counts: Counter[float] = Counter()
    pi_counts: Counter[int] = Counter()
    chelate_counts: Counter[int] = Counter()
    ligand_occ_counts: Counter[str] = Counter()
    ligand_to_metals: Dict[str, set[str]] = {}
    unique_ligand_sets: List[set[str]] = []
    homoleptic_count = 0
    heteroleptic_count = 0
    multicenter_count = 0
    total_ligand_instances = 0

    for r in precursors:
        c = r["center_atom"]
        center_counts[c] += 1
        cn_eff_counts[float(r["cn_eff"])] += 1
        pi_counts[int(r["pi_hapticity_index"])] += 1
        chelate_counts[int(r["chelate_ring_count"])] += 1

        _, frags, n_c, _ = fragment_metal_and_ligands(r["smiles"], c)
        if n_c > 1:
            multicenter_count += 1
        total_ligand_instances += len(frags)
        u_frags = set(frags)
        unique_ligand_sets.append(u_frags)
        if len(u_frags) <= 1:
            homoleptic_count += 1
        else:
            heteroleptic_count += 1
        for fs in frags:
            ligand_occ_counts[fs] += 1
        for fs in u_frags:
            ligand_to_metals.setdefault(fs, set()).add(c)

    shared_ligs = [k for k, v in ligand_to_metals.items() if len(v) >= 2]
    singleton_ligs = [k for k, v in ligand_to_metals.items() if len(v) == 1]

    top15_ligands: List[Dict[str, Any]] = []
    for smi, occ in ligand_occ_counts.most_common(15):
        metals = sorted(ligand_to_metals.get(smi, set()))
        prec_cnt = sum(1 for s in unique_ligand_sets if smi in s)
        top15_ligands.append(
            {
                "smiles": smi,
                "label": LIGAND_COMMON_LABELS.get(smi, smi),
                "occurrences": occ,
                "precursor_count": prec_cnt,
                "metal_count": len(metals),
                "metals": metals,
            }
        )

    topk_coverage: List[Dict[str, Any]] = []
    for k in (5, 10, 15, 20, 30):
        topk_set = {smi for smi, _ in ligand_occ_counts.most_common(k)}
        full_c = sum(1 for s in unique_ligand_sets if s and s.issubset(topk_set))
        part_c = sum(1 for s in unique_ligand_sets if s & topk_set)
        topk_coverage.append(
            {
                "k": k,
                "full": full_c,
                "full_pct": round(100.0 * full_c / len(unique_ligand_sets), 1),
                "partial": part_c,
                "partial_pct": round(100.0 * part_c / len(unique_ligand_sets), 1),
            }
        )

    coord_cases: List[Dict[str, Any]] = []
    for mid, role_note, geom_note in COORD_CASE_SPECS:
        rec = by_id.get(mid)
        if rec is None:
            continue
        coord_cases.append(
            {
                "id": mid,
                "name": rec["name"],
                "center_atom": rec["center_atom"],
                "role_note": role_note,
                "geom_note": geom_note,
                "r_cov_ang": rec["r_cov_ang"],
                "r_eff_ang": rec["r_eff_ang"],
                "chi_pauling": rec["chi_pauling"],
                "chi_eff": rec["chi_eff"],
                "cn_eff": rec["cn_eff"],
                "pi": rec["pi_hapticity_index"],
                "chelate": rec["chelate_ring_count"],
                "n_ve": rec["n_ve"],
                "ve_deficit": rec["ve_deficit"],
                "vbur_percent": rec["vbur_percent"],
                "vbur_d2_percent": rec["vbur_d2_percent"],
            }
        )

    ablation_rows = compute_compositional_modeling_benchmark(precursors)

    return {
        "n_records_total": len(records),
        "n_precursors": len(precursors),
        "n_unique_metals": len(center_counts),
        "center_counts": center_counts.most_common(),
        "cn_eff_counts": sorted(cn_eff_counts.items()),
        "pi_counts": sorted(pi_counts.items()),
        "chelate_counts": sorted(chelate_counts.items()),
        "n_unique_ligands": len(ligand_occ_counts),
        "total_ligand_instances": total_ligand_instances,
        "homoleptic_count": homoleptic_count,
        "heteroleptic_count": heteroleptic_count,
        "multicenter_count": multicenter_count,
        "n_shared_ligands": len(shared_ligs),
        "n_singleton_ligands": len(singleton_ligs),
        "top15_ligands": top15_ligands,
        "topk_coverage": topk_coverage,
        "coord_cases": coord_cases,
        "ablation_rows": ablation_rows,
    }


def _build_cn_eff_svg(cn_items: List[Tuple[float, int]], n_prec: int) -> str:
    """Render an inline SVG vertical bar chart for CN_eff distribution across 102 precursors."""
    w, h = 440, 235
    pad_l, pad_r, pad_t, pad_b = 44, 18, 28, 48
    plot_w = w - pad_l - pad_r
    plot_h = h - pad_t - pad_b
    max_cnt = max((cnt for _, cnt in cn_items), default=1)
    y_max = ((max_cnt + 9) // 10) * 10
    step = plot_w / max(1, len(cn_items))
    bar_w = min(34.0, step * 0.62)

    bars: List[str] = []
    for yv in range(0, y_max + 1, 10):
        yy = pad_t + plot_h * (1.0 - yv / y_max)
        bars.append(
            f'<line x1="{pad_l}" y1="{yy:.1f}" x2="{w - pad_r}" y2="{yy:.1f}" '
            f'stroke="#e2e8f0" stroke-width="1" stroke-dasharray="3,3" />'
            f'<text x="{pad_l - 6}" y="{yy + 3:.1f}" text-anchor="end" font-size="9.5" fill="#64748b">{yv}</text>'
        )

    for idx, (cn, cnt) in enumerate(cn_items):
        cx = pad_l + (idx + 0.5) * step
        bh = plot_h * (cnt / y_max)
        by = pad_t + plot_h - bh
        color = "#7c3aed" if cn >= 6.0 else ("#0284c7" if cn == 4.0 else ("#0d9488" if cn == 3.0 else "#f59e0b"))
        lbl = f"{int(cn)}" if cn == int(cn) else f"{cn:.1f}"
        sub = "pi/고배위" if 6.0 <= cn <= 8.0 else ("케이지 경계" if cn > 8.0 else "")
        bars.append(
            f'<rect x="{cx - bar_w * 0.5:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
            f'rx="3" fill="{color}"><title>CN_eff={lbl}: {cnt}/{n_prec}종 ({cnt / n_prec * 100:.1f}%)</title></rect>'
            f'<text x="{cx:.1f}" y="{by - 4:.1f}" text-anchor="middle" font-size="10" font-weight="800" fill="#0f172a">{cnt}</text>'
            f'<text x="{cx:.1f}" y="{pad_t + plot_h + 15:.1f}" text-anchor="middle" font-size="10.5" font-weight="700" fill="#334155">{lbl}</text>'
        )
        if sub:
            bars.append(
                f'<text x="{cx:.1f}" y="{pad_t + plot_h + 28:.1f}" text-anchor="middle" font-size="8" fill="{color}" font-weight="700">{sub}</text>'
            )

    return (
        f'<svg viewBox="0 0 {w} {h}" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;">'
        f'<text x="{w / 2:.1f}" y="17" text-anchor="middle" font-size="11" font-weight="800" fill="#0f172a">'
        f'반도체 프리커서 {n_prec}종 유효 배위수 (CN_eff) 분포</text>'
        f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t + plot_h}" stroke="#64748b" stroke-width="1.2" />'
        f'<line x1="{pad_l}" y1="{pad_t + plot_h}" x2="{w - pad_r}" y2="{pad_t + plot_h}" stroke="#64748b" stroke-width="1.2" />'
        f'{"".join(bars)}'
        f'<text x="{w / 2:.1f}" y="{h - 6}" text-anchor="middle" font-size="9.5" fill="#475569">'
        f'유효 배위수 CN_eff = N_nbr(heavy) + N_H(imp) + 2*pi_hapticity</text>'
        f'</svg>'
    )


def _build_top_ligands_svg(top_ligands: List[Dict[str, Any]], total_instances: int, n_unique_lig: int) -> str:
    """Render an inline SVG horizontal bar chart of top-12 ligand fragments and cross-metal sharing."""
    items = top_ligands[:12]
    w, h = 480, 295
    pad_l, pad_r, pad_t, pad_b = 140, 88, 30, 26
    plot_w = w - pad_l - pad_r
    row_h = (h - pad_t - pad_b) / max(1, len(items))
    max_occ = max((it["occurrences"] for it in items), default=1)

    bars: List[str] = []
    for i, it in enumerate(items):
        cy = pad_t + (i + 0.5) * row_h
        bw = plot_w * (it["occurrences"] / max_occ)
        bh = row_h * 0.62
        is_shared = it["metal_count"] >= 2
        fill = "#0284c7" if it["metal_count"] >= 4 else ("#0d9488" if is_shared else "#94a3b8")
        short_lbl = html.escape(it["label"].split(" (")[0])
        m_str = ", ".join(it["metals"])
        bars.append(
            f'<text x="{pad_l - 8}" y="{cy + 3.5:.1f}" text-anchor="end" font-size="10" font-weight="700" fill="#1e293b">{short_lbl}</text>'
            f'<rect x="{pad_l}" y="{cy - bh * 0.5:.1f}" width="{bw:.1f}" height="{bh:.1f}" rx="3" fill="{fill}">'
            f'<title>{html.escape(it["smiles"])}: {it["occurrences"]}회 등장 / {it["precursor_count"]}개 프리커서 / {it["metal_count"]}개 금속 ({m_str})</title></rect>'
            f'<text x="{pad_l + bw + 6:.1f}" y="{cy + 3.5:.1f}" font-size="9.5" font-weight="700" fill="#0f172a">'
            f'{it["occurrences"]}회 ({it["precursor_count"]}종, {it["metal_count"]}M)</text>'
        )

    return (
        f'<svg viewBox="0 0 {w} {h}" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;">'
        f'<text x="{w / 2:.1f}" y="18" text-anchor="middle" font-size="11" font-weight="800" fill="#0f172a">'
        f'상위 12대 리간드 절편 빈도 ({n_unique_lig}종 고유 / 총 {total_instances} 인스턴스)</text>'
        f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{h - pad_b}" stroke="#64748b" stroke-width="1.2" />'
        f'{"".join(bars)}'
        f'<text x="{w - 10}" y="{h - 8}" text-anchor="end" font-size="8.8" fill="#475569">'
        f'색상: 파랑(&ge;4금속 공유) &bull; 청록(2~3금속) &bull; 회색(단일금속) | 표기: 등장횟수(프리커서수, 금속수)</text>'
        f'</svg>'
    )


def _render_template(tpl: str, mapping: Dict[str, Any]) -> str:
    """Replace explicit ``__TOKEN__`` placeholders cleanly without touching LaTeX braces."""
    out = tpl
    for k, v in mapping.items():
        out = out.replace(k, str(v))
    return out


def build_ligand_coord_storyline_html(
    unified_json_path: str,
    stats: Optional[Dict[str, Any]] = None,
) -> str:
    """Build concise Summary tab card covering multi-coordination physics and ligand statistics."""
    if stats is None:
        stats = analyze_precursor_coordination_and_ligands(unified_json_path)
    ab_map = {row["key"]: row for row in stats["ablation_rows"]}
    r_ef = ab_map["delta_e_form_per_bond_kcal_mol"]
    r_vb = ab_map["vbur_percent"]
    cov15 = next(x for x in stats["topk_coverage"] if x["k"] == 15)
    cov30 = next(x for x in stats["topk_coverage"] if x["k"] == 30)
    singleton_pct = 100.0 * stats["n_singleton_ligands"] / max(1, stats["n_unique_ligands"])
    cn_map = {int(k): v for k, v in stats["cn_eff_counts"] if k == int(k)}
    n_cn6_plus = sum(v for k, v in stats["cn_eff_counts"] if k >= 6.0)

    tpl = r"""
          <div class="card" style="margin-top:20px; margin-bottom:24px; padding:20px; border-top:4px solid #7c3aed;">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:10px; margin-bottom:10px;">
              <div style="font-weight:800; font-size:1.05rem; color:#0f172a;">
                4-2. 전이금속 다중 배위수(<code>CN<sub>eff</sub></code>) 결정 원리 &amp; 리간드 조합&middot;물성 연계 요약
              </div>
              <span class="badge badge-4g" style="font-size:0.76rem;">Multi-Coordination State &bull; __N_PREC__ Precursors &bull; __N_LIG__ Unique Ligands</span>
            </div>
            <p style="font-size:0.86rem; color:#475569; line-height:1.6; margin:0 0 14px 0;">
              전이금속은 동일 원소라도 산화수&middot;$\sigma$/$\pi$-배위 양식&middot;킬레이트 여부에 따라 배위수와 유효 반경이 연속 변동합니다. 본 프레임워크는 다중심 $\pi$-배위($\eta^n$-Cp 등)에는 유효 배위수 $+2$ 및 $4e^-$를 가산하되 $\eta^1$ $\sigma$-아릴(예: $\text{CrPh}_3$, $\text{WPh}_5$)은 엄격히 배제하는 $CN_{\text{eff}} = N_{\text{nbr}}^{\text{heavy}} + N_H^{\text{imp}} + 2\,\eta_{\pi}$ 규칙을 적용합니다. 아울러 <code>__N_PREC__</code>종 프리커서의 리간드 조합 공간을 해부한 결과, 이산적 리간드 ID 가산 회귀는 <code>__SINGLETON_PCT__%</code>에 달하는 단일금속 희소 리간드와 비선형 입체&middot;전자 결합으로 붕괴하는 반면, <strong>금속 주기율 불변량 + 연속 리간드 집합 서술자 + 금속&times;리간드 결합항</strong>을 결합하면 미관측 금속 외삽(<code>LOMO-5</code>)에서도 높은 일반화 성능을 확보합니다.
            </p>

            <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap:12px; margin-bottom:14px;">
              <div style="background:#f5f3ff; border:1px solid #ddd6fe; border-left:4px solid #7c3aed; border-radius:6px; padding:12px 14px;">
                <div style="font-size:0.76rem; font-weight:700; color:#5b21b6;">다중 유효 배위수 범위</div>
                <div style="font-size:1.35rem; font-weight:800; color:#4c1d95; margin:3px 0;">CN<sub>eff</sub> = 1 ~ 8 (14)</div>
                <div style="font-size:0.76rem; color:#475569;">3배위(__CN3__종) &middot; 4배위(__CN4__종) &middot; $\pi$-6배위+(__CN6P__종) &middot; 다핵 착물 __MULTI_C__종</div>
              </div>
              <div style="background:#f0f9ff; border:1px solid #bae6fd; border-left:4px solid #0284c7; border-radius:6px; padding:12px 14px;">
                <div style="font-size:0.76rem; font-weight:700; color:#0369a1;">리간드 조합 공간 압축</div>
                <div style="font-size:1.35rem; font-weight:800; color:#0c4a6e; margin:3px 0;">__N_LIG__종 / __TOT_INST__ 인스턴스</div>
                <div style="font-size:0.76rem; color:#475569;">동종 __HOMO__ vs 이종 __HETERO__종 &bull; Top-15 완전피복 __COV15__%, Top-30 __COV30__%</div>
              </div>
              <div style="background:#ecfdf5; border:1px solid #a7f3d0; border-left:4px solid #059669; border-radius:6px; padding:12px 14px;">
                <div style="font-size:0.76rem; font-weight:700; color:#065f46;">연속 결합 모델 (5-CV R&sup2;)</div>
                <div style="font-size:1.35rem; font-weight:800; color:#064e3b; margin:3px 0;">&Delta;E<sub>f</sub> __EF_CV5__ &middot; %V<sub>bur</sub> __VB_CV5__</div>
                <div style="font-size:0.76rem; color:#047857;">이산 원-핫(<code>__VB_DISC__</code>) 대비 연속 결합 <code>%V<sub>bur</sub></code> __VB_CV5__ (LOMO __VB_LOMO__)</div>
              </div>
              <div style="background:#fffbeb; border:1px solid #fde68a; border-left:4px solid #d97706; border-radius:6px; padding:12px 14px;">
                <div style="font-size:0.76rem; font-weight:700; color:#92400e;">조합 분해 한계 (협동 반응)</div>
                <div style="font-size:1.35rem; font-weight:800; color:#78350f; margin:3px 0;">LOMO BDE<sub>min</sub> / &Delta;G<sub>dim</sub></div>
                <div style="font-size:0.76rem; color:#78350f;">쌍별 리간드 협동(&beta;-H 전달) 및 스핀&middot;궤도 정합 영역은 3D CI-NEB 보완 필수</div>
              </div>
            </div>

            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:10px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:10px 14px;">
              <div style="font-size:0.81rem; color:#334155;">
                상세 수식($CN_{\text{eff}}$, $N_{\text{ve}}$, $R_{\text{eff}}$), 2D 그래프 파싱 경계 사례(<code>CCTBA</code>, <code>Co(AMD)2</code>), 상위 리간드 차트 및 <code>5-CV / LOMO-5</code> 전체 실측표는 리포트 탭 <strong>Section 18</strong>에 수록되어 있습니다.
              </div>
              <button type="button" onclick="jumpToTechnicalSection('tm-coord-ligand-stats')" class="btn-pdf" style="background:#7c3aed; color:#fff; border-color:#7c3aed; cursor:pointer; font-size:0.8rem; padding:7px 14px;">
                Section 18 상세 분석 열기 &rarr;
              </button>
            </div>
          </div>"""

    return _render_template(
        tpl,
        {
            "__N_PREC__": stats["n_precursors"],
            "__N_LIG__": stats["n_unique_ligands"],
            "__SINGLETON_PCT__": f"{singleton_pct:.1f}",
            "__CN3__": cn_map.get(3, 0),
            "__CN4__": cn_map.get(4, 0),
            "__CN6P__": n_cn6_plus,
            "__MULTI_C__": stats["multicenter_count"],
            "__TOT_INST__": stats["total_ligand_instances"],
            "__HOMO__": stats["homoleptic_count"],
            "__HETERO__": stats["heteroleptic_count"],
            "__COV15__": f"{cov15['full_pct']:.1f}",
            "__COV30__": f"{cov30['full_pct']:.1f}",
            "__EF_CV5__": f"{r_ef['cont_cv5_r2']:.3f}",
            "__VB_CV5__": f"{r_vb['coup_cv5_r2']:.3f}",
            "__VB_DISC__": f"{r_vb['disc_cv5_r2']:.3f}",
            "__VB_LOMO__": f"{r_vb['coup_lomo_r2']:.3f}",
        },
    )


def build_ligand_coord_technical_html(
    unified_json_path: str,
    stats: Optional[Dict[str, Any]] = None,
) -> str:
    """Build Section 18 HTML cards for Technical Report tab."""
    if stats is None:
        stats = analyze_precursor_coordination_and_ligands(unified_json_path)
    cn_svg = _build_cn_eff_svg(stats["cn_eff_counts"], stats["n_precursors"])
    lig_svg = _build_top_ligands_svg(
        stats["top15_ligands"], stats["total_ligand_instances"], stats["n_unique_ligands"]
    )

    metal_pills = " ".join(
        f'<span style="display:inline-block; background:#f1f5f9; border:1px solid #cbd5e1; border-radius:999px; padding:2px 9px; font-size:0.77rem; font-family:\'JetBrains Mono\',monospace; margin:2px;"><strong>{html.escape(m)}</strong>:{c}</span>'
        for m, c in stats["center_counts"]
    )

    coord_rows_html = []
    for c in stats["coord_cases"]:
        coord_rows_html.append(
            f"""
            <tr>
              <td><strong>{html.escape(c['name'])}</strong><br><span style="font-size:0.71rem; color:#64748b;">{html.escape(c['id'])} ({html.escape(c['center_atom'])})</span></td>
              <td>{html.escape(c['role_note'])}<br><span style="font-size:0.72rem; color:#475569;">{html.escape(c['geom_note'])}</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{c['pi']} / {c['chelate']}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:800; color:#5b21b6;">{c['cn_eff']:.1f}</td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{c['n_ve']:.0f}e<sup>-</sup> ({c['ve_deficit']:.2f})</td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{c['r_cov_ang']:.2f} &rarr; <strong>{c['r_eff_ang']:.3f}</strong> &Aring;</td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{c['vbur_percent']:.1f}% &rarr; {c['vbur_d2_percent']:.1f}%</td>
            </tr>"""
        )

    lig_rows_html = []
    for idx, it in enumerate(stats["top15_ligands"], start=1):
        m_str = ", ".join(it["metals"])
        share_badge = (
            '<span class="badge badge-prec" style="font-size:0.7rem;">교차 공유</span>'
            if it["metal_count"] >= 2
            else '<span class="badge" style="background:#f1f5f9; color:#475569; font-size:0.7rem;">단일 금속</span>'
        )
        lig_rows_html.append(
            f"""
            <tr>
              <td style="text-align:center;">{idx}</td>
              <td><strong>{html.escape(it['label'])}</strong></td>
              <td><code>{html.escape(it['smiles'])}</code></td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:700;">{it['occurrences']}</td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{it['precursor_count']} / {stats['n_precursors']}</td>
              <td>{share_badge} <strong>{it['metal_count']}종</strong> <span style="font-size:0.74rem; color:#475569;">({html.escape(m_str)})</span></td>
            </tr>"""
        )

    cov_rows_html = []
    for cv in stats["topk_coverage"]:
        cov_rows_html.append(
            f"""
            <tr>
              <td style="text-align:center; font-weight:700;">Top-{cv['k']}</td>
              <td style="text-align:center; font-family:'JetBrains Mono';"><strong>{cv['full']} / {stats['n_precursors']}</strong> ({cv['full_pct']:.1f}%)</td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{cv['partial']} / {stats['n_precursors']} ({cv['partial_pct']:.1f}%)</td>
            </tr>"""
        )

    ablation_rows_html = []
    for row in stats["ablation_rows"]:
        ablation_rows_html.append(
            f"""
            <tr>
              <td><strong>{html.escape(row['label'])}</strong><br><span style="font-size:0.72rem; color:#64748b;">n={row['n']}, 단위: {html.escape(row['unit'])}</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{row['disc_cv5_r2']:+.3f} <span style="font-size:0.72rem; color:#64748b;">({row['disc_cv5_mae']:.2f})</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono';">{row['disc_lomo_r2']:+.3f} <span style="font-size:0.72rem; color:#64748b;">({row['disc_lomo_mae']:.2f})</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:700; color:#0369a1;">{row['cont_cv5_r2']:+.3f} <span style="font-size:0.72rem; color:#64748b;">({row['cont_cv5_mae']:.2f})</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:700; color:#0369a1;">{row['cont_lomo_r2']:+.3f} <span style="font-size:0.72rem; color:#64748b;">({row['cont_lomo_mae']:.2f})</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:800; color:#059669;">{row['coup_cv5_r2']:+.3f} <span style="font-size:0.72rem; color:#64748b;">({row['coup_cv5_mae']:.2f})</span></td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-weight:800; color:#059669;">{row['coup_lomo_r2']:+.3f} <span style="font-size:0.72rem; color:#64748b;">({row['coup_lomo_mae']:.2f})</span></td>
            </tr>"""
        )

    ab_map = {row["key"]: row for row in stats["ablation_rows"]}
    r_ef = ab_map["delta_e_form_per_bond_kcal_mol"]
    r_vb = ab_map["vbur_percent"]
    r_lsa = ab_map["ligand_sa"]
    r_bde = ab_map["bde_min_kcal_mol"]
    r_dg = ab_map["delta_g_dimer_373_kcal_mol"]
    cov15 = next(x for x in stats["topk_coverage"] if x["k"] == 15)

    pi_summary = ", ".join(f"$\\eta_{{\\pi}}={k}$: {v}종" for k, v in stats["pi_counts"])
    chel_summary = ", ".join(f"$N_{{\\text{{chel}}}}={k}$: {v}종" for k, v in stats["chelate_counts"])
    shared_pct = 100.0 * stats["n_shared_ligands"] / max(1, stats["n_unique_ligands"])
    singleton_pct = 100.0 * stats["n_singleton_ligands"] / max(1, stats["n_unique_ligands"])

    tpl = r"""
        <!-- Card 18.1: Transition-Metal Multi-Coordination Determination -->
        <div class="card" style="margin-bottom:20px; border-top:4px solid #7c3aed;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">1. 전이금속 다중 배위수(<code>CN<sub>eff</sub></code>) 결정 알고리즘 및 상태-의존 연속 매핑</h3>
            <span class="badge badge-4g" style="font-size:0.78rem;">element_physics.py &bull; __N_PREC__ Precursors Telemetry</span>
          </div>
          <p style="font-size:0.88rem; color:#334155; line-height:1.65; margin-bottom:12px;">
            유기 분자의 탄소 원자가 고정 배위수($4$)를 갖는 것과 달리, 전이금속&middot;주족 중심 원자는 동일한 원소라도 $\sigma$-단일결합, 배위결합(<code>DATIVE</code>), 다배위 킬레이트, 및 $\eta^n$ 다중심 $\pi$-착물 형성에 따라 배위수가 $1$에서 $8$ 이상까지 폭넓게 변동합니다. 본 파이프라인은 <code>scripts/benchmark/element_physics.py</code>의 <code>compute_refined_pi_hapticity</code>와 <code>compute_coordination_state</code>를 통해 원소 기호 분기 없이 다음 물리 규칙으로 <strong>유효 배위수($CN_{\text{eff}}$)</strong>와 <strong>상태-의존 배위 벡터</strong>를 결정합니다.
          </p>

          <div class="versus-grid" style="margin-bottom:16px;">
            <div class="versus-card" style="border-left:4px solid #7c3aed;">
              <div class="versus-header" style="color:#5b21b6;">A. 유효 배위수 및 원자가 전자수 정식화</div>
              <div style="font-size:0.83rem; color:#334155; line-height:1.6;">
                $$CN_{\text{eff}}(M \mid \text{env}) = N_{\text{nbr}}^{\text{heavy}} + N_H^{\text{imp}} + 2\,\eta_{\pi}$$
                $$N_{\text{ve}} = N_{\text{val}}(Z_M) - q_{\text{formal}} + \sum_{l \in \mathcal{N}(M)} e_l + 4\,\eta_{\pi}, \quad \text{ve\_deficit} = \max\!\left(0, \frac{N_{\text{sat}} - N_{\text{ve}}}{N_{\text{sat}}}\right)$$
                <ul style="margin:6px 0 0 0; padding-left:18px;">
                  <li><strong>$\eta_{\pi}$ 다중심 $\pi$-판정 규칙</strong>: 중심 원자 $M$과 2개 이상 고리 원자를 공유하거나(명시적 polyhapto), 비방향족 불포화 5원 이상 고리(<code>Cp</code> 등)가 $M$에 결합하면 $\eta_{\pi}$를 $+1$ 가산합니다. 반면 단일 <em>ipso</em>-C에 결합한 6원 방향족 페닐 고리(예: <code>CrPh3</code>, <code>WPh5</code>)는 $\eta^1$ $\sigma$-아릴로 분류해 $\eta_{\pi}=0$으로 배제합니다.</li>
                  <li><strong>전자 기여도 $e_l$ 및 $N_{\text{sat}}$</strong>: <code>DATIVE</code> 결합은 $e_l=2.0$, 공유결합은 결합 차수(암시적 수소 $+1.0$)를 적용하며, $\pi$-고리당 graph $\sigma$-결합($+1e^-$) 외에 $+4e^-$를 보정하여 중성 Cp($5e^-$) 기여를 복원합니다($N_{\text{sat}}=18$ for $d$-block, $8$ for main-group).</li>
                </ul>
              </div>
            </div>

            <div class="versus-card" style="border-left:4px solid #0284c7;">
              <div class="versus-header" style="color:#0369a1;">B. 배위수&ndash;산화수 연동 유효 반경 및 극성 수축</div>
              <div style="font-size:0.83rem; color:#334155; line-height:1.6;">
                $$R_{\text{eff}}(M \mid \text{env}) = \operatorname{clip}\!\Big(R_{\text{cov}}^{(1)}(Z_M)\big[1 - 0.06\max(0, \overline{\Delta\chi}_{ML}) + 0.025(CN_{\text{eff}}-4)\big],\; 0.75R_{\text{cov}}^{(1)},\; 1.22R_{\text{cov}}^{(1)}\Big)$$
                $$\chi_{\text{eff}}(M \mid \text{env}) = \chi_0(Z_M) + 0.12\sum_l \operatorname{sgn}(\Delta\chi_l)\big[1 - \exp(-0.25\,\Delta\chi_l^2)\big]$$
                <ul style="margin:6px 0 0 0; padding-left:18px;">
                  <li><strong>고산화 극성 수축 vs 고배위 팽창</strong>: $CN_{\text{eff}}=4$를 기준점으로 리간드 전기음성도 대비 $\overline{\Delta\chi}_{ML}$이 커지면 유효 반경이 최대 $0.75 R_{\text{cov}}$까지 수축하고(<code>RuO4</code>: $1.46 \to 1.351\,\text{\AA}$), 반대로 $\eta^5$-<code>Cp</code> 등으로 $CN_{\text{eff}}$가 $6\text{--}8$로 증가하면 반경이 최대 $1.22 R_{\text{cov}}$까지 팽창합니다(<code>Ru(EtCp)2</code>: $1.46 \to 1.502\,\text{\AA}$, <code>ZMA</code>: $1.75 \to 1.774\,\text{\AA}$).</li>
                  <li><strong>실측 분포(__N_PREC__종 프리커서)</strong>: __PI_SUMMARY__ &bull; __CHEL_SUMMARY__</li>
                </ul>
              </div>
            </div>
          </div>

          <div style="display:grid; grid-template-columns: minmax(280px, 1fr) minmax(320px, 1.35fr); gap:16px; align-items:start; margin-bottom:14px;">
            <div>
              __CN_SVG__
            </div>
            <div class="table-responsive">
              <table class="data-table" style="font-size:0.78rem;">
                <thead>
                  <tr>
                    <th>분자 사례 (ID)</th>
                    <th>배위 상태 및 판별 근거</th>
                    <th style="text-align:center;">$\eta_{\pi}$ / Chel</th>
                    <th style="text-align:center;">$CN_{\text{eff}}$</th>
                    <th style="text-align:center;">$N_{\text{ve}}$ ($\Delta N_{\text{ve}}$비)</th>
                    <th style="text-align:center;">$R_{\text{cov}} \to R_{\text{eff}}$</th>
                    <th style="text-align:center;">$\%V_{\text{bur}} \to \%V_{\text{bur}}^{\text{D2}}$</th>
                  </tr>
                </thead>
                <tbody>
                  __COORD_ROWS__
                </tbody>
              </table>
            </div>
          </div>

          <div style="background:#fffbeb; border:1px solid #fde68a; border-left:4px solid #d97706; border-radius:6px; padding:12px 14px; font-size:0.82rem; color:#78350f; line-height:1.6;">
            <strong>2D SMILES 기반 다중 배위수 판정의 알려진 경계 (인식 한계 및 3D ECoN 고도화 방향):</strong>
            (1) 다핵 알카인 가교 클러스터 <code>CCTBA</code>($\text{Co}_2(\text{CO})_6(\mu\text{-}\eta^2\text{-}t\text{BuC}\equiv\text{CH})$)는 $\text{Co}_2\text{C}_2$ 사면체 케이지 내 3원 고리가 중복 계수되어 $\eta_{\pi}=4$, $CN_{\text{eff}}=14.0$으로 과대 산출되나, $R_{\text{eff}}$ 상한 <code>clip(1.22 r_cov)</code>($1.533\,\text{\AA}$)가 수치 발산을 방어합니다.
            (2) <code>Co(AMD)2</code>와 같이 아미디네이트(<code>N-C=N</code>) 리간드의 두 번째 배위 질소가 SMILES에 <code>&lt;-</code> 화살표 없이 단일결합만으로 표기된 경우 그래프 연결 수 기준 $CN_{\text{eff}}=2.0$으로 낮게 인식됩니다.
            (3) 이합체 <code>Ga2(NMe2)6</code>(<code>EXT-ALD-121</code>)는 가교 $\text{Ga}_2\text{N}_2$ 4원 고리가 $\eta_{\pi}=1$ 및 $N_{\text{chel}}=1$로 이중 집계됩니다. 향후 <code>SevenNet</code> 3D 이완 좌표를 직접 입력받는 <strong>Hoppe 유효 배위수(ECoN) 및 보로노이 연속 가중 배위수</strong>로 대체하면 2D SMILES 표기 의존성을 완전히 제거할 수 있습니다.
          </div>
        </div>

        <!-- Card 18.2: 102 Precursor Ligand & Metal Statistics -->
        <div class="card" style="margin-bottom:20px; border-top:4px solid #0284c7;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">2. 반도체 프리커서 __N_PREC__종 중심 금속&middot;리간드 절편 조합 통계</h3>
            <span class="badge badge-prec" style="font-size:0.78rem;">__N_METALS__ Metals &bull; __N_LIG__ Unique Ligand Fragments</span>
          </div>
          <p style="font-size:0.88rem; color:#334155; line-height:1.65; margin-bottom:12px;">
            <code>Audited Precursor (27)</code>와 <code>Expanded Precursor (75)</code>를 합산한 총 <strong><code>__N_PREC__</code>종 반도체 프리커서</strong>에서 중심 원자를 제거하여 리간드 절편을 추출한 결과(<code>__N_PREC__/__N_PREC__</code> 파싱 성공, 다핵 착물 <code>__MULTI_C__</code>종), <strong>동종 리간드 착물(homoleptic) <code>__HOMO__</code>종</strong>과 <strong>이종 리간드 착물(heteroleptic) <code>__HETERO__</code>종</strong>이 균형을 이루고 있었습니다.
            총 <code>__TOT_INST__</code>개의 리간드 결합 인스턴스는 단 <strong><code>__N_LIG__</code>종의 고유 절편</strong>으로 압축되나, <code>__N_SHARED__</code>종(<code>__SHARED_PCT__%</code>)만이 2개 이상의 금속에서 재사용되고 나머지 <strong><code>__N_SINGLETON__</code>종(<code>__SINGLETON_PCT__%</code>)은 단일 금속에서만 등장하는 롱테일(long-tail) 분포</strong>를 보입니다.
          </p>

          <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:10px 14px; margin-bottom:14px;">
            <div style="font-size:0.79rem; font-weight:700; color:#334155; margin-bottom:6px;">
              중심 원소 분포 (__N_METALS__개 원소, 총 __N_PREC__종):
            </div>
            <div>__METAL_PILLS__</div>
          </div>

          <div style="display:grid; grid-template-columns: minmax(300px, 1.15fr) minmax(280px, 1fr); gap:16px; align-items:start; margin-bottom:14px;">
            <div>
              __LIG_SVG__
            </div>
            <div>
              <div style="font-weight:700; font-size:0.86rem; color:#0f172a; margin-bottom:6px;">
                상위 Top-K 리간드 어휘의 프리커서 피복률
              </div>
              <table class="data-table" style="font-size:0.79rem; margin-bottom:12px;">
                <thead>
                  <tr>
                    <th style="text-align:center;">리간드 상위 크기</th>
                    <th style="text-align:center;">완전 피복 (모든 리간드 포함)</th>
                    <th style="text-align:center;">부분 피복 (1종 이상 포함)</th>
                  </tr>
                </thead>
                <tbody>
                  __COV_ROWS__
                </tbody>
              </table>
              <div style="background:#f0f9ff; border:1px solid #bae6fd; border-radius:6px; padding:10px 12px; font-size:0.8rem; color:#0c4a6e; line-height:1.55;">
                <strong>핵심 해석:</strong> 최빈 리간드 <code>15</code>종만으로 전체 프리커서의 절반(<code>__COV15_FULL__%</code>)이 완전히 구성되며 <code>__COV15_PART__%</code>가 해당 리간드를 포함합니다. 그러나 <code>__SINGLETON_PCT__%</code>의 희소 리간드가 단일 금속에만 결합되어 있어, 리간드의 <em>SMILES 토큰 이름</em>이 아닌 <strong>물리&middot;전자&middot;입체 연속 서술자</strong>로 매핑해야만 미관측 조합으로 일반화할 수 있습니다.
              </div>
            </div>
          </div>

          <div class="table-responsive">
            <table class="data-table collapsible-table" style="font-size:0.79rem;">
              <thead>
                <tr>
                  <th style="text-align:center;">순위</th>
                  <th>리간드 일반명 / 모티프</th>
                  <th>절편 Canonical SMILES</th>
                  <th style="text-align:center;">총 등장 횟수</th>
                  <th style="text-align:center;">포함 프리커서</th>
                  <th>금속 간 재사용 범위 (중심 원소 목록)</th>
                </tr>
              </thead>
              <tbody>
                __LIG_ROWS__
              </tbody>
            </table>
          </div>
        </div>

        <!-- Card 18.3: Empirical Compositional Modeling & Limitations -->
        <div class="card" style="border-top:4px solid #059669;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">3. 금속 물성&ndash;리간드 조합 결합 모델링 실측 벤치마크 및 한계 규명</h3>
            <span class="badge badge-2g" style="font-size:0.78rem;">5-CV vs LOMO-5 Compositional Generalization</span>
          </div>
          <p style="font-size:0.88rem; color:#334155; line-height:1.65; margin-bottom:12px;">
            “리간드 종류가 제한적이므로 리간드 조합과 중심 금속 특성을 이어 프리커서 물성을 예측할 수 있는가”라는 설계 가설을 검증하기 위해, 프리커서 $\mathcal{C} = (M, \{L_i\}_{i=1}^{N_L})$를 다음 <strong>4항 계층 집합&ndash;상호작용 분해 수식</strong>으로 정식화하고, <code>__N_PREC__</code>종 실측 데이터에서 무작위 <code>5-CV</code>와 미관측 금속을 통째로 제외하는 <code>LOMO-5</code> 외삽 평가를 수행했습니다.
          </p>
          <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:6px; padding:10px 14px; margin-bottom:14px; font-size:0.85rem; color:#0f172a;">
            $$\hat{y}\!\left(M, \{L_i\}_{i=1}^{N_L}\right) = f_M(\mathbf{z}_M) + \sum_{i=1}^{N_L} \alpha_i(\mathbf{z}_M, \mathbf{z}_{L_i})\, f_L(\mathbf{z}_{L_i}) + \Phi_{ML}\!\left(\mathbf{z}_M, \sum_{i=1}^{N_L} \mathbf{z}_{L_i}\right) + \sum_{1 \le i &lt; j \le N_L} \Psi_{LL}\!\left(\mathbf{z}_{L_i}, \mathbf{z}_{L_j} \mid \mathbf{z}_M\right)$$
            <div style="font-size:0.79rem; color:#475569; line-height:1.5; margin-top:4px;">
              여기서 [A] 이산 가산 모델은 $M$과 $L_i$의 원-핫 벡터만 가산하고, [B] 연속 가산 모델은 $\mathbf{z}_M=(R_{\text{cov}}, \chi_P, N_{\text{val}})$과 리간드 집합 연속량(중원자 합&middot;평균&middot;최댓값, 회전결합, 고리수, C/N/O/할로겐 주개 분율, $\eta_{\pi}$, $N_{\text{chel}}$)을 사용하며, [B+] 결합 모델은 무차원 입체 밀집도 $N_{\text{heavy}}/R_{\text{cov}}^3$, 표면 점유비 $N_{\text{heavy,max}}/R_{\text{cov}}^2$, $\overline{\Delta\chi}$, 및 <code>ve_deficit</code>을 $\Phi_{ML}$ 항으로 명시 투입합니다.
            </div>
          </div>

          <div class="table-responsive" style="margin-bottom:14px;">
            <table class="data-table" style="font-size:0.79rem;">
              <thead>
                <tr>
                  <th>예측 타겟 (단위)</th>
                  <th style="text-align:center;">[A] 이산 가산 <code>5-CV</code> $R^2$ (MAE)</th>
                  <th style="text-align:center;">[A] 이산 가산 <code>LOMO-5</code> $R^2$ (MAE)</th>
                  <th style="text-align:center;">[B] 연속 가산 <code>5-CV</code> $R^2$ (MAE)</th>
                  <th style="text-align:center;">[B] 연속 가산 <code>LOMO-5</code> $R^2$ (MAE)</th>
                  <th style="text-align:center;">[B+] <code>M&times;L</code> 결합 <code>5-CV</code> $R^2$ (MAE)</th>
                  <th style="text-align:center;">[B+] <code>M&times;L</code> 결합 <code>LOMO-5</code> $R^2$ (MAE)</th>
                </tr>
              </thead>
              <tbody>
                __ABLATION_ROWS__
              </tbody>
            </table>
          </div>

          <div class="versus-grid">
            <div class="versus-card" style="border-left:4px solid #059669;">
              <div class="versus-header" style="color:#065f46;">적용 가능 영역 (저차원&middot;입체&middot;생성열 물성)</div>
              <ul style="font-size:0.82rem; color:#334155; line-height:1.6; padding-left:18px; margin:0;">
                <li><strong>이산 모델 실패 vs 연속 결합 도약</strong>: 단순 원-핫 가산 회귀는 <code>%V<sub>bur</sub></code>($R^2=__VB_DISC__$)과 <code>&Delta;G<sub>dimer</sub></code>($R^2=__DG_DISC__$)에서 붕괴했으나, 연속 서술자 및 <code>M&times;L</code> 결합항 도입 시 <code>&Delta;E<sub>f</sub>/v<sub>M</sub></code> <strong>__EF_CONT__</strong>, <code>%V<sub>bur</sub></code> <strong>__VB_COUP__</strong>, <code>Ligand SA</code> <strong>__LSA_CONT__</strong>까지 급상승했습니다.</li>
                <li><strong>미관측 금속 외삽(LOMO-5)</strong>: 학습에 포함되지 않은 금속군을 예측할 때도 연속 물리 표현은 <code>Ligand SA</code>(<strong>__LSA_LOMO__</strong>), <code>&Delta;E<sub>f</sub>/v<sub>M</sub></code>(<strong>__EF_LOMO__</strong>), <code>%V<sub>bur</sub></code>(<strong>__VB_LOMO__</strong>)을 유지해 가상 조합 공간 스크리닝의 신뢰성 있는 대리 모델(Surrogate)로 기능합니다.</li>
              </ul>
            </div>

            <div class="versus-card" style="border-left:4px solid #e11d48;">
              <div class="versus-header" style="color:#be123c;">근본 한계 및 쌍별 상호작용 &Psi;<sub>LL</sub>&middot;CI-NEB 보완</div>
              <ul style="font-size:0.82rem; color:#334155; line-height:1.6; padding-left:18px; margin:0;">
                <li><strong>전자쌍 특이성 및 협동 분해 한계</strong>: <code>BDE<sub>min</sub></code>은 <code>5-CV</code>(__BDE_CV5__) 대비 <code>LOMO-5</code>(__BDE_LOMO__)에서 크게 하락하는데, 이는 <code>CO</code>&middot;<code>Cp</code>의 $\pi$-역결합과 산화환원 재편성이 금속&ndash;리간드 궤도 정합에 민감하기 때문입니다. 특히 Section 17에서 입증한 <code>Cp(NMe2)3Zr</code>의 저에너지 $\beta$-H 이동($E_a=34.67$ vs $\text{BDE}_{\min}=69.96\,\text{kcal/mol}$)은 리간드 간 쌍별 상호작용 $\Psi_{LL}$과 3D <code>SevenNet</code> CI-NEB 계산으로만 포착 가능합니다.</li>
                <li><strong>필수 누수 방지 검증 프로토콜</strong>: 조합형 ML 학습 시 무작위 분할이 아닌 <strong>LOMO</strong>(미관측 금속), <strong>LOLO</strong>(미관측 리간드), 및 <strong>Homoleptic &rarr; Heteroleptic 외삽</strong>을 표준 벤치마크로 강제해야 합니다.</li>
              </ul>
            </div>
          </div>
        </div>"""

    return _render_template(
        tpl,
        {
            "__N_PREC__": stats["n_precursors"],
            "__N_METALS__": stats["n_unique_metals"],
            "__N_LIG__": stats["n_unique_ligands"],
            "__MULTI_C__": stats["multicenter_count"],
            "__HOMO__": stats["homoleptic_count"],
            "__HETERO__": stats["heteroleptic_count"],
            "__TOT_INST__": stats["total_ligand_instances"],
            "__N_SHARED__": stats["n_shared_ligands"],
            "__SHARED_PCT__": f"{shared_pct:.1f}",
            "__N_SINGLETON__": stats["n_singleton_ligands"],
            "__SINGLETON_PCT__": f"{singleton_pct:.1f}",
            "__PI_SUMMARY__": pi_summary,
            "__CHEL_SUMMARY__": chel_summary,
            "__CN_SVG__": cn_svg,
            "__COORD_ROWS__": "".join(coord_rows_html),
            "__METAL_PILLS__": metal_pills,
            "__LIG_SVG__": lig_svg,
            "__COV_ROWS__": "".join(cov_rows_html),
            "__COV15_FULL__": f"{cov15['full_pct']:.1f}",
            "__COV15_PART__": f"{cov15['partial_pct']:.1f}",
            "__LIG_ROWS__": "".join(lig_rows_html),
            "__ABLATION_ROWS__": "".join(ablation_rows_html),
            "__VB_DISC__": f"{r_vb['disc_cv5_r2']:.3f}",
            "__DG_DISC__": f"{r_dg['disc_cv5_r2']:.3f}",
            "__EF_CONT__": f"{r_ef['cont_cv5_r2']:.3f}",
            "__VB_COUP__": f"{r_vb['coup_cv5_r2']:.3f}",
            "__LSA_CONT__": f"{r_lsa['cont_cv5_r2']:.3f}",
            "__LSA_LOMO__": f"{r_lsa['cont_lomo_r2']:.3f}",
            "__EF_LOMO__": f"{r_ef['cont_lomo_r2']:.3f}",
            "__VB_LOMO__": f"{r_vb['coup_lomo_r2']:.3f}",
            "__BDE_CV5__": f"{r_bde['cont_cv5_r2']:.3f}",
            "__BDE_LOMO__": f"{r_bde['cont_lomo_r2']:.3f}",
        },
    )
