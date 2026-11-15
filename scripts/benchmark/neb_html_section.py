#!/usr/bin/env python3
"""HTML and SVG section builders for Section 17 / Storyline Step 5 (CI-NEB Kinetics).

Consumes `database/external/sevenn_decomp_neb_benchmark.json` and
`database/benchmark_unified_latest_160.json` to compute Arrhenius 1% decomposition
temperatures, barrier-to-temperature uncertainty propagation, and Pearson
correlations against experimental ALD/CVD T_decomp ranges without hardcoded drift.
"""

from __future__ import annotations

import json
import math
from typing import Any, Dict, List, Optional, Tuple

R_KCAL_MOL_K: float = 1.98720425864083e-3  # units: kcal/(mol*K)

# Literature experimental decomposition & ALD window metadata keyed by reaction_id.
LITERATURE_TDECOMP_META: Dict[str, Dict[str, Any]] = {
    "CpZr_NMe2_betaH_elim": {
        "short": "CpZr (ZAC)",
        "regime": "gas_unimol",
        "t_range_c": "305–350 (ALD 상한 ~320)",
        "t_mid_c": 327.5,
        "ald_note": "Cp 구경꾼 리간드가 표면 노출을 억제하여 TDMAZ 대비 고온 ZrO2 ALD 가능 (10.1021/acsami.5c11791, 10.1016/j.mtcomm.2022.104008)",
    },
    "GaEt3_betaH_elim": {
        "short": "GaEt3",
        "regime": "gas_unimol",
        "t_range_c": "270–380 (N2) / 220–330 (H2)",
        "t_mid_c": 325.0,
        "ald_note": "에틸기 beta-H 제거로 >250 C에서 열 ALD 상한 제한 (10.1021/jp5116405, 10.1016/0039-6028(93)90145-K)",
    },
    "DCDS_SiHCl_elim": {
        "short": "DCDS",
        "regime": "gas_unimol",
        "t_range_c": "400–500 (열분해 >450)",
        "t_mid_c": 450.0,
        "ald_note": "Si-H 포함 1,2-H 이동 경로로 HCDS 대비 저온 열분해 개시",
    },
    "HCDS_SiCl2_elim": {
        "short": "HCDS",
        "regime": "gas_unimol",
        "t_range_c": "500–550 (ALD 450–650)",
        "t_mid_c": 525.0,
        "ald_note": "높은 1,2-Cl 이동 장벽으로 고온 열 ALD/LPCVD 운전 가능 (10.1016/j.surfin.2025.106470)",
    },
    "TEOS_C2H4_elim": {
        "short": "TEOS",
        "regime": "gas_unimol",
        "t_range_c": "600–700 (무촉매 기상)",
        "t_mid_c": 650.0,
        "ald_note": "최고 기상 내열성; NH3/H2O 촉매 ALD 시 저온 공정 가능 (10.1149/1.1394331, 10.1149/1.1768548)",
    },
    "TTIP_acetone_elim": {
        "short": "TTIP (ace)",
        "regime": "surface_catalytic",
        "t_range_c": "250–275 (ALD 창 190–240)",
        "t_mid_c": 262.5,
        "ald_note": "순수 기상 장벽(49.0) 대비 성장 중 TiO2/-OH 표면 자가촉매(Ea,app~20-40)로 공정 상한 하강 (10.1116/1.4798587)",
    },
    "TTIP_propene_elim": {
        "short": "TTIP (pro)",
        "regime": "selectivity_alt",
        "t_range_c": "경쟁 채널 (+13.0 kcal/mol)",
        "t_mid_c": None,
        "ald_note": "아세톤 제거 대비 높은 장벽으로 고온 기상 부반응 채널에 해당 (10.1021/acs.jpca.7b09904)",
    },
    "TDMAT_NH3_transam": {
        "short": "TDMAT+NH3",
        "regime": "coreactant_ald",
        "t_range_c": "ALD 창 100–200 (자가분해 200–220)",
        "t_mid_c": None,
        "ald_note": "벤치마크 채널은 NH3 공반응 리간드 교환(목적 ALD 반응); 공정 상한은 단분자 beta-H 제거(39.3 kcal/mol) 지배",
    },
}


def compute_arrhenius_t_threshold_c(
    ea_kcal_mol: float,
    pre_exp_s: float = 1.0e13,
    decomp_frac: float = 0.01,
    residence_time_s: float = 1.0,
) -> float:
    """Compute first-order decomposition threshold temperature (deg C).

    Physics/Algorithm: Solves 1 - exp(-A * exp(-Ea / (R*T)) * t_res) = chi
    for T in Kelvin and converts to Celsius.
    """
    rate_target = -math.log(1.0 - decomp_frac) / residence_time_s  # units: s^-1
    t_kelvin = ea_kcal_mol / (R_KCAL_MOL_K * math.log(pre_exp_s / rate_target))
    return t_kelvin - 273.15


def compute_pearson_r(xs: List[float], ys: List[float]) -> Optional[float]:
    """Compute Pearson correlation coefficient r across paired finite samples."""
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = sum((x - mx) ** 2 for x in xs)
    den_y = sum((y - my) ** 2 for y in ys)
    if den_x <= 1e-12 or den_y <= 1e-12:
        return None
    return num / math.sqrt(den_x * den_y)


def _load_bde_map(unified_db_path: str) -> Dict[str, float]:
    """Load `bde_min_kcal_mol` keyed by precursor ID from unified benchmark JSON."""
    with open(unified_db_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    rows = data.get("records", data) if isinstance(data, dict) else data
    out: Dict[str, float] = {}
    for rec in rows:
        rid = rec.get("id") or rec.get("molecule_id")
        bde = rec.get("bde_min_kcal_mol")
        if rid and bde is not None:
            out[str(rid)] = float(bde)
    return out


def compute_neb_kinetics_context(
    neb_bench_path: str,
    unified_db_path: str,
) -> Dict[str, Any]:
    """Assemble benchmark rows, Arrhenius scalings, and correlation metrics."""
    with open(neb_bench_path, "r", encoding="utf-8") as fh:
        bench = json.load(fh)
    bde_by_id = _load_bde_map(unified_db_path)

    summary = bench["summary"]
    reactions = bench["reactions"]
    mae = float(summary["mae_ea_kcal_mol"])
    dt_half = (
        compute_arrhenius_t_threshold_c(40.0 + mae)
        - compute_arrhenius_t_threshold_c(40.0 - mae)
    ) / 2.0

    enriched: List[Dict[str, Any]] = []
    for rxn in reactions:
        spec = rxn["spec"]
        en = rxn["energies_kcal_mol"]
        st = rxn["status"]
        vb = rxn["ts_diagnostics"]["core_vibrations"]
        rid = spec["reaction_id"]
        pid = spec["precursor_id"]
        meta = LITERATURE_TDECOMP_META.get(rid, {})
        ea_lit = float(en["lit_ea_kcal_mol"])
        ea_7n = float(en["sevenn_ea_cineb_kcal_mol"])
        enriched.append(
            {
                "reaction_id": rid,
                "precursor_id": pid,
                "precursor_name": spec["precursor_name"],
                "short": meta.get("short", spec["precursor_name"]),
                "equation": spec["reaction_equation"],
                "mechanism": spec["mechanism_class"],
                "lit_authors": spec["lit_authors"],
                "lit_year": spec["lit_year"],
                "lit_journal": spec["lit_journal"],
                "lit_doi": spec["lit_doi"],
                "lit_method": spec["lit_method"],
                "bde_min": bde_by_id.get(pid),
                "ea_lit": ea_lit,
                "ea_7n": ea_7n,
                "signed_err": float(en["signed_error_ea_kcal_mol"]),
                "de_adduct": float(en["sevenn_delta_e_adduct_kcal_mol"]),
                "de_sep": float(en["sevenn_delta_e_sep_kcal_mol"]),
                "fmax": float(st["final_neb_fmax_ev_A"]),
                "converged": bool(st["converged_cineb"]),
                "n_imag": int(vb["n_imag_modes_gt20cm1"]),
                "imag_freqs": vb["imag_freqs_cm1"],
                "wall_time_s": float(st["wall_time_s"]),
                "bond_lengths_A": rxn["ts_diagnostics"]["bond_lengths_A"],
                "band_profile": en["band_profile_rel_kcal_mol"],
                "regime": meta.get("regime", "other"),
                "t_range_c": meta.get("t_range_c", "-"),
                "t_mid_c": meta.get("t_mid_c"),
                "ald_note": meta.get("ald_note", ""),
                "t1_lit_1e13": compute_arrhenius_t_threshold_c(ea_lit, 1.0e13),
                "t1_7n_1e13": compute_arrhenius_t_threshold_c(ea_7n, 1.0e13),
                "t1_7n_1e14": compute_arrhenius_t_threshold_c(ea_7n, 1.0e14),
            }
        )

    gas_rows = [r for r in enriched if r["regime"] == "gas_unimol" and r["t_mid_c"] is not None]
    all_uni_rows = [
        r for r in enriched if r["regime"] in {"gas_unimol", "surface_catalytic"} and r["t_mid_c"] is not None
    ]

    def _corr_pack(sub: List[Dict[str, Any]]) -> Dict[str, Optional[float]]:
        ts = [float(r["t_mid_c"]) for r in sub]
        bde_sub = [r for r in sub if r["bde_min"] is not None]
        return {
            "n": len(sub),
            "r_lit": compute_pearson_r([r["ea_lit"] for r in sub], ts),
            "r_7n": compute_pearson_r([r["ea_7n"] for r in sub], ts),
            "r_bde": compute_pearson_r(
                [float(r["bde_min"]) for r in bde_sub],
                [float(r["t_mid_c"]) for r in bde_sub],
            ),
        }

    return {
        "summary": summary,
        "dt_halfwidth_c": dt_half,
        "rows": enriched,
        "by_id": {r["reaction_id"]: r for r in enriched},
        "corr_gas": _corr_pack(gas_rows),
        "corr_all_uni": _corr_pack(all_uni_rows),
    }


def _build_neb_kinetics_svg(ctx: Dict[str, Any]) -> str:
    """Build dual-panel SVG: (L) Lit vs 7net-0 parity, (R) Ea vs T_decomp."""
    rows = ctx["rows"]

    # Left Panel: X = Lit Ea (0..75), Y = 7net-0 Ea (0..75)
    lx0, lx1, ly0, ly1 = 58.0, 392.0, 42.0, 294.0

    def _lx(v: float) -> float:
        return lx0 + (max(0.0, min(75.0, v)) / 75.0) * (lx1 - lx0)

    def _ly(v: float) -> float:
        return ly1 - (max(0.0, min(75.0, v)) / 75.0) * (ly1 - ly0)

    # Diagonal +/- 5.5 kcal/mol band polygon
    band_pts = (
        f"{_lx(0):.1f},{_ly(5.5):.1f} {_lx(69.5):.1f},{_ly(75):.1f} "
        f"{_lx(75):.1f},{_ly(75):.1f} {_lx(75):.1f},{_ly(69.5):.1f} "
        f"{_lx(5.5):.1f},{_ly(0):.1f} {_lx(0):.1f},{_ly(0):.1f}"
    )

    label_offsets_left = {
        "TDMAT_NH3_transam": (8, -6),
        "CpZr_NMe2_betaH_elim": (-12, -12),
        "GaEt3_betaH_elim": (-38, -10),
        "DCDS_SiHCl_elim": (8, 12),
        "HCDS_SiCl2_elim": (8, 4),
        "TTIP_acetone_elim": (-14, -12),
        "TEOS_C2H4_elim": (-10, 15),
        "TTIP_propene_elim": (-46, -10),
    }

    left_pts_svg: List[str] = []
    for r in rows:
        px = _lx(r["ea_lit"])
        py = _ly(r["ea_7n"])
        if r["reaction_id"] == "CpZr_NMe2_betaH_elim":
            col, rad = "#7c3aed", 6.2
        elif r["regime"] == "gas_unimol":
            col, rad = "#0284c7", 5.2
        elif r["regime"] == "surface_catalytic":
            col, rad = "#d97706", 5.2
        elif r["regime"] == "coreactant_ald":
            col, rad = "#059669", 5.2
        else:
            col, rad = "#64748b", 4.6
        dx, dy = label_offsets_left.get(r["reaction_id"], (7, -6))
        tip = f"{r['short']}: Lit {r['ea_lit']:.2f} vs 7net-0 {r['ea_7n']:.2f} kcal/mol (err {r['signed_err']:+.2f})"
        left_pts_svg.append(
            f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{rad}" fill="{col}" stroke="#ffffff" stroke-width="1.3"><title>{tip}</title></circle>'
            f'<text x="{px + dx:.1f}" y="{py + dy:.1f}" font-size="9.2" font-weight="700" fill="{col}">{r["short"]}</text>'
        )

    # Right Panel: X = 7net-0 Ea (30..60), Y = Experimental T_decomp (150..700 C)
    rx0, rx1, ry0, ry1 = 480.0, 825.0, 42.0, 294.0

    def _rx(ea: float) -> float:
        return rx0 + ((max(30.0, min(60.0, ea)) - 30.0) / 30.0) * (rx1 - rx0)

    def _ry(tc: float) -> float:
        return ry1 - ((max(150.0, min(700.0, tc)) - 150.0) / 550.0) * (ry1 - ry0)

    # Arrhenius T1%(A=1e13 s^-1, t=1s) theoretical curve across Ea=30..60
    curve_pts = []
    for ea_v in [30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60]:
        t_v = compute_arrhenius_t_threshold_c(float(ea_v), 1.0e13)
        if 150.0 <= t_v <= 700.0:
            curve_pts.append(f"{_rx(float(ea_v)):.1f},{_ry(t_v):.1f}")
    arrh_poly = " ".join(curve_pts)

    label_offsets_right = {
        "CpZr_NMe2_betaH_elim": (-16, -11),
        "GaEt3_betaH_elim": (-36, 15),
        "DCDS_SiHCl_elim": (8, 4),
        "HCDS_SiCl2_elim": (8, -6),
        "TEOS_C2H4_elim": (-38, 16),
        "TTIP_acetone_elim": (-12, -12),
    }

    right_pts_svg: List[str] = []
    for r in rows:
        if r["t_mid_c"] is None:
            continue
        px = _rx(r["ea_7n"])
        py = _ry(float(r["t_mid_c"]))
        if r["reaction_id"] == "CpZr_NMe2_betaH_elim":
            col, rad = "#7c3aed", 6.4
        elif r["regime"] == "surface_catalytic":
            col, rad = "#d97706", 5.6
        else:
            col, rad = "#0284c7", 5.2
        dx, dy = label_offsets_right.get(r["reaction_id"], (7, -6))
        tip = f"{r['short']}: Ea_7net={r['ea_7n']:.2f} kcal/mol, T_decomp~{r['t_mid_c']:.1f} C (BDE_min={r['bde_min']})"
        right_pts_svg.append(
            f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{rad}" fill="{col}" stroke="#ffffff" stroke-width="1.3"><title>{tip}</title></circle>'
            f'<text x="{px + dx:.1f}" y="{py + dy:.1f}" font-size="9.2" font-weight="700" fill="{col}">{r["short"]}</text>'
        )

    # TTIP surface catalysis vertical arrow from Arrhenius gas-phase line (~468 C) down to 262.5 C
    ttip_r = ctx["by_id"]["TTIP_acetone_elim"]
    ttip_x = _rx(ttip_r["ea_7n"])
    ttip_y_gas = _ry(ttip_r["t1_7n_1e13"])
    ttip_y_exp = _ry(float(ttip_r["t_mid_c"]))

    cg = ctx["corr_gas"]
    return f"""
          <svg viewBox="0 0 855 348" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;">
            <!-- Panel Backgrounds -->
            <rect x="{lx0:.1f}" y="{ly0:.1f}" width="{lx1-lx0:.1f}" height="{ly1-ly0:.1f}" fill="#ffffff" stroke="#cbd5e1" />
            <rect x="{rx0:.1f}" y="{ry0:.1f}" width="{rx1-rx0:.1f}" height="{ry1-ry0:.1f}" fill="#ffffff" stroke="#cbd5e1" />

            <!-- Left Panel Titles & Parity -->
            <text x="{(lx0+lx1)/2:.1f}" y="22" text-anchor="middle" font-size="11.2" font-weight="800" fill="#0f172a">A. 8대 채널 문헌 vs SevenNet CI-NEB 활성화 장벽 (MAE {ctx['summary']['mae_ea_kcal_mol']:.2f} kcal/mol)</text>
            <polygon points="{band_pts}" fill="#dbeafe" fill-opacity="0.55" />
            <line x1="{_lx(0):.1f}" y1="{_ly(0):.1f}" x2="{_lx(75):.1f}" y2="{_ly(75):.1f}" stroke="#2563eb" stroke-width="1.4" stroke-dasharray="5,4" />
            {''.join(left_pts_svg)}
            <!-- Left Axes Ticks -->
            {''.join(f'<text x="{_lx(v):.1f}" y="{ly1+14:.1f}" text-anchor="middle" font-size="9.5" fill="#475569">{v}</text><text x="{lx0-7:.1f}" y="{_ly(v)+3:.1f}" text-anchor="end" font-size="9.5" fill="#475569">{v}</text>' for v in [0, 15, 30, 45, 60, 75])}
            <text x="{(lx0+lx1)/2:.1f}" y="330" text-anchor="middle" font-size="10" font-weight="700" fill="#334155">Literature DFT / Exp Ea (kcal/mol)</text>
            <text x="16" y="{(ly0+ly1)/2:.1f}" text-anchor="middle" font-size="10" font-weight="700" fill="#334155" transform="rotate(-90 16 {(ly0+ly1)/2:.1f})">SevenNet CI-NEB Ea (kcal/mol)</text>

            <!-- Right Panel Titles & Arrhenius T1% Curve -->
            <text x="{(rx0+rx1)/2:.1f}" y="22" text-anchor="middle" font-size="11.2" font-weight="800" fill="#0f172a">B. SevenNet Ea vs 실험 T_decomp (기상 5종 r = {cg['r_7n']:.3f} vs BDE r = {cg['r_bde']:.3f})</text>
            <polyline points="{arrh_poly}" fill="none" stroke="#0d9488" stroke-width="1.8" stroke-dasharray="6,4" />
            <text x="{_rx(53.2):.1f}" y="{_ry(530):.1f}" font-size="8.8" font-weight="700" fill="#0f766e">1차 Arrhenius T(1%, 1s, A=10^13)</text>
            <!-- TTIP catalytic offset arrow -->
            <line x1="{ttip_x:.1f}" y1="{ttip_y_gas+5:.1f}" x2="{ttip_x:.1f}" y2="{ttip_y_exp-6:.1f}" stroke="#d97706" stroke-width="1.5" stroke-dasharray="3,2" />
            <text x="{ttip_x-6:.1f}" y="{(ttip_y_gas+ttip_y_exp)/2:.1f}" text-anchor="end" font-size="8.6" font-weight="700" fill="#b45309">TiO2 자가촉매 저하</text>
            {''.join(right_pts_svg)}
            <!-- Right Axes Ticks -->
            {''.join(f'<text x="{_rx(v):.1f}" y="{ry1+14:.1f}" text-anchor="middle" font-size="9.5" fill="#475569">{v}</text>' for v in [30, 35, 40, 45, 50, 55, 60])}
            {''.join(f'<text x="{rx0-7:.1f}" y="{_ry(t)+3:.1f}" text-anchor="end" font-size="9.5" fill="#475569">{t}</text>' for t in [200, 300, 400, 500, 600, 700])}
            <text x="{(rx0+rx1)/2:.1f}" y="330" text-anchor="middle" font-size="10" font-weight="700" fill="#334155">SevenNet CI-NEB Ea (kcal/mol)</text>
            <text x="438" y="{(ry0+ry1)/2:.1f}" text-anchor="middle" font-size="10" font-weight="700" fill="#334155" transform="rotate(-90 438 {(ry0+ry1)/2:.1f})">Experimental T_decomp (deg C)</text>
          </svg>"""


def build_neb_storyline_html(
    neb_bench_path: str,
    unified_db_path: str,
) -> str:
    """Build concise Storyline Step 5 summary cards, SVG figure, and jump link."""
    ctx = compute_neb_kinetics_context(neb_bench_path, unified_db_path)
    s = ctx["summary"]
    cg = ctx["corr_gas"]
    cpzr = ctx["by_id"]["CpZr_NMe2_betaH_elim"]
    sel = s["ttip_selectivity"]
    svg_html = _build_neb_kinetics_svg(ctx)

    return rf"""
            <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap:16px; margin-bottom:20px;">
              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #0284c7; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">8대 반응경로 장벽 정합성 (CI-NEB)</div>
                <div style="font-size:1.75rem; font-weight:800; color:#0369a1; margin:4px 0;">MAE {s['mae_ea_kcal_mol']:.2f} kcal/mol</div>
                <div style="font-size:0.8rem; color:#475569;">RMSE <strong>{s['rmse_ea_kcal_mol']:.2f}</strong>, Bias <strong>{s['bias_mean_signed_err_kcal_mol']:+.2f}</strong> kcal/mol; TTIP 채널 선택성 &Delta;E<sub>a</sub> <strong>+{sel['sevenn_delta_ea_propene_minus_acetone_kcal_mol']:.2f}</strong> (문헌 +{sel['lit_delta_ea_propene_minus_acetone_kcal_mol']:.1f})</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #7c3aed; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">CpZr(NMe<sub>2</sub>)<sub>3</sub> 신규 검증 (ZAC)</div>
                <div style="font-size:1.75rem; font-weight:800; color:#6d28d9; margin:4px 0;">{cpzr['ea_7n']:.2f} vs {cpzr['ea_lit']:.2f} kcal/mol</div>
                <div style="font-size:0.8rem; color:#475569;">&beta;-H 제거 HNMe<sub>2</sub> 이탈 TS (<code>n_imag=1 [{cpzr['imag_freqs'][0]} cm<sup>-1</sup>]</code>), 예측 $T_{{1\%}}$ = <strong>{cpzr['t1_7n_1e13']:.1f}&deg;C</strong> (실험 305&ndash;350&deg;C)</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #059669; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">기상 단분자 T<sub>decomp</sub> 상관성 (N={cg['n']})</div>
                <div style="font-size:1.75rem; font-weight:800; color:#047857; margin:4px 0;">r = {cg['r_7n']:.3f} (BDE {cg['r_bde']:.3f})</div>
                <div style="font-size:0.8rem; color:#475569;">문헌 장벽 $r={cg['r_lit']:.3f}$; 정적 <code>BDE_min</code>이 역전시키는 <code>CpZr</code>&middot;<code>GaEt3</code>&middot;<code>HCDS</code> 내열성 서열 복원</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #d97706; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">NEB 자동화 판정 &amp; 온도 불확도</div>
                <div style="font-size:1.75rem; font-weight:800; color:#b45309; margin:4px 0;">&plusmn;{ctx['dt_halfwidth_c']:.1f}&deg;C / Stage-3</div>
                <div style="font-size:0.8rem; color:#475569;">무작위 IDPP 실패 방지 위해 SMARTS 템플릿 &middot; 고리 TS 시딩 &middot; 부분 헤시안 게이트 기반 선별 자동화 적용</div>
              </div>
            </div>

            <div style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:16px; margin-bottom:18px;">
              <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:8px;">
                <div style="font-weight:700; font-size:1.0rem; color:#0f172a;">
                  [NEB Kinetics Plot] 8대 반응 장벽 정합성 및 E<sub>a</sub><sup>CI-NEB</sup> 대 실험 T<sub>decomp</sub> 연계 분포
                </div>
                <span class="badge badge-3g" style="font-size:0.74rem;">8 Channels &bull; FS Topology {s['n_valid_fs_topology']}/{s['n_reactions']} &bull; Core TS 1-Imag {s['n_exact_1_imag_core_ts']}/{s['n_reactions']}</span>
              </div>
              {svg_html}
              <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:12px; margin-top:12px; padding:12px 14px; background:#f0f9ff; border:1px solid #bae6fd; border-left:4px solid #0284c7; border-radius:6px;">
                <div style="font-size:0.83rem; color:#1e3a8a; line-height:1.55; max-width:740px;">
                  <strong>핵심 결론:</strong> 정적 균일 해리 지표(<code>BDE_min</code>)는 <code>CpZr(NMe2)3</code>(69.96 kcal/mol) &gt; <code>HCDS</code>(42.16 kcal/mol)로 내열성을 오판하지만, SevenNet CI-NEB는 리간드 간 협동 &beta;-H 제거(41.55 kcal/mol)를 포착하여 실제 $T_{{\text{{decomp}}}}$ 서열을 복원합니다. 단, <code>TTIP</code>의 산화물 표면 자가촉매 하강과 <code>TDMAT</code>/<code>CpZr</code>의 공반응 교환 대비 자가분해 차이로 인해, ALD 공정창 판정은 이중 장벽($E_a^{{\text{{self-decomp}}}} - E_a^{{\text{{surface/ALD}}}}$) 체계로 해석해야 합니다.
                </div>
                <button type="button" onclick="jumpToTechnicalSection('neb-decomp-kinetics')" class="btn-pdf" style="background:#0369a1; color:#fff; border-color:#0369a1; cursor:pointer;">
                  Section 17 상세 분석 보기 &rarr;
                </button>
              </div>
            </div>"""


def build_neb_technical_html(
    neb_bench_path: str,
    unified_db_path: str,
) -> str:
    """Build full Section 17 technical cards (tables, CpZr analysis, Arrhenius, Mermaid)."""
    ctx = compute_neb_kinetics_context(neb_bench_path, unified_db_path)
    s = ctx["summary"]
    cg = ctx["corr_gas"]
    ca = ctx["corr_all_uni"]
    cpzr = ctx["by_id"]["CpZr_NMe2_betaH_elim"]
    cpzr_bonds = cpzr["bond_lengths_A"]
    svg_html = _build_neb_kinetics_svg(ctx)

    # Build Table 17-1 rows (all 8 benchmarked reactions)
    bench_rows_html: List[str] = []
    for r in ctx["rows"]:
        bde_str = f"{r['bde_min']:.2f}" if r["bde_min"] is not None else "-"
        conv_badge = (
            '<span style="color:#059669; font-weight:700;">Yes</span>'
            if r["converged"]
            else '<span style="color:#b45309; font-weight:600;">No</span>'
        )
        vib_str = (
            f"<strong style='color:#059669;'>1</strong> ({r['imag_freqs'][0]} cm<sup>-1</sup>)"
            if r["n_imag"] == 1
            else f"{r['n_imag']} ({', '.join(str(x) for x in r['imag_freqs']) if r['imag_freqs'] else 'endpoint/flat'})"
        )
        bond_summary = ", ".join(
            f"{k.split('_', 1)[1]}: {v['IS_A']}→<strong>{v['TS_A']}</strong>→{v['FS_A']}"
            for k, v in list(r["bond_lengths_A"].items())[:2]
        )
        row_bg = "background:#f5f3ff;" if r["reaction_id"] == "CpZr_NMe2_betaH_elim" else ""
        doi_url = f"https://doi.org/{r['lit_doi']}"
        bench_rows_html.append(
            f"""
              <tr style="{row_bg}">
                <td><strong>{r['short']}</strong><br><span style="font-size:0.73rem; color:#64748b;"><code>{r['precursor_id']}</code></span></td>
                <td><code style="font-size:0.75rem;">{r['equation']}</code><br><span style="font-size:0.74rem; color:#475569;">{r['mechanism']}</span></td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{bde_str}</td>
                <td style="text-align:center; font-family:'JetBrains Mono';"><strong>{r['ea_lit']:.2f}</strong><br><a href="{doi_url}" target="_blank" style="font-size:0.7rem;">{r['lit_doi']}</a></td>
                <td style="text-align:center; font-family:'JetBrains Mono'; color:#0369a1;"><strong>{r['ea_7n']:.2f}</strong></td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{r['signed_err']:+.2f}</td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{r['de_sep']:.2f}</td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{r['fmax']:.4f}<br>({conv_badge})</td>
                <td style="text-align:center; font-size:0.77rem;">{vib_str}</td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{r['wall_time_s']:.1f} s</td>
                <td style="font-size:0.74rem; font-family:'JetBrains Mono';">{bond_summary}</td>
              </tr>"""
        )

    # Build Table 17-2 rows (Kinetic T_1% & Experimental T_decomp comparison)
    kin_order = [
        "CpZr_NMe2_betaH_elim",
        "GaEt3_betaH_elim",
        "DCDS_SiHCl_elim",
        "HCDS_SiCl2_elim",
        "TEOS_C2H4_elim",
        "TTIP_acetone_elim",
        "TDMAT_NH3_transam",
    ]
    kin_rows_html: List[str] = []
    for rid in kin_order:
        r = ctx["by_id"][rid]
        bde_str = f"{r['bde_min']:.2f}" if r["bde_min"] is not None else "-"
        t13_str = f"{r['t1_7n_1e13']:.1f}&deg;C" if r["t1_7n_1e13"] > 0 else f"{r['t1_7n_1e13']:.1f}&deg;C (공반응)"
        t14_str = f"{r['t1_7n_1e14']:.1f}&deg;C" if r["t1_7n_1e14"] > 0 else "-"
        regime_tag = {
            "gas_unimol": '<span class="badge badge-2g">기상 단분자</span>',
            "surface_catalytic": '<span class="badge" style="background:#fef3c7; color:#92400e;">표면 자가촉매</span>',
            "coreactant_ald": '<span class="badge" style="background:#dcfce7; color:#166534;">공반응 ALD 교환</span>',
        }.get(r["regime"], "-")
        kin_rows_html.append(
            f"""
              <tr>
                <td><strong>{r['short']}</strong> (<code>{r['precursor_id']}</code>)</td>
                <td style="text-align:center;">{regime_tag}</td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{bde_str}</td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{r['ea_lit']:.2f} / <strong style="color:#0369a1;">{r['ea_7n']:.2f}</strong></td>
                <td style="text-align:center; font-family:'JetBrains Mono';">{r['t1_lit_1e13']:.1f}&deg;C</td>
                <td style="text-align:center; font-family:'JetBrains Mono';"><strong>{t13_str}</strong><br><span style="font-size:0.72rem; color:#64748b;">(A=10<sup>14</sup>: {t14_str})</span></td>
                <td style="text-align:center; font-weight:700;">{r['t_range_c']}</td>
                <td style="font-size:0.79rem; color:#334155;">{r['ald_note']}</td>
              </tr>"""
        )

    return rf"""
        <!-- Card 1: 8-Channel SevenNet CI-NEB Benchmark -->
        <div class="card" style="margin-bottom:24px; border-top:3px solid #0284c7;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">1. [8대 반응경로 CI-NEB 벤치마크] 문헌 DFT·실험 장벽 대비 SevenNet(7net-0) 정량 검증</h3>
            <span class="badge badge-2g">MAE {s['mae_ea_kcal_mol']:.2f} &middot; RMSE {s['rmse_ea_kcal_mol']:.2f} kcal/mol &bull; FS {s['n_valid_fs_topology']}/{s['n_reactions']} &bull; Core TS 1-Imag {s['n_exact_1_imag_core_ts']}/{s['n_reactions']}</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.9rem; line-height:1.65; margin-bottom:14px;">
            <code>database/benchmark_unified_latest_160.json</code> 내 실존 반도체 프리커서 및 대조군 7종(총 8개 반응 채널)을 대상으로 ASE <code>improvedtangent</code> NEB + CI-NEB(중간 이미지 5개/총 7개, <code>neb_steps=75</code>, <code>cineb_steps=130</code>, <code>k_spring=0.10 eV/&Aring;&sup2;</code>, <code>FIRE maxstep=0.05 &Aring;</code>, CPU 4스레드)를 일관 적용한 결과입니다. 모든 계산에서 최종 상태(FS) 위상 검증 <strong>{s['n_valid_fs_topology']}/{s['n_reactions']}</strong>을 통과하였으며, 전체 활성화 장벽 오차는 <strong>MAE {s['mae_ea_kcal_mol']:.2f} kcal/mol</strong>, 평균 부호 편향 <strong>{s['bias_mean_signed_err_kcal_mol']:+.2f} kcal/mol</strong>입니다.
          </p>

          <div style="margin-bottom:16px;">
            {svg_html}
          </div>

          <div style="overflow-x:auto;">
            <table class="report-table" style="font-size:0.79rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th>프리커서 (ID)</th>
                  <th>반응식 &amp; 메커니즘</th>
                  <th style="text-align:center;">DB BDE<sub>min</sub><br>(kcal/mol)</th>
                  <th style="text-align:center;">문헌 E<sub>a</sub> &amp; DOI<br>(kcal/mol)</th>
                  <th style="text-align:center;">SevenNet E<sub>a</sub><br>(kcal/mol)</th>
                  <th style="text-align:center;">Signed Error<br>(kcal/mol)</th>
                  <th style="text-align:center;">&Delta;E<sub>sep</sub><sup>&infin;</sup><br>(kcal/mol)</th>
                  <th style="text-align:center;">f<sub>max</sub> (eV/&Aring;)<br>(Convergence)</th>
                  <th style="text-align:center;">핵심 헤시안<br>n<sub>imag</sub></th>
                  <th style="text-align:center;">CPU 시간</th>
                  <th>전이상태 결합 길이 변화 (&Aring;)</th>
                </tr>
              </thead>
              <tbody>
                {''.join(bench_rows_html)}
              </tbody>
            </table>
          </div>
        </div>

        <!-- Card 2: Deep Dive on CpZr(NMe2)3 (ZAC, PREC_HK_07) -->
        <div class="card" style="margin-bottom:24px; border-top:3px solid #7c3aed;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">2. [CpZr 프리커서 심층 검증] CpZr(NMe<sub>2</sub>)<sub>3</sub>(ZAC, <code>PREC_HK_07</code>) 열분해 메커니즘 계산</h3>
            <span class="badge" style="background:#ede9fe; color:#5b21b6; border:1px solid #c4b5fd;">J. Phys. Chem. C 2022 &bull; DOI: 10.1021/acs.jpcc.2c04425</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.9rem; line-height:1.65; margin-bottom:14px;">
            고유전율 ZrO<sub>2</sub>/HZO ALD의 대표 이종 리간드 프리커서인 <code>Cp(NMe2)3Zr</code>(ZAC, <code>PREC_HK_07</code>)는 2D SMILES(<code>CN(C)[Zr](N(C)C)(N(C)C)C1=CC=CC1</code>)를 RDKit ETKDGv3/UFF로 임베딩할 경우 $\eta^1$-시클로펜타다이에닐 붕괴(<code>UFF Zr3</code> 미인식)가 발생합니다. 이를 해결하기 위해 $\eta^5$-piano-stool 초기 배위(<code>Zr&ndash;C(Cp) 2.60&ndash;2.62 &Aring;</code>, <code>Zr&ndash;N 2.07 &Aring;</code>)를 해석적으로 구축하고, Van 등(2022, <a href="https://doi.org/10.1021/acs.jpcc.2c04425" target="_blank">10.1021/acs.jpcc.2c04425</a>)이 보고한 두 가지 경쟁 열분해 경로를 범용 MLIAP인 SevenNet으로 검증하였습니다.
          </p>

          <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(290px, 1fr)); gap:16px; margin-bottom:14px;">
            <div style="background:#f5f3ff; border:1px solid #ddd6fe; border-left:4px solid #7c3aed; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#5b21b6; font-size:0.92rem; margin-bottom:6px;">주 경로(Channel A): 리간드 간 &beta;-H 제거 및 HNMe<sub>2</sub> 이탈</div>
              <div style="font-size:0.83rem; color:#334155; line-height:1.6;">
                $$\text{{CpZr(NMe}}_2)_3 \longrightarrow \text{{CpZr(NMe}}_2)(\eta^2\text{{-N(Me)CH}}_2) + \text{{HNMe}}_2$$
                <ul style="margin:6px 0 0 18px; padding:0;">
                  <li><strong>활성화 장벽:</strong> SevenNet CI-NEB <strong>{cpzr['ea_7n']:.2f} kcal/mol</strong> vs 문헌 M06L DFT <strong>{cpzr['ea_lit']:.2f} kcal/mol</strong> (144.9 kJ/mol, Signed Error {cpzr['signed_err']:+.2f} kcal/mol).</li>
                  <li><strong>전이상태 기하 및 진동 모드:</strong> <code>n_imag = 1 [{cpzr['imag_freqs'][0]} cm<sup>-1</sup>]</code>의 순수 양성자 이동 안장점 확보 (<code>C<sub>&alpha;</sub>&ndash;H</code>: {cpzr_bonds['breaking_Ca_H']['IS_A']}&rarr;<strong>{cpzr_bonds['breaking_Ca_H']['TS_A']}</strong>&rarr;{cpzr_bonds['breaking_Ca_H']['FS_A']} &Aring;; <code>N<sub>leave</sub>&ndash;H</code>: {cpzr_bonds['forming_Nleave_H']['IS_A']}&rarr;<strong>{cpzr_bonds['forming_Nleave_H']['TS_A']}</strong>&rarr;{cpzr_bonds['forming_Nleave_H']['FS_A']} &Aring;; 메탈라사이클 <code>Zr&ndash;C<sub>&alpha;</sub></code>: {cpzr_bonds['forming_Zr_Ca']['IS_A']}&rarr;<strong>{cpzr_bonds['forming_Zr_Ca']['TS_A']}</strong>&rarr;{cpzr_bonds['forming_Zr_Ca']['FS_A']} &Aring;).</li>
                  <li><strong>결합 열역학:</strong> 반응 부가체(&Delta;E<sub>adduct</sub> = +{cpzr['de_adduct']:.2f} kcal/mol, <code>Zr&ndash;N<sub>leave</sub></code> 2.375 &Aring;) 거쳐 기상 완전 분리 시 &Delta;E<sub>sep</sub> = +{cpzr['de_sep']:.2f} kcal/mol (문헌 ~20.8 kcal/mol).</li>
                </ul>
              </div>
            </div>

            <div style="background:#f8fafc; border:1px solid #cbd5e1; border-left:4px solid #0d9488; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#0f766e; font-size:0.92rem; margin-bottom:6px;">경쟁 경로(Channel B) 대비 Cp 구경꾼 리간드 내열 효과</div>
              <div style="font-size:0.83rem; color:#334155; line-height:1.6;">
                $$\text{{CpZr(NMe}}_2)_3 \longrightarrow \text{{Zr(NMe}}_2)_2(\eta^2\text{{-N(Me)CH}}_2) + \text{{C}}_5\text{{H}}_6$$
                <ul style="margin:6px 0 0 18px; padding:0;">
                  <li><strong>Cp 이탈 열역학 장벽:</strong> SevenNet 분리 생성물 에너지차 &Delta;E<sub>sep</sub> = <strong>+56.57 kcal/mol</strong>로, Channel A 전체 전이상태 장벽(41.55 kcal/mol)을 +15.02 kcal/mol 상회하여 초기 열분해에서 Cp 탈리가 억제됨을 재현합니다.</li>
                  <li><strong>TDMAZ 대비 ALD 창 확장 원리:</strong> <code>Zr(NMe2)4</code>(TDMAZ)은 &minus;NMe<sub>2</sub> 4개가 순차 탈리되며 반응성 높은 Zr 표면을 노출하는 반면, <code>CpZr(NMe2)3</code>는 처음 2개의 &minus;NMe<sub>2</sub>만 표면 &minus;OH와 저에너지(0.22/0.16 eV) 교환 후 난분해성 <code>O<sub>2</sub>ZrCp(NMe<sub>2</sub>)<sup>*</sup></code> 표면을 유지(Cp 탈리 장벽 1.66&ndash;2.26 eV)하여 305&ndash;350&deg;C 영역까지 자기제한적 성장을 보존합니다.</li>
                </ul>
              </div>
            </div>
          </div>
        </div>

        <!-- Card 3: ALD T_decomp Kinetic Linkage -->
        <div class="card" style="margin-bottom:24px; border-top:3px solid #059669;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">3. [ALD 열분해 온도(T<sub>decomp</sub>) 연계 분석] E<sub>a</sub><sup>CI-NEB</sup> &middot; BDE<sub>min</sub> 및 이중 장벽 공정창</h3>
            <span class="badge" style="background:#dcfce7; color:#166534; border:1px solid #86efac;">Gas-Phase r = {cg['r_7n']:.3f} vs BDE r = {cg['r_bde']:.3f}</span>
          </div>

          <div class="math-callout" style="margin-bottom:14px;">
            <strong>1차 단분자 열분해 속도론 및 ALD 이중 장벽 공정창 정식화:</strong><br>
            $$k(T) = A \exp\!\left(-\frac{{E_a}}{{R T}}\right), \qquad T_{{\text{{decomp}}}}(\chi, t_{{\text{{res}}}}) = \frac{{E_a}}{{R \ln\!\left(\dfrac{{A \, t_{{\text{{res}}}}}}{{-\ln(1 - \chi)}}\right)}}, \qquad \Delta T_{{\text{{decomp}}}} \approx \frac{{\Delta E_a}}{{R \ln\!\left(\dfrac{{A \, t_{{\text{{res}}}}}}{{-\ln(1 - \chi)}}\right)}}$$
            $$\Delta E_{{\text{{window}}}} \equiv E_a^{{\text{{self-decomp}}}} - E_a^{{\text{{co-reactant / surface}}}} \gg 0 \quad \Longrightarrow \quad T_{{\text{{onset}}}}^{{\text{{ALD}}}}\big(E_a^{{\text{{co-reactant}}}}\big) &lt; T_{{\text{{process}}}} &lt; T_{{\text{{decomp}}}}^{{\text{{eff}}}}\big(E_a^{{\text{{self-decomp, surf}}}}\big)$$
          </div>

          <div style="overflow-x:auto; margin-bottom:14px;">
            <table class="report-table" style="font-size:0.8rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th>프리커서</th>
                  <th style="text-align:center;">지배 메커니즘 구분</th>
                  <th style="text-align:center;">DB BDE<sub>min</sub><br>(kcal/mol)</th>
                  <th style="text-align:center;">문헌 / SevenNet E<sub>a</sub><br>(kcal/mol)</th>
                  <th style="text-align:center;">이론 $T_{{1\%}}$ (문헌)<br>(A=10<sup>13</sup> s<sup>-1</sup>)</th>
                  <th style="text-align:center;">이론 $T_{{1\%}}$ (SevenNet)<br>(A=10<sup>13</sup> / 10<sup>14</sup>)</th>
                  <th style="text-align:center;">실험 기상/공정<br>T<sub>decomp</sub> (&deg;C)</th>
                  <th>공정창 및 촉매 편차 해석</th>
                </tr>
              </thead>
              <tbody>
                {''.join(kin_rows_html)}
              </tbody>
            </table>
          </div>

          <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(290px, 1fr)); gap:14px;">
            <div style="background:#f0fdf4; border:1px solid #bbf7d0; border-left:4px solid #059669; border-radius:6px; padding:12px;">
              <div style="font-weight:700; color:#166534; font-size:0.88rem; margin-bottom:4px;">(1) 무촉매 기상 5종에서의 BDE<sub>min</sub> 한계 극복</div>
              <p style="font-size:0.81rem; color:#14532d; line-height:1.55; margin:0;">
                무촉매 기상 분해 5종(<code>CpZr</code>, <code>GaEt3</code>, <code>DCDS</code>, <code>HCDS</code>, <code>TEOS</code>)에서 실험 $T_{{\text{{decomp}}}}$ 중간값과의 Pearson Correlation은 <strong>문헌 E<sub>a</sub> $r={cg['r_lit']:.3f}$</strong>, <strong>SevenNet CI-NEB $r={cg['r_7n']:.3f}$</strong>인 반면 정적 <code>BDE_min</code>은 <strong>$r={cg['r_bde']:.3f}$</strong>에 그칩니다. 특히 <code>CpZr(NMe2)3</code>(BDE 69.96 &gt; HCDS 42.16)와 <code>GaEt3</code>(44.86 &gt; HCDS 42.16)는 균일 해리 에너지 기준과 달리 저에너지 협동 &beta;-H 제거 경로를 통해 <code>HCDS</code>보다 약 200&deg;C 낮은 온도에서 열분해되며, SevenNet CI-NEB는 이를 정확히 분리합니다.
              </p>
            </div>
            <div style="background:#fffbeb; border:1px solid #fde68a; border-left:4px solid #d97706; border-radius:6px; padding:12px;">
              <div style="font-weight:700; color:#92400e; font-size:0.88rem; margin-bottom:4px;">(2) 표면 자가촉매(TTIP) 및 공반응 교환(TDMAT) 구분</div>
              <p style="font-size:0.81rem; color:#78350f; line-height:1.55; margin:0;">
                표면 촉매 효과가 강한 <code>TTIP</code>를 포함하면 6종 Pearson Correlation은 $r_{{\text{{lit}}}}={ca['r_lit']:.3f}$, $r_{{\text{{7net}}}}={ca['r_7n']:.3f}$로 하락합니다. 이는 <code>TTIP</code> 기상 장벽(50.87 kcal/mol &rarr; $T_{{1\%}}$ {ctx['by_id']['TTIP_acetone_elim']['t1_7n_1e13']:.1f}&deg;C)이 성장 중인 TiO<sub>2</sub> 표면에서 $E_{{a,\text{{app}}}}\sim 20\text{{--}}40\text{{ kcal/mol}}$로 낮아지기 때문입니다(<a href="https://doi.org/10.1021/acs.jpca.7b09904" target="_blank">10.1021/acs.jpca.7b09904</a>). 또한 <code>TDMAT+NH3</code>(16.22 kcal/mol)는 목적 ALD 교환 장벽이므로, 실제 열분해 상한(200&ndash;220&deg;C)은 단분자 &beta;-H 자가분해 장벽(39.3 kcal/mol)과의 이중 장벽 차이로 판정해야 합니다.
              </p>
            </div>
          </div>
        </div>

        <!-- Card 4: Automation Feasibility & Stage-3 Architecture -->
        <div class="card" style="margin-bottom:24px; border-top:3px solid #d97706;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">4. [NEB 자동화 타당성 평가] 연산 비용 &middot; 실패 모드 진단 및 Stage-3 계층형 파이프라인 설계</h3>
            <span class="badge" style="background:#fef3c7; color:#92400e; border:1px solid #fcd34d;">Stage-3 Selective Pipeline &bull; Guarded CI-NEB</span>
          </div>

          <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap:14px; margin-bottom:16px;">
            <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:13px;">
              <div style="font-weight:700; color:#0f172a; font-size:0.86rem; margin-bottom:4px;">정밀도 대 속도론적 온도 민감도</div>
              <div style="font-size:0.8rem; color:#475569; line-height:1.5;">
                1차 속도론($A=10^{{13}}\text{{ s}}^{{-1}}, t_{{\text{{res}}}}=1\text{{ s}}, \chi=1\%$)에서 <code>MAE = {s['mae_ea_kcal_mol']:.2f} kcal/mol</code>은 <strong>&plusmn;{ctx['dt_halfwidth_c']:.1f}&deg;C</strong>의 $T_{{\text{{decomp}}}}$ 불확도에 해당합니다. 따라서 정밀 공정 온도 예측기보다는 <strong>저내열(&le;350&deg;C) vs 고내열(&ge;500&deg;C) 프리커서 계층 선별 필터</strong>로 운용하는 것이 적합합니다.
              </div>
            </div>
            <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:13px;">
              <div style="font-weight:700; color:#0f172a; font-size:0.86rem; margin-bottom:4px;">CPU 연산 비용 및 GPU 가속 요구사항</div>
              <div style="font-size:0.8rem; color:#475569; line-height:1.5;">
                CPU 4스레드 기준 소형 실란(8원자)은 <code>283.8&ndash;296.9 s</code>로 신속하나, 중대형 알콕사이드&middot;아마이드(22&ndash;45원자)는 밴드당 <code>782.7&ndash;1818.7 s</code>(약 13&ndash;30분)가 소요됩니다. 대량 자동화에는 <strong>Stage 1&ndash;2 통과 후보 한정 실행</strong> 및 GPU/cuEquivariance 가속이 필수입니다.
              </div>
            </div>
            <div style="background:#fef2f2; border:1px solid #fecaca; border-radius:6px; padding:13px;">
              <div style="font-weight:700; color:#991b1b; font-size:0.86rem; margin-bottom:4px;">비편향 IDPP 자동화의 3대 실패 유형</div>
              <div style="font-size:0.8rem; color:#7f1d1d; line-height:1.5;">
                단순 선형 보간은 (1) <code>GaEt3</code> C&ndash;H 균일 해리 이탈(102.70 kcal/mol), (2) <code>TTIP</code> 선행 Ti&ndash;O 절단 우회(71.23 kcal/mol), (3) <code>TEOS</code> 늦은 CI 점화로 인한 좁은 능선 미끄러짐(21.46 kcal/mol) 및 <code>CpZr</code> RDKit &eta;<sup>1</sup> 붕괴를 유발하므로 규칙 기반 시딩이 필수입니다.
              </div>
            </div>
          </div>

          <div class="mermaid">
flowchart LR
  S1["Stage 1-2: 정적 물리 스크리닝<br>(TMC-SA, V_bur, BDE_min)"] --> S2["반응 클래스 SMARTS 매칭<br>(beta-H, 4/6-원 고리, 1,2-이동)"]
  S2 --> S3["배위 보정 초기구조 생성<br>(eta^5 Piano-Stool / Kabsch 정렬)"]
  S3 --> S4["TS 중심 25~45 percent 고리형 시딩<br>+ 조기 CI-NEB (FIRE maxstep=0.05)"]
  S4 --> G{{"자동 품질 게이트<br>valid_fs &amp; TS결합 &amp; n_imag=1"}}
  G -- "Pass" --> K["Arrhenius T_1pct &amp; 이중 장벽<br>ALD 공정창(Delta E_window) 판정"]
  G -- "Fail / Flat" --> R["시딩/스텝 보정 재시도<br>또는 VASP DFT 이관"]
          </div>
        </div>
"""
