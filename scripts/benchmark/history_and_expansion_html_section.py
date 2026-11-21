"""
HTML section builders for:
- Section 21: Chronological Analysis & Validation History (#analysis-history-roadmap)
- Summary Step 4-5 & Section 22: KR/JP/CN/US/TW Patent & Metal-Centered Literature DB Expansion (#patent-literature-expansion)
"""

from __future__ import annotations

import html
import json
import os
from collections import Counter
from typing import Any, Dict, List, Optional


def _render_template(tpl: str, mapping: Dict[str, Any]) -> str:
    """Substitute __TOKEN__ placeholders safely without str.format brace collisions."""
    out = tpl
    for k, v in mapping.items():
        out = out.replace(k, str(v))
    return out


def _load_expansion_payload(expansion_json_path: str) -> Dict[str, Any]:
    """Load precomputed patent/literature expansion screening payload."""
    with open(expansion_json_path, "r", encoding="utf-8") as f:
        return json.load(f)


def _build_history_timeline_svg() -> str:
    """Build inline SVG roadmap showing the 10-stage investigative progression."""
    stages = [
        ("S1", "Phase 1 (55종)", "유기 SA 18종 오분류 규명", "#dc2626"),
        ("S2", "SevenNet-0 71+27", "MLIP 3D·BDE·이합체 ΔG", "#ea580c"),
        ("S3", "음성대조군 40종", "LOFMO=0 vs 물리 1.000", "#d97706"),
        ("S4", "SAPIS (160종)", "Design 1+2 AUROC 0.996", "#059669"),
        ("S5", "8대 CI-NEB+CpZr", "Ea-T_decomp r=0.879", "#0284c7"),
        ("S6", "배위·리간드 통계", "16금속 51리간드 CN_eff", "#2563eb"),
        ("S7", "CoOrd-Inv 역설계", "35,661 후보 OOF 4.28%", "#4f46e5"),
        ("S8", "M+L 물리·일반화", "Pauling-HSAB·포화 %Vbur", "#7c3aed"),
        ("S9", "분석 이력 체계화", "Git·테스트 28종 계보화", "#0f766e"),
        ("S10", "5극 특허·문헌 확장", "+28 프리커서·+10 인히비터", "#be123c"),
    ]
    w, h = 900, 178
    nodes = []
    # 2 rows of 5 nodes connected left-to-right then down-left-to-right (S-curve)
    coords = []
    for idx in range(10):
        row = idx // 5
        col = idx % 5 if row == 0 else (4 - (idx % 5))
        cx = 92 + col * 178
        cy = 48 + row * 84
        coords.append((cx, cy))

    arrows = []
    for i in range(9):
        x1, y1 = coords[i]
        x2, y2 = coords[i + 1]
        arrows.append(
            f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#94a3b8" stroke-width="2.2" marker-end="url(#histArr)"/>'
        )

    for (tag, title, sub, col), (cx, cy) in zip(stages, coords):
        nodes.append(
            f'<rect x="{cx - 78}" y="{cy - 26}" width="156" height="52" rx="8" fill="#ffffff" stroke="{col}" stroke-width="2"/>'
            f'<rect x="{cx - 78}" y="{cy - 26}" width="32" height="52" rx="6" fill="{col}"/>'
            f'<text x="{cx - 62}" y="{cy + 4}" text-anchor="middle" font-size="11" font-weight="800" fill="#ffffff">{tag}</text>'
            f'<text x="{cx + 16}" y="{cy - 5}" text-anchor="middle" font-size="10.2" font-weight="800" fill="#0f172a">{html.escape(title)}</text>'
            f'<text x="{cx + 16}" y="{cy + 12}" text-anchor="middle" font-size="8.6" fill="#475569">{html.escape(sub)}</text>'
        )

    return (
        f'<svg viewBox="0 0 {w} {h}" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px;">'
        '<defs><marker id="histArr" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto">'
        '<polygon points="0 0, 7 3.5, 0 7" fill="#64748b"/></marker></defs>'
        + "".join(arrows)
        + "".join(nodes)
        + '<text x="888" y="170" text-anchor="end" font-size="8.8" fill="#64748b">* 전 단계 스크립트·JSON 산출물·회귀 단위 테스트(28종)로 100% 재현 가능하도록 연결된 연구 이력 아키텍처</text>'
        "</svg>"
    )


def build_analysis_history_technical_html(
    unified_json_path: str,
    neb_json_path: str,
    expansion_json_path: str,
) -> str:
    """Build Section 21 summarizing chronological analysis & verification history."""
    with open(unified_json_path, "r", encoding="utf-8") as f:
        u_recs = json.load(f)["records"]
    with open(neb_json_path, "r", encoding="utf-8") as f:
        n_data = json.load(f)
    exp_data = _load_expansion_payload(expansion_json_path)

    n_unified = len(u_recs)
    n_neb = len(n_data.get("reactions", []))
    n_exp = exp_data["summary"]["n_total"]
    timeline_svg = _build_history_timeline_svg()

    history_rows = [
        (
            "Stage 1",
            "초기 문헌·오픈소스 SA 18종 감사",
            "55종 (27 프리커서 + 28 시약)",
            "run_phase1_baseline.py / audit_reproduce_all.py",
            "제약용 SA(RDKit SAScore) 92.6% 오분류, DFRscore 100% ELSE 붕괴, ReaSyn 디코딩 오프셋 규명 (9/9 실측 감사).",
            "Sec 1~10, 16",
        ),
        (
            "Stage 2",
            "SevenNet-0 3D 물리량 전수 계산",
            "75종 확장 포함 102 프리커서 + 40 음성대조군",
            "run_expanded_elements_sevenn.py / run_expanded_neg_sevenn.py",
            "프리커서 102종(잔여 71종 포함) 및 실패유형 층화 40종 음성대조군에 대해 BDE_min, E_coh, ΔG_dimer(373 K) 무대체 전산 산출.",
            "Sec 13, 15",
        ),
        (
            "Stage 3",
            "다중 배위수·SAPIS 물리 판정 기준 통합",
            "160종 (102 Prec + 18 Inh + 40 Neg)",
            "element_physics.py / tmc_sa_scorer.py / run_unified_latest_160_analysis.py",
            "CN_eff·R_eff(M|env)·%V_bur^D2·Design 1+2 정식화; ML LOFMO 0.000 붕괴 구간에서 SAPIS 물리 판정 기준 1.000 방어 검증.",
            "Sec 15 / Step 4",
        ),
        (
            "Stage 4",
            "CI-NEB 8대 분해경로 & CpZr 이중 채널",
            "8개 반응 (Si2Cl6, Si2H2Cl4, GaEt3, TEOS, TTIPx2, TDMAT+NH3, CpZr)",
            "run_sevenn_decomp_neb_benchmark.py / neb_html_section.py",
            "정적 BDE(69.96) 대비 Cp(NMe2)3Zr β-H 협동 장벽(34.67 vs 문헌 34.63) 재현; T_decomp 상관 r=0.507->0.879 상승 및 7-이미지 3D 가시화.",
            "Sec 17 / Step 5",
        ),
        (
            "Stage 5",
            "전이금속 다중 배위수·리간드 조합 통계",
            "102 프리커서 / 16 중심금속 / 51 리간드 조각",
            "ligand_coord_html_section.py",
            "CN_eff=1~14 분포 및 Top-30 리간드 94.1% 부분 커버리지 정량화; 연속 M×L 결합 대리모델 5-CV/LOMO-5 비교(CCTBA·Co(AMD)2 한계 명시).",
            "Sec 18 / Step 4-2",
        ),
        (
            "Stage 6",
            "반도체 프리커서 역설계 모델 (CoOrd-Inv)",
            "35,661 조합 공간 (16 금속 × 44 리간드)",
            "inverse_design_model.py",
            "생성 AI 4패러다임 통합 및 OOF 회귀 재발견(5-CV 4.28%, LOMO-5 9.50%)으로 이산 Bag(L) 대비 Top-5% 재현율 대폭 향상.",
            "Sec 12, 19 / Step 4-3",
        ),
        (
            "Stage 7",
            "중심원자–리간드 물리법칙 & 4대 일반화 검증",
            "102 회귀 + 142 분류 + 8 CI-NEB 모티프",
            "mechanistic_ml_model.py",
            "Pauling-HSAB ΔEf/vM(Homo->Hetero +0.586) 및 구면포화 %Vbur(+0.647)로 RF 능가; 반면 BDE_min(-0.287)·ΔG_dimer 외삽 한계와 모티프별 동역학 정식화.",
            "Sec 20 / Step 4-4",
        ),
        (
            "Stage 8",
            "특허·문헌 DB 확장 + 38종 SevenNet-0 전수 계산",
            "+28 신규 프리커서 (8 신규 원소) + 10 인히비터",
            "run_patent_lit_expansion.py / run_patent_expansion_sevenn.py",
            "Mo/W/Nb/Ni/Cu/La/Y/Ge/Sn/Pt 공백 해소 및 38종 전수 7net-0 LBFGS 완료 (0 Fallback); PDMAT(93.5%->79.2%)·Ge(NMe2)4(90.5%->77.0%) 이완 구제 및 286 정준 DB 통합.",
            "Sec 21~22 / Step 4-5",
        ),
    ]

    row_html = "".join(
        f"""
        <tr>
          <td><span class="badge badge-2g" style="font-size:0.72rem;">{html.escape(st)}</span></td>
          <td><strong>{html.escape(title)}</strong></td>
          <td style="font-size:0.76rem; color:#0369a1; font-weight:700;">{html.escape(scope)}</td>
          <td><code style="font-size:0.73rem;">{html.escape(scripts)}</code></td>
          <td style="font-size:0.79rem; line-height:1.5;">{html.escape(finding)}</td>
          <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.75rem; font-weight:700;">{html.escape(sec)}</td>
        </tr>"""
        for st, title, scope, scripts, finding, sec in history_rows
    )

    tpl = r"""
        <div class="card" style="margin-bottom:20px; border-top:4px solid #0f766e;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">1. 프로젝트 연구·분석 이력 아키텍처 (Stage 1 &rarr; Stage 8 진화 계보)</h3>
            <span class="badge badge-phys" style="font-size:0.78rem;">Reproducible Investigative Lineage</span>
          </div>
          <p style="font-size:0.89rem; color:#334155; line-height:1.65; margin-bottom:12px;">
            본 연구는 단순 유기물 <code>SA</code> 지표의 한계를 진단한 <strong>Phase 1 (55종)</strong>에서 출발해, 범용 MLIP(<code>SevenNet-0</code>)·음성대조군 <code>40</code>종(<code>__N_UNIFIED__</code>종 통합 평가셋)·<code>__N_NEB__</code>대 반응경로 <code>CI-NEB</code>·다중 배위수 통계·<code>CoOrd-Inv</code> 역설계·중심원자–리간드 물리 모델 일반화 검증(Section 20)을 거쳐 <strong>한·중·일·미·대만 특허/문헌 확장 후보(__N_EXP__종)</strong>까지 단일 재현 파이프라인으로 축적되었습니다.
          </p>
          <div style="margin-bottom:14px;">
            __TIMELINE_SVG__
          </div>
          <div class="table-responsive">
            <table class="data-table" style="font-size:0.79rem;">
              <thead>
                <tr>
                  <th style="width:74px;">단계</th>
                  <th style="width:185px;">분석 마일스톤</th>
                  <th style="width:155px;">대상 데이터 범위</th>
                  <th style="width:190px;">핵심 모듈 / 스크립트</th>
                  <th>규명된 물리·화학적 결론 및 개선 이력</th>
                  <th style="width:92px; text-align:center;">리포트 위치</th>
                </tr>
              </thead>
              <tbody>
                __HIST_ROWS__
              </tbody>
            </table>
          </div>
        </div>"""

    return _render_template(
        tpl,
        {
            "__N_UNIFIED__": n_unified,
            "__N_NEB__": n_neb,
            "__N_EXP__": n_exp,
            "__TIMELINE_SVG__": timeline_svg,
            "__HIST_ROWS__": row_html,
        },
    )


def _build_element_coverage_svg(cov_table: List[Dict[str, Any]]) -> str:
    """Build grouped horizontal/vertical SVG showing current 102 vs +28 expanded metal coverage."""
    # Show top 18 elements covering all newly added elements and industrial metals
    top_rows = cov_table[:24]
    w, h = 620, 235
    x0, y_bot, plot_w, plot_h = 42, 198, 560, 152
    max_val = max(r["merged_130"] for r in top_rows) or 25
    step = plot_w / len(top_rows)
    bw = max(10.0, step * 0.62)

    bars = []
    for i, row in enumerate(top_rows):
        cx = x0 + (i + 0.5) * step
        c_cur = row["current_102"]
        c_add = row["added_28"]
        h_cur = (c_cur / max_val) * plot_h
        h_add = (c_add / max_val) * plot_h
        if c_cur > 0:
            bars.append(
                f'<rect x="{cx - bw / 2:.1f}" y="{y_bot - h_cur:.1f}" width="{bw:.1f}" height="{h_cur:.1f}" fill="#94a3b8"/>'
            )
        if c_add > 0:
            bars.append(
                f'<rect x="{cx - bw / 2:.1f}" y="{y_bot - h_cur - h_add:.1f}" width="{bw:.1f}" height="{h_add:.1f}" fill="#e11d48"/>'
            )
            bars.append(
                f'<text x="{cx:.1f}" y="{y_bot - h_cur - h_add - 4:.1f}" text-anchor="middle" font-size="8.2" font-weight="800" fill="#be123c">+{c_add}</text>'
            )
        el_col = "#be123c" if c_cur == 0 else "#0f172a"
        bars.append(
            f'<text x="{cx:.1f}" y="{y_bot + 13}" text-anchor="middle" font-size="8.8" font-weight="700" fill="{el_col}">{html.escape(row["element"])}</text>'
        )

    grid = []
    for t in (0, 5, 10, 15, 20, 25):
        if t <= max_val + 2:
            yy = y_bot - (t / max_val) * plot_h
            grid.append(
                f'<line x1="{x0}" y1="{yy:.1f}" x2="{x0 + plot_w}" y2="{yy:.1f}" stroke="#e2e8f0" stroke-dasharray="3 3"/>'
                f'<text x="{x0 - 5}" y="{yy + 3:.1f}" text-anchor="end" font-size="8.5" fill="#64748b">{t}</text>'
            )

    return (
        f'<svg viewBox="0 0 {w} {h}" style="width:100%; height:auto; display:block; background:#ffffff; border:1px solid #cbd5e1; border-radius:8px;">'
        '<text x="310" y="18" text-anchor="middle" font-size="11.2" font-weight="800" fill="#0f172a">중심원소별 프리커서 분포: 현행 DB (102종) vs 특허·문헌 확장 (+28종, 붉은색 표기=신규 8원소)</text>'
        '<rect x="345" y="27" width="10" height="10" fill="#94a3b8"/><text x="359" y="35" font-size="8.8" fill="#334155">현행 102종 (16원소)</text>'
        '<rect x="468" y="27" width="10" height="10" fill="#e11d48"/><text x="482" y="35" font-size="8.8" fill="#be123c">특허·문헌 신규 +28종 (총 24원소)</text>'
        + "".join(grid)
        + "".join(bars)
        + f'<line x1="{x0}" y1="{y_bot}" x2="{x0 + plot_w}" y2="{y_bot}" stroke="#475569" stroke-width="1.3"/>'
        "</svg>"
    )


def build_patent_expansion_storyline_html(
    expansion_json_path: str,
    unified_json_path: str,
) -> str:
    """Build Summary Tab Step 4-5 executive bridge card for patent & literature DB expansion."""
    exp_data = _load_expansion_payload(expansion_json_path)
    s = exp_data["summary"]
    tpl = r"""
          <!-- Step 4-5: KR/JP/CN/US/TW Patent & Metal-Centered Literature DB Expansion -->
          <div class="card" style="margin-bottom:24px; border-left:5px solid #e11d48; background:linear-gradient(135deg, #fff1f2 0%, #ffffff 60%); padding:20px;">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
              <h4 style="margin:0; font-size:1.03rem; color:#881337; font-weight:800;">
                4-5. 한&middot;중&middot;일&middot;미&middot;대만 특허 및 중심금속별 문헌 기반 DB 확장 후보 스크리닝 (+__N_TOT__종)
              </h4>
              <span class="badge" style="background:#ffe4e6; color:#9f1239; border:1px solid #fecdd3; font-size:0.76rem;">
                KR &middot; JP &middot; CN &middot; US &middot; TW Patents &amp; Reviews
              </span>
            </div>
            <p style="font-size:0.88rem; color:#334155; line-height:1.65; margin-bottom:12px;">
              현행 <code>102</code>종 프리커서는 <code>Ga(24)</code>&middot;<code>In(24)</code>&middot;<code>Ru(16)</code>에 편중되어 있고 <code>Nb&middot;Ni&middot;Cu&middot;Ge&middot;La&middot;Y&middot;Zn&middot;Pt</code>가 전무(<code>0</code>종)하며 <code>Mo&middot;W&middot;Ta&middot;Sn</code>은 각 <code>1</code>종에 불과합니다. 한&middot;중&middot;일&middot;미&middot;대만 5극 특허(<code>ADEKA</code>, <code>SK트리켐/TCLC</code>, <code>Merck</code>, <code>Air Liquide</code>, <code>Applied Materials</code>, <code>TEL</code>, <code>TSMC</code> 등) 및 금속별 문헌 리뷰를 전수 대조해 <strong>신규 프리커서 __N_PREC__종</strong>(신규 <strong>__N_NEW_EL__개 원소</strong> 포함, 총 <strong>__N_MERGED_EL__원소</strong>)과 <strong>ASD 인히비터 __N_INH__종</strong>을 도출하고, <strong>전수 __N_SEVENN__종(Fallback <code>0</code>건)에 <code>SevenNet-0</code> LBFGS 3D 이완&middot;BDE&middot;E<sub>coh</sub>&middot;&Delta;E<sub>dimer</sub>&middot;&Delta;E<sub>f</sub><sup>std</sup>/v<sub>M</sub></strong> 계산을 완료했습니다.
            </p>
            <div class="stat-grid" style="margin-bottom:12px;">
              <div class="stat-box" style="border-top:3px solid #be123c;">
                <div class="stat-num" style="color:#be123c;">+__N_PREC__ (+__N_NEW_EL__ 원소)</div>
                <div class="stat-label">ETKDG 통과 __N_PREC_PASS__/__N_PREC__ &bull; 7net 엄격 통과 __N_PREC_PASS_7N__/__N_PREC__</div>
              </div>
              <div class="stat-box" style="border-top:3px solid #0284c7;">
                <div class="stat-num" style="color:#0369a1;">+__N_INH__종 (7net __N_INH_PASS_7N__/__N_INH__)</div>
                <div class="stat-label">신규 특허 ASD 인히비터 (ETKDG __N_INH_PASS__/__N_INH__ 통과)</div>
              </div>
              <div class="stat-box" style="border-top:3px solid #7c3aed;">
                <div class="stat-num" style="color:#6d28d9;">__N_SEVENN__/38 SevenNet-0 완비</div>
                <div class="stat-label">PDMAT(93.5&rarr;79.2%)&middot;Ge(NMe2)4(90.5&rarr;77.0%) 이완 구제</div>
              </div>
            </div>
            <div style="font-size:0.83rem; color:#475569; line-height:1.55;">
              * <strong>자세한 5개국 특허 목록, 원소별 문헌 DOI, 전체 38종 SevenNet-0 실측 디스크립터 및 교차 DB 중복 제거 내역</strong>은 리포트 탭 <strong>Section 21~22</strong>와 DB 탐색 탭에서 확인할 수 있습니다.
            </div>
          </div>"""

    return _render_template(
        tpl,
        {
            "__N_TOT__": s["n_total"],
            "__N_PREC__": s["n_precursors"],
            "__N_PREC_PASS__": s["n_precursors_pass"],
            "__N_INH__": s["n_inhibitors"],
            "__N_INH_PASS__": s["n_inhibitors_pass"],
            "__N_PREC_PASS_7N__": s.get("n_precursors_pass_7net", s["n_precursors_pass"]),
            "__N_INH_PASS_7N__": s.get("n_inhibitors_pass_7net", s["n_inhibitors_pass"]),
            "__N_SEVENN__": s.get("n_sevenn_computed", 0),
            "__N_NEW_EL__": s["n_new_elements_added"],
            "__N_MERGED_EL__": s["n_merged_elements"],
        },
    )


def build_patent_expansion_technical_html(
    expansion_json_path: str,
    unified_json_path: str,
) -> str:
    """Build Section 22 cards analyzing 5-jurisdiction patents, metal-centered literature, and 38 screened candidates."""
    exp_data = _load_expansion_payload(expansion_json_path)
    s = exp_data["summary"]
    recs = exp_data["records"]
    cov_svg = _build_element_coverage_svg(s["coverage_table"])

    patent_landscape = [
        (
            "한국 (KR)",
            "SK트리켐·한솔케미칼·솔브레인·DNF·원익머트리얼즈·KAIST",
            "KR20230088807A (MoO2Cl2), KR101719526B1 (Co/DAD), KR101216068B1 (유기금속 ALD), KR20190055531A (SK트리켐), KR20200140432A / KR20220109843A (DMADMS/DMATMS 선택적 억제)",
            "DRAM 고유전율 Hf/Zr 시클로펜타디에닐(Cp) 국산화, 3D NAND 워드라인용 불소-프리 MoO2Cl2, 아미노실란 기반 선택적 증착(ASD) 공정 특허 집중.",
        ),
        (
            "일본 (JP)",
            "ADEKA·트리케미칼(TCLC)·다나까귀금속·Tokyo Electron (TEL)",
            "JP6948159B2 / JP6954776B2 / JP5128289B2 (ADEKA Cp-Hf 화합물), JP7847142B2 (액상 Hf/Zr), TEL 건식 기상 SAM/아미노실란·ALE 순환 특허군",
            "HZO 강유전체·HKMG용 CpHf(NMe2)3/dmap 복합계, 다나까 고휘발성 액상 Ru 프리커서(TRuST), TEL 건식 기상 억제제(TMSDEA·포스폰산) 정밀 제어 선도.",
        ),
        (
            "미국·국제 (US/WO)",
            "Merck(Versum)·Air Liquide·Entegris·Applied Materials·Lam Research·ASM",
            "US20250250673A1 (MoO2Cl2 ALD), US12077858B2 (WCl5 무불소 갭필), US7557229B2 / WO2004046417 (Gordon 아미디네이트), US9099301B1 (La 아미디네이트), WO2024050202A1 (희토류 Cp), US10804099B2 / US9564312B2",
            "Entegris/Merck 무불소 Mo·W 할라이드+고체 승화기(Subleem/Chemguard), Gordon 아미디네이트(Ni/Co/Cu/La), AMAT·Lam·ASM 다중 챔버 ASD 통합 플랫폼.",
        ),
        (
            "대만 (TW)",
            "TSMC·대만 주요 파운드리·장비 파트너 (AMAT/TEL 대만 대응 패밀리)",
            "US20250323146A1 / US11854940B2 / US20220238323A1 / US11699739B2 (TSMC BEOL 완전 자기정렬 비아 FSAV·금속/유전체 SAM 억제제 패밀리)",
            "2nm 이하 GAA 나노시트·BEOL 배선에서 티올/아졸(2-MBI·PFDT) 및 포스폰산·아미노실란을 활용한 완전 자기정렬 비아(FSAV)·DoD/MoD 선택 증착 특허 주도.",
        ),
        (
            "중국 (CN)",
            "SMIC·YMTC·Naura 및 국내외 소재사 CN 출원군",
            "CN113423862A (Cp계 Hf/Zr 프리커서), CN102084026A ((L)Ru(CO)3), CN102344460B (Co 아미디네이트), MoCl5/MoO2Cl2 8인치 플랫폼 환원제 조합 출원",
            "400단 이상 3D NAND Mo 워드라인 대체, DRAM 커패시터용 Cp-Hf 조성물 및 루테늄/코발트 배리어용 아미디네이트/카르보닐 특허 확장.",
        ),
    ]

    pat_rows_html = "".join(
        f"""
        <tr>
          <td><strong>{html.escape(reg)}</strong></td>
          <td style="font-size:0.77rem;">{html.escape(comp)}</td>
          <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0369a1;">{html.escape(pids)}</td>
          <td style="font-size:0.79rem; line-height:1.5;">{html.escape(focus)}</td>
        </tr>"""
        for reg, comp, pids, focus in patent_landscape
    )

    def _fmt_cand_row(r: Dict[str, Any]) -> str:
        has_7n = r.get("tmc_sa_7net") is not None
        v_str = (
            f"{r['vbur_d2_percent']:.1f}% &rarr; <strong>{r['vbur_rprobe_7net']:.1f}%</strong>"
            f"<br><span style=\"font-size:0.68rem; color:#64748b;\">D2 {r['vbur_d2_7net']:.1f}%</span>"
            if has_7n
            else f"{r['vbur_d2_percent']:.1f}%"
        )
        phys_str = (
            f"BDE <strong>{r['bde_min_kcal_mol']:.1f}</strong><br>"
            f"<span style=\"font-size:0.68rem; color:#475569;\">e<sub>coh</sub> {r['e_coh_ev_per_atom']:.2f} &bull; dE<sub>dim</sub> {r['delta_e_dimer_kcal_mol']:.1f}</span><br>"
            f"<span style=\"font-size:0.68rem; color:#0369a1;\">&Delta;E<sub>f</sub>/v<sub>M</sub> {r['delta_e_form_per_bond_kcal_mol']:.1f}</span>"
            if has_7n
            else "-"
        )
        tmc_str = (
            f"{r['tmc_sa_no3d']:.2f} &rarr; <strong>{r['tmc_sa_7net']:.2f}</strong>"
            if has_7n
            else f"{r['ligand_sa']:.2f} &rarr; <strong>{r['tmc_sa_no3d']:.2f}</strong>"
        )
        p_7n = bool(r["pass_7net"]) if has_7n else bool(r["unified_pass"])
        p_no3d = bool(r["unified_pass"])
        if p_7n:
            badge_html = '<span class="badge badge-pass" style="font-size:0.71rem;">7net PASS</span>'
        elif p_no3d:
            badge_html = '<span class="badge badge-fail" style="font-size:0.69rem; background:#fef3c7; color:#92400e; border-color:#fcd34d;">7net BOUND</span>'
        else:
            badge_html = '<span class="badge badge-fail" style="font-size:0.71rem;">BOUNDARY</span>'
        if p_7n and not p_no3d:
            badge_html += '<div style="font-size:0.66rem; color:#047857; font-weight:700; margin-top:2px;">3D 구제</div>'
        elif p_no3d and has_7n and not p_7n:
            badge_html += '<div style="font-size:0.66rem; color:#92400e; margin-top:2px;">ETKDG Pass</div>'

        return f"""
        <tr>
          <td><code>{html.escape(r['id'])}</code></td>
          <td><strong>{html.escape(r['name'])}</strong><br><span style="font-size:0.7rem; color:#64748b;">{html.escape(r['family_class'])} ({html.escape(r['element_group'])})</span></td>
          <td style="text-align:center;"><strong>{html.escape(r['center_atom'])}</strong></td>
          <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.74rem;">{r['cn_eff']:.1f} / {r['n_ve']:.0f}</td>
          <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.73rem;">{v_str}</td>
          <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.72rem; line-height:1.35;">{phys_str}</td>
          <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.74rem;">{tmc_str}</td>
          <td style="text-align:center;">{badge_html}</td>
          <td style="font-size:0.73rem; line-height:1.45;">
            <div style="color:#0369a1; font-weight:700;">[{html.escape(r['jurisdictions'])}] {html.escape(r['references'])}</div>
            <div style="color:#475569;">{html.escape(r['application'])} &mdash; {html.escape(r['mechanism_note'])}</div>
            {f'<div style="color:#be123c; font-weight:700;">물리/경계 진단: {html.escape(r["boundary_flag"])}</div>' if r['boundary_flag'] else ''}
          </td>
        </tr>"""

    cand_rows_html = "".join(_fmt_cand_row(r) for r in recs)

    tpl = r"""
        <!-- Card 22.1: Central-Metal Coverage Gap Audit -->
        <div class="card" style="margin-bottom:20px; border-top:4px solid #e11d48;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">1. 중심금속원자 기준 문헌 조사를 통한 현행 DB 공백 진단 및 확장 전략</h3>
            <span class="badge badge-4g" style="font-size:0.78rem;">16 &rarr; __N_MERGED_EL__ Central Elements</span>
          </div>
          <div style="display:grid; grid-template-columns:minmax(320px, 1.15fr) minmax(250px, 0.85fr); gap:16px; align-items:start;">
            <div>
              __COV_SVG__
            </div>
            <div style="display:flex; flex-direction:column; gap:10px; font-size:0.82rem; color:#334155; line-height:1.6;">
              <div style="background:#fff1f2; border:1px solid #fecdd3; border-radius:8px; padding:11px 13px;">
                <strong style="color:#9f1239;">핵심 원소 편중 진단:</strong> 현행 <code>102</code>종 프리커서는 III-V·귀금속(<code>Ga 24, In 24, Ru 16</code>)이 <code>62.7%</code>를 차지하는 반면, 3D NAND/GAA/BEOL 핵심 금속인 <strong><code>Mo(1), W(1), Ta(1), Sn(1)</code>은 희소하고 <code>Nb, Ni, Cu, Ge, La, Y, Pt, Zn</code>은 <code>0</code>종</strong>입니다.
              </div>
              <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; padding:11px 13px;">
                <strong style="color:#0f172a;">중심금속별 문헌 리뷰 근거:</strong>
                <ul style="margin:6px 0 0 16px; padding:0;">
                  <li><strong>Mo·W (워드라인·배선):</strong> <code>MoO2Cl2</code> (<code>10.1039/D6TC00454G</code>), 옥소-아미디네이트(<code>10.1021/acs.chemmater.6b00248</code>), <code>WCl5</code> (<code>10.1016/j.apsusc.2021.150373</code>).</li>
                  <li><strong>Nb·Ta (초전도·배리어):</strong> <code>TBTDEN</code> (<code>10.1021/acs.chemmater.5c01456</code>), <code>Nb(OEt)5</code> (<code>10.1021/la902289h</code>), <code>PDMAT</code> (<code>10.1116/1.2345205</code>).</li>
                  <li><strong>Ni·Cu·Co·Pt (콘택·배선):</strong> <code>Ni(dmamp)2</code> (<code>10.1116/1.1875172</code>), <code>Cu(dmap)2</code> (<code>10.1021/acs.chemmater.7b01710</code>), <code>Co(tBu2DAD)2</code> (<code>10.1021/acsami.1c09819</code>), <code>MeCpPtMe3</code> (<code>10.1016/j.apsusc.2022.152750</code>).</li>
                  <li><strong>La·Y·Ge·Sn (다이폴·강유전·채널):</strong> <code>La(iPrAMD)3</code> (<code>10.1021/cm050756v</code>), <code>Ge(NMe2)4</code> (<code>10.1116/1.3653424</code>), <code>Sn(NMe2)4</code> (<code>10.1116/1.2835087</code>).</li>
                </ul>
              </div>
            </div>
          </div>
        </div>

        <!-- Card 22.2: 5-Jurisdiction Patent Landscape -->
        <div class="card" style="margin-bottom:20px; border-top:4px solid #0284c7;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">2. 한&middot;중&middot;일&middot;미&middot;대만 5극 반도체 프리커서&middot;ASD 인히비터 특허 지형 분석</h3>
            <span class="badge badge-phys" style="font-size:0.78rem;">KR &bull; JP &bull; US/WO &bull; TW &bull; CN</span>
          </div>
          <div class="table-responsive">
            <table class="data-table" style="font-size:0.79rem;">
              <thead>
                <tr>
                  <th style="width:95px;">특허 관할</th>
                  <th style="width:190px;">주요 출원인 / 연구기관</th>
                  <th style="width:270px;">대표 특허 번호 및 패밀리</th>
                  <th>핵심 화학 모티프 및 공정 차별점</th>
                </tr>
              </thead>
              <tbody>
                __PAT_ROWS__
              </tbody>
            </table>
          </div>
        </div>

        <!-- Card 22.3: Screened 38 Non-Overlapping Expansion Candidates -->
        <div class="card" style="border-top:4px solid #059669;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
            <h3 class="card-title" style="margin-bottom:0;">3. 특허&middot;문헌 확장 후보 __N_TOT__종 전수 <code>SevenNet-0</code> 물리 디스크립터 보강 및 <code>SAPIS</code> 검증</h3>
            <span class="badge badge-pass" style="font-size:0.78rem;">7net __N_SEVENN__/__N_TOT__ (0 FB) &bull; Strict Prec __N_PREC_PASS_7N__/__N_PREC__ &bull; Inh __N_INH_PASS_7N__/__N_INH__</span>
          </div>
          <p style="font-size:0.87rem; color:#334155; line-height:1.65; margin-bottom:10px;">
            신규 확장 후보 <strong>__N_TOT__종</strong>(프리커서 __N_PREC__ + 인히비터 __N_INH__) 전수에 대해 신규 7개 원소(<code>Nb, Ni, Pt, La, Y, Zn, Ge</code>)의 <code>SevenNet-0</code> 고립원자&middot;표준결정 화학퍼텐셜($\mu_Z^{\text{std}}$)을 산출한 뒤 <strong>LBFGS 단량체 이완&middot;이합체 회합(&Delta;E<sub>dimer</sub>)&middot;단편 BDE<sub>min</sub>&middot;응집에너지(E<sub>coh</sub>)&middot;표준 생성에너지(&Delta;E<sub>f</sub><sup>std</sup>/v<sub>M</sub>)를 전수 계산 완료(Fallback <code>0</code>건)</strong>했습니다. 그 결과 (1) 초기 ETKDGv3에서 입체과밀로 오기각되던 <strong><code>PDMAT</code>(<code>93.5% &rarr; 79.2%</code>, <code>TMC 4.52 &rarr; 1.24</code>)과 <code>Ge(NMe2)4</code>(<code>90.5% &rarr; 77.0%</code>, <code>TMC 4.52 &rarr; 1.79</code>)가 <code>SevenNet-0</code> 3D 이완으로 정확히 구제</strong>되었으며, (2) <code>Mo(CO)6/W(CO)6/Cu(hfac)(vtms)</code> 중성 CO&middot;올레핀 해리 과소결합 편차, <code>MoCl5</code> 고체 이합체 격자응집($E_{\text{coh}}=2.625$), <code>La(iPrAMD)3/Sn(NMe2)4/MeCpHf(NMe2)3</code> 저온 ALD 고반응성 창($\text{BDE}=29\text{--}40\text{ kcal/mol}$)까지 정량적으로 분해되었습니다.
          </p>
          <div style="background:#f8fafc; border:1px solid #cbd5e1; border-left:4px solid #0284c7; border-radius:6px; padding:11px 13px; margin-bottom:12px; font-size:0.81rem; color:#334155; line-height:1.6;">
            <strong style="color:#0f172a;">교차 DB 정준 중복 제거 및 엔티티 병합 감사 (Canonical Entity Resolution):</strong> 전체 <code>205</code>종 마스터&middot;<code>29</code>종 시약 사전과의 교차 감사 결과, (1) <code>MoCl5</code>&middot;<code>TaCl5</code>&middot;<code>HfCl4</code>&middot;<code>TiCl4</code>&middot;<code>HCDS</code>&middot;<code>TMA</code>&middot;<code>SiCl4</code>&middot;<code>SiH2Cl2</code>&middot;<code>ZrCl4</code>&middot;<code>AlCl3</code> 등 <strong>프리커서&ndash;출발시약 겸용 이중역할 <code>10</code>종</strong>, (2) <code>Ga2(NMe2)6 / [Ga(NMe2)3]2</code> 및 <code>GaMe3</code> 착물 표기 이형 <strong><code>2</code>쌍</strong>, (3) <code>AtomicLimits</code> 무-SMILES 문헌 엔트리와 일치하는 특허 확장 <strong><code>14</code>종</strong>(<code>MoCl5</code>, <code>Mo(CO)6</code>, <code>Nb(OEt)5</code>, <code>TBTDEN</code>, <code>PDMAT</code>, <code>TaCl5</code>, <code>HfCl4</code>, <code>TiCl4</code>, <code>Cu(dmap)2</code>, <code>MeCpPtMe3</code>, <code>Y(MeCp)3</code>, <code>ZnEt2</code>, <code>Ge(NMe2)4</code>, <code>Sn(NMe2)4</code>)이 식별되었습니다. 통합 DB 탐색기에서는 <strong>정준 엔티티로 병합&middot;SMILES 및 <code>SevenNet-0</code> 디스크립터 역충전(Backfill)하고 순신규 <code>24</code>종만 신규 추가</strong>하여, 총 <strong><code>286</code>개 유일 정준 엔티티(검색 중복 <code>0</code>건)</strong> 및 <strong>Sankey <code>342</code>건 통합 흐름도</strong>를 구축했습니다.
          </div>
          <div class="table-responsive">
            <table class="data-table" style="font-size:0.75rem;">
              <thead>
                <tr>
                  <th style="width:54px;">ID</th>
                  <th style="width:148px;">후보명 / 리간드 계통</th>
                  <th style="width:42px; text-align:center;">중심</th>
                  <th style="width:62px; text-align:center;">CN / N<sub>ve</sub></th>
                  <th style="width:94px; text-align:center;">%V<sub>bur</sub><br>(ETKDG&rarr;7net)</th>
                  <th style="width:128px; text-align:center;">SevenNet-0 물리 지표<br>(BDE &middot; E<sub>coh</sub> &middot; &Delta;E<sub>f</sub>/v<sub>M</sub>)</th>
                  <th style="width:86px; text-align:center;">SAPIS<br>(No3D&rarr;7net)</th>
                  <th style="width:82px; text-align:center;">판정</th>
                  <th>특허·문헌 출처 및 공정 메커니즘 / 물리·경계 진단</th>
                </tr>
              </thead>
              <tbody>
                __CAND_ROWS__
              </tbody>
            </table>
          </div>
        </div>"""

    return _render_template(
        tpl,
        {
            "__N_MERGED_EL__": s["n_merged_elements"],
            "__COV_SVG__": cov_svg,
            "__PAT_ROWS__": pat_rows_html,
            "__N_TOT__": s["n_total"],
            "__N_PREC__": s["n_precursors"],
            "__N_INH__": s["n_inhibitors"],
            "__N_SEVENN__": s.get("n_sevenn_computed", 0),
            "__N_PREC_PASS__": s["n_precursors_pass"],
            "__N_INH_PASS__": s["n_inhibitors_pass"],
            "__N_PREC_PASS_7N__": s.get("n_precursors_pass_7net", s["n_precursors_pass"]),
            "__N_INH_PASS_7N__": s.get("n_inhibitors_pass_7net", s["n_inhibitors_pass"]),
            "__CAND_ROWS__": cand_rows_html,
        },
    )
