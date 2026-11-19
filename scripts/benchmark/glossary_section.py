"""Glossary Tab Builder (``용어 사전``) for Semiconductor SA & Precursor Modeling.

Generates the standalone HTML content for the 5th top-level view tab
(``#tab-content-glossary``, button label ``용어 사전``) covering foundational
terminology across SA models, SMILES/cheminformatics, coordination chemistry,
and semiconductor ALD/CVD precursor physics.
"""

from __future__ import annotations

import html
from typing import Dict, List, Tuple

# (id, category_key, category_label, term, full_name, summary, detail)
# Ordered alphabetically (symbols/A-Z first, then Korean ㄱ-ㅎ).
GLOSSARY_ENTRIES: List[Tuple[str, str, str, str, str, str, str]] = [
    (
        "gl-phys",
        "phys",
        "반도체 공정·물리 지표",
        "%V_bur / BDE_min / E_coh",
        "Buried Volume / Minimum Bond Dissociation Energy / Cohesive Energy",
        "금속 중심 주변의 입체적 보호 정도(%V_bur), 가장 약한 금속–리간드 결합 세기(BDE_min), 분자 간 응집 에너지(E_coh)를 뜻하는 3D 물리 지표입니다.",
        "유기 SA 모델이 구별하지 못하는 금속 착체의 열분해 안정성, 이합체화(고분자화) 위험, 휘발성을 정량 검증하기 위해 SevenNet-0 구조 최적화로부터 직접 계산합니다.",
    ),
    (
        "gl-addhs-chiral",
        "chem",
        "배위화학·SMILES",
        "AddHs 및 카이랄 (Chiral) 중심",
        "Explicit Hydrogen Addition & Chirality Center",
        "AddHs는 생략되어 있던 수소 원자(H)를 분자 그래프에 명시적 노드로 추가하는 연산이며, 카이랄 중심은 거울상 이성질체를 만드는 비대칭 원자입니다.",
        "SAScore는 수소가 생략된 heavy-atom 그래프를 전제로 설계되었기 때문에 AddHs를 제거하지 않으면 원자 수(n_At)와 지문 조각이 늘어나 점수가 인위적으로 급등합니다. "
        "카이랄 중심 개수는 복잡도 페널티(Eq. 2)의 log10(n_Chiral + 1) 항에 직접 가산됩니다.",
    ),
    (
        "gl-process",
        "phys",
        "반도체 공정·물리 지표",
        "ALD / CVD / AS-ALD / Negative Control",
        "Atomic Layer / Chemical Vapor Deposition & Negative Control",
        "반도체 초미세 박막 증착 공정(ALD/CVD), 특정 표면만 선택적으로 막는 차단제(AS-ALD 인히비터), 그리고 모델 오탐 검증용 불안정 대조군(Negative Control)을 뜻합니다.",
        "Negative Control은 TiEt4, ZrMe4처럼 실제 공정에서 불안정하거나 합성이 불가능한 40종 분자로 구성되어, SA 모델이 이를 무분별하게 통과(False Positive)시키지 않는지 검증하는 기준선 역할을 합니다.",
    ),
    (
        "gl-brsascore",
        "sa",
        "SA 모델·해싱",
        "BR-SAScore (BScore / RScore)",
        "Building-Block & Reaction Synthetic Accessibility Score",
        "상용 출발물질(Building Block, B)과 특허 합성 반응(Reaction, R) 데이터를 분리해 평가하는 3세대 하이브리드 SA 모델입니다.",
        "eMolecules 2.31억 종 기반 BScore와 USPTO 특허 반응 ~100만 건 기반 RScore 중 더 유리한 값(max)을 채택합니다. "
        "특히 0 미만인 병목 조각만 평균(Eq. 12)내어, 흔한 치환기가 많아도 어려운 핵심 조각이 있으면 합성 난해로 민감하게 포착합니다.",
    ),
    (
        "gl-morgan",
        "sa",
        "SA 모델·해싱",
        "Morgan Fingerprint",
        "Circular Morgan Substructure Fingerprint (r = 0~2)",
        "분자 내 각 원자를 중심으로 결합 거리 반경 r = 0, 1, 2 이내 이웃 원자의 화학적 환경을 부호화해 고유 정수 해시 ID(bitId)로 변환하는 원형 분자 지문입니다.",
        "Radius 0은 중심 원자 자체, Radius 1은 직접 결합한 1차 이웃, Radius 2는 2결합 떨어진 고리·치환기 맥락을 포착합니다. "
        "SAScore와 BR-SAScore는 이 bitId가 학습 사전에서 얼마나 자주 등장했는지를 조회(lookup)하여 합성 난이도를 계산합니다.",
    ),
    (
        "gl-oov",
        "sa",
        "SA 모델·해싱",
        "OOV (사전 미등재)",
        "Out-Of-Vocabulary",
        "학습 데이터셋에 한 번도 등장하지 않아 모델의 기준 사전(어휘집)에 등록되어 있지 않은 분자 조각이나 원소를 뜻합니다.",
        "자연어 처리에서 어휘집에 없는 단어를 의미하던 용어가 화학 AI 모델로 확장되었습니다. "
        "SAScore (2009)와 BR-SAScore (2024)는 유기물 데이터셋 기반이므로 금속–리간드 결합 조각이 OOV로 분류되어 각각 바닥값(-4.0, -6.0)을 받습니다. "
        "GASA 같은 그래프 신경망에서는 미지원 원소 입력 시 에러를 내지 않고 특성 벡터를 전부 0(Zero-Vector)으로 처리하는 문제가 발생합니다.",
    ),
    (
        "gl-smiles",
        "chem",
        "배위화학·SMILES",
        "SMILES",
        "Simplified Molecular-Input Line-Entry System",
        "2차원 분자 구조 그래프를 알파벳·숫자·기호로 이루어진 1줄 텍스트 문자열로 변환해 컴퓨터가 읽을 수 있게 만든 표준 표기법입니다.",
        "원자 기호(C, N, O, Hf), 단일/이중 결합(=, #), 고리 연결 번호(1, 2), 가지치기 괄호(()), 배위 결합 화살표(->) 등을 조합해 분자를 표현합니다. "
        "동일한 화학 물질이라도 고리 닫힘이나 배위 결합 표기 방식에 따라 여러 SMILES 문자열이 만들어질 수 있어, 그래프 기반 SA 점수가 달라지는 원인이 됩니다.",
    ),
    (
        "gl-tmcsa",
        "phys",
        "반도체 공정·물리 지표",
        "SAPIS & ligand_sa",
        "Synthetic Accessibility of Precursors & Inhibitors for Semiconductors",
        "본 연구에서 제안한 반도체 프리커서 및 AS-ALD 인히비터 통합 2트랙 평가 체계로, 분리된 프로리간드/유기물의 합성성(ligand_sa)과 3D 배위·열역학 안정성을 결합한 지표입니다.",
        "Track 1(ligand_sa)은 프로리간드별 min(BR-SAScore, SAScore) 중 최댓값으로 리간드 합성 난이도를 측정하고, "
        "Track 2는 SevenNet-0 기계학습 포텐셜로 계산한 입체 차폐율(%V_bur), 결합 해리 에너지(BDE_min), 응집 에너지(E_coh) 게이트를 결합해 최종 합성 가능성을 판정합니다.",
    ),
    (
        "gl-hapticity",
        "chem",
        "배위화학·SMILES",
        "η¹ / η⁵ 배위 (Hapticity)",
        "Hapticity (η^n, 햅티시티)",
        "그리스 문자 에타(η)와 위첨자 n은 하나의 리간드 안에서 연속된 n개의 원자가 동시에 중심 금속에 결합하고 있음을 나타냅니다.",
        "시클로펜타디에닐(Cp) 리간드는 고리 내 탄소 1개만 금속과 시그마 결합을 이루는 η¹ 방식과, 5원 방향족 고리의 π 전자 구름 전체가 면으로 배위하는 η⁵ 방식으로 나뉩니다. "
        "고배위(η⁵)로 그래프를 재구성할수록 3중 이상 다중 고리 복잡도가 추가되어 바닐라 SAScore 점수가 급격히 상승합니다.",
    ),
    (
        "gl-unleaked",
        "sa",
        "SA 모델·해싱",
        "데이터 누설 없는 확장",
        "Unleaked Domain Expansion (max(existing, -0.47))",
        "모델의 성능을 공정하게 검증하기 위해 평가 대상(정답 프리커서)의 조각 정보는 사전에 절대 넣지 않고, 실제로 구매 가능한 출발 시약의 조각만 사전에 반영하는 원칙입니다.",
        "만약 평가할 프리커서 데이터를 사전에 미리 학습시키면 점수가 인위적으로 좋아지는 데이터 누설(Data Leakage)이 발생합니다. "
        "또한 시약 조각을 0.0으로 덮어쓰면 BR-SAScore의 음수 필터(Score < 0)에서 탈락해 분모가 줄어드는 부작용이 생기므로, 음수 영역 하한(-0.47)으로만 상향 보정합니다.",
    ),
    (
        "gl-chelate",
        "chem",
        "배위화학·SMILES",
        "킬레이트 (Chelate) & 배위 결합 (->)",
        "Chelate Ring / Dative Coordination Bond",
        "하나의 리간드가 두 개 이상의 비공유 전자쌍 주개 원자로 중심 금속에 동시에 결합해 고리 구조를 형성하는 배위 형태입니다.",
        "대표 예로 아세틸아세토네이트(acac)가 산소 2개로 Hf에 결합해 6원 고리를 만드는 구조가 있습니다. "
        "SMILES에서 이를 단좌 개방형으로 둘지 dative bond(->)로 고리를 닫을지에 따라 고리 수와 Morgan Radius 1~2 지문이 바뀌어 SAScore가 크게 달라집니다.",
    ),
    (
        "gl-proligand",
        "chem",
        "배위화학·SMILES",
        "프로리간드 (Pro-ligand) 및 분리",
        "Neutral Protonated Pro-ligand (HL_i) & Decoupling",
        "금속에 배위된 음이온 리간드에 수소(H)를 붙여 중성 유기 분자 형태로 되돌린 리간드 원형(예: -NMe2 -> Me2NH, -Cp -> HCp)입니다.",
        "실제 공정용 프리커서 합성에서는 금속 할라이드/알킬 출발물질과 이 중성 프로리간드를 반응시켜 프리커서를 만듭니다. "
        "따라서 금속 중심을 제거하고 프로리간드만 분리해 평가(compute_ligand_sa)하면 금속 OOV 붕괴 없이 실제 유기 리간드 합성 난이도를 정확히 측정할 수 있습니다.",
    ),
]


def build_glossary_html() -> str:
    """Build the searchable ``용어 사전`` tab HTML fragment."""
    cards_html: List[str] = []
    for gid, cat_key, cat_lbl, term, full_name, summary, detail in GLOSSARY_ENTRIES:
        badge_style = (
            "background:#eff6ff; color:#1d4ed8; border:1px solid #bfdbfe;"
            if cat_key == "sa"
            else (
                "background:#f0fdfa; color:#0f766e; border:1px solid #99f6e4;"
                if cat_key == "chem"
                else "background:#fef3c7; color:#92400e; border:1px solid #fde68a;"
            )
        )
        search_blob = html.escape(f"{term} {full_name} {cat_lbl} {summary} {detail}".lower())
        cards_html.append(
            f"""
          <article id="{gid}" class="card glossary-card" data-gcat="{cat_key}" data-gtext="{search_blob}" style="margin:0; display:flex; flex-direction:column; gap:8px;">
            <div style="display:flex; justify-content:space-between; align-items:flex-start; gap:8px; flex-wrap:wrap;">
              <div>
                <h3 style="margin:0; font-size:1.02rem; font-weight:800; color:#0f172a;">{html.escape(term)}</h3>
                <div style="font-size:0.78rem; font-family:'JetBrains Mono', monospace; color:#0369a1; margin-top:2px;">{html.escape(full_name)}</div>
              </div>
              <span style="font-size:0.72rem; font-weight:700; padding:2px 8px; border-radius:999px; {badge_style}">{html.escape(cat_lbl)}</span>
            </div>
            <p style="margin:0; font-size:0.86rem; font-weight:600; color:#1e293b; line-height:1.6;">{html.escape(summary)}</p>
            <p style="margin:0; font-size:0.82rem; color:#475569; line-height:1.65;">{html.escape(detail)}</p>
          </article>"""
        )

    return f"""
      <div class="db-controls-card" style="margin-bottom:16px;">
        <div class="db-search-row" style="margin-bottom:10px;">
          <input type="text" id="glossarySearchInput" class="db-search-input" placeholder="용어 검색 (예: OOV, SMILES, 킬레이트, eta, 프로리간드, 누설, 카이랄...)" oninput="filterGlossary()">
          <button type="button" class="btn-toc-utility" style="width:auto; padding:9px 14px;" onclick="resetGlossaryFilter()">초기화</button>
          <span id="glossaryCountBadge" style="font-family:'JetBrains Mono', monospace; font-weight:700; color:#0369a1; font-size:0.84rem;">총 {len(GLOSSARY_ENTRIES)}개 용어</span>
        </div>
        <div class="db-pill-bar" id="glossaryCatPills" style="margin-bottom:0;">
          <button type="button" class="db-pill active" data-gcat="all" onclick="setGlossaryCategory('all', this)">전체 ({len(GLOSSARY_ENTRIES)})</button>
          <button type="button" class="db-pill" data-gcat="sa" onclick="setGlossaryCategory('sa', this)">SA 모델·해싱</button>
          <button type="button" class="db-pill" data-gcat="chem" onclick="setGlossaryCategory('chem', this)">배위화학·SMILES</button>
          <button type="button" class="db-pill" data-gcat="phys" onclick="setGlossaryCategory('phys', this)">반도체 공정·물리 지표</button>
        </div>
      </div>

      <div id="glossaryEmptyState" class="card" style="display:none; text-align:center; padding:32px; color:#64748b;">
        일치하는 용어가 없습니다. 검색어나 분류 필터를 확인해 주십시오.
      </div>

      <div id="glossaryCardGrid" class="grid-2" style="gap:14px;">
        {''.join(cards_html)}
      </div>
    """
