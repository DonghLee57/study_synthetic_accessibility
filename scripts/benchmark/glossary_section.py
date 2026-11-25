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
        "상용 출발물질과 특허 합성 반응 데이터를 분리해 평가하는 3세대 하이브리드 합성 가능성 평가 모델입니다.",
        "eMolecules 2.31억 종 기반 BScore와 USPTO 특허 반응 약 100만 건 기반 RScore 중 더 유리한 최댓값을 채택합니다. "
        "특히 점수가 0 미만인 합성 병목 조각만 선별하여 평균을 산출하므로, 흔한 치환기가 많더라도 난해한 핵심 조각이 포함되어 있으면 합성 난이도 상승을 민감하게 포착합니다.",
    ),
    (
        "gl-cv5",
        "sa",
        "SA 모델·해싱",
        "CV-5",
        "Stratified 5-Fold Cross-Validation",
        "전체 데이터셋을 양성군과 대조군의 비율이 동일하게 유지되도록 5개 부분집합으로 나눈 뒤, 4개로 학습하고 나머지 1개로 검증하는 과정을 5회 반복하여 평균 성능을 산출하는 층화 교차검증 기법입니다.",
        "단일 분할 평가에서 특정 표본 구성에 따라 성능이 과대 추정되는 현상을 방지합니다. "
        "다만 학습 폴드와 검증 폴드가 동일한 실패 메커니즘 분포를 공유하므로, 학습 단계에서 경험하지 못한 새로운 화학적 붕괴 메커니즘에 대한 일반화 능력은 별도의 LOFMO 외삽 평가로 검증해야 합니다.",
    ),
    (
        "gl-lofmo",
        "sa",
        "SA 모델·해싱",
        "LOFMO",
        "Leave-One-Failure-Mechanism-Out Extrapolation",
        "합성 불가 대조군의 6대 물리·화학적 실패 메커니즘 중 특정 메커니즘 하나를 학습 과정에서 완전히 배제한 뒤, 처음 마주하는 해당 실패 메커니즘 분자들을 모델이 정확히 기각하는지 측정하는 외삽 검증 프로토콜입니다.",
        "무작위 교차검증에서는 높은 정확도를 보이던 통계적 기계학습 분류기라도 학습 데이터에 없던 새로운 실패 원리를 만나면 기각률이 급락하는 한계가 나타납니다. "
        "본 연구에서는 결합 해리 에너지와 입체 차폐율 등 물리 법칙에 기반한 SAPIS 판정 기준이 미지의 실패 메커니즘에서도 일관된 차단 성능을 유지함을 입증하는 핵심 기준으로 사용합니다.",
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
        "gl-else-node",
        "sa",
        "SA 모델·해싱",
        "단일 미지 노드 (ELSE)",
        "Unknown Atom Bucket / ELSE Node",
        "원자 종류를 고정된 길이의 벡터로 입력받는 그래프 신경망 모델에서 사전 정의된 주요 유기 원소 외의 모든 원소를 하나로 묶어 배정하는 기타 범주 노드를 뜻합니다.",
        "대표적으로 DFRscore 등 유기합성 기반 그래프 신경망은 탄소, 질소, 산소, 황, 인, 할로겐 등 약 10여 종의 유기 원소만 개별 항목으로 구분하고 나머지 원소는 모두 단일 ELSE 노드로 처리합니다. "
        "이 구조에 반도체 프리커서를 입력하면 하프늄, 지르코늄, 티타늄, 루테늄, 텅스텐처럼 산화수와 배위 성질이 전혀 다른 전이금속들이 모두 동일한 ELSE 노드로 치환되어 중심 금속 간 변별력이 완전히 상실됩니다.",
    ),
    (
        "gl-unleaked",
        "sa",
        "SA 모델·해싱",
        "데이터 누설 없는 확장",
        "Unleaked Domain Expansion (max(existing, -0.47))",
        "모델의 성능을 공정하게 검증하기 위해 평가 대상인 정답 프리커서의 조각 정보는 사전에 넣지 않고, 실제로 구매 가능한 출발 시약의 조각만 사전에 반영하는 원칙입니다.",
        "평가 대상 프리커서를 사전에 미리 등록하면 점수가 인위적으로 개선되는 데이터 누설이 발생하여 신규 후보 물질에 일반화할 수 없습니다. "
        "또한 시약 조각 점수를 0.0으로 덮어쓰면 BR-SAScore의 음수 조각 평균 계산에서 해당 조각이 제외되어 분모가 줄어드는 부작용이 발생하므로, 음수 영역 내 하한값인 -0.47로만 상향 보정합니다.",
    ),
    (
        "gl-logodds",
        "sa",
        "SA 모델·해싱",
        "로그오즈 및 로그오즈비",
        "Log-Odds & Bayesian Log-Odds Ratio",
        "어떤 사건이 일어날 확률을 일어나지 않을 확률로 나눈 승산 비율에 로그를 취한 값이며, 두 집단 간 출현 확률 비율에 로그를 취한 값을 로그오즈비라고 합니다.",
        "확률 p에 대해 오즈는 p / (1 - p)로 정의되고, 여기에 로그를 취한 로그오즈 log(p / (1 - p))는 확률 0.5를 기준으로 0이 되며 유불리에 따라 양수와 음수로 대칭 변환됩니다. "
        "1세대 베이지안 모델인 SYBA는 특정 부분구조 조각이 합성 용이 분자군에서 나타날 확률을 비합성 분자군에서 나타날 확률로 나눈 로그오즈비 log(P(f | 합성 용이) / P(f | 비합성))를 조각마다 미리 계산해 두고, 분자에 포함된 모든 조각의 로그오즈비를 단순 합산하여 양수이면 합성 가능, 음수이면 합성 난해로 판정합니다.",
    ),
    (
        "gl-rxn-corpus",
        "sa",
        "SA 모델·해싱",
        "반응 코퍼스",
        "Chemical Reaction Corpus",
        "특허 문헌과 학술 논문에서 수집한 수많은 화학 반응의 출발물질, 시약, 생성물 구조 변화를 기계가 학습할 수 있도록 정제해 구축한 대규모 화학 반응 기록 집합입니다.",
        "언어학에서 대규모 말뭉치를 뜻하던 개념을 화학 정보학에 도입한 용어로, 미국 특허청 반응 데이터베이스인 USPTO와 상용 문헌 데이터베이스인 Reaxys 반응 집합이 대표적입니다. "
        "2세대 이후의 합성 가능성 평가 모델과 역합성 탐색 엔진은 정적인 단일 분자 목록 대신 이 반응 코퍼스를 학습하여 실제 결합 형성 규칙과 단계별 구조 복잡도 변화를 반영합니다.",
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
          <input type="text" id="glossarySearchInput" class="db-search-input" placeholder="용어 검색 (예: 로그오즈, ELSE, 반응 코퍼스, LOFMO, CV-5, OOV, SMILES...)" oninput="filterGlossary()">
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
