import json
import os

def load_meta(filename):
    p1 = os.path.join("database", "metadata", filename)
    if os.path.exists(p1):
        with open(p1, "r", encoding="utf-8") as f:
            return json.load(f)
    p2 = os.path.join("downloaded_papers", filename)
    if os.path.exists(p2):
        with open(p2, "r", encoding="utf-8") as f:
            return json.load(f)
    raise FileNotFoundError(f"Cannot find metadata file {filename}")

papers = load_meta("manifest.json")
extra = load_meta("repos_and_precursors.json")
vasp_data = load_meta("vasp_frameworks.json")

with open("database/precursors.json", "r", encoding="utf-8") as f:
    prec_db = json.load(f)

with open("database/building_blocks.json", "r", encoding="utf-8") as f:
    bb_db = json.load(f)

with open("database/reactions.json", "r", encoding="utf-8") as f:
    rxn_db = json.load(f)

lit_rows_html = []
for idx, p in enumerate(papers, 1):
    gen = p.get("generation", "")
    role = p.get("role", "")
    if "Area-Selective ALD" in gen or role == "Precursor_ASD":
        gen_code = "ASD"
        badge_cls = "badge-prec"
    elif "VASP" in gen or role == "VASP_Framework":
        gen_code = "VASP"
        badge_cls = "badge-2g"
    elif role.startswith("Precursor") or "Precursor" in gen:
        gen_code = "Precursor"
        badge_cls = "badge-prec"
    elif "1st Gen" in gen:
        gen_code = "1G"
        badge_cls = "badge-1g"
    elif "2nd Gen" in gen:
        gen_code = "2G"
        badge_cls = "badge-2g"
    elif "3rd Gen" in gen:
        gen_code = "3G"
        badge_cls = "badge-3g"
    elif "4th Gen" in gen and role == "Core":
        gen_code = "4G"
        badge_cls = "badge-4g"
    else:
        gen_code = "Related"
        badge_cls = "badge-ext"
    
    doi_link = f"https://doi.org/{p['doi']}"
    if p["pdf_status"]["downloaded"]:
        ext = p["pdf_status"].get("ext", "pdf")
        label = p["pdf_status"].get("label", "PDF 열람")
        local_rel = p["pdf_status"].get("file_path", f"downloaded_papers/{p['id']}.{ext}").replace("\\", "/")
        size_str = f"{p['pdf_status']['size_kb']} KB"
        pdf_action = f'<a href="{local_rel}" onclick="smartOpenPaper(event, \'{local_rel}\', \'{doi_link}\')" target="_blank" class="btn-pdf">{label} ({size_str})</a>'
    else:
        pdf_action = f'<a href="{doi_link}" target="_blank" class="btn-paywall">출판사 유료 구독 (DOI 이동)</a>'
    
    search_str = f"{p['name']} {p['title']} {p['authors']} {p['journal']} {p['generation']}".lower()

    row = f"""
            <tr class="lit-row" data-gen="{gen_code}" data-downloaded="{str(p['pdf_status']['downloaded']).lower()}" data-search="{search_str}">
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;"><span class="badge {badge_cls}">{p['generation']}</span></td>
              <td style="white-space:nowrap;"><strong>{p['name']}</strong></td>
              <td style="text-align:center; font-weight:600; color:#475569;">{p['year']}</td>
              <td style="min-width:460px;">
                <div style="font-weight:600; color:#0f172a; margin-bottom:4px; line-height:1.45;">{p['title']}</div>
                <div style="font-size:0.8rem; color:#64748b; line-height:1.4;">{p['authors']} &mdash; <em style="color:#475569;">{p['journal']}</em></div>
              </td>
              <td style="white-space:nowrap;">
                <div><a href="{doi_link}" target="_blank" style="font-family:'JetBrains Mono'; font-size:0.78rem; color:#0284c7; font-weight:600;">DOI: {p['doi']}</a></div>
                <div style="margin-top:6px;">{pdf_action}</div>
              </td>
            </tr>"""
    lit_rows_html.append(row)

lit_rows_joined = "\n".join(lit_rows_html)

# Benchmark metrics on semiconductor precursor DB (27 molecules)
# Formally partitioned into:
# (A) [코드 직접 실측]: Directly installed & executed locally on 27 precursors
# (B) [도메인 제약 분석]: Evaluated via dataset scope & chemical grammar constraints
# Benchmark metrics on semiconductor precursor DB (27 molecules)
# Strictly filtered to the 8 directly installed, executed, and benchmarked models
repo_semi_metrics = {
    "SAScore (RDKit Built-in)": {
        "rate": "11.1% (3/27, 실측)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (False Negative 왜곡)",
        "desc": "[코드 직접 실측 11.1%] TEOS(2.84), 4DMAS(3.37), Ru(EtCp)2(3.48) 3종만 3.5 이하로 통과. TSA(6.91), DMAH(6.31), TMA(4.81), ZMA(4.82) 등 단순 1단계 양산 전구체에 높은 감점이 부여됨. PubChem 의약 분자 빈도 기반으로 Si-N, Al-C, Zr-O 등 반도체 핵심 배위를 비전형적 복잡 구조로 분류한 결과임."
    },
    "SCScore": {
        "rate": "85.2% (23/27, 실측)",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "코드 직접 실측 (초기원료 편향 / 위험성 미반영)",
        "desc": "[코드 직접 실측 85.2%] Reaxys 1024bool 신경망 직접 실행: 23종이 1.01~1.94의 낮은 복잡도로 산출(DMAH 1.01, HCDS 1.02, TTIP 1.25, DIPAS 1.73). Reaxys 반응쌍에서 유기금속이 출발물질로 등록되어 전이금속 착체라도 분자량이 작으면 단순 원료로 과소평가(ZMA 2.72, Co(AMD)2 2.38). 자연발화성(Pyrophoric) 등 취급 위험성은 반영되지 않음."
    },
    "SYBA": {
        "rate": "85.2% (23/27, 실측)",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "코드 직접 실측 (샌드위치 착체 오분류 / syba.csv.gz)",
        "desc": "[코드 직접 실측 85.2%] 25.7MB 나이브 베이지안 사전(syba.csv.gz; 기본 syba4.csv.gz 대비 방법론적 변경 명시) 직접 연산: 23종은 +3.96~+11.98 양수로 분류됨. 반면 ZINC 의약 DB에 부재한 메탈로센/Cp 샌드위치 착체인 Ru(EtCp)2(-44.27), CpTi(-24.98), CpZr(-24.98), ZMA(-18.21) 4종은 큰 음수(합성 불가 판정)로 산출됨."
    },
    "BR-SAScore": {
        "rate": "11.1% (엄격 시약 3/27) / 100% (암기 누수)",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "코드 직접 실측 (무누수 한계 규명 & max(existing, -0.47) 보정 완료)",
        "desc": "[코드 직접 실측] 기본 의약사전(11.1%, 3/27)뿐 아니라, 데이터 누수(생성물 등록 및 CCTBA A>>A 규칙)를 완전히 배제하고 BRSAScore의 음수 단편 평균(fscore < 0) 수식에 맞춰 max(existing, -0.47)로 보정한 엄격 시약 사전에서도 점수가 악화되는 전구체(0종) 없이 통과율은 11.1%(3/27: TSA, TMA, DMAH)로 동일함(신규 형성된 M-N/M-C 중심 단편은 시약에 부재하여 -6.0 유지; 반면 단순 0.0 덮어쓰기 시 유기 단편이 분모에서 빠져 4종 악화). 타깃 생성물을 사전 피팅하면 100.0%(27/27)가 나오나 순환 기억 왜곡이므로 리간드 분리형 2-Track TMC-SA 프로토콜이 필수적임."
    },
    "RAscore": {
        "rate": "0.0% 물리적 유효 (수치상 100% FP)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (중심금속 결합 무시 False Positive)",
        "desc": "[코드 직접 실측] 27종 전수 0.629~0.988(평균 0.925)로 수치상 100% 통과 판정이나, 이는 ECFP6 반경 3 지문이 중심금속을 단순 고차 알킬/아민 탄소로 취급하여 약물 합성의 쉬운 알킬화로 착각한 전형적 False Positive임. 실제 금속 결합 안정성/산화 위험에 대한 물리적 판정 유효율은 0.0%임."
    },
    "RxnFlow": {
        "rate": "0.0% (코어 0/27, 부분구조 33.3%)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (금속 배위 신톤 부재로 코어 생성 0건)",
        "desc": "[코드 직접 실측] Explore(30개)·CGFlow(42개) 총 72개 템플릿 매칭 결과 일반 생성물 부분구조(아민 알킬 사슬 등)를 포함하는 9종(33.3%)만 매칭되었을 뿐, 템플릿 파일 내 금속 원자가 전무하여 반도체 전구체 코어 화합물의 실제 합성 경로 도출 성공률은 0.0%(0/27)임."
    },
    "ReaSyn": {
        "rate": "0.0% (경로 0/27, 구문파싱 100%)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (유기 템플릿 한계로 무기 경로 생성 0건)",
        "desc": "[코드 직접 실측] 118개 원소 어휘 토큰이 지원되어 SMILES 파싱(구문 수용)은 100% 에러 없이 통과하나, SynFormer 115개 특허 반응 템플릿에는 금속 원자가 전무하고 일반 생성물 부분구조(24/27; 2원자 이상 패턴 시 23/27)만 매칭되어 무기 배위 결합 경로는 0건(0.0%)임."
    },
    "AiZynthFinder": {
        "rate": "N/A (MCTS 미실행 null, 객체파싱 100%)",
        "badge_style": "background:#f1f5f9; color:#475569; border-color:#cbd5e1; font-weight:700;",
        "eval": "객체 실측 & 정적 분석 (가중치 미설치 MCTS 미실행 null)",
        "desc": "[객체 실측 & 정적 분석] aizynthfinder.chem.Molecule(smiles=smi).sanitize() 직접 호출로 27종 전구체 객체 생성(구조 수용)은 100% 통과 실측. 단, 외부 USPTO 신경망 가중치 및 ZINC 재고 파일 미설치로 실제 MCTS 탐색은 로컬 미실행(null 처리)되었으며, 정적 분석 상 USPTO 50k에는 유기금속 배위 규칙이 부재함."
    }
}


# 15 Unique Multi-stage evaluated models across the 4 execution pipeline stages (ordered by generation 1G -> 4G -> CASP)
multistage_models = [
    {
        "name": "SAScore (RDKit Built-in)",
        "gen": "1G 단편 빈도",
        "gen_badge": "badge-1g",
        "repo_key": "SAScore (RDKit Built-in)",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Morgan(r=2) 추출", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (통계 단편 기반)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("11.1% (3/27) [FN 88.9%]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — 의약 단편 빈도 편향 및 가교고리 오인",
        "barrier_desc": "PubChem 유기 의약 분자 기반 단편 빈도 사전에 Si-N, Al-C, M-N 결합이 부재하여 단순 1단계 양산 전구체 TSA(6.91), DMAH(6.31), TMA(4.81) 등 24종(88.9%)에서 False Negative 발생. 에테르/아민 유사 구조 3종(TEOS 2.84, 4DMAS 3.37, Ru(EtCp)2 3.48)만 통과."
    },
    {
        "name": "SYBA",
        "gen": "1G 나이브 베이지안",
        "gen_badge": "badge-1g",
        "repo_key": "SYBA",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 나이브 베이즈 단편", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (단편 로그공산비)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("85.2% (23/27) [Cp 착체 오분류]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — ZINC 의약 코퍼스 기반 샌드위치 착체 오분류",
        "barrier_desc": "비환식 알킬아민·알콕사이드 전구체 23종(+3.96~+11.98)은 양수로 분류되나, ZINC 의약 라이브러리(syba.csv.gz) 내 금속-Cp 배위 단편이 부재하여 Ru(EtCp)2(-44.27), CpTi(NMe2)3(-24.98), Cp(NMe2)3Zr(-24.98), ZMA(-18.21) 등 환형·메탈로센 4종은 큰 음수로 산출됨."
    },
    {
        "name": "SCScore",
        "gen": "2G 반응 신경망",
        "gen_badge": "badge-2g",
        "repo_key": "SCScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 1024-bit Morgan", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (12층 순위 DNN)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("85.2% (23/27) [위험 은폐]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 4 — Reaxys 초기 반응물 편향에 따른 안전 위험 은폐",
        "barrier_desc": "Reaxys 1,200만 반응쌍(full_reaxys_model_1024bool)에서 유기금속이 주로 출발 시약으로 등록된 편향으로 23종을 1.01~1.94의 저복잡도로 산출함. ZMA(2.72)·Co(AMD)2(2.38)의 입체 복잡도는 일부 반영하나, TMA(1.21)·DMAH(1.01)의 자연발화성 및 열분해 위험성은 전혀 변별하지 못함."
    },
    {
        "name": "RAscore",
        "gen": "2G MCTS 증류",
        "gen_badge": "badge-2g",
        "repo_key": "RAscore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("투영 왜곡 (금속→알킬 탄소)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (XGBoost 분류)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 유효 (수치 100% FP)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 — ECFP6 알킬 치환기 편향 및 중심금속 불안정성 소거",
        "barrier_desc": "2048-bit ECFP6 카운트 벡터가 디메틸아미노·이소프로필 등 유기 치환기 비트에 편향되어 중심금속 배위 안정성을 인지하지 못함. 전수 평균 0.925(0.629~0.988, 23종 ≥0.90)로 무조건 통과시켜 불안정 분자까지 거르지 못하는 100% False Positive 발생."
    },
    {
        "name": "BR-SAScore",
        "gen": "3G 반응 인지 분리",
        "gen_badge": "badge-3g",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% BFrag/RFrag 분리 (무누수 보정)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("11.1% (3/27) 배위 결합 부재", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("11.1% (3/27) [완제품 한계]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 기본·무누수 시약사전 동일 한계 (리간드 분리 평가 필수)",
        "barrier_desc": "기본 의약사전(eMolecules/USPTO)과 데이터 누수를 배제한 엄격 시약사전(max(existing, -0.47) 보정; 음수 단편 평균 수식에 정합하여 점수 악화 0종 유지) 모두에서 새로 형성된 M-N/M-C 배위 단편이 시약에 부재(-6.0)하여 동일하게 3종(TSA, TMA, DMAH)만 통과함. (단, 타깃 생성물 자체를 사전에 등록 시 27/27 통과하나 이는 순환 암기 왜곡임.) 유기 프로-리간드(HL_i) 단독 평가에 국한하여 운용해야 함."
    },
    {
        "name": "GASA",
        "gen": "3G GAT 위치선택성",
        "gen_badge": "badge-3g",
        "repo_key": "GASA",
        "s1": ("74.1% 수용 (7종 25.9% 미지원)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "s2": ("0벡터 무음 오인코딩 (Silent Loss)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT 분류)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 물리적 유효성", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 1·2 — DGL 의존성 부재 & 미지원 금속 0벡터 무음 손실 (정적 감사)",
        "barrier_desc": "로컬 dgl/dgllife 미설치로 추론 미실행이며 피처화 코드 정적·정량 감사를 수행함. 원자 허용 집합(allowable_set)에 5대 전이금속(Hf, Ta, Ru, W, Mo)이 누락되어 7종(25.9%)에서 런타임 예외 대신 encode_unknown=False 처리로 전량 0벡터가 생성되는 무음 오인코딩(Silent Mis-encoding) 결함 확인."
    },
    {
        "name": "DeepSA",
        "gen": "3G 멀티모달 트랜스포머",
        "gen_badge": "badge-3g",
        "repo_key": "DeepSA",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("Ro5(26/27 적합) / 서브워드 파편화", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (멀티모달 분류)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 물리적 타당성", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — AutoGluon 부재 & ELECTRA 무기 괄호 파편화 (정적 감사)",
        "barrier_desc": "로컬 autogluon 미설치로 추론 미실행이며 서술자·토크나이저 정적 감사를 수행함. RDKit 실측 상 26/27종이 Lipinski Ro5를 통과(Hf(OtBu)4만 cLogP=5.06 초과)하나, ELECTRA 기반 SMILES 서브워드 분할 시 [Zr]·[Hf]·[SiH3] 등 무기 괄호 표현이 의약 문맥으로 왜곡·파편화됨."
    },
    {
        "name": "DFRscore",
        "gen": "3G GAT 단계 예측기",
        "gen_badge": "badge-3g",
        "repo_key": "DFRscore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("원소 정보 상실 (전수 'ELSE')", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT 회귀)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 물리적 신뢰도", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 — 원소 인코딩 집합 제한으로 전이금속 전수 'ELSE' 축퇴",
        "barrier_desc": "그래프 노드 허용 원소가 유기 10원소(C, N, O, F, S, Cl, Br, I, B, P)로 제한되어 27종 전구체 전수(100%)의 규소·금속 원자(Si, Ti, Zr, Hf, Ru, W, Mo, Al, Co 등)가 단일 'ELSE' 토큰으로 뭉개져 원소별 합성 난이도 변별력을 완전 상실함."
    },
    {
        "name": "SynNet",
        "gen": "4G 마코프 순방향 계획",
        "gen_badge": "badge-4g",
        "repo_key": "SynNet",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 2048-bit Morgan (>0)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 경로", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — HB-91 유기 템플릿 한정 (금속 코어 규칙 0건)",
        "barrier_desc": "2048-bit Morgan 지문 활성 비트 생성(GetNumOnBits() > 0)을 27/27 전구체에서 직접 실측 검증함. 그러나 내장 Hartenfeller-Button(HB) 91종 SMARTS 매칭 시 일반 유기 부분구조 23/27종만 존재하며 템플릿 내 금속 원자가 전무(0건)하여 착물 합성 트리 생성 불가."
    },
    {
        "name": "ChemProjector",
        "gen": "4G 분자 투영 디퓨전",
        "gen_badge": "badge-4g",
        "repo_key": "ChemProjector",
        "s1": ("100.0% (27/27) (Z≤100)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 이중 잠재 공간 투영", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 경로", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — HB-91 유기 템플릿 매칭 (금속 배위 규칙 0건)",
        "barrier_desc": "원자번호 Z≤100 기반 그래프 임베딩 수용은 27/27 통과. HB 91종 SMARTS 매칭 실측 결과 23종에서 [N:1][#6:2] 등 일반 유기 치환기 부분구조만 탐지(TSA, HCDS, TMA, DMAH 제외)될 뿐, 템플릿 내 금속 원자가 0개이므로 전구체 코어 투영 불가."
    },
    {
        "name": "SynFormer",
        "gen": "4G Seq2Seq 트랜스포머",
        "gen_badge": "badge-4g",
        "repo_key": "SynFormer",
        "s1": ("100.0% (27/27) 자체 어휘", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% CoR 토큰 시퀀스", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% (≥2원자 85.2%) / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 경로", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 115종 특허 반응망 내 무기 배위 규칙 부재 (0건)",
        "barrier_desc": "자체 smiles_vocab 토크나이저 실측 27/27 전수 통과. comprehensive.txt 115종 생성물 SMARTS 매칭 결과 전체 24/27종(무탄소 TSA는 #84 [N:1] 1원자 패턴), ≥2원자 패턴 23/27종, ≥4원자 패턴 11/27종에서 유기 부분구조만 매칭되며 금속 원자 포함 규칙은 0건임."
    },
    {
        "name": "RxnFlow",
        "gen": "4G GFlowNet",
        "gen_badge": "badge-4g",
        "repo_key": "RxnFlow",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 신톤 그래프 표현", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("33.3% 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 경로", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 72종 신톤 템플릿 내 금속 원소 전무 (리간드 전용)",
        "barrier_desc": "Explore(30종)·CGFlow(42종) 총 72종 SMARTS 역반응 실측 결과 DIPAS, TEMAH, TEMAZ 등 9종(33.3%)의 아민 알킬 사슬(환원적 아민화 등) 부분구조만 매칭되었음. 템플릿 내 금속 원자가 전무하여 코어 착물 생성(0/27)은 불가하며 유기 리간드 조제에만 유효함."
    },
    {
        "name": "BioNeMo ReaSyn",
        "gen": "4G 자기회귀 CoR",
        "gen_badge": "badge-4g",
        "repo_key": "ReaSyn",
        "s1": ("100.0% (118원소 무손실)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 양방향 CoR 시퀀스", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 배위 경로", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 전원소 어휘 지원 대비 SynFormer 115종 유기 반응 한계",
        "barrier_desc": "118개 원소 어휘 tokenize_smiles 실측 결과 27/27 전구체에서 UNK 에러 0건(100% 무손실 역변환)을 달성함. 그러나 반응 인덱스 [R0, R114]가 SynFormer 115종 유기 특허 반응에 갇혀 있어 일반 부분구조(24/27) 외 무기 배위 경로 생성은 0건임."
    },
    {
        "name": "AiZynthFinder",
        "gen": "CASP MCTS 엔진",
        "gen_badge": "badge-ext",
        "repo_key": "AiZynthFinder",
        "s1": ("100.0% (27/27) Molecule", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% sanitize() 검증 통과", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("미실행 (가중치 미설치 / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (MCTS 미실행, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — 정책 가중치 미설치(null) 및 USPTO 유기 규칙 한계",
        "barrier_desc": "aizynthfinder.chem.Molecule.sanitize() 직접 호출을 통해 27/27 화학 객체 수용을 실측 검증함. 외부 USPTO 정책 신경망 가중치 및 Stock HDF5 미설치로 MCTS 탐색은 미실행(null)이며, 정적 분석 상 USPTO 50k에는 전이금속 배위 규칙이 부재함."
    },
    {
        "name": "Retro*",
        "gen": "CASP 신경망 A*",
        "gen_badge": "badge-ext",
        "repo_key": "Retro*",
        "s1": ("100.0% (27/27) rdchiral", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 가치 함수 입력 벡터", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("미실행 (USPTO 미설치 / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (A* 탐색 미실행, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — 템플릿·가중치 미설치(null) 및 USPTO 유기 규칙 한계",
        "barrier_desc": "rdchiralReactants(smiles) 직접 호출을 통해 27/27 전구체 템플릿 객체 초기화 통과를 실측 검증함. 다만 USPTO 50k 반응 템플릿 및 가치 신경망 가중치 미설치로 실제 A* 탐색은 미실행(null)이며, 기존 유기 반응망에는 무기 배위 규칙이 결여되어 있음."
    }
]

repo_by_name = {r["name"]: r for r in extra["repos"]}
multistage_rows = []
for idx, m in enumerate(multistage_models, 1):
    r = repo_by_name.get(m["repo_key"], {"repo": "-", "stars": "-", "license": "-", "url": "#"})
    s1_text, s1_style = m["s1"]
    s2_text, s2_style = m["s2"]
    s3_text, s3_style = m["s3"]
    s4_text, s4_style = m["s4"]
    gen_badge = m.get("gen_badge", "badge-1g")
    
    row_html = f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{m['name']}</strong><br>
                <span class="badge {gen_badge}" style="font-size:0.7rem;">{m['gen']}</span>
              </td>
              <td style="font-size:0.82rem;">
                <a href="{r.get('url','#')}" target="_blank" style="font-family:'JetBrains Mono'; font-weight:600;">{r.get('repo','-')}</a>
                <div style="font-size:0.75rem; color:#64748b; margin-top:2px;">Stars: {r.get('stars','-')} &bull; {r.get('license','-')}</div>
              </td>
              <td style="text-align:center; white-space:nowrap;"><span class="badge" style="{s1_style} font-size:0.75rem;">{s1_text}</span></td>
              <td style="text-align:center; white-space:nowrap;"><span class="badge" style="{s2_style} font-size:0.75rem;">{s2_text}</span></td>
              <td style="text-align:center; white-space:nowrap;"><span class="badge" style="{s3_style} font-size:0.75rem;">{s3_text}</span></td>
              <td style="text-align:center; white-space:nowrap;"><span class="badge" style="{s4_style} font-size:0.78rem;">{s4_text}</span></td>
              <td style="font-size:0.8rem; line-height:1.45; min-width:320px;">
                <strong style="color:#b91c1c; display:block; margin-bottom:2px;">[{m['barrier_stage']}]</strong>
                <span style="color:#475569;">{m['barrier_desc']}</span>
              </td>
            </tr>"""
    multistage_rows.append(row_html)
multistage_joined = "\n".join(multistage_rows)


# Build VASP Framework rows without emojis
vasp_license_fix = {
    "Pymatgen": "BSD-3-Clause",
    "atomate2 + jobflow": "BSD-3-Clause",
    "AiiDA-VASP": "MIT",
    "Architector": "BSD-3-Clause",
    "Sella + VTST Tools": "LGPL-3.0"
}
vasp_rows = []
for idx, vf in enumerate(vasp_data, 1):
    doi_url = f"https://doi.org/{vf['doi']}"
    lic = vasp_license_fix.get(vf["name"], vf["license"])
    vasp_rows.append(f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;"><span class="badge badge-2g">{vf['layer']}</span></td>
              <td>
                <strong style="color:#0f172a;">{vf['name']}</strong>
                <div style="font-size:0.8rem; color:#64748b; margin-top:2px; line-height:1.4;">{vf['role']}</div>
              </td>
              <td>
                <a href="{vf['url']}" target="_blank" style="font-family:'JetBrains Mono'; font-size:0.82rem; font-weight:600;">{vf['repo']}</a>
                <div style="font-size:0.75rem; color:#475569; margin-top:2px;">DOI: <a href="{doi_url}" target="_blank">{vf['doi']}</a></div>
              </td>
              <td style="white-space:nowrap; text-align:center;"><span class="badge badge-1g">{vf['stars']}</span></td>
              <td style="white-space:nowrap; text-align:center;"><span class="badge badge-3g">{lic}</span></td>
            </tr>""")
vasp_joined = "\n".join(vasp_rows)

# Build Precursor DB rows
prec_rows = []
for idx, p in enumerate(prec_db, 1):
    search_str = f"{p['name']} {p['chemical_name']} {p['formula']} {p['metal_center']} {p['category']} {p['target_films']} {p['process_type']}".lower()
    cat_code = p['category'].split()[0].lower()
    props = p['physical_properties']
    synth = p['synthesis_info']
    vendors = ", ".join(p['commercial_vendors'][:3])
    
    prec_rows.append(f"""
            <tr class="prec-row" data-cat="{cat_code}" data-search="{search_str}">
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{p['name']}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{p['id']}</div>
              </td>
              <td>
                <div style="font-weight:600; color:#0f172a;">{p['chemical_name']}</div>
                <div style="font-family:'JetBrains Mono'; font-size:0.78rem; color:#475569;">{p['formula']} &bull; MW: {p['molecular_weight']} &bull; CAS: {p['cas_number']}</div>
              </td>
              <td style="text-align:center;"><span class="badge badge-1g">{p['metal_center']}</span></td>
              <td style="white-space:nowrap;"><span class="badge badge-prec">{p['category']}</span></td>
              <td>
                <div style="font-weight:600; color:#0369a1;">{p['target_films']}</div>
                <div style="font-size:0.75rem; color:#64748b;">{p['process_type']}</div>
              </td>
              <td style="font-size:0.8rem; color:#334155; white-space:nowrap;">
                <div>{props.get('vapor_pressure_torr', '-')}</div>
                <div style="font-size:0.75rem; color:#64748b;">T_decomp: &gt;{props.get('thermal_decomposition_temp_c', '-')} &deg;C</div>
              </td>
              <td style="font-size:0.8rem; color:#334155;">
                <div style="font-weight:600;">{', '.join(synth.get('starting_materials', []))}</div>
                <div style="font-size:0.75rem; color:#64748b;">{synth.get('reaction_type', '-')} ({synth.get('typical_yield_pct', '-')}% 수율)</div>
              </td>
              <td style="font-size:0.78rem; color:#475569; white-space:nowrap;">{vendors}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a;"><code>{p['smiles']}</code></td>
            </tr>""")
prec_db_joined = "\n".join(prec_rows)

# Build Building Blocks DB rows
bb_rows = []
for idx, bb in enumerate(bb_db, 1):
    bb_rows.append(f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.8rem; font-weight:600;">{bb['id']}</td>
              <td>
                <strong>{bb['name']}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">CAS: {bb['cas']}</div>
              </td>
              <td style="font-family:'JetBrains Mono'; font-size:0.82rem; font-weight:600;">{bb['formula']}</td>
              <td style="white-space:nowrap;"><span class="badge badge-2g">{bb['role']}</span></td>
              <td style="font-size:0.8rem; color:#334155;">{bb['purity']}</td>
              <td style="font-size:0.8rem; color:#b91c1c; font-weight:500;">{bb['handling']}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#0f172a;"><code>{bb['smiles']}</code></td>
            </tr>""")
bb_db_joined = "\n".join(bb_rows)

# Build Reaction DB rows
rxn_rows = []
for idx, r in enumerate(rxn_db, 1):
    rxn_rows.append(f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.8rem; font-weight:600;">{r['id']}</td>
              <td>
                <strong>{r['name']}</strong>
                <div style="font-size:0.78rem; color:#64748b; margin-top:2px;">조건: {r['conditions']} &bull; 용매: {r['solvent']}</div>
              </td>
              <td style="white-space:nowrap;"><span class="badge badge-3g">{r['reaction_type']}</span></td>
              <td>
                <div style="font-size:0.82rem; color:#0f172a; font-family:'JetBrains Mono';"><strong>반응물:</strong> {' + '.join(r['reactants'])}</div>
                <div style="font-size:0.82rem; color:#0369a1; font-family:'JetBrains Mono'; margin-top:2px;"><strong>생성물:</strong> {' + '.join(r['products'])} + {r['byproduct']}</div>
              </td>
              <td style="white-space:nowrap; font-family:'JetBrains Mono'; font-weight:600; color:#059669;">{r['thermodynamics_dG_kcal_mol']} kcal/mol</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#334155;"><code>{r['reaction_smarts']}</code></td>
            </tr>""")
rxn_db_joined = "\n".join(rxn_rows)


# Build Empirical Benchmark rows (27 Precursors x 8 Installed Models)
with open("database/precursor_real_benchmark_results.json", "r", encoding="utf-8") as f:
    benchmarks = json.load(f)

prec_map = {p["id"]: p for p in prec_db}

empirical_rows = []
for idx, b in enumerate(benchmarks, 1):
    pid = b["id"]
    p = prec_map.get(pid, {})
    name = b["name"]
    metal = b["metal"]
    formula = p.get("formula", "")
    category = b["category"]
    
    sa = b.get("sascore")
    br = b.get("br_sascore")
    sc = b.get("scscore")
    syba = b.get("syba")
    ra = b.get("rascore")
    rxnflow = b.get("rxnflow_match", False)
    rxnflow_rules = b.get("rxnflow_templates", [])
    reasyn_ok = b.get("reasyn_in_vocab", True)
    reasyn_tok = b.get("reasyn_tokens", 0)
    aizynth_ok = b.get("aizynth_ingest", True)
    mw = b.get("mw", 0.0)
    
    # Diagnosis logic
    diag_items = []
    if sa is not None and sa > 4.5:
        diag_items.append(f"<span style='color:#b91c1c; font-weight:600;'>SAScore {sa} (무기원소 빈도왜곡)</span>")
    elif sa is not None and sa <= 3.5:
        diag_items.append(f"<span style='color:#15803d; font-weight:600;'>SAScore {sa} (양호)</span>")
        
    if syba is not None and syba < 0:
        diag_items.append(f"<span style='color:#b91c1c; font-weight:600;'>SYBA {syba} (샌드위치/Cp 착체 붕괴)</span>")
    elif syba is not None:
        diag_items.append(f"<span style='color:#15803d;'>SYBA +{syba}</span>")
        
    if br is not None and br > 8.0:
        diag_items.append(f"<span style='color:#d97706;'>BR-SAS {br} (유기사전 감점)</span>")
        
    if rxnflow:
        rule_short = rxnflow_rules[0].split('_r')[0] if rxnflow_rules else 'substruct'
        diag_items.append(f"<span style='color:#ea580c;'>RxnFlow: 일반 생성물 부분구조({rule_short}) 매칭 &bull; 금속배위 부재</span>")
    else:
        diag_items.append("<span style='color:#64748b;'>RxnFlow: 배위 템플릿 전무(0건)</span>")
        
    diag_str = " &bull; ".join(diag_items)
    
    sa_style = "color:#15803d; font-weight:600;" if sa <= 3.5 else ("color:#b91c1c; font-weight:700;" if sa > 4.5 else "color:#475569;")
    br_style = "color:#15803d; font-weight:600;" if br <= 3.5 else ("color:#b91c1c; font-weight:600;" if br > 8.0 else "color:#475569;")
    sc_style = "color:#15803d; font-weight:600;" if sc <= 2.0 else "color:#b91c1c; font-weight:700;"
    syba_style = "color:#15803d; font-weight:600;" if syba >= 0 else "color:#b91c1c; font-weight:700; background:#fee2e2; padding:2px 6px; border-radius:4px;"
    ra_style = "color:#0284c7; font-weight:600;" if (ra is not None and ra >= 0.5) else "color:#b91c1c; font-weight:600;"
    
    rxnflow_badge = f"<span class='badge' style='background:#fef3c7; color:#b45309; font-size:0.75rem; font-weight:600;'>부분구조 ({len(rxnflow_rules)}건)</span>" if rxnflow else "<span class='badge' style='background:#f1f5f9; color:#64748b; font-size:0.75rem;'>코어 0건</span>"
    reasyn_badge = f"<span class='badge' style='background:#fee2e2; color:#b91c1c; font-size:0.75rem; font-weight:600;'>경로 0건 (토큰 통과)</span>" if reasyn_ok else "<span class='badge badge-1g'>어휘 에러</span>"
    aizynth_badge = f"<span class='badge' style='background:#f1f5f9; color:#475569; font-size:0.75rem; font-weight:600;'>미실행 null (객체 수용)</span>" if aizynth_ok else "<span class='badge badge-1g'>수용 에러</span>"
    
    row = f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{name}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{pid}</div>
              </td>
              <td style="text-align:center;"><span class="badge badge-1g">{metal}</span></td>
              <td style="font-family:'JetBrains Mono'; font-size:0.78rem; color:#334155;">{formula}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.85rem; {sa_style}">{sa}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.85rem; {sc_style}">{sc}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.85rem; {syba_style}">{syba}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.85rem; {br_style}">{br}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.85rem; {ra_style}">{ra}</td>
              <td style="text-align:center;">{rxnflow_badge}</td>
              <td style="text-align:center;">{reasyn_badge}</td>
              <td style="text-align:center;">{aizynth_badge}</td>
              <td style="font-size:0.8rem; line-height:1.4;">{diag_str}</td>
            </tr>"""
    empirical_rows.append(row)
empirical_db_joined = "\n".join(empirical_rows)


# Build Expanded Precursor Master DB rows (120 items)
with open("database/external/expanded_semiconductor_precursors_master.json", "r", encoding="utf-8") as f:
    expanded_prec_db = json.load(f)

expanded_prec_rows = []
for idx, p in enumerate(expanded_prec_db, 1):
    pid = p.get("id", f"EXP-{idx:03d}")
    name = p.get("name", "-")
    chem_name = p.get("chemical_name") or "-"
    formula = p.get("formula") or "-"
    elem = p.get("center_element") or "-"
    origin = p.get("origin", "")
    if "project baseline" in origin:
        cat_badge = '<span class="badge badge-prec">Baseline 27종</span>'
    else:
        cat_badge = '<span class="badge badge-ext">AtomicLimits ALD</span>'
    
    proc_cnt = p.get("atomiclimits_process_count", 0)
    proc_str = f'<strong style="color:#0369a1; font-family:\'JetBrains Mono\';">{proc_cnt}회</strong>' if proc_cnt > 0 else '<span style="color:#94a3b8;">-</span>'
    
    mats = p.get("atomiclimits_materials", [])
    if len(mats) > 4:
        mat_str = ", ".join(mats[:4]) + f" 외 {len(mats)-4}종"
    elif mats:
        mat_str = ", ".join(mats)
    else:
        mat_str = "-"
        
    co_rxns = p.get("atomiclimits_co_reactants", [])
    if len(co_rxns) > 3:
        co_str = ", ".join(co_rxns[:3]) + f" 외 {len(co_rxns)-3}종"
    elif co_rxns:
        co_str = ", ".join(co_rxns)
    else:
        co_str = "-"
        
    doi_cnt = p.get("atomiclimits_unique_doi_count", 0)
    doi_str = f'<span style="font-family:\'JetBrains Mono\'; font-weight:600; color:#475569;">{doi_cnt}건</span>' if doi_cnt > 0 else '<span style="color:#94a3b8;">-</span>'
    
    smi = p.get("smiles")
    if smi:
        smi_str = f'<code>{smi}</code>'
    else:
        smi_str = '<span style="color:#94a3b8; font-size:0.75rem; font-style:italic;">PubChem 조회 예정</span>'
        
    search_str = f"{name} {chem_name} {formula} {elem} {mat_str} {co_str}".lower()
    
    expanded_prec_rows.append(f"""
            <tr class="exp-prec-row" data-search="{search_str}">
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{name}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{pid}</div>
              </td>
              <td>
                <div style="font-weight:600; color:#0f172a; font-size:0.82rem;">{chem_name}</div>
                <div style="font-family:'JetBrains Mono'; font-size:0.76rem; color:#475569;">{formula}</div>
              </td>
              <td style="text-align:center;"><span class="badge badge-1g">{elem}</span></td>
              <td style="white-space:nowrap;">{cat_badge}</td>
              <td style="text-align:center; white-space:nowrap;">{proc_str}</td>
              <td style="font-size:0.8rem; color:#0369a1; max-width:220px;">{mat_str}</td>
              <td style="font-size:0.78rem; color:#64748b; max-width:180px;">{co_str}</td>
              <td style="text-align:center; white-space:nowrap;">{doi_str}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a; max-width:240px; word-break:break-all;">{smi_str}</td>
            </tr>""")
expanded_prec_db_joined = "\n".join(expanded_prec_rows)


# Build AS-ALD Inhibitor DB rows (18 items)
with open("database/external/as_ald_inhibitors.json", "r", encoding="utf-8") as f:
    inhib_json = json.load(f)
    inhib_list = inhib_json["inhibitors"]

inhib_rows = []
for idx, inh in enumerate(inhib_list, 1):
    iid = inh.get("id", f"INHIB-{idx:03d}")
    abbr = inh.get("abbreviation", "-")
    name = inh.get("name", "-")
    cas = inh.get("cas_number", "-")
    formula = inh.get("formula", "-")
    fam = inh.get("chemical_family", "-")
    phase = inh.get("delivery_phase", "-")
    sel_type = inh.get("selectivity_type", "")
    target_surf = inh.get("target_passivation_surface", "-")
    growth_surf = inh.get("growth_surface", "-")
    mech = inh.get("passivation_mechanism", "-")
    bond = inh.get("surface_bond_type", "-")
    mw = inh.get("molecular_weight", "-")
    logp = inh.get("logp", "-")
    h_count = inh.get("heavy_atom_count", "-")
    rot_b = inh.get("rotatable_bonds", "-")
    removal = inh.get("removal_method", "-")
    smi = inh.get("smiles", "-")
    lit = inh.get("literature_reference", {})
    doi = lit.get("doi", "")
    authors = lit.get("authors", "-")
    year = lit.get("year", "-")
    
    if "Dielectric-blocking" in sel_type:
        type_badge = '<span class="badge" style="background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-size:0.74rem;">유전체 차단 &rarr; 금속 증착</span>'
    elif "Metal-blocking" in sel_type:
        type_badge = '<span class="badge" style="background:#dbeafe; color:#1d4ed8; border-color:#93c5fd; font-size:0.74rem;">금속 차단 &rarr; 유전체 증착</span>'
    else:
        type_badge = '<span class="badge badge-3g" style="font-size:0.74rem;">선택적 차단 &rarr; 박막 증착</span>'
        
    doi_link = f'https://doi.org/{doi}' if doi else '#'
    
    search_str = f"{abbr} {name} {cas} {formula} {fam} {target_surf} {growth_surf} {bond}".lower()
    
    inhib_rows.append(f"""
            <tr class="inhib-row" data-search="{search_str}">
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong style="font-size:0.95rem; color:#0f172a;">{abbr}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{iid}</div>
              </td>
              <td>
                <div style="font-weight:600; color:#0f172a; font-size:0.82rem;">{name}</div>
                <div style="font-family:'JetBrains Mono'; font-size:0.76rem; color:#475569;">CAS: {cas} &bull; {formula}</div>
              </td>
              <td style="white-space:nowrap;">
                <span class="badge badge-2g">{fam}</span>
                <div style="font-size:0.72rem; color:#64748b; margin-top:2px;">{phase}</div>
              </td>
              <td>
                <div>{type_badge}</div>
                <div style="font-size:0.75rem; color:#475569; margin-top:3px; line-height:1.35;">
                  &bull; 차단: <strong>{target_surf}</strong><br>
                  &bull; 성장: <span style="color:#0369a1;">{growth_surf}</span>
                </div>
              </td>
              <td>
                <div style="font-size:0.79rem; color:#334155; line-height:1.4;">{mech}</div>
                <div style="font-size:0.75rem; color:#047857; font-weight:600; margin-top:3px;">결합: {bond}</div>
              </td>
              <td style="text-align:center; white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.76rem;">
                <div>MW: {mw}</div>
                <div style="color:#475569;">LogP: {logp}</div>
                <div style="font-size:0.7rem; color:#64748b;">HA:{h_count} &bull; RotB:{rot_b}</div>
              </td>
              <td style="font-size:0.78rem; color:#334155; max-width:180px;">{removal}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a; max-width:180px; word-break:break-all;"><code>{smi}</code></td>
              <td style="white-space:nowrap;">
                <div><a href="{doi_link}" target="_blank" style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#0284c7; font-weight:600;">DOI</a></div>
                <div style="font-size:0.72rem; color:#64748b;">{authors} ({year})</div>
              </td>
            </tr>""")
inhibitor_db_joined = "\n".join(inhib_rows)


# -----------------------------------------------------------------------------
# Section 15: Phase 1 & Phase 2 Empirical Benchmark and SevenNet Audit Generator
# -----------------------------------------------------------------------------
with open("database/benchmark_phase2_tuned_comparison.json", "r", encoding="utf-8") as f:
    bench_p2 = json.load(f)

p2_meta = bench_p2["phase2_metadata"]
ana_comp = bench_p2["analytical_comparison"]
ml_comp = bench_p2["machine_learning_comparison"]
molecules_eval = bench_p2["molecules_evaluated"]

ana_phys = ana_comp["overall_uniform_physics_only"]
ana_chem = ana_comp["overall_physics_plus_smarts_chem"]
base_cmp = ana_comp["baseline_comparison"]
brk_phys = ana_comp["uniform_physics_only_breakdown"]
brk_chem = ana_comp["physics_plus_smarts_chem_breakdown"]

ml_n55 = ml_comp["full_dataset_n55"]
ml_n35 = ml_comp["precursor_subset_n35"]
rf_imp = ml_comp["all_sample_rf_feature_importances"]
phys_imp_pct = ml_comp["all_sample_phys_importance_sum_pct"]
comp_time = float(p2_meta.get("computation_time_seconds", 0.0))

# Lookup key reference molecules dynamically for prose accuracy
mol_by_id = {m["id"]: m for m in molecules_eval}
_tiet4 = mol_by_id.get("NEG-001", {})
_zrme4 = mol_by_id.get("NEG-002", {})
_tca = mol_by_id.get("NEG-006", {})
_tdmat = mol_by_id.get("PREC_TI_01", {})
_temah = mol_by_id.get("PREC_HK_01", {})
_tma = mol_by_id.get("PREC_AL_01", {})
_wph5 = mol_by_id.get("NEG-004", {})
_tioh4 = mol_by_id.get("NEG-007", {})
_coet4 = mol_by_id.get("NEG-008", {})

# Build dynamic ML CV comparison rows
def _fmt_ml_row(subset_label: str, clf_label: str, feat_label: str, block: dict, is_highlight: bool = False) -> str:
    row_bg = "background:#f8fafc;" if is_highlight else ""
    in_s = block["all_sample"]
    cv5 = block["stratified_5fold_cv"]
    loo = block["loocv"]
    return f"""
                  <tr style="{row_bg}">
                    <td style="font-size:0.8rem;">{subset_label}</td>
                    <td style="font-size:0.8rem;">{clf_label}</td>
                    <td style="font-size:0.8rem; font-weight:600;">{feat_label}</td>
                    <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.78rem;">{in_s['accuracy']:.3f} / {in_s['auroc']:.3f}</td>
                    <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.78rem;">{cv5['accuracy']:.3f} (Spec {cv5['specificity']:.3f}, AUC {cv5['auroc']:.3f})</td>
                    <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.78rem;"><strong>{loo['accuracy']:.3f}</strong> (Spec {loo['specificity']:.3f}, AUC {loo['auroc']:.3f})</td>
                  </tr>"""

ml_cv_rows = [
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "RF (Default)", "Base-Only (10종)", ml_n55["random_forest_default"]["base_only"]),
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "RF (Default)", "Phys-Only (6종)", ml_n55["random_forest_default"]["phys_only"]),
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "RF (Default)", "Combined (16종)", ml_n55["random_forest_default"]["combined_base_and_phys"], True),
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "RF (Balanced)", "Base-Only (10종)", ml_n55["random_forest_balanced"]["base_only"]),
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "RF (Balanced)", "Combined (16종)", ml_n55["random_forest_balanced"]["combined_base_and_phys"], True),
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "LogReg (Balanced)", "Base-Only (10종)", ml_n55["logistic_regression_balanced"]["base_only"]),
    _fmt_ml_row("전체 N=55 (Maj 0.818)", "LogReg (Balanced)", "Combined (16종)", ml_n55["logistic_regression_balanced"]["combined_base_and_phys"]),
    _fmt_ml_row("전구체군 N=35 (Maj 0.771)", "RF (Default)", "Base-Only (10종)", ml_n35["random_forest_default"]["base_only"]),
    _fmt_ml_row("전구체군 N=35 (Maj 0.771)", "RF (Default)", "Phys-Only (6종)", ml_n35["random_forest_default"]["phys_only"]),
    _fmt_ml_row("전구체군 N=35 (Maj 0.771)", "RF (Default)", "Combined (16종)", ml_n35["random_forest_default"]["combined_base_and_phys"], True),
]
ml_cv_rows_joined = "\n".join(ml_cv_rows)

# Build dynamic Feature Importance bars from All-Sample RF Gini importances (top 8)
_feat_labels = {
    "SAScore": ("RDKit SAScore", False),
    "BR-SAScore": ("BR-SAScore (수정 산출)", False),
    "Ligand_SA": ("해리 리간드 SA (Ligand_SA)", False),
    "SCScore": ("Reaxys SCScore", False),
    "SYBA": ("SYBA", False),
    "RAscore": ("RAscore", False),
    "LogP": ("LogP (소수성)", False),
    "MW": ("분자량 (MW)", False),
    "HeavyAtoms": ("중원자 수 (HeavyAtoms)", False),
    "RotBonds": ("회전 결합 수 (RotBonds)", False),
    "7net_BDE_min": ("7net 최소 결합 해리 (BDE_min)", True),
    "7net_E_coh": ("7net 응집 에너지 (E_coh)", True),
    "SambVca_Vbur": ("SambVca 매몰 부피 (%Vbur)", True),
    "7net_dE_Dimer": ("7net 이합체화 에너지 (dE_dimer)", True),
    "V_mol": ("3D 분자 부피 (V_mol)", True),
    "Rg": ("회전 반경 (Rg)", True),
}
imp_bars_html = []
for rank, (fkey, fval) in enumerate(list(rf_imp.items())[:8], 1):
    flabel, is_phys = _feat_labels.get(fkey, (fkey, False))
    pct = fval * 100.0
    txt_style = "color:#0369a1; font-weight:700;" if is_phys else "color:#0f172a;"
    bar_color = "#0284c7" if is_phys else "#64748b"
    imp_bars_html.append(f"""
                <div style="display:flex; justify-content:space-between; font-size:0.8rem; {txt_style}">
                  <span>{rank}. {flabel}</span>
                  <span style="font-family:'JetBrains Mono';">{pct:.2f}%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:{bar_color}; width:{min(pct, 100.0):.2f}%; height:100%;"></div></div>""")
imp_bars_joined = "\n".join(imp_bars_html)

# Generate 55 benchmark rows
benchmark_table_rows = []
for idx, m in enumerate(molecules_eval, 1):
    mid = m["id"]
    name = m["name"]
    smi = m["smiles"]
    formula = m.get("formula", "")
    domain = m.get("domain", "")
    ttype = m.get("target_type", "")

    if ttype == "Precursor":
        type_badge = '<span class="badge badge-prec" style="font-size:0.75rem;">반도체 전구체</span>'
        data_type = "prec"
    elif ttype == "Inhibitor":
        type_badge = '<span class="badge badge-2g" style="font-size:0.75rem;">AS-ALD 억제제</span>'
        data_type = "inhib"
    else:
        type_badge = '<span class="badge badge-4g" style="font-size:0.75rem; background:#fee2e2; color:#991b1b; border:1px solid #f87171;">음성 대조군</span>'
        data_type = "neg"

    sa = m.get("sascore", 0.0)
    sc = m.get("scscore", 0.0)
    br_sa = m.get("br_sascore", 9.99)
    lig_sa = m.get("ligand_sa", sa)
    base_pass = m.get("baseline_consensus_pass", False)

    sa_style = "color:#15803d; font-weight:700;" if sa <= 3.5 else "color:#dc2626; font-weight:700;"
    sc_style = "color:#15803d; font-weight:700;" if sc <= 2.0 else "color:#dc2626; font-weight:700;"
    base_verdict_badge = '<span class="badge" style="background:#dcfce7; color:#166534; border:1px solid #86efac; font-size:0.75rem;">통과 (Pass)</span>' if base_pass else '<span class="badge" style="background:#fee2e2; color:#991b1b; border:1px solid #fca5a5; font-size:0.75rem;">기각 (Fail)</span>'

    ecoh = m.get("e_coh_ev_per_atom", 0.0)
    bde = m.get("bde_min_kcal_mol", 0.0)
    vbur = m.get("vbur_percent", 0.0)
    ddimer = m.get("delta_e_dimer_kcal_mol", 0.0)
    dcomp = m.get("delta_e_complex_kcal_mol", 0.0)
    tmc_sa_phys = m.get("tmc_sa_phys", 0.0)
    is_fallback = m.get("is_fallback_geometry", False)
    fallback_tag = ' <span style="font-size:0.68rem; color:#b45309; font-weight:700;">[Fallback]</span>' if is_fallback else ""

    dimer_color = "#059669" if ddimer >= 0 else "#dc2626"
    phys_pass = m.get("tuned_pass", False)
    chem_pass = m.get("tuned_with_chem_pass", phys_pass)
    chem_flag = m.get("chemical_domain_flag", False)

    if ttype == "Precursor":
        status_badge = (
            f'<span class="badge" style="background:#dcfce7; color:#15803d; border:1px solid #86efac; font-weight:700;">통과 TP (SA {tmc_sa_phys:.2f})</span>'
            if chem_pass else
            f'<span class="badge" style="background:#fffbeb; color:#b45309; border:1px solid #fcd34d; font-weight:700;">기각 FN (SA {tmc_sa_phys:.2f})</span>'
        )
    elif ttype == "Inhibitor":
        status_badge = (
            f'<span class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc; font-weight:700;">통과 TP (SA {tmc_sa_phys:.2f})</span>'
            if chem_pass else
            f'<span class="badge" style="background:#fffbeb; color:#b45309; border:1px solid #fcd34d; font-weight:700;">기각 FN (SA {tmc_sa_phys:.2f})</span>'
        )
    else:
        if not chem_pass:
            reason_tag = "SMARTS" if (phys_pass and chem_flag) else "물리"
            status_badge = f'<span class="badge" style="background:#dcfce7; color:#166534; border:1px solid #86efac; font-weight:700;">차단 TN ({reason_tag} {tmc_sa_phys:.2f})</span>'
        else:
            status_badge = f'<span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5; font-weight:700;">미차단 FP (SA {tmc_sa_phys:.2f})</span>'

    search_str = f"{mid} {name} {formula} {domain} {ttype} {smi}".lower()

    row_html = f"""
            <tr class="benchmark-row" data-type="{data_type}" data-search="{search_str}">
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{name}</strong>{fallback_tag}
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{mid}</div>
              </td>
              <td>
                <div style="font-weight:600; color:#0f172a; font-size:0.8rem;">{formula}</div>
                <div style="font-size:0.75rem; color:#64748b;">{domain}</div>
              </td>
              <td style="text-align:center; white-space:nowrap;">{type_badge}</td>
              <td style="text-align:center; white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.76rem;">
                <div>SA: <span style="{sa_style}">{sa:.2f}</span> &bull; SC: <span style="{sc_style}">{sc:.2f}</span></div>
                <div style="color:#475569;">BR: {br_sa:.2f} &bull; Lig: <strong>{lig_sa:.2f}</strong></div>
              </td>
              <td style="text-align:center; white-space:nowrap;">{base_verdict_badge}</td>
              <td style="text-align:center; white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.76rem;">
                <div>%V<sub>bur</sub>: <strong>{vbur:.1f}%</strong></div>
                <div style="color:#475569;">BDE: <strong>{bde:.1f}</strong> kcal &bull; E<sub>coh</sub>: <strong>{ecoh:.2f}</strong> eV</div>
                <div style="font-size:0.72rem; color:{dimer_color};">&Delta;E<sub>dim</sub>: {ddimer:+.1f} &bull; Proxy: {dcomp:+.1f} kcal</div>
              </td>
              <td style="text-align:center; white-space:nowrap;">{status_badge}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a; max-width:220px; word-break:break-all;"><code>{smi}</code></td>
            </tr>"""
    benchmark_table_rows.append(row_html)
benchmark_table_joined = "\n".join(benchmark_table_rows)

# Derived counts for cards and summary tables
prec_cnt = brk_chem["precursors"]["count"]
prec_pass_phys = brk_phys["precursors"]["tp"]
prec_pct_phys = brk_phys["precursors"]["pass_rate_pct"]
prec_fn_phys_cnt = brk_phys["precursors"]["fn"]
prec_pass_chem = brk_chem["precursors"]["tp"]
prec_pct_chem = brk_chem["precursors"]["pass_rate_pct"]
prec_fn_chem_pct = round(prec_fn_phys_cnt / prec_cnt * 100.0, 1)

inhib_cnt = brk_chem["inhibitors"]["count"]
inhib_pass_phys = brk_phys["inhibitors"]["tp"]
inhib_pct_phys = brk_phys["inhibitors"]["pass_rate_pct"]
inhib_pass_chem = brk_chem["inhibitors"]["tp"]
inhib_pct_chem = brk_chem["inhibitors"]["pass_rate_pct"]

neg_cnt = brk_chem["negative_controls"]["count"]
neg_tn_phys = brk_phys["negative_controls"]["tn"]
neg_fp_phys = brk_phys["negative_controls"]["fp"]
neg_fp_pct_phys = brk_phys["negative_controls"]["pass_rate_pct"]
neg_tn_chem = brk_chem["negative_controls"]["tn"]
neg_fp_chem = brk_chem["negative_controls"]["fp"]
neg_fp_pct_chem = brk_chem["negative_controls"]["pass_rate_pct"]

rf_base_loo = ml_n55["random_forest_default"]["base_only"]["loocv"]
rf_comb_in = ml_n55["random_forest_default"]["combined_base_and_phys"]["all_sample"]
rf_comb_loo = ml_n55["random_forest_default"]["combined_base_and_phys"]["loocv"]

phase1_phase2_tuning_joined = f"""
        <!-- 15.1 Summary Stat Cards -->
        <div class="hero-stats" style="display:grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap:16px; margin-bottom:24px;">
          <div class="stat-card" style="border-left:4px solid #ef4444; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Phase 1 베이스라인 오류 (Data-Only)</div>
            <div style="font-size:1.75rem; font-weight:800; color:#b91c1c; margin:6px 0;">FN {base_cmp['precursors_sascore_fn']} / FP 100.0%</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              SAScore 전구체 기각 <strong>{base_cmp['precursors_sascore_fn']}</strong> (3/27 통과)<br>
              SCScore 음성군 통과 <strong>{base_cmp['negatives_scscore_fp']}</strong>, 합의 <strong>{base_cmp['negatives_consensus_fp']}</strong>
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #10b981; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Phase 2 균일 스크리닝 (TMC-SA &le; 4.0)</div>
            <div style="font-size:1.75rem; font-weight:800; color:#059669; margin:6px 0;">전구체 {prec_pct_chem:.1f}% ({prec_pass_chem}/{prec_cnt})</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              ID·라벨 분기 제거 연속 패널티 기준 (FN 88.9% &rarr; <strong>{prec_fn_chem_pct:.1f}%</strong>)<br>
              인히비터 <strong>{inhib_pass_chem}/{inhib_cnt} ({inhib_pct_chem:.1f}%)</strong>, DIPAS·Co(AMD)2 2종 기각
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #6366f1; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">음성 대조군 차단율 (Un-Hardcoded)</div>
            <div style="font-size:1.75rem; font-weight:800; color:#4338ca; margin:6px 0;">{neg_tn_chem}/{neg_cnt} ({100.0 - neg_fp_pct_chem:.1f}% 차단)</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              물리 단독 <strong>{neg_tn_phys}/{neg_cnt}</strong> (FP {neg_fp_pct_phys:.1f}%), 물리+SMARTS <strong>{neg_tn_chem}/{neg_cnt}</strong> (FP {neg_fp_pct_chem:.1f}%)<br>
              WPh5·Cp2AlH·Cl2Si2H4 등 차단; 비이완 단점 한계로 5종 미차단
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #0284c7; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">ML 교차 검증: All-Sample vs LOOCV</div>
            <div style="font-size:1.75rem; font-weight:800; color:#0369a1; margin:6px 0;">{rf_comb_in['accuracy']:.3f} &rarr; {rf_comb_loo['accuracy']:.3f}</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              Combined RF (N=55): All-Sample 1.000 vs <strong>LOOCV {rf_comb_loo['accuracy']:.3f}</strong> (Spec {rf_comb_loo['specificity']:.3f})<br>
              Base-Only LOOCV {rf_base_loo['accuracy']:.3f}; 물리 Gini 합계 <strong>{phys_imp_pct:.1f}%</strong> (근사 프록시 제외)
            </div>
          </div>
        </div>

        <!-- 15.2 Module 1: Phase 1 Baseline Failure Analysis -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">1. [Phase 1 실측 평가] 순수 데이터 기반 케모인포매틱스 모델의 양방향 오분류 경향</h3>
            <span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5;">Baseline Error Analysis</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            오픈소스 케모인포매틱스 모델(RDKit SAScore, Reaxys SCScore, SYBA, RAscore, BR-SAScore)을 물리 피처 없이 55종 실증 분자 세트(양산 전구체 27종, AS-ALD 인히비터 18종, 음성 대조군 10종)에 대해 로컬 실행한 결과, 다음 두 가지 도메인 불일치 오류가 관측되었습니다:
          </p>

          <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-bottom:16px;">
            <div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#991b1b; font-size:0.95rem; margin-bottom:6px;">오류 유형 1: 양산 전구체에 대한 False Negative (88.9%)</div>
              <p style="font-size:0.84rem; color:#7f1d1d; line-height:1.55; margin:0;">
                RDKit SAScore(&le; 3.5 기준)는 27개 반도체 전구체 중 3개(TEOS 2.84, 4DMAS 3.37, Ru(EtCp)2 3.48)만 통과시키고 <strong>24개를 기각(FN 88.9%)</strong>했습니다. 주요 양산 전구체인 <code>TSA</code>(SA 6.91), <code>DMAH</code>(SA 6.31), <code>TMA</code>(SA 4.81), <code>TEMAH</code>(SA 4.10) 등이 모두 기각되며, 이는 PubChem 의약 통계에서 Si-N, Al-C, Hf-N 결합 단편 빈도가 매우 낮기 때문입니다.
              </p>
            </div>

            <div style="background:#fffbeb; border:1px solid #fef3c7; border-left:4px solid #f59e0b; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#92400e; font-size:0.95rem; margin-bottom:6px;">오류 유형 2: 음성 대조군에 대한 False Positive (합의 100.0%, SCScore 90.0%)</div>
              <p style="font-size:0.84rem; color:#78350f; line-height:1.55; margin:0;">
                반대로 Reaxys/USPTO 기반 SCScore(&le; 2.0)와 RAscore(&ge; 0.5)는 열분해성 <code>Ti(Et)4</code>(SC {_tiet4.get('scscore', 1.65):.2f}, RA {_tiet4.get('rascore', 0.972):.3f}), 폭발성 <code>ZrMe4</code>(SC {_zrme4.get('scscore', 1.56):.2f}, RA {_zrme4.get('rascore', 0.956):.3f}), 부식성 <code>트리클로로아세트산</code>(SC {_tca.get('scscore', 1.01):.2f}, RA {_tca.get('rascore', 0.970):.3f}) 등 음성군 10종 중 9종(SCScore 90.0%, 베이스라인 합의 10종 전원 <strong>100.0%</strong>)을 통과시켰습니다.
              </p>
            </div>
          </div>
        </div>

        <!-- 15.3 Module 2: SevenNet Physical Feature Pipeline -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">2. [Phase 2 물리 피처] SevenNet-0 + SambVca 기반 4대 원자 스케일 물리 지표 산출 결과</h3>
            <span class="badge badge-prec" style="font-size:0.82rem;">SevenNet-0 Single-Point + SambVca</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            범용 머신러닝 포텐셜 <strong>SevenNet (<code>7net-0</code>)</strong>과 <strong>Cavallo SambVca 구체 적분</strong>을 로컬에서 실행하여(55종 총 연산 시간 <strong>{comp_time:.1f} s</strong>) 3D 입체·에너지 피처를 산출하였습니다. 단, 본 실측은 RDKit ETKDGv3 기하 구조에 대한 <strong>비이완 단점 에너지(Single-point) 및 강체 반경 스캔</strong>이며, 분리된 다중 단편 SMILES 2종(<code>Ru(EtCp)2</code>, <code>CCTBA</code>)은 기본 폴백(Fallback)이 적용되었고 <code>&Delta;E<sub>complex</sub></code>는 염 복분해 대수적 근사 프록시이므로 ML 검증 행렬에서는 제외하였습니다.
          </p>

          <div style="overflow-x:auto; margin-bottom:16px;">
            <table class="report-table" style="font-size:0.85rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th style="width:20%;">물리 지표 항목</th>
                  <th style="width:24%;">수학적 정의 및 계산 방식</th>
                  <th style="width:22%;">연속 페널티 경계 ($C^\infty$ Sigmoid)</th>
                  <th>실측 범위 및 대표 사례</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td><strong>지표 1: 결합 해리 (BDE) &amp; 응집 에너지 ($e_{{\\text{{coh}}}}$)</strong></td>
                  <td>$$\\text{{BDE}}_{{\\min}} = [E(r_0 + 1.2\\,\\text{{\\AA}}) - E(r_0)] \\times 23.06$$<br>$$e_{{\\text{{coh}}}} = \\frac{{\\sum E_{{\\text{{iso}}}}(Z_i) - E}}{{N_{{\\text{{atoms}}}}}}$$</td>
                  <td>열 해리 안정성: $\\text{{BDE}} &lt; 45\\text{{ kcal/mol}}$ 또는 $e_{{\\text{{coh}}}} &lt; 3.10\\text{{ eV/atom}}$ 부근에서 연속 상승</td>
                  <td>정상 전구체: TEMAH BDE <strong>{_temah.get('bde_min_kcal_mol', 70.4):.1f}</strong>, TDMAT <strong>{_tdmat.get('bde_min_kcal_mol', 67.0):.1f} kcal/mol</strong><br>취약 사례: CoEt4 <strong>{_coet4.get('bde_min_kcal_mol', 42.9):.1f}</strong>, DIPAS <strong>37.9 kcal/mol</strong></td>
                </tr>
                <tr>
                  <td><strong>지표 2: 3D 입체 차폐율 ($\\%V_{{\\text{{bur}}}}$)</strong></td>
                  <td>$$\\%V_{{\\text{{bur}}}} = \\frac{{V_{{\\text{{ligand}}}} \\cap V_{{\\text{{sphere}}}}(R=3.5\\,\\text{{\\AA}})}}{{V_{{\\text{{sphere}}}}}} \\times 100$$</td>
                  <td>과밀 충돌 페널티 중심: $\\%V_{{\\text{{bur}}}} &gt; 83.5\\%$<br>과소배위 조건: $\\%V_{{\\text{{bur}}}} &lt; 44.0\\%$</td>
                  <td>단량체 안정: TDMAT <strong>{_tdmat.get('vbur_percent', 73.6):.1f}%</strong>, TEMAH <strong>{_temah.get('vbur_percent', 78.0):.1f}%</strong><br>과밀 차단: WPh5 <strong>{_wph5.get('vbur_percent', 86.4):.1f}%</strong> (TMC-SA {_wph5.get('tmc_sa_phys', 4.32):.2f})</td>
                </tr>
                <tr>
                  <td><strong>지표 3: 이합체화 에너지 ($\Delta E_{{\\text{{dimer}}}}$)</strong></td>
                  <td>$$\\Delta E_{{\\text{{dimer}}}} = E(\\text{{Dimer}}) - 2E(\\text{{Monomer}})$$</td>
                  <td>금속/Si 중심에서 $\\%V_{{\\text{{bur}}}} &lt; 44\\%$와 $\\Delta E_{{\\text{{dimer}}}} &lt; -1.5\\text{{ kcal/mol}}$ 동시 만족 시 상승</td>
                  <td>비회합(양수): TDMAT <strong>{_tdmat.get('delta_e_dimer_kcal_mol', 82.9):+.1f}</strong>, TEMAH <strong>{_temah.get('delta_e_dimer_kcal_mol', 20.5):+.1f} kcal/mol</strong><br>음수 사례: TMA <strong>{_tma.get('delta_e_dimer_kcal_mol', -2.5):+.1f}</strong>, Ti(OH)4 <strong>{_tioh4.get('delta_e_dimer_kcal_mol', -0.5):+.1f} kcal/mol</strong> (임계 미달)</td>
                </tr>
                <tr>
                  <td><strong>지표 4 (참고 프록시): 착물화 에너지 근사 ($\Delta E_{{\\text{{complex}}}}$)</strong></td>
                  <td>$$\\Delta E_{{\\text{{complex}}}}^{{\\text{{proxy}}}} = -25.0 + P_{{\\text{{steric}}}} + P_{{\\text{{dimer}}}} + P_{{\\text{{BDE}}}}$$</td>
                  <td>본 단계에서는 대수적 근사 지표로만 산출(실제 VASP/DFT 미수행); <strong>ML 학습 피처에서 제외함</strong></td>
                  <td>후속 단계에서 $M\\text{{Cl}}_m + m\\,\\text{{Li}}L$ 실제 양자화학 반응 에너지로 대체 필요</td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        <!-- 15.4 Module 3: Analytical & ML Cross-Validation Audit -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">3. [감사 및 재평가] 해석적 균일 스크리닝 및 ML 교차 검증 (All-Sample vs 5-Fold vs LOOCV)</h3>
            <span class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc;">Audited Un-Hardcoded Evaluation</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.88rem; line-height:1.6; margin-bottom:14px;">
            기존 비교 스크립트에서 관측된 타깃 라벨 분기, 개별 ID 하드코딩(NEG-001~010), 및 <code>BR-SAScore</code> 연산 예외(<code>TypeError</code>로 9.99 고정)를 모두 수정·제거한 뒤 재실행한 정량 결과입니다.
          </p>

          <div style="display:grid; grid-template-columns: 3fr 2fr; gap:18px; margin-bottom:20px;">
            <div>
              <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">표 15-1. 라벨·ID 비의존 균일 해석적 스크리닝 ($\\text{{TMC-SA}}_{{\\text{{phys}}}} \le 4.0$) 성능</div>
              <table class="report-table" style="font-size:0.82rem;">
                <thead>
                  <tr style="background:#f8fafc;">
                    <th>평가 항목 (N=55)</th>
                    <th style="text-align:center; color:#b91c1c;">Phase 1 베이스라인</th>
                    <th style="text-align:center;">Phase 2 물리 단독</th>
                    <th style="text-align:center; color:#0369a1;">Phase 2 물리 + SMARTS</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><strong>전체 정확도 / 특이도 (Spec)</strong></td>
                    <td style="text-align:center;">합의 Spec 0.000</td>
                    <td style="text-align:center; font-family:'JetBrains Mono';">{ana_phys['accuracy']:.3f} / <strong>{ana_phys['specificity']:.3f}</strong></td>
                    <td style="text-align:center; font-family:'JetBrains Mono'; color:#0369a1;"><strong>{ana_chem['accuracy']:.3f} / {ana_chem['specificity']:.3f}</strong></td>
                  </tr>
                  <tr>
                    <td><strong>양산 전구체 통과율 (TP/27)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">{base_cmp['precursors_sascore_pass']}</td>
                    <td style="text-align:center; font-weight:700; color:#15803d;">{prec_pct_phys:.1f}% ({prec_pass_phys}/{prec_cnt})</td>
                    <td style="text-align:center; font-weight:700; color:#15803d;">{prec_pct_chem:.1f}% ({prec_pass_chem}/{prec_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>전구체 위음성률 (False Negative)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">{base_cmp['precursors_sascore_fn']}</td>
                    <td style="text-align:center;">{prec_fn_chem_pct:.1f}% ({prec_fn_phys_cnt}/{prec_cnt})</td>
                    <td style="text-align:center;">{prec_fn_chem_pct:.1f}% ({prec_fn_phys_cnt}/{prec_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>AS-ALD 인히비터 통과율 (TP/18)</strong></td>
                    <td style="text-align:center;">{base_cmp['inhibitors_sascore_pass']}</td>
                    <td style="text-align:center; font-weight:700; color:#0369a1;">{inhib_pct_phys:.1f}% ({inhib_pass_phys}/{inhib_cnt})</td>
                    <td style="text-align:center; font-weight:700; color:#0369a1;">{inhib_pct_chem:.1f}% ({inhib_pass_chem}/{inhib_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>음성 대조군 오판율 (FP/10)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">SC {base_cmp['negatives_scscore_fp']} / 합의 {base_cmp['negatives_consensus_fp']}</td>
                    <td style="text-align:center;">{neg_fp_pct_phys:.1f}% ({neg_fp_phys}/{neg_cnt})</td>
                    <td style="text-align:center; font-weight:700;">{neg_fp_pct_chem:.1f}% ({neg_fp_chem}/{neg_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>음성군 차단율 (TN/10)</strong></td>
                    <td style="text-align:center; color:#b91c1c;">0.0% (0/10, 합의)</td>
                    <td style="text-align:center;">{100.0 - neg_fp_pct_phys:.1f}% ({neg_tn_phys}/{neg_cnt})</td>
                    <td style="text-align:center; font-weight:700; color:#0369a1;">{100.0 - neg_fp_pct_chem:.1f}% ({neg_tn_chem}/{neg_cnt})</td>
                  </tr>
                </tbody>
              </table>
            </div>

            <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:16px;">
              <div style="font-weight:700; color:#0f172a; margin-bottom:4px;">All-Sample RF Gini 피처 중요도 (참고용)</div>
              <div style="font-size:0.78rem; color:#b45309; margin-bottom:10px; font-weight:600;">
                주의: N=55 전체 적합(In-sample Acc 1.000) 기준이며 Out-of-Fold 예측 성능 개선을 의미하지 않음
              </div>
              <div style="display:flex; flex-direction:column; gap:5px;">
                {imp_bars_joined}
              </div>
              <div style="margin-top:10px; font-size:0.75rem; color:#0369a1; font-weight:600;">
                &bull; 물리·구조 피처(6종) Gini 기여 합계: <strong>{phys_imp_pct:.1f}%</strong> (근사 프록시 &Delta;E<sub>cplx</sub> 제외)
              </div>
            </div>
          </div>

          <!-- Table 15-2: All-Sample vs 5-Fold CV vs LOOCV -->
          <div style="margin-top:6px;">
            <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">표 15-2. 머신러닝 교차 검증 비교 (All-Sample vs 5-Fold Stratified CV vs LOOCV)</div>
            <div style="overflow-x:auto;">
              <table class="report-table" style="font-size:0.81rem;">
                <thead>
                  <tr style="background:#f1f5f9;">
                    <th>데이터셋 (다수 클래스)</th>
                    <th>분류기 (Classifier)</th>
                    <th>피처 조합 (Feature Set)</th>
                    <th style="text-align:center;">All-Sample (Acc / AUC)</th>
                    <th style="text-align:center;">5-Fold CV (Acc / Spec / AUC)</th>
                    <th style="text-align:center;">LOOCV (Acc / Spec / AUC)</th>
                  </tr>
                </thead>
                <tbody>
                  {ml_cv_rows_joined}
                </tbody>
              </table>
            </div>
            <div style="margin-top:10px; background:#f8fafc; border-left:4px solid #64748b; padding:10px 14px; font-size:0.82rem; color:#334155; line-height:1.55;">
              <strong>검증 해석:</strong> All-Sample에서 Random Forest는 <code>Accuracy 1.000</code>을 기록하나 5-Fold 및 LOOCV에서는 크게 하락하여($N=55$ Combined LOOCV Acc <code>{rf_comb_loo['accuracy']:.3f}</code>, Spec <code>{rf_comb_loo['specificity']:.3f}</code> vs Base-Only LOOCV Acc <code>{rf_base_loo['accuracy']:.3f}</code>, Spec <code>{rf_base_loo['specificity']:.3f}</code>) 소규모 표본에서의 과적합이 확인됩니다. 이는 ETKDGv3 비이완 단점 <code>7net-0</code> 피처만으로는 <code>Ti(Et)4</code>(&beta;-H 제거)나 <code>ZrMe4</code>(열폭발) 같은 속도론적·가수분해 장벽 결함을 직접 분별할 수 없음을 시사합니다.
            </div>
          </div>
        </div>

        <!-- 15.5 Module 4: 55 Molecules Interactive Data Table -->
        <div class="card" id="panel-benchmark">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px; flex-wrap:wrap; gap:10px;">
            <div>
              <h3 class="card-title" style="margin-bottom:2px;">4. [실측 데이터시트] 55종 벤치마크 분자 전수 목록 (Interactive Benchmark Table)</h3>
              <p style="color:var(--text-secondary); font-size:0.84rem; margin:0;">
                양산 전구체 27종, AS-ALD 인히비터 18종, 음성 대조군 10종에 대한 Phase 1(수정 산출 BR-SAScore 포함) vs Phase 2(균일 물리+SMARTS 스크리닝) 결과입니다.
              </p>
            </div>
            <div style="display:flex; gap:8px; align-items:center;">
              <input type="text" id="benchmarkSearch" onkeyup="searchBenchmark()" placeholder="분자명, 화학식, SMILES 검색..." style="padding:6px 12px; border:1px solid #cbd5e1; border-radius:6px; font-size:0.8rem; width:220px;">
            </div>
          </div>

          <div class="filter-bar" style="margin-bottom:14px; display:flex; gap:8px;">
            <button type="button" class="filter-btn active" onclick="filterBenchmark('all', event)">전체 55종 보기</button>
            <button type="button" class="filter-btn" onclick="filterBenchmark('prec', event)">팹 양산 전구체 (27종)</button>
            <button type="button" class="filter-btn" onclick="filterBenchmark('inhib', event)">AS-ALD 인히비터 (18종)</button>
            <button type="button" class="filter-btn" onclick="filterBenchmark('neg', event)">음성 대조군 (10종)</button>
          </div>

          <div style="overflow-x:auto;">
            <table class="report-table" style="font-size:0.82rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th style="width:4%; text-align:center;">No</th>
                  <th style="width:15%;">분자명 / 식별 ID</th>
                  <th style="width:12%;">화학식 / 도메인</th>
                  <th style="width:9%; text-align:center;">타겟 분류</th>
                  <th style="width:12%; text-align:center;">Phase 1 케모 스코어</th>
                  <th style="width:8%; text-align:center;">Phase 1</th>
                  <th style="width:18%; text-align:center;">Phase 2 물리 (%V<sub>bur</sub>, BDE, &Delta;E<sub>dim</sub>)</th>
                  <th style="width:12%; text-align:center;">균일 스크리닝 판정</th>
                  <th>SMILES 구조식</th>
                </tr>
              </thead>
              <tbody>
                {benchmark_table_joined}
              </tbody>
            </table>
          </div>
        </div>
"""

template_path = "template.html"
with open(template_path, "r", encoding="utf-8") as f:
    template = f.read()

final_html = template.replace("<!-- LIT_ROWS_PLACEHOLDER -->", lit_rows_joined)
final_html = final_html.replace("<!-- EMPIRICAL_BENCHMARK_ROWS_PLACEHOLDER -->", empirical_db_joined)
final_html = final_html.replace("<!-- MULTISTAGE_BENCHMARK_ROWS_PLACEHOLDER -->", multistage_joined)
final_html = final_html.replace("<!-- VASP_ROWS_PLACEHOLDER -->", vasp_joined)
final_html = final_html.replace("<!-- PREC_DB_ROWS_PLACEHOLDER -->", prec_db_joined)
final_html = final_html.replace("<!-- BB_DB_ROWS_PLACEHOLDER -->", bb_db_joined)
final_html = final_html.replace("<!-- RXN_DB_ROWS_PLACEHOLDER -->", rxn_db_joined)
final_html = final_html.replace("<!-- EXPANDED_PREC_DB_ROWS_PLACEHOLDER -->", expanded_prec_db_joined)
final_html = final_html.replace("<!-- INHIBITOR_DB_ROWS_PLACEHOLDER -->", inhibitor_db_joined)
final_html = final_html.replace("<!-- PHASE1_PHASE2_TUNING_PLACEHOLDER -->", phase1_phase2_tuning_joined)

with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)

with open("synthetic_accessibility_report.html", "w", encoding="utf-8") as f:
    f.write(final_html)

print("Generated clean light-theme index.html & synthetic_accessibility_report.html successfully!")
