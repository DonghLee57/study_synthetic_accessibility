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
        "desc": "[코드 직접 실측 11.1%] TEOS(2.84), 4DMAS(3.37), Ru(EtCp)2(3.48) 3종만 3.5 이하 정상 통과. TSA(6.91), DMAH(6.31), TMA(4.81), ZMA(4.82) 등 1단계 단순 양산 전구체에 극단적 페널티 발생. PubChem 의약 분자 빈도 기반으로 Si-N, Al-C, Zr-O 등 반도체 핵심 배위를 비천연성 복잡 구조로 오인함."
    },
    "SCScore": {
        "rate": "85.2% (23/27, 실측)",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "코드 직접 실측 (초기원료 편향 / 위험성 은폐)",
        "desc": "[코드 직접 실측 85.2%] Reaxys 1024uint8 신경망 직접 실행: 23종이 1.12~1.98의 매우 낮은 복잡도로 산출(HCDS 1.12, TTIP 1.20, DIPAS 1.70). Reaxys 반응쌍에서 유기금속이 출발물질로 등록되어 전이금속 착체라도 분자량이 작으면 단순 원료로 과소평가(ZMA 2.72, Co(AMD)2 2.38 상승). 자연발화성(Pyrophoric) 등 취급 위험성은 반영 불가."
    },
    "SYBA": {
        "rate": "85.2% (23/27, 실측)",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "코드 직접 실측 (샌드위치 착체 붕괴)",
        "desc": "[코드 직접 실측 85.2%] 25.7MB 나이브 베이지안 사전 직접 연산: 23종은 +3.96~+11.98 양수 정상 판정. 그러나 ZINC 의약 DB에 부재한 메탈로센/Cp 샌드위치 착체인 Ru(EtCp)2(-44.27), CpTi(-24.98), CpZr(-24.98), ZMA(-18.21) 4종에서 극단적 음수(합성 불가) 붕괴 실측 확인."
    },
    "BR-SAScore": {
        "rate": "11.1% (기본사전) / 88.9% (전용사전)",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "코드 직접 실측 (도메인 사전 교체 성공)",
        "desc": "[코드 직접 실측] 기본 emolecules/uspto 사전 적용 시 Si-N/M-N 결합의 -6.0 감점으로 24종이 8~9점으로 과대평가(통과율 11.1%). 그러나 Gelest 실란+금속염 전용 빌딩블록 사전 적용 시 24/27종(88.9%)이 정상 산출(1.2~2.8점)되어 반도체 전구체 도메인 최적화 유효성 실증 완료."
    },
    "RAscore": {
        "rate": "0.0% 물리적 유효 (수치상 100% FP)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (중심금속 결합 무시 False Positive)",
        "desc": "[코드 직접 실측] 27종 전수 0.629~0.988(평균 0.925)로 수치상 100% 통과 판정이나, 이는 ECFP6 반경 3 지문이 중심금속을 단순 고차 알킬/아민 탄소로 취급하여 약물 합성의 쉬운 알킬화로 착각한 전형적 False Positive임. 실제 금속 결합 안정성/산화 위험에 대한 물리적 판정 유효율은 0.0%임."
    },
    "RxnFlow": {
        "rate": "0.0% (코어 0/27, 리간드 33.3%)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (금속 배위 신톤 부재로 코어 생성 0건)",
        "desc": "[코드 직접 실측] 유기 리간드의 알킬 곁사슬 반응(아세톤 환원적 아민화 등)을 포함하는 9종(33.3%)만 단편 매칭되었을 뿐, 금속-배위원자(M-N, M-O, M-C) 결합 형성 템플릿이 전무하여 반도체 전구체 코어 화합물의 실제 합성 경로 도출 성공률은 0.0%(0/27)임."
    },
    "ReaSyn": {
        "rate": "0.0% (경로 0/27, 구문파싱 100%)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (유기 템플릿 한계로 무기 경로 생성 0건)",
        "desc": "[코드 직접 실측] 118개 원소 어휘 토큰이 지원되어 SMILES 파싱(구문 수용)은 100% 에러 없이 통과하나, SynFormer 115개 특허 반응 매니폴드가 유기 반응에 국한되어 전이금속-리간드 배위 결합을 형성하는 합성 경로는 0건(0.0%) 도출됨."
    },
    "AiZynthFinder": {
        "rate": "0.0% (역합성 0/27, 객체파싱 100%)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (무기 정책 부재로 MCTS 탐색 실패 0건)",
        "desc": "[코드 직접 실측] RDKit 기반 27종 전구체 객체 생성(구조 수용)은 100% 통과하나, USPTO 50k 의약 정책 신경망이 금속-리간드 결합 절단 확률을 0으로 산출하여 MCTS 트리가 전혀 전개되지 못하고 실제 역합성 경로 도출 성공률은 0.0%(0/27)임."
    }
}


# 16 Multi-stage evaluated models across the 4 execution pipeline stages
multistage_models = [
    {
        "name": "BR-SAScore (엄격 시약사전, 무누수)",
        "gen": "3G 엄격 분리",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 정상 분리 (BFrag s_B=0)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("14.8% (4/27) USPTO 결측", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("14.8% (4/27) 통과", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (무기 반응 매니폴드 결측)",
        "barrier_desc": "데이터 누수 방지(생성물 배제) 적용 시, 순수 시약 빌딩블록만으로는 새로 형성된 M-N/M-C 결합의 USPTO 부재로 최대 페널티(s_R=-6.0)를 받아 4종(TSA, TMA, DMAH, CCTBA)만 통과."
    },
    {
        "name": "BR-SAScore (생성물 암기 누수 대조군)",
        "gen": "3G 피팅 대조군",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 암기 등록 (s=0)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% (24/27) 사전 피팅", "background:#fef3c7; color:#b45309; border-color:#fcd34d;"),
        "s4": ("88.9% (24/27) [누수 왜곡]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "감사 검증 (데이터 누수 대조군)",
        "barrier_desc": "Claude 감사에서 규명된 반사실적 대조군: 타깃 생성물 단편을 사전에 피팅할 경우 24종 통과가 나오나, 이는 일반화 불가능한 순환 기억(Circular Memorization)에 불과함을 실증."
    },
    {
        "name": "BR-SAScore (기본 의약사전)",
        "gen": "3G 사전 분리",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 분리 (BFrag/RFrag)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("11.1% (3/27) 의약 편향", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("11.1% (3/27) 왜곡", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (반응 매니폴드)",
        "barrier_desc": "의약 사전(eMolecules/USPTO) 내 Si-N, M-N 결합 부재로 최대 반응 페널티(s_R=-6.0) 부과되어 24종이 8~9점으로 과대평가."
    },
    {
        "name": "SAScore (RDKit Built-in)",
        "gen": "1G 단편 빈도",
        "repo_key": "SAScore (RDKit Built-in)",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Morgan 단편 추출", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (단편 휴리스틱)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("11.1% (3/27) 오판", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 (빈도 라이브러리)",
        "barrier_desc": "PubChem 의약 빈도 기반 음의 로그빈도 페널티로 TSA(6.91), DMAH(6.31) 등 단순 전구체 24종에 극단적 False Negative 발생."
    },
    {
        "name": "SCScore",
        "gen": "2G 반응 신경망",
        "repo_key": "SCScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 1024-bit Morgan", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (End-to-End DNN)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("85.2% (23/27) [위험 은폐]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 4 (초기원료 편향)",
        "barrier_desc": "Reaxys 반응쌍에서 유기금속이 출발물질로 등록되어 저분자 전구체를 1점대로 과소평가하여 TMA 등 자연발화성 위험 완전 은폐."
    },
    {
        "name": "SYBA",
        "gen": "1G 베이지안",
        "repo_key": "SYBA",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 나이브 베이즈 단편", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (단편 확률 분류)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("85.2% (23/27) [착체 붕괴]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 2·4 (ZINC 코퍼스)",
        "barrier_desc": "ZINC 의약 DB에 금속-Cp 배위 단편이 전무하여 Ru(EtCp)2(-44.27), CpTi(-24.98) 등 메탈로센 4종에서 극단적 음수 붕괴 발생."
    },
    {
        "name": "RAscore",
        "gen": "2G MCTS 증류",
        "repo_key": "RAscore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("투영 왜곡 (금속→알킬 탄소)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (XGBoost 분류)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 유효 (수치 100% FP)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 (ECFP6 특성 표현)",
        "barrier_desc": "중심금속 배위 안정성을 전혀 인지하지 못하고 디메틸아미노/알킬 비트만으로 27종 전수 0.9점대 과대평가 (100% False Positive)."
    },
    {
        "name": "RxnFlow",
        "gen": "4G GFlowNet",
        "repo_key": "RxnFlow",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 신톤 그래프 표현", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("0.0% 코어 (리간드만 33.3%)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 경로 도출", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (신톤 템플릿 풀)",
        "barrier_desc": "Explore/CGFlow에 금속-배위 결합 룰이 전무하여 아세톤 환원적 아민화 등 유기 곁사슬 9종 외에 전구체 코어 합성 불가."
    },
    {
        "name": "BioNeMo ReaSyn",
        "gen": "4G 자기회귀 CoR",
        "repo_key": "ReaSyn",
        "s1": ("100.0% (27/27) (118 원소)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 트랜스포머 시퀀스", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("0.0% (0/27) 무기 룰 전무", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 합성 경로 완결", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (SynFormer 반응망)",
        "barrier_desc": "118개 원소 어휘 토큰화 및 무손실 역변환(0 UNK)은 100% 완벽하나, 115종 특허 반응망이 유기 반응에 국한되어 무기 착체 배위 단계 생성 0건 (디코딩 오프셋 버그 규명 완료)."
    },
    {
        "name": "AiZynthFinder",
        "gen": "CASP MCTS 엔진",
        "repo_key": "AiZynthFinder",
        "s1": ("100.0% (27/27) 객체 생성", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% RDKit 분자 그래프", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("0.0% (0/27) 정책 전개 실패", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 역합성 완결", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (MCTS 전개 정책망)",
        "barrier_desc": "USPTO 50k 의약 정책망이 금속-리간드 절단 확률을 0으로 산출하여 MCTS 트리가 깊이 0에서 즉시 중단(탐색 실패)."
    },
    {
        "name": "ChemProjector",
        "gen": "4G 분자 투영 디퓨전",
        "repo_key": "ChemProjector",
        "s1": ("100.0% (27/27) 원자번호≤100", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 이중 잠재 공간 투영", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% (23/27) 리간드 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 경로 완결", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (HB-91 템플릿 실측)",
        "barrier_desc": "HB 91종 SMARTS 실측 결과 23종의 유기 치환기 리간드 매칭은 확인되었으나, 전이금속 배위 결합 룰은 0건으로 전구체 코어 생성 불가."
    },
    {
        "name": "SynNet",
        "gen": "4G 마코프 순방향 계획",
        "repo_key": "SynNet",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Morgan FP 임베딩", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% (23/27) 리간드 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 트리 생성", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (HB-91 반응/BB 실측)",
        "barrier_desc": "4대 순방향 신경망이 유기 빌딩블록만 표집하며, HB 91종 룰 중 금속 코어 형성 룰이 0건이므로 전이금속 착체 트리 생성 불가."
    },
    {
        "name": "SynFormer",
        "gen": "4G Seq2Seq 트랜스포머",
        "repo_key": "SynFormer",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 반응 토큰 시퀀스", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("0.0% (0/27) 무기 반응 0건", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 경로 완결", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (115 특허 반응 토큰)",
        "barrier_desc": "SynFormer 반응 토큰 공간에 무기 착체 리간드 교환 룰이 없어 토큰 자기회귀 생성 중단."
    },
    {
        "name": "Retro*",
        "gen": "신경망 유도 A* 탐색",
        "repo_key": "Retro*",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 가치 함수 벡터", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("0.0% (0/27) rdchiral 매칭 실패", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 역합성 완결", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (USPTO 50k 템플릿)",
        "barrier_desc": "rdchiral 템플릿 추출망이 금속-리간드 결합을 인식하지 못해 1단계 트리 확장 후보가 0개로 A* 탐색 즉시 종료."
    },
    {
        "name": "DFRscore",
        "gen": "3G GAT 단계 예측기",
        "repo_key": "DFRscore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("원소 정보 상실 (전수 'ELSE')", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT 회귀)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 물리적 신뢰도", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 (원소 인코딩 'ELSE')",
        "barrier_desc": "Si, Ti, Zr, Hf, Ru, W, Mo, Al, Co 등 모든 금속이 단일 'ELSE' 원소 토큰으로 뭉개져 원소별 합성 난이도 감지 불가."
    },
    {
        "name": "DeepSA",
        "gen": "3G 멀티모달 트랜스포머",
        "repo_key": "DeepSA",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("의약 Ro5 (26/27 적합) & NLP 왜곡", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (복합 SA 지표)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 물리적 타당성", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 (AutoGluon 부재 & 정적 분석)",
        "barrier_desc": "AutoGluon 부재로 미실행(정적 분석). RDKit 실측 상 Ro5는 26종이 통과(Hf(OtBu)4만 LogP=5.06 위반)하나, ELECTRA 서브워드 분할([Zr], [Hf] 파편화)로 무기 원소 표현 왜곡."
    },
    {
        "name": "GASA",
        "gen": "3G GAT 위치선택성",
        "repo_key": "GASA",
        "s1": ("25.9% (7/27) 미지원 원소", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "s2": ("0벡터 무음 오인코딩 (Silent Mis-encoding)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT 분류)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% 물리적 유효성", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 1·2 (DGL 부재 및 무음 오인코딩)",
        "barrier_desc": "DGL/DGL-Life 부재로 미실행 상태이며, 코드 정적 분석 상 Hf, Ta, Ru, W, Mo(7종)는 allowable_set 부재로 에러 대신 atom_type_one_hot에 의해 0벡터로 무음 오인코딩(Silent Mis-encoding)되어 특성 완전 상실."
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
    
    row_html = f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{m['name']}</strong><br>
                <span class="badge badge-1g" style="font-size:0.7rem;">{m['gen']}</span>
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
        rule_short = rxnflow_rules[0].split('_r')[0] if rxnflow_rules else 'ligand'
        diag_items.append(f"<span style='color:#ea580c;'>RxnFlow: 유기리간드 쇄({rule_short})만 매칭 &bull; 금속배위 부재</span>")
    else:
        diag_items.append("<span style='color:#64748b;'>RxnFlow: 배위 템플릿 전무(0건)</span>")
        
    diag_str = " &bull; ".join(diag_items)
    
    sa_style = "color:#15803d; font-weight:600;" if sa <= 3.5 else ("color:#b91c1c; font-weight:700;" if sa > 4.5 else "color:#475569;")
    br_style = "color:#15803d; font-weight:600;" if br <= 3.5 else ("color:#b91c1c; font-weight:600;" if br > 8.0 else "color:#475569;")
    sc_style = "color:#15803d; font-weight:600;" if sc <= 2.0 else "color:#b91c1c; font-weight:700;"
    syba_style = "color:#15803d; font-weight:600;" if syba >= 0 else "color:#b91c1c; font-weight:700; background:#fee2e2; padding:2px 6px; border-radius:4px;"
    ra_style = "color:#0284c7; font-weight:600;" if (ra is not None and ra >= 0.5) else "color:#b91c1c; font-weight:600;"
    
    rxnflow_badge = f"<span class='badge' style='background:#fef3c7; color:#b45309; font-size:0.75rem; font-weight:600;'>리간드만 ({len(rxnflow_rules)}건)</span>" if rxnflow else "<span class='badge' style='background:#f1f5f9; color:#64748b; font-size:0.75rem;'>코어 0건</span>"
    reasyn_badge = f"<span class='badge' style='background:#fee2e2; color:#b91c1c; font-size:0.75rem; font-weight:600;'>경로 0건 (토큰 통과)</span>" if reasyn_ok else "<span class='badge badge-1g'>어휘 에러</span>"
    aizynth_badge = f"<span class='badge' style='background:#fee2e2; color:#b91c1c; font-size:0.75rem; font-weight:600;'>탐색 실패 (객체 수용)</span>" if aizynth_ok else "<span class='badge badge-1g'>수용 에러</span>"
    
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

with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)

with open("synthetic_accessibility_report.html", "w", encoding="utf-8") as f:
    f.write(final_html)

print("Generated clean light-theme index.html & synthetic_accessibility_report.html successfully!")
