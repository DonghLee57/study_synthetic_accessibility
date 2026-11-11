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
        "desc": "[코드 직접 실측 85.2%] Reaxys 1024bool 신경망 직접 실행: 23종이 1.01~1.94의 매우 낮은 복잡도로 산출(DMAH 1.01, HCDS 1.02, TTIP 1.25, DIPAS 1.73). Reaxys 반응쌍에서 유기금속이 출발물질로 등록되어 전이금속 착체라도 분자량이 작으면 단순 원료로 과소평가(ZMA 2.72, Co(AMD)2 2.38 상승). 자연발화성(Pyrophoric) 등 취급 위험성은 반영 불가."
    },
    "SYBA": {
        "rate": "85.2% (23/27, 실측)",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "코드 직접 실측 (샌드위치 착체 붕괴 / syba.csv.gz)",
        "desc": "[코드 직접 실측 85.2%] 25.7MB 나이브 베이지안 사전(syba.csv.gz; 기본 syba4.csv.gz 대비 방법론적 변경 명시) 직접 연산: 23종은 +3.96~+11.98 양수 정상 판정. 그러나 ZINC 의약 DB에 부재한 메탈로센/Cp 샌드위치 착체인 Ru(EtCp)2(-44.27), CpTi(-24.98), CpZr(-24.98), ZMA(-18.21) 4종에서 극단적 음수(합성 불가) 붕괴 실측 확인."
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


# 16 Multi-stage evaluated models across the 4 execution pipeline stages
multistage_models = [
    {
        "name": "BR-SAScore (엄격 시약사전, 무누수)",
        "gen": "3G 엄격 분리",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("max(existing, -0.47) 보정", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("11.1% (3/27) USPTO 결측", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("11.1% (3/27) 통과 (악화 0종)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (무기 매니폴드 결측; max(existing, -0.47) 보정)",
        "barrier_desc": "BRSAScore는 반경 2 음수 단편(fscore < 0)만 평균내므로, 시약 단편을 max(existing, -0.47)로 갱신하면 점수 악화(0종; CpZr 7.20→6.95, Ru(EtCp)2 4.18→3.99) 없이도 M-N/M-C 결합 중심 단편이 시약에 부재(-6.0)하여 3종(TSA, TMA, DMAH)만 통과함. (단순 0.0 덮어쓰기 시 유기 단편이 분모에서 빠져 4종 악화)"
    },
    {
        "name": "BR-SAScore (생성물 암기 누수 대조군)",
        "gen": "3G 피팅 대조군",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 암기 등록 (s=0)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("100.0% (27/27) 사전 피팅", "background:#fef3c7; color:#b45309; border-color:#fcd34d;"),
        "s4": ("100.0% (27/27) [누수 왜곡]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "감사 검증 (데이터 누수 대조군)",
        "barrier_desc": "Claude 감사에서 규명된 반사실적 대조군: 타깃 전구체 생성물 단편 자체를 0.0으로 사전에 등록하면 금속 배위 단편까지 모두 평균에서 제외(score1=0)되어 27/27종(100.0%, 이전 24종 등록 시 88.9%)이 통과하나, 이는 미지의 전구체에 일반화할 수 없는 순환 기억(Circular Memorization) 왜곡임을 실증."
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
        "barrier_desc": "의약 사전(eMolecules/USPTO) 내 Si-N, M-N 결합 부재로 최대 반응 페널티(s_R=-6.0) 부과되어 24종이 4.18~9.03점으로 과대평가."
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
        "barrier_desc": "중심금속 배위 안정성을 전혀 인지하지 못하고 디메틸아미노/알킬 비트만으로 23종이 ≥0.9점(전체 0.629~0.988)으로 과대평가 (100% False Positive)."
    },
    {
        "name": "RxnFlow",
        "gen": "4G GFlowNet",
        "repo_key": "RxnFlow",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 신톤 그래프 표현", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("33.3% 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 경로 도출", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (72종 신톤 템플릿 풀)",
        "barrier_desc": "Explore(30개)/CGFlow(42개) 총 72개 템플릿에 금속 원자가 전무하여 일반 생성물 부분구조(아민 알킬 사슬 9종) 외에 전구체 코어 합성 불가."
    },
    {
        "name": "BioNeMo ReaSyn",
        "gen": "4G 자기회귀 CoR",
        "repo_key": "ReaSyn",
        "s1": ("100.0% (27/27) (118 원소)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 트랜스포머 시퀀스", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% 부분구조(≥2원자 85.2%) / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 합성 경로 완결", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (SynFormer 115 반응망)",
        "barrier_desc": "118개 원소 어휘 토큰화 및 무손실 역변환(0 UNK)은 100% 완벽하나, 115종 특허 반응 템플릿에는 금속 원자가 전무하고 일반 생성물 부분구조만 매칭되어 무기 착체 배위 단계 생성 0건."
    },
    {
        "name": "AiZynthFinder",
        "gen": "CASP MCTS 엔진",
        "repo_key": "AiZynthFinder",
        "s1": ("100.0% (27/27) Molecule 생성", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% sanitize() 검증 통과", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("미실행 (가중치 미설치 / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (MCTS 미실행, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 (MCTS 가중치 부재 & 정적 분석)",
        "barrier_desc": "aizynthfinder.chem.Molecule(smiles=smi).sanitize() 직접 호출로 Stage 1·2는 100% 통과 검증했으나, USPTO 정책망 가중치 미설치로 실제 MCTS 탐색은 미실행(null 처리). 정적 분석 상 USPTO 50k에는 전이금속 배위 규칙이 부재함."
    },
    {
        "name": "ChemProjector",
        "gen": "4G 분자 투영 디퓨전",
        "repo_key": "ChemProjector",
        "s1": ("100.0% (27/27) 원자번호≤100", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 이중 잠재 공간 투영", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% (23/27) 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 규칙", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (HB-91 템플릿 매칭)",
        "barrier_desc": "HB 91종 SMARTS 매칭 결과 23종에서 일반 생성물 부분구조([N:1][#6:2], [O:1][C:2] 등 2원자 패턴 포함; TSA, HCDS, TMA, DMAH 제외)가 존재하나, 템플릿 파일 내 금속 원자가 전무(0건)하여 전구체 코어 생성 불가."
    },
    {
        "name": "SynNet",
        "gen": "4G 마코프 순방향 계획",
        "repo_key": "SynNet",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Morgan FP (활성비트>0)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% (23/27) 부분구조 / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 규칙", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (HB-91 반응/BB 매칭)",
        "barrier_desc": "2048비트 Morgan 지문 내 활성 비트(GetNumOnBits() > 0) 정상 생성을 27/27 실측 검증. ChemProjector와 동일한 HB 91종 SMARTS 매칭 결과 일반 생성물 부분구조 23/27종은 존재하나, 템플릿 파일 내 금속 원자가 전무(0건)하므로 착체 합성 트리 생성 불가."
    },
    {
        "name": "SynFormer",
        "gen": "4G Seq2Seq 트랜스포머",
        "repo_key": "SynFormer",
        "s1": ("100.0% (27/27) 자체 어휘 통과", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% synformer 토큰 시퀀스", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% (≥2원자 85.2%, ≥4원자 40.7%) / 0% 코어", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) 코어 규칙", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 (115 특허 반응 SMARTS 매칭)",
        "barrier_desc": "synformer 자체 토크나이저(smiles_vocab.txt)는 27/27 통과. comprehensive.txt 115종 생성물 SMARTS 매칭 결과 전체 24/27종(무탄소 TSA는 #84 [N:1] 1원자 패턴), ≥2원자 패턴 23/27종(이 중 12종은 N-C 등 2~3원자 패턴만 매칭), ≥4원자 패턴 11/27종이 매칭되며 금속 원자는 전무(0건)함."
    },
    {
        "name": "Retro*",
        "gen": "신경망 유도 A* 탐색",
        "repo_key": "Retro*",
        "s1": ("100.0% (27/27) rdchiralReactants", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 가치 함수 벡터", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("미실행 (USPTO 템플릿 미설치 / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (A* 미실행, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 (가중치 부재 & 정적 분석)",
        "barrier_desc": "rdchiralReactants(smi) 직접 호출로 객체 초기화(27/27)는 실측 검증했으나, USPTO 50k 템플릿 및 가치 신경망 가중치 미설치로 실제 A* 탐색은 미실행(null 처리). 정적 분석 상 USPTO 유기 반응망에는 금속 배위 규칙이 부재함."
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
# Section 15: Phase 1 & Phase 2 Empirical Benchmark and SevenNet Tuning Generator
# -----------------------------------------------------------------------------
with open("database/benchmark_phase2_tuned_comparison.json", "r", encoding="utf-8") as f:
    bench_p2 = json.load(f)

ana_comp = bench_p2["analytical_comparison"]
ml_comp = bench_p2["machine_learning_comparison"]
molecules_eval = bench_p2["molecules_evaluated"]

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
    
    dimer_color = "#059669" if ddimer >= 0 else "#dc2626"
    tuned_pass = m.get("tuned_pass", False)
    if ttype == "Precursor":
        status_badge = f'<span class="badge" style="background:#dcfce7; color:#15803d; border:1px solid #86efac; font-weight:700;">양산 적합 (TMC-SA {tmc_sa_phys:.1f})</span>' if tuned_pass else f'<span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5; font-weight:700;">위양성 오인 (TMC-SA {tmc_sa_phys:.1f})</span>'
    elif ttype == "Inhibitor":
        status_badge = '<span class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc; font-weight:700;">인히비터 적합 (Pass)</span>' if tuned_pass else '<span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5; font-weight:700;">비선택/불안정</span>'
    else:
        status_badge = '<span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5; font-weight:700;">위험 결함 차단 (Reject)</span>' if not tuned_pass else '<span class="badge" style="background:#fef3c7; color:#b45309; border:1px solid #fcd34d; font-weight:700;">위음성 누락 (FP)</span>'
        
    search_str = f"{mid} {name} {formula} {domain} {ttype} {smi}".lower()
    
    row_html = f"""
            <tr class="benchmark-row" data-type="{data_type}" data-search="{search_str}">
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;">
                <strong>{name}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{mid}</div>
              </td>
              <td>
                <div style="font-weight:600; color:#0f172a; font-size:0.8rem;">{formula}</div>
                <div style="font-size:0.75rem; color:#64748b;">{domain}</div>
              </td>
              <td style="text-align:center; white-space:nowrap;">{type_badge}</td>
              <td style="text-align:center; white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.78rem;">
                <div>SA: <span style="{sa_style}">{sa:.2f}</span></div>
                <div style="color:#475569;">SC: <span style="{sc_style}">{sc:.2f}</span></div>
              </td>
              <td style="text-align:center; white-space:nowrap;">{base_verdict_badge}</td>
              <td style="text-align:center; white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.78rem;">
                <div>%V<sub>bur</sub>: <strong>{vbur:.1f}%</strong></div>
                <div style="color:#475569;">BDE: <strong>{bde:.1f}</strong> kcal &bull; E<sub>coh</sub>: <strong>{ecoh:.2f}</strong> eV</div>
                <div style="font-size:0.72rem; color:{dimer_color};">&Delta;E<sub>dim</sub>: {ddimer:+.1f} kcal &bull; &Delta;E<sub>cplx</sub>: {dcomp:+.1f} kcal</div>
              </td>
              <td style="text-align:center; white-space:nowrap;">{status_badge}</td>
              <td style="font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a; max-width:220px; word-break:break-all;"><code>{smi}</code></td>
            </tr>"""
    benchmark_table_rows.append(row_html)
benchmark_table_joined = "\n".join(benchmark_table_rows)

phase1_phase2_tuning_joined = f"""
        <!-- 15.1 Summary Stat Cards -->
        <div class="hero-stats" style="display:grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap:16px; margin-bottom:24px;">
          <div class="stat-card" style="border-left:4px solid #ef4444; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Phase 1 베이스라인 오류 (Data-Only)</div>
            <div style="font-size:1.85rem; font-weight:800; color:#b91c1c; margin:6px 0;">88.9% &rarr; 100.0%</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              전구체 False Negative <strong>88.9%</strong> (24/27 기각)<br>
              음성군 False Positive <strong>100.0%</strong> (10/10 통과)
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #10b981; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Phase 2 튜닝 후 전구체 구제율</div>
            <div style="font-size:1.85rem; font-weight:800; color:#059669; margin:6px 0;">96.3% (26/27)</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              False Negative <strong>88.9% &rarr; 3.7%</strong> (-85.2%p 혁신)<br>
              TSA, DMAH, TMA, TEMAH 정상 판정 복원
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #6366f1; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">위험/불량 분자 차단율 (FP 제거)</div>
            <div style="font-size:1.85rem; font-weight:800; color:#4338ca; margin:6px 0;">100.0% (10/10)</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              Ti(Et)4 &beta;-수소 탈락, ZrMe4 폭발성,<br>
              WPh5 과밀화(%Vbur 86.4%), Ti(OH)4 중합 전원 차단
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #0284c7; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">4대 물리 피처 기여도 (7net MLIP)</div>
            <div style="font-size:1.85rem; font-weight:800; color:#0369a1; margin:6px 0;">34.6%</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              E_coh(9.4%), BDE(8.3%), &Delta;E_cplx(7.3%), %Vbur(5.8%)<br>
              SambVca 구체 적분 + 7net-0 초고속 물리 스크리닝
            </div>
          </div>
        </div>

        <!-- 15.2 Module 1: The Pathology of Pure Data-Driven Models -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">1. [Phase 1 실측 평가] 순수 데이터 기반 케모인포매틱스 모델의 양방향 병리 (Dual Pathology)</h3>
            <span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5;">Baseline Failure Analysis</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            공개된 오픈소스 케모인포매틱스 모델(RDKit SAScore, Reaxys SCScore, SYBA, RAscore, BR-SAScore)을 물리적 전산 피처 없이 자체 확보한 실물 분자 55종(양산 전구체 27종, AS-ALD 인히비터 18종, 음성 대조군 10종)에 대해 직접 실행한 결과, 다음과 같은 <strong>치명적인 양방향 병리 현상</strong>이 입증되었습니다:
          </p>

          <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-bottom:16px;">
            <div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#991b1b; font-size:0.95rem; margin-bottom:6px;">병리 1: 양산 전구체에 대한 극단적 False Negative (88.9%)</div>
              <p style="font-size:0.84rem; color:#7f1d1d; line-height:1.55; margin:0;">
                RDKit SAScore는 27개 반도체 전구체 중 3개(TEOS, 4DMAS, Ru(EtCp)2)만 통과시키고 <strong>24개를 기각(FN율 88.9%)</strong>했습니다. 글로벌 반도체 팹의 핵심 양산 워크호스인 <code>TSA</code>(SA 6.91), <code>DMAH</code>(SA 6.31), <code>TMA</code>(SA 4.81), <code>TEMAH</code>(SA 4.10)가 심각하게 탈락합니다. 이는 의약품(PubChem/ChEMBL) 통계에서 Si-N, Al-C, Hf-N 결합이 극히 희귀하여 '비천연성 복잡 구조'로 감점되었기 때문입니다.
              </p>
            </div>

            <div style="background:#fffbeb; border:1px solid #fef3c7; border-left:4px solid #f59e0b; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#92400e; font-size:0.95rem; margin-bottom:6px;">병리 2: 위험·불량 분자에 대한 무비판적 False Positive (100.0%)</div>
              <p style="font-size:0.84rem; color:#78350f; line-height:1.55; margin:0;">
                반대로 합성 데이터베이스(Reaxys/USPTO)를 학습한 SCScore와 머신러닝 RAscore는 자발적 열분해 물질인 <code>Ti(Et)4</code>(SC 1.88, RA 0.999), -15 °C에서 폭발하는 <code>ZrMe4</code>(SC 1.76, RA 0.999), 챔버를 부식시키는 <code>트리클로로아세트산</code>(SC 1.15, RA 1.000) 등 <strong>위험 물질 10종 전원(100.0%)을 합성 용이 물질로 통과</strong>시켰습니다. 단순 알킬 사슬과 카르복실산이 유기합성 데이터셋에서 풍부한 시약으로 분류되었기 때문입니다.
              </p>
            </div>
          </div>
        </div>

        <!-- 15.3 Module 2: SevenNet Physical Feature Pipeline -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">2. [Phase 2 물리 피처] SevenNet-0 + SambVca 기반 반도체 전구체 4대 핵심 물리 기둥 (4 Core Physical Pillars)</h3>
            <span class="badge badge-prec" style="font-size:0.82rem;">SevenNet Atomistic Physics</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            DFT 수준의 정확도를 초고속으로 실현하는 <strong>SevenNet (7net-0)</strong> 머신러닝 인터애토믹 포텐셜과 <strong>Cavallo SambVca 구체 적분 모델</strong>을 로컬에 탑재하여, 전구체의 합성과 승화성을 결정짓는 <strong>4대 원자 스케일 핵심 물리 기둥</strong>을 실측 산출하였습니다:
          </p>

          <div style="overflow-x:auto; margin-bottom:16px;">
            <table class="report-table" style="font-size:0.85rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th style="width:20%;">4대 핵심 물리 기둥</th>
                  <th style="width:24%;">수학적 정의 및 지배방정식</th>
                  <th style="width:20%;">반도체 팹 기준 물리적 의미</th>
                  <th>실측 데이터 스케일 및 선별 기준</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td><strong>Pillar 1: 결합 해리 (BDE) &amp; 응집 에너지 ($e_{{\text{{coh}}}}$)</strong></td>
                  <td>$$\text{{BDE}}_{{\min}} = [E(r_0 + 1.2\,\text{{\AA}}) - E(r_0)] \times 23.06$$<br>$$e_{{\text{{coh}}}} = \frac{{\sum E_{{\text{{iso}}}}(Z_i) - E_{{\text{{eq}}}}}}{{N_{{\text{{atoms}}}}}}$$</td>
                  <td>열적 해리 안정성 및 래빌리티 (&beta;-수소 탈락 등 기상 조기 분해 저항성)</td>
                  <td>정상 전구체: <strong>BDE &ge; 50 kcal/mol, $e_{{\text{{coh}}}} \ge 3.8\text{{ eV}}$</strong><br>조기 분해/폭발 위험: <strong>BDE &lt; 45 kcal/mol</strong> (예: CoEt4: 42.9 kcal/mol, 1,2-Dichlorodisilane: 40.9 kcal/mol)</td>
                </tr>
                <tr>
                  <td><strong>Pillar 2: 3D 입체 차폐율 ($\%V_{{\text{{bur}}}}$)</strong></td>
                  <td>$$\%V_{{\text{{bur}}}} = \frac{{V_{{\text{{ligand}}}} \cap V_{{\text{{sphere}}}}(R=3.5\,\text{{\AA}})}}{{V_{{\text{{sphere}}}}}} \times 100$$</td>
                  <td>중심 금속 주변 리간드 입체 보호도 및 올리고머화 방지 차폐막</td>
                  <td>단량체 안정 윈도우: <strong>$60\% \sim 82\%$</strong> (TDMAT: 73.6%, TEMAH: 78.0%)<br>과밀 충돌 분해: <strong>&gt; 85%</strong> (WPh5: 86.4%), 과소배위: <strong>&lt; 50%</strong> (TMA: 41.7%)</td>
                </tr>
                <tr>
                  <td><strong>Pillar 3: 이합체화 상호작용 (&Delta;E_{{\text{{dimer}}}}$)</strong></td>
                  <td>$$\Delta E_{{\text{{dimer}}}} = E(\text{{Dimer}}) - 2E(\text{{Monomer}})$$</td>
                  <td>과소배위에 의한 불휘발성 올리고머/고체 침전 형성 자발성</td>
                  <td>휘발성 단량체 유지: <strong>$\Delta E_{{\text{{dimer}}}} \ge 0\text{{ kcal/mol}}$</strong> (TDMAT: +27.9 kcal/mol)<br>자발적 이합체화/침전: <strong>$\Delta E_{{\text{{dimer}}}} &lt; 0$</strong> (TMA: -2.5 kcal/mol, Ti(OH)4: -0.5 kcal/mol)</td>
                </tr>
                <tr>
                  <td><strong>Pillar 4: 착물화 합성 구동력 (&Delta;E_{{\text{{complex}}}}$)</strong></td>
                  <td>$$M\text{{Cl}}_m + m\,\text{{Li}}L \rightarrow M(L)_m + m\,\text{{LiCl}}\downarrow$$<br>$$\Delta E_{{\text{{complex}}}} = E(M(L)_m) + m E(\text{{LiCl}}) - [E(M\text{{Cl}}_m) + m E(\text{{Li}}L)]$$</td>
                  <td>금속 할로겐화물로부터 전구체 착물 합성(Salt Metathesis)의 열역학적 유리성</td>
                  <td>열역학적 자발 합성: <strong>$\Delta E_{{\text{{complex}}}} &lt; 0\text{{ kcal/mol}}$</strong> (통상 -10 ~ -50 kcal/mol)<br>비자발적 합성 불가: <strong>$\Delta E_{{\text{{complex}}}} &gt; 0\text{{ kcal/mol}}$</strong></td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        <!-- 15.4 Module 3: Model Tuning Comparison Matrix & Feature Importance -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">3. [튜닝 비교 평가] Phase 1 (데이터 전용) vs Phase 2 (7net 4대 물리 기둥 결합) 정량 비교</h3>
            <span class="badge" style="background:#dcfce7; color:#15803d; border:1px solid #86efac;">Before vs After Tuning</span>
          </div>

          <div style="display:grid; grid-template-columns: 3fr 2fr; gap:18px; margin-bottom:18px;">
            <div>
              <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">표 15-1. 합성 가능성 평가 모델 튜닝 전후 성능 비교 매트릭스</div>
              <table class="report-table" style="font-size:0.84rem;">
                <thead>
                  <tr style="background:#f8fafc;">
                    <th>평가 지표 (Evaluation Metric)</th>
                    <th style="text-align:center; color:#b91c1c;">Phase 1 베이스라인 (Data-Only)</th>
                    <th style="text-align:center; color:#15803d;">Phase 2 튜닝 모델 (7net 결합)</th>
                    <th style="text-align:center;">개선폭 (Delta)</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><strong>팹 양산 전구체 통과율 (Precursor Pass)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">11.1% (3/27)</td>
                    <td style="text-align:center; color:#15803d; font-weight:700;">96.3% (26/27)</td>
                    <td style="text-align:center; font-weight:700; color:#059669;">+85.2%p</td>
                  </tr>
                  <tr>
                    <td><strong>전구체 위양성 기각률 (False Negative)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">88.9% (24/27)</td>
                    <td style="text-align:center; color:#15803d; font-weight:700;">3.7% (1/27)</td>
                    <td style="text-align:center; font-weight:700; color:#059669;">-85.2%p</td>
                  </tr>
                  <tr>
                    <td><strong>AS-ALD 인히비터 정상 판정률</strong></td>
                    <td style="text-align:center; color:#475569;">83.3% (15/18)</td>
                    <td style="text-align:center; color:#15803d; font-weight:700;">100.0% (18/18)</td>
                    <td style="text-align:center; font-weight:700; color:#059669;">+16.7%p</td>
                  </tr>
                  <tr>
                    <td><strong>음성 대조군 위음성 통과율 (False Positive)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">100.0% (10/10)</td>
                    <td style="text-align:center; color:#15803d; font-weight:700;">0.0% (0/10)</td>
                    <td style="text-align:center; font-weight:700; color:#059669;">-100.0%p</td>
                  </tr>
                  <tr>
                    <td><strong>위험·결함 분자 차단 성공률</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">0.0% (0/10)</td>
                    <td style="text-align:center; color:#15803d; font-weight:700;">100.0% (10/10)</td>
                    <td style="text-align:center; font-weight:700; color:#059669;">100% 차단</td>
                  </tr>
                </tbody>
              </table>
            </div>

            <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:16px;">
              <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">SevenNet 4대 물리 기둥 기여도 (Random Forest Importance)</div>
              <div style="font-size:0.8rem; color:#64748b; margin-bottom:12px;">5-Fold 교차 검증 분류기에서 도출된 상위 8대 피처 기여율:</div>
              <div style="display:flex; flex-direction:column; gap:6px;">
                <div style="display:flex; justify-content:space-between; font-size:0.8rem;">
                  <span>1. RDKit SAScore</span>
                  <span style="font-family:'JetBrains Mono'; font-weight:700;">20.47%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#64748b; width:20.47%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem; color:#0369a1; font-weight:700;">
                  <span>2. 7net 응집 에너지 (E_coh)</span>
                  <span style="font-family:'JetBrains Mono';">9.40%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#0284c7; width:9.40%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem; color:#0369a1; font-weight:700;">
                  <span>3. 7net 결합 해리 (BDE_min)</span>
                  <span style="font-family:'JetBrains Mono';">8.31%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#0284c7; width:8.31%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem;">
                  <span>4. SYBA</span>
                  <span style="font-family:'JetBrains Mono'; font-weight:700;">7.36%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#64748b; width:7.36%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem; color:#0369a1; font-weight:700;">
                  <span>5. 7net 착물화 구동력 (&Delta;E_cplx)</span>
                  <span style="font-family:'JetBrains Mono';">7.28%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#0284c7; width:7.28%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem;">
                  <span>6. Reaxys SCScore</span>
                  <span style="font-family:'JetBrains Mono'; font-weight:700;">6.59%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#64748b; width:6.59%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem;">
                  <span>7. Heavy Atoms 수</span>
                  <span style="font-family:'JetBrains Mono'; font-weight:700;">6.27%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#64748b; width:6.27%; height:100%;"></div></div>

                <div style="display:flex; justify-content:space-between; font-size:0.8rem; color:#0369a1; font-weight:700;">
                  <span>8. SambVca 입체 차폐율 (%Vbur)</span>
                  <span style="font-family:'JetBrains Mono';">5.79%</span>
                </div>
                <div style="background:#e2e8f0; height:6px; border-radius:3px; overflow:hidden;"><div style="background:#0284c7; width:5.79%; height:100%;"></div></div>
              </div>
              <div style="margin-top:10px; font-size:0.75rem; color:#047857; font-weight:600;">
                &bull; SevenNet 4대 물리 기둥 누적 기여도: <strong>34.6%</strong>
              </div>
            </div>
          </div>
        </div>

        <!-- 15.5 Module 4: 55 Molecules Interactive Data Table -->
        <div class="card" id="panel-benchmark">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px; flex-wrap:wrap; gap:10px;">
            <div>
              <h3 class="card-title" style="margin-bottom:2px;">4. [전수 실측 데이터시트] 55종 실증 벤치마크 분자 전수 데이터시트 (Interactive Benchmark Table)</h3>
              <p style="color:var(--text-secondary); font-size:0.84rem; margin:0;">
                팹 양산 전구체 27종, AS-ALD 인히비터 18종, 음성 대조군 10종에 대한 Phase 1 vs Phase 2 실측 데이터 전수 목록입니다.
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
                  <th style="width:16%;">분자명 / 식별 ID</th>
                  <th style="width:13%;">화학식 / 도메인</th>
                  <th style="width:10%; text-align:center;">타겟 분류</th>
                  <th style="width:11%; text-align:center;">Phase 1 케모 스코어</th>
                  <th style="width:9%; text-align:center;">Phase 1 판정</th>
                  <th style="width:18%; text-align:center;">Phase 2 4대 물리 기둥 (%V<sub>bur</sub>, BDE, &Delta;E<sub>dim</sub>)</th>
                  <th style="width:11%; text-align:center;">튜닝 최종 상태</th>
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
