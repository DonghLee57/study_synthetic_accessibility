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
        "rate": "7.4% (2/27)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "Direct Execution (False Negative Bias)",
        "desc": "[Empirical 7.4%] TEOS(2.84), 4DMAS(3.37) 2종만 threshold 3.5 이하로 통과함(연결형 Ru(EtCp)2 보정 시 4.70). TSA(6.91), DMAH(6.31), TMA(4.81), ZMA(4.82) 등 단순 1단계 양산 전구체에 높은 감점이 부여됨. PubChem 의약 Fragment Frequency 기반으로 Si-N, Al-C, Zr-O 배위를 희귀 구조로 분류한 결과임."
    },
    "SCScore": {
        "rate": "77.8% (21/27)",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "Direct Execution (Early-Reactant Bias)",
        "desc": "[Empirical 77.8%] Reaxys 1024bool DNN 직접 실행 결과 21종이 1.01~1.89의 낮은 점수로 산출됨(DMAH 1.01, HCDS 1.02, TTIP 1.25, DIPAS 1.73). Reaxys 반응쌍에서 유기금속이 출발물질로 등록되어 저분자량 착체가 낮은 complexity로 평가되며(ZMA 2.72, Co(AMD)2 2.38, CCTBA 2.32, Ru(EtCp)2 2.27), 자연발화성(Pyrophoric) 위험은 반영되지 않음."
    },
    "SYBA": {
        "rate": "85.2% (23/27)",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "Direct Execution (Cp Collapse / syba.csv.gz)",
        "desc": "[Empirical 85.2%] 25.7MB Naive Bayes Lookup Table(syba.csv.gz; 기본 syba4.csv.gz 대비 방법론적 변경 명시) 직접 연산 결과 23종은 +2.28~+11.98 양수로 분류됨. 반면 ZINC DB에 부재한 메탈로센/Cp 착체인 CpTi(-24.98), CpZr(-24.98), ZMA(-18.21), Ru(EtCp)2(-17.31) 4종은 큰 음수로 산출됨."
    },
    "BR-SAScore": {
        "rate": "11.1% (Strict Reagent 3/27) / 100% (Leaked)",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "Direct Execution (Unleaked max(existing, -0.47))",
        "desc": "[Empirical] 기본 사전(11.1%, 3/27)과 Data Leakage(생성물 등록 및 CCTBA A>>A 규칙)를 배제한 Strict Reagent Dictionary(음수 fragment 평균 수식 fscore < 0에 맞춰 max(existing, -0.47) 보정) 모두 점수 악화 0종에 Pass Rate 11.1%(3/27: TSA, TMA, DMAH)로 동일함(신규 형성 M-N/M-C 단편 -6.0 유지; 0.0 덮어쓰기 시 유기 단편 제외로 3종 악화). 생성물 사전 등록 시 100.0%(27/27)는 Circular Memorization 결과이므로 Pro-ligand 분리형 2-Track TMC-SA 적용이 필요함."
    },
    "RAscore": {
        "rate": "0.0% Physical Validity (100% FP)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "Direct Execution (False Positive Distortion)",
        "desc": "[Empirical] 27종 전수 0.602~0.988(Mean 0.908)로 수치상 100% Pass 판정이나, ECFP6(radius 3) Fingerprint Count Vector가 중심금속을 일반 알킬/아민 탄소 환경으로 투영한 결과임. 실제 금속 결합 안정성 및 열분해 위험에 대한 물리적 변별력은 0.0%임."
    },
    "RxnFlow": {
        "rate": "0.0% Core (33.3% Substructure)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "Direct Execution (0 Metal Templates)",
        "desc": "[Empirical] Explore(30개)·CGFlow(42개) 총 72개 Template SMARTS 매칭 결과 일반 생성물 부분구조(아민 알킬 사슬 등)를 포함하는 9종(33.3%)만 매칭되었으며, Template 내 금속 원자가 전무하여 전구체 코어 합성 경로 도출률은 0.0%(0/27)임."
    },
    "ReaSyn": {
        "rate": "0.0% Pathway (100% Tokenization)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "Direct Execution (118-Elem Vocab vs Organic Manifold)",
        "desc": "[Empirical] 118원소 Vocabulary 지원으로 SMILES Tokenization 및 무손실 복원은 100%(27/27) 통과하나, SynFormer 115개 반응 Template에는 금속 원자가 없고 일반 유기 부분구조(24/27; ≥2-atom 시 23/27)만 매칭되어 무기 배위 경로는 0건(0.0%)임."
    },
    "AiZynthFinder": {
        "rate": "N/A (MCTS Unexecuted null, 100% Sanitize)",
        "badge_style": "background:#f1f5f9; color:#475569; border-color:#cbd5e1; font-weight:700;",
        "eval": "Object Sanitize + Static Analysis (null)",
        "desc": "[Object Check & Static Analysis] aizynthfinder.chem.Molecule(smiles=smi).sanitize() 직접 호출 시 27종 전구체 객체 생성은 100% 통과함. 외부 USPTO Policy Neural Network Weights 및 Stock HDF5 미설치로 MCTS 탐색은 미실행(null)이며, 정적 분석 상 USPTO 50k에는 유기금속 배위 규칙이 부재함."
    }
}


# 15 Unique Multi-stage evaluated models across the 4 execution pipeline stages (ordered by generation 1G -> 4G -> CASP)
multistage_models = [
    {
        "name": "SAScore (RDKit Built-in)",
        "gen": "1G Fragment Frequency",
        "gen_badge": "badge-1g",
        "repo_key": "SAScore (RDKit Built-in)",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Morgan(r=2)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (Statistical Lookup)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("7.4% (2/27) [FN 92.6%]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — Fragment Frequency Bias & Ring Penalty",
        "barrier_desc": "PubChem 의약 DB 기반 Fragment Frequency 사전에 Si-N, Al-C, M-N 결합이 부재하여 1단계 양산 전구체 TSA(6.91), DMAH(6.31), TMA(4.81) 등 25종(92.6%)에서 False Negative 발생. 에테르/아민 유사 구조 2종(TEOS 2.84, 4DMAS 3.37)만 통과."
    },
    {
        "name": "SYBA",
        "gen": "1G Naive Bayes",
        "gen_badge": "badge-1g",
        "repo_key": "SYBA",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Bayesian Subtrees", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (Log-Odds Ratio)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("85.2% (23/27) [Cp Collapse]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — ZINC Fragment Bias & Metallocene Collapse",
        "barrier_desc": "비환식 알킬아민·알콕사이드 전구체 23종(+2.28~+11.98)은 양수로 분류되나, ZINC 라이브러리(syba.csv.gz) 내 금속-Cp 단편 부재로 CpTi(NMe2)3(-24.98), Cp(NMe2)3Zr(-24.98), ZMA(-18.21), Ru(EtCp)2(-17.31) 4종이 큰 음수로 산출됨."
    },
    {
        "name": "SCScore",
        "gen": "2G Reaction DNN",
        "gen_badge": "badge-2g",
        "repo_key": "SCScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 1024-bit Morgan", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (12-Layer Ranking DNN)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("77.8% (21/27) [Hazard Blind]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 4 — Reaxys Early-Stage Reactant Bias",
        "barrier_desc": "Reaxys 1,200만 반응쌍(full_reaxys_model_1024bool)에서 유기금속이 출발 시약으로 등록된 편향으로 21종이 1.01~1.89로 산출됨. ZMA(2.72)·Co(AMD)2(2.38)·CCTBA(2.32)·Ru(EtCp)2(2.27)의 복잡도는 구분하나 TMA(1.21)·DMAH(1.01)의 자연발화성 및 열분해 위험성은 변별하지 못함."
    },
    {
        "name": "RAscore",
        "gen": "2G MCTS Distillation",
        "gen_badge": "badge-2g",
        "repo_key": "RAscore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("ECFP6 Distortion", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (XGBoost Classifier)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Valid (100% FP)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 — ECFP6 Alkyl Bias & Metal Instability Blindness",
        "barrier_desc": "2048-bit ECFP6 Count Vector가 디메틸아미노·이소프로필 등 유기 치환기 비트에 편향되어 중심금속 배위 안정성을 인지하지 못함. 전수 Mean 0.908(0.602~0.988, 22종 ≥0.90)로 통과하여 불안정 분자를 걸러내지 못하는 100% False Positive 발생."
    },
    {
        "name": "BR-SAScore",
        "gen": "3G Decoupled BFrag/RFrag",
        "gen_badge": "badge-3g",
        "repo_key": "BR-SAScore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% BFrag/RFrag (Unleaked)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("11.1% (3/27) Coordination Gap", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("11.1% (3/27) [Intact Limit]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — Strict Reagent Dictionary Limit (Pro-Ligand Split Required)",
        "barrier_desc": "기본 사전(eMolecules/USPTO)과 Data Leakage를 배제한 Strict Reagent Dictionary(max(existing, -0.47) 보정; 음수 fragment 평균 수식 정합으로 악화 0종) 모두 신규 M-N/M-C 단편 부재(-6.0)로 동일하게 3종(TSA, TMA, DMAH)만 통과함(생성물 사전 등록 시 27/27은 Circular Memorization 결과). 유기 Pro-ligand(HL_i) 분리 평가가 요구됨."
    },
    {
        "name": "GASA",
        "gen": "3G GAT Classifier",
        "gen_badge": "badge-3g",
        "repo_key": "GASA",
        "s1": ("74.1% Supported (7 OOV Metals)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "s2": ("All-Zero Silent Mis-encoding", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT Classifier)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Physical Validity", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 1·2 — DGL Absent & Silent All-Zero Vector Loss (Static Audit)",
        "barrier_desc": "로컬 dgl/dgllife 미설치로 추론 미실행이며 Featurization 코드 정적·정량 감사를 수행함. Allowable Atom Set(allowable_set)에 5대 전이금속(Hf, Ta, Ru, W, Mo)이 누락되어 7종(25.9%)에서 예외 대신 encode_unknown=False로 All-Zero Vector가 생성되는 Silent Mis-encoding 결함 확인."
    },
    {
        "name": "DeepSA",
        "gen": "3G Multimodal Transformer",
        "gen_badge": "badge-3g",
        "repo_key": "DeepSA",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("Ro5(26/27) / Subword Split", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (Multimodal Classifier)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Physical Validity", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — AutoGluon Absent & ELECTRA Bracket Splitting (Static Audit)",
        "barrier_desc": "로컬 autogluon 미설치로 추론 미실행이며 Descriptor·Tokenizer 정적 감사를 수행함. RDKit 실측 상 26/27종이 Lipinski Ro5를 통과(Hf(OtBu)4만 cLogP=5.06 초과)하나, ELECTRA 기반 SMILES Subword Splitting 시 [Zr]·[Hf]·[SiH3] 등 무기 괄호 표현이 파편화됨."
    },
    {
        "name": "DFRscore",
        "gen": "3G GAT Step Regressor",
        "gen_badge": "badge-3g",
        "repo_key": "DFRscore",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("Element Loss ('ELSE' Collapse)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT Regression)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Physical Grounding", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 — Organic-Only Atom Vocabulary ('ELSE' Collapse)",
        "barrier_desc": "Graph Node 허용 원소가 유기 10원소(C, N, O, F, S, Cl, Br, I, B, P)로 제한되어 27종 전구체 전수(100%)의 규소·금속 원자(Si, Ti, Zr, Hf, Ru, W, Mo, Al, Co 등)가 단일 'ELSE' Token으로 축퇴되어 원소별 변별력을 상실함."
    },
    {
        "name": "SynNet",
        "gen": "4G Markov Forward Plan",
        "gen_badge": "badge-4g",
        "repo_key": "SynNet",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 2048-bit Morgan (>0)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% Substruct / 0% Core", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) Core Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — HB-91 Organic Templates Only (0 Metal Rules)",
        "barrier_desc": "2048-bit Morgan Fingerprint 활성 비트(GetNumOnBits() > 0)를 27/27 전구체에서 검증함. 그러나 Hartenfeller-Button(HB) 91종 SMARTS 매칭 시 일반 유기 부분구조 23/27종만 존재하며 Template 내 금속 원자가 전무(0건)하여 착물 합성 트리 생성 불가."
    },
    {
        "name": "ChemProjector",
        "gen": "4G Dual-Diffusion Projection",
        "gen_badge": "badge-4g",
        "repo_key": "ChemProjector",
        "s1": ("100.0% (27/27) (Z≤100)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Latent Embedding", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("85.2% Substruct / 0% Core", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) Core Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — HB-91 Substructure Only (0 Metal Rules)",
        "barrier_desc": "원자번호 Z≤100 기반 Graph Embedding 수용은 27/27 통과. HB 91종 SMARTS 매칭 결과 23종에서 [N:1][#6:2] 등 일반 유기 치환기 부분구조만 탐지(TSA, HCDS, TMA, DMAH 제외)될 뿐, Template 내 금속 원자가 0개이므로 전구체 코어 Projection 불가."
    },
    {
        "name": "SynFormer",
        "gen": "4G Seq2Seq Transformer",
        "gen_badge": "badge-4g",
        "repo_key": "SynFormer",
        "s1": ("100.0% (27/27) Native Vocab", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% CoR Token Sequence", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% (≥2-atom 85.2%) / 0% Core", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) Core Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 115 Patent Templates Without Coordination Rules",
        "barrier_desc": "자체 smiles_vocab Tokenizer 실측 27/27 전수 통과. comprehensive.txt 115종 생성물 SMARTS 매칭 결과 전체 24/27종(무탄소 TSA는 #84 [N:1] 1원자 패턴), ≥2원자 패턴 23/27종, ≥4원자 패턴 11/27종에서 유기 부분구조만 매칭되며 금속 원자 포함 규칙은 0건임."
    },
    {
        "name": "RxnFlow",
        "gen": "4G GFlowNet",
        "gen_badge": "badge-4g",
        "repo_key": "RxnFlow",
        "s1": ("100.0% (27/27)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Synthon Graph", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("33.3% Substruct / 0% Core", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) Core Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 72 Synthon Templates Lack Metal Elements",
        "barrier_desc": "Explore(30종)·CGFlow(42종) 총 72종 SMARTS 역반응 실측 결과 DIPAS, TEMAH, TEMAZ 등 9종(33.3%)의 아민 알킬 사슬(환원적 아민화 등) 부분구조만 매칭됨. Template 내 금속 원자가 전무하여 코어 착물 생성(0/27)은 불가하며 유기 Pro-ligand 조제에만 유효함."
    },
    {
        "name": "BioNeMo ReaSyn",
        "gen": "4G Autoregressive CoR",
        "gen_badge": "badge-4g",
        "repo_key": "ReaSyn",
        "s1": ("100.0% (118-Elem Lossless)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Bidirectional CoR Sequence", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("88.9% Substruct / 0% Core", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": ("0.0% (0/27) Coordination Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — Universal Vocab vs SynFormer-115 Organic Manifold",
        "barrier_desc": "118개 원소 Vocabulary 기반 tokenize_smiles 실측 결과 27/27 전구체에서 OOV(UNK) 0건(100% 무손실 역변환)을 달성함. 그러나 Reaction Manifold [R0, R114]가 SynFormer 115종 유기 특허 반응에 제한되어 일반 부분구조(24/27) 외 무기 배위 경로 생성은 0건임."
    },
    {
        "name": "AiZynthFinder",
        "gen": "CASP MCTS Engine",
        "gen_badge": "badge-ext",
        "repo_key": "AiZynthFinder",
        "s1": ("100.0% (27/27) Molecule", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% sanitize() Validated", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("Unexecuted (Weights Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (MCTS Unexecuted, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — Policy Weights Absent (null) & USPTO Organic Rules",
        "barrier_desc": "aizynthfinder.chem.Molecule.sanitize() 직접 호출을 통해 27/27 분자 객체 수용을 검증함. 외부 USPTO Policy Neural Network Weights 및 Stock HDF5 미설치로 MCTS 탐색은 미실행(null)이며, 정적 분석 상 USPTO 50k에는 전이금속 배위 규칙이 부재함."
    },
    {
        "name": "Retro*",
        "gen": "CASP Neural A*",
        "gen_badge": "badge-ext",
        "repo_key": "Retro*",
        "s1": ("100.0% (27/27) rdchiral", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Value Input Validated", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("Unexecuted (USPTO Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (A* Unexecuted, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — Templates/Weights Absent (null) & USPTO Organic Rules",
        "barrier_desc": "rdchiralReactants(smiles) 직접 호출을 통해 27/27 전구체 Template 객체 초기화 통과를 검증함. 다만 USPTO 50k 반응 Template 및 Value Neural Network Weights 미설치로 A* Search는 미실행(null)이며, 기존 유기 반응망에는 무기 배위 규칙이 결여되어 있음."
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
    
    rxnflow_badge = f"<span class='badge' style='background:#fef3c7; color:#b45309; font-size:0.75rem; font-weight:600;'>Substructure ({len(rxnflow_rules)})</span>" if rxnflow else "<span class='badge' style='background:#f1f5f9; color:#64748b; font-size:0.75rem;'>0 Core</span>"
    reasyn_badge = f"<span class='badge' style='background:#fee2e2; color:#b91c1c; font-size:0.75rem; font-weight:600;'>0 Route (Token OK)</span>" if reasyn_ok else "<span class='badge badge-1g'>Vocab Error</span>"
    aizynth_badge = f"<span class='badge' style='background:#f1f5f9; color:#475569; font-size:0.75rem; font-weight:600;'>Unexecuted null (Sanitize OK)</span>" if aizynth_ok else "<span class='badge badge-1g'>Sanitize Error</span>"
    
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
fallback_cnt = int(p2_meta.get("fallback_count", sum(1 for m in molecules_eval if m.get("is_fallback_geometry"))))

# Lookup key reference molecules dynamically for prose accuracy
mol_by_id = {m["id"]: m for m in molecules_eval}
_tiet4 = mol_by_id.get("NEG-001", {})
_zrme4 = mol_by_id.get("NEG-002", {})
_silcb = mol_by_id.get("NEG-003", {})
_tca = mol_by_id.get("NEG-006", {})
_tdmat = mol_by_id.get("PREC_TI_01", {})
_temah = mol_by_id.get("PREC_HK_01", {})
_dipas = mol_by_id.get("PREC_SI_01", {})
_tma = mol_by_id.get("PREC_AL_01", {})
_ruetc = mol_by_id.get("PREC_MET_03", {})
_cctba = mol_by_id.get("PREC_MET_05", {})
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
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Random Forest (Default)", "Base-Only (10)", ml_n55["random_forest_default"]["base_only"]),
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Random Forest (Default)", "Phys-Only (6)", ml_n55["random_forest_default"]["phys_only"]),
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Random Forest (Default)", "Combined (16)", ml_n55["random_forest_default"]["combined_base_and_phys"], True),
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Random Forest (Balanced)", "Base-Only (10)", ml_n55["random_forest_balanced"]["base_only"]),
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Random Forest (Balanced)", "Combined (16)", ml_n55["random_forest_balanced"]["combined_base_and_phys"], True),
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Logistic Regression (Balanced)", "Base-Only (10)", ml_n55["logistic_regression_balanced"]["base_only"]),
    _fmt_ml_row("Full Dataset N=55 (Majority 0.818)", "Logistic Regression (Balanced)", "Combined (16)", ml_n55["logistic_regression_balanced"]["combined_base_and_phys"]),
    _fmt_ml_row("Precursor Subset N=35 (Majority 0.771)", "Random Forest (Default)", "Base-Only (10)", ml_n35["random_forest_default"]["base_only"]),
    _fmt_ml_row("Precursor Subset N=35 (Majority 0.771)", "Random Forest (Default)", "Phys-Only (6)", ml_n35["random_forest_default"]["phys_only"]),
    _fmt_ml_row("Precursor Subset N=35 (Majority 0.771)", "Random Forest (Default)", "Combined (16)", ml_n35["random_forest_default"]["combined_base_and_phys"], True),
]
ml_cv_rows_joined = "\n".join(ml_cv_rows)

# Build dynamic Feature Importance bars from All-Sample RF Gini importances (top 8)
_feat_labels = {
    "SAScore": ("RDKit SAScore", False),
    "BR-SAScore": ("BR-SAScore", False),
    "Ligand_SA": ("Ligand SA", False),
    "SCScore": ("Reaxys SCScore", False),
    "SYBA": ("SYBA", False),
    "RAscore": ("RAscore", False),
    "LogP": ("LogP", False),
    "MW": ("MW", False),
    "HeavyAtoms": ("HeavyAtoms", False),
    "RotBonds": ("RotBonds", False),
    "7net_BDE_min": ("7net BDE_min", True),
    "7net_E_coh": ("7net E_coh", True),
    "SambVca_Vbur": ("SambVca %V_bur", True),
    "7net_dE_Dimer": ("7net dE_dimer", True),
    "V_mol": ("V_mol", True),
    "Rg": ("Rg", True),
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
        type_badge = '<span class="badge badge-prec" style="font-size:0.75rem;">Precursor</span>'
        data_type = "prec"
    elif ttype == "Inhibitor":
        type_badge = '<span class="badge badge-2g" style="font-size:0.75rem;">Inhibitor</span>'
        data_type = "inhib"
    else:
        type_badge = '<span class="badge badge-4g" style="font-size:0.75rem; background:#fee2e2; color:#991b1b; border:1px solid #f87171;">Negative Control</span>'
        data_type = "neg"

    sa = m.get("sascore", 0.0)
    sc = m.get("scscore", 0.0)
    br_sa = m.get("br_sascore", 9.99)
    lig_sa = m.get("ligand_sa", sa)
    base_pass = m.get("baseline_consensus_pass", False)

    sa_style = "color:#15803d; font-weight:700;" if sa <= 3.5 else "color:#dc2626; font-weight:700;"
    sc_style = "color:#15803d; font-weight:700;" if sc <= 2.0 else "color:#dc2626; font-weight:700;"
    base_verdict_badge = '<span class="badge" style="background:#dcfce7; color:#166534; border:1px solid #86efac; font-size:0.75rem;">Pass</span>' if base_pass else '<span class="badge" style="background:#fee2e2; color:#991b1b; border:1px solid #fca5a5; font-size:0.75rem;">Fail</span>'

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

    if ttype in ("Precursor", "Inhibitor"):
        status_badge = (
            f'<span class="badge" style="background:#dcfce7; color:#15803d; border:1px solid #86efac; font-weight:700;">TP (SA {tmc_sa_phys:.2f})</span>'
            if chem_pass else
            f'<span class="badge" style="background:#fffbeb; color:#b45309; border:1px solid #fcd34d; font-weight:700;">FN (SA {tmc_sa_phys:.2f})</span>'
        )
    else:
        if not chem_pass:
            reason_tag = "SMARTS" if (phys_pass and chem_flag) else "Physics"
            status_badge = f'<span class="badge" style="background:#dcfce7; color:#166534; border:1px solid #86efac; font-weight:700;">TN ({reason_tag} {tmc_sa_phys:.2f})</span>'
        else:
            status_badge = f'<span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5; font-weight:700;">FP (SA {tmc_sa_phys:.2f})</span>'

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

prec_sa_passers = [b for b in benchmarks if b["sascore"] <= 3.5]
prec_sa_pass_names = ", ".join(f"{b['name']} {b['sascore']:.2f}" for b in prec_sa_passers)
prec_sa_fn_cnt = len(benchmarks) - len(prec_sa_passers)

neg_cons_tn = int(base_cmp.get("negatives_consensus_tn", neg_cnt - int(str(base_cmp["negatives_consensus_fp"]).split("/")[0])))
neg_cons_fp_cnt = neg_cnt - neg_cons_tn
neg_cons_fp_pct = round(neg_cons_fp_cnt / neg_cnt * 100.0, 1)
neg_cons_spec_pct = round(100.0 - neg_cons_fp_pct, 1)

prec_failed_names = [mol_by_id[i]["name"].split(" (")[0] for i in brk_chem["precursors"]["failed_ids"]]
prec_failed_note = (", ".join(prec_failed_names) + f" {len(prec_failed_names)}종 FN") if prec_failed_names else "전구체 0종 FN"
neg_rejected_names = [m["name"].split(" (")[0] for m in molecules_eval if m["target_type"] == "Negative Control" and not m.get("tuned_with_chem_pass", True)]
neg_missed_names = [mol_by_id[i]["name"].split(" (")[0] for i in brk_chem["negative_controls"]["missed_negative_ids"]]
neg_rejected_note = ", ".join(neg_rejected_names[:3]) + (" 등" if len(neg_rejected_names) > 3 else "")

phase2_banner_metrics = (
    f"Precursor Pass Rate {prec_pct_chem:.1f}% [{prec_pass_chem}/{prec_cnt}], "
    f"Negative Control Specificity {100.0 - neg_fp_pct_chem:.1f}% [{neg_tn_chem}/{neg_cnt}]"
)


def _step4_sample_row(mid: str, note: str) -> str:
    m = mol_by_id.get(mid, {})
    ttype = m.get("target_type", "Precursor")
    if ttype == "Precursor":
        tbadge = '<span class="badge badge-prec">Precursor</span>'
    else:
        tbadge = '<span class="badge badge-4g">Negative Control</span>'
    pass_flag = bool(m.get("tuned_with_chem_pass", False))
    tmc_sa = float(m.get("tmc_sa_phys", 0.0))
    if ttype == "Precursor":
        vbadge = f'<span class="badge badge-success">TP ({tmc_sa:.2f})</span>' if pass_flag else f'<span class="badge badge-warning">FN ({tmc_sa:.2f})</span>'
    else:
        vbadge = f'<span class="badge badge-success">TN ({tmc_sa:.2f})</span>' if not pass_flag else f'<span class="badge badge-4g">FP ({tmc_sa:.2f})</span>'
    return f"""
                  <tr>
                    <td><strong>{m.get('name', mid).split(' (')[0]}</strong></td>
                    <td>{tbadge}</td>
                    <td style="text-align:center; font-weight:700;">{float(m.get('vbur_percent', 0.0)):.1f}%</td>
                    <td style="text-align:center; font-weight:700;">{float(m.get('delta_e_dimer_kcal_mol', 0.0)):+.1f} kcal/mol</td>
                    <td style="text-align:center; font-weight:700;">{float(m.get('bde_min_kcal_mol', 0.0)):.1f} kcal/mol</td>
                    <td>{vbadge} {note}</td>
                  </tr>"""


_step4_rows = "\n".join([
    _step4_sample_row("PREC_TI_01", "SevenNet-0 LBFGS 이완 후 단량체 입체 차폐 및 BDE 창 충족"),
    _step4_sample_row("PREC_SI_01", "분리 단편 BDE(55.5 kcal/mol)로 단일 원자 충돌 저평가 해소; Si-H %V_bur < 44% 일괄 Penalty로 FN"),
    _step4_sample_row("PREC_AL_01", f"과소배위(%V_bur {float(_tma.get('vbur_percent', 0.0)):.1f}%) 조건에서의 이합체화 경향 및 낮은 Ligand SA 반영"),
    _step4_sample_row("PREC_MET_03", "연결형 SMILES 보정으로 NaN Fallback 제거 및 LBFGS 이완(%V_bur 36.9% 일괄 Penalty로 FN)"),
    _step4_sample_row("PREC_MET_05", "Co2(CO)6 카보닐 배위 정정(C12H10Co2O6); Co-CO 해리 반영 (T_decomp 160 C)"),
    _step4_sample_row("NEG-003", f"이완 고리 신장 BDE {float(_silcb.get('bde_min_kcal_mol', 0.0)):.1f} kcal/mol 산출로 반방향족 고리 변형 TN 차단"),
    _step4_sample_row("NEG-004", f"비중첩 간격({float(_wph5.get('dimer_init_gap_ang', 2.8)):.2f} A) 적용으로 이합체 충돌 해소; LBFGS 이완 후 %V_bur {float(_wph5.get('vbur_percent', 0.0)):.1f}%로 Threshold 미만(FP)"),
    _step4_sample_row("NEG-008", f"분리 단편 BDE {float(_coet4.get('bde_min_kcal_mol', 0.0)):.1f} kcal/mol 및 Ligand SA 합산으로 TN 차단"),
])

phase2_step4_summary_joined = f"""
            <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap:16px; margin-bottom:20px;">
              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #10b981; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">Precursor False Negative Rate (FNR)</div>
                <div style="font-size:1.8rem; font-weight:800; color:#059669; margin:4px 0;">{base_cmp['precursors_sascore_fn']} &rarr; {prec_fn_chem_pct:.1f}%</div>
                <div style="font-size:0.8rem; color:#475569;">Precursor {prec_pass_chem}/{prec_cnt} Pass ({prec_pct_chem:.1f}%), Inhibitor {inhib_pass_chem}/{inhib_cnt} Pass ({inhib_pct_chem:.1f}%)</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #6366f1; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">Negative Control FPR (Un-Hardcoded)</div>
                <div style="font-size:1.8rem; font-weight:800; color:#4338ca; margin:4px 0;">{neg_cons_fp_pct:.1f}% &rarr; {neg_fp_pct_chem:.1f}%</div>
                <div style="font-size:0.8rem; color:#475569;">Physics Only TN {neg_tn_phys}/{neg_cnt} &rarr; +SMARTS TN {neg_tn_chem}/{neg_cnt} ({neg_fp_chem}종 FP 잔존)</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #0284c7; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">Random Forest Validation: All-Sample vs LOOCV</div>
                <div style="font-size:1.8rem; font-weight:800; color:#0369a1; margin:4px 0;">{rf_comb_in['accuracy']:.3f} &rarr; {rf_comb_loo['accuracy']:.3f}</div>
                <div style="font-size:0.8rem; color:#475569;">Combined LOOCV Accuracy {rf_comb_loo['accuracy']:.3f} (Specificity {rf_comb_loo['specificity']:.3f}) vs Base-Only {rf_base_loo['accuracy']:.3f}</div>
              </div>
            </div>

            <div style="font-weight:700; font-size:1.05rem; color:#0f172a; margin-bottom:10px;">
              대표 분자별 SevenNet-0 LBFGS 이완 실측치 및 Uniform Screening 판정 ({comp_time:.1f}s 실행, Fallback {fallback_cnt}/{len(molecules_eval)}종)
            </div>
            <div style="overflow-x:auto; margin-bottom:20px;">
              <table class="report-table" style="font-size:0.84rem;">
                <thead>
                  <tr style="background:#f1f5f9;">
                    <th>Molecule</th>
                    <th>Target Class</th>
                    <th style="text-align:center;">$\\%V_{{\\text{{bur}}}}$</th>
                    <th style="text-align:center;">$\\Delta E_{{\\text{{dimer}}}}$</th>
                    <th style="text-align:center;">$\\text{{BDE}}_{{\\min}}$</th>
                    <th>Uniform Verdict (TMC-SA<sub>phys</sub> &le; 4.0) &amp; Physical Diagnosis</th>
                  </tr>
                </thead>
                <tbody>
{_step4_rows}
                </tbody>
              </table>
            </div>"""

phase1_phase2_tuning_joined = f"""
        <!-- 15.1 Summary Stat Cards -->
        <div class="hero-stats" style="display:grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap:16px; margin-bottom:24px;">
          <div class="stat-card" style="border-left:4px solid #ef4444; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Phase 1 Baseline Error (Data-Only)</div>
            <div style="font-size:1.75rem; font-weight:800; color:#b91c1c; margin:6px 0;">FN {base_cmp['precursors_sascore_fn']} / FP {neg_cons_fp_pct:.1f}%</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              SAScore Precursor FNR <strong>{base_cmp['precursors_sascore_fn']}</strong> ({base_cmp['precursors_sascore_pass']} Pass)<br>
              SCScore Negative FPR <strong>{base_cmp['negatives_scscore_fp']}</strong>, Consensus <strong>{base_cmp['negatives_consensus_fp']}</strong>
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #10b981; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Phase 2 Uniform Screening (TMC-SA &le; 4.0)</div>
            <div style="font-size:1.75rem; font-weight:800; color:#059669; margin:6px 0;">Precursor {prec_pct_chem:.1f}% ({prec_pass_chem}/{prec_cnt})</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              Label/ID-free Continuous Penalty (FNR {base_cmp['precursors_sascore_fn']} &rarr; <strong>{prec_fn_chem_pct:.1f}%</strong>)<br>
              Inhibitor Pass Rate <strong>{inhib_pass_chem}/{inhib_cnt} ({inhib_pct_chem:.1f}%)</strong>, {prec_failed_note}
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #6366f1; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">Negative Control Specificity (Un-Hardcoded)</div>
            <div style="font-size:1.75rem; font-weight:800; color:#4338ca; margin:6px 0;">{neg_tn_chem}/{neg_cnt} (Spec {100.0 - neg_fp_pct_chem:.1f}%)</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              Physics Only TN <strong>{neg_tn_phys}/{neg_cnt}</strong> (FPR {neg_fp_pct_phys:.1f}%), +SMARTS TN <strong>{neg_tn_chem}/{neg_cnt}</strong> (FPR {neg_fp_pct_chem:.1f}%)<br>
              {neg_rejected_note} 차단; Ground-State 한계로 {neg_fp_chem}종 FP 잔존 ({', '.join(neg_missed_names) if neg_missed_names else '없음'})
            </div>
          </div>

          <div class="stat-card" style="border-left:4px solid #0284c7; padding:18px;">
            <div style="font-size:0.82rem; font-weight:700; color:#64748b; text-transform:uppercase;">ML Cross-Validation: All-Sample vs LOOCV</div>
            <div style="font-size:1.75rem; font-weight:800; color:#0369a1; margin:6px 0;">{rf_comb_in['accuracy']:.3f} &rarr; {rf_comb_loo['accuracy']:.3f}</div>
            <div style="font-size:0.82rem; color:#475569; line-height:1.5;">
              Combined Random Forest (N=55): In-Sample 1.000 vs <strong>LOOCV Accuracy {rf_comb_loo['accuracy']:.3f}</strong> (Specificity {rf_comb_loo['specificity']:.3f})<br>
              Base-Only LOOCV {rf_base_loo['accuracy']:.3f}; Physical Gini Importance Sum <strong>{phys_imp_pct:.1f}%</strong> (Proxy 제외)
            </div>
          </div>
        </div>

        <!-- 15.2 Module 1: Phase 1 Baseline Failure Analysis -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">1. [Phase 1 Benchmark] Data-Driven Base Model의 Dual Misclassification 분석</h3>
            <span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5;">Baseline Error Analysis</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            오픈소스 케모인포매틱스 Base Model(RDKit SAScore, Reaxys SCScore, SYBA, RAscore, BR-SAScore)을 물리 Descriptor 없이 55종 Dataset(전구체 27종, 인히비터 18종, Negative Control 10종)에 대해 로컬 실행한 결과, 다음 두 가지 도메인 불일치 오류가 관측되었습니다.
          </p>

          <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-bottom:16px;">
            <div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#991b1b; font-size:0.95rem; margin-bottom:6px;">Error Type 1: Precursor False Negative ({base_cmp['precursors_sascore_fn']})</div>
              <p style="font-size:0.84rem; color:#7f1d1d; line-height:1.55; margin:0;">
                RDKit SAScore(&le; 3.5)는 27개 반도체 전구체 중 {len(prec_sa_passers)}개({prec_sa_pass_names})만 통과시키고 <strong>{prec_sa_fn_cnt}개를 기각(FNR {base_cmp['precursors_sascore_fn']})</strong>합니다. 주요 양산 전구체인 <code>TSA</code>(SA 6.91), <code>DMAH</code>(SA 6.31), <code>TMA</code>(SA 4.81), <code>TEMAH</code>(SA 4.10)가 모두 탈락하며, 이는 PubChem 의약 통계 내 Si-N, Al-C, Hf-N Fragment Frequency가 매우 낮기 때문입니다.
              </p>
            </div>

            <div style="background:#fffbeb; border:1px solid #fef3c7; border-left:4px solid #f59e0b; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#92400e; font-size:0.95rem; margin-bottom:6px;">Error Type 2: Negative Control False Positive (Consensus {neg_cons_fp_pct:.1f}%, SCScore {base_cmp['negatives_scscore_fp']})</div>
              <p style="font-size:0.84rem; color:#78350f; line-height:1.55; margin:0;">
                Reaxys/USPTO 기반 SCScore(&le; 2.0)와 RAscore(&ge; 0.5)는 열분해성 <code>Ti(Et)4</code>(SC {_tiet4.get('scscore', 1.65):.2f}, RA {_tiet4.get('rascore', 0.972):.3f}), 폭발성 <code>ZrMe4</code>(SC {_zrme4.get('scscore', 1.56):.2f}, RA {_zrme4.get('rascore', 0.956):.3f}), 부식성 <code>트리클로로아세트산</code>(SC {_tca.get('scscore', 1.01):.2f}, RA {_tca.get('rascore', 0.970):.3f}) 등 Negative Control 10종 중 SCScore {base_cmp['negatives_scscore_fp']} 및 Baseline Consensus <strong>{base_cmp['negatives_consensus_fp']}을 Pass로 판정</strong>합니다.
              </p>
            </div>
          </div>
        </div>

        <!-- 15.3 Module 2: SevenNet Physical Feature Pipeline (Upgraded Protocol) -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">2. [Phase 2 Physical Features] SevenNet-0 LBFGS 이완 및 비중첩 이합체·분리 단편 BDE</h3>
            <span class="badge badge-prec" style="font-size:0.82rem;">SevenNet-0 LBFGS + SambVca (Fallback {fallback_cnt}/55)</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            구조 정합성 감사에서 확인된 <code>CCTBA</code>(<code>[Co](=O)(=O)=O</code> &rarr; hexacarbonyl <code>C12H10Co2O6</code>), <code>Ru(EtCp)2</code>(분리형 &rarr; 연결형 <code>C14H18Ru</code>), <code>CpZr</code>/<code>CpTi</code>(화학식 <code>H23</code> 정합), 및 <code>Cp2AlH</code> SMILES를 정정하여 Multi-Fragment 3D 중첩 NaN Fallback을 <strong>{fallback_cnt}/55종(전수 정상 계산)</strong>으로 제거하였습니다. 모든 55종 분자는 <strong>SevenNet (<code>7net-0</code>) ASE LBFGS 구조 이완</strong>(단량체 <code>fmax=0.04 eV/&Aring;</code>), 비중첩 주축 정렬 이합체 간격($\\Delta z_{{\\text{{shift}}}} = z_{{\\text{{span}}}} + 2.8\\,\\text{{\\AA}}$ 후 이합체 이완), 및 분리 단편 해리 BDE(비고리 결합 <code>E(ML_{{n-1}})+E(L)-E(ML_n)</code>; 고리 결합 <code>FixAtoms</code> 외환 이동 신장)로 재연산되었으며 총 실행 시간은 <strong>{comp_time:.1f} s</strong>입니다. (<code>&Delta;E<sub>complex</sub></code>는 Algebraic Proxy이므로 ML Feature Matrix에서 제외됨)
          </p>

          <div style="overflow-x:auto; margin-bottom:16px;">
            <table class="report-table" style="font-size:0.85rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th style="width:20%;">Physical Descriptor</th>
                  <th style="width:24%;">수학적 정의 및 계산 방식</th>
                  <th style="width:22%;">Continuous Penalty 경계 ($C^\\infty$ Fermi-Dirac Sigmoid)</th>
                  <th>실측 범위 및 대표 사례</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td><strong>Descriptor 1: 결합 해리 (BDE) &amp; 응집 에너지 ($e_{{\\text{{coh}}}}$)</strong></td>
                  <td>$$\\text{{BDE}}_{{\\min}} = \\min_b [E_{{\\text{{rel}}}}(ML_{{n-1}}) + E_{{\\text{{rel}}}}(L) - E_{{\\text{{rel}}}}(ML_n)] \\times 23.06$$<br>$$e_{{\\text{{coh}}}} = \\frac{{\\sum E_{{\\text{{iso}}}}(Z_i) - E_{{\\text{{rel}}}}}}{{N_{{\\text{{atoms}}}}}}$$</td>
                  <td>열 해리 안정성: $\\text{{BDE}} &lt; 45\\text{{ kcal/mol}}$ 또는 $e_{{\\text{{coh}}}} &lt; 3.10\\text{{ eV/atom}}$ 부근에서 연속 상승</td>
                  <td>정상 전구체: TEMAH BDE <strong>{_temah.get('bde_min_kcal_mol', 0.0):.1f}</strong>, TDMAT <strong>{_tdmat.get('bde_min_kcal_mol', 0.0):.1f}</strong>, DIPAS <strong>{_dipas.get('bde_min_kcal_mol', 0.0):.1f} kcal/mol</strong><br>고리 변형/저안정: Silacyclobutadiene <strong>{_silcb.get('bde_min_kcal_mol', 0.0):.1f}</strong>, CoEt4 <strong>{_coet4.get('bde_min_kcal_mol', 0.0):.1f}</strong>, CCTBA <strong>{_cctba.get('bde_min_kcal_mol', 0.0):.1f} kcal/mol</strong></td>
                </tr>
                <tr>
                  <td><strong>Descriptor 2: 3D 입체 매몰 부피 ($\\%V_{{\\text{{bur}}}}$)</strong></td>
                  <td>$$\\%V_{{\\text{{bur}}}} = \\frac{{V_{{\\text{{ligand}}}} \\cap V_{{\\text{{sphere}}}}(R=3.5\\,\\text{{\\AA}})}}{{V_{{\\text{{sphere}}}}}} \\times 100$$</td>
                  <td>과밀 충돌 Penalty 중심: $\\%V_{{\\text{{bur}}}} &gt; 83.5\\%$<br>과소배위 조건: $\\%V_{{\\text{{bur}}}} &lt; 44.0\\%$</td>
                  <td>LBFGS 이완 단량체: TDMAT <strong>{_tdmat.get('vbur_percent', 0.0):.1f}%</strong>, TEMAH <strong>{_temah.get('vbur_percent', 0.0):.1f}%</strong>, Ru(EtCp)2 <strong>{_ruetc.get('vbur_percent', 0.0):.1f}%</strong><br>고차폐 사례: WPh5 <strong>{_wph5.get('vbur_percent', 0.0):.1f}%</strong> (TMC-SA {_wph5.get('tmc_sa_phys', 0.0):.2f}), TBTDET <strong>{float(mol_by_id.get('PREC_MET_06', {}).get('vbur_percent', 0.0)):.1f}%</strong></td>
                </tr>
                <tr>
                  <td><strong>Descriptor 3: 비중첩 이합체화 에너지 ($\\Delta E_{{\\text{{dimer}}}}$)</strong></td>
                  <td>$$\\Delta E_{{\\text{{dimer}}}} = [E_{{\\text{{rel}}}}(\\text{{Dimer}}; \\Delta z \\ge z_{{\\text{{span}}}}+2.8\\text{{\\AA}}) - 2E_{{\\text{{rel}}}}(\\text{{Monomer}})] \\times 23.06$$</td>
                  <td>금속/Si 중심에서 $\\%V_{{\\text{{bur}}}} &lt; 44\\%$와 $\\Delta E_{{\\text{{dimer}}}} &lt; -1.5\\text{{ kcal/mol}}$ 동시 만족 시 상승</td>
                  <td>초기 간격 $\\ge 2.8\\text{{\\AA}}$ 보장 (WPh5 gap <strong>{_wph5.get('dimer_init_gap_ang', 0.0):.2f} &Aring;</strong>, $\\Delta E_{{\\text{{dim}}}}$ <strong>{_wph5.get('delta_e_dimer_kcal_mol', 0.0):+.1f} kcal/mol</strong>)<br>주요 사례: TDMAT <strong>{_tdmat.get('delta_e_dimer_kcal_mol', 0.0):+.1f}</strong>, TMA <strong>{_tma.get('delta_e_dimer_kcal_mol', 0.0):+.1f}</strong>, Ti(OH)4 <strong>{_tioh4.get('delta_e_dimer_kcal_mol', 0.0):+.1f} kcal/mol</strong></td>
                </tr>
                <tr>
                  <td><strong>Descriptor 4 (Proxy): 착물화 에너지 근사 ($\\Delta E_{{\\text{{complex}}}}$)</strong></td>
                  <td>$$\\Delta E_{{\\text{{complex}}}}^{{\\text{{proxy}}}} = -25.0 + P_{{\\text{{steric}}}} + P_{{\\text{{dimer}}}} + P_{{\\text{{BDE}}}}$$</td>
                  <td>본 단계에서는 Algebraic Proxy로만 산출(실제 VASP/DFT 미수행); <strong>ML Feature Set에서 제외됨</strong></td>
                  <td>후속 단계에서 $M\\text{{Cl}}_m + m\\,\\text{{Li}}L$ 실제 양자화학 반응 에너지로 대체 필요</td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        <!-- 15.4 Module 3: Analytical & ML Cross-Validation Audit -->
        <div class="card" style="margin-bottom:24px;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">3. [Audit &amp; Re-evaluation] Uniform Analytical Screening 및 ML Cross-Validation (All-Sample vs Stratified 5-Fold CV vs LOOCV)</h3>
            <span class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc;">Audited Un-Hardcoded Evaluation</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.88rem; line-height:1.6; margin-bottom:14px;">
            SMILES 정정, LBFGS 이완·비중첩 이합체·분리 단편 BDE 재계산, Target Label Branching 제거, 개별 ID Hardcoding 제거, 및 <code>BR-SAScore</code> 연산 정상화를 모두 반영하여 재실행한 정량 결과입니다.
          </p>

          <div style="display:grid; grid-template-columns: 3fr 2fr; gap:18px; margin-bottom:20px;">
            <div>
              <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">표 15-1. Label/ID-Free Uniform Analytical Screening ($\\text{{TMC-SA}}_{{\\text{{phys}}}} \\le 4.0$) 성능</div>
              <table class="report-table" style="font-size:0.82rem;">
                <thead>
                  <tr style="background:#f8fafc;">
                    <th>Evaluation Metric (N=55)</th>
                    <th style="text-align:center; color:#b91c1c;">Phase 1 Baseline</th>
                    <th style="text-align:center;">Phase 2 Physics Only</th>
                    <th style="text-align:center; color:#0369a1;">Phase 2 Physics + SMARTS</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><strong>Accuracy / Specificity</strong></td>
                    <td style="text-align:center;">Consensus Spec {base_cmp.get('negatives_consensus_specificity', neg_cons_tn / neg_cnt):.3f}</td>
                    <td style="text-align:center; font-family:'JetBrains Mono';">{ana_phys['accuracy']:.3f} / <strong>{ana_phys['specificity']:.3f}</strong></td>
                    <td style="text-align:center; font-family:'JetBrains Mono'; color:#0369a1;"><strong>{ana_chem['accuracy']:.3f} / {ana_chem['specificity']:.3f}</strong></td>
                  </tr>
                  <tr>
                    <td><strong>Precursor Pass Rate (TP/27)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">{base_cmp['precursors_sascore_pass']}</td>
                    <td style="text-align:center; font-weight:700; color:#15803d;">{prec_pct_phys:.1f}% ({prec_pass_phys}/{prec_cnt})</td>
                    <td style="text-align:center; font-weight:700; color:#15803d;">{prec_pct_chem:.1f}% ({prec_pass_chem}/{prec_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>Precursor FNR (FN/27)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">{base_cmp['precursors_sascore_fn']}</td>
                    <td style="text-align:center;">{prec_fn_chem_pct:.1f}% ({prec_fn_phys_cnt}/{prec_cnt})</td>
                    <td style="text-align:center;">{prec_fn_chem_pct:.1f}% ({prec_fn_phys_cnt}/{prec_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>Inhibitor Pass Rate (TP/18)</strong></td>
                    <td style="text-align:center;">{base_cmp['inhibitors_sascore_pass']}</td>
                    <td style="text-align:center; font-weight:700; color:#0369a1;">{inhib_pct_phys:.1f}% ({inhib_pass_phys}/{inhib_cnt})</td>
                    <td style="text-align:center; font-weight:700; color:#0369a1;">{inhib_pct_chem:.1f}% ({inhib_pass_chem}/{inhib_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>Negative Control FPR (FP/10)</strong></td>
                    <td style="text-align:center; color:#b91c1c; font-weight:700;">SC {base_cmp['negatives_scscore_fp']} / Cons {base_cmp['negatives_consensus_fp']}</td>
                    <td style="text-align:center;">{neg_fp_pct_phys:.1f}% ({neg_fp_phys}/{neg_cnt})</td>
                    <td style="text-align:center; font-weight:700;">{neg_fp_pct_chem:.1f}% ({neg_fp_chem}/{neg_cnt})</td>
                  </tr>
                  <tr>
                    <td><strong>Negative Control Specificity (TN/10)</strong></td>
                    <td style="text-align:center; color:#b91c1c;">{neg_cons_spec_pct:.1f}% ({neg_cons_tn}/{neg_cnt}, Cons)</td>
                    <td style="text-align:center;">{100.0 - neg_fp_pct_phys:.1f}% ({neg_tn_phys}/{neg_cnt})</td>
                    <td style="text-align:center; font-weight:700; color:#0369a1;">{100.0 - neg_fp_pct_chem:.1f}% ({neg_tn_chem}/{neg_cnt})</td>
                  </tr>
                </tbody>
              </table>
            </div>

            <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:16px;">
              <div style="font-weight:700; color:#0f172a; margin-bottom:4px;">All-Sample Random Forest Gini Feature Importance</div>
              <div style="font-size:0.78rem; color:#b45309; margin-bottom:10px; font-weight:600;">
                N=55 In-Sample Fit (Accuracy 1.000) 기준이며 Out-of-Fold 일반화 성능 향상과 무관함
              </div>
              <div style="display:flex; flex-direction:column; gap:5px;">
                {imp_bars_joined}
              </div>
              <div style="margin-top:10px; font-size:0.75rem; color:#0369a1; font-weight:600;">
                &bull; Physical Feature (6종) Gini Importance 합계: <strong>{phys_imp_pct:.1f}%</strong> (Proxy &Delta;E<sub>cplx</sub> 제외)
              </div>
            </div>
          </div>

          <!-- Table 15-2: All-Sample vs 5-Fold CV vs LOOCV -->
          <div style="margin-top:6px;">
            <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">표 15-2. ML Cross-Validation 비교 (All-Sample vs Stratified 5-Fold CV vs LOOCV)</div>
            <div style="overflow-x:auto;">
              <table class="report-table" style="font-size:0.81rem;">
                <thead>
                  <tr style="background:#f1f5f9;">
                    <th>Dataset (Majority Baseline)</th>
                    <th>Classifier</th>
                    <th>Feature Set</th>
                    <th style="text-align:center;">All-Sample (Acc / AUROC)</th>
                    <th style="text-align:center;">Stratified 5-Fold CV (Acc / Spec / AUROC)</th>
                    <th style="text-align:center;">LOOCV (Acc / Spec / AUROC)</th>
                  </tr>
                </thead>
                <tbody>
                  {ml_cv_rows_joined}
                </tbody>
              </table>
            </div>
            <div style="margin-top:10px; background:#f8fafc; border-left:4px solid #64748b; padding:10px 14px; font-size:0.82rem; color:#334155; line-height:1.55;">
              <strong>Validation Analysis:</strong> In-Sample에서 Random Forest는 <code>Accuracy 1.000</code>을 기록하나 Stratified 5-Fold CV 및 LOOCV에서 하락하여($N=55$ Combined LOOCV Accuracy <code>{rf_comb_loo['accuracy']:.3f}</code>, Specificity <code>{rf_comb_loo['specificity']:.3f}</code> vs Base-Only LOOCV Accuracy <code>{rf_base_loo['accuracy']:.3f}</code>, Specificity <code>{rf_base_loo['specificity']:.3f}</code>) 소규모 Sample에서의 Overfitting이 확인됩니다. 다만 LBFGS 이완·비중첩 이합체·분리 단편 BDE 적용 후 Combined LOOCV AUROC는 <strong><code>{rf_comb_loo['auroc']:.3f}</code></strong>로 Base-Only(<code>{rf_base_loo['auroc']:.3f}</code>)를 상회하며, Ground-State <code>7net-0</code> Feature만으로는 <code>Ti(Et)4</code>(&beta;-H 제거)나 <code>ZrMe4</code>(열폭발) 같은 활성화 장벽($\\Delta G^\\ddagger$) 의존 반응 속도론적 결함을 직접 분별할 수 없으므로 원소군별 Guardrail로 병행 운용합니다(<code>tests/test_benchmark_reproducibility.py</code> 12건 단위 테스트로 Phase 1·SevenNet·LOOCV 결정론적 재현성 검증 완료).
            </div>
          </div>
        </div>

        <!-- 15.5 Module 4: 55 Molecules Interactive Data Table -->
        <div class="card" id="panel-benchmark">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:12px; flex-wrap:wrap; gap:10px;">
            <div>
              <h3 class="card-title" style="margin-bottom:2px;">4. [Benchmark Datasheet] 55종 분자 전수 실측 테이블 (Interactive Benchmark Table)</h3>
              <p style="color:var(--text-secondary); font-size:0.84rem; margin:0;">
                양산 전구체 27종, AS-ALD 인히비터 18종, Negative Control 10종에 대한 Phase 1(정정 SMILES &amp; BR-SAScore 포함) vs Phase 2(Relaxed LBFGS Physics + SMARTS Screening) 결과입니다.
              </p>
            </div>
            <div style="display:flex; gap:8px; align-items:center;">
              <input type="text" id="benchmarkSearch" onkeyup="searchBenchmark()" placeholder="분자명, 화학식, SMILES 검색..." style="padding:6px 12px; border:1px solid #cbd5e1; border-radius:6px; font-size:0.8rem; width:220px;">
            </div>
          </div>

          <div class="filter-bar" style="margin-bottom:14px; display:flex; gap:8px;">
            <button type="button" class="filter-btn active" onclick="filterBenchmark('all', event)">전체 55종 보기</button>
            <button type="button" class="filter-btn" onclick="filterBenchmark('prec', event)">Precursor (27종)</button>
            <button type="button" class="filter-btn" onclick="filterBenchmark('inhib', event)">Inhibitor (18종)</button>
            <button type="button" class="filter-btn" onclick="filterBenchmark('neg', event)">Negative Control (10종)</button>
          </div>

          <div style="overflow-x:auto;">
            <table class="report-table" style="font-size:0.82rem;">
              <thead>
                <tr style="background:#f1f5f9;">
                  <th style="width:4%; text-align:center;">No</th>
                  <th style="width:15%;">Molecule / ID</th>
                  <th style="width:12%;">Formula / Domain</th>
                  <th style="width:9%; text-align:center;">Target Class</th>
                  <th style="width:12%; text-align:center;">Phase 1 SA Scores</th>
                  <th style="width:8%; text-align:center;">Phase 1</th>
                  <th style="width:18%; text-align:center;">Phase 2 Physics (%V<sub>bur</sub>, BDE, &Delta;E<sub>dim</sub>)</th>
                  <th style="width:12%; text-align:center;">Uniform Verdict</th>
                  <th>SMILES</th>
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
final_html = final_html.replace("<!-- PHASE2_BANNER_METRICS -->", phase2_banner_metrics)
final_html = final_html.replace("<!-- PHASE2_STEP4_SUMMARY_PLACEHOLDER -->", phase2_step4_summary_joined)
final_html = final_html.replace("<!-- PHASE1_PHASE2_TUNING_PLACEHOLDER -->", phase1_phase2_tuning_joined)

with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)

with open("synthetic_accessibility_report.html", "w", encoding="utf-8") as f:
    f.write(final_html)

print("Generated clean light-theme index.html & synthetic_accessibility_report.html successfully!")
