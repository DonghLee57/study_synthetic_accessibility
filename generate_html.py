import json
import os
import re

from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

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
# 18-Model & Expanded Semiconductor DB metrics are loaded dynamically below from
# database/benchmark_expanded_models_db.json.


# Load expanded 18-model & expanded semiconductor DB benchmark results
expanded_db_eval_path = os.path.join("database", "benchmark_expanded_models_db.json")
with open(expanded_db_eval_path, "r", encoding="utf-8") as f:
    expanded_db_eval = json.load(f)
exp_p = expanded_db_eval["precursors_102"]["counts"]
exp_p_dist = expanded_db_eval["precursors_102"]["distributions"]
exp_i = expanded_db_eval["inhibitors_18"]["counts"]
exp_r = expanded_db_eval["reagents_29"]["counts"]
N_PREC = exp_p["total"]  # 102 SMILES-verified precursors (Master 205 중 SMILES 보유 전수)

# 18 Unique Multi-stage evaluated models across the 4 execution pipeline stages (ordered by generation 1G -> 4G)
multistage_models = [
    {
        "name": "SAScore (RDKit Built-in)",
        "gen": "1G Heuristic",
        "gen_badge": "badge-1g",
        "repo_key": "SAScore (RDKit Built-in)",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Morgan(r=2)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (Statistical Lookup)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": (f"{exp_p['sascore_pass']/N_PREC*100:.1f}% ({exp_p['sascore_pass']}/{N_PREC}) [FN {(N_PREC-exp_p['sascore_pass'])/N_PREC*100:.1f}%]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — Fragment Frequency Bias & Ring Penalty",
        "barrier_desc": f"PubChem 의약 DB 기반 Fragment Frequency 사전에서 Si-N, Al-C, M-N, Ir/Rh/Ga/In/As/Ru 배위 결합이 부재하여 102종 반도체 프리커서 중 {N_PREC-exp_p['sascore_pass']}종(Mean {exp_p_dist['sascore_mean']})에서 False Negative가 발생함. 에테르·아민 유사 Si 프리커서 2종(TEOS, 4DMAS)만 통과하며 AS-ALD 인히비터는 {exp_i['sascore_pass']}/{exp_i['total']} 통과함."
    },
    {
        "name": "SYBA",
        "gen": "1G Naive Bayes",
        "gen_badge": "badge-1g",
        "repo_key": "SYBA",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Bayesian Subtrees", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (Log-Odds Ratio)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": (f"{exp_p['syba_pass']/N_PREC*100:.1f}% ({exp_p['syba_pass']}/{N_PREC}) [Cp Collapse]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — ZINC Fragment Bias & Metallocene Collapse",
        "barrier_desc": f"102종 프리커서 중 비환식 알킬·알콕사이드 {exp_p['syba_pass']}종(Mean {exp_p_dist['syba_mean']})은 양수로 분류되나, ZINC 라이브러리(syba.csv.gz) 내 금속-Cp·다핵 착물 단편 부재로 {N_PREC-exp_p['syba_pass']}종이 음수로 붕괴함."
    },
    {
        "name": "SCScore",
        "gen": "2G Reaction DNN",
        "gen_badge": "badge-2g",
        "repo_key": "SCScore",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% 1024-bit Morgan", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("N/A (12-Layer Ranking DNN)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": (f"{exp_p['scscore_pass']/N_PREC*100:.1f}% ({exp_p['scscore_pass']}/{N_PREC}) [Hazard Blind]", "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;"),
        "barrier_stage": "Stage 4 — Reaxys Early-Stage Reactant Bias",
        "barrier_desc": f"Reaxys 1,200만 반응쌍(full_reaxys_model_1024bool)에서 단순 유기금속·할로겐화물이 출발 시약으로 등록된 편향으로 {exp_p['scscore_pass']}/{N_PREC}종(Mean {exp_p_dist['scscore_mean']})이 2.0 이하로 통과함. 알킬·할로겐화물 프리커서의 자연발화성 및 열분해 위험성은 변별하지 못함."
    },
    {
        "name": "RAscore",
        "gen": "2G MCTS Distillation",
        "gen_badge": "badge-2g",
        "repo_key": "RAscore",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("ECFP6 Alkyl Bias", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (XGBoost Classifier)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": (f"{exp_p['rascore_pass']/N_PREC*100:.1f}% ({exp_p['rascore_pass']}/{N_PREC}) [FP Bias]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 — ECFP6 Alkyl Bias & Metal Instability Blindness",
        "barrier_desc": f"2048-bit ECFP6 Count Vector가 알킬·아민 치환기 비트에 편향되어 102종 중 {exp_p['rascore_pass']}종({exp_p['rascore_pass']/N_PREC*100:.1f}%, Mean {exp_p_dist['rascore_mean']}, {exp_p_dist['rascore_min']}~{exp_p_dist['rascore_max']})이 0.5 이상으로 통과하여 불안정 착물을 걸러내지 못함."
    },
    {
        "name": "RetroGNN",
        "gen": "2G MPNN Distillation",
        "gen_badge": "badge-2g",
        "repo_key": "RetroGNN",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Graph / ECFP Input", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("Unexecuted (Weights Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (USPTO-Distilled, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — Weights Absent (null) & USPTO Organic Distillation",
        "barrier_desc": "RDKit 기반 분자 그래프 및 2048-bit ECFP 입력(102/102) 수용은 검증됨. 다만 사전학습 MPNN 가중치 미제공으로 추론은 미실행(null)이며, USPTO 역합성 경로 길이 증류 기반이므로 무기 배위 결합 합성 난이도를 반영하지 못함."
    },
    {
        "name": "MPScore",
        "gen": "2G Domain ML (Cage)",
        "gen_badge": "badge-2g",
        "repo_key": "MPScore (Materials Precursor)",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": (f"0.0% ({exp_p['mpscore_organic_in_domain']}/{N_PREC}) Metal OOD", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("Unexecuted (Zenodo Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (Organic Cage OOD, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 2·3 — Zenodo Weights Absent (null) & 100% Metal OOD",
        "barrier_desc": f"ECFP4 활성 비트(102/102)는 산출되나 유기 다공성 케이지(POC) 프리커서(C,H,N,O,F,S,B,P,Cl,Br,I) 영역에 해당되는 프리커서는 0/102(0.0%)임(AS-ALD 인히비터는 {exp_i['mpscore_organic_in_domain']}/{exp_i['total']} 해당). Zenodo 가중치 미설치로 추론은 미실행(null)임."
    },
    {
        "name": "BR-SAScore",
        "gen": "3G Decoupled BFrag/RFrag",
        "gen_badge": "badge-3g",
        "repo_key": "BR-SAScore",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% BFrag/RFrag (Unleaked)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": (f"{exp_p['br_sascore_strict_pass']/N_PREC*100:.1f}% ({exp_p['br_sascore_strict_pass']}/{N_PREC}) Coordination Gap", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": (f"{exp_p['br_sascore_strict_pass']/N_PREC*100:.1f}% ({exp_p['br_sascore_strict_pass']}/{N_PREC}) [Intact Limit]", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — Strict Reagent Dictionary Limit (Pro-Ligand Split Required)",
        "barrier_desc": f"기본 사전과 Data Leakage를 배제한 Strict Reagent Dictionary(max(existing, -0.47) 보정) 모두 102종 중 {exp_p['br_sascore_strict_pass']}종({exp_p['br_sascore_strict_pass']/N_PREC*100:.1f}%)만 통과함(상용 시약 29종은 {exp_r['br_sascore_strict_pass']}/{exp_r['total']} 통과). 유기 Pro-ligand(HL_i) 분리 평가가 요구됨."
    },
    {
        "name": "GASA",
        "gen": "3G GAT Classifier",
        "gen_badge": "badge-3g",
        "repo_key": "GASA",
        "s1": (f"{exp_p['gasa_supported']/N_PREC*100:.1f}% Supported ({exp_p['gasa_silent_misencoded']} OOV)", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "s2": ("All-Zero Silent Mis-encoding", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT Classifier)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Physical Validity", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 1·2 — DGL Absent & 52.0% Silent All-Zero Vector Loss",
        "barrier_desc": f"로컬 dgl/dgllife 미설치로 Featurization 정적·정량 감사를 수행함. Allowable Atom Set(43원소)에서 Ga, Ir, Rh, Ru, Hf, Ta, W, Mo 8원소가 누락되어 102종 중 {exp_p['gasa_silent_misencoded']}종({exp_p['gasa_silent_misencoded']/N_PREC*100:.1f}%)에서 encode_unknown=False로 All-Zero Vector가 생성됨."
    },
    {
        "name": "DeepSA",
        "gen": "3G Multimodal Transformer",
        "gen_badge": "badge-3g",
        "repo_key": "DeepSA",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": (f"Ro5({exp_p['deepsa_ro5_pass']}/{N_PREC}) / Bracket Split({exp_p['deepsa_bracket_fragmented']})", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (Multimodal Classifier)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Physical Validity", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2·4 — AutoGluon Absent & ELECTRA Bracket Splitting",
        "barrier_desc": f"로컬 autogluon 미설치로 Descriptor·Tokenizer 정적 감사를 수행함. 102종 중 {exp_p['deepsa_ro5_pass']}종이 Lipinski Ro5를 통과하나 전수 {exp_p['deepsa_bracket_fragmented']}/{N_PREC}종(100.0%)이 무기 대괄호 표기([Ga], [Ir], [In] 등)를 포함하여 Subword 파편화가 발생함."
    },
    {
        "name": "DFRscore",
        "gen": "3G GAT Step Regressor",
        "gen_badge": "badge-3g",
        "repo_key": "DFRscore",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": (f"100% Loss ({exp_p['dfrscore_else_collapsed']}/{N_PREC} 'ELSE')", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s3": ("N/A (GAT Regression)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("0.0% Physical Grounding", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 2 — Organic-Only Atom Vocabulary ('ELSE' Collapse)",
        "barrier_desc": f"Graph Node 허용 원소가 유기 10원소(C,N,O,F,S,Cl,Br,I,B,P)로 제한되어 102종 프리커서 전수({exp_p['dfrscore_else_collapsed']}/{N_PREC}, 100.0%)의 중심 원자가 단일 'ELSE' Token으로 축퇴됨(AS-ALD 인히비터만 {exp_i['dfrscore_organic_only']}/{exp_i['total']} 유기 원소)."
    },
    {
        "name": "SynNet",
        "gen": "4G Markov Forward Plan",
        "gen_badge": "badge-4g",
        "repo_key": "SynNet",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": (f"100% Morgan ({exp_p['synnet_morgan_active']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": (f"{exp_p['synnet_hb_substruct']/N_PREC*100:.1f}% Substruct / 0% M-L Rule", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": (f"0.0% (0/{N_PREC}) M-L Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — HB-91 Organic Templates Only (0 Metal Rules)",
        "barrier_desc": f"2048-bit Morgan 활성 비트는 {exp_p['synnet_morgan_active']}/{N_PREC} 전수 검증됨. HB-91 SMARTS 매칭 시 유기 부분구조 {exp_p['synnet_hb_substruct']}/{N_PREC}종({exp_p['synnet_hb_substruct']/N_PREC*100:.1f}%)만 탐지되며 금속–리간드(M–L) 배위 반응 규칙은 0건임."
    },
    {
        "name": "ChemProjector",
        "gen": "4G Dual-Diffusion Projection",
        "gen_badge": "badge-4g",
        "repo_key": "ChemProjector",
        "s1": (f"100.0% ({exp_p['chemprojector_z100_pass']}/{N_PREC}) (Z<=100)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Latent Embedding", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": (f"{exp_p['chemprojector_hb_substruct']/N_PREC*100:.1f}% Substruct / 0% M-L Rule", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": (f"0.0% (0/{N_PREC}) M-L Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — HB-91 Substructure Only (0 Metal Rules)",
        "barrier_desc": f"Z<=100 Graph 수용은 {exp_p['chemprojector_z100_pass']}/{N_PREC} 통과. HB-91 SMARTS 매칭 결과 {exp_p['chemprojector_hb_substruct']}/{N_PREC}종에서 유기 치환기 부분구조만 탐지되며, Template 내 금속 원자가 0개이므로 금속–리간드(M–L) 배위 결합 경로 Projection 불가."
    },
    {
        "name": "SynFormer",
        "gen": "4G Seq2Seq Transformer",
        "gen_badge": "badge-4g",
        "repo_key": "SynFormer",
        "s1": (f"100.0% ({exp_p['synformer_vocab_pass']}/{N_PREC}) Native Vocab", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% CoR Token Sequence", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": (f"{exp_p['synformer_substruct']/N_PREC*100:.1f}% (>={exp_p['synformer_ge2_substruct']} 2-atom) / 0% M-L Rule", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": (f"0.0% (0/{N_PREC}) M-L Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 115 Patent Templates Without Coordination Rules",
        "barrier_desc": f"자체 smiles_vocab Tokenizer 실측 {exp_p['synformer_vocab_pass']}/{N_PREC} 전수 통과. comprehensive.txt 115종 매칭 결과 {exp_p['synformer_substruct']}/{N_PREC}종(>=2원자 {exp_p['synformer_ge2_substruct']}/{N_PREC}종)에서 유기 부분구조만 매칭되며 금속–리간드(M–L) 배위 반응 규칙은 0건임."
    },
    {
        "name": "SynFlowNet",
        "gen": "4G GFlowNet (HB Space)",
        "gen_badge": "badge-4g",
        "repo_key": "SynFlowNet",
        "s1": (f"100.0% ({exp_p['synflownet_rdkit_ingest_pass']}/{N_PREC}) RDKit Mol", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Graph Action Space", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("Unexecuted (Weights Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (HB Organic Rules, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — Weights Absent (null) & HB Organic Action Space",
        "barrier_desc": "RDKit 분자 객체 수용(102/102)은 검증됨. 로컬 저장소 및 체크포인트 미설치로 GFlowNet 샘플링은 미실행(null)이며, HB 기반 유기 빌딩블록 반응 공간만 정의되어 금속 착물 합성은 불가함."
    },
    {
        "name": "RxnFlow",
        "gen": "4G GFlowNet (72 Rules)",
        "gen_badge": "badge-4g",
        "repo_key": "RxnFlow",
        "s1": (f"100.0% ({exp_p['rdkit_valid']}/{N_PREC})", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Synthon Graph", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": (f"{exp_p['rxnflow_substruct']/N_PREC*100:.1f}% ({exp_p['rxnflow_substruct']}/{N_PREC}) / 0% M-L Rule", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": (f"0.0% (0/{N_PREC}) M-L Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — 72 Synthon Templates Lack Metal Elements",
        "barrier_desc": f"Explore(30종)·CGFlow(42종) 총 72종 SMARTS 역반응 실측 결과 102종 중 {exp_p['rxnflow_substruct']}종({exp_p['rxnflow_substruct']/N_PREC*100:.1f}%)의 유기 치환기 부분구조만 매칭됨. Template 내 금속 원자가 전무하여 금속–리간드(M–L) 배위 결합 경로 생성(0/{N_PREC})은 불가함."
    },
    {
        "name": "BioNeMo ReaSyn",
        "gen": "4G Autoregressive CoR",
        "gen_badge": "badge-4g",
        "repo_key": "ReaSyn",
        "s1": (f"100.0% Vocab ({exp_p['reasyn_lossless_roundtrip']}/{N_PREC} Exact)", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Bidirectional CoR Sequence", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": (f"{exp_p['synformer_substruct']/N_PREC*100:.1f}% Substruct / 0% M-L Rule", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5;"),
        "s4": (f"0.0% (0/{N_PREC}) M-L Route", "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;"),
        "barrier_stage": "Stage 3 — Universal Vocab vs SynFormer-115 Organic Manifold",
        "barrier_desc": f"118개 원소 Vocabulary 실측 결과 {exp_p['reasyn_vocab_pass']}/{N_PREC} OOV(UNK) 0건이며 {exp_p['reasyn_lossless_roundtrip']}/{N_PREC}종({exp_p['reasyn_lossless_roundtrip']/N_PREC*100:.1f}%)이 무손실 역변환됨. 다만 SynFormer-115 유기 반응 공간에 제한되어 무기 배위 경로 생성은 0건임."
    },
    {
        "name": "AiZynthFinder",
        "gen": "4G CASP MCTS Engine",
        "gen_badge": "badge-ext",
        "repo_key": "AiZynthFinder",
        "s1": (f"100.0% ({exp_p['aizynth_ingest_pass']}/{N_PREC}) Molecule", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% sanitize() Validated", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("Unexecuted (Weights Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (MCTS Unexecuted, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — Policy Weights Absent (null) & USPTO Organic Rules",
        "barrier_desc": f"aizynthfinder.chem.Molecule.sanitize() 직접 호출을 통해 {exp_p['aizynth_ingest_pass']}/{N_PREC} 분자 객체 수용을 검증함. 외부 USPTO Policy 가중치 및 Stock HDF5 미설치로 MCTS 탐색은 미실행(null)이며, USPTO 50k에는 전이금속 배위 규칙이 부재함."
    },
    {
        "name": "Retro*",
        "gen": "4G CASP Neural A*",
        "gen_badge": "badge-ext",
        "repo_key": "Retro*",
        "s1": (f"100.0% ({exp_p['retro_star_rdchiral_ingest_pass']}/{N_PREC}) rdchiral", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s2": ("100% Value Input Validated", "background:#dcfce7; color:#15803d; border-color:#86efac;"),
        "s3": ("Unexecuted (USPTO Absent / null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1;"),
        "s4": ("N/A (A* Unexecuted, null)", "background:#f1f5f9; color:#64748b; border-color:#cbd5e1; font-weight:700;"),
        "barrier_stage": "Stage 3 — Templates/Weights Absent (null) & USPTO Organic Rules",
        "barrier_desc": f"rdchiralReactants(smiles) 직접 호출을 통해 {exp_p['retro_star_rdchiral_ingest_pass']}/{N_PREC} 프리커서 Template 객체 초기화 통과를 검증함. 다만 USPTO 50k Template 및 Value Neural Network 가중치 미설치로 A* Search는 미실행(null)임."
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
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;">{idx:02d}</td>
              <td style="word-break:keep-all; overflow-wrap:anywhere;">
                <strong style="font-size:0.81rem;">{m['name']}</strong><br>
                <span class="badge {gen_badge}" style="font-size:0.67rem; margin-top:2px;">{m['gen']}</span>
              </td>
              <td style="font-size:0.76rem; overflow-wrap:anywhere;">
                <a href="{r.get('url','#')}" target="_blank" style="font-family:'JetBrains Mono'; font-weight:600; font-size:0.73rem;">{r.get('repo','-')}</a>
                <div style="font-size:0.70rem; color:#64748b; margin-top:2px;">Stars: {r.get('stars','-')} &bull; {r.get('license','-')}</div>
              </td>
              <td style="text-align:center;"><span class="badge" style="{s1_style}">{s1_text}</span></td>
              <td style="text-align:center;"><span class="badge" style="{s2_style}">{s2_text}</span></td>
              <td style="text-align:center;"><span class="badge" style="{s3_style}">{s3_text}</span></td>
              <td style="text-align:center;"><span class="badge" style="{s4_style}">{s4_text}</span></td>
              <td style="font-size:0.77rem; line-height:1.45; word-break:keep-all;">
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
    
    rxnflow_badge = f"<span class='badge' style='background:#fef3c7; color:#b45309; font-weight:600;'>Substruct ({len(rxnflow_rules)})</span>" if rxnflow else "<span class='badge' style='background:#f1f5f9; color:#64748b;'>0 Core</span>"
    reasyn_badge = f"<span class='badge' style='background:#fee2e2; color:#b91c1c; font-weight:600;'>0 Route (Token OK)</span>" if reasyn_ok else "<span class='badge badge-1g'>Vocab Error</span>"
    aizynth_badge = f"<span class='badge' style='background:#f1f5f9; color:#475569; font-weight:600;'>Unexecuted null (Sanitize OK)</span>" if aizynth_ok else "<span class='badge badge-1g'>Sanitize Error</span>"
    
    row = f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.72rem; color:#64748b;">{idx:02d}</td>
              <td style="word-break:keep-all; overflow-wrap:anywhere;">
                <strong style="font-size:0.77rem;">{name}</strong>
                <div style="font-family:'JetBrains Mono'; font-size:0.68rem; color:#64748b;">{pid}</div>
              </td>
              <td style="text-align:center;">
                <span class="badge badge-1g" style="margin-bottom:2px;">{metal}</span>
                <div style="font-family:'JetBrains Mono'; font-size:0.67rem; color:#334155; overflow-wrap:anywhere;">{formula}</div>
              </td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.75rem; {sa_style}">{sa}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.75rem; {sc_style}">{sc}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.74rem; {syba_style}">{syba}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.75rem; {br_style}">{br}</td>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.75rem; {ra_style}">{ra}</td>
              <td style="text-align:center;">{rxnflow_badge}</td>
              <td style="text-align:center;">{reasyn_badge}</td>
              <td style="text-align:center;">{aizynth_badge}</td>
              <td style="font-size:0.74rem; line-height:1.4; word-break:keep-all;">{diag_str}</td>
            </tr>"""
    empirical_rows.append(row)
empirical_db_joined = "\n".join(empirical_rows)


# Build Expanded Precursor Master DB rows (Top 120 + Focus Elements)
with open("database/external/expanded_semiconductor_precursors_master.json", "r", encoding="utf-8") as f:
    expanded_prec_db = json.load(f)

exp_sevenn_path = "database/external/expanded_elements_sevenn_features.json"
if os.path.exists(exp_sevenn_path):
    with open(exp_sevenn_path, "r", encoding="utf-8") as f:
        exp_sevenn_by_id = {r["id"]: r for r in json.load(f).get("features", []) if r.get("id")}
else:
    exp_sevenn_by_id = {}

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
prec_failed_note = (", ".join(prec_failed_names) + f" {len(prec_failed_names)}종 FN") if prec_failed_names else "프리커서 0종 FN"
neg_rejected_names = [m["name"].split(" (")[0] for m in molecules_eval if m["target_type"] == "Negative Control" and not m.get("tuned_with_chem_pass", True)]
neg_missed_names = [mol_by_id[i]["name"].split(" (")[0] for i in brk_chem["negative_controls"]["missed_negative_ids"]]
neg_rejected_note = ", ".join(neg_rejected_names[:3]) + (" 등" if len(neg_rejected_names) > 3 else "")

# Load unified latest 160-molecule analysis & refinement reaction energy check for Step 4
unified_path = os.path.join("database", "benchmark_unified_latest_160.json")
with open(unified_path, "r", encoding="utf-8") as f:
    unified_160 = json.load(f)

refine_path = os.path.join("database", "benchmark_refinement_and_rxn_check.json")
with open(refine_path, "r", encoding="utf-8") as f:
    refine_res = json.load(f)

_ref_v2 = refine_res["items_1_and_4"]["cohort_55_v2"]
_ref_exp75 = refine_res["items_1_and_4"]["expanded_75_topology"]
_ref_rxns = refine_res["item_2_reaction_energy_check"]["reactions"]
_rxn_teos = _ref_rxns[0]["delta_e_kcal_mol"]["SevenNet-0"]
_rxn_tdmat_g = _ref_rxns[1]["delta_e_kcal_mol"]["SevenNet-0"]
_rxn_tdmat_s = _ref_rxns[1]["delta_e_with_solid_licl_kcal_mol"]["SevenNet-0"]
_rxn_temah_g = _ref_rxns[2]["delta_e_kcal_mol"]["SevenNet-0"]
_rxn_temah_s = _ref_rxns[2]["delta_e_with_solid_licl_kcal_mol"]["SevenNet-0"]

_u_ds = unified_160["dataset_summary"]
_u_sum = unified_160["unified_v2_summary"]
_u_lat = unified_160["latent_space_160"]
_u_mech = unified_160["by_failure_mechanism"]
_u_ml = unified_160["ml_cv5_160"]
_u_lofmo = unified_160["lofmo_extrapolation_160"]
_u_records = unified_160["records"]
_u_by_id = {r["id"]: r for r in _u_records}

phase2_banner_metrics = (
    f"전체 102종 프리커서 TPR {_u_sum['precursors_all_tpr_pct']:.1f}% [{_u_sum['precursors_all_tp']}], "
    f"40종 Negative Control TNR {_u_sum['negative_all_tnr_pct']:.1f}% [{_u_sum['negative_all_tn']}] "
    f"(감사 55종 기준: 프리커서 25/27 [92.6%], Neg TNR 10/10 [100.0%])"
)


def _svg_x(vbur: float) -> float:
    return 70.0 + (max(0.0, min(100.0, vbur)) / 100.0) * 690.0


def _svg_y_risk(score: float) -> float:
    clamped = max(0.0, min(10.0, score))
    return 340.0 - (clamped / 10.0) * 310.0


_plot160_points_svg = []
for _r in _u_records:
    _grp = _r["group"]
    _vbur = float(_r["vbur_percent"])
    _vmin = float(_r["v_min_adaptive"])
    _sa_v2 = float(_r["tmc_sa_v2"])
    _crisk = int(_r["coord_risk_count"])
    _y_score = min(9.8, _sa_v2 + 2.5 * _crisk)
    _cx, _cy = _svg_x(_vbur), _svg_y_risk(_y_score)
    _nm = str(_r["name"]).replace('"', '&quot;')
    _vsrc = _r["vbur_source"]
    _bde_str = f" | BDE={_r['bde_min_kcal_mol']:.1f} kcal/mol" if _r["bde_min_kcal_mol"] is not None else ""
    _verdict = "PASS" if _r["unified_pass"] else f"REJECT ({', '.join(_r['gate_reasons']) or 'TMC-SA>4.0'})"
    _tip = f"{_nm} [{_r['id']}] | {_grp} | %V_bur={_vbur:.1f}% (V_min={_vmin:.0f}%, {_vsrc}) | TMC-SA_v2={_sa_v2:.2f} | {_verdict}{_bde_str}"

    if "Audited Precursor" in _grp:
        _plot160_points_svg.append(
            f'<circle cx="{_cx:.1f}" cy="{_cy:.1f}" r="5.0" fill="#059669" fill-opacity="0.86" stroke="#fff" stroke-width="1" '
            f'data-tip="{_tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{_tip}</title></circle>'
        )
    elif "Expanded Precursor" in _grp:
        _plot160_points_svg.append(
            f'<polygon points="{_cx:.1f},{_cy-5.2:.1f} {_cx+5.2:.1f},{_cy:.1f} {_cx:.1f},{_cy+5.2:.1f} {_cx-5.2:.1f},{_cy:.1f}" '
            f'fill="#7c3aed" fill-opacity="0.84" stroke="#fff" stroke-width="0.9" '
            f'data-tip="{_tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{_tip}</title></polygon>'
        )
    elif "Inhibitor" in _grp:
        _plot160_points_svg.append(
            f'<rect x="{_cx-4.2:.1f}" y="{_cy-4.2:.1f}" width="8.4" height="8.4" fill="#0284c7" fill-opacity="0.86" stroke="#fff" stroke-width="0.9" '
            f'data-tip="{_tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{_tip}</title></rect>'
        )
    else:
        _plot160_points_svg.append(
            f'<polygon points="{_cx:.1f},{_cy-5.8:.1f} {_cx-5.2:.1f},{_cy+4.2:.1f} {_cx+5.2:.1f},{_cy+4.2:.1f}" '
            f'fill="#dc2626" fill-opacity="0.92" stroke="#fff" stroke-width="0.9" '
            f'data-tip="{_tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{_tip}</title></polygon>'
        )


def _build_latent_svg_160(points: list, x_key: str, y_key: str, x_label: str, y_label: str) -> str:
    xs = [float(p[x_key]) for p in points]
    ys = [float(p[y_key]) for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xpad = max(0.5, (xmax - xmin) * 0.08)
    ypad = max(0.5, (ymax - ymin) * 0.08)
    xmin -= xpad; xmax += xpad
    ymin -= ypad; ymax += ypad
    def _mx(v: float) -> float:
        return 50.0 + ((v - xmin) / (xmax - xmin)) * 340.0
    def _my(v: float) -> float:
        return 250.0 - ((v - ymin) / (ymax - ymin)) * 220.0
    elems = []
    for p in points:
        cx, cy = _mx(float(p[x_key])), _my(float(p[y_key]))
        grp = p["group"]
        nm = str(p["name"]).replace('"', '&quot;')
        sub_lbl = p["failure_mechanism"] if p["label"] == 0 else p["center_atom"]
        tip = f"{nm} [{p['id']} | {sub_lbl}] | {x_label}={p[x_key]:.2f}, {y_label}={p[y_key]:.2f}"
        if "Audited Precursor" in grp:
            elems.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.3" fill="#059669" fill-opacity="0.86" stroke="#fff" stroke-width="0.9" data-tip="{tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{tip}</title></circle>')
        elif "Expanded Precursor" in grp:
            elems.append(f'<polygon points="{cx:.1f},{cy-4.6:.1f} {cx+4.6:.1f},{cy:.1f} {cx:.1f},{cy+4.6:.1f} {cx-4.6:.1f},{cy:.1f}" fill="#7c3aed" fill-opacity="0.82" stroke="#fff" stroke-width="0.8" data-tip="{tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{tip}</title></polygon>')
        elif "Inhibitor" in grp:
            elems.append(f'<rect x="{cx-3.8:.1f}" y="{cy-3.8:.1f}" width="7.6" height="7.6" fill="#0284c7" fill-opacity="0.86" stroke="#fff" stroke-width="0.8" data-tip="{tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{tip}</title></rect>')
        else:
            elems.append(f'<polygon points="{cx:.1f},{cy-5.2:.1f} {cx-4.7:.1f},{cy+3.8:.1f} {cx+4.7:.1f},{cy+3.8:.1f}" fill="#dc2626" fill-opacity="0.92" stroke="#fff" stroke-width="0.8" data-tip="{tip}" onmousemove="showStep4Tip(evt, this)" onmouseleave="hideStep4Tip()"><title>{tip}</title></polygon>')
    pts_str = "\n                  ".join(elems)
    return f"""<svg viewBox="0 0 410 290" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px;">
                  <line x1="50" y1="250" x2="390" y2="250" stroke="#64748b" stroke-width="1.2"/>
                  <line x1="50" y1="30" x2="50" y2="250" stroke="#64748b" stroke-width="1.2"/>
                  <text x="220" y="278" font-size="11" font-weight="700" text-anchor="middle" fill="#0f172a">{x_label}</text>
                  <text x="18" y="140" font-size="11" font-weight="700" text-anchor="middle" fill="#0f172a" transform="rotate(-90 18 140)">{y_label}</text>
                  {pts_str}
                </svg>"""


_pca_160_svg = _build_latent_svg_160(
    _u_records,
    "pc1",
    "pc2",
    f"PC1 ({_u_lat['pca_explained_variance_ratio'][0]*100:.1f}%)",
    f"PC2 ({_u_lat['pca_explained_variance_ratio'][1]*100:.1f}%)",
)
_tsne_160_svg = _build_latent_svg_160(
    _u_records,
    "tsne1",
    "tsne2",
    "t-SNE Dim 1",
    "t-SNE Dim 2",
)


def _step4_unified_row(mid: str, note: str) -> str:
    r = _u_by_id.get(mid)
    if not r:
        return ""
    grp = r["group"]
    if "Audited Precursor" in grp:
        tbadge = '<span class="badge badge-prec">Audited Precursor</span>'
    elif "Expanded Precursor" in grp:
        tbadge = f'<span class="badge badge-2g">Expanded ({r["center_atom"]})</span>'
    else:
        tbadge = '<span class="badge badge-4g">Negative Control</span>'

    if r["label"] == 1:
        vbadge = f'<span class="badge badge-success">TP (SA_v2 {r["tmc_sa_v2"]:.2f})</span>' if r["unified_pass"] else f'<span class="badge badge-warning">FN (SA_v2 {r["tmc_sa_v2"]:.2f})</span>'
    else:
        if not r["unified_pass"]:
            reasons = ", ".join(r["gate_reasons"]) or "Physics"
            vbadge = f'<span class="badge badge-success">TN ({reasons})</span>'
        else:
            vbadge = f'<span class="badge badge-4g">FP (SA_v2 {r["tmc_sa_v2"]:.2f})</span>'

    bde_txt = f"{r['bde_min_kcal_mol']:.1f} kcal/mol" if r["bde_min_kcal_mol"] is not None else "N/A (ETKDGv3)"
    return f"""
                  <tr>
                    <td><strong>{r['name'].split(' (')[0]}</strong> <span style="font-size:0.74rem; color:#64748b;">({r['id']})</span></td>
                    <td>{tbadge}</td>
                    <td style="text-align:center; font-weight:700;">{r['vbur_percent']:.1f}% <span style="font-size:0.72rem; color:#64748b;">(V<sub>min</sub> {r['v_min_adaptive']:.0f}%)</span></td>
                    <td style="text-align:center;">{bde_txt}</td>
                    <td style="text-align:center;">&pi;={r['pi_hapticity_index']} / Ch={r['chelate_ring_count']}</td>
                    <td>{vbadge} {note}</td>
                  </tr>"""


_step4_rep_rows = "\n".join([
    _step4_unified_row("PREC_SI_01", "Si 중심 적응형 V<sub>min</sub>(34%) 적용으로 기존 FN(4.97)에서 TP(2.26)로 구제"),
    _step4_unified_row("PREC_SI_05", "Si 중심 적응형 V<sub>min</sub>(34%) + 분리 단편 BDE 반영으로 TP(3.48) 구제"),
    _step4_unified_row("PREC_MET_03", "&pi;-햅티시티 적응형 V<sub>min</sub>(34%) 적용으로 기존 FN(6.03)에서 TP(3.71)로 구제"),
    _step4_unified_row("EXT-ALD-099", "고리 &eta;-디엔 위상 인식으로 sp3 알킬 오탐 해소 (100% 정상 통과)"),
    _step4_unified_row("EXT-ALD-100", "Ru(EBCHD) &pi;-햅티시티 위상 게이트 인식 (오탐 제거)"),
    _step4_unified_row("NEG-001", "동종 sp3 티타늄 알킬(&beta;-H 제거) 배위 게이트 TN 차단"),
    _step4_unified_row("NEG-004", "비중첩 이합체 이완 및 과밀 다중-아릴 게이트 TN 차단"),
    _step4_unified_row("NEG-036", "불안정 카복실산(TFA) 부식 게이트 TN 차단"),
    _step4_unified_row("NEG-019", "미탐 잔여 3건 중 1건: 반도체 16원소 외곽 Cr 착물 (주기율표 블록 일반화 필요)"),
    _step4_unified_row("NEG-039", "미탐 잔여 3건 중 2건: Cl<sub>2</sub>SiH-SiHCl<sub>2</sub> (Si-Si 불균등화·SevenNet BDE 미연산)"),
])

_mech_labels = [
    ("unstable_homoleptic_tm_alkyl", "불안정 동종 sp3 전이금속 알킬 (&beta;-H 제거)"),
    ("nonvolatile_poly_aryl", "비휘발성 과포화 다중-아릴 금속 착물"),
    ("polyhydroxide_solgel_condensation", "다중 수산화물 (Sol-Gel 즉시 축합)"),
    ("corrosive_carboxylic_acid", "부식성 카복실산 (프로톤 분해·식각)"),
    ("inert_nonfunctional_alkane", "비활성 무작용기 알케인 (화학흡착 불가)"),
    ("strained_or_pyrophoric_hydride_halide", "고변형 실라사이클·자연발화성 수소화물/할로겐화물"),
]
_mech_rows_html = "\n".join(
    f"<tr><td>{label}</td>"
    f"<td>{_u_mech[k]['total']}</td>"
    f"<td>{_u_mech[k]['sascore_rejected']}/{_u_mech[k]['total']} ({_u_mech[k]['sascore_tnr_pct']:.1f}%)</td>"
    f"<td>{_u_mech[k]['syba_rejected']}/{_u_mech[k]['total']} ({_u_mech[k]['syba_tnr_pct']:.1f}%)</td>"
    f"<td><strong>{_u_mech[k]['unified_rejected']}/{_u_mech[k]['total']} ({_u_mech[k]['unified_tnr_pct']:.1f}%)</strong></td></tr>"
    for k, label in _mech_labels
)

_opt_x_adapt = _svg_x(34.0)
_opt_x1, _opt_x2 = _svg_x(44.0), _svg_x(83.5)
_opt_y_top = _svg_y_risk(10.0)
_opt_y_4 = _svg_y_risk(4.0)
_opt_y_bot = _svg_y_risk(0.0)
_plot160_svg_joined = "\n                ".join(_plot160_points_svg)

phase2_step4_summary_joined = f"""
            <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap:16px; margin-bottom:20px;">
              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #10b981; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">전체 Precursor TPR (102종, TMC-SA<sub>v2</sub>)</div>
                <div style="font-size:1.75rem; font-weight:800; color:#059669; margin:4px 0;">{_u_sum['baseline_sascore_precursors_tpr_pct']:.1f}% &rarr; {_u_sum['precursors_all_tpr_pct']:.1f}%</div>
                <div style="font-size:0.8rem; color:#475569;"><strong>{_u_sum['precursors_all_tp']} Pass</strong> (감사 {_u_sum['audited_precursors_tp']} [{_u_sum['audited_precursors_tpr_pct']:.1f}%] &middot; 확장 {_u_sum['expanded_precursors_tp']} [{_u_sum['expanded_precursors_tpr_pct']:.1f}%]; FNR {_u_sum['precursors_all_fnr_pct']:.1f}%)</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #6366f1; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">6대 실패 메커니즘 Negative Control TNR (40종)</div>
                <div style="font-size:1.75rem; font-weight:800; color:#4338ca; margin:4px 0;">{_u_sum['baseline_sascore_neg_tnr_pct']:.1f}% &rarr; {_u_sum['negative_all_tnr_pct']:.1f}%</div>
                <div style="font-size:0.8rem; color:#475569;"><strong>{_u_sum['negative_all_tn']} 차단</strong> (FPR {_u_sum['negative_all_fpr_pct']:.1f}%; SYBA 10.0% 대비 +82.5%p &middot; 미탐 3종: CrPh<sub>3</sub>, Cl<sub>2</sub>SiH-SiHCl<sub>2</sub>, Cp<sub>2</sub>GaH)</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #0284c7; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">160종 통합 CV-5 vs LOFMO 외삽 TNR</div>
                <div style="font-size:1.75rem; font-weight:800; color:#0369a1; margin:4px 0;">LR {_u_ml['LogisticRegression']['cv5_balanced_acc']:.3f} / Gate {_u_lofmo['Unified_TMC_SA_v2_Gate']['mean_unseen_tnr']:.3f}</div>
                <div style="font-size:0.8rem; color:#475569;">CV-5 Bal.Acc (LR {_u_ml['LogisticRegression']['cv5_balanced_acc']:.3f}, SVC {_u_ml['SVC_RBF']['cv5_balanced_acc']:.3f}) vs LOFMO 미지 메커니즘 (SVC 0.000 &middot; OCSVM {_u_lofmo['OneClassSVM_RBF']['mean_unseen_tnr']:.3f} &middot; 물리 게이트 {_u_lofmo['Unified_TMC_SA_v2_Gate']['mean_unseen_tnr']:.3f})</div>
              </div>

              <div style="background:#f8fafc; border:1px solid #e2e8f0; border-left:4px solid #7c3aed; border-radius:6px; padding:16px;">
                <div style="font-size:0.8rem; font-weight:700; color:#64748b;">3D %V<sub>bur</sub> 계측 &amp; SevenNet-0 &Delta;E<sub>syn</sub></div>
                <div style="font-size:1.75rem; font-weight:800; color:#6d28d9; margin:4px 0;">{_u_ds['total_molecules']}/{_u_ds['total_molecules']} (100%)</div>
                <div style="font-size:0.8rem; color:#475569;">SevenNet-0 LBFGS {_u_ds['sevenn_relaxed_count']}종 + RDKit ETKDGv3 {_u_ds['etkdg_3d_count']}종 &middot; TEOS {_rxn_teos:.2f} / TDMAT {_rxn_tdmat_s:.2f} kcal/mol</div>
              </div>
            </div>

            <!-- View 1: Unified 160-Molecule Evaluation Space (%V_bur vs TMC-SA_v2 + Coord Risk) -->
            <div style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:16px; margin-bottom:20px; position:relative;">
              <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:8px;">
                <div style="font-weight:700; font-size:1.02rem; color:#0f172a;">
                  [Unified Plot 1] 최신 통합 160종 물리-배위 위험 공간 (적응형 %V<sub>bur</sub> 창 vs TMC-SA<sub>v2</sub> + Coord Risk)
                </div>
                <span class="badge badge-2g" style="font-size:0.75rem;">102 Precursors (27+75) &bull; 18 Inhibitors &bull; 40 Negative Controls</span>
              </div>
              <p style="font-size:0.82rem; color:#475569; margin:0 0 10px 0; line-height:1.5;">
                업데이트된 DB 전체 <strong>160종 분자</strong>를 최신 모델(<code>TMC-SA<sub>v2</sub></code> 적응형 하한 $V_{{\min}}=34\\%$ [Si·&pi;-haptic] / $44\\%$ [일반 TM] + 6대 실패 메커니즘 배위 게이트)로 단일 재평가한 결과입니다. 하단 녹색·청록색 허용 창(Composite Risk &le; <code>4.0</code>) 내에 <strong>Precursor {_u_sum['precursors_all_tp']} ({_u_sum['precursors_all_tpr_pct']:.1f}%)</strong>이 수렴하고, 붉은 삼각형 Negative Control은 <strong>{_u_sum['negative_all_tn']} ({_u_sum['negative_all_tnr_pct']:.1f}%)</strong> 상단 거부 영역으로 분리됩니다.
              </p>
              <div id="step4-plot-tooltip" style="display:none; position:fixed; z-index:9999; background:#0f172a; color:#f8fafc; padding:6px 10px; border-radius:6px; font-size:0.76rem; font-family:'JetBrains Mono', monospace; pointer-events:none; box-shadow:0 4px 12px rgba(0,0,0,0.2);"></div>
              <svg viewBox="0 0 820 400" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px;">
                <!-- Rejection zone above Y = 4.0 -->
                <rect x="70" y="{_opt_y_top:.1f}" width="690" height="{_opt_y_4 - _opt_y_top:.1f}" fill="#fee2e2" fill-opacity="0.32" />
                <!-- Adaptive Si / pi-haptic window extension (34% to 44%, Y <= 4.0) -->
                <rect x="{_opt_x_adapt:.1f}" y="{_opt_y_4:.1f}" width="{_opt_x1 - _opt_x_adapt:.1f}" height="{_opt_y_bot - _opt_y_4:.1f}" fill="#ccfbf1" fill-opacity="0.65" />
                <!-- Standard TM Optimal window shading (44% to 83.5%, Y <= 4.0) -->
                <rect x="{_opt_x1:.1f}" y="{_opt_y_4:.1f}" width="{_opt_x2 - _opt_x1:.1f}" height="{_opt_y_bot - _opt_y_4:.1f}" fill="#dcfce7" fill-opacity="0.62" />
                <text x="80" y="{_opt_y_top + 18:.1f}" font-size="11" font-weight="700" fill="#b91c1c">Rejection Zone (TMC-SA_v2 &gt; 4.0 or Coordination/Reactivity Gate Triggered)</text>
                <text x="{_opt_x_adapt + 4:.1f}" y="{_opt_y_bot - 8:.1f}" font-size="10" font-weight="700" fill="#0f766e">Si/pi V_min 34%</text>
                <text x="{_opt_x1 + 8:.1f}" y="{_opt_y_bot - 8:.1f}" font-size="10.5" font-weight="700" fill="#15803d">Pass Window (44~83.5%, Composite Risk &lt;= 4.0)</text>
                <!-- Grid & Axes -->
                <line x1="70" y1="340" x2="760" y2="340" stroke="#64748b" stroke-width="1.5"/>
                <line x1="70" y1="30" x2="70" y2="340" stroke="#64748b" stroke-width="1.5"/>
                <!-- Threshold dashed lines -->
                <line x1="{_opt_x_adapt:.1f}" y1="30" x2="{_opt_x_adapt:.1f}" y2="340" stroke="#0d9488" stroke-width="1.2" stroke-dasharray="4,4"/>
                <line x1="{_opt_x1:.1f}" y1="30" x2="{_opt_x1:.1f}" y2="340" stroke="#059669" stroke-width="1.2" stroke-dasharray="4,4"/>
                <line x1="{_opt_x2:.1f}" y1="30" x2="{_opt_x2:.1f}" y2="340" stroke="#059669" stroke-width="1.2" stroke-dasharray="4,4"/>
                <line x1="70" y1="{_opt_y_4:.1f}" x2="760" y2="{_opt_y_4:.1f}" stroke="#dc2626" stroke-width="1.4" stroke-dasharray="5,4"/>
                <text x="764" y="{_opt_y_4 + 4:.1f}" font-size="10.5" font-weight="700" fill="#dc2626">Threshold 4.0</text>
                <!-- X ticks -->
                <text x="70" y="358" font-size="11" text-anchor="middle" fill="#475569">0%</text>
                <text x="{_svg_x(25):.1f}" y="358" font-size="11" text-anchor="middle" fill="#475569">25%</text>
                <text x="{_opt_x_adapt:.1f}" y="358" font-size="11" text-anchor="middle" fill="#0d9488" font-weight="700">34%</text>
                <text x="{_opt_x1:.1f}" y="358" font-size="11" text-anchor="middle" fill="#059669" font-weight="700">44%</text>
                <text x="{_svg_x(60):.1f}" y="358" font-size="11" text-anchor="middle" fill="#475569">60%</text>
                <text x="{_opt_x2:.1f}" y="358" font-size="11" text-anchor="middle" fill="#059669" font-weight="700">83.5%</text>
                <text x="760" y="358" font-size="11" text-anchor="middle" fill="#475569">100%</text>
                <text x="415" y="382" font-size="12" font-weight="700" text-anchor="middle" fill="#0f172a">3D Percent Buried Volume (%V_bur, R = 3.5 A; 85 SevenNet-0 + 75 ETKDGv3)</text>
                <!-- Y ticks -->
                <text x="62" y="{_svg_y_risk(10)+4:.1f}" font-size="11" text-anchor="end" fill="#475569">10.0</text>
                <text x="62" y="{_svg_y_risk(7)+4:.1f}" font-size="11" text-anchor="end" fill="#475569">7.0</text>
                <text x="62" y="{_svg_y_risk(4)+4:.1f}" font-size="11" text-anchor="end" fill="#dc2626" font-weight="700">4.0</text>
                <text x="62" y="{_svg_y_risk(2)+4:.1f}" font-size="11" text-anchor="end" fill="#475569">2.0</text>
                <text x="62" y="{_svg_y_risk(0)+4:.1f}" font-size="11" text-anchor="end" fill="#475569">0.0</text>
                <text x="20" y="185" font-size="12" font-weight="700" text-anchor="middle" fill="#0f172a" transform="rotate(-90 20 185)">Composite Risk (TMC-SA_v2 + Coord Gate)</text>
                <!-- Data points (All 160 molecules) -->
                {_plot160_svg_joined}
              </svg>
              <div style="display:flex; gap:18px; flex-wrap:wrap; margin-top:10px; font-size:0.8rem; color:#334155;">
                <span><span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:#059669; margin-right:5px;"></span>Audited Precursor (27)</span>
                <span><span style="display:inline-block; width:10px; height:10px; transform:rotate(45deg); background:#7c3aed; margin-right:6px;"></span>Expanded Ir/Rh/As/Ga/In/Ru Precursor (75)</span>
                <span><span style="display:inline-block; width:10px; height:10px; background:#0284c7; margin-right:5px;"></span>AS-ALD Inhibitor (18)</span>
                <span><span style="display:inline-block; width:0; height:0; border-left:6px solid transparent; border-right:6px solid transparent; border-bottom:10px solid #dc2626; margin-right:5px;"></span>Negative Control across 6 Failure Classes (40)</span>
              </div>
              <script>
                function showStep4Tip(evt, el) {{
                  var tip = document.getElementById('step4-plot-tooltip');
                  if (!tip) return;
                  tip.textContent = el.getAttribute('data-tip');
                  tip.style.display = 'block';
                  tip.style.left = (evt.clientX + 12) + 'px';
                  tip.style.top = (evt.clientY - 28) + 'px';
                }}
                function hideStep4Tip() {{
                  var tip = document.getElementById('step4-plot-tooltip');
                  if (tip) tip.style.display = 'none';
                }}
              </script>
            </div>

            <!-- View 2: Unified 160-Molecule Latent Space Projections (PCA & t-SNE) -->
            <div style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:16px; margin-bottom:20px;">
              <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:8px;">
                <div style="font-weight:700; font-size:1.02rem; color:#0f172a;">
                  [Unified Plot 2] 통합 160종 잠재 공간 재정의 (24-Descriptor PCA &amp; t-SNE &bull; 120 Pos vs 40 Neg)
                </div>
                <span class="badge badge-4g" style="font-size:0.75rem;">24 Features (%V_bur + V_min + pi_hapticity + Chelate + Coord/Vol)</span>
              </div>
              <p style="font-size:0.82rem; color:#475569; margin:0 0 12px 0; line-height:1.55;">
                최신 <strong>160종 통합 평가 데이터셋</strong>({_u_ds['positive_total']} Positive vs {_u_ds['negative_controls_total']} Negative across 6 failure classes)에 24개 물리·위상·배위 디스크립터를 적용한 <code>PCA</code>(선형 2D 분산 설명률 {(_u_lat['pca_explained_variance_ratio'][0]+_u_lat['pca_explained_variance_ratio'][1])*100:.1f}%, Silhouette <strong>{_u_lat['pca_silhouette']:.3f}</strong>)와 <code>t-SNE</code>(비선형 국소 매니폴드, Silhouette <strong>{_u_lat['tsne_silhouette']:.3f}</strong>) 투영입니다. 2차원 선형 투영만으로는 다중 실패 메커니즘 분리(<code>PCA 2D Bal.Acc {_u_lat['pca_2d_balanced_acc']:.3f}</code>)가 불충분하지만, 24차원 전체 비선형 공간에서는 <code>SVC_RBF</code> Balanced Accuracy가 <strong>{_u_lat['full_space_svc_balanced_acc']:.3f}</strong>까지 상승합니다.
              </p>
              <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap:14px; margin-bottom:12px;">
                <div>{_pca_160_svg}</div>
                <div>{_tsne_160_svg}</div>
              </div>
              <div style="display:flex; gap:18px; flex-wrap:wrap; font-size:0.78rem; color:#334155;">
                <span><span style="display:inline-block; width:10px; height:10px; border-radius:50%; background:#059669; margin-right:5px;"></span>Audited Precursor (27)</span>
                <span><span style="display:inline-block; width:10px; height:10px; transform:rotate(45deg); background:#7c3aed; margin-right:6px;"></span>Expanded Ir/Rh/As/Ga/In/Ru (75)</span>
                <span><span style="display:inline-block; width:10px; height:10px; background:#0284c7; margin-right:5px;"></span>AS-ALD Inhibitor (18)</span>
                <span><span style="display:inline-block; width:0; height:0; border-left:6px solid transparent; border-right:6px solid transparent; border-bottom:10px solid #dc2626; margin-right:5px;"></span>Negative Control (40)</span>
              </div>
            </div>

            <!-- View 3: Unified 160-Molecule Model Comparison (CV-5 vs LOFMO Extrapolation) -->
            <div style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:16px; margin-bottom:20px;">
              <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:8px;">
                <div style="font-weight:700; font-size:1.02rem; color:#0f172a;">
                  [Unified Benchmark] 최신 160종 통합 평가 데이터셋 모델 성능 (CV-5 내삽 vs LOFMO 미지 실패 메커니즘 외삽)
                </div>
                <span class="badge badge-3g" style="font-size:0.75rem;">24 Features &bull; N=160 (120 Pos : 40 Neg)</span>
              </div>
              <div style="overflow-x:auto; margin-bottom:14px;">
                <table class="report-table" style="font-size:0.81rem;">
                  <thead>
                    <tr style="background:#f1f5f9;">
                      <th>Model / Rule</th>
                      <th>유형 및 학습 가정</th>
                      <th>Precursor TPR (102종) / CV-5 TPR</th>
                      <th>Neg Control TNR (40종) / CV-5 TNR</th>
                      <th>CV-5 Bal.Acc (AUC)</th>
                      <th>LOFMO 미지 메커니즘 평균 TNR</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td><strong>RDKit SAScore (&le; 3.5 / &gt; 4.0)</strong></td>
                      <td>1G Organic Fragment Statistics</td>
                      <td>TPR {_u_sum['baseline_sascore_precursors_tpr_pct']:.1f}% ({_u_sum['baseline_sascore_precursors_tp']})</td>
                      <td>TNR {_u_sum['baseline_sascore_neg_tnr_pct']:.1f}% ({_u_sum['baseline_sascore_neg_tn']})</td>
                      <td>N/A (Unsupervised)</td>
                      <td>{_u_sum['baseline_sascore_neg_tnr_pct']/100.0:.3f} (도메인 불일치)</td>
                    </tr>
                    <tr style="background:#f0fdf4;">
                      <td><strong>TMC-SA<sub>v2</sub> + Coord/Topology Gate</strong></td>
                      <td>Pro-ligand SA + Adaptive V<sub>min</sub> + 물리·위상 규칙</td>
                      <td><strong>TPR {_u_sum['precursors_all_tpr_pct']:.1f}% ({_u_sum['precursors_all_tp']})</strong></td>
                      <td><strong>TNR {_u_sum['negative_all_tnr_pct']:.1f}% ({_u_sum['negative_all_tn']})</strong></td>
                      <td>N/A (Rule + Physics)</td>
                      <td><strong>{_u_lofmo['Unified_TMC_SA_v2_Gate']['mean_unseen_tnr']:.3f}</strong> (규칙 기반 외삽 강건)</td>
                    </tr>
                    <tr>
                      <td><strong>OneClassSVM (RBF)</strong></td>
                      <td>Positive-Only Anomaly Detection</td>
                      <td>CV-5 TPR {_u_ml['OneClassSVM_RBF']['cv5_tpr_sensitivity']:.3f}</td>
                      <td>CV-5 TNR {_u_ml['OneClassSVM_RBF']['cv5_tnr_specificity']:.3f}</td>
                      <td>{_u_ml['OneClassSVM_RBF']['cv5_balanced_acc']:.3f} ({_u_ml['OneClassSVM_RBF']['cv5_roc_auc']:.3f})</td>
                      <td><strong>{_u_lofmo['OneClassSVM_RBF']['mean_unseen_tnr']:.3f}</strong> (ML 최고 외삽성)</td>
                    </tr>
                    <tr>
                      <td><strong>LogisticRegression (Balanced)</strong></td>
                      <td>Supervised Linear (24 Features)</td>
                      <td>CV-5 TPR {_u_ml['LogisticRegression']['cv5_tpr_sensitivity']:.3f}</td>
                      <td>CV-5 TNR {_u_ml['LogisticRegression']['cv5_tnr_specificity']:.3f}</td>
                      <td><strong>{_u_ml['LogisticRegression']['cv5_balanced_acc']:.3f}</strong> ({_u_ml['LogisticRegression']['cv5_roc_auc']:.3f})</td>
                      <td style="color:#b91c1c; font-weight:700;">{_u_lofmo['LogisticRegression']['mean_unseen_tnr']:.3f} (미지 메커니즘 붕괴)</td>
                    </tr>
                    <tr>
                      <td><strong>SVC (RBF, Balanced)</strong></td>
                      <td>Supervised Non-Linear Kernel</td>
                      <td>CV-5 TPR {_u_ml['SVC_RBF']['cv5_tpr_sensitivity']:.3f}</td>
                      <td>CV-5 TNR {_u_ml['SVC_RBF']['cv5_tnr_specificity']:.3f}</td>
                      <td>{_u_ml['SVC_RBF']['cv5_balanced_acc']:.3f} ({_u_ml['SVC_RBF']['cv5_roc_auc']:.3f})</td>
                      <td style="color:#b91c1c; font-weight:700;">{_u_lofmo['SVC_RBF']['mean_unseen_tnr']:.3f} (미지 메커니즘 붕괴)</td>
                    </tr>
                    <tr>
                      <td><strong>RandomForest / XGBoost (Weighted)</strong></td>
                      <td>Supervised Tree Ensemble</td>
                      <td>RF {_u_ml['RandomForest']['cv5_tpr_sensitivity']:.3f} / XGB {_u_ml['XGBoost_Weighted']['cv5_tpr_sensitivity']:.3f}</td>
                      <td>RF {_u_ml['RandomForest']['cv5_tnr_specificity']:.3f} / XGB {_u_ml['XGBoost_Weighted']['cv5_tnr_specificity']:.3f}</td>
                      <td>RF {_u_ml['RandomForest']['cv5_balanced_acc']:.3f} / XGB {_u_ml['XGBoost_Weighted']['cv5_balanced_acc']:.3f}</td>
                      <td style="color:#b91c1c; font-weight:700;">RF {_u_lofmo['RandomForest']['mean_unseen_tnr']:.3f} (학습 외 실패 미탐)</td>
                    </tr>
                  </tbody>
                </table>
              </div>

              <div style="font-weight:700; font-size:0.95rem; color:#0f172a; margin-bottom:8px;">
                6대 실패 메커니즘별 40종 Negative Control 정밀 차단률 (TNR Breakdown)
              </div>
              <div style="overflow-x:auto;">
                <table class="data-table">
                  <thead>
                    <tr>
                      <th>실패 메커니즘 (6개 층화 범주)</th>
                      <th>표본 수</th>
                      <th>SAScore 차단 (TNR)</th>
                      <th>SYBA 차단 (TNR)</th>
                      <th>TMC-SA<sub>v2</sub> + Coord Gate 차단 (TNR)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {_mech_rows_html}
                  </tbody>
                </table>
              </div>
            </div>

            <!-- View 4: Representative Molecules Diagnostic Summary & Item 2 SevenNet-0 Reaction Energy -->
            <div style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:16px; margin-bottom:20px;">
              <div style="font-weight:700; font-size:1.0rem; color:#0f172a; margin-bottom:8px;">
                [Diagnostic Samples &amp; Item 2 Reaction Energy] 주요 구제·차단 사례 및 SevenNet-0 &Delta;E<sub>syn</sub> 정합성
              </div>
              <div style="overflow-x:auto; margin-bottom:12px;">
                <table class="report-table" style="font-size:0.82rem;">
                  <thead>
                    <tr style="background:#f1f5f9;">
                      <th>Molecule (ID)</th>
                      <th>Classification</th>
                      <th style="text-align:center;">%V<sub>bur</sub> (Adaptive V<sub>min</sub>)</th>
                      <th style="text-align:center;">SevenNet-0 BDE<sub>min</sub></th>
                      <th style="text-align:center;">Topology (&pi; / Chelate)</th>
                      <th>Latest Verdict &amp; Diagnostic Note</th>
                    </tr>
                  </thead>
                  <tbody>
{_step4_rep_rows}
                  </tbody>
                </table>
              </div>
              <div style="background:#faf5ff; border:1px solid #e9d5ff; border-radius:6px; padding:11px 14px; font-size:0.81rem; color:#581c87; line-height:1.55;">
                <strong>Item 2 (SevenNet-0 화학양론적 합성 반응 에너지 &Delta;E<sub>syn</sub> 검증):</strong>
                알코올 분해(<code>SYN_TEOS</code>: <strong>{_rxn_teos:.2f} kcal/mol</strong>)는 기체상에서도 발열성 정합을 보이며, 염 복분해 경로는 암염 결정 <code>LiCl(s)</code> 격자 기준을 적용할 때 실제 합성 구동력이 정확히 재현됩니다(<code>SYN_TDMAT</code> {_rxn_tdmat_g:+.2f} &rarr; <strong>{_rxn_tdmat_s:.2f} kcal/mol</strong>, <code>SYN_TEMAH</code> {_rxn_temah_g:+.2f} &rarr; <strong>{_rxn_temah_s:.2f} kcal/mol</strong>; 캐시된 <code>MACE-MP0</code>은 로컬 <code>e3nn</code> 직렬화 비호환으로 미구동).
              </div>
            </div>

            <!-- View 5: Physical Rationale of Selected Descriptors & 4 Extrapolation Improvements -->
            <div class="card" style="border-left:4px solid #0284c7; margin-top:20px;">
              <h4 style="margin-top:0; color:#0f172a;">선택된 계산 물성(Descriptor)의 물리적 배경·공정 중요성 및 외삽(Extrapolation) 개선 방향</h4>
              <p style="font-size:0.84rem; color:#475569; margin-bottom:12px; line-height:1.6;">
                본 프레임워크의 물리·배위 Descriptor는 단순 통계 상관이 아니라 <strong>반도체 기상 공정(ALD/CVD) 프리커서의 단량체 휘발성·열적 배관 안정성·합성 분리 가능성</strong>을 결정짓는 미시 물리 법칙에 근거해 선정되었습니다.
              </p>
              <div style="overflow-x:auto; margin-bottom:14px;">
                <table class="report-table" style="font-size:0.81rem;">
                  <thead>
                    <tr style="background:#f1f5f9;">
                      <th>물성 그룹 및 핵심 Descriptor</th>
                      <th>물리·화학적 선정 배경 (기구적 의미)</th>
                      <th>반도체 ALD/CVD 공정 중요성 및 임계 경계</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td><strong>입체 배위 차폐</strong><br><code>%V<sub>bur</sub> (R=3.5 A)</code><br>적응형 하한 <code>V<sub>min</sub></code></td>
                      <td>중심 원자 반경 $R=3.5\text{{ Å}}$ 구 내 리간드 점유 부피 비율(Falivene et al., DOI <code>10.1038/s41557-019-0319-5</code>). 과소배위 시 이합체화·자연발화가 발생하며, 과밀 배위 시 결합 약화·승화 저해가 유발됨.</td>
                      <td><strong>Si·&pi;-햅틱 $V_{{\min}}=34\%$ / 일반 TM $44\% \sim 83.5\%$</strong>: <code>DIPAS</code>(41.2%)·<code>TSA</code>(35.8%)·<code>Ru(EtCp)2</code>(36.9%)의 정상 단량체성을 구제하고 <code>WPh5</code>류 과밀 아릴 비휘발군을 구별함.</td>
                    </tr>
                    <tr>
                      <td><strong>열분해 내성</strong><br><code>BDE<sub>min</sub></code>, <code>&Delta;E<sub>dimer</sub></code></td>
                      <td><code>SevenNet-0</code> 3D 이완 후 최약 리간드 해리 에너지($\text{{BDE}}_{{\min}}$) 및 비중첩 이합체(&ge;2.8 A) 응집 에너지($\Delta E_{{\text{{dimer}}}}$). 기체상 전달 중 최약 결합 절단이 열분해 온도 $T_{{\text{{decomp}}}}$를 지배함.</td>
                      <td><strong>$\text{{BDE}}_{{\min}} \ge 45\text{{ kcal/mol}}$</strong>: 버블러 가열 구간 내 조기 열분해(<code>CCTBA</code> CO 이탈, <code>TiEt4</code> &beta;-H 제거)와 수산화물 솔-겔 중합(<code>M(OH)n</code>)을 구분하는 열역학 축.</td>
                    </tr>
                    <tr>
                      <td><strong>배위·위상 안정성</strong><br><code>pi_hapticity_index</code>,<br><code>chelate_ring_count</code>,<br><code>is_homoleptic_tm_sp3</code></td>
                      <td>$\pi$-역공여(&eta;-Cp/diene) 및 킬레이트 고리 효과를 동종 $sp^3$ &sigma;-알킬(&beta;-수소 제거 취약)과 위상적으로 구분함.</td>
                      <td>SMILES에서 $sp^3$로 표기된 $\eta$-디엔 착물(<code>Ru(IMBCH)</code>, <code>Ru(EBCHD)</code>) 오탐 <code>0/75</code> 해소 및 불안정 호모렙틱 알킬(<code>TiEt4</code>, <code>ZrMe4</code>, <code>HfEt4</code>, <code>TaMe5</code>) <code>100%</code> 차단.</td>
                    </tr>
                    <tr>
                      <td><strong>휘발성·표면 반응성</strong><br><code>metal_oh_count</code>, <code>cooh_count</code>,<br><code>is_inert_alkane</code>, <code>strained_si_ring</code>,<br><code>TPSA / HBD / MW</code></td>
                      <td>분자 간 수소결합·옥솔레이션 중합(<code>M-OH</code>), 강프로톤 산 식각(<code>-COOH</code>), 비활성 알케인(흡착점 부재), 3/4각 실라사이클 고리 변형과 분자량·극성 표면적을 정량화함.</td>
                      <td>ALD 기화 불가 솔-겔(<code>Ti/Zr/Hf/Si(OH)4</code>)과 무반응 용매(<code>Decane</code>, <code>Cyclohexane</code>) 및 부식성 산(<code>TFA</code>, <code>TCA</code>)의 공정 부적합성을 차단함.</td>
                    </tr>
                  </tbody>
                </table>
              </div>

              <div style="background:#fff7ed; border:1px solid #fed7aa; border-radius:6px; padding:14px;">
                <div style="font-weight:700; font-size:0.9rem; color:#9a3412; margin-bottom:8px;">
                  통합 LOFMO 외삽 감사 결과 및 차기 4대 개선 과제
                </div>
                <p style="font-size:0.82rem; color:#7c2d12; margin:0 0 10px 0; line-height:1.55;">
                  통합 160종 데이터셋의 무작위 <code>CV-5</code>에서는 지도학습(<code>LogisticRegression</code> Bal.Acc <strong>{_u_ml['LogisticRegression']['cv5_balanced_acc']:.3f}</strong>, <code>SVC_RBF</code> <strong>{_u_ml['SVC_RBF']['cv5_balanced_acc']:.3f}</strong>)이 높게 보이나, 신규 실패 메커니즘 1종을 통째로 제외하는 <strong>LOFMO 외삽</strong>에서는 지도학습이 <code>0.000~0.024</code>로 붕괴합니다. 반면 양성(Positive) 분포만 학습하는 <strong><code>OneClassSVM_RBF</code>({_u_lofmo['OneClassSVM_RBF']['mean_unseen_tnr']:.3f})</strong>과 <strong>물리·위상 게이트({_u_lofmo['Unified_TMC_SA_v2_Gate']['mean_unseen_tnr']:.3f})</strong>가 미지 실패군을 방어합니다. Pro-ligand SA 단일 기준 적용 시 확인된 3건의 미탐(<code>NEG-019 CrPh3</code>, <code>NEG-039 Cl2SiH-SiHCl2</code>, <code>NEG-040 Cp2GaH</code>) 해소를 위해 다음 4가지 업그레이드가 요구됩니다:
                </p>
                <ul style="font-size:0.81rem; color:#7c2d12; margin:0; padding-left:18px; line-height:1.6;">
                  <li><strong>주기율표 블록 일반화 및 Si-Si/M-H 하이드라이드 규칙 추가</strong>: 고정 16원소 목록 외곽(<code>CrPh3</code>)과 불균등화 디실레인(<code>Cl2SiH-SiHCl2</code>)·불안정 갈륨 하이드라이드(<code>Cp2GaH</code>)를 포괄하도록 d/f/p-블록 산화수·<code>M-H</code>/<code>Si-Si</code> 규칙 및 <code>SevenNet-0</code> 물리 연산을 전수 적용해야 합니다.</li>
                  <li><strong>스핀·전하 인지 MLIP / DFT 결합</strong>: <code>SevenNet-0</code>(MPtrj 무기 결정 학습)은 중성 <code>CO</code>·할라이드·아렌 라디칼 해리에서 음의 BDE(<code>Ir(CO)3(tBu3Cyp) -8.8</code>, <code>GaCl3 -3.5 kcal/mol</code>)를 산출하므로, 분자 스핀·전하를 학습한 <code>MACE-OMol</code> 또는 <code>g-xTB</code>/$\omega\text{{B97M-V}}$로 BDE 산출층을 교체해야 합니다.</li>
                  <li><strong>응축상 &Delta;E<sub>syn</sub> 및 정량 증기압(P<sub>vap</sub>) 통합</strong>: <code>LiCl(s)</code> 암염 기준을 자동 반영한 합성 반응 에너지(&Delta;E<sub>syn</sub>)와 Clausius-Clapeyron 기화 창($\Delta T = T_{{\text{{decomp}}}} - T_{{\text{{vap}}}}$)을 점수식에 직접 결합해야 합니다.</li>
                  <li><strong>외삽 표준 프로토콜 상시화</strong>: 무작위 CV-5가 아닌 <strong>Leave-One-Element-Out (LOEO)</strong> 및 <strong>LOFMO</strong>를 정규 회귀 테스트 기준으로 채택해야 합니다.</li>
                </ul>
              </div>
            </div>"""

phase1_phase2_tuning_joined = f"""
        <!-- Latest Unified 160-Molecule Upgrade Callout -->
        <div style="background:#f0fdf4; border:1px solid #bbf7d0; border-left:4px solid #059669; border-radius:8px; padding:14px 18px; margin-bottom:20px; font-size:0.84rem; color:#14532d; line-height:1.6;">
          <strong>최신 버전 모델(<code>TMC-SA<sub>v2</sub></code> + Coord/Topology Gate) 통합 업그레이드 요약 (160종 전체 DB 기준):</strong>
          아래 표·데이터시트는 <code>SevenNet-0</code> 정밀 LBFGS 이완·이합체·분리 단편 BDE 기준선(<code>N=55</code> 감사 세트)을 기록하며, 최신 버전에서는 <strong>Si·&pi;-햅티시티 적응형 $V_{{\min}}=34\%$</strong> 및 <strong>6대 실패 메커니즘 배위 게이트</strong>를 결합하여 감사 세트 <strong>Precursor 25/27 (92.6%), Neg TNR 10/10 (100.0%)</strong> 및 전체 업데이트 <strong>160종 평가 데이터셋 Precursor {_u_sum['precursors_all_tp']} ({_u_sum['precursors_all_tpr_pct']:.1f}%), Neg TNR {_u_sum['negative_all_tn']} ({_u_sum['negative_all_tnr_pct']:.1f}%)</strong>을 달성하였습니다.
        </div>

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
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">1. [Phase 1 Benchmark] Data-Driven Base Model의 Dual Misclassification 분석</h3>
            <span class="badge" style="background:#fee2e2; color:#b91c1c; border:1px solid #fca5a5;">Baseline Error Analysis</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.92rem; line-height:1.65; margin-bottom:14px;">
            오픈소스 케모인포매틱스 Base Model(RDKit SAScore, Reaxys SCScore, SYBA, RAscore, BR-SAScore)을 물리 Descriptor 없이 55종 Dataset(프리커서 27종, 인히비터 18종, Negative Control 10종)에 대해 로컬 실행한 결과, 다음 두 가지 도메인 불일치 오류가 관측되었습니다.
          </p>

          <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-bottom:16px;">
            <div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:6px; padding:14px;">
              <div style="font-weight:700; color:#991b1b; font-size:0.95rem; margin-bottom:6px;">Error Type 1: Precursor False Negative ({base_cmp['precursors_sascore_fn']})</div>
              <p style="font-size:0.84rem; color:#7f1d1d; line-height:1.55; margin:0;">
                RDKit SAScore(&le; 3.5)는 27개 반도체 프리커서 중 {len(prec_sa_passers)}개({prec_sa_pass_names})만 통과시키고 <strong>{prec_sa_fn_cnt}개를 기각(FNR {base_cmp['precursors_sascore_fn']})</strong>합니다. 주요 양산 프리커서인 <code>TSA</code>(SA 6.91), <code>DMAH</code>(SA 6.31), <code>TMA</code>(SA 4.81), <code>TEMAH</code>(SA 4.10)가 모두 탈락하며, 이는 PubChem 의약 통계 내 Si-N, Al-C, Hf-N Fragment Frequency가 매우 낮기 때문입니다.
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
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
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
                  <td>정상 프리커서: TEMAH BDE <strong>{_temah.get('bde_min_kcal_mol', 0.0):.1f}</strong>, TDMAT <strong>{_tdmat.get('bde_min_kcal_mol', 0.0):.1f}</strong>, DIPAS <strong>{_dipas.get('bde_min_kcal_mol', 0.0):.1f} kcal/mol</strong><br>고리 변형/저안정: Silacyclobutadiene <strong>{_silcb.get('bde_min_kcal_mol', 0.0):.1f}</strong>, CoEt4 <strong>{_coet4.get('bde_min_kcal_mol', 0.0):.1f}</strong>, CCTBA <strong>{_cctba.get('bde_min_kcal_mol', 0.0):.1f} kcal/mol</strong></td>
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
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <h3 class="card-title" style="margin-bottom:0;">3. [Audit &amp; Re-evaluation] Uniform Analytical Screening 및 ML Cross-Validation (All-Sample vs Stratified 5-Fold CV vs LOOCV)</h3>
            <span class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc;">Audited Un-Hardcoded Evaluation</span>
          </div>
          <p style="color:var(--text-secondary); font-size:0.88rem; line-height:1.6; margin-bottom:14px;">
            SMILES 정정, LBFGS 이완·비중첩 이합체·분리 단편 BDE 재계산, Target Label Branching 제거, 개별 ID Hardcoding 제거, 및 <code>BR-SAScore</code> 연산 정상화를 모두 반영하여 재실행한 정량 결과입니다.
          </p>

          <div style="display:grid; grid-template-columns: 3fr 2fr; gap:18px; margin-bottom:20px;">
            <div>
              <div style="font-weight:700; color:#0f172a; margin-bottom:8px;">표 15-1. Label/ID-Free Uniform Analytical Screening ($\\text{{TMC-SA}}_{{\\text{{phys}}}} \\le 4.0$) 성능</div>
              <div style="overflow-x:auto;">
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
                양산 프리커서 27종, AS-ALD 인히비터 18종, Negative Control 10종에 대한 Phase 1(정정 SMILES &amp; BR-SAScore 포함) vs Phase 2(Relaxed LBFGS Physics + SMARTS Screening) 결과입니다.
              </p>
            </div>
            <div style="display:flex; gap:8px; align-items:center; width:100%; max-width:260px;">
              <input type="text" id="benchmarkSearch" onkeyup="searchBenchmark()" placeholder="분자명, 화학식, SMILES 검색..." style="padding:6px 12px; border:1px solid #cbd5e1; border-radius:6px; font-size:0.8rem; width:100%;">
            </div>
          </div>

          <div class="filter-bar" style="margin-bottom:14px; display:flex; flex-wrap:wrap; gap:8px;">
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

# -----------------------------------------------------------------------------
# DB Explorer Payload Builder (RDKit 2D SVGs + Unified Searchable DB Records)
# -----------------------------------------------------------------------------
def _render_mol2d_svg(smi: str) -> str:
    """Render compact standalone RDKit 2D SVG for a SMILES string."""
    if not smi:
        return ""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            return ""
        mol.UpdatePropertyCache(strict=False)
    try:
        rdDepictor.Compute2DCoords(mol)
    except Exception:
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(240, 165)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.08
    opts.bondLineWidth = 1.6
    opts.minFontSize = 11
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg.split("?>", 1)[-1]
    svg = re.sub(r"<!--.*?-->", "", svg, flags=re.S)
    svg = re.sub(r"\s+width=['\"][^'\"]*['\"]", "", svg, count=1)
    svg = re.sub(r"\s+height=['\"][^'\"]*['\"]", "", svg, count=1)
    svg = re.sub(r"\s+", " ", svg).strip()
    return svg


expanded_by_id = {x.get("id"): x for x in expanded_prec_db}
bench_by_id = {b["id"]: b for b in benchmarks}
svg_dict = {}


def _norm_elem(val: str) -> str:
    v = (val or "").strip()
    if not v or v in {"C", "N", "O", "S", "Cl", "F"}:
        return "Organic"
    return v


def _ensure_svg(smi: str) -> bool:
    if not smi:
        return False
    if smi not in svg_dict:
        rendered = _render_mol2d_svg(smi)
        if rendered:
            svg_dict[smi] = rendered
    return smi in svg_dict


db_records = []

# 1. Baseline 27 Precursors (enriched with AtomicLimits, SA, and SevenNet-0 physics)
for p in prec_db:
    pid = p["id"]
    smi = p.get("smiles") or ""
    _ensure_svg(smi)
    exp = expanded_by_id.get(pid, {})
    ev = mol_by_id.get(pid, {})
    bm = bench_by_id.get(pid, {})
    props = p.get("physical_properties", {})
    synth = p.get("synthesis_info", {})
    db_records.append({
        "id": pid,
        "db_type": "precursor",
        "db_label": "Baseline Precursor",
        "is_expanded_master": True,
        "name": p["name"],
        "chemical_name": p.get("chemical_name", "-"),
        "formula": p.get("formula", "-"),
        "mw": p.get("molecular_weight"),
        "cas": p.get("cas_number", "-"),
        "element": p.get("metal_center", "-"),
        "category": p.get("category", "-"),
        "smiles": smi,
        "has_2d": bool(smi and smi in svg_dict),
        "target_films": p.get("target_films", "-"),
        "process_type": p.get("process_type", "-"),
        "vapor_pressure": props.get("vapor_pressure_torr", "-"),
        "t_decomp_c": props.get("thermal_decomposition_temp_c"),
        "state_at_rt": props.get("state_at_rt", "-"),
        "bp_c": props.get("boiling_point_c"),
        "mp_c": props.get("melting_point_c"),
        "starting_materials": synth.get("starting_materials", []),
        "reaction_type": synth.get("reaction_type", "-"),
        "yield_pct": synth.get("typical_yield_pct"),
        "byproduct": synth.get("byproduct", "-"),
        "vendors": p.get("commercial_vendors", []),
        "doi": p.get("reference_doi", ""),
        "atomiclimits_name": exp.get("atomiclimits_name", ""),
        "atomiclimits_process_count": exp.get("atomiclimits_process_count", 0),
        "atomiclimits_materials": exp.get("atomiclimits_materials", []),
        "atomiclimits_co_reactants": exp.get("atomiclimits_co_reactants", []),
        "atomiclimits_doi_count": exp.get("atomiclimits_unique_doi_count", 0),
        "sascore": ev.get("sascore"),
        "scscore": ev.get("scscore"),
        "syba": ev.get("syba"),
        "rascore": ev.get("rascore"),
        "br_sascore": ev.get("br_sascore"),
        "ligand_sa": ev.get("ligand_sa"),
        "rxnflow_match": bm.get("rxnflow_match"),
        "rxnflow_rules": bm.get("rxnflow_templates", []),
        "reasyn_in_vocab": bm.get("reasyn_in_vocab"),
        "bde_min": ev.get("bde_min_kcal_mol"),
        "vbur": ev.get("vbur_percent"),
        "de_dimer": ev.get("delta_e_dimer_kcal_mol"),
        "e_coh": ev.get("e_coh_ev_per_atom"),
        "tmc_sa_phys": ev.get("tmc_sa_phys"),
        "phase2_pass": ev.get("tuned_with_chem_pass"),
        "domain_flag": ev.get("chemical_domain_flag", "-"),
    })

# 2. Expanded AtomicLimits Precursors (Top-120 + Ir/Rh/As/Ga/In/Ru expansion)
for idx, exp in enumerate(expanded_prec_db, 1):
    pid = exp.get("id", f"EXP-{idx:03d}")
    if pid in prec_map:
        continue
    smi = exp.get("smiles") or ""
    _ensure_svg(smi)
    ev7 = exp_sevenn_by_id.get(pid, {})
    db_records.append({
        "id": pid,
        "db_type": "expanded",
        "db_label": "AtomicLimits ALD",
        "is_expanded_master": True,
        "name": exp.get("name", "-"),
        "chemical_name": exp.get("chemical_name") or exp.get("atomiclimits_name") or "-",
        "formula": exp.get("formula") or "-",
        "mw": exp.get("molecular_weight"),
        "cas": exp.get("cas_number") or "-",
        "element": exp.get("center_element") or "-",
        "category": exp.get("category") or exp.get("origin", "AtomicLimits ALD Database"),
        "smiles": smi,
        "has_2d": bool(smi and smi in svg_dict),
        "target_films": ", ".join(exp.get("atomiclimits_materials", [])) or "-",
        "process_type": "ALD",
        "atomiclimits_name": exp.get("atomiclimits_name", ""),
        "atomiclimits_process_count": exp.get("atomiclimits_process_count", 0),
        "atomiclimits_materials": exp.get("atomiclimits_materials", []),
        "atomiclimits_co_reactants": exp.get("atomiclimits_co_reactants", []),
        "atomiclimits_doi_count": exp.get("atomiclimits_unique_doi_count", 0),
        "bde_min": ev7.get("bde_min_kcal_mol"),
        "vbur": ev7.get("vbur_percent"),
        "de_dimer": ev7.get("delta_e_dimer_kcal_mol"),
        "e_coh": ev7.get("e_coh_ev_per_atom"),
    })

# 3. AS-ALD Inhibitors (18 entries)
for inh in inhib_list:
    iid = inh.get("id", "-")
    smi = inh.get("smiles") or ""
    _ensure_svg(smi)
    ev = mol_by_id.get(iid, {})
    lit = inh.get("literature_reference", {})
    m_obj = Chem.MolFromSmiles(smi) if smi else None
    calc_elem = "Si" if (m_obj and any(a.GetSymbol() == "Si" for a in m_obj.GetAtoms())) else "Organic"
    db_records.append({
        "id": iid,
        "db_type": "inhibitor",
        "db_label": "AS-ALD Inhibitor",
        "is_expanded_master": False,
        "name": inh.get("abbreviation", "-"),
        "chemical_name": inh.get("name", "-"),
        "formula": inh.get("formula", "-"),
        "mw": inh.get("molecular_weight"),
        "cas": inh.get("cas_number", "-"),
        "element": _norm_elem(ev.get("center_atom") or calc_elem),
        "category": inh.get("chemical_family", "-"),
        "smiles": smi,
        "has_2d": bool(smi and smi in svg_dict),
        "selectivity_type": inh.get("selectivity_type", "-"),
        "target_surface": inh.get("target_passivation_surface", "-"),
        "growth_surface": inh.get("growth_surface", "-"),
        "passivation_mechanism": inh.get("passivation_mechanism", "-"),
        "surface_bond": inh.get("surface_bond_type", "-"),
        "delivery_phase": inh.get("delivery_phase", "-"),
        "removal_method": inh.get("removal_method", "-"),
        "logp": inh.get("logp"),
        "doi": lit.get("doi", ""),
        "authors": f"{lit.get('authors', '-')} ({lit.get('year', '-')})",
        "sascore": ev.get("sascore"),
        "scscore": ev.get("scscore"),
        "syba": ev.get("syba"),
        "rascore": ev.get("rascore"),
        "br_sascore": ev.get("br_sascore"),
        "ligand_sa": ev.get("ligand_sa"),
        "bde_min": ev.get("bde_min_kcal_mol"),
        "vbur": ev.get("vbur_percent"),
        "de_dimer": ev.get("delta_e_dimer_kcal_mol"),
        "e_coh": ev.get("e_coh_ev_per_atom"),
        "tmc_sa_phys": ev.get("tmc_sa_phys"),
        "phase2_pass": ev.get("tuned_with_chem_pass"),
        "domain_flag": ev.get("chemical_domain_flag", "-"),
    })

# 4. Building Blocks / Commercial Reagents (29 entries)
for bb in bb_db:
    smi = bb.get("smiles") or ""
    _ensure_svg(smi)
    m_obj = Chem.MolFromSmiles(smi) if smi else None
    mw_val = round(Descriptors.MolWt(m_obj), 2) if m_obj else None
    elem = "-"
    if m_obj:
        metals = [a.GetSymbol() for a in m_obj.GetAtoms() if a.GetSymbol() not in {"C", "H", "N", "O", "F", "Cl", "Br", "I", "S", "P", "B"}]
        elem = metals[0] if metals else "Organic"
    db_records.append({
        "id": bb["id"],
        "db_type": "reagent",
        "db_label": "Reagent / BB",
        "is_expanded_master": False,
        "name": bb["name"],
        "chemical_name": bb["name"],
        "formula": bb.get("formula", "-"),
        "mw": mw_val,
        "cas": bb.get("cas", "-"),
        "element": elem,
        "category": bb.get("role", "-"),
        "smiles": smi,
        "has_2d": bool(smi and smi in svg_dict),
        "purity": bb.get("purity", "-"),
        "handling": bb.get("handling", "-"),
    })

# 5. Accredited Synthesis Reactions (12 entries)
for r in rxn_db:
    reactants = r.get("reactants", [])
    products = r.get("products", [])
    for s in reactants + products:
        _ensure_svg(s)
    elem = "-"
    for s in products:
        m_obj = Chem.MolFromSmiles(s)
        if m_obj:
            metals = [a.GetSymbol() for a in m_obj.GetAtoms() if a.GetSymbol() not in {"C", "H", "N", "O", "F", "Cl", "Br", "I", "S", "P", "B"}]
            if metals:
                elem = metals[0]
                break
    db_records.append({
        "id": r["id"],
        "db_type": "reaction",
        "db_label": "Reaction Template",
        "is_expanded_master": False,
        "name": r["name"],
        "chemical_name": r.get("reaction_type", "-"),
        "formula": " + ".join(reactants) + " -> " + " + ".join(products),
        "mw": None,
        "cas": "-",
        "element": elem,
        "category": r.get("reaction_type", "-"),
        "smiles": r.get("reaction_smarts", ""),
        "reactants": reactants,
        "products": products,
        "byproduct": r.get("byproduct", "-"),
        "solvent": r.get("solvent", "-"),
        "conditions": r.get("conditions", "-"),
        "dg_kcal_mol": r.get("thermodynamics_dG_kcal_mol"),
        "has_2d": any(s in svg_dict for s in reactants + products),
    })

# 6. Negative Controls (10 entries)
for ev in molecules_eval:
    if not ev["id"].startswith("NEG-"):
        continue
    smi = ev.get("smiles") or ""
    _ensure_svg(smi)
    db_records.append({
        "id": ev["id"],
        "db_type": "negative",
        "db_label": "Negative Control",
        "is_expanded_master": False,
        "name": ev["name"],
        "chemical_name": ev.get("domain", "Negative Control"),
        "formula": ev.get("formula", "-"),
        "mw": ev.get("mw"),
        "cas": "-",
        "element": _norm_elem(ev.get("center_atom", "-")),
        "category": ev.get("domain", "-"),
        "smiles": smi,
        "has_2d": bool(smi and smi in svg_dict),
        "sascore": ev.get("sascore"),
        "scscore": ev.get("scscore"),
        "syba": ev.get("syba"),
        "rascore": ev.get("rascore"),
        "br_sascore": ev.get("br_sascore"),
        "ligand_sa": ev.get("ligand_sa"),
        "bde_min": ev.get("bde_min_kcal_mol"),
        "vbur": ev.get("vbur_percent"),
        "de_dimer": ev.get("delta_e_dimer_kcal_mol"),
        "e_coh": ev.get("e_coh_ev_per_atom"),
        "tmc_sa_phys": ev.get("tmc_sa_phys"),
        "phase2_pass": ev.get("tuned_with_chem_pass"),
        "domain_flag": ev.get("chemical_domain_flag", "-"),
        "phys_penalty_detail": ev.get("phys_penalty_detail", ""),
    })

db_explorer_payload_json = json.dumps(
    {"svgs": svg_dict, "records": db_records},
    ensure_ascii=False,
).replace("</", "<\\/")


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
final_html = final_html.replace("<!-- DB_EXPLORER_DATA_PLACEHOLDER -->", db_explorer_payload_json)

with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)

with open("synthetic_accessibility_report.html", "w", encoding="utf-8") as f:
    f.write(final_html)

print("Generated clean light-theme index.html & synthetic_accessibility_report.html successfully!")
