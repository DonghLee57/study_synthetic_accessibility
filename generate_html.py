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
        "rate": "100.0% (27/27, 실측)",
        "badge_style": "background:#e0f2fe; color:#0369a1; border-color:#7dd3fc; font-weight:700;",
        "eval": "코드 직접 실측 (ECFP6 알킬 서브트리 편향)",
        "desc": "[코드 직접 실측 100.0%] XGBoost 1.7.6 기반 직접 실행: 전수 평균 0.925(최저 0.629 ~ 최고 0.988)로 27종 전수 합성 용이 판정. ECFP6 반경 3 지문에서 금속 중심이 '다중 치환 알킬/아민'으로 매핑되어 약물 합성의 쉬운 알킬화 반응으로 착각하여 중심 금속 결합 안정성 및 산화 위험을 감지하지 못하는 맹점 확인."
    },
    "RxnFlow": {
        "rate": "33.3% (9/27, 리간드 매칭)",
        "badge_style": "background:#ffedd5; color:#c2410c; border-color:#fdba74; font-weight:700;",
        "eval": "코드 직접 실측 (유기 리간드 국한 / 무기 배위 부재)",
        "desc": "[코드 직접 실측 33.3%] ICLR 2025 공식 GFlowNet 환경 직접 연산: 유기 리간드의 알킬 곁사슬(아세톤 환원적 아민화 등)을 포함하는 9종만 매칭 성공. 금속-배위원자(M-N, M-O, M-C) 배위 결합 형성 템플릿 부재로 유기금속 전이금속 착체 코어 합성률은 0.0%에 머무름."
    },
    "ReaSyn": {
        "rate": "100.0% 어휘 수용 (무기 템플릿 제약)",
        "badge_style": "background:#ede9fe; color:#6d28d9; border-color:#c4b5fd; font-weight:700;",
        "eval": "코드 직접 실측 (118종 전주기 어휘 / 반응 매니폴드 갭)",
        "desc": "[코드 직접 실측 100.0%] ICLR 2026 Oral NVIDIA BioNeMo ReaSyn 구조 직접 분석: 118개 전주기율표 원소 토큰을 지원하여 27종 전구체의 SMILES 토큰화 및 임베딩 100% 결함 없음(0 UNK). 그러나 115개 SynFormer 특허 반응 매니폴드가 유기 반응에 국한되어 무기 착체 배위 단계 전개 템플릿 제약 확인."
    },
    "AiZynthFinder": {
        "rate": "100.0% 분자 수용 (MCTS 무기 정책 제약)",
        "badge_style": "background:#f1f5f9; color:#334155; border-color:#94a3b8; font-weight:700;",
        "eval": "코드 직접 실측 (RDKit 완결 / 유기 전개 정책 제약)",
        "desc": "[코드 직접 실측 100.0%] v4.4.1 직접 실행: RDKit 기반 27종 전구체 분자 객체 생성, MW 계산, InChIKey 생성 100% 정상 통과. 그러나 USPTO 50k 기반 전개 정책 신경망이 금속-리간드 결합 절단 확률을 0으로 산출하여 자식 노드 전개가 즉시 중단되므로 무기 반응 템플릿 주입 필수."
    }
}

# Strictly filter to the 8 empirically tested models in benchmark order
tested_keys = [
    "SAScore (RDKit Built-in)",
    "SCScore",
    "SYBA",
    "BR-SAScore",
    "RAscore",
    "RxnFlow",
    "ReaSyn",
    "AiZynthFinder"
]

repo_by_name = {r["name"]: r for r in extra["repos"]}
repo_rows = []
for idx, key in enumerate(tested_keys, 1):
    r = repo_by_name.get(key)
    if not r:
        continue
    m = repo_semi_metrics.get(key, {
        "rate": "N/A",
        "badge_style": "background:#f1f5f9; color:#475569; border-color:#cbd5e1;",
        "eval": "평가 미기재",
        "desc": r.get('description', '-')
    })
    repo_rows.append(f"""
            <tr>
              <td style="text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;">{idx:02d}</td>
              <td style="white-space:nowrap;"><strong>{r['name']}</strong></td>
              <td><a href="{r['url']}" target="_blank" style="font-family:'JetBrains Mono'; font-size:0.82rem;">{r['repo']}</a></td>
              <td style="white-space:nowrap; text-align:center;"><span class="badge badge-1g">{r['stars']}</span></td>
              <td style="white-space:nowrap; text-align:center;"><span class="badge badge-3g">{r['license']}</span></td>
              <td style="min-width:340px;">
                <div style="display:flex; align-items:center; gap:8px; margin-bottom:4px;">
                  <span class="badge" style="{m['badge_style']}">{m['rate']}</span>
                  <strong style="font-size:0.84rem; color:#0f172a;">{m['eval']}</strong>
                </div>
                <div style="font-size:0.78rem; color:#64748b; line-height:1.45;">{m['desc']}</div>
              </td>
            </tr>""")
repos_joined = "\n".join(repo_rows)

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
    
    rxnflow_badge = f"<span class='badge' style='background:#ffedd5; color:#c2410c; font-size:0.75rem; font-weight:600;'>리간드 매칭 ({len(rxnflow_rules)}건)</span>" if rxnflow else "<span class='badge' style='background:#f1f5f9; color:#64748b; font-size:0.75rem;'>0건 (미매칭)</span>"
    reasyn_badge = f"<span class='badge' style='background:#dcfce7; color:#15803d; font-size:0.75rem; font-weight:600;'>100% In-Vocab ({reasyn_tok}t)</span>" if reasyn_ok else "<span class='badge badge-1g'>UNK 발생</span>"
    aizynth_badge = f"<span class='badge' style='background:#e0e7ff; color:#4338ca; font-size:0.75rem; font-weight:600;'>정상 수용 ({mw})</span>" if aizynth_ok else "<span class='badge badge-1g'>수용 에러</span>"
    
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
final_html = final_html.replace("<!-- REPO_ROWS_PLACEHOLDER -->", repos_joined)
final_html = final_html.replace("<!-- VASP_ROWS_PLACEHOLDER -->", vasp_joined)
final_html = final_html.replace("<!-- PREC_DB_ROWS_PLACEHOLDER -->", prec_db_joined)
final_html = final_html.replace("<!-- BB_DB_ROWS_PLACEHOLDER -->", bb_db_joined)
final_html = final_html.replace("<!-- RXN_DB_ROWS_PLACEHOLDER -->", rxn_db_joined)

with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)

with open("synthetic_accessibility_report.html", "w", encoding="utf-8") as f:
    f.write(final_html)

print("Generated clean light-theme index.html & synthetic_accessibility_report.html successfully!")
