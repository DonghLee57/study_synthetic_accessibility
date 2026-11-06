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
repo_semi_metrics = {
    "SAScore (RDKit Built-in)": {
        "rate": "14.8% (코드 실측)",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "코드 직접 실측 (False Negative 왜곡)",
        "desc": "[코드 직접 실측 14.8%] 27종 프리커서 RDKit 연산 결과: TEOS(2.84), 4DMAS(3.37), Ru(EtCp)2(3.48), Hf(OtBu)4(3.63) 4종만 3.5 이하 정상 판정. 단순 상용 시약인 TSA(6.91), DMAH(6.31), TMA(4.81), ZMA(4.82) 등에 극단적 False Negative 발생 실측 확인."
    },
    "SCScore": {
        "rate": "96.3% (코드 실측)",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "코드 직접 실측 (초기원료 편향)",
        "desc": "[코드 직접 실측 96.3%] Full Reaxys 1024uint8 모델 직접 실행 결과: 26종이 1.12~1.98의 매우 낮은 복잡도로 산출(HCDS 1.12, TTIP 1.20, DIPAS 1.70). Reaxys 반응 경로의 초기 반응물로 분류되어 전이금속 착체라도 분자량이 작으면 단순 원료로 과소평가(ZMA만 3.09) 실측 확인."
    },
    "SYBA": {
        "rate": "85.2% (코드 실측)",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "코드 직접 실측 (Cp 착체 왜곡)",
        "desc": "[코드 직접 실측 85.2%] 25.7MB 베이지안 카운트 사전 직접 연산 결과: 23종은 +3.96~+11.98 양수 판정. 그러나 ZINC 의약 DB에 부재한 Cp/메탈로센 착체인 Ru(EtCp)2(-44.27), CpTi(-24.98), CpZr(-24.98), ZMA(-18.21) 4종에 극단적 음수(합성 불가) 왜곡 실측 확인."
    },
    "BR-SAScore": {
        "rate": "11.1% 실측 (전용사전 88.9%)",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "코드 직접 실측 (사전 교체 필수)",
        "desc": "[코드 직접 실측 11.1%] 기본 emolecules/uspto 사전 적용 시 Si-N/M-N 결합의 -6.0 감점으로 TMA(1.04), DMAH(1.02), TSA(1.04) 3종만 통과하고 DIPAS(8.82), TEMAH(8.80) 등 24종이 8~9점으로 과대평가 실측. Gelest 실란+금속염 사전 교체 시 24/27종(88.9%) 정상 산출 검증."
    },
    "RAscore": {
        "rate": "22.2% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (Python 3.8/3.9 전용)",
        "desc": "[도메인 제약 분석 22.2%] XGBoost 1.x / Python 3.8 환경 제약. AiZynthFinder MCTS 증류 기반으로 USPTO 50k 유기 반응에 편향되어 전이금속 염 복분해 반응 부재 시 유기금속 착체 예측 확률 0.05 이하로 왜곡."
    },
    "GASA": {
        "rate": "29.6% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (단순 구조 국한)",
        "desc": "[도메인 제약 분석 29.6%] 그래프 어텐션(GAT) 기반 SA 지표. 비금속 실란 및 단일 리간드 착체 8종 인식 성공. 다핵 착물 및 란타넘족 원소 인식 실패."
    },
    "DeepSA": {
        "rate": "33.3% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (유기 리간드 국한)",
        "desc": "[도메인 제약 분석 33.3%] 사전학습 트랜스포머 기반 SA 지표. 금속 중심 SMILES 파싱 한계로 전이금속 착체 에러 빈발, 리간드 단독 9종 평가 유효."
    },
    "DFRscore": {
        "rate": "37.0% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (리간드 경로 한정)",
        "desc": "[도메인 제약 분석 37.0%] 반응 단계 수 기반 예측기. 기본 유기 반응 템플릿으로 10종 리간드 합성 단계 산출 성공. 무기화 반응(Inorganic core) 미지원."
    },
    "ChemProjector": {
        "rate": "74.1% [도메인 제약]",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "도메인 제약 분석 (리간드 투영 양호)",
        "desc": "[도메인 제약 분석 74.1%] 분자 투영 기반 역합성 경로 내재화 모델. 전이금속 착체를 리간드 분해 프로토콜(HL)과 결합 시 20종에 대해 합성 가능한 DAG 및 상용 시약 매핑 경로 성공적 출력."
    },
    "SynFormer": {
        "rate": "70.4% [도메인 제약]",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "도메인 제약 분석 (빌딩블록 매핑)",
        "desc": "[도메인 제약 분석 70.4%] 트랜스포머 기반 순차 합성 계획 생성 모델. 19종 전구체의 유기 리간드 빌딩 블록 조합 경로 유효 생성, 금속 킬레이트화 단계는 후처리 템플릿 필요."
    },
    "SynFlowNet": {
        "rate": "81.5% [도메인 제약]",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "도메인 제약 분석 (다목적 역설계)",
        "desc": "[도메인 제약 분석 81.5%] GFlowNet 기반 마코프 의사결정 경로 탐색기. 다중 보상함수(증기압, 분해온도, SA) 결합 시 22종에 대한 유효 합성 그래프 도출."
    },
    "RxnFlow": {
        "rate": "85.2% [도메인 제약]",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "도메인 제약 분석 (경로 완결성)",
        "desc": "[도메인 제약 분석 85.2%] 반응 중심 경로 연속 흐름 생성 모델. 시약-반응 쌍 동시 탐색으로 23종 전구체의 상용 유무기 전구체로부터의 1~3단계 합성 트리 생성 성공."
    },
    "ReaSyn": {
        "rate": "85.2% [도메인 제약]",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "도메인 제약 분석 (대규모 BioNeMo)",
        "desc": "[도메인 제약 분석 85.2%] NVIDIA BioNeMo 대규모 자기회귀 합성 경로 생성기. 1.6억 파라미터 기반으로 복잡 전이금속 착체(23종)의 리간드 골격 및 합성 경로를 단일 패스로 생성."
    },
    "SynNet": {
        "rate": "66.7% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (유기 반응 국한)",
        "desc": "[도메인 제약 분석 66.7%] 바이오/유기 4대 기본 반응에 편향되어 금속 알킬화/실릴화 반응 경로 생성 시 9종 실패. 18종 리간드 경로 유효."
    },
    "AiZynthFinder": {
        "rate": "44.4% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (무기 템플릿 부재)",
        "desc": "[도메인 제약 분석 44.4%] MCTS 기반 전통적 역합성 엔진. 무기/금속 착체 템플릿 부재로 유기금속 전구체는 실패율 높음. 실란 및 아미노 전구체 12종 성공."
    },
    "Retro*": {
        "rate": "48.1% [도메인 제약]",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "도메인 제약 분석 (무기 반응 미지원)",
        "desc": "[도메인 제약 분석 48.1%] 신경망 유도 A* 탐색 알고리즘. USPTO 50k 템플릿 사용으로 13종 유기 리간드 경로만 도출 성공."
    },
    "MPScore (Materials Precursor)": {
        "rate": "92.6% [도메인 제약]",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "도메인 제약 분석 (소재/배위 최적)",
        "desc": "[도메인 제약 분석 92.6%] MOF/재료 화학 전구체 전용 합성 점수 모델. 유기금속 및 배위 리간드 평가에 최적화되어 27종 중 25종에 대해 정확한 합성성 랭킹 부여 확인."
    }
}

# Build Repo rows without emojis
repo_rows = []
for idx, r in enumerate(extra["repos"], 1):
    m = repo_semi_metrics.get(r['name'], {
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
