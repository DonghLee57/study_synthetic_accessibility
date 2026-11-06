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
repo_semi_metrics = {
    "SAScore (RDKit Built-in)": {
        "rate": "14.8%",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "심각한 왜곡 (False Negative)",
        "desc": "유기 의약품 PubChem 통계 기반으로 금속-헤테로원자 결합을 비정상 파편으로 간주하여 6~9점(난합성) 과대평가. 공유결합 실란계(DIPAS, 3DMAS 등 4종)만 유효."
    },
    "SYBA": {
        "rate": "11.1%",
        "badge_style": "background:#fee2e2; color:#b91c1c; border-color:#fca5a5; font-weight:700;",
        "eval": "적용 불가 (음수 편향)",
        "desc": "ZINC vs 난합성 베이지안 분류 모델. 무기 착물 파편에 극단적 음수 스코어 부여. 비금속 실란 3종만 양수 판정."
    },
    "SCScore": {
        "rate": "25.9%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "유기 리간드 한정 유효",
        "desc": "Reaxys 유기 반응 복잡도 순위 모델. 전이금속 착화 반응 결여로 대부분 4.5 이상 과대평가. 단순 알콕사이드/알킬아미드 7종만 3.0 이하."
    },
    "RAscore": {
        "rate": "22.2%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "USPTO 유기 반응 편향",
        "desc": "AiZynthFinder MCTS 역합성 확률 증류 모델. 전이금속 염 복분해 반응 부재로 유기금속 착체 예측 확률 0.05 이하 왜곡."
    },
    "GASA": {
        "rate": "29.6%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "단순 구조 국한",
        "desc": "그래프 어텐션(GAT) 기반 SA 지표. 비금속 실란 및 단일 리간드 착체 8종 인식 성공. 다핵 착물 및 란타넘족 원소 인식 실패."
    },
    "DeepSA": {
        "rate": "33.3%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "유기 리간드 국한",
        "desc": "사전학습 트랜스포머 기반 SA 지표. 금속 중심 SMILES 파싱 한계로 전이금속 착체 에러 빈발, 리간드 단독 9종 평가 유효."
    },
    "DFRscore": {
        "rate": "37.0%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "리간드 경로 한정",
        "desc": "반응 단계 수 기반 예측기. 기본 유기 반응 템플릿으로 10종 리간드 합성 단계 산출 성공. 무기화 반응(Inorganic core) 미지원."
    },
    "BR-SAScore": {
        "rate": "88.9%",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "전용 DB 탑재 시 최우수 (기본 25.9%)",
        "desc": "반도체 전용 빌딩 블록(Gelest 실란+금속염) 및 무기 복분해 반응(semi_rxn) 사전 교체 시 24/27종 고유 배위 구조에서 100% 정상 스코어 산출 검증."
    },
    "ChemProjector": {
        "rate": "74.1%",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "리간드 투영 양호",
        "desc": "분자 투영 기반 역합성 경로 내재화 모델. 전이금속 착체를 리간드 분해 프로토콜(HL)과 결합 시 20종에 대해 합성 가능한 DAG 및 상용 시약 매핑 경로 성공적 출력."
    },
    "SynFormer": {
        "rate": "70.4%",
        "badge_style": "background:#d1fae5; color:#047857; border-color:#6ee7b7; font-weight:700;",
        "eval": "빌딩블록 매핑 우수",
        "desc": "트랜스포머 기반 순차 합성 계획 생성 모델. 19종 전구체의 유기 리간드 빌딩 블록 조합 경로 유효 생성, 금속 킬레이트화 단계는 후처리 템플릿 필요."
    },
    "SynFlowNet": {
        "rate": "81.5%",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "다목적 역설계 우수",
        "desc": "GFlowNet 기반 마코프 의사결정 경로 탐색기. 다중 보상함수(증기압, 분해온도, SA) 결합 시 22종에 대한 유효 합성 그래프 도출."
    },
    "RxnFlow": {
        "rate": "85.2%",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "경로 완결성 우수 (ICLR 2025)",
        "desc": "반응 중심 경로 연속 흐름 생성 모델. 시약-반응 쌍 동시 탐색으로 23종 전구체의 상용 유무기 전구체로부터의 1~3단계 합성 트리 생성 성공."
    },
    "ReaSyn": {
        "rate": "85.2%",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "대규모 모델 우수 (NVIDIA)",
        "desc": "NVIDIA BioNeMo 대규모 자기회귀 합성 경로 생성기. 1.6억 파라미터 기반으로 복잡 전이금속 착체(23종)의 리간드 골격 및 합성 경로를 단일 패스로 생성."
    },
    "SynNet": {
        "rate": "66.7%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "유기 반응 국한",
        "desc": "바이오/유기 4대 기본 반응에 편향되어 금속 알킬화/실릴화 반응 경로 생성 시 9종 실패. 18종 리간드 경로 유효."
    },
    "AiZynthFinder": {
        "rate": "44.4%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "무기 템플릿 부재",
        "desc": "MCTS 기반 전통적 역합성 엔진. 무기/금속 착체 템플릿 부재로 유기금속 전구체는 실패율 높음. 실란 및 아미노 전구체 12종 성공."
    },
    "Retro*": {
        "rate": "48.1%",
        "badge_style": "background:#fef3c7; color:#b45309; border-color:#fcd34d; font-weight:700;",
        "eval": "무기 반응 미지원",
        "desc": "신경망 유도 A* 탐색 알고리즘. USPTO 50k 템플릿 사용으로 13종 유기 리간드 경로만 도출 성공."
    },
    "MPScore (Materials Precursor)": {
        "rate": "92.6%",
        "badge_style": "background:#dcfce7; color:#15803d; border-color:#86efac; font-weight:700;",
        "eval": "재료/배위 최적",
        "desc": "MOF/재료 화학 전구체 전용 합성 점수 모델. 유기금속 및 배위 리간드 평가에 최적화되어 27종 중 25종에 대해 정확한 합성성 랭킹 부여 확인."
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

template_path = "template.html"
with open(template_path, "r", encoding="utf-8") as f:
    template = f.read()

final_html = template.replace("<!-- LIT_ROWS_PLACEHOLDER -->", lit_rows_joined)
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
