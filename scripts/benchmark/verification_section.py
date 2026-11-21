"""Dataset Provenance Verification Tab Builder (``출처 검증``) for the 304-Entity Dataset.

Generates the standalone HTML content for the provenance verification tab
(``#tab-content-verification``, button label ``출처 검증``), enumerating all
304 baseline dataset entries (205 ALD/CVD precursors, 18 AS-ALD inhibitors,
29 reagents, 12 synthesis reactions, and 40 negative controls) with source
provenance, clickable DOI / AtomicLimits / PubChem verification links, and
persistent checkbox + note tracking via ``localStorage`` and JSON export/import.
"""

from __future__ import annotations

import html
import json
import os
import re
import urllib.parse
from typing import Any, Dict, List


def _load_json(rel_path: str) -> Any:
    with open(rel_path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _extract_list(obj: Any, *keys: str) -> List[Dict[str, Any]]:
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        for k in keys:
            if k in obj and isinstance(obj[k], list):
                return obj[k]
    return []


def _pubchem_url(query: str) -> str:
    q = (query or "").strip()
    if not q or q == "-":
        return "https://pubchem.ncbi.nlm.nih.gov/"
    return f"https://pubchem.ncbi.nlm.nih.gov/#query={urllib.parse.quote(q)}"


def _doi_link(doi: str, label: str = "") -> str:
    d = re.sub(r"^(https?://)?(dx\.)?doi\.org/", "", (doi or "").strip(), flags=re.IGNORECASE)
    if not d:
        return ""
    disp = html.escape(label or f"DOI: {d}")
    href = f"https://doi.org/{urllib.parse.quote(d, safe='/')}"
    return (
        f'<a href="{href}" target="_blank" rel="noopener noreferrer" '
        f'style="color:#0369a1; font-weight:600; text-decoration:none; border-bottom:1px dotted #7dd3fc;">{disp}</a>'
    )


def _ext_link(url: str, label: str, color: str = "#0f766e") -> str:
    return (
        f'<a href="{html.escape(url)}" target="_blank" rel="noopener noreferrer" '
        f'style="color:{color}; font-weight:600; text-decoration:none; border-bottom:1px dotted {color};">{html.escape(label)}</a>'
    )


def build_verification_html() -> str:
    """Build the interactive 304-entity ``출처 검증`` tab HTML fragment."""
    prec27_list = _extract_list(_load_json(os.path.join("database", "precursors.json")), "precursors")
    prec27_by_id = {p.get("id"): p for p in prec27_list}
    prec27_by_name = {(p.get("name") or "").lower(): p for p in prec27_list}

    exp_master = _extract_list(
        _load_json(os.path.join("database", "external", "expanded_semiconductor_precursors_master.json")),
        "precursors",
        "records",
    )
    al_list = _extract_list(
        _load_json(os.path.join("database", "external", "atomiclimits_precursors_extracted.json")),
        "precursors",
    )
    al_by_name: Dict[str, Dict[str, Any]] = {}
    for rec in al_list:
        pname = (rec.get("precursor") or "").strip()
        if pname:
            al_by_name[pname.lower()] = rec
        for alias in rec.get("aliases_in_source") or []:
            al_by_name[str(alias).strip().lower()] = rec

    pat_list = _extract_list(
        _load_json(os.path.join("database", "external", "patent_lit_expansion_candidates.json")),
        "records",
        "candidates",
    )
    pat_doi_by_name: Dict[str, str] = {}
    for pr in pat_list:
        refs = str(pr.get("references") or "")
        m_doi = re.search(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", refs)
        if m_doi:
            pat_doi_by_name[(pr.get("name") or "").strip().lower()] = m_doi.group(0).rstrip(".,;)")

    inh_list = _extract_list(
        _load_json(os.path.join("database", "external", "as_ald_inhibitors.json")),
        "inhibitors",
    )
    bb_list = _extract_list(_load_json(os.path.join("database", "building_blocks.json")), "building_blocks")
    rxn_list = _extract_list(_load_json(os.path.join("database", "reactions.json")), "reactions")
    neg_list = _extract_list(
        _load_json(os.path.join("database", "expanded_negative_controls.json")),
        "negative_controls",
        "records",
    )

    rows_data: List[Dict[str, str]] = []

    # 1. 205 ALD/CVD Precursors (27 Baseline + 178 Expanded)
    for r in exp_master:
        rid = str(r.get("id") or "")
        name = str(r.get("name") or "")
        elem = str(r.get("center_element") or "-")
        cas = str(r.get("cas_number") or "").strip()
        origin = str(r.get("origin") or r.get("property_source") or "AtomicLimits ALD DB")
        al_name = str(r.get("atomiclimits_name") or name).strip()
        al_rec = al_by_name.get(al_name.lower()) or al_by_name.get(name.lower())
        p27 = prec27_by_id.get(rid) or prec27_by_name.get(name.lower())

        links: List[str] = []
        source_desc = origin
        if p27 and p27.get("reference_doi"):
            doi_val = str(p27["reference_doi"]).strip()
            links.append(_doi_link(doi_val))
            source_desc = f"Baseline 원저 문헌 ({origin})"
        elif name.lower() in pat_doi_by_name:
            links.append(_doi_link(pat_doi_by_name[name.lower()]))
            source_desc = "AtomicLimits ALD DB + 특허·문헌 검증"
        elif al_rec and al_rec.get("dois"):
            first_doi = str(al_rec["dois"][0]).strip()
            if first_doi.startswith("10."):
                links.append(_doi_link(first_doi, f"DOI: {first_doi}"))
            proc_n = al_rec.get("process_count", 0)
            doi_n = al_rec.get("unique_doi_count", len(al_rec.get("dois") or []))
            source_desc = f"AtomicLimits ALD DB (공정 {proc_n}건 · 문헌 {doi_n}편)"

        if al_rec is not None:
            links.append(_ext_link("https://www.atomiclimits.com/alddatabase/", "AtomicLimits DB", "#0f766e"))
        pc_cid = r.get("pubchem_cid") or (p27.get("pubchem_cid") if p27 else None)
        if pc_cid:
            cas_sfx = f" · CAS {cas}" if (cas and cas != "-") else ""
            links.append(_ext_link(f"https://pubchem.ncbi.nlm.nih.gov/compound/{pc_cid}", f"PubChem CID {pc_cid}{cas_sfx}", "#475569"))
        elif cas and cas != "-":
            links.append(_ext_link(_pubchem_url(cas), f"PubChem ({cas})", "#475569"))
        elif not links:
            links.append(_ext_link(_pubchem_url(name), "PubChem 검색", "#475569"))

        rows_data.append({
            "id": rid,
            "cat_key": "precursor",
            "cat_label": "프리커서",
            "name": name,
            "elem": elem,
            "source": source_desc,
            "links_html": " &middot; ".join(links),
        })

    # 2. 18 AS-ALD Inhibitors
    for r in inh_list:
        rid = str(r.get("id") or "")
        abbr = str(r.get("abbreviation") or "")
        name = str(r.get("name") or "")
        disp_name = f"{abbr} ({name})" if abbr and abbr != name else name
        cas = str(r.get("cas_number") or r.get("cas") or "").strip()
        lit = r.get("literature_reference") or {}
        doi_val = str(lit.get("doi") or "").strip()
        links = []
        if doi_val:
            links.append(_doi_link(doi_val))
            src_desc = f"AS-ALD 문헌 ({lit.get('authors', '')} {lit.get('year', '')})".strip()
        else:
            src_desc = f"AS-ALD 문헌 · {r.get('category', 'Inhibitor')}"
        links.append(_ext_link(_pubchem_url(cas if (cas and cas != "-") else name), f"PubChem ({cas})" if (cas and cas != "-") else "PubChem 검색", "#475569"))
        rows_data.append({
            "id": rid,
            "cat_key": "inhibitor",
            "cat_label": "인히비터",
            "name": disp_name,
            "elem": "Si" if "silane" in name.lower() or "silazane" in name.lower() else "Organic",
            "source": src_desc,
            "links_html": " &middot; ".join(links),
        })

    # 3. 29 Reagents / Building Blocks
    for r in bb_list:
        rid = str(r.get("id") or "")
        name = str(r.get("name") or "")
        cas = str(r.get("cas_number") or r.get("cas") or "").strip()
        cat = str(r.get("category") or "상용 시약")
        purity = str(r.get("purity") or "").strip()
        src_desc = f"상용 시약 카탈로그 ({cat}{', ' + purity if purity else ''})"
        links = [_ext_link(_pubchem_url(cas if (cas and cas != "-") else name), f"PubChem ({cas})" if (cas and cas != "-") else "PubChem 검색", "#475569")]
        rows_data.append({
            "id": rid,
            "cat_key": "reagent",
            "cat_label": "시약",
            "name": name,
            "elem": str(r.get("formula") or "-"),
            "source": src_desc,
            "links_html": " &middot; ".join(links),
        })

    # 4. 12 Synthesis Reactions
    for r in rxn_list:
        rid = str(r.get("id") or "")
        name = str(r.get("name") or "")
        target_p = str(r.get("target_precursor") or r.get("product_name") or "")
        p27 = prec27_by_name.get(target_p.lower())
        if not p27:
            for p_item in prec27_list:
                if (p_item.get("name") or "").lower() in name.lower():
                    p27 = p_item
                    break
        links = []
        if p27 and p27.get("reference_doi"):
            links.append(_doi_link(str(p27["reference_doi"])))
        links.append(_ext_link(_pubchem_url(target_p or name), "PubChem 검색", "#475569"))
        src_desc = f"공인 합성 반응 ({r.get('reaction_type', 'Synthesis')})"
        rows_data.append({
            "id": rid,
            "cat_key": "reaction",
            "cat_label": "합성 반응",
            "name": name,
            "elem": str(r.get("solvent") or "-"),
            "source": src_desc,
            "links_html": " &middot; ".join(links),
        })

    # 5. 40 Negative Controls
    for r in neg_list:
        rid = str(r.get("id") or "")
        name = str(r.get("name") or "")
        mech = str(r.get("failure_mechanism") or "Negative Control")
        orig = str(r.get("origin") or "Curated Negative Control")
        src_desc = f"{orig} ({mech})"
        links = [_ext_link(_pubchem_url(name.split("(")[0].strip()), "PubChem 검색", "#475569")]
        rows_data.append({
            "id": rid,
            "cat_key": "negative",
            "cat_label": "Negative",
            "name": name,
            "elem": "-",
            "source": src_desc,
            "links_html": " &middot; ".join(links),
        })

    total_n = len(rows_data)
    cat_counts: Dict[str, int] = {}
    for row in rows_data:
        cat_counts[row["cat_key"]] = cat_counts.get(row["cat_key"], 0) + 1

    cat_badge_styles = {
        "precursor": "background:#eff6ff; color:#1d4ed8; border:1px solid #bfdbfe;",
        "inhibitor": "background:#f0fdfa; color:#0f766e; border:1px solid #99f6e4;",
        "reagent": "background:#fefce8; color:#a16207; border:1px solid #fde68a;",
        "reaction": "background:#f5f3ff; color:#6d28d9; border:1px solid #ddd6fe;",
        "negative": "background:#fef2f2; color:#b91c1c; border:1px solid #fecaca;",
    }

    trs_html: List[str] = []
    for row in rows_data:
        rid_esc = html.escape(row["id"])
        search_blob = html.escape(
            f"{row['id']} {row['cat_label']} {row['name']} {row['elem']} {row['source']}".lower()
        )
        b_style = cat_badge_styles.get(row["cat_key"], "background:#f1f5f9; color:#334155;")
        trs_html.append(
            f'<tr class="ver-row" data-vid="{rid_esc}" data-vcat="{row["cat_key"]}" data-vtext="{search_blob}">'
            f'<td style="text-align:center;"><input type="checkbox" class="ver-chk" data-vid="{rid_esc}" onchange="toggleVerifyCheck(\'{rid_esc}\', this.checked)" aria-label="{rid_esc} 확인"></td>'
            f'<td style="font-family:\'JetBrains Mono\',monospace; font-size:0.76rem; color:#475569; white-space:nowrap;">{rid_esc}</td>'
            f'<td style="white-space:nowrap;"><span style="font-size:0.7rem; font-weight:700; padding:2px 7px; border-radius:999px; {b_style}">{html.escape(row["cat_label"])}</span></td>'
            f'<td style="font-weight:700; color:#0f172a;">{html.escape(row["name"])} <span style="font-weight:500; color:#64748b; font-size:0.76rem;">({html.escape(row["elem"])})</span></td>'
            f'<td style="font-size:0.8rem; color:#334155;">{html.escape(row["source"])}</td>'
            f'<td style="font-size:0.78rem; line-height:1.55;">{row["links_html"]}</td>'
            f'<td><input type="text" class="ver-note" data-vid="{rid_esc}" placeholder="검증 메모..." oninput="saveVerifyNote(\'{rid_esc}\', this.value)" style="width:100%; min-width:120px; padding:4px 7px; border:1px solid #cbd5e1; border-radius:5px; font-size:0.76rem; background:#ffffff;"></td>'
            f'</tr>'
        )

    return f"""
      <div class="db-controls-card" style="margin-bottom:14px;">
        <div style="display:flex; flex-wrap:wrap; justify-content:space-between; align-items:center; gap:10px; margin-bottom:10px;">
          <div>
            <h2 style="margin:0; font-size:1.25rem; font-weight:800; color:#0f172a;">데이터셋 출처 전수 검증 ({total_n}종)</h2>
            <p style="margin:3px 0 0 0; font-size:0.82rem; color:#475569;">
              각 분자의 출처 확보 경로와 원문(DOI · AtomicLimits · PubChem) 링크를 직접 확인하고 체크박스와 메모로 기록할 수 있습니다. (브라우저 자동 저장 및 JSON 백업 지원)
            </p>
          </div>
          <div style="display:flex; align-items:center; gap:8px; flex-wrap:wrap;">
            <span id="verifyProgressBadge" style="font-family:'JetBrains Mono',monospace; font-weight:800; font-size:0.84rem; color:#0f766e; background:#f0fdfa; border:1px solid #99f6e4; padding:5px 11px; border-radius:999px;">
              확인 완료: 0 / {total_n} (0.0%)
            </span>
            <button type="button" class="db-pill" onclick="exportVerifyState()" style="font-size:0.76rem; padding:5px 10px;">JSON 내보내기</button>
            <label class="db-pill" style="font-size:0.76rem; padding:5px 10px; cursor:pointer; margin:0;">
              JSON 불러오기
              <input type="file" accept=".json" style="display:none;" onchange="importVerifyState(event)">
            </label>
            <button type="button" class="btn-toc-utility" style="width:auto; padding:5px 10px; font-size:0.76rem;" onclick="resetVerifyChecks()">체크 초기화</button>
          </div>
        </div>

        <div style="height:7px; background:#e2e8f0; border-radius:999px; overflow:hidden; margin-bottom:12px;">
          <div id="verifyProgressBar" style="height:100%; width:0%; background:linear-gradient(90deg, #0d9488, #0284c7); transition:width 0.2s ease;"></div>
        </div>

        <div class="db-search-row" style="margin-bottom:10px;">
          <input type="text" id="verifySearchInput" class="db-search-input" placeholder="ID, 분자명, 원소, 출처, DOI 검색..." oninput="filterVerifyTable()">
          <div class="db-pill-bar" id="verifyStatusPills" style="margin:0;">
            <button type="button" class="db-pill active" data-vstat="all" onclick="setVerifyStatus('all', this)">전체 상태</button>
            <button type="button" class="db-pill" data-vstat="unchecked" onclick="setVerifyStatus('unchecked', this)">미확인만</button>
            <button type="button" class="db-pill" data-vstat="checked" onclick="setVerifyStatus('checked', this)">확인 완료만</button>
          </div>
          <span id="verifyShownBadge" style="font-family:'JetBrains Mono',monospace; font-size:0.8rem; font-weight:700; color:#0369a1;">표시: {total_n}건</span>
        </div>

        <div class="db-pill-bar" id="verifyCatPills" style="margin-bottom:0;">
          <button type="button" class="db-pill active" data-vcat="all" onclick="setVerifyCategory('all', this)">전체 ({total_n})</button>
          <button type="button" class="db-pill" data-vcat="precursor" onclick="setVerifyCategory('precursor', this)">프리커서 ({cat_counts.get('precursor', 0)})</button>
          <button type="button" class="db-pill" data-vcat="inhibitor" onclick="setVerifyCategory('inhibitor', this)">인히비터 ({cat_counts.get('inhibitor', 0)})</button>
          <button type="button" class="db-pill" data-vcat="reagent" onclick="setVerifyCategory('reagent', this)">시약 ({cat_counts.get('reagent', 0)})</button>
          <button type="button" class="db-pill" data-vcat="reaction" onclick="setVerifyCategory('reaction', this)">합성 반응 ({cat_counts.get('reaction', 0)})</button>
          <button type="button" class="db-pill" data-vcat="negative" onclick="setVerifyCategory('negative', this)">Negative Control ({cat_counts.get('negative', 0)})</button>
        </div>
      </div>

      <div class="card" style="padding:0; overflow:hidden;">
        <div style="overflow-x:auto;">
          <table class="data-table" id="verifyTable" style="width:100%; border-collapse:collapse; font-size:0.82rem; margin:0;">
            <thead>
              <tr>
                <th style="width:52px; text-align:center;">확인</th>
                <th style="width:108px;">ID</th>
                <th style="width:84px;">분류</th>
                <th style="width:220px;">분자명 (원소)</th>
                <th style="width:240px;">출처 확보 경로</th>
                <th style="min-width:220px;">검증 링크</th>
                <th style="width:170px;">확인 메모</th>
              </tr>
            </thead>
            <tbody id="verifyTableBody">
              {''.join(trs_html)}
            </tbody>
          </table>
        </div>
      </div>
    """
