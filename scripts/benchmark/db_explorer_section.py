"""Domain DB Explorer (Tab 3) & Dataset Architecture Sankey SVG builder.

Extracted from ``generate_html.py`` to modularize the 286-entity canonical
deduplicated database payload, RDKit 2D structure dictionary, and proportional
Sankey diagram generation.
"""

from __future__ import annotations

import json
import os
import re
from typing import Any, Dict

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from html_common import (
    canon_smi as _canon_smi,
    norm_chem_key as _norm_chem_key,
    norm_elem as _norm_elem,
    render_mol2d_svg as _render_mol2d_svg,
)


def build_db_explorer_and_sankey(shared_ctx: Dict[str, Any]) -> Dict[str, str]:
    """Build ``DB_EXPLORER_DATA_PLACEHOLDER`` and ``DATASET_SANKEY_PLACEHOLDER``.

    Physics/Algorithm:
        Normalizes and deduplicates 342 curated precursor, inhibitor, reagent,
        reaction, negative-control, and patent-expansion records into 286 unique
        chemical entities with minified RDKit 2D SVGs.

    Args:
        shared_ctx: Context dictionary returned by ``build_main_report_sections()``.

    Returns:
        Dictionary mapping placeholder tokens to rendered strings.
    """
    prec_db = shared_ctx["prec_db"]
    bb_db = shared_ctx["bb_db"]
    rxn_db = shared_ctx["rxn_db"]
    benchmarks = shared_ctx["benchmarks"]
    expanded_prec_db = shared_ctx["expanded_prec_db"]
    exp_sevenn_by_id = shared_ctx["exp_sevenn_by_id"]
    inhib_list = shared_ctx["inhib_list"]
    mol_by_id = shared_ctx["mol_by_id"]
    prec_map = shared_ctx["prec_map"]
    molecules_eval = shared_ctx["molecules_eval"]
    _u_ds = shared_ctx["_u_ds"]
    _u_by_id = shared_ctx["_u_by_id"]

    expanded_by_id = {x.get("id"): x for x in expanded_prec_db}
    bench_by_id = {b["id"]: b for b in benchmarks}
    svg_dict = {}



    def _ensure_svg(smi: str) -> bool:
        if not smi:
            return False
        if smi not in svg_dict:
            rendered = _render_mol2d_svg(smi)
            if rendered:
                svg_dict[smi] = rendered
        return smi in svg_dict


    expansion_cand_path = os.path.join("database", "external", "patent_lit_expansion_candidates.json")
    with open(expansion_cand_path, "r", encoding="utf-8") as f:
        patent_exp_payload = json.load(f)
    patent_exp_records = patent_exp_payload.get("records", [])

    db_records = []
    rec_by_id = {}
    rec_by_csmi = {}
    rec_by_norm_key = {}




    def _register_db_rec(rec: dict) -> None:
        rec.setdefault("alt_db_types", [])
        rec.setdefault("aliases", [])
        rec.setdefault("secondary_ids", [])
        rec.setdefault("role_badges", [])
        db_records.append(rec)
        rec_by_id[rec["id"]] = rec
        cs = _canon_smi(rec.get("smiles", ""))
        if cs and rec.get("db_type") != "reaction" and cs not in rec_by_csmi:
            rec_by_csmi[cs] = rec
        if rec.get("db_type") in {"precursor", "expanded"}:
            for k_raw in (rec.get("name"), rec.get("formula"), rec.get("atomiclimits_name")):
                nk = _norm_chem_key(k_raw)
                if nk and nk not in rec_by_norm_key:
                    rec_by_norm_key[nk] = rec


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
        _register_db_rec({
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

    # 2. Expanded AtomicLimits Precursors (deduplicating pure notation-alias pairs inside master)
    NOTATION_ALIAS_TARGET = {
        "EXT-ALD-133": "EXT-ALD-121",  # [Ga(NMe2)3]2 == Ga2(NMe2)6
        "EXT-ALD-179": "EXT-ALD-178",  # GaMe3[MeO(CH2)2NHtBu] == GaMe3(CH3OCH2CH2NHtBu)
    }
    pending_notation_aliases = []

    for idx, exp in enumerate(expanded_prec_db, 1):
        pid = exp.get("id", f"EXP-{idx:03d}")
        if pid in prec_map:
            continue
        if pid in NOTATION_ALIAS_TARGET:
            pending_notation_aliases.append((pid, exp))
            continue
        smi = exp.get("smiles") or ""
        _ensure_svg(smi)
        ev7 = exp_sevenn_by_id.get(pid, {})
        _register_db_rec({
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

    for alias_id, exp in pending_notation_aliases:
        target = rec_by_id.get(NOTATION_ALIAS_TARGET[alias_id])
        if target:
            aname = exp.get("name") or exp.get("atomiclimits_name") or ""
            if aname and aname not in target["aliases"] and aname != target["name"]:
                target["aliases"].append(aname)
            if alias_id not in target["secondary_ids"]:
                target["secondary_ids"].append(alias_id)
            for mat in exp.get("atomiclimits_materials", []):
                if mat not in target["atomiclimits_materials"]:
                    target["atomiclimits_materials"].append(mat)
            target["atomiclimits_process_count"] = (target.get("atomiclimits_process_count") or 0) + (exp.get("atomiclimits_process_count") or 0)

    # 3. AS-ALD Inhibitors (18 entries)
    for inh in inhib_list:
        iid = inh.get("id", "-")
        smi = inh.get("smiles") or ""
        _ensure_svg(smi)
        ev = mol_by_id.get(iid, {})
        lit = inh.get("literature_reference", {})
        m_obj = Chem.MolFromSmiles(smi) if smi else None
        calc_elem = "Si" if (m_obj and any(a.GetSymbol() == "Si" for a in m_obj.GetAtoms())) else "Organic"
        _register_db_rec({
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

    # 4. Building Blocks / Commercial Reagents (29 entries; merge 10 dual-role reagents into canonical precursor records)
    DUAL_ROLE_BB_EXPLICIT = {
        "BB_02": "EXT-ALD-117",  # Dichlorosilane (SiH2Cl2 == SiCl2H2)
    }
    n_merged_dual_role = 0
    for bb in bb_db:
        smi = bb.get("smiles") or ""
        _ensure_svg(smi)
        m_obj = Chem.MolFromSmiles(smi) if smi else None
        mw_val = round(Descriptors.MolWt(m_obj), 2) if m_obj else None
        cs = _canon_smi(smi)
        target = None
        if bb["id"] in DUAL_ROLE_BB_EXPLICIT:
            target = rec_by_id.get(DUAL_ROLE_BB_EXPLICIT[bb["id"]])
        if not target and cs:
            target = rec_by_csmi.get(cs)
        if not target:
            target = rec_by_norm_key.get(_norm_chem_key(bb.get("formula"))) or rec_by_norm_key.get(_norm_chem_key(bb.get("name")))

        if target:
            n_merged_dual_role += 1
            if "reagent" not in target["alt_db_types"]:
                target["alt_db_types"].append("reagent")
            if "Dual-Role: Reagent" not in target["role_badges"]:
                target["role_badges"].append("Dual-Role: Reagent")
            if bb["id"] not in target["secondary_ids"]:
                target["secondary_ids"].append(bb["id"])
            if bb["name"] not in target["aliases"] and bb["name"] != target["name"]:
                target["aliases"].append(bb["name"])
            if (not target.get("smiles")) and smi:
                target["smiles"] = smi
                target["has_2d"] = bool(smi in svg_dict)
                if cs and cs not in rec_by_csmi:
                    rec_by_csmi[cs] = target
            if (not target.get("formula") or target.get("formula") == "-") and bb.get("formula"):
                target["formula"] = bb["formula"]
            if target.get("mw") is None and mw_val is not None:
                target["mw"] = mw_val
            if (not target.get("cas") or target.get("cas") == "-") and bb.get("cas"):
                target["cas"] = bb["cas"]
            if not target.get("purity"):
                target["purity"] = bb.get("purity", "-")
            if not target.get("handling"):
                target["handling"] = bb.get("handling", "-")
            continue

        elem = "-"
        if m_obj:
            metals = [a.GetSymbol() for a in m_obj.GetAtoms() if a.GetSymbol() not in {"C", "H", "N", "O", "F", "Cl", "Br", "I", "S", "P", "B"}]
            elem = metals[0] if metals else "Organic"
        _register_db_rec({
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
        _register_db_rec({
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
        _register_db_rec({
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

    # 7. KR/JP/CN/US/TW Patent & Metal-Centered Literature Expansion (38 entries: backfill 14 AtomicLimits stubs + add 24 net-new)
    PATENT_BACKFILL_TARGET = {
        "EXP_02": "EXT-ALD-058",  # MoCl5
        "EXP_03": "EXT-ALD-066",  # Mo(CO)6
        "EXP_07": "EXT-ALD-045",  # Nb(OEt)5
        "EXP_08": "EXT-ALD-113",  # TBTDEN == Nb(NtBu)(NEt2)3
        "EXP_10": "EXT-ALD-120",  # PDMAT == Ta(NMe2)5
        "EXP_11": "EXT-ALD-055",  # TaCl5
        "EXP_14": "EXT-ALD-041",  # HfCl4
        "EXP_15": "EXT-ALD-029",  # TiCl4
        "EXP_20": "EXT-ALD-072",  # Cu(dmap)2
        "EXP_22": "EXT-ALD-067",  # MeCpPtMe3 == Pt(CpMe)Me3
        "EXP_24": "EXT-ALD-081",  # Y(MeCp)3 == Y(CpMe)3
        "EXP_25": "EXT-ALD-028",  # ZnEt2 (DEZ)
        "EXP_26": "EXT-ALD-085",  # Ge(NMe2)4
        "EXP_27": "EXT-ALD-031",  # Sn(NMe2)4 (TDMASn)
    }
    n_backfilled_al = 0
    n_net_new_patent = 0

    for exp_rec in patent_exp_records:
        eid = exp_rec["id"]
        smi = exp_rec.get("smiles") or ""
        _ensure_svg(smi)
        cs = _canon_smi(smi)
        m_obj = Chem.MolFromSmiles(smi) if smi else None
        formula_str = rdMolDescriptors.CalcMolFormula(m_obj) if m_obj else "-"
        mw_val = round(Descriptors.MolWt(m_obj), 2) if m_obj else None
        elem = _norm_elem(exp_rec.get("center_atom", "-"))
        diag_note = exp_rec.get("boundary_flag") or (
            f"All SevenNet-0 & SAPIS physical/chemical gates passed (BDE={exp_rec.get('bde_min_kcal_mol')} kcal/mol, TMC={exp_rec.get('tmc_sa_7net')})"
            if exp_rec.get("pass_7net")
            else ("All ETKDG SAPIS physical & chemical gates passed" if exp_rec.get("unified_pass") else "")
        )
        has_7n = exp_rec.get("tmc_sa_7net") is not None
        vbur_best = exp_rec.get("vbur_rprobe_7net") if has_7n else exp_rec.get("vbur_d2_percent")
        tmc_best = exp_rec.get("tmc_sa_7net") if has_7n else exp_rec.get("tmc_sa_no3d")
        pass_best = exp_rec.get("pass_7net") if has_7n else exp_rec.get("unified_pass")

        target = rec_by_id.get(PATENT_BACKFILL_TARGET.get(eid, ""))
        if not target and cs:
            target = rec_by_csmi.get(cs)

        if target:
            n_backfilled_al += 1
            if "patent_expansion" not in target["alt_db_types"]:
                target["alt_db_types"].append("patent_expansion")
            if "Patent/Lit Backfilled" not in target["role_badges"]:
                target["role_badges"].append("Patent/Lit Backfilled")
            if has_7n and "SevenNet-0 3D" not in target["role_badges"]:
                target["role_badges"].append("SevenNet-0 3D")
            if eid not in target["secondary_ids"]:
                target["secondary_ids"].append(eid)
            if exp_rec["name"] not in target["aliases"] and exp_rec["name"] != target["name"]:
                target["aliases"].append(exp_rec["name"])
            if (not target.get("smiles")) and smi:
                target["smiles"] = smi
                target["has_2d"] = bool(smi in svg_dict)
                if cs and cs not in rec_by_csmi:
                    rec_by_csmi[cs] = target
            if (not target.get("formula") or target.get("formula") == "-") and formula_str != "-":
                target["formula"] = formula_str
            if target.get("mw") is None and mw_val is not None:
                target["mw"] = mw_val
            if (not target.get("element") or target.get("element") == "-") and elem:
                target["element"] = elem
            if not target.get("target_films") or target.get("target_films") == "-":
                target["target_films"] = exp_rec.get("application", "-")
            else:
                target["target_films"] = f"{target['target_films']} | {exp_rec.get('application', '')}".strip(" |")
            target["sascore"] = target.get("sascore") if target.get("sascore") is not None else exp_rec.get("sascore")
            target["br_sascore"] = target.get("br_sascore") if target.get("br_sascore") is not None else exp_rec.get("br_sascore")
            target["ligand_sa"] = target.get("ligand_sa") if target.get("ligand_sa") is not None else exp_rec.get("ligand_sa")
            target["vbur"] = vbur_best if vbur_best is not None else target.get("vbur")
            target["vbur_d2_7net"] = exp_rec.get("vbur_d2_7net")
            target["vbur_etkdg"] = exp_rec.get("vbur_d2_percent")
            target["bde_min"] = exp_rec.get("bde_min_kcal_mol") if target.get("bde_min") is None else target.get("bde_min")
            target["e_coh"] = exp_rec.get("e_coh_ev_per_atom") if target.get("e_coh") is None else target.get("e_coh")
            target["de_dimer"] = exp_rec.get("delta_e_dimer_kcal_mol") if target.get("de_dimer") is None else target.get("de_dimer")
            target["dg_dimer_373"] = exp_rec.get("delta_g_dimer_373_kcal_mol")
            target["de_form_bond"] = exp_rec.get("delta_e_form_per_bond_kcal_mol")
            target["tmc_sa_no3d"] = exp_rec.get("tmc_sa_no3d")
            target["tmc_sa_phys"] = tmc_best if tmc_best is not None else target.get("tmc_sa_phys")
            target["phase2_pass"] = pass_best if pass_best is not None else target.get("phase2_pass")
            target["coord_number"] = exp_rec.get("cn_eff")
            target["valence_electrons"] = exp_rec.get("n_ve")
            target["jurisdictions"] = exp_rec.get("jurisdictions", "")
            target["patent_references"] = exp_rec.get("references", "")
            target["mechanism_note"] = exp_rec.get("mechanism_note", "")
            target["diag_note"] = diag_note
            continue

        n_net_new_patent += 1
        is_inh = exp_rec.get("target_type") == "inhibitor"
        badges = ["Net-New Expansion"]
        if has_7n:
            badges.append("SevenNet-0 3D")
        _register_db_rec({
            "id": eid,
            "db_type": "patent_expansion",
            "db_label": "Patent/Lit Inhibitor" if is_inh else "Patent/Lit Precursor",
            "is_expanded_master": False,
            "alt_db_types": ["inhibitor"] if is_inh else [],
            "role_badges": badges,
            "name": exp_rec["name"],
            "chemical_name": f"{exp_rec.get('family_class', '-')} ({exp_rec.get('element_group', '-')})",
            "formula": formula_str,
            "mw": mw_val,
            "cas": "-",
            "element": elem,
            "category": f"{exp_rec.get('element_group', '-')} | {exp_rec.get('family_class', '-')}",
            "smiles": smi,
            "has_2d": bool(smi and smi in svg_dict),
            "target_films": "-" if is_inh else exp_rec.get("application", "-"),
            "target_surface": exp_rec.get("application", "-") if is_inh else "",
            "process_type": "AS-ALD / SAM Passivation" if is_inh else "ALD / CVD (Patent & Lit)",
            "jurisdictions": exp_rec.get("jurisdictions", ""),
            "patent_references": exp_rec.get("references", ""),
            "mechanism_note": exp_rec.get("mechanism_note", ""),
            "diag_note": diag_note,
            "coord_number": exp_rec.get("cn_eff"),
            "valence_electrons": exp_rec.get("n_ve"),
            "sascore": exp_rec.get("sascore"),
            "br_sascore": exp_rec.get("br_sascore"),
            "ligand_sa": exp_rec.get("ligand_sa"),
            "vbur": vbur_best,
            "vbur_d2_7net": exp_rec.get("vbur_d2_7net"),
            "vbur_etkdg": exp_rec.get("vbur_d2_percent"),
            "bde_min": exp_rec.get("bde_min_kcal_mol"),
            "e_coh": exp_rec.get("e_coh_ev_per_atom"),
            "de_dimer": exp_rec.get("delta_e_dimer_kcal_mol"),
            "dg_dimer_373": exp_rec.get("delta_g_dimer_373_kcal_mol"),
            "de_form_bond": exp_rec.get("delta_e_form_per_bond_kcal_mol"),
            "tmc_sa_no3d": exp_rec.get("tmc_sa_no3d"),
            "tmc_sa_phys": tmc_best,
            "phase2_pass": pass_best,
            "domain_flag": "boundary_diagnostic" if exp_rec.get("boundary_flag") else "patent_lit_pass",
        })

    _THD = "O/C(=C/C(=O)C(C)(C)C)C(C)(C)C"
    _ACAC = "O/C(C)=C\\C(C)=O"
    _HFAC = "O/C(C(F)(F)F)=C\\C(=O)C(F)(F)F"
    _FOD = "O/C(C(C)(C)C)=C\\C(=O)C(F)(F)C(F)(F)C(F)(F)F"
    _OD = "O/C(CC)=C\\C(=O)CC"
    _DMAMP = "OC(C)(C)CN(C)C"
    _DMAMB = "OC(C)(CC)CN(C)C"
    _CP = "C1C=CC=C1"
    _ETCP = "CCC1=CC=CC1"
    _IPRCP = "CC(C)C1=CC=CC1"
    _CP_IPR3 = "CC(C)C1=C(C(C)C)C=C(C(C)C)C1"
    _CURATED_DISPLAY_2D_SMILES: Dict[str, str] = {
        "LiOtBu": "[Li]OC(C)(C)C",
        "La(thd)3": f"[La]({_THD})({_THD}){_THD}",
        "Ta(OEt)5": "CCO[Ta](OCC)(OCC)(OCC)OCC",
        "Te(SiMe3)2": "C[Si](C)(C)[Te][Si](C)(C)C",
        "Mn(thd)3": f"[Mn]({_THD})({_THD}){_THD}",
        "Hf(NEt2)4": "CCN(CC)[Hf](N(CC)CC)(N(CC)CC)N(CC)CC",
        "LiN(SiMe3)2": "[Li]N([Si](C)(C)C)[Si](C)(C)C",
        "FeCp2": f"[Fe]({_CP}){_CP}",
        "Sn(dmamp)2": f"[Sn]({_DMAMP}){_DMAMP}",
        "Ca(thd)2": f"[Ca]({_THD}){_THD}",
        "Sr(thd)2": f"[Sr]({_THD}){_THD}",
        "WF6": "F[W](F)(F)(F)(F)F",
        "Y(thd)3": f"[Y]({_THD})({_THD}){_THD}",
        "Ag(fod)(PEt3)": f"CCP(CC)CC->[Ag]{_FOD}",
        "Si2H6": "[SiH3][SiH3]",
        "Te(OEt)4": "CCO[Te](OCC)(OCC)OCC",
        "CoCp2": f"[Co]({_CP}){_CP}",
        "Ni(acac)2": f"[Ni]({_ACAC}){_ACAC}",
        "NiCp2": f"[Ni]({_CP}){_CP}",
        "Yb(thd)3": f"[Yb]({_THD})({_THD}){_THD}",
        "Cu(acac)2": f"[Cu]({_ACAC}){_ACAC}",
        "NaOtBu": "[Na]OC(C)(C)C",
        "Te(SiEt3)2": "CC[Si](CC)(CC)[Te][Si](CC)(CC)CC",
        "Er(thd)3": f"[Er]({_THD})({_THD}){_THD}",
        "Se(SiMe3)2": "C[Si](C)(C)[Se][Si](C)(C)C",
        "SnCl4": "Cl[Sn](Cl)(Cl)Cl",
        "Sr(CpiPr3)2": f"[Sr]({_CP_IPR3}){_CP_IPR3}",
        "TiF4": "F[Ti](F)(F)F",
        "Fe(thd)3": f"[Fe]({_THD})({_THD}){_THD}",
        "Pb(N(SiMe3)2)2": "C[Si](C)(C)N([Si](C)(C)C)[Pb]N([Si](C)(C)C)[Si](C)(C)C",
        "Sb(OEt)3": "CCO[Sb](OCC)OCC",
        "Sc(thd)3": f"[Sc]({_THD})({_THD}){_THD}",
        "SnI4": "I[Sn](I)(I)I",
        "Ti(OMe)4": "CO[Ti](OC)(OC)OC",
        "V(NEtMe)4": "CCN(C)[V](N(C)CC)(N(C)CC)N(C)CC",
        "ZnMe2": "C[Zn]C",
        "[Cu(sBuAMD)]2": "CCC(C)N=C(C)N(C(C)CC)[Cu]",
        "Co(thd)2": f"[Co]({_THD}){_THD}",
        "KOtBu": "[K]OC(C)(C)C",
        "Ni(dmamb)2": f"[Ni]({_DMAMB}){_DMAMB}",
        "Ni(tBu2AMD)2": "CC(=NC(C)(C)C)N(C(C)(C)C)[Ni]N(C(C)(C)C)C(C)=NC(C)(C)C",
        "Si(OtBu)3OH": "CC(C)(C)O[Si](O)(OC(C)(C)C)OC(C)(C)C",
        "TaF5": "F[Ta](F)(F)(F)F",
        "BH3(NHMe2)": "C[NH+](C)[BH3-]",
        "CdMe2": "C[Cd]C",
        "Cu(hfac)2": f"[Cu]({_HFAC}){_HFAC}",
        "La(iPrFMD)3": "CC(C)N=CN(C(C)C)[La](N(C(C)C)C=NC(C)C)N(C(C)C)C=NC(C)C",
        "Mg(thd)2": f"[Mg]({_THD}){_THD}",
        "Pb(dmamp)2": f"[Pb]({_DMAMP}){_DMAMP}",
        "Pb(thd)2": f"[Pb]({_THD}){_THD}",
        "Pd(hfac)2": f"[Pd]({_HFAC}){_HFAC}",
        "Sb(NMe2)3": "CN(C)[Sb](N(C)C)N(C)C",
        "Se(SiEt3)2": "CC[Si](CC)(CC)[Se][Si](CC)(CC)CC",
        "VO(OiPr)3": "CC(C)O[V](=O)(OC(C)C)OC(C)C",
        "ZnCl2": "Cl[Zn]Cl",
        "Zr(OtBu)4": "CC(C)(C)O[Zr](OC(C)(C)C)(OC(C)(C)C)OC(C)(C)C",
        "Cu(thd)2": f"[Cu]({_THD}){_THD}",
        "Gd(thd)3": f"[Gd]({_THD})({_THD}){_THD}",
        "Ge(OEt)4": "CCO[Ge](OCC)(OCC)OCC",
        "La(CpiPr)3": f"[La]({_IPRCP})({_IPRCP}){_IPRCP}",
        "Mg(EtCp)2": f"[Mg]({_ETCP}){_ETCP}",
        "Mn(EtCp)2": f"[Mn]({_ETCP}){_ETCP}",
        "SbCl3": "Cl[Sb](Cl)Cl",
        "SnEt4": "CC[Sn](CC)(CC)CC",
        "GaNp3": "CC(C)(C)C[Ga](CC(C)(C)C)CC(C)(C)C",
        "In(triaz)3": "CN=NN(C)[In](N(C)N=NC)N(C)N=NC",
        "Ga(triaz)3": "CN=NN(C)[Ga](N(C)N=NC)N(C)N=NC",
        "In(dpguan)3": "CN(C)C(=NC(C)C)N(C(C)C)[In](N(C(C)C)C(=NC(C)C)N(C)C)N(C(C)C)C(=NC(C)C)N(C)C",
        "InMe2(NtBu(CH2)2(CH3)2OMe)": "C[In](C)N(C(C)(C)C)CC(C)(C)OC",
        "Ru(CpEtMe)2": "CCC1=C(C)C=CC1[Ru]C1C=CC(C)=C1CC",
        "Ru(Me2pyr)2": "CC1=CC(C)=C([Ru]C2=C(C)C=C(C)N2)N1",
        "Ru(od)3": f"[Ru]({_OD})({_OD}){_OD}",
        "Ru(thd)2(cod)": f"[Ru]({_THD})({_THD})C1=CCCC=CCC1",
        "RuCp(CO)2((CH2)3N(CH3)2)": f"CN(C)CCC[Ru]({_CP})(<-[C-]#[O+])<-[C-]#[O+]",
        "RuCp(CpCH(Me)(NMe2))": f"CC(N(C)C)C1=CC=CC1[Ru]{_CP}",
    }

    gate_reason_labels: Dict[str, str] = {
        "sa_gt_4": "합성 난이도 초과 (>4.0)",
        "homoleptic_tm_sp3": "불안정 금속 알킬 (β-H 제거/열분해)",
        "overcrowded_poly_aryl_tm": "다중 아릴 입체 과밀",
        "poly_hydroxide": "다중 수산화기·응집 위험",
        "carboxylic_acid": "카복실산 부반응·비휘발",
        "inert_alkane": "비활성 알케인 (배위 불가)",
        "strained_si_ring": "고리 변형 실란 불안정",
        "pi_conjugated_main_group_hydride": "공액 하이드라이드 불안정",
        "catenated_main_group_hydride": "연쇄 하이드라이드 불안정",
    }

    ens_by_id: Dict[str, Dict[str, Any]] = {}
    ens_path = os.path.join("database", "external", "umlip_ensemble_198.json")
    if os.path.exists(ens_path):
        with open(ens_path, "r", encoding="utf-8") as fh:
            ens_by_id = json.load(fh).get("records", {})

    tdec_by_id: Dict[str, Dict[str, Any]] = {}
    tdec_path = os.path.join("database", "external", "thermal_decomposition_db.json")
    if os.path.exists(tdec_path):
        with open(tdec_path, "r", encoding="utf-8") as fh:
            tdec_by_id = json.load(fh).get("records", {})

    for rec in db_records:
        # Backfill display-only 2D RDKit structure for known AtomicLimits precursors
        if not rec.get("has_2d"):
            disp_smi = _CURATED_DISPLAY_2D_SMILES.get(rec.get("name", ""))
            if disp_smi and _ensure_svg(disp_smi):
                rec["smiles"] = disp_smi
                rec["has_2d"] = True
                m_disp = Chem.MolFromSmiles(disp_smi)
                if m_disp is not None:
                    if not rec.get("formula") or rec.get("formula") == "-":
                        rec["formula"] = rdMolDescriptors.CalcMolFormula(m_disp)
                    if rec.get("mw") is None:
                        rec["mw"] = round(Descriptors.MolWt(m_disp), 2)

        u_rec = _u_by_id.get(rec.get("id"))
        if u_rec is not None:
            if u_rec.get("tmc_sa_v2") is not None:
                rec["tmc_sa_phys"] = round(float(u_rec["tmc_sa_v2"]), 2)
            if u_rec.get("unified_pass") is not None:
                rec["phase2_pass"] = bool(u_rec["unified_pass"])
            rec["gate_reasons"] = list(u_rec.get("gate_reasons") or [])

        sapis_val = rec.get("tmc_sa_phys")
        if sapis_val is None:
            sapis_val = rec.get("tmc_sa_no3d")
        rec["sapis_score"] = round(float(sapis_val), 2) if sapis_val is not None else None
        rec["sapis_pass"] = rec.get("phase2_pass")

        # Attach 4-uMLIP robust ensemble results when available
        ens_rec = ens_by_id.get(rec.get("id"))
        if ens_rec is None:
            for sid in rec.get("secondary_ids") or []:
                if sid in ens_by_id:
                    ens_rec = ens_by_id[sid]
                    break
        if ens_rec is not None and ens_rec.get("n_models", 0) >= 2:
            rec["ens_sapis"] = round(float(ens_rec["sapis_median"]), 2)
            rec["ens_min"] = round(float(ens_rec["sapis_min"]), 2)
            rec["ens_max"] = round(float(ens_rec["sapis_max"]), 2)
            rec["ens_sigma_bde"] = round(float(ens_rec["sigma_bde_kcal_mol"]), 1)
            rec["ens_iqr_bde"] = round(float(ens_rec["iqr_bde_kcal_mol"]), 1)
            rec["ens_grade"] = str(ens_rec["confidence_grade"])
            rec["ens_dft_flag"] = bool(ens_rec["dft_verification_flag"])
            rec["ens_agree"] = str(ens_rec.get("verdict_agreement") or "")
            rec["ens_spread"] = round(float(ens_rec.get("sapis_spread", float(ens_rec["sapis_max"]) - float(ens_rec["sapis_min"]))), 2)
            rec["ens_per_model"] = ens_rec.get("per_model_sapis") or {}
            rec["ens_bde_med"] = ens_rec.get("bde_min_median_kcal_mol")
            rec["ens_n_models"] = int(ens_rec["n_models"])
            # Use robust 4-uMLIP median as the primary displayed SAPIS score
            rec["sapis_score"] = rec["ens_sapis"]

        # Attach experimental thermal decomposition temperature (T_dec) & molecular mechanism when available
        tdec_rec = tdec_by_id.get(rec.get("id"))
        if tdec_rec is None:
            for sid in rec.get("secondary_ids") or []:
                if sid in tdec_by_id:
                    tdec_rec = tdec_by_id[sid]
                    break
        if tdec_rec is not None:
            rec["t_dec"] = {
                "min_c": tdec_rec.get("t_dec_c_min"),
                "max_c": tdec_rec.get("t_dec_c_max"),
                "display": tdec_rec.get("t_dec_display"),
                "metric": tdec_rec.get("t_dec_metric"),
                "data_type": tdec_rec.get("data_type"),
                "cond": tdec_rec.get("measurement_conditions") or {},
                "ea": tdec_rec.get("activation_energy") or {},
                "mech": tdec_rec.get("molecular_mechanism") or {},
                "refs": tdec_rec.get("references") or [],
            }

        # Uniform concise source label
        doi_str = (rec.get("doi") or "").strip()
        pat_str = (rec.get("patent_references") or "").strip()
        al_cnt = rec.get("atomiclimits_process_count") or 0
        dt = rec.get("db_type")
        if doi_str and al_cnt:
            rec["card_source"] = "DOI · AtomicLimits"
        elif doi_str:
            rec["card_source"] = "DOI 문헌"
        elif pat_str:
            rec["card_source"] = "특허 · DOI" if "DOI:" in pat_str.upper() else "특허 문헌"
        elif al_cnt:
            rec["card_source"] = f"AtomicLimits ({al_cnt}건)"
        elif dt == "inhibitor":
            rec["card_source"] = "AS-ALD 문헌"
        elif dt == "reagent":
            rec["card_source"] = "상용 시약 DB"
        elif dt == "reaction":
            rec["card_source"] = "공인 합성 반응"
        elif dt == "negative":
            rec["card_source"] = "Negative Control"
        else:
            rec["card_source"] = rec.get("db_label") or "-"

        # Detailed risk / rejection reason label
        g_reasons = rec.get("gate_reasons") or []
        diag_note_str = (rec.get("diag_note") or "").strip()
        bde_v = rec.get("bde_min")
        lig_v = rec.get("ligand_sa")
        if u_rec is not None:
            if bde_v is None:
                bde_v = u_rec.get("bde_min_kcal_mol")
            if lig_v is None:
                lig_v = u_rec.get("ligand_sa")

        if rec.get("sapis_pass") is False:
            if diag_note_str:
                rec["card_risk"] = diag_note_str
            elif g_reasons:
                parts = []
                for g in g_reasons[:2]:
                    if g == "sa_gt_4" and rec["sapis_score"] is not None:
                        parts.append(f"합성 난이도 초과 (SAPIS {rec['sapis_score']:.2f} > 4.0)")
                    else:
                        parts.append(gate_reason_labels.get(g, g))
                rec["card_risk"] = ", ".join(parts)
            elif bde_v is not None and float(bde_v) < 20.0 and (lig_v is None or float(lig_v) <= 3.5):
                rec["card_risk"] = f"7net 저결합 해리 페널티 (BDE={float(bde_v):.1f} kcal/mol, ligand_sa={float(lig_v):.2f})"
            elif dt == "inhibitor" and rec["sapis_score"] is not None and lig_v is not None and float(rec["sapis_score"]) > float(lig_v) + 0.5:
                rec["card_risk"] = f"BR-SAScore 사전 미등재 조각 패널티 (BR={rec['sapis_score']:.2f}, ligand_sa={float(lig_v):.2f})"
            elif rec["sapis_score"] is not None and rec["sapis_score"] > 4.0:
                rec["card_risk"] = f"합성·물리 점수 초과 (SAPIS {rec['sapis_score']:.2f} > 4.0)"
            else:
                rec["card_risk"] = "합성·물리 판정 기준 차단"
        elif rec.get("sapis_pass") is True:
            rec["card_risk"] = "낮음 (판정 기준 통과)"
        elif dt == "reagent":
            rec["card_risk"] = rec.get("handling") or "일반 취급"
        elif dt == "reaction":
            rec["card_risk"] = f"ΔG {rec['dg_kcal_mol']} kcal/mol" if rec.get("dg_kcal_mol") is not None else "-"
        else:
            rec["card_risk"] = "미평가 (문헌 참고용)"

    db_meta_summary = {
        "n_unique_entities": len(db_records),
        "n_merged_dual_role_reagents": n_merged_dual_role,
        "n_merged_notation_aliases": len(pending_notation_aliases),
        "n_backfilled_atomiclimits": n_backfilled_al,
        "n_net_new_patent": n_net_new_patent,
        "n_patent_total": len(patent_exp_records),
    }

    db_explorer_payload_json = json.dumps(
        {"svgs": svg_dict, "records": db_records, "meta": db_meta_summary},
        ensure_ascii=False,
        separators=(",", ":"),
    ).replace("</", "<\\/")


    def _build_dataset_sankey_svg() -> str:
        """Build a proportional-width SVG Sankey diagram integrating 304 Baseline + 38 Patent/Lit Expansion = 342 Curated."""
        n_prec_cat = len(expanded_prec_db)
        n_prec_eval = int(_u_ds["precursors_total"])
        n_prec_ref = n_prec_cat - n_prec_eval
        n_inhib = int(_u_ds["inhibitors_total"])
        n_neg = int(_u_ds["negative_controls_total"])
        n_reag_rxn = len(bb_db) + len(rxn_db)
        n_eval = int(_u_ds["total_molecules"])
        n_excl = n_prec_ref + n_reag_rxn
        n_svn = int(_u_ds["sevenn_relaxed_count"])
        n_etk_base = int(_u_ds["etkdg_3d_count"])
        n_pos = int(_u_ds["positive_total"])

        pat_s = patent_exp_payload.get("summary", {})
        n_pat_tot = int(pat_s.get("n_total", len(patent_exp_records)))
        n_pat_prec = int(pat_s.get("n_precursors", 0))
        n_pat_inh = int(pat_s.get("n_inhibitors", 0))
        n_pat_7n = int(pat_s.get("n_sevenn_computed", 0))
        n_pat_7n_pass = int(pat_s.get("n_precursors_pass_7net", 0)) + int(pat_s.get("n_inhibitors_pass_7net", 0))
        n_pat_no3d_pass = int(pat_s.get("n_precursors_pass", 0)) + int(pat_s.get("n_inhibitors_pass", 0))
        n_pat_pass = n_pat_7n_pass if n_pat_7n == n_pat_tot and n_pat_tot > 0 else n_pat_no3d_pass
        n_pat_bnd = n_pat_tot - n_pat_pass
        n_total_curated = n_prec_cat + n_inhib + n_neg + n_reag_rxn + n_pat_tot

        n_svn_pos = sum(1 for r in _u_by_id.values() if r["label"] == 1 and "SevenNet" in str(r.get("vbur_source", "")))
        n_svn_neg = sum(1 for r in _u_by_id.values() if r["label"] == 0 and "SevenNet" in str(r.get("vbur_source", "")))
        n_etk_pos = n_pos - n_svn_pos
        n_etk_neg = n_neg - n_svn_neg

        scale = 0.70
        bar_w = 14.0
        x0, x1, x2, x3 = 205.0, 415.0, 605.0, 760.0

        def _layout(items, y_start=56.0, gap=14.0):
            pos = {}
            y = y_start
            for key, val, color, title in items:
                h = max(16.0, val * scale)
                pos[key] = {"x": 0.0, "y": y, "h": h, "val": val, "color": color, "title": title, "in_off": y, "out_off": y}
                y += h + gap
            return pos

        c0 = _layout([
            ("cat_prec", n_prec_cat, "#06b6d4", "프리커서 카탈로그"),
            ("cat_inh", n_inhib, "#0284c7", "선택적 증착 인히비터"),
            ("cat_pat", n_pat_tot, "#f59e0b", "특허·문헌 확장 후보"),
            ("cat_neg", n_neg, "#ef4444", "합성 불가 대조군"),
            ("cat_rr", n_reag_rxn, "#64748b", "출발 시약 및 합성 반응"),
        ], y_start=56.0, gap=12.0)
        for d in c0.values():
            d["x"] = x0

        y_eval = 56.0
        y_pat1 = y_eval + n_eval * scale + 24.0
        y_ref1 = y_pat1 + n_pat_tot * scale + 28.0
        c1 = {
            "eval_160": {"x": x1, "y": y_eval, "h": n_eval * scale, "val": n_eval, "color": "#0f766e", "title": "통합 분류 평가셋", "in_off": y_eval, "out_off": y_eval},
            "pat_38": {"x": x1, "y": y_pat1, "h": n_pat_tot * scale, "val": n_pat_tot, "color": "#d97706", "title": "특허·문헌 스크리닝셋", "in_off": y_pat1, "out_off": y_pat1},
            "ref_144": {"x": x1, "y": y_ref1, "h": n_excl * scale, "val": n_excl, "color": "#94a3b8", "title": "참조·사전 전용 데이터", "in_off": y_ref1, "out_off": y_ref1},
        }

        c2_items = [("svn", n_svn, "#10b981", "벤치마크 3D 최적화")]
        if n_etk_base > 0:
            c2_items.append(("etk_base", n_etk_base, "#8b5cf6", "초기 3D 구조 생성"))
        c2_items.append(("etk_pat", n_pat_tot, "#059669" if (n_pat_7n == n_pat_tot and n_pat_tot > 0) else "#f59e0b", "확장 후보 3D 최적화"))
        c2 = _layout(c2_items, y_start=56.0, gap=24.0)
        for d in c2.values():
            d["x"] = x2

        c3 = _layout([
            ("pos_120", n_pos, "#16a34a", "합성 가능 양성군"),
            ("neg_40", n_neg, "#dc2626", "합성 불가 대조군"),
            ("pat_out", n_pat_tot, "#0284c7", "확장 후보 SAPIS 판정"),
        ], y_start=56.0, gap=24.0)
        for d in c3.values():
            d["x"] = x3

        flows = [
            (c0["cat_prec"], c1["eval_160"], n_prec_eval, "#06b6d4", f"프리커서 카탈로그 → 통합 분류 평가셋: {n_prec_eval}종"),
            (c0["cat_inh"], c1["eval_160"], n_inhib, "#0284c7", f"선택적 증착 인히비터 → 통합 분류 평가셋: {n_inhib}종"),
            (c0["cat_neg"], c1["eval_160"], n_neg, "#ef4444", f"합성 불가 대조군 → 통합 분류 평가셋: {n_neg}종"),
            (c0["cat_pat"], c1["pat_38"], n_pat_tot, "#f59e0b", f"특허·문헌 확장 후보 → 특허·문헌 스크리닝셋: {n_pat_tot}종"),
            (c0["cat_prec"], c1["ref_144"], n_prec_ref, "#94a3b8", f"프리커서 카탈로그 → 참조·사전 전용 데이터: {n_prec_ref}건"),
            (c0["cat_rr"], c1["ref_144"], n_reag_rxn, "#64748b", f"출발 시약 및 합성 반응 → 참조·사전 전용 데이터: {n_reag_rxn}건"),
            (c1["eval_160"], c2["svn"], n_svn, "#10b981", f"통합 분류 평가셋 → 벤치마크 3D 최적화: {n_svn}종"),
        ]
        if n_etk_base > 0:
            flows.append((c1["eval_160"], c2["etk_base"], n_etk_base, "#8b5cf6", f"통합 분류 평가셋 → 초기 3D 구조 생성: {n_etk_base}종"))
        flows.extend([
            (c1["pat_38"], c2["etk_pat"], n_pat_tot, "#059669" if n_pat_7n == n_pat_tot else "#f59e0b", f"특허·문헌 스크리닝셋 → 확장 후보 3D 최적화: {n_pat_tot}종"),
            (c2["svn"], c3["pos_120"], n_svn_pos, "#16a34a", f"벤치마크 3D 최적화 → 합성 가능 양성군: {n_svn_pos}종"),
            (c2["svn"], c3["neg_40"], n_svn_neg, "#059669", f"벤치마크 3D 최적화 → 합성 불가 대조군: {n_svn_neg}종"),
            (c2["etk_pat"], c3["pat_out"], n_pat_tot, "#0284c7", f"확장 후보 3D 최적화 → 확장 후보 SAPIS 판정: {n_pat_tot}종"),
        ])

        ribbon_svg = []
        for src, dst, val, color, tip in flows:
            dy = val * scale
            y0a = src["out_off"]
            y0b = y0a + dy
            y1a = dst["in_off"]
            y1b = y1a + dy
            src["out_off"] = y0b
            dst["in_off"] = y1b
            sx = src["x"] + bar_w
            dx = dst["x"]
            mx = 0.5 * (sx + dx)
            d_path = (
                f"M {sx:.1f},{y0a:.1f} "
                f"C {mx:.1f},{y0a:.1f} {mx:.1f},{y1a:.1f} {dx:.1f},{y1a:.1f} "
                f"L {dx:.1f},{y1b:.1f} "
                f"C {mx:.1f},{y1b:.1f} {mx:.1f},{y0b:.1f} {sx:.1f},{y0b:.1f} Z"
            )
            ribbon_svg.append(
                f'<path d="{d_path}" fill="{color}" class="sankey-flow"><title>{tip}</title></path>'
            )

        nodes_svg = []
        for d in c0.values():
            ym = d["y"] + d["h"] * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]} ({d["val"]}건)</title></rect>'
                f'<text x="{d["x"] - 10:.1f}" y="{ym + 4:.1f}" text-anchor="end" font-size="11.2" font-weight="700" fill="#0f172a">{d["title"]} ({d["val"]})</text>'
            )

        for k, d in c1.items():
            ym = d["y"] + d["h"] * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]} ({d["val"]}건)</title></rect>'
            )
            if k in {"eval_160", "pat_38"}:
                bw = 152 if k == "pat_38" else 136
                bx = d["x"] + 7 - bw * 0.5
                tcol = "#92400e" if k == "pat_38" else "#0f766e"
                nodes_svg.append(
                    f'<rect x="{bx:.1f}" y="{ym - 12:.1f}" width="{bw}" height="24" rx="5" fill="#ffffff" fill-opacity="0.95" stroke="#cbd5e1" stroke-width="0.9" />'
                    f'<text x="{d["x"] + 7:.1f}" y="{ym + 4:.1f}" text-anchor="middle" font-size="10.8" font-weight="700" fill="{tcol}">{d["title"]} ({d["val"]})</text>'
                )
            else:
                nodes_svg.append(
                    f'<text x="{d["x"] + bar_w + 10:.1f}" y="{ym + 4:.1f}" text-anchor="start" font-size="11.2" font-weight="700" fill="#475569">{d["title"]} ({d["val"]})</text>'
                )

        for d in c2.values():
            ym = d["y"] + d["h"] * 0.5
            bw = 146
            bx = d["x"] + 7 - bw * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]} ({d["val"]}종)</title></rect>'
                f'<rect x="{bx:.1f}" y="{ym - 12:.1f}" width="{bw}" height="24" rx="5" fill="#ffffff" fill-opacity="0.95" stroke="#cbd5e1" stroke-width="0.9" />'
                f'<text x="{d["x"] + 7:.1f}" y="{ym + 4:.1f}" text-anchor="middle" font-size="10.8" font-weight="700" fill="#0f172a">{d["title"]} ({d["val"]})</text>'
            )

        for d in c3.values():
            ym = d["y"] + d["h"] * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]} ({d["val"]}종)</title></rect>'
                f'<text x="{d["x"] + bar_w + 10:.1f}" y="{ym + 4:.1f}" text-anchor="start" font-size="11.2" font-weight="700" fill="#0f172a">{d["title"]} ({d["val"]})</text>'
            )

        headers_svg = (
            f'<text x="{x0 + 7:.1f}" y="25" text-anchor="middle" font-size="11.5" font-weight="800" fill="#0369a1">1. 원본 데이터 수집 ({n_total_curated})</text>'
            f'<text x="{x1 + 7:.1f}" y="25" text-anchor="middle" font-size="11.5" font-weight="800" fill="#0f766e">2. 평가·스크리닝 분할 ({n_total_curated})</text>'
            f'<text x="{x2 + 7:.1f}" y="25" text-anchor="middle" font-size="11.5" font-weight="800" fill="#059669">3. 3D 물리 구조 최적화 ({n_eval + n_pat_tot})</text>'
            f'<text x="{x3 + 7:.1f}" y="25" text-anchor="middle" font-size="11.5" font-weight="800" fill="#334155">4. 최종 검증 및 판정 ({n_eval + n_pat_tot})</text>'
            f'<line x1="18" y1="38" x2="942" y2="38" stroke="#e2e8f0" stroke-width="1" />'
        )

        return f"""
              <div style="margin-bottom:14px;">
                <style>.sankey-flow {{ fill-opacity: 0.42; transition: fill-opacity 0.15s ease; cursor: pointer; }} .sankey-flow:hover {{ fill-opacity: 0.72; }}</style>
                <svg viewBox="0 0 960 368" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;">
                  {headers_svg}
                  {''.join(ribbon_svg)}
                  {''.join(nodes_svg)}
                </svg>
              </div>"""


    dataset_sankey_html = _build_dataset_sankey_svg()

    return {
        "<!-- DB_EXPLORER_DATA_PLACEHOLDER -->": db_explorer_payload_json,
        "<!-- DATASET_SANKEY_PLACEHOLDER -->": dataset_sankey_html,
    }
