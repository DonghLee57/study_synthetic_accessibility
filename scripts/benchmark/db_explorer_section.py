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
            f"All SevenNet-0 & TMC-SA physical/chemical gates passed (BDE={exp_rec.get('bde_min_kcal_mol')} kcal/mol, TMC={exp_rec.get('tmc_sa_7net')})"
            if exp_rec.get("pass_7net")
            else ("All ETKDG TMC-SA physical & chemical gates passed" if exp_rec.get("unified_pass") else "")
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

        scale = 0.72
        bar_w = 14.0
        x0, x1, x2, x3 = 198.0, 412.0, 602.0, 772.0

        def _layout(items, y_start=54.0, gap=12.0):
            pos = {}
            y = y_start
            for key, val, color, title, sub in items:
                h = val * scale
                pos[key] = {"x": 0.0, "y": y, "h": h, "val": val, "color": color, "title": title, "sub": sub, "in_off": y, "out_off": y}
                y += h + gap
            return pos

        c0 = _layout([
            ("cat_prec", n_prec_cat, "#06b6d4", "Precursor Catalog", "Core 27 · AL 93 · Exp 85"),
            ("cat_inh", n_inhib, "#0284c7", "AS-ALD Inhibitors", "Baseline SMILES 18/18"),
            ("cat_pat", n_pat_tot, "#f59e0b", "Patent & Lit Expansion", f"+{n_pat_prec} Prec · +{n_pat_inh} Inh ({n_backfilled_al} AL-BF + {n_net_new_patent} New)"),
            ("cat_neg", n_neg, "#ef4444", "Negative Controls", "Core 10 · 6-Mech 30"),
            ("cat_rr", n_reag_rxn, "#64748b", "Reagents & Reactions", f"{len(bb_db)} Reagents (10 Dual-Role) + {len(rxn_db)} Rxns"),
        ], y_start=54.0, gap=11.0)
        for d in c0.values():
            d["x"] = x0

        y_eval = 54.0
        y_pat1 = y_eval + n_eval * scale + 20.0
        y_ref1 = y_pat1 + n_pat_tot * scale + 28.0
        c1 = {
            "eval_160": {"x": x1, "y": y_eval, "h": n_eval * scale, "val": n_eval, "color": "#0f766e", "title": "Pinned Benchmark", "sub": f"SMILES Valid ({n_pos}:{n_neg})", "in_off": y_eval, "out_off": y_eval},
            "pat_38": {"x": x1, "y": y_pat1, "h": n_pat_tot * scale, "val": n_pat_tot, "color": "#d97706", "title": "Patent Screening", "sub": f"+8 Metals (16→{pat_s.get('n_merged_elements', 24)}) · Deduplicated", "in_off": y_pat1, "out_off": y_pat1},
            "ref_144": {"x": x1, "y": y_ref1, "h": n_excl * scale, "val": n_excl, "color": "#94a3b8", "title": "Reference / Dict Only", "sub": f"AL {n_prec_ref - 10} (incl. {n_backfilled_al} BF) · Meta 10 · Dict {len(bb_db)} · Rxn {len(rxn_db)}", "in_off": y_ref1, "out_off": y_ref1},
        }

        c2_items = [("svn", n_svn, "#10b981", "SevenNet LBFGS 3D", "Core 55 · Exp 75 · Neg 30 · 100%")]
        if n_etk_base > 0:
            c2_items.append(("etk_base", n_etk_base, "#8b5cf6", "ETKDGv3 Guarded", "Benchmark Subset"))
        if n_pat_7n == n_pat_tot and n_pat_tot > 0:
            c2_items.append(("etk_pat", n_pat_tot, "#059669", "Exp SevenNet LBFGS", f"{n_pat_prec} Prec + {n_pat_inh} Inh · 0 Fallback"))
        else:
            c2_items.append(("etk_pat", n_pat_tot, "#f59e0b", "ETKDGv3 3D + TMC-SA", f"{n_pat_prec} Prec + {n_pat_inh} Inh · 100% Embed"))
        c2 = _layout(c2_items, y_start=54.0, gap=22.0)
        for d in c2.values():
            d["x"] = x2

        pat_out_sub = (
            f"7net {n_pat_pass}/{n_pat_tot} Pass (PDMAT/Ge Rescued) · {n_pat_bnd} Diag"
            if (n_pat_7n == n_pat_tot and n_pat_tot > 0)
            else f"{n_pat_pass} Pass ({int(pat_s.get('n_precursors_pass', 0))}+{int(pat_s.get('n_inhibitors_pass', 0))}) · {n_pat_bnd} Boundary Flag"
        )
        c3 = _layout([
            ("pos_120", n_pos, "#16a34a", "Benchmark Positive", f"{n_prec_eval}+{n_inhib} · OCSVM · CV-5"),
            ("neg_40", n_neg, "#dc2626", "Negative Control", "LOFMO Test · CV-5 · TMC-SA"),
            ("pat_out", n_pat_tot, "#0284c7", "Expansion 7net TMC-SA", pat_out_sub),
        ], y_start=54.0, gap=22.0)
        for d in c3.values():
            d["x"] = x3

        flows = [
            (c0["cat_prec"], c1["eval_160"], n_prec_eval, "#06b6d4", f"Precursor Catalog -> Pinned Eval Set: {n_prec_eval} (27 Core + 75 Expanded)"),
            (c0["cat_inh"], c1["eval_160"], n_inhib, "#0284c7", f"AS-ALD Inhibitors -> Pinned Eval Set: {n_inhib}"),
            (c0["cat_neg"], c1["eval_160"], n_neg, "#ef4444", f"Negative Controls -> Pinned Eval Set: {n_neg} (10 Core + 30 Stratified)"),
            (c0["cat_pat"], c1["pat_38"], n_pat_tot, "#f59e0b", f"KR/JP/CN/US/TW Patent & Lit Expansion -> TMC-SA Screen & DB Backfill: {n_pat_tot} ({n_backfilled_al} AL Backfill + {n_net_new_patent} Net-New)"),
            (c0["cat_prec"], c1["ref_144"], n_prec_ref, "#94a3b8", f"Precursor Catalog -> Reference Only: {n_prec_ref} (AtomicLimits + Meta 10; {n_backfilled_al} stubs SMILES-backfilled in DB Explorer)"),
            (c0["cat_rr"], c1["ref_144"], n_reag_rxn, "#64748b", f"Reagents ({len(bb_db)}, 10 Dual-Role Merged in DB Search) + Reactions ({len(rxn_db)}) -> Dict / Ref: {n_reag_rxn}"),
            (c1["eval_160"], c2["svn"], n_svn, "#10b981", f"Pinned Eval Set -> SevenNet LBFGS 3D: {n_svn} ({n_svn_pos} Pos + {n_svn_neg} Neg)"),
        ]
        if n_etk_base > 0:
            flows.append((c1["eval_160"], c2["etk_base"], n_etk_base, "#8b5cf6", f"Eval -> ETKDGv3: {n_etk_base}"))
        flows.extend([
            (c1["pat_38"], c2["etk_pat"], n_pat_tot, "#059669" if n_pat_7n == n_pat_tot else "#f59e0b", f"Patent Expansion -> SevenNet-0 LBFGS 3D + BDE/Dimer/ΔEf: {n_pat_tot} (0 Fallbacks)"),
            (c2["svn"], c3["pos_120"], n_svn_pos, "#16a34a", f"SevenNet LBFGS -> Benchmark Positive: {n_svn_pos}"),
            (c2["svn"], c3["neg_40"], n_svn_neg, "#059669", f"SevenNet LBFGS -> Negative Control Validation: {n_svn_neg}"),
            (c2["etk_pat"], c3["pat_out"], n_pat_tot, "#0284c7", f"Expansion 7net TMC-SA -> {n_pat_pass}/{n_pat_tot} Strict Pass (ETKDG {n_pat_no3d_pass}/{n_pat_tot}) + {n_pat_bnd} Boundary Diagnostics"),
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
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]}: {d["val"]}</title></rect>'
            )
            nodes_svg.append(
                f'<text x="{d["x"] - 8:.1f}" y="{ym - 1:.1f}" text-anchor="end" font-size="11" font-weight="800" fill="#0f172a">{d["val"]} {d["title"]}</text>'
                f'<text x="{d["x"] - 8:.1f}" y="{ym + 11:.1f}" text-anchor="end" font-size="9" fill="#64748b">{d["sub"]}</text>'
            )

        for k, d in c1.items():
            ym = d["y"] + d["h"] * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]}: {d["val"]}</title></rect>'
            )
            if k in {"eval_160", "pat_38"}:
                bw = 158 if k == "pat_38" else 142
                bx = d["x"] + 7 - bw * 0.5
                tcol = "#92400e" if k == "pat_38" else "#0f766e"
                nodes_svg.append(
                    f'<rect x="{bx:.1f}" y="{ym - 16:.1f}" width="{bw}" height="31" rx="5" fill="#ffffff" fill-opacity="0.93" stroke="#cbd5e1" stroke-width="0.8" />'
                    f'<text x="{d["x"] + 7:.1f}" y="{ym - 2:.1f}" text-anchor="middle" font-size="10.8" font-weight="800" fill="{tcol}">{d["val"]} {d["title"]}</text>'
                    f'<text x="{d["x"] + 7:.1f}" y="{ym + 10:.1f}" text-anchor="middle" font-size="8.6" fill="#475569">{d["sub"]}</text>'
                )
            else:
                nodes_svg.append(
                    f'<text x="{d["x"] + bar_w + 10:.1f}" y="{ym - 2:.1f}" text-anchor="start" font-size="11" font-weight="800" fill="#475569">{d["val"]} {d["title"]}</text>'
                    f'<text x="{d["x"] + bar_w + 10:.1f}" y="{ym + 11:.1f}" text-anchor="start" font-size="9" fill="#64748b">{d["sub"]} (분류 제외)</text>'
                )

        for d in c2.values():
            ym = d["y"] + d["h"] * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]}: {d["val"]}</title></rect>'
                f'<rect x="{d["x"] - 64:.1f}" y="{ym - 16:.1f}" width="142" height="31" rx="5" fill="#ffffff" fill-opacity="0.93" stroke="#cbd5e1" stroke-width="0.8" />'
                f'<text x="{d["x"] + 7:.1f}" y="{ym - 2:.1f}" text-anchor="middle" font-size="10.6" font-weight="800" fill="#0f172a">{d["val"]} {d["title"]}</text>'
                f'<text x="{d["x"] + 7:.1f}" y="{ym + 10:.1f}" text-anchor="middle" font-size="8.6" fill="#475569">{d["sub"]}</text>'
            )

        for d in c3.values():
            ym = d["y"] + d["h"] * 0.5
            nodes_svg.append(
                f'<rect x="{d["x"]:.1f}" y="{d["y"]:.1f}" width="{bar_w}" height="{d["h"]:.1f}" rx="2" fill="{d["color"]}"><title>{d["title"]}: {d["val"]}</title></rect>'
                f'<text x="{d["x"] + bar_w + 8:.1f}" y="{ym - 1:.1f}" text-anchor="start" font-size="11" font-weight="800" fill="#0f172a">{d["val"]} {d["title"]}</text>'
                f'<text x="{d["x"] + bar_w + 8:.1f}" y="{ym + 11:.1f}" text-anchor="start" font-size="8.8" fill="#475569">{d["sub"]}</text>'
            )

        headers_svg = (
            f'<text x="{x0 + 7:.1f}" y="24" text-anchor="middle" font-size="11" font-weight="800" fill="#0369a1">1. Curated Sources ({n_total_curated})</text>'
            f'<text x="{x1 + 7:.1f}" y="24" text-anchor="middle" font-size="11" font-weight="800" fill="#0f766e">2. Benchmark &amp; Patent Split ({n_total_curated})</text>'
            f'<text x="{x2 + 7:.1f}" y="24" text-anchor="middle" font-size="11" font-weight="800" fill="#059669">3. 3D Structure Engine ({n_eval + n_pat_tot})</text>'
            f'<text x="{x3 + 7:.1f}" y="24" text-anchor="middle" font-size="11" font-weight="800" fill="#334155">4. Eval &amp; Screening ({n_eval + n_pat_tot})</text>'
            f'<line x1="18" y1="36" x2="942" y2="36" stroke="#e2e8f0" stroke-width="1" />'
        )

        return f"""
              <div style="margin-bottom:14px;">
                <style>.sankey-flow {{ fill-opacity: 0.44; transition: fill-opacity 0.15s ease; cursor: pointer; }} .sankey-flow:hover {{ fill-opacity: 0.72; }}</style>
                <svg viewBox="0 0 960 392" style="width:100%; height:auto; display:block; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;">
                  {headers_svg}
                  {''.join(ribbon_svg)}
                  {''.join(nodes_svg)}
                  <text x="944" y="382" text-anchor="end" font-size="9.3" fill="#64748b">* Raw Curated {n_total_curated} ({n_eval + n_excl} Baseline + {n_pat_tot} Patent/Lit) &bull; DB Explorer Canonical Unique = {len(db_records)} ({n_merged_dual_role} Dual-Role Reagents &amp; {n_backfilled_al} AL Backfills Merged)</text>
                </svg>
              </div>"""


    dataset_sankey_html = _build_dataset_sankey_svg()

    return {
        "<!-- DB_EXPLORER_DATA_PLACEHOLDER -->": db_explorer_payload_json,
        "<!-- DATASET_SANKEY_PLACEHOLDER -->": dataset_sankey_html,
    }
