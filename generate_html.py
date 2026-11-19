"""Top-level HTML report orchestrator for Synthetic Accessibility & Precursor DB.

Assembles the four modular view tabs:
  1. ``#tab-content-storyline``: Executive Summary (Steps 1-5)
  2. ``#tab-content-technical``: Deep Technical Report (Sections 1-22)
  3. ``#tab-content-database``: Unified 286-Entity 2D Structure DB Explorer
  4. ``#tab-content-brtheory``: Foundational Theory (2009 SAScore & 2024 BR-SAScore)

Applies lossless SVG/HTML compression via ``scripts/benchmark/html_common.py``
and writes ``index.html`` and ``synthetic_accessibility_report.html``.
"""

from __future__ import annotations

import os
import sys
from typing import Dict

BENCH_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts", "benchmark")
if BENCH_DIR not in sys.path:
    sys.path.insert(0, BENCH_DIR)

from br_sascore_theory_section import build_br_sascore_theory_html
from glossary_section import build_glossary_html
from intro_section import build_intro_html
from db_explorer_section import build_db_explorer_and_sankey
from history_and_expansion_html_section import (
    build_analysis_history_technical_html,
    build_patent_expansion_storyline_html,
    build_patent_expansion_technical_html,
)
from html_common import COMPACT_REPORT_CSS, compress_html_report
from inverse_design_model import (
    build_inverse_design_storyline_html,
    build_inverse_design_technical_html,
    compute_inverse_design_analysis,
)
from ligand_coord_html_section import (
    analyze_precursor_coordination_and_ligands,
    build_ligand_coord_storyline_html,
    build_ligand_coord_technical_html,
)
from main_report_sections import build_main_report_sections
from mechanistic_ml_model import (
    build_mechanistic_storyline_html,
    build_mechanistic_technical_html,
    evaluate_mechanistic_ml_models,
)
from neb_html_section import (
    build_neb_storyline_html,
    build_neb_technical_html,
    build_neb_trajectory_payload_json,
)


def _read_template_asset(filename: str) -> str:
    """Read a modular template asset from ``scripts/benchmark/templates``."""
    asset_path = os.path.join(BENCH_DIR, "templates", filename)
    with open(asset_path, "r", encoding="utf-8") as fh:
        return fh.read()


def build_full_report_html() -> str:
    """Assemble and compress the complete standalone HTML report.

    Physics/Algorithm:
        Executes each domain section builder, injects the modular tab templates
        into ``template.html``, replaces all dynamic data placeholders, and runs
        :func:`compress_html_report` for lossless SVG/table minification.

    Returns:
        Final compressed HTML document string.
    """
    main_placeholders, shared_ctx = build_main_report_sections()
    db_placeholders = build_db_explorer_and_sankey(shared_ctx)

    neb_bench_path = os.path.join("database", "external", "sevenn_decomp_neb_benchmark.json")
    unified_bench_path = os.path.join("database", "benchmark_unified_latest_160.json")
    expansion_cand_path = os.path.join("database", "external", "patent_lit_expansion_candidates.json")

    neb_decomp_storyline_html = build_neb_storyline_html(neb_bench_path, unified_bench_path)
    neb_decomp_technical_html = build_neb_technical_html(neb_bench_path, unified_bench_path)
    neb_trajectory_json = build_neb_trajectory_payload_json(neb_bench_path, unified_bench_path)

    ligand_coord_stats = analyze_precursor_coordination_and_ligands(unified_bench_path)
    ligand_coord_storyline_html = build_ligand_coord_storyline_html(
        unified_bench_path, stats=ligand_coord_stats
    )
    ligand_coord_technical_html = build_ligand_coord_technical_html(
        unified_bench_path, stats=ligand_coord_stats
    )

    inverse_design_stats = compute_inverse_design_analysis(unified_bench_path)
    inverse_design_storyline_html = build_inverse_design_storyline_html(
        unified_bench_path, inv_stats=inverse_design_stats
    )
    inverse_design_technical_html = build_inverse_design_technical_html(
        unified_bench_path, inv_stats=inverse_design_stats
    )

    mech_ml_stats = evaluate_mechanistic_ml_models(unified_bench_path, neb_bench_path)
    mechanistic_ml_storyline_html = build_mechanistic_storyline_html(
        unified_bench_path, neb_bench_path, mech_stats=mech_ml_stats
    )
    mechanistic_ml_technical_html = build_mechanistic_technical_html(
        unified_bench_path, neb_bench_path, mech_stats=mech_ml_stats
    )

    analysis_history_technical_html = build_analysis_history_technical_html(
        unified_bench_path, neb_bench_path, expansion_cand_path
    )
    patent_expansion_storyline_html = build_patent_expansion_storyline_html(
        expansion_cand_path, unified_bench_path
    )
    patent_expansion_technical_html = build_patent_expansion_technical_html(
        expansion_cand_path, unified_bench_path
    )
    br_sascore_theory_html = build_br_sascore_theory_html()
    glossary_html = build_glossary_html()
    intro_html = build_intro_html()

    with open("template.html", "r", encoding="utf-8") as fh:
        html_doc = fh.read()

    shell_replacements: Dict[str, str] = {
        "<!-- REPORT_STYLES_PLACEHOLDER -->": _read_template_asset("report_styles.css") + COMPACT_REPORT_CSS,
        "<!-- TAB_STORYLINE_SHELL_PLACEHOLDER -->": _read_template_asset("tab_storyline.html"),
        "<!-- TAB_TECHNICAL_SHELL_PLACEHOLDER -->": _read_template_asset("tab_technical.html"),
        "<!-- TAB_DATABASE_SHELL_PLACEHOLDER -->": _read_template_asset("tab_database.html"),
        "<!-- REPORT_MAIN_JS_PLACEHOLDER -->": _read_template_asset("report_main.js"),
        "<!-- REPORT_NEB_JS_PLACEHOLDER -->": _read_template_asset("report_neb.js"),
    }
    for token, content in shell_replacements.items():
        html_doc = html_doc.replace(token, content)

    section_replacements: Dict[str, str] = {
        **main_placeholders,
        **db_placeholders,
        "<!-- NEB_DECOMP_STORYLINE_PLACEHOLDER -->": neb_decomp_storyline_html,
        "<!-- NEB_DECOMP_TECHNICAL_PLACEHOLDER -->": neb_decomp_technical_html,
        "<!-- LIGAND_COORD_STORYLINE_PLACEHOLDER -->": ligand_coord_storyline_html,
        "<!-- LIGAND_COORD_TECHNICAL_PLACEHOLDER -->": ligand_coord_technical_html,
        "<!-- INVERSE_DESIGN_STORYLINE_PLACEHOLDER -->": inverse_design_storyline_html,
        "<!-- INVERSE_DESIGN_TECHNICAL_PLACEHOLDER -->": inverse_design_technical_html,
        "<!-- MECHANISTIC_ML_STORYLINE_PLACEHOLDER -->": mechanistic_ml_storyline_html,
        "<!-- MECHANISTIC_ML_TECHNICAL_PLACEHOLDER -->": mechanistic_ml_technical_html,
        "<!-- ANALYSIS_HISTORY_TECHNICAL_PLACEHOLDER -->": analysis_history_technical_html,
        "<!-- PATENT_EXPANSION_STORYLINE_PLACEHOLDER -->": patent_expansion_storyline_html,
        "<!-- PATENT_EXPANSION_TECHNICAL_PLACEHOLDER -->": patent_expansion_technical_html,
        "<!-- BR_SASCORE_THEORY_PLACEHOLDER -->": br_sascore_theory_html,
        "<!-- GLOSSARY_PLACEHOLDER -->": glossary_html,
        "<!-- INTRO_PLACEHOLDER -->": intro_html,
        "__NEB_TRAJECTORY_JSON_PLACEHOLDER__": neb_trajectory_json,
    }
    for token, content in section_replacements.items():
        html_doc = html_doc.replace(token, content)

    return compress_html_report(html_doc)


def main() -> None:
    """Entry point for generating ``index.html`` and ``synthetic_accessibility_report.html``."""
    final_html = build_full_report_html()
    for out_name in ("index.html", "synthetic_accessibility_report.html"):
        with open(out_name, "w", encoding="utf-8") as fh:
            fh.write(final_html)
    print("Generated clean light-theme index.html & synthetic_accessibility_report.html successfully!")


if __name__ == "__main__":
    main()
