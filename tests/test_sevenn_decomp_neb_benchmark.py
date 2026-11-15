#!/usr/bin/env python3
"""Unit tests for the SevenNet-0 decomposition NEB benchmark schema and helpers."""

from __future__ import annotations

import json
import os
import sys
import unittest

import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts", "benchmark"))

from run_sevenn_decomp_neb_benchmark import (  # noqa: E402
    BENCHMARK_SPECS,
    REACTION_BUILDERS,
    kabsch_align,
    summarize_results,
)
from neb_html_section import (  # noqa: E402
    build_neb_storyline_html,
    build_neb_technical_html,
    compute_arrhenius_t_threshold_c,
    compute_neb_kinetics_context,
)


class TestSevenNetDecompNebBenchmark(unittest.TestCase):
    """Verify literature benchmark metadata, Kabsch alignment, and aggregation logic."""

    def test_specs_match_builders_and_have_dois(self) -> None:
        self.assertGreaterEqual(len(BENCHMARK_SPECS), 7)
        spec_ids = [s.reaction_id for s in BENCHMARK_SPECS]
        self.assertEqual(len(spec_ids), len(set(spec_ids)))
        for spec in BENCHMARK_SPECS:
            self.assertIn(spec.reaction_id, REACTION_BUILDERS)
            self.assertTrue(spec.lit_doi.startswith("10."))
            self.assertGreater(spec.lit_ea_kcal_mol, 0.0)

    def test_kabsch_align_recovers_rigid_rotation_and_translation(self) -> None:
        rng = np.random.default_rng(7)
        ref = rng.normal(size=(8, 3))
        ref -= np.mean(ref, axis=0)
        theta = 0.65
        rot_true = np.array(
            [
                [np.cos(theta), -np.sin(theta), 0.0],
                [np.sin(theta), np.cos(theta), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        mobile = (ref @ rot_true) + np.array([2.5, -1.2, 0.8])
        aligned = kabsch_align(mobile, ref, anchor_indices=[0, 1, 2, 3, 4])
        self.assertLess(float(np.max(np.abs(aligned - ref))), 1e-10)

    def test_summarize_results_computes_mae_and_ttip_selectivity(self) -> None:
        mock = [
            {
                "spec": {"reaction_id": "TTIP_acetone_elim", "lit_ea_kcal_mol": 49.0},
                "status": {"valid_fs_topology": True, "converged_cineb": True},
                "energies_kcal_mol": {
                    "sevenn_ea_cineb_kcal_mol": 44.0,
                    "signed_error_ea_kcal_mol": -5.0,
                    "abs_error_ea_kcal_mol": 5.0,
                },
                "ts_diagnostics": {"core_vibrations": {"n_imag_modes_gt20cm1": 1}},
            },
            {
                "spec": {"reaction_id": "TTIP_propene_elim", "lit_ea_kcal_mol": 62.0},
                "status": {"valid_fs_topology": True, "converged_cineb": True},
                "energies_kcal_mol": {
                    "sevenn_ea_cineb_kcal_mol": 65.0,
                    "signed_error_ea_kcal_mol": 3.0,
                    "abs_error_ea_kcal_mol": 3.0,
                },
                "ts_diagnostics": {"core_vibrations": {"n_imag_modes_gt20cm1": 1}},
            },
        ]
        summary = summarize_results(mock)
        self.assertEqual(summary["n_reactions"], 2)
        self.assertEqual(summary["mae_ea_kcal_mol"], 4.0)
        self.assertEqual(summary["bias_mean_signed_err_kcal_mol"], -1.0)
        self.assertTrue(summary["ttip_selectivity"]["correct_channel_order"])

    def test_saved_benchmark_artifact_integrity(self) -> None:
        bench_path = os.path.join(
            PROJECT_ROOT, "database", "external", "sevenn_decomp_neb_benchmark.json"
        )
        self.assertTrue(os.path.exists(bench_path), f"Missing {bench_path}")
        with open(bench_path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        summary = data["summary"]
        self.assertEqual(summary["n_reactions"], 8)
        self.assertEqual(summary["n_valid_fs_topology"], 8)
        self.assertLess(summary["mae_ea_kcal_mol"], 6.0)
        self.assertTrue(summary["ttip_selectivity"]["correct_channel_order"])
        for rxn in data["reactions"]:
            self.assertTrue(rxn["status"]["valid_fs_topology"])
            self.assertEqual(len(rxn["energies_kcal_mol"]["band_profile_rel_kcal_mol"]), 7)
            self.assertTrue(rxn["spec"]["lit_doi"].startswith("10."))

    def test_neb_html_section_kinetics_and_render(self) -> None:
        neb_path = os.path.join(
            PROJECT_ROOT, "database", "external", "sevenn_decomp_neb_benchmark.json"
        )
        unified_path = os.path.join(
            PROJECT_ROOT, "database", "benchmark_unified_latest_160.json"
        )
        ctx = compute_neb_kinetics_context(neb_path, unified_path)
        self.assertEqual(ctx["corr_gas"]["n"], 5)
        self.assertGreater(ctx["corr_gas"]["r_lit"], 0.95)
        self.assertGreater(ctx["corr_gas"]["r_7n"], 0.80)
        self.assertGreater(ctx["corr_gas"]["r_7n"], ctx["corr_gas"]["r_bde"])
        self.assertAlmostEqual(
            compute_arrhenius_t_threshold_c(41.55, pre_exp_s=1.0e13),
            332.3,
            places=1,
        )

        story_html = build_neb_storyline_html(neb_path, unified_path)
        tech_html = build_neb_technical_html(neb_path, unified_path)
        for chunk in (story_html, tech_html):
            self.assertIn("CpZr", chunk)
            self.assertNotIn("\x07", chunk)
            self.assertNotIn("\x0c", chunk)
            self.assertNotIn("\\x08", chunk)
        self.assertIn("T_{\\text{decomp}}", tech_html)
        self.assertIn('G{"자동 품질 게이트', tech_html)


if __name__ == "__main__":
    unittest.main()
