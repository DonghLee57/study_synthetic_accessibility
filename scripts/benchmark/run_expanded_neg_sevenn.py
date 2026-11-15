#!/usr/bin/env python3
"""
SevenNet-0 Physical Feature Evaluation for Critical Negative Controls
---------------------------------------------------------------------
Computes relaxed SevenNet-0 (7net-0) LBFGS geometry, cohesive energy, isolated-fragment
BDE_min, buried volume, and non-overlapping dimerization energy for key negative controls
lacking prior MLIP calculations (NEG-019, NEG-038, NEG-039, NEG-040) by reusing
``compute_precursor_physical_features`` from ``run_phase2_sevenn_features.py``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Set

import torch
from sevenn.calculator import SevenNetCalculator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from run_phase2_sevenn_features import compute_precursor_physical_features  # noqa: E402

BASE_SEVENN_PATH = os.path.join(PROJECT_ROOT, "database", "external", "sevenn_physical_features.json")
NEG_DB_PATH = os.path.join(PROJECT_ROOT, "database", "expanded_negative_controls.json")
OUT_PATH = os.path.join(PROJECT_ROOT, "database", "external", "expanded_neg_sevenn_features.json")


def _write_payload(results: List[dict], elapsed_s: float) -> None:
    payload = {
        "metadata": {
            "framework": "SevenNet (7net-0 LBFGS Relaxed) + Cavallo SambVca Sphere",
            "selection_mode": "all_expanded_negative_controls",
            "target_ids": [r["id"] for r in results],
            "total_molecules_computed": len(results),
            "fallback_count": sum(1 for r in results if r.get("is_fallback")),
            "computation_time_seconds": round(elapsed_s, 1),
        },
        "features": results,
    }
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
        f.write("\n")


def main(num_threads: int = 4, resume: bool = True) -> None:
    """Compute SevenNet LBFGS physical features for all expanded negative controls."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    torch.set_num_threads(num_threads)

    with open(NEG_DB_PATH, "r", encoding="utf-8") as f:
        neg_records = json.load(f)["records"]

    base_ids: Set[str] = set()
    if os.path.exists(BASE_SEVENN_PATH):
        with open(BASE_SEVENN_PATH, "r", encoding="utf-8") as f:
            base_ids.update(r["id"] for r in json.load(f).get("features", []) if r.get("id"))

    cached_by_id: Dict[str, dict] = {}
    if resume and os.path.exists(OUT_PATH):
        with open(OUT_PATH, "r", encoding="utf-8") as f:
            for rec in json.load(f).get("features", []):
                if rec.get("id") and not rec.get("is_fallback"):
                    cached_by_id[rec["id"]] = rec

    selected = [r for r in neg_records if r["id"] not in base_ids]
    print(f"Selected {len(selected)} expanded negative controls (cached={len(cached_by_id)}, to_compute={len(selected) - len(cached_by_id)}).", flush=True)

    calc = SevenNetCalculator(model="7net-0", device="cpu")
    frag_cache: Dict[str, float] = {}
    results: List[dict] = []
    t_start = time.time()

    for idx, item in enumerate(selected, 1):
        nid = item["id"]
        if nid in cached_by_id:
            results.append(cached_by_id[nid])
            print(f"[{idx:02d}/{len(selected):02d}] {nid} | {item['name'][:24]:24s} | cached", flush=True)
            continue

        t0 = time.time()
        p_data = compute_precursor_physical_features(item["smiles"], calc, frag_cache=frag_cache)
        dt = time.time() - t0
        rec = {
            "id": nid,
            "name": item["name"],
            "smiles": item["smiles"],
            "failure_mechanism": item.get("failure_mechanism", ""),
            "center_atom": p_data["center_atom"],
            "is_fallback": p_data["is_fallback"],
            "mono_relax_steps": p_data["mono_relax_steps"],
            "mono_fmax_ev_ang": p_data["mono_fmax_ev_ang"],  # units: eV/A
            "dimer_init_gap_ang": p_data["dimer_init_gap_ang"],  # units: A
            "dimer_relax_steps": p_data["dimer_relax_steps"],
            "bde_mode": p_data["bde_mode"],
            "e_pot_ev": p_data["e_eq"],  # units: eV
            "e_coh_ev_per_atom": p_data["e_coh"],  # units: eV/atom
            "bde_min_kcal_mol": p_data["bde_min"],  # units: kcal/mol
            "vbur_percent": p_data["vbur_percent"],  # units: %
            "delta_e_dimer_kcal_mol": p_data["delta_e_dimer_kcal_mol"],  # units: kcal/mol
            "delta_e_complex_kcal_mol": p_data["delta_e_complex_kcal_mol"],  # units: kcal/mol
            "radius_of_gyration_ang": p_data["rg"],  # units: A
            "mol_volume_ang3": p_data["vol_mol"],  # units: A^3
        }
        results.append(rec)
        print(
            f"[{idx:02d}/{len(selected):02d}] {nid} | {item['name'][:24]:24s} | "
            f"%Vbur={p_data['vbur_percent']:.1f}% | BDE={p_data['bde_min']:.2f} kcal/mol | "
            f"ecoh={p_data['e_coh']:.3f} eV | dEdim={p_data['delta_e_dimer_kcal_mol']:+.2f} kcal/mol | {dt:.1f}s",
            flush=True,
        )
        _write_payload(results, time.time() - t_start)

    _write_payload(results, time.time() - t_start)
    print(f"Saved {len(results)} -> {OUT_PATH}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--threads", type=int, default=4, help="PyTorch CPU threads")
    parser.add_argument("--no-resume", action="store_true", help="Recompute all entries without cache")
    cli_args = parser.parse_args()
    main(num_threads=cli_args.threads, resume=not cli_args.no_resume)
