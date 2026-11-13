#!/usr/bin/env python3
"""SevenNet-0 Physical Feature Evaluation for Expanded Focus-Element Precursors.

Reads ``database/external/expanded_semiconductor_precursors_master.json``, selects
curated precursors carrying validated ``SMILES`` for target elements
(``Ir, Rh, As, Ga, In, Ru``) that are not already in ``sevenn_physical_features.json``,
and computes the relaxed SevenNet-0 physical descriptors by reusing
``compute_precursor_physical_features`` from ``run_phase2_sevenn_features.py``.
"""

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

FOCUS_ELEMENTS: tuple = ("Ir", "Rh", "As", "Ga", "In", "Ru")
MASTER_PATH = os.path.join(PROJECT_ROOT, "database", "external", "expanded_semiconductor_precursors_master.json")
BASE_SEVENN_PATH = os.path.join(PROJECT_ROOT, "database", "external", "sevenn_physical_features.json")
OUT_PATH = os.path.join(PROJECT_ROOT, "database", "external", "expanded_elements_sevenn_features.json")


def select_target_precursors(
    master_records: List[dict],
    existing_ids: Set[str],
    max_per_element: int = 0,
) -> List[dict]:
    """Select curated focus-element precursors with valid SMILES sorted by ALD frequency.

    Physics/Algorithm:
        Prioritizes non-baseline focus-element records (``Ir, Rh, As, Ga, In, Ru``)
        ordered by descending ``atomiclimits_process_count`` so high-relevance ALD
        precursors are evaluated first under a deterministic budget.
    """
    by_el: Dict[str, List[dict]] = {el: [] for el in FOCUS_ELEMENTS}
    for rec in master_records:
        el = rec.get("center_element")
        if el not in by_el:
            continue
        if rec.get("id") in existing_ids:
            continue
        if not rec.get("smiles"):
            continue
        by_el[el].append(rec)

    selected: List[dict] = []
    for el in FOCUS_ELEMENTS:
        rows = sorted(
            by_el[el],
            key=lambda r: (-int(r.get("atomiclimits_process_count") or 0), r.get("name", "")),
        )
        if max_per_element > 0:
            rows = rows[:max_per_element]
        selected.extend(rows)
    return selected


def run_expanded_calculation(max_per_element: int = 6, num_threads: int = 4, resume: bool = True) -> None:
    """Compute SevenNet-0 LBFGS physical features for curated focus-element precursors.

    Physics/Algorithm:
        Evaluates cohesive energy (``e_coh_ev_per_atom``, units: eV/atom), relaxed
        fragment/stretch dissociation energy (``bde_min_kcal_mol``, units: kcal/mol),
        Cavallo buried volume (``vbur_percent``, units: %), non-overlapping dimerization
        energy (``delta_e_dimer_kcal_mol``, units: kcal/mol), and algebraic complexation
        proxy (``delta_e_complex_kcal_mol``, units: kcal/mol) via ``SevenNet-0``.
    """
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    torch.set_num_threads(num_threads)

    with open(MASTER_PATH, "r", encoding="utf-8") as f:
        master_records = json.load(f)

    existing_ids: Set[str] = set()
    if os.path.exists(BASE_SEVENN_PATH):
        with open(BASE_SEVENN_PATH, "r", encoding="utf-8") as f:
            base_data = json.load(f)
        existing_ids.update(r["id"] for r in base_data.get("features", []) if r.get("id"))

    cached_by_id: Dict[str, dict] = {}
    if resume and os.path.exists(OUT_PATH):
        with open(OUT_PATH, "r", encoding="utf-8") as f:
            prev = json.load(f)
        for rec in prev.get("features", []):
            if rec.get("id") and not rec.get("is_fallback"):
                cached_by_id[rec["id"]] = rec

    targets = select_target_precursors(master_records, existing_ids, max_per_element=max_per_element)
    print(f"Selected {len(targets)} focus-element precursors (max_per_element={max_per_element}, cached={len(cached_by_id)}).")

    calc = SevenNetCalculator(model="7net-0", device="cpu")
    frag_cache: dict = {}
    results: List[dict] = []
    t_start = time.time()

    for idx, item in enumerate(targets, 1):
        pid = item["id"]
        if pid in cached_by_id:
            results.append(cached_by_id[pid])
            print(f"[{idx:02d}/{len(targets):02d}] {pid:11s} | {item['name']:20s} | cached", flush=True)
            continue

        t0 = time.time()
        p_data = compute_precursor_physical_features(item["smiles"], calc, frag_cache=frag_cache)
        dt = time.time() - t0
        is_monomer_volatile = bool((60.0 <= p_data["vbur_percent"] <= 82.0) and (p_data["delta_e_dimer_kcal_mol"] >= -5.0))
        is_thermally_labile = bool(p_data["bde_min"] < 40.0 or p_data["e_coh"] < 3.0)

        rec = {
            "id": pid,
            "name": item["name"],
            "smiles": item["smiles"],
            "formula": item.get("formula") or "-",
            "center_atom": p_data["center_atom"],
            "center_element_expected": item.get("center_element"),
            "atomiclimits_process_count": item.get("atomiclimits_process_count", 0),
            "is_fallback": p_data["is_fallback"],
            "mono_relax_steps": p_data["mono_relax_steps"],
            "mono_fmax_ev_ang": p_data["mono_fmax_ev_ang"],
            "dimer_init_gap_ang": p_data["dimer_init_gap_ang"],
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
            "is_monomer_volatile": is_monomer_volatile,
            "is_thermally_labile": is_thermally_labile,
        }
        results.append(rec)
        print(
            f"[{idx:02d}/{len(targets):02d}] {pid:11s} | {item['name']:20s} | "
            f"%Vbur:{p_data['vbur_percent']:5.1f}% | BDE:{p_data['bde_min']:6.1f} | "
            f"dEdim:{p_data['delta_e_dimer_kcal_mol']:6.1f} | ecoh:{p_data['e_coh']:5.3f} | "
            f"fb:{p_data['is_fallback']} | {dt:.1f}s",
            flush=True,
        )

        _write_payload(results, time.time() - t_start, max_per_element)

    _write_payload(results, time.time() - t_start, max_per_element)
    print(f"Completed {len(results)} focus-element calculations -> {OUT_PATH}")


def _write_payload(results: List[dict], elapsed_s: float, max_per_element: int) -> None:
    by_el: Dict[str, int] = {}
    for r in results:
        el = r.get("center_element_expected") or r.get("center_atom") or "-"
        by_el[el] = by_el.get(el, 0) + 1
    payload = {
        "metadata": {
            "framework": "SevenNet (7net-0 LBFGS Relaxed) + Cavallo SambVca Sphere",
            "focus_elements": list(FOCUS_ELEMENTS),
            "max_per_element": max_per_element,
            "total_molecules_computed": len(results),
            "counts_by_element": by_el,
            "fallback_count": sum(1 for r in results if r.get("is_fallback")),
            "computation_time_seconds": round(elapsed_s, 1),
        },
        "features": results,
    }
    with open(OUT_PATH, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
        f.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--max-per-element", type=int, default=6, help="Max curated entries per element (0 = all)")
    parser.add_argument("--threads", type=int, default=4, help="PyTorch CPU threads")
    parser.add_argument("--no-resume", action="store_true", help="Recompute all entries without cache")
    args = parser.parse_args()
    run_expanded_calculation(
        max_per_element=args.max_per_element,
        num_threads=args.threads,
        resume=not args.no_resume,
    )


if __name__ == "__main__":
    main()
