#!/usr/bin/env python3
"""
Unified Semiconductor Precursor SA & Coordination-Physics Analysis (N = 160)
----------------------------------------------------------------------------
Orchestrates the end-to-end 160-molecule evaluation (27 Audited Precursors +
75 Expanded Precursors + 18 AS-ALD Inhibitors + 40 Negative Controls across
6 failure mechanisms) by composing shared benchmark modules:
  - ``element_physics``: Multi-state periodic coordination & thermodynamic engine
  - ``tmc_sa_scorer``: 2D/SMARTS coordination descriptors, gates, and TMC-SA scorer
  - ``ml_evaluation``: Stratified CV-5, LOFMO extrapolation, Descriptor 4 audit, PCA/t-SNE

Outputs:
  - database/benchmark_unified_latest_160.json
"""

from __future__ import annotations

import json
import os
import sys
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import RDConfig, rdMolDescriptors

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(RDConfig.RDContribDir, "SA_Score"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts", "benchmark"))

import sascorer  # type: ignore
from syba.syba import SybaClassifier  # type: ignore

from element_physics import (
    DELTA_R_SHELL_ANG,
    compute_formation_energy_metrics,
    compute_intrinsic_topology,
    universal_vbur_lower_bound,
)
from ml_evaluation import (
    compute_latent_projections,
    evaluate_cv5_and_lofmo,
    evaluate_descriptor4_audit,
    load_json,
    rate_pct,
    save_json,
)
from tmc_sa_scorer import (
    compute_chemical_coord_descriptors,
    compute_coord_risk_count,
    compute_design12_geometry_metrics,
    compute_ligand_sa,
    compute_tmc_sa_score,
    evaluate_coord_chem_gate,
)

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

try:
    from BRSAScore import SAScorer as BRSAScorer  # type: ignore
except Exception:
    BRSAScorer = None


def _build_unified_entries(
    p2_rows: List[Dict[str, Any]],
    exp_db: Dict[str, Any],
    neg_db: Dict[str, Any],
    exp_7n_by_id: Dict[str, Dict[str, Any]],
    neg_7n_by_id: Dict[str, Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Assemble standardized 160-molecule metadata table and SevenNet energy lookup."""
    phys_by_id: Dict[str, Dict[str, Any]] = {r["id"]: r for r in p2_rows}
    phys_by_id.update(exp_7n_by_id)
    phys_by_id.update(neg_7n_by_id)

    entries: List[Dict[str, Any]] = []
    seen_ids = set()

    for r in p2_rows:
        if r["target_type"] in ("Precursor", "Inhibitor"):
            grp = "Audited Precursor (27)" if r["target_type"] == "Precursor" else "AS-ALD Inhibitor (18)"
            entries.append({
                "id": r["id"],
                "name": r["name"],
                "smiles": r["smiles"],
                "group": grp,
                "target_type": r["target_type"],
                "label": 1,
                "failure_mechanism": "pos",
                "origin": "audited_55",
            })
            seen_ids.add(r["id"])

    for r in exp_db["precursors_102"]["records"]:
        if r["id"] in seen_ids or r.get("is_master_precursor"):
            continue
        entries.append({
            "id": r["id"],
            "name": r["name"],
            "smiles": r["smiles"],
            "group": "Expanded Precursor (75)",
            "target_type": "Precursor",
            "label": 1,
            "failure_mechanism": "pos",
            "origin": "expanded_75",
        })
        seen_ids.add(r["id"])

    for r in neg_db["records"]:
        entries.append({
            "id": r["id"],
            "name": r["name"],
            "smiles": r["smiles"],
            "group": "Negative Control (40)",
            "target_type": "Negative Control",
            "label": 0,
            "failure_mechanism": r["failure_mechanism"],
            "origin": r.get("origin", "expanded_neg_30"),
        })

    assert len(entries) == 160, f"Expected 160 molecules, got {len(entries)}"
    return entries, phys_by_id


def main() -> None:
    p2_rows = load_json(os.path.join(PROJECT_ROOT, "database", "benchmark_phase2_tuned_comparison.json"))["molecules_evaluated"]
    exp_db = load_json(os.path.join(PROJECT_ROOT, "database", "benchmark_expanded_models_db.json"))
    neg_db = load_json(os.path.join(PROJECT_ROOT, "database", "expanded_negative_controls.json"))
    exp_7n = load_json(os.path.join(PROJECT_ROOT, "database", "external", "expanded_elements_sevenn_features.json"))["features"]

    neg_7n_path = os.path.join(PROJECT_ROOT, "database", "external", "expanded_neg_sevenn_features.json")
    neg_7n_list = load_json(neg_7n_path)["features"] if os.path.exists(neg_7n_path) else []

    entries, phys_by_id = _build_unified_entries(
        p2_rows, exp_db, neg_db,
        {r["id"]: r for r in exp_7n},
        {r["id"]: r for r in neg_7n_list},
    )

    syba_clf = SybaClassifier()
    syba_clf.fitDefaultScore()
    br_scorer = BRSAScorer() if BRSAScorer is not None else None

    records: List[Dict[str, Any]] = []
    mech_stats: Dict[str, Dict[str, int]] = defaultdict(lambda: {"total": 0, "sascore_rejected": 0, "syba_rejected": 0, "unified_rejected": 0})

    for entry in entries:
        mid, smi = entry["id"], entry["smiles"]
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES: {smi}"
        formula = rdMolDescriptors.CalcMolFormula(mol)

        desc = compute_chemical_coord_descriptors(smi)
        topo = compute_intrinsic_topology(smi)
        geom = compute_design12_geometry_metrics(smi)
        center_sym = str(geom["center_sym"])
        r_cov_ang = float(geom["r_cov_ang"])
        r_probe_ang = float(geom["r_probe_ang"])
        n_val = float(geom["n_val"])
        chi_pauling = float(geom["chi_pauling"])
        is_coord_center = bool(topo["is_coord_center"])

        sa_val = float(sascorer.calculateScore(mol))
        syba_val = float(syba_clf.predict(mol=mol))
        br_val = float(br_scorer.calculateScore(smi)[0]) if br_scorer is not None else sa_val

        phys_rec = phys_by_id.get(mid)
        if phys_rec is not None and "ligand_sa" in phys_rec:
            lig_sa = float(phys_rec["ligand_sa"])
        else:
            lig_sa = compute_ligand_sa(smi, sa_val, br_val, br_scorer=br_scorer)

        gate = evaluate_coord_chem_gate(smi, lig_sa, desc, topo, include_design2=True)
        v_min = universal_vbur_lower_bound(int(topo["pi_hapticity_index"]))

        if phys_rec is not None:
            vbur_r35 = float(phys_rec["vbur_percent"])
            radial_shift = (float(geom["vbur_rprobe"]) - float(geom["vbur_r35"])) if geom["embed_ok"] else 0.0
            vbur_val = round(max(0.0, min(100.0, vbur_r35 + radial_shift)), 2)
            vbur_d2_val = round(max(0.0, min(100.0, vbur_r35 + (float(geom["vbur_d2_percent"]) - float(geom["vbur_r35"])))), 2) if geom["embed_ok"] else vbur_val
            bde_val: Optional[float] = float(phys_rec["bde_min_kcal_mol"])
            ddim_val: Optional[float] = float(phys_rec["delta_e_dimer_kcal_mol"])
            ecoh_val: Optional[float] = float(phys_rec["e_coh_ev_per_atom"])
            epot_val: Optional[float] = float(phys_rec["e_pot_ev"])
            dcomp_legacy_val: Optional[float] = float(phys_rec["delta_e_complex_kcal_mol"])
            vbur_source = "SevenNet-0-LBFGS+Design2State"
            syms_mol = geom.get("syms") or [a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms()]
            form_metrics = compute_formation_energy_metrics(epot_val, syms_mol, int(geom["center_degree"]))
            de_form_atom: Optional[float] = float(form_metrics["delta_e_form_ev_per_atom"])
            de_form_bond: Optional[float] = float(form_metrics["delta_e_form_per_bond_kcal_mol"])
            dg_dimer_373: Optional[float] = round(ddim_val + float(geom["td_s_dimer_373"]), 2)
        else:
            vbur_r35 = float(geom["vbur_r35"])
            vbur_val = float(geom["vbur_rprobe"])
            vbur_d2_val = float(geom["vbur_d2_percent"])
            bde_val = ddim_val = ecoh_val = epot_val = dcomp_legacy_val = de_form_atom = de_form_bond = dg_dimer_373 = None
            vbur_source = "RDKit-ETKDGv3-Design2State"

        if phys_rec is not None and "ligand_sa" in phys_rec:
            base_sa = float(phys_rec["ligand_sa"])
            if entry["target_type"] == "Inhibitor" and not is_coord_center:
                base_sa = float(phys_rec.get("br_sascore", base_sa))
            tmc_sa = compute_tmc_sa_score(base_sa, vbur_val, v_min, is_coord_center, bde_val=bde_val, ecoh_val=ecoh_val, ddim_val=ddim_val)
        else:
            center_h = mol.GetAtomWithIdx(int(topo["center_idx"])).GetTotalNumHs()
            tmc_sa = compute_tmc_sa_score(
                lig_sa, vbur_val, v_min, is_coord_center, bde_val=bde_val,
                is_d_block=bool(topo["is_d_block"]), ve_deficit=float(geom.get("ve_deficit", 0.0)),
                chi_pauling=chi_pauling, center_h=center_h,
            )

        tmc_pass = (tmc_sa <= 4.0) and gate["passed"]
        coord_risk_count = compute_coord_risk_count(desc, topo, is_coord_center)
        composite_risk = round(min(9.8, tmc_sa + 4.0 * coord_risk_count), 3)
        risk_index = round(
            max(0.0, (sa_val - 4.0) * 0.5) + 4.0 * coord_risk_count + (1.2 if (vbur_val > 83.5 and is_coord_center) else 0.0),
            3,
        )

        mech = entry["failure_mechanism"]
        if entry["label"] == 0:
            mech_stats[mech]["total"] += 1
            if sa_val > 4.0:
                mech_stats[mech]["sascore_rejected"] += 1
            if syba_val < 0.0:
                mech_stats[mech]["syba_rejected"] += 1
            if not tmc_pass:
                mech_stats[mech]["unified_rejected"] += 1

        records.append({
            **entry,
            "formula": formula,
            "center_atom": center_sym,
            "r_cov_ang": r_cov_ang,
            "r_probe_ang": r_probe_ang,
            "r_eff_ang": float(geom["r_eff_ang"]),
            "r_probe_d2_ang": float(geom["r_probe_d2_ang"]),
            "n_val": n_val,
            "chi_pauling": chi_pauling,
            "chi_eff": float(geom["chi_eff"]),
            "delta_chi_mean": float(geom["delta_chi_mean"]),
            "ionicity_mean": float(geom["ionicity_mean"]),
            "cn_eff": float(geom["cn_eff"]),
            "n_ve": float(geom["n_ve"]),
            "ve_deficit": float(geom["ve_deficit"]),
            "coord_propensity": float(geom["coord_propensity"]),
            "sascore": round(sa_val, 3),
            "syba": round(syba_val, 3),
            "br_sascore": round(br_val, 3),
            "ligand_sa": round(lig_sa, 3),
            "vbur_r35_percent": round(vbur_r35, 2),
            "vbur_percent": round(vbur_val, 2),
            "vbur_d2_percent": round(vbur_d2_val, 2),
            "void_cone_deg": float(geom["void_cone_deg"]),
            "v_min_adaptive": round(v_min, 1),
            "vbur_margin": round(vbur_val - v_min, 2),
            "vbur_source": vbur_source,
            "bde_min_kcal_mol": round(bde_val, 2) if bde_val is not None else None,
            "delta_e_dimer_kcal_mol": round(ddim_val, 2) if ddim_val is not None else None,
            "delta_g_dimer_373_kcal_mol": dg_dimer_373,
            "e_coh_ev_per_atom": round(ecoh_val, 3) if ecoh_val is not None else None,
            "delta_e_form_ev_per_atom": de_form_atom,
            "delta_e_form_per_bond_kcal_mol": de_form_bond,
            "delta_e_complex_legacy_kcal_mol": round(dcomp_legacy_val, 2) if dcomp_legacy_val is not None else None,
            "pi_hapticity_index": int(topo["pi_hapticity_index"]),
            "chelate_ring_count": int(topo["chelate_ring_count"]),
            "is_homoleptic_tm_sp3_refined": bool(topo["is_homoleptic_tm_sp3_refined"]),
            "catenated_main_group_hydride": bool(topo.get("catenated_main_group_hydride")),
            "pi_conjugated_main_group_hydride": bool(topo.get("pi_conjugated_main_group_hydride")),
            "tmc_sa_v2": round(tmc_sa, 3),
            "unified_pass": bool(tmc_pass),
            "gate_reasons": gate["reasons"],
            "coord_risk_count": int(coord_risk_count),
            "composite_risk": composite_risk,
            "risk_index": risk_index,
            "descriptors": desc,
        })

    # Summary & statistical evaluations
    precursors_all = [r for r in records if r["target_type"] == "Precursor"]
    prec_aud = [r for r in records if r["group"] == "Audited Precursor (27)"]
    prec_exp = [r for r in records if r["group"] == "Expanded Precursor (75)"]
    inhib_all = [r for r in records if r["target_type"] == "Inhibitor"]
    pos_all = [r for r in records if r["label"] == 1]
    neg_all = [r for r in records if r["label"] == 0]

    prec_tp = sum(1 for r in precursors_all if r["unified_pass"])
    prec_aud_tp = sum(1 for r in prec_aud if r["unified_pass"])
    prec_exp_tp = sum(1 for r in prec_exp if r["unified_pass"])
    inhib_tp = sum(1 for r in inhib_all if r["unified_pass"])
    pos_tp = sum(1 for r in pos_all if r["unified_pass"])
    neg_tn = sum(1 for r in neg_all if not r["unified_pass"])
    neg_fp = [r["id"] for r in neg_all if r["unified_pass"]]

    neg_sas_tn = sum(1 for r in neg_all if r["sascore"] > 4.0)
    neg_syba_tn = sum(1 for r in neg_all if r["syba"] < 0.0)
    prec_sas_tp = sum(1 for r in precursors_all if r["sascore"] <= 3.5)

    desc_keys = sorted(records[0]["descriptors"].keys())
    scalar_cols = [
        "sascore", "syba", "vbur_percent", "vbur_d2_percent", "void_cone_deg",
        "r_cov_ang", "r_probe_ang", "r_eff_ang", "r_probe_d2_ang",
        "n_val", "chi_pauling", "chi_eff", "delta_chi_mean", "ionicity_mean",
        "cn_eff", "n_ve", "ve_deficit", "coord_propensity",
        "pi_hapticity_index", "chelate_ring_count",
    ]
    feature_names = scalar_cols + desc_keys
    X = np.array([[float(r[c]) for c in scalar_cols] + [float(r["descriptors"][k]) for k in desc_keys] for r in records], dtype=float)
    y = np.array([r["label"] for r in records], dtype=int)

    X_pca, X_tsne, latent_stats = compute_latent_projections(X, y)
    for i, r in enumerate(records):
        r["pc1"] = round(float(X_pca[i, 0]), 4)
        r["pc2"] = round(float(X_pca[i, 1]), 4)
        r["tsne1"] = round(float(X_tsne[i, 0]), 4)
        r["tsne2"] = round(float(X_tsne[i, 1]), 4)

    ml_cv5, lofmo = evaluate_cv5_and_lofmo(X, y, [r["failure_mechanism"] for r in records])
    descriptor4_audit_85 = evaluate_descriptor4_audit(records)

    mech_summary = {
        mech: {
            **st,
            "sascore_tnr_pct": rate_pct(st["sascore_rejected"], st["total"]),
            "syba_tnr_pct": rate_pct(st["syba_rejected"], st["total"]),
            "unified_tnr_pct": rate_pct(st["unified_rejected"], st["total"]),
        }
        for mech, st in mech_stats.items()
    }
    unique_mechs = sorted({r["failure_mechanism"] for r in records if r["failure_mechanism"] != "pos"})
    lofmo["Unified_TMC_SA_v2_Gate"] = {
        "mean_unseen_tnr": round(float(np.mean([mech_summary[m]["unified_tnr_pct"] / 100.0 for m in unique_mechs])), 4),
        "per_unseen_mechanism_tnr": {m: round(mech_summary[m]["unified_tnr_pct"] / 100.0, 4) for m in unique_mechs},
    }

    refine_res = load_json(os.path.join(PROJECT_ROOT, "database", "benchmark_refinement_and_rxn_check.json"))
    out = {
        "dataset_summary": {
            "total_molecules": len(records),
            "positive_total": len(pos_all),
            "precursors_total": len(precursors_all),
            "audited_precursors": len(prec_aud),
            "expanded_precursors": len(prec_exp),
            "inhibitors_total": len(inhib_all),
            "negative_controls_total": len(neg_all),
            "sevenn_relaxed_count": sum(1 for r in records if r["vbur_source"].startswith("SevenNet-0")),
            "etkdg_3d_count": sum(1 for r in records if r["vbur_source"].startswith("RDKit-ETKDGv3")),
            "feature_count": len(feature_names),
            "feature_names": feature_names,
        },
        "design1_model_spec": {
            "mode": "Element-Property-Driven Design 1 (No Element-Symbol Branching)",
            "probe_radius_formula": f"R_probe(Z_M) = R_cov(Z_M) + {DELTA_R_SHELL_ANG:.2f} A (Cordero 2008)",
            "v_min_formula": "V_min^univ(pi_true) = 44.0 - 5.0 * min(2, pi_true)",
            "descriptor4_replacement": "delta_e_form_per_bond_kcal_mol = (E_pot - sum n_Z*mu_Z^std)*23.0605 / deg(M)",
            "descriptor3_entropy_complement": "delta_g_dimer_373_kcal_mol = delta_e_dimer - T*DeltaS_dimer^trans+rot(373.15 K)",
            "center_atom_rule": "Pauling chi_P <= 2.40 outside organic nonmetals + d-block Z range (Cr auto-detected)",
        },
        "design2_model_spec": {
            "mode": "State-Dependent Multi-Property Coordination & Polarity-Weighted %V_bur^D2 (Design 2)",
            "state_effective_radius": "R_eff(M|env) = clip(R_cov*(1 - 0.06*max(0, mean_DeltaChi) + 0.025*(CN_eff - 4)), 0.75*R_cov, 1.22*R_cov)",
            "state_effective_chi": "chi_eff(M|env) = chi_0 + 0.12*sum(sgn(DeltaChi_l)*(1 - exp(-0.25*DeltaChi_l^2)))",
            "tm_coordination_vector": "[delta_chi_mean, chi_eff, r_eff_ang, cn_eff, n_ve, ve_deficit, coord_propensity]",
            "weighted_vbur_d2": "w_MH = 0.0 (terminal M-H) and w_j = clip(1.0 + 0.20*(chi_P(j) - 2.55), 0.75, 1.25)",
            "open_void_cone": "200-pt Fibonacci S^2 largest open ray half-angle theta_max_void (deg)",
            "hydride_gates": [
                "catenated_main_group_hydride: non-d-block hydride with M-M bond & B/C/Al/Si/Ga/Ge/Sn hydride",
                "pi_conjugated_main_group_hydride: non-d-block hydride on Cp/pi-conjugated ligand (resolves NEG-040 Cp2GaH)",
                "sevenn_labile_hydride_pen: non-d-block electron-deficient (ve_deficit > 0, chi_P <= 2.05) hydride with 0 <= SevenNet-0 BDE_min < 30 kcal/mol",
            ],
        },
        "unified_v2_summary": {
            "precursors_all_tp": f"{prec_tp}/{len(precursors_all)}",
            "precursors_all_tpr_pct": rate_pct(prec_tp, len(precursors_all)),
            "precursors_all_fnr_pct": round(100.0 - rate_pct(prec_tp, len(precursors_all)), 1),
            "audited_precursors_tp": f"{prec_aud_tp}/{len(prec_aud)}",
            "audited_precursors_tpr_pct": rate_pct(prec_aud_tp, len(prec_aud)),
            "expanded_precursors_tp": f"{prec_exp_tp}/{len(prec_exp)}",
            "expanded_precursors_tpr_pct": rate_pct(prec_exp_tp, len(prec_exp)),
            "inhibitors_tp": f"{inhib_tp}/{len(inhib_all)}",
            "inhibitors_tpr_pct": rate_pct(inhib_tp, len(inhib_all)),
            "positive_all_tp": f"{pos_tp}/{len(pos_all)}",
            "positive_all_tpr_pct": rate_pct(pos_tp, len(pos_all)),
            "negative_all_tn": f"{neg_tn}/{len(neg_all)}",
            "negative_all_tnr_pct": rate_pct(neg_tn, len(neg_all)),
            "negative_all_fpr_pct": round(100.0 - rate_pct(neg_tn, len(neg_all)), 1),
            "negative_fp_ids": neg_fp,
            "baseline_sascore_precursors_tp": f"{prec_sas_tp}/{len(precursors_all)}",
            "baseline_sascore_precursors_tpr_pct": rate_pct(prec_sas_tp, len(precursors_all)),
            "baseline_sascore_neg_tn": f"{neg_sas_tn}/{len(neg_all)}",
            "baseline_sascore_neg_tnr_pct": rate_pct(neg_sas_tn, len(neg_all)),
            "baseline_syba_neg_tn": f"{neg_syba_tn}/{len(neg_all)}",
            "baseline_syba_neg_tnr_pct": rate_pct(neg_syba_tn, len(neg_all)),
        },
        "descriptor4_audit_85": descriptor4_audit_85,
        "latent_space_160": latent_stats,
        "by_failure_mechanism": mech_summary,
        "ml_cv5_160": ml_cv5,
        "lofmo_extrapolation_160": lofmo,
        "item2_reaction_energy_summary": refine_res["item_2_reaction_energy_check"]["reactions"],
        "records": records,
    }

    out_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
    save_json(out_path, out)
    print(f"Saved unified 160-molecule analysis -> {out_path}")
    print(json.dumps({k: v for k, v in out.items() if k not in ("records", "item2_reaction_energy_summary")}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
