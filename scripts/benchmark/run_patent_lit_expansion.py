"""
Generate and screen a curated, non-overlapping semiconductor precursor and ASD
inhibitor expansion candidate set sourced from KR/JP/CN/US/TW patents and
metal-centered peer-reviewed literature.
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from typing import Any, Dict, List

from rdkit import Chem

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(CURRENT_DIR))
EXTERNAL_DIR = os.path.join(PROJECT_ROOT, "external")
for p in (CURRENT_DIR, os.path.join(EXTERNAL_DIR, "br_sascore")):
    if p not in sys.path:
        sys.path.insert(0, p)

from element_physics import (  # noqa: E402
    MAIN_NONMETALS,
    compute_intrinsic_topology,
    universal_vbur_lower_bound,
)
from tmc_sa_scorer import (  # noqa: E402
    compute_chemical_coord_descriptors,
    compute_design12_geometry_metrics,
    compute_ligand_sa,
    compute_tmc_sa_score,
    evaluate_coord_chem_gate,
    sascorer,
)

try:
    from BRSAScore import SAScorer as BRSAScorer  # type: ignore
except Exception:
    BRSAScorer = None


# (name, smiles, group_label, family_class, target_application, jurisdictions, refs, mechanism_note)
EXPANSION_SPECS: List[tuple] = [
    # --- Group 6 (Mo, W) ---
    (
        "MoO2Cl2",
        "O=[Mo](=O)(Cl)Cl",
        "Group 6 (Mo/W)",
        "Oxohalide",
        "F-free Mo word-line / DRAM contact ALD",
        "US, KR, EP, JP",
        "US20250250673A1; US20240239684A1; KR20230088807A; EP4019470A1; DOI:10.1039/D6TC00454G",
        "Volatile Mo(VI) oxochloride replacing WF6 to suppress F-induced SiO2 damage; solid delivery via Subleem/Chemguard.",
    ),
    (
        "MoCl5",
        "Cl[Mo](Cl)(Cl)(Cl)Cl",
        "Group 6 (Mo/W)",
        "Pentahalide",
        "Mo / MoCx / TMD ALD & self-etching AS-ALD",
        "US, CN, KR",
        "Merck/EMD; Entegris; CN11... MoCl5+DHP; DOI:10.1021/acs.chemmater.6b00248",
        "Lowers O-impurity vs MoO2Cl2 when paired with H2 or cyclohexadiene/dihydropyrazine reductants; shows selective oxide self-etching.",
    ),
    (
        "Mo(CO)6",
        "[C-]#[O+]->[Mo](<-[C-]#[O+])(<-[C-]#[O+])(<-[C-]#[O+])(<-[C-]#[O+])<-[C-]#[O+]",
        "Group 6 (Mo/W)",
        "Zero-valent Carbonyl",
        "Low-T Mo / MoN / MoS2 ALD & PEALD",
        "US, KR, JP",
        "PatSnap Mo ALD review; DOI:10.1021/cr900056b",
        "18-electron octahedral Mo(0) complex enabling halide-free ligand dissociation below 200 C.",
    ),
    (
        "MoO2(thd)2",
        "CC(C)(C)C(=O)C=C(C(C)(C)C)O[Mo](=O)(=O)OC(=CC(=O)C(C)(C)C)C(C)(C)C",
        "Group 6 (Mo/W)",
        "Oxo-beta-diketonate",
        "MoOx / Mo-doped oxide ALD",
        "US, EP",
        "DOI:10.1021/acs.chemmater.6b00248",
        "Halide-free hexavalent Mo(VI) chelate with bulky tBu groups suppressing dimerization.",
    ),
    (
        "WCl5",
        "Cl[W](Cl)(Cl)(Cl)Cl",
        "Group 6 (Mo/W)",
        "Pentahalide",
        "F-free W contact & gate-fill ALD",
        "US, KR, TW",
        "US12077858B2; US12002679B2; DOI:10.1016/j.apsusc.2021.150373",
        "Fluorine-free tungsten halide for high-aspect-ratio gap fill without WF6 fluorine diffusion.",
    ),
    (
        "W(CO)6",
        "[C-]#[O+]->[W](<-[C-]#[O+])(<-[C-]#[O+])(<-[C-]#[O+])(<-[C-]#[O+])<-[C-]#[O+]",
        "Group 6 (Mo/W)",
        "Zero-valent Carbonyl",
        "Low-T W / WN / WCx ALD",
        "US, JP",
        "DOI:10.1021/cr900056b",
        "18-electron octahedral W(0) carbonyl benchmark for halide-free low-T seed deposition.",
    ),
    # --- Group 5 (Nb, Ta) ---
    (
        "Nb(OEt)5",
        "CCO[Nb](OCC)(OCC)(OCC)OCC",
        "Group 5 (Nb/Ta)",
        "Alkoxide",
        "Nb2O5 / LiNbO3 / high-k capacitor ALD",
        "US, JP, KR",
        "DOI:10.1021/la902289h",
        "Liquid Nb(V) ethoxide with clean beta-H/hydrolysis reactivity for conformal Nb2O5.",
    ),
    (
        "TBTDEN",
        "CCN(CC)[Nb](=NC(C)(C)C)(N(CC)CC)N(CC)CC",
        "Group 5 (Nb/Ta)",
        "Imido-Amide",
        "Superconducting NbCN & NbN barrier ALD",
        "US, JP",
        "DOI:10.1021/acs.chemmater.5c01456",
        "Nb analogue of TBTDET combining a strong Nb=NtBu anchor with labile NEt2 ligands.",
    ),
    (
        "NbCl5",
        "Cl[Nb](Cl)(Cl)(Cl)Cl",
        "Group 5 (Nb/Ta)",
        "Pentahalide",
        "NbN / Nb2O5 / ferroelectric doping ALD",
        "US, CN",
        "DOI:10.1021/cr900056b; DOI:10.1021/acs.chemmater.8b03082",
        "Carbon-free inorganic Nb halide delivering high thermal stability in >350 C ALD.",
    ),
    (
        "PDMAT",
        "CN(C)[Ta](N(C)C)(N(C)C)(N(C)C)N(C)C",
        "Group 5 (Nb/Ta)",
        "Alkylamide",
        "Ta2O5 / TaN barrier ALD",
        "US, KR, JP",
        "DOI:10.1116/1.2345205",
        "Pentakis(dimethylamido)tantalum(V); small Ta-N radius yields %Vbur^D2 = 93.5% (high steric crowding boundary).",
    ),
    (
        "TaCl5",
        "Cl[Ta](Cl)(Cl)(Cl)Cl",
        "Group 5 (Nb/Ta)",
        "Pentahalide",
        "High-T TaN / Ta2O5 ALD",
        "US, KR",
        "DOI:10.1021/cr900056b",
        "Classic carbon-free Ta(V) halide standard.",
    ),
    # --- Group 4 (Hf, Ti) ---
    (
        "CpHf(NMe2)3 (HyALD)",
        "CN(C)[Hf](N(C)C)(N(C)C)C1=CC=CC1",
        "Group 4 (Hf/Ti)",
        "Cp-Alkylamide Heteroleptic",
        "HZO ferroelectric & DRAM high-k HfO2 ALD",
        "JP, KR, US, CN",
        "JP6948159B2; JP6954776B2; JP5128289B2; US8003814B2; DOI:10.1116/1.4972210",
        "ADEKA/SK Trichem benchmark; fills the missing Hf analogue of CpZr(NMe2)3 (ZAC) and CpTi(NMe2)3.",
    ),
    (
        "MeCpHf(NMe2)3",
        "CC1=CC=CC1[Hf](N(C)C)(N(C)C)N(C)C",
        "Group 4 (Hf/Ti)",
        "Alkyl-Cp Heteroleptic",
        "Low-melting liquid DRAM HfO2 ALD",
        "CN, KR, JP, US",
        "CN113423862A; CN101310037B; SK Trichem / ADEKA JP7847142B2 family",
        "Methyl substitution on Cp depresses melting point while keeping %Vbur^D2 = 81.2% within pass band.",
    ),
    (
        "HfCl4",
        "Cl[Hf](Cl)(Cl)Cl",
        "Group 4 (Hf/Ti)",
        "Tetrahalide",
        "Carbon-free HKMG HfO2 ALD",
        "US, KR, TW",
        "DOI:10.1021/cr900056b",
        "Inorganic HKMG benchmark offering zero carbon incorporation.",
    ),
    (
        "TiCl4",
        "Cl[Ti](Cl)(Cl)Cl",
        "Group 4 (Hf/Ti)",
        "Tetrahalide",
        "TiN / TiO2 high-aspect-ratio ALD",
        "US, KR, TW, JP",
        "DOI:10.1021/cr900056b",
        "Industry workhorse for ultra-conformal TiN contact barriers.",
    ),
    # --- Interconnect & Noble Metals (Co, Ni, Cu, Pt) ---
    (
        "CpCo(CO)2",
        "C1=CC=CC1[Co](<-[C-]#[O+])<-[C-]#[O+]",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Half-sandwich Carbonyl",
        "Low-T Co cap & CoSi2 contact ALD",
        "US, JP, KR",
        "DOI:10.1021/cr900056b; Gelest Co interconnect line",
        "18-electron liquid Co(I) precursor bridging CCTBA and Co(AMD)2.",
    ),
    (
        "Co(tBu2DAD)2",
        "CC(C)(C)N1=CC=N(C(C)(C)C)[Co]12N(C(C)(C)C)=CC=N2C(C)(C)C",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Diazadienyl (DAD)",
        "Area-selective Co interconnect ALD",
        "KR, US",
        "KR101719526B1; DOI:10.1021/acsami.1c09819",
        "Redox-active 1,4-diazadiene chelate reacting cleanly with HCOOH or tBuNH2.",
    ),
    (
        "Ni(dmamp)2",
        "CC(C)(CN(C)C)O[Ni]OC(C)(C)CN(C)C",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Aminoalkoxide Chelate",
        "Ni / NiO / NiSi contact ALD",
        "KR, US, JP",
        "DOI:10.1116/1.1875172",
        "First Ni precursor in repo; bidentate aminoalkoxide suppresses oligomerization.",
    ),
    (
        "Ni(iPr-AMD)2",
        "CC(C)N=C(C)N(C(C)C)[Ni]N(C(C)C)C(C)=NC(C)C",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Amidinate",
        "Ni metal / Ni3N / RRAM ALD",
        "US, KR, CN",
        "US7557229B2; WO2004046417",
        "Gordon-type volatile Ni(II) amidinate benchmark.",
    ),
    (
        "Cu(dmap)2",
        "CC(O[Cu]OC(C)CN(C)C)CN(C)C",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Aminoalkoxide Chelate",
        "Low-T low-resistivity Cu interconnect ALD",
        "KR, US",
        "DOI:10.1021/acs.chemmater.7b01710; DOI:10.1063/1.4932057",
        "First Cu precursor in repo; reduces cleanly with tBuNHNH2 or Et2Zn.",
    ),
    (
        "Cu(hfac)(vtms)",
        "C=C[Si](C)(C)C->[Cu]1OC(=CC(=O1)C(F)(F)F)C(F)(F)F",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Beta-diketonate Pi-Adduct",
        "Cu seed CVD / selective inhibitor CVD",
        "US, KR, TW, JP",
        "US20160148839A1; DOI:10.1021/acs.chemmater.0c04771",
        "CupraSelect liquid Cu(I) disproportionation precursor.",
    ),
    (
        "MeCpPtMe3",
        "CC1=CC=CC1[Pt](C)(C)C",
        "Interconnect/Noble (Co/Ni/Cu/Pt)",
        "Cp-Alkyl Mixed",
        "Pt electrode & selective catalytic seed ALD",
        "US, JP, KR",
        "DOI:10.1016/j.apsusc.2022.152750; DOI:10.1116/6.0001352",
        "18-electron mixed Cp-methyl Pt(IV) standard; demonstrates that Me ligands + Cp stabilize Pt against homoleptic alkane failure.",
    ),
    # --- Rare-Earth & Group 12-14 (La, Y, Zn, Ge, Sn, Al) ---
    (
        "La(iPrAMD)3",
        "CC(C)N=C(C)N(C(C)C)[La](N(C(C)C)C(C)=NC(C)C)N(C(C)C)C(C)=NC(C)C",
        "Rare-Earth / Main-Group",
        "Tris-Amidinate",
        "La2O3 / LaSiOx HKMG dipole V_th tuning ALD",
        "US, KR, TW",
        "US9099301B1; DOI:10.1021/cm050756v",
        "First f-block lanthanide in repo for GAA/HKMG dipole threshold voltage control.",
    ),
    (
        "Y(MeCp)3",
        "CC1=CC=CC1[Y](C1=CC=CC1C)C1=CC=CC1C",
        "Rare-Earth / Main-Group",
        "Tris-AlkylCp",
        "Y2O3 / Y-doped HfO2 ferroelectric ALD",
        "WO, US, KR",
        "WO2024050202A1; DOI:10.1021/cm050756v",
        "First Group 3 rare-earth Cp precursor for ferroelectric HZO stabilization.",
    ),
    (
        "ZnEt2 (DEZ)",
        "CC[Zn]CC",
        "Rare-Earth / Main-Group",
        "Dialkyl d10",
        "ZnO / IGZO / AZO channel ALD",
        "US, KR, TW, CN",
        "US20210228302A1; DOI:10.1021/cr900056b",
        "Standard IGZO/ZnO precursor; exposes d10 closed-shell refinement need in homoleptic_tm_sp3 gate.",
    ),
    (
        "Ge(NMe2)4",
        "CN(C)[Ge](N(C)C)(N(C)C)N(C)C",
        "Rare-Earth / Main-Group",
        "Tetra-Alkylamide",
        "GeO2 / GST phase-change / SiGe ALD",
        "US, KR, JP",
        "DOI:10.1116/1.3653424",
        "First Ge precursor in repo; small Ge covalent radius causes %Vbur^D2 = 90.5% steric gate boundary.",
    ),
    (
        "Sn(NMe2)4 (TDMASn)",
        "CN(C)[Sn](N(C)C)(N(C)C)N(C)C",
        "Rare-Earth / Main-Group",
        "Tetra-Alkylamide",
        "SnO2 / EUV resist / OS semiconductor ALD",
        "US, KR, TW",
        "DOI:10.1116/1.2835087",
        "Pairs with Ge(NMe2)4 to map Group 14 ionic-radius relaxation (%Vbur^D2 = 79.2%, Clean PASS).",
    ),
    (
        "TiBA (Al(iBu)3)",
        "CC(C)C[Al](CC(C)C)CC(C)C",
        "Rare-Earth / Main-Group",
        "Trialkyl Al",
        "Low-O Al metal / AlN / AlOx ALD",
        "US, JP",
        "DOI:10.1021/cr900056b",
        "Beta-H-bearing Al precursor complementing TMA by enabling cleaner low-carbon Al metal CVD/ALD.",
    ),
    # --- Patented ASD Inhibitors (10) ---
    (
        "TMSDEA",
        "CCN(CC)[Si](C)(C)C",
        "AS-ALD Inhibitor",
        "Bulky Aminosilane SMI",
        "Dielectric-OH silylation blocking",
        "JP, KR, US, TW",
        "Tokyo Electron / KAIST; KR20200140432A family",
        "Bulkier NEt2 leaving group than DMATMS; tunes kinetic steric coverage on SiO2.",
    ),
    (
        "MSTFA",
        "CN(C(=O)C(F)(F)F)[Si](C)(C)C",
        "AS-ALD Inhibitor",
        "Fluoroacetamide Silane SMI",
        "Rapid low-T OH-group silylation",
        "US, KR, TW",
        "Silylation SMI literature; Applied Materials/Lam ASD context",
        "Electron-withdrawing CF3CO leaves fast at low substrate temperature.",
    ),
    (
        "BSTFA",
        "FC(F)(F)C(=N[Si](C)(C)C)O[Si](C)(C)C",
        "AS-ALD Inhibitor",
        "Dual-TMS Trifluoroacetamide",
        "Ultra-fast oxide passivation SMI",
        "US, JP",
        "Vapor silylation SMI review",
        "Delivers two TMS capping units per reagent with volatile neutral byproduct.",
    ),
    (
        "TBDMSCl",
        "CC(C)(C)[Si](C)(C)Cl",
        "AS-ALD Inhibitor",
        "Bulky Chlorosilane SMI",
        "Hydrophobic bulky SiO2 cap",
        "US, KR",
        "US9564312B2 / US10804099B2 aminosilane/chlorosilane context",
        "tBuMe2Si umbrella provides stronger nucleation delay than Me3SiCl (TMSCl).",
    ),
    (
        "HMDSO",
        "C[Si](C)(C)O[Si](C)(C)C",
        "AS-ALD Inhibitor",
        "Disiloxane Passivant",
        "Amine/Cl-free plasma/vapor SiO2 passivation",
        "US, KR, TW",
        "Clean vapor-phase Si-CH3 passivant",
        "Eliminates HCl and HNMe2 corrosive byproducts during in-situ chamber passivation.",
    ),
    (
        "PPA",
        "O=P(O)(O)c1ccccc1",
        "AS-ALD Inhibitor",
        "Aromatic Phosphonic Acid",
        "Dielectric-on-Dielectric (TiN/Al2O3 vs SiO2) SAM",
        "US, TW, KR, JP",
        "IMEC / TSMC DoD PA SAM family",
        "Pi-stacking phenyl ring increases lateral packing density over alkyl PA.",
    ),
    (
        "MPA",
        "CP(=O)(O)O",
        "AS-ALD Inhibitor",
        "Vapor Alkylphosphonic Acid",
        "All-dry vapor-dosed metal-oxide blocker",
        "US, KR, JP",
        "TEL / Applied Materials CVT PA passivation",
        "Low MW enables dry vapor-phase pulsing in standard ALD manifolds without solvent dips.",
    ),
    (
        "Hhfac",
        "O=C(CC(=O)C(F)(F)F)C(F)(F)F",
        "AS-ALD Inhibitor",
        "Fluorinated Beta-Diketone",
        "Selective blocking & hybrid AS-ALD/ALE",
        "US, JP, TW",
        "Merck / Colorado ALE-ASD hybrid",
        "Chelates and volatilizes stray metal nuclei via volatile M(hfac)x byproducts.",
    ),
    (
        "2-MBI",
        "Sc1nc2ccccc2[nH]1",
        "AS-ALD Inhibitor",
        "Mercapto-Azole SAM",
        "Cu / Co / Ru BEOL FSAV selective blocking",
        "TW, US, KR",
        "TSMC US20250323146A1 / US11854940B2 BEOL FSAV context",
        "Combines thiol-metal anchor + benzimidazole pi-ring for dense Cu/Ru barrier passivation.",
    ),
    (
        "PFDT",
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCS",
        "AS-ALD Inhibitor",
        "Perfluoroalkanethiol SAM",
        "Superhydrophobic Cu/Co metal-vs-oxide blocker",
        "TW, US, JP, KR",
        "TSMC US20250323146A1 / SKKU fluorothiol SAM",
        "Perfluoroalkyl chain maximizes water contact angle and steric precursor repulsion.",
    ),
]


def _ligand_sa_score(
    smi: str,
    target_type: str,
    sa_val: float,
    br_val: float,
    br_scorer: Any,
) -> float:
    """Apply uniform repo SA scoring rule: ligand fragmentation for complexes, BR-SA for small-Si SMIs."""
    mol = Chem.MolFromSmiles(smi)
    if target_type == "inhibitor":
        # Non-metal/silicon-organic inhibitors use min(BR-SA, SA) consistent with DB records (DMATMS/HMDS)
        cand = [sa_val]
        if br_scorer is not None:
            try:
                bs, _ = br_scorer.calculateScore(smi)
                cand.append(float(bs))
            except Exception:
                pass
        return round(min(cand), 3)
    return compute_ligand_sa(smi, sa_val, br_val, br_scorer=br_scorer)


def build_expansion_dataset(
    unified_json_path: str,
    output_json_path: str,
) -> Dict[str, Any]:
    """Validate, screen, deduplicate, and save patent/literature expansion candidates."""
    with open(unified_json_path, "r", encoding="utf-8") as f:
        existing_recs = json.load(f)["records"]

    ex_cans = set()
    ex_names = set()
    ex_prec_el = Counter()
    for r in existing_recs:
        ex_names.add(str(r.get("name", "")).lower())
        m = Chem.MolFromSmiles(r["smiles"])
        if m:
            ex_cans.add(Chem.MolToSmiles(m))
        if r.get("group") in {"Audited Precursor (27)", "Expanded Precursor (75)"}:
            c_el = r.get("center_atom")
            if c_el:
                ex_prec_el[c_el] += 1

    br_scorer = BRSAScorer() if BRSAScorer is not None else None
    records: List[Dict[str, Any]] = []

    for idx, spec in enumerate(EXPANSION_SPECS, start=1):
        name, smi, grp_lbl, fam_cls, app, juris, refs, note = spec
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES in expansion spec: {name} -> {smi}"
        can = Chem.MolToSmiles(mol)
        assert can not in ex_cans and name.lower() not in ex_names, f"Duplicate candidate: {name}"

        ttype = "inhibitor" if grp_lbl == "AS-ALD Inhibitor" else "precursor"
        sa_val = round(float(sascorer.calculateScore(mol)), 3)
        br_val = round(float(br_scorer.calculateScore(smi)[0]), 3) if br_scorer is not None else sa_val
        lig_sa = _ligand_sa_score(smi, ttype, sa_val, br_val, br_scorer)

        topo = compute_intrinsic_topology(smi)
        g = compute_design12_geometry_metrics(smi)
        desc = compute_chemical_coord_descriptors(smi)
        vbur_d2 = float(g["vbur_d2_percent"])
        v_min = universal_vbur_lower_bound(int(topo["pi_hapticity_index"]))
        center_h = mol.GetAtomWithIdx(int(topo["center_idx"])).GetTotalNumHs()

        tmc_sa = compute_tmc_sa_score(
            lig_sa,
            vbur_d2,
            v_min,
            bool(topo["is_coord_center"]),
            is_d_block=bool(topo["is_d_block"]),
            ve_deficit=float(g["ve_deficit"]),
            chi_pauling=float(g["chi_pauling"]),
            center_h=int(center_h),
        )
        gate = evaluate_coord_chem_gate(smi, lig_sa, desc, topo, include_design2=True)
        unified_pass = bool(gate["passed"] and tmc_sa <= 4.0)

        # Boundary / refinement note for honest diagnosis
        boundary_flag = ""
        if topo["center_atom"] == "Zn" and "homoleptic_tm_sp3" in gate["reasons"]:
            boundary_flag = "Zn(II) d10 폐각 착물 — 현 SMARTS 게이트 정교화 대상 (DEZ)"
        elif vbur_d2 > 83.5 and ttype == "precursor":
            boundary_flag = f"소형 중심원자 사면체/5배위 입체밀집 (%Vbur={vbur_d2:.1f}% > 83.5%)"

        records.append(
            {
                "id": f"EXP_{idx:02d}",
                "name": name,
                "smiles": smi,
                "canonical_smiles": can,
                "target_type": ttype,
                "element_group": grp_lbl,
                "family_class": fam_cls,
                "application": app,
                "jurisdictions": juris,
                "references": refs,
                "mechanism_note": note,
                "center_atom": str(topo["center_atom"]),
                "is_d_block": bool(topo["is_d_block"]),
                "cn_eff": round(float(g["cn_eff"]), 1),
                "n_ve": round(float(g["n_ve"]), 1),
                "pi_hapticity": int(topo["pi_hapticity_index"]),
                "chelate_rings": int(topo["chelate_ring_count"]),
                "embed_ok": bool(g["embed_ok"]),
                "vbur_d2_percent": round(vbur_d2, 2),
                "sascore": sa_val,
                "br_sascore": br_val,
                "ligand_sa": lig_sa,
                "tmc_sa_no3d": tmc_sa,
                "gate_passed": bool(gate["passed"]),
                "gate_reasons": gate["reasons"],
                "unified_pass": unified_pass,
                "boundary_flag": boundary_flag,
            }
        )

    exp_prec_el = Counter(r["center_atom"] for r in records if r["target_type"] == "precursor")
    all_els = sorted(set(ex_prec_el.keys()) | set(exp_prec_el.keys()))
    coverage_table = [
        {
            "element": el,
            "current_102": int(ex_prec_el.get(el, 0)),
            "added_28": int(exp_prec_el.get(el, 0)),
            "merged_130": int(ex_prec_el.get(el, 0) + exp_prec_el.get(el, 0)),
        }
        for el in sorted(all_els, key=lambda e: (-(ex_prec_el.get(e, 0) + exp_prec_el.get(e, 0)), e))
    ]

    summary = {
        "n_total": len(records),
        "n_precursors": sum(1 for r in records if r["target_type"] == "precursor"),
        "n_inhibitors": sum(1 for r in records if r["target_type"] == "inhibitor"),
        "n_precursors_pass": sum(1 for r in records if r["target_type"] == "precursor" and r["unified_pass"]),
        "n_inhibitors_pass": sum(1 for r in records if r["target_type"] == "inhibitor" and r["unified_pass"]),
        "n_new_elements_added": sum(1 for el in exp_prec_el if ex_prec_el.get(el, 0) == 0),
        "new_elements": sorted(el for el in exp_prec_el if ex_prec_el.get(el, 0) == 0),
        "n_merged_elements": len(all_els),
        "coverage_table": coverage_table,
    }

    payload = {"summary": summary, "records": records}
    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    with open(output_json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return payload


if __name__ == "__main__":
    u_path = os.path.join(PROJECT_ROOT, "database", "benchmark_unified_latest_160.json")
    o_path = os.path.join(PROJECT_ROOT, "database", "external", "patent_lit_expansion_candidates.json")
    res = build_expansion_dataset(u_path, o_path)
    print("SUMMARY:", json.dumps(res["summary"], ensure_ascii=False))
    for r in res["records"]:
        print(
            f"{r['id']} | {r['target_type']:9s} | {r['name']:20s} | c={r['center_atom']:2s} | "
            f"CN={r['cn_eff']:4.1f} | ve={r['n_ve']:4.1f} | Vbur={r['vbur_d2_percent']:5.1f} | "
            f"ligSA={r['ligand_sa']:4.2f} | TMC={r['tmc_sa_no3d']:4.2f} | pass={r['unified_pass']} | {r['boundary_flag']}"
        )
