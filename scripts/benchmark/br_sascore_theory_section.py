#!/usr/bin/env python3
"""
Foundational Background Theory Tab Builder (``배경이론``)
---------------------------------------------------------
Generates the standalone HTML content for the 4th top-level view tab
(``tab-content-brtheory``, button label ``배경이론``) organized into two
cohesive sections:
  - Section 1 (``#brt-sec1``): Vanilla 2009 SAScore fundamentals, Hf complex vs.
    pro-ligand evaluations, and SMILES graph reconstruction sensitivity
    explaining literature variance (e.g., 5.179 / 5.491 / 6.74).
  - Section 2 (``#brt-sec2``): Unified BR-SAScore deep dive covering origin,
    Eq. (1)-(13) & interactive simulator, eMolecules (231M) + USPTO (~1M)
    331,332-fragment dictionary audit, BRSAScore.py implementation & 2D atom
    contribution maps, Retro* TS1-TS3 benchmarks, and SAPIS domain extension.
"""

from __future__ import annotations

import gzip
import html
import math
import os
import pickle
import re
import sys
from typing import Any, Dict, List, Tuple

import BRSAScore
from rdkit import Chem
from rdkit.Chem import RDConfig, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

sys.path.insert(0, os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore

# Exact literature citations & dataset constants (Chen & Jung, J. Cheminform. 2024, 16:83)
BR_PAPER_DOI = "10.1186/s13321-024-00879-0"
SA_2009_DOI = "10.1186/1758-2946-1-8"
DEFAULT_FRAG_PENALTY = -6.0  # units: dimensionless log-frequency floor
DEFAULT_PENALTY_BUFFER = 1.0  # units: dimensionless normalization buffer
RADIUS_MORGAN = 2  # units: topological bond radius (r = 2)


def compute_vanilla_sascore_issues() -> Dict[str, Any]:
    """
    Evaluate 2009 RDKit vanilla ``sascorer`` across SMILES graph reconstruction
    variants of ``Hf(OMe)3(acac)``, ``Hf(acac)4``, ``HfCp(OMe)(acac)2`` and their
    neutral free pro-ligands (``MeOH``, ``Hacac``, ``HCp``).

    Physics/Algorithm:
        Demonstrates why identical transition-metal complexes yield wide SA
        spreads (0.5-2.5+ units) under identical SMILES vector chemistry depending
        on dative ring closure, eta5-Cp topology, explicit hydrogens (AddHs),
        zwitterionic/ionic charge states, and stereochemical center assignment.
    """
    complex_specs: List[Tuple[str, str, str, str, bool]] = [
        (
            "Hf(OMe)3(acac)",
            "이온 해리형 (Ionic)",
            "CC(=O)[CH-]C(C)=O.C[O-].C[O-].C[O-].[Hf+4]",
            "금속-리간드 공유/배위 결합이 끊어져 유기 이온 빈도만 반영",
            False,
        ),
        (
            "Hf(OMe)3(acac)",
            "개방형 단좌 에놀레이트 (Open)",
            "CO[Hf](OC)(OC)OC(C)=CC(C)=O",
            "acac를 비고리형 단좌(Hf-O 1개)로 표기; 킬레이트 고리 미반영",
            False,
        ),
        (
            "Hf(OMe)3(acac)",
            "6원 킬레이트 배위형 (Chelate ->)",
            "CO[Hf]1(OC)(OC)OC(C)=CC(C)=O->1",
            "카보닐 O->Hf 배위 닫힘으로 6원 금속고리 희귀 조각(+감점) 반영",
            False,
        ),
        (
            "Hf(OMe)3(acac)",
            "킬레이트 + AddHs 미제거 상태",
            "CO[Hf]1(OC)(OC)OC(C)=CC(C)=O->1",
            "명시적 H가 nAtoms 및 Morgan 지문 조각에 합산되어 점수 급등 (+1.62)",
            True,
        ),
        (
            "Hf(acac)4",
            "개방형 4-단좌 (Open)",
            "CC(=O)C=C(C)O[Hf](OC(C)=CC(C)=O)(OC(C)=CC(C)=O)OC(C)=CC(C)=O",
            "비킬레이트 표기로 acac 유기 조각 비중이 커져 금속 패널티 희석",
            False,
        ),
        (
            "Hf(acac)4",
            "8배위 4중 킬레이트 (Chelate ->)",
            "CC1=CC(C)=O->[Hf]234(O1)(OC(C)=CC(C)=O->2)(OC(C)=CC(C)=O->3)OC(C)=CC(C)=O->4",
            "4개 6원 킬레이트 고리 형성 (외부 발표치 5.491과 근사한 5.4854)",
            False,
        ),
        (
            "Hf(acac)4",
            "4중 킬레이트 양쪽성이온 ([Hf-4]/[O+])",
            "CC1=CC(C)=[O+][Hf-4]234(O1)(OC(C)=CC(C)=[O+]2)(OC(C)=CC(C)=[O+]3)OC(C)=CC(C)=[O+]4",
            "화살표(->) 대신 형식전하([O+]/[Hf-4]) 재구성 시 추가 감점",
            False,
        ),
        (
            "HfCp(OMe)(acac)2",
            "개방형 + 시그마(eta1)-Cp",
            "CO[Hf](C1=CC=CC1)(OC(C)=CC(=O)C)OC(C)=CC(=O)C",
            "Cp가 1개의 Hf-C 단일결합으로만 연결되고 acac 킬레이트 미닫힘",
            False,
        ),
        (
            "HfCp(OMe)(acac)2",
            "2중 킬레이트 + eta1-Cp",
            "CO[Hf]12(C3=CC=CC3)(OC(C)=CC(C)=O->1)OC(C)=CC(C)=O->2",
            "acac 2개는 6원 킬레이트로 닫고 Cp는 eta1 단일결합으로 연결",
            False,
        ),
        (
            "HfCp(OMe)(acac)2",
            "2중 킬레이트 + eta5 다중배위 Cp",
            "[CH]12[CH]3[CH]4[CH]5[CH]1[Hf]234567(OC)(OC(C)=CC(C)=O->6)OC(C)=CC(C)=O->7",
            "Cp 5탄소 eta5 결합 재구성 시 브리지/다중고리 조각 증가로 6점대 진입",
            False,
        ),
    ]

    complex_rows: List[Dict[str, Any]] = []
    for comp, rep, smi, note, use_addhs in complex_specs:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES for {comp} ({rep}): {smi}"
        if use_addhs:
            mol = Chem.AddHs(mol)
        sa_val = float(sascorer.calculateScore(mol))  # units: 1-10 SA scale
        n_at = mol.GetNumAtoms()
        size_pen = math.pow(n_at, 1.005) - n_at if n_at > 0 else 0.0
        n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        complex_rows.append(
            {
                "complex": comp,
                "rep": rep,
                "smiles": smi,
                "n_atoms": n_at,
                "size_pen": round(size_pen, 3),
                "n_chiral": n_chiral,
                "sa_2009": round(sa_val, 4),
                "note": note,
            }
        )

    ligand_specs: List[Tuple[str, str, str, str]] = [
        ("OMe- (Methoxide)", "MeOH (메탄올)", "CO", "알콕사이드 출발물질; 원자 2개 소분자"),
        ("acac- (Acetylacetonate)", "Hacac (keto형)", "CC(=O)CC(C)=O", "베타-디케톤 케토 토토머"),
        ("acac- (Acetylacetonate)", "Hacac (enol형)", "CC(=O)/C=C(/C)O", "배위 직전 시스-에놀형 구조"),
        ("Cp- (Cyclopentadienyl)", "HCp (시클로펜타디엔)", "C1=CCC=C1", "공액 5원 고리 리간드 전구 프로톤화물"),
    ]
    ligand_rows: List[Dict[str, Any]] = []
    for ion, pro_name, smi, note in ligand_specs:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid ligand SMILES: {smi}"
        sa_val = float(sascorer.calculateScore(mol))  # units: 1-10 SA scale
        ligand_rows.append(
            {
                "ion": ion,
                "pro_name": pro_name,
                "smiles": smi,
                "n_atoms": mol.GetNumAtoms(),
                "sa_2009": round(sa_val, 4),
                "note": note,
            }
        )

    return {"complex_rows": complex_rows, "ligand_rows": ligand_rows}


def _load_br_pickle_dicts() -> Tuple[Dict[int, float], Dict[int, float], Dict[int, float]]:
    """
    Load the three gzipped fragment score dictionaries shipped with ``BRSAScore``.

    Physics/Algorithm:
        ``BScores_emolecules.pkl.gz`` stores scaled building-block log-frequencies
        (Eq. 9), ``RScores_uspto.pkl.gz`` stores 2^-d reaction-center weighted
        USPTO log-frequencies (Eq. 10), and ``BRScores_uspto_emolecules.pkl.gz``
        stores their coordinate-wise envelope max(BScore_i, RScore_i) (Eq. 11).
    """
    pdir = os.path.join(os.path.dirname(BRSAScore.__file__), "pickle")
    with gzip.open(os.path.join(pdir, "BScores_emolecules.pkl.gz"), "rb") as fb:
        b_scores: Dict[int, float] = pickle.load(fb)
    with gzip.open(os.path.join(pdir, "RScores_uspto.pkl.gz"), "rb") as fr:
        r_scores: Dict[int, float] = pickle.load(fr)
    with gzip.open(os.path.join(pdir, "BRScores_uspto_emolecules.pkl.gz"), "rb") as fbr:
        br_scores: Dict[int, float] = pickle.load(fbr)
    return b_scores, r_scores, br_scores


def _summarize_score_dict(scores: Dict[int, float]) -> Dict[str, Any]:
    """Compute exact key count, extrema, mean, negative ratio, and score bins."""
    vals = list(scores.values())
    n_total = len(vals)
    neg_cnt = sum(1 for v in vals if v < 0.0)
    bins = [0, 0, 0, 0, 0, 0]
    for v in vals:
        if v < -2.0:
            bins[0] += 1
        elif v < -1.0:
            bins[1] += 1
        elif v < 0.0:
            bins[2] += 1
        elif v < 1.0:
            bins[3] += 1
        elif v < 2.0:
            bins[4] += 1
        else:
            bins[5] += 1
    return {
        "count": n_total,
        "min": round(min(vals), 3),
        "max": round(max(vals), 3),
        "mean": round(sum(vals) / n_total, 3),
        "neg_count": neg_cnt,
        "non_neg_count": n_total - neg_cnt,
        "neg_pct": round(100.0 * neg_cnt / n_total, 2),
        "bins": bins,
    }


def compute_br_pickle_statistics() -> Dict[str, Any]:
    """
    Compute empirical statistics and Venn overlap across ``BRSAScore`` pickles.

    Physics/Algorithm:
        Measures exact cardinality and disjoint vs. shared fragment environments
        between eMolecules building blocks (BFrags) and USPTO reaction-center
        neighborhoods (RFrags).
    """
    b_scores, r_scores, br_scores = _load_br_pickle_dicts()
    b_keys = set(b_scores.keys())
    r_keys = set(r_scores.keys())
    br_keys = set(br_scores.keys())

    both_keys = b_keys & r_keys
    b_only_keys = b_keys - r_keys
    r_only_keys = r_keys - b_keys

    b_wins = sum(1 for k in both_keys if b_scores[k] > r_scores[k])
    r_wins = sum(1 for k in both_keys if r_scores[k] > b_scores[k])
    ties = len(both_keys) - b_wins - r_wins

    return {
        "B": _summarize_score_dict(b_scores),
        "R": _summarize_score_dict(r_scores),
        "BR": _summarize_score_dict(br_scores),
        "overlap": {
            "b_only": len(b_only_keys),
            "r_only": len(r_only_keys),
            "both": len(both_keys),
            "union": len(br_keys),
            "b_wins_in_overlap": b_wins,
            "r_wins_in_overlap": r_wins,
            "ties_in_overlap": ties,
        },
    }


def compute_br_worked_examples() -> List[Dict[str, Any]]:
    """
    Run live ``BRSAScore.SAScorer`` and 2009 ``sascorer`` on canonical organic,
    pro-ligand, and semiconductor organometallic molecules with atom SVGs.

    Physics/Algorithm:
        Reproduces the paper's Aspirin walkthrough (Fig. 2) alongside semiconductor
        complexes and pro-ligands, capturing RDKit ``contribution_to_svg`` atom
        heatmaps that expose missing ``-6.0`` radius-2 fragments.
    """
    scorer = BRSAScore.SAScorer()
    specs: List[Tuple[str, str, str, str]] = [
        (
            "Aspirin",
            "CC(=O)Oc1ccccc1C(=O)O",
            "Paper Fig. 2 기준 예시",
            "방향족 치환 3개 환경만 R=2로 평가되며 모두 상용·반응 사전에 존재해 2.09 산출.",
        ),
        (
            "MeOH (Methanol)",
            "CO",
            "단순 알콕사이드 출발물질",
            "원자 2개 소형 분자로 Radius-2 환경이 0개이므로 BR-fragmentScore=0 -> 1.01.",
        ),
        (
            "Hacac (아세틸아세톤, keto)",
            "CC(=O)CC(=O)C",
            "베타-디케톤 프로리간드 (케토)",
            "상용 시약임에도 중앙 메틸렌(C3) 1개가 사전 미등재(-6.0)로 판정되어 5.20으로 과대 패널티(SAPIS에서는 min(BR, SA)=2.17로 복구).",
        ),
        (
            "HCp (Cyclopentadiene)",
            "C1=CC=CC1",
            "시클로펜타디엔 프로리간드",
            "공액 디엔 고리 R=2 환경 5개(-0.93 ~ -2.50) 평균 + 고리 복잡도 합산 -> 3.33.",
        ),
        (
            "TEOS (Si(OEt)4)",
            "CCO[Si](OCC)(OCC)OCC",
            "공유결합형 SiO2 프리커서",
            "Si-O-Et 9개 환경이 USPTO/eMolecules에 일부 등재(-1.60 ~ -2.37) -> 3.82.",
        ),
        (
            "TDMAH (Hf(NMe2)4)",
            "CN(C)[Hf](N(C)C)(N(C)C)N(C)C",
            "대표 Hf-아미드 ALD 프리커서",
            "모든 Hf-N 주변 R=2 환경 5개가 사전 미등재(-6.0)에 걸려 8.62로 오분류.",
        ),
        (
            "Hf(acac)4",
            "CC1=CC(C)=O[Hf]234(O1)(OC(C)=CC(C)=O2)(OC(C)=CC(C)=O3)OC(C)=CC(C)=O4",
            "8배위 킬레이트 Hf 프리커서",
            "배위 R=2 환경 21개 전부가 -6.0 및 다중고리 복잡도 가중 -> 9.04 (Synthesis Hard 오판).",
        ),
    ]

    results: List[Dict[str, Any]] = []
    for name, smi, role_tag, note in specs:
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"Invalid SMILES for {name}: {smi}"
        sa_2009 = float(sascorer.calculateScore(mol))  # units: 1-10 SA scale
        br_score, contribs = scorer.calculateScore(smi)
        br_val = float(br_score)  # units: 1-10 BR-SA scale

        bi: Dict[int, List[Tuple[int, int]]] = {}
        fps = rdMolDescriptors.GetMorganFingerprint(mol, RADIUS_MORGAN, useChirality=True, bitInfo=bi)
        r2_bits = [b for b, vs in bi.items() if vs[0][1] == RADIUS_MORGAN]
        unseen_cnt = sum(1 for val in contribs.values() if float(val) <= -5.999)

        svg_raw = scorer.contribution_to_svg(smi, contribs)
        svg_start = svg_raw.find("<svg")
        svg_clean = svg_raw[svg_start:] if svg_start >= 0 else svg_raw
        svg_clean = (
            svg_clean.replace("\t", " ")
            .replace("width='400px' height='370px'", "style='width:210px; height:155px; max-width:100%; display:block;'")
        )
        svg_clean = re.sub(r"<!--.*?-->", "", svg_clean, flags=re.S)
        svg_clean = re.sub(r"\s+class=['\"](?:bond|atom)-[^'\"]*['\"]", "", svg_clean)
        svg_clean = re.sub(r"\s+", " ", svg_clean).strip()

        results.append(
            {
                "name": name,
                "smiles": smi,
                "role_tag": role_tag,
                "note": note,
                "sa_2009": round(sa_2009, 2),
                "br_sascore": round(br_val, 2),
                "n_atoms": mol.GetNumAtoms(),
                "n_fps_all": len(fps.GetNonzeroElements()),
                "n_r2_envs": len(r2_bits),
                "n_contribs_neg": len(contribs),
                "n_unseen_minus6": unseen_cnt,
                "min_contrib": round(min(contribs.values()), 2) if contribs else 0.0,
                "svg": svg_clean,
            }
        )
    return results


def _build_morgan_fingerprint_2d_svg() -> str:
    """
    Generate a 2D RDKit structural diagram of ``Hf(OMe)3(acac)`` highlighting
    Morgan Fingerprint circular shells around the central ``Hf`` atom
    (``r=0`` center, ``r=1`` 1-bond shell, ``r=2`` 2-bond shell).
    """
    smi = "CO[Hf]1(OC)(OC)OC(C)=CC(C)=O->1"
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    rdDepictor.Compute2DCoords(mol)

    center_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Hf")
    r1_atoms = {nb.GetIdx() for nb in mol.GetAtomWithIdx(center_idx).GetNeighbors()}
    r2_atoms = {
        nb2.GetIdx()
        for n1 in r1_atoms
        for nb2 in mol.GetAtomWithIdx(n1).GetNeighbors()
    } - r1_atoms - {center_idx}

    atom_colors: Dict[int, Tuple[float, float, float, float]] = {
        center_idx: (0.90, 0.22, 0.27, 0.55),
    }
    for idx in r1_atoms:
        atom_colors[idx] = (0.96, 0.62, 0.04, 0.45)
    for idx in r2_atoms:
        atom_colors[idx] = (0.01, 0.52, 0.78, 0.30)

    hl_atoms = [center_idx] + sorted(r1_atoms) + sorted(r2_atoms)
    hl_bonds: List[int] = []
    bond_colors: Dict[int, Tuple[float, float, float, float]] = {}
    for b in mol.GetBonds():
        u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if u == center_idx or v == center_idx:
            hl_bonds.append(b.GetIdx())
            bond_colors[b.GetIdx()] = (0.96, 0.62, 0.04, 0.50)
        elif (u in r1_atoms and v in r2_atoms) or (v in r1_atoms and u in r2_atoms):
            hl_bonds.append(b.GetIdx())
            bond_colors[b.GetIdx()] = (0.01, 0.52, 0.78, 0.40)

    drawer = rdMolDraw2D.MolDraw2DSVG(320, 210)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.5
    opts.highlightRadius = 0.38
    prep = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer.DrawMolecule(
        prep,
        highlightAtoms=hl_atoms,
        highlightBonds=hl_bonds,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    drawer.FinishDrawing()
    svg_raw = drawer.GetDrawingText().replace("svg:", "")
    svg_start = svg_raw.find("<svg")
    svg_clean = svg_raw[svg_start:] if svg_start >= 0 else svg_raw
    svg_clean = (
        svg_clean.replace("\t", " ")
        .replace("width='320px' height='210px'", "style='width:100%; max-width:320px; height:auto; display:block; margin:0 auto;'")
    )
    svg_clean = re.sub(r"<!--.*?-->", "", svg_clean, flags=re.S)
    svg_clean = re.sub(r"\s+class=['\"](?:bond|atom)-[^'\"]*['\"]", "", svg_clean)
    return re.sub(r"\s+", " ", svg_clean).strip()


def _build_cpti_extrapolation_2d_svg() -> str:
    """
    Generate a 2D RDKit structural diagram of ``CpTi(NMe2)3`` color-coded by
    live ``BRSAScore`` atom contributions:
      - Crimson: metal-adjacent unseen environments (``-6.0``)
      - Blue: Cp ring environments extensible by reagent DB (``-1.016 / -1.324 -> -0.47``)
    """
    smi = "CN(C)[Ti](N(C)C)(N(C)C)C1=CC=CC1"
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    rdDepictor.Compute2DCoords(mol)

    scorer = BRSAScore.SAScorer()
    _, contrib = scorer.calculateScore(smi)

    unseen_atoms = {idx for idx, val in contrib.items() if val <= -5.999}
    reagent_atoms = {idx for idx, val in contrib.items() if -5.999 < val < 0.0}

    atom_colors: Dict[int, Tuple[float, float, float, float]] = {}
    for idx in unseen_atoms:
        atom_colors[idx] = (0.90, 0.22, 0.27, 0.48)
    for idx in reagent_atoms:
        atom_colors[idx] = (0.01, 0.52, 0.78, 0.42)

    hl_atoms = sorted(unseen_atoms | reagent_atoms)
    hl_bonds: List[int] = []
    bond_colors: Dict[int, Tuple[float, float, float, float]] = {}
    for b in mol.GetBonds():
        u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if u in reagent_atoms and v in reagent_atoms:
            hl_bonds.append(b.GetIdx())
            bond_colors[b.GetIdx()] = (0.01, 0.52, 0.78, 0.45)
        elif u in unseen_atoms and v in unseen_atoms:
            hl_bonds.append(b.GetIdx())
            bond_colors[b.GetIdx()] = (0.90, 0.22, 0.27, 0.35)

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 210)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.4
    opts.highlightRadius = 0.38
    prep = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer.DrawMolecule(
        prep,
        highlightAtoms=hl_atoms,
        highlightBonds=hl_bonds,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    drawer.FinishDrawing()
    svg_raw = drawer.GetDrawingText().replace("svg:", "")
    svg_start = svg_raw.find("<svg")
    svg_clean = svg_raw[svg_start:] if svg_start >= 0 else svg_raw
    svg_clean = (
        svg_clean.replace("\t", " ")
        .replace("width='300px' height='210px'", "style='width:100%; max-width:300px; height:auto; display:block; margin:0 auto;'")
    )
    svg_clean = re.sub(r"<!--.*?-->", "", svg_clean, flags=re.S)
    svg_clean = re.sub(r"\s+class=['\"](?:bond|atom)-[^'\"]*['\"]", "", svg_clean)
    return re.sub(r"\s+", " ", svg_clean).strip()


def _render_mini_mol_svg(smi: str, highlight_metal_bonds: bool = False) -> str:
    """Render a compact RDKit 2D SVG for the Pro-ligand decoupling diagram."""
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None, f"Invalid SMILES: {smi}"
    rdDepictor.Compute2DCoords(mol)
    hl_atoms: List[int] = []
    hl_bonds: List[int] = []
    atom_colors: Dict[int, Tuple[float, float, float, float]] = {}
    bond_colors: Dict[int, Tuple[float, float, float, float]] = {}
    if highlight_metal_bonds:
        ti_idxs = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Ti"}
        for idx in ti_idxs:
            hl_atoms.append(idx)
            atom_colors[idx] = (0.90, 0.22, 0.27, 0.50)
        for b in mol.GetBonds():
            if b.GetBeginAtomIdx() in ti_idxs or b.GetEndAtomIdx() in ti_idxs:
                hl_bonds.append(b.GetIdx())
                bond_colors[b.GetIdx()] = (0.96, 0.62, 0.04, 0.55)

    drawer = rdMolDraw2D.MolDraw2DSVG(180, 125)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.2
    opts.highlightRadius = 0.38
    prep = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer.DrawMolecule(
        prep,
        highlightAtoms=hl_atoms,
        highlightBonds=hl_bonds,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    drawer.FinishDrawing()
    svg_raw = drawer.GetDrawingText().replace("svg:", "")
    svg_start = svg_raw.find("<svg")
    svg_clean = svg_raw[svg_start:] if svg_start >= 0 else svg_raw
    svg_clean = (
        svg_clean.replace("\t", " ")
        .replace("width='180px' height='125px'", "style='width:100%; max-width:165px; height:auto; display:block; margin:0 auto;'")
    )
    svg_clean = re.sub(r"<!--.*?-->", "", svg_clean, flags=re.S)
    svg_clean = re.sub(r"\s+class=['\"](?:bond|atom)-[^'\"]*['\"]", "", svg_clean)
    return re.sub(r"\s+", " ", svg_clean).strip()


def _build_proligand_decoupling_diagram_html() -> str:
    """Build the 2D structural cleavage diagram for *Pro-ligand decoupling."""
    complex_svg = _render_mini_mol_svg("CN(C)[Ti](N(C)C)(N(C)C)C1=CC=CC1", highlight_metal_bonds=True)
    me2nh_svg = _render_mini_mol_svg("CNC", highlight_metal_bonds=False)
    hcp_svg = _render_mini_mol_svg("C1=CC=CC1", highlight_metal_bonds=False)
    return f"""
      <div style="display:flex; align-items:center; justify-content:center; gap:8px; flex-wrap:wrap;">
        <div style="flex:1 1 135px; text-align:center;">
          {complex_svg}
          <div style="font-size:0.72rem; font-weight:700; color:#0f172a; margin-top:2px;">CpTi(NMe<sub>2</sub>)<sub>3</sub></div>
          <div style="font-size:0.68rem; color:#b91c1c;">금속(Ti) &middot; M-L 결합 절단</div>
        </div>
        <div style="font-size:1.1rem; font-weight:800; color:#0284c7; padding:0 2px;">&rarr;</div>
        <div style="flex:1 1 190px; display:grid; grid-template-columns:1fr auto 1fr; align-items:center; gap:4px;">
          <div style="text-align:center; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:4px;">
            {me2nh_svg}
            <div style="font-size:0.71rem; font-weight:700; color:#0f172a;">Me<sub>2</sub>NH</div>
            <div style="font-size:0.67rem; color:#0f766e;">min(1.02, 2.22) = <strong>1.02</strong></div>
          </div>
          <div style="font-size:0.95rem; font-weight:800; color:#64748b;">+</div>
          <div style="text-align:center; background:#f0fdfa; border:1px solid #5eead4; border-radius:6px; padding:4px;">
            {hcp_svg}
            <div style="font-size:0.71rem; font-weight:700; color:#0f172a;">HCp</div>
            <div style="font-size:0.67rem; color:#0f766e;">min(3.33, 3.31) = <strong>3.31</strong></div>
          </div>
        </div>
      </div>
    """


def _build_vanilla_reconstruction_rows_html(vanilla_data: Dict[str, Any]) -> Tuple[str, str]:
    """Render HTML table rows for Section 1 (Hf complex variants and free pro-ligands)."""
    c_rows: List[str] = []
    for row in vanilla_data["complex_rows"]:
        is_chelate = "5.2957" if row["sa_2009"] == 5.2957 else ("5.4854" if row["sa_2009"] == 5.4854 else ("6.4305" if row["sa_2009"] == 6.4305 else ""))
        row_style = ' style="background:#f0fdf4;"' if is_chelate else (' style="background:#fef2f2;"' if row["sa_2009"] > 6.5 else "")
        c_rows.append(
            f"""
            <tr{row_style}>
              <td><strong>{html.escape(row["complex"])}</strong></td>
              <td>{html.escape(row["rep"])}</td>
              <td style="text-align:center;">{row["n_atoms"]}</td>
              <td style="text-align:right; font-family:'JetBrains Mono', monospace;">{row["size_pen"]:.3f}</td>
              <td style="text-align:right; font-weight:800; color:#0f172a;">{row["sa_2009"]:.3f}</td>
              <td style="font-size:0.8rem; color:#475569;">{html.escape(row["note"])}</td>
            </tr>"""
        )

    l_rows: List[str] = []
    for row in vanilla_data["ligand_rows"]:
        l_rows.append(
            f"""
            <tr>
              <td><strong>{html.escape(row["ion"])}</strong></td>
              <td>{html.escape(row["pro_name"])}</td>
              <td><code>{html.escape(row["smiles"])}</code></td>
              <td style="text-align:center;">{row["n_atoms"]}</td>
              <td style="text-align:right; font-weight:800; color:#0f766e;">{row["sa_2009"]:.3f}</td>
              <td style="font-size:0.8rem; color:#475569;">{html.escape(row["note"])}</td>
            </tr>"""
        )
    return "\n".join(c_rows), "\n".join(l_rows)


def _build_fragment_distribution_svg(stats: Dict[str, Any]) -> str:
    """
    Render a compact vector SVG comparing the 6 histogram bins and Venn overlap
    across ``BScores``, ``RScores``, and ``BRScores``.
    """
    bins_labels = ["[-3, -2)", "[-2, -1)", "[-1, 0)", "[0, 1)", "[1, 2)", "[2, 3]"]
    b_bins = stats["B"]["bins"]
    r_bins = stats["R"]["bins"]
    br_bins = stats["BR"]["bins"]
    max_val = max(max(b_bins), max(r_bins), max(br_bins), 1)

    bars_svg: List[str] = []
    base_y = 195.0
    chart_h = 155.0
    for idx, label in enumerate(bins_labels):
        gx = 66.0 + idx * 94.0
        for s_idx, (cnt, color) in enumerate(
            [(b_bins[idx], "#0284c7"), (r_bins[idx], "#d97706"), (br_bins[idx], "#0f766e")]
        ):
            h = (cnt / max_val) * chart_h
            bx = gx + s_idx * 22.0
            by = base_y - h
            bars_svg.append(
                f'<rect x="{bx:.1f}" y="{by:.1f}" width="18" height="{max(h, 1.5):.1f}" '
                f'rx="2" fill="{color}" opacity="0.88">'
                f'<title>{label}: {cnt:,}</title></rect>'
            )
            if cnt >= 1000:
                lbl_k = f"{cnt / 1000:.0f}k"
                bars_svg.append(
                    f'<text x="{bx + 9:.1f}" y="{by - 4:.1f}" text-anchor="middle" '
                    f'font-size="8.5" font-weight="700" fill="#334155">{lbl_k}</text>'
                )
        bars_svg.append(
            f'<text x="{gx + 29:.1f}" y="213" text-anchor="middle" font-size="10" '
            f'font-weight="600" fill="#475569">{html.escape(label)}</text>'
        )

    ov = stats["overlap"]
    total_u = max(ov["union"], 1)
    w_bar = 540.0
    w_bo = w_bar * ov["b_only"] / total_u
    w_both = w_bar * ov["both"] / total_u
    w_ro = w_bar * ov["r_only"] / total_u
    x0 = 66.0
    y_ov = 258.0

    venn_svg = (
        f'<text x="{x0:.1f}" y="{y_ov - 10:.1f}" font-size="10.5" font-weight="700" fill="#0f172a">'
        f'Morgan Fingerprint (r=2) 사전 집합 분할: 총 {ov["union"]:,}개 고유 조각 = B 전용 + B∩R + R 전용</text>'
        f'<rect x="{x0:.1f}" y="{y_ov:.1f}" width="{w_bo:.1f}" height="24" fill="#0284c7" />'
        f'<rect x="{x0 + w_bo:.1f}" y="{y_ov:.1f}" width="{w_both:.1f}" height="24" fill="#7c3aed" />'
        f'<rect x="{x0 + w_bo + w_both:.1f}" y="{y_ov:.1f}" width="{w_ro:.1f}" height="24" fill="#d97706" />'
        f'<text x="{x0 + w_bo / 2:.1f}" y="{y_ov + 15.5:.1f}" text-anchor="middle" font-size="9.5" font-weight="700" fill="#fff">'
        f'B 전용 {ov["b_only"]:,} ({100 * ov["b_only"] / total_u:.1f}%)</text>'
        f'<text x="{x0 + w_bo + w_both / 2:.1f}" y="{y_ov + 15.5:.1f}" text-anchor="middle" font-size="9.5" font-weight="700" fill="#fff">'
        f'B∩R {ov["both"]:,} ({100 * ov["both"] / total_u:.1f}%)</text>'
        f'<text x="{x0 + w_bo + w_both + w_ro / 2:.1f}" y="{y_ov + 15.5:.1f}" text-anchor="middle" font-size="9.5" font-weight="700" fill="#fff">'
        f'R 전용 {ov["r_only"]:,} ({100 * ov["r_only"] / total_u:.1f}%)</text>'
    )

    return f"""
      <svg viewBox="0 0 650 306" style="width:100%; max-width:680px; height:auto; display:block; margin:12px auto 4px auto; background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px;">
        <text x="66" y="22" font-size="11.5" font-weight="800" fill="#0f172a">BRSAScore==0.1.1 내장 조각 사전 점수 분포 (실측 키 개수 기준)</text>
        <circle cx="410" cy="18" r="5" fill="#0284c7"/><text x="419" y="21.5" font-size="9.5" font-weight="700" fill="#0369a1">BScore ({stats["B"]["count"]:,})</text>
        <circle cx="495" cy="18" r="5" fill="#d97706"/><text x="504" y="21.5" font-size="9.5" font-weight="700" fill="#b45309">RScore ({stats["R"]["count"]:,})</text>
        <circle cx="576" cy="18" r="5" fill="#0f766e"/><text x="585" y="21.5" font-size="9.5" font-weight="700" fill="#0f766e">BR ({stats["BR"]["count"]:,})</text>
        <line x1="60" y1="40" x2="60" y2="195" stroke="#94a3b8" stroke-width="1"/>
        <line x1="60" y1="195" x2="625" y2="195" stroke="#64748b" stroke-width="1.2"/>
        <text x="54" y="45" text-anchor="end" font-size="9" fill="#64748b">{max_val // 1000}k</text>
        <text x="54" y="120" text-anchor="end" font-size="9" fill="#64748b">{max_val // 2000}k</text>
        <text x="54" y="198" text-anchor="end" font-size="9" fill="#64748b">0</text>
        {''.join(bars_svg)}
        {venn_svg}
      </svg>"""


def _build_worked_example_cards_html(examples: List[Dict[str, Any]]) -> str:
    """Render RDKit atom-contribution cards generated by ``contribution_to_svg``."""
    cards: List[str] = []
    for ex in examples:
        delta = round(ex["br_sascore"] - ex["sa_2009"], 2)
        delta_sign = "+" if delta >= 0 else ""
        badge_bg = "#dcfce7; color:#166534;" if ex["br_sascore"] <= 3.5 else (
            "#fef3c7; color:#92400e;" if ex["br_sascore"] <= 6.0 else "#fee2e2; color:#991b1b;"
        )
        cards.append(
            f"""
            <div style="background:#fff; border:1px solid #cbd5e1; border-radius:10px; padding:14px; display:flex; flex-direction:column; gap:8px;">
              <div style="display:flex; justify-content:space-between; align-items:flex-start; gap:8px;">
                <div>
                  <div style="font-weight:800; font-size:0.94rem; color:#0f172a;">{html.escape(ex["name"])}</div>
                  <div style="font-size:0.74rem; color:#0369a1; font-weight:600;">{html.escape(ex["role_tag"])}</div>
                </div>
                <span style="padding:3px 8px; border-radius:999px; font-size:0.74rem; font-weight:800; {badge_bg}">
                  BR: {ex["br_sascore"]:.2f}
                </span>
              </div>
              <div style="display:flex; justify-content:center; align-items:center; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; min-height:152px; overflow:hidden;">
                {ex["svg"]}
              </div>
              <div style="display:grid; grid-template-columns:repeat(3,1fr); gap:6px; font-size:0.75rem; text-align:center;">
                <div style="background:#f1f5f9; padding:5px; border-radius:6px;">
                  <div style="color:#64748b;">SAScore (2009)</div>
                  <strong style="color:#0f172a;">{ex["sa_2009"]:.2f}</strong>
                </div>
                <div style="background:#eff6ff; padding:5px; border-radius:6px;">
                  <div style="color:#1d4ed8;">BR - (2009)</div>
                  <strong style="color:#1e40af;">{delta_sign}{delta:.2f}</strong>
                </div>
                <div style="background:#fef2f2; padding:5px; border-radius:6px;">
                  <div style="color:#b91c1c;">Unseen (-6)</div>
                  <strong style="color:#991b1b;">{ex["n_unseen_minus6"]}/{ex["n_contribs_neg"]}</strong>
                </div>
              </div>
              <div style="font-size:0.74rem; color:#475569; line-height:1.45;">
                {html.escape(ex["note"])}
              </div>
              <div style="font-family:'JetBrains Mono', monospace; font-size:0.71rem; color:#64748b; background:#f8fafc; padding:5px 7px; border-radius:5px; word-break:break-all;">
                {html.escape(ex["smiles"])}
              </div>
            </div>"""
        )
    return "\n".join(cards)


def build_br_sascore_theory_html() -> str:
    """
    Build the complete 2-section ``배경이론`` HTML fragment combining:
      1) Section 1 (``#brt-sec1``): 2009 vanilla SAScore principles, Hf-complex vs.
         pro-ligand evaluation, and SMILES graph reconstruction sensitivity.
      2) Section 2 (``#brt-sec2``): Unified single-section BR-SAScore deep dive
         (equations, dataset pickles, code audit, benchmarks, SAPIS).

    Physics/Algorithm:
        All reconstruction scores, pickle dictionary counts/bins, and atom
        contribution heatmaps are computed live from RDKit ``sascorer`` and
        ``BRSAScore==0.1.1`` on build.
    """
    vanilla_data = compute_vanilla_sascore_issues()
    vanilla_complex_rows_html, vanilla_ligand_rows_html = _build_vanilla_reconstruction_rows_html(vanilla_data)
    morgan_2d_svg = _build_morgan_fingerprint_2d_svg()
    cpti_2d_svg = _build_cpti_extrapolation_2d_svg()
    proligand_diagram_html = _build_proligand_decoupling_diagram_html()

    stats = compute_br_pickle_statistics()
    examples = compute_br_worked_examples()
    dist_svg = _build_fragment_distribution_svg(stats)
    cards_html = _build_worked_example_cards_html(examples)

    b_stat = stats["B"]
    r_stat = stats["R"]
    br_stat = stats["BR"]
    ov = stats["overlap"]

    html_template = r"""
      <!-- SECTION 2: SASCORE (2009) PRINCIPLES, STRUCTURAL LIMITS & SMILES RECONSTRUCTION -->
      <section id="brt-sec1">
        <div class="section-header">
          <div class="section-tag">Chapter 02 &bull; Base Model Deep Dive I</div>
          <h2 class="section-title">2. SAScore (2009) 모델의 원리, 구조적 한계 및 SMILES 재구성 이슈</h2>
          <p class="section-desc" style="margin-top:4px;">반도체 SA 모델을 이해하기 위한 기초 이론 가이드입니다. <code>SAScore (2009)</code>의 수식 규격과 일반 유기합성 및 분자 생성 모델에서 나타나는 구조적 한계, 그리고 유기금속 SMILES 그래프 재구성에 따른 점수 편차 원인을 통합 검증합니다.</p>
        </div>

        <p style="line-height:1.72; color:#334155;">
          Ertl &amp; Schuffenhauer가 제안한 <code>SAScore (2009)</code> 모델은 PubChem 대표 1% 샘플인 934,046종 화합물에서 추출한 Morgan Fingerprint(r&le;2) 출현 빈도 사전 <code>fpscores.pkl.gz</code>와 위상 복잡도 패널티의 선형 결합을 기반으로 합니다. 하지만 정적인 부분구조 출현 친숙도만을 집계하는 설계 특성상 일반 유기 분자 설계는 물론, 이를 반도체용 프리커서에 그대로 적용하면 다음과 같은 세 가지 근본 한계가 드러납니다.
        </p>
        <ul style="margin:6px 0 12px 20px; line-height:1.68; font-size:0.86rem; color:#334155;">
          <li><strong>① 정적 빈도·위상 패널티와 실제 역합성 경로의 괴리:</strong> 상용 카이랄 출발물질 유래 분자에 과도한 복잡도 감점을 부여하는 반면, 분자 생성 모델의 단독 목적 함수로 사용할 때는 단순 방향족·알킬 사슬만 반복 연결하는 보상 해킹이 발생하여 실제 역합성 경로 존재율과 큰 괴리를 보입니다.</li>
          <li><strong>② 사전 미등재 및 바닥값 붕괴:</strong> 유기물 중심의 PubChem 사전에 금속&middot;배위 결합 조각이 등재되어 있지 않아 다수의 분자 조각이 최하점인 <code>-4.0</code>을 부여받으며, 실제보다 합성 난이도가 과대 평가됩니다.</li>
          <li><strong>③ SMILES 그래프 재구성 민감도:</strong> 동일한 프리커서라도 킬레이트 고리 닫힘, $\eta^5$ 배위, <code>AddHs</code> 처리 등 그래프 표현 방식에 따라 추출되는 지문 조각과 고리 패널티가 달라져 산출 점수가 크게 요동칩니다.</li>
        </ul>

        <div class="callout-info" style="margin:10px 0 14px 0; background:#f8fafc; border-color:#cbd5e1;">
          <strong style="font-size:0.94rem; color:#0f172a;">*Morgan Fingerprint</strong>
          <div class="brt-split-grid">
            <ul style="margin:0 0 0 18px; line-height:1.68; font-size:0.84rem; color:#334155;">
              <li><strong>원형 부분구조 해싱 원리:</strong> 분자 그래프에서 각 원자를 중심으로 결합 거리 반경 $r = 0, 1, 2$ 이내 이웃 원자들의 원소 종류&middot;원자가&middot;결합 차수&middot;방향족성&middot;고리 포함 여부를 단계적으로 부호화해 고유 정수 해시 식별자인 <code>bitId</code>로 변환하는 원형 분자 지문입니다.</li>
              <li><strong>반경별 의미와 재구성 민감도:</strong> 우측 그림에서 적색으로 표시된 <code>Radius 0</code>은 중심 원자 자체를 나타내고, 주황색 <code>Radius 1</code>은 직접 결합한 1차 이웃 원자를 나타내며, 청색 <code>Radius 2</code>는 두 결합 떨어진 고리 및 치환기 맥락을 포착합니다. 따라서 동일한 프리커서라도 킬레이트 닫힘이나 $\eta^5$ 배위처럼 그래프 연결 방식을 재구성하면 생성되는 <code>Radius 1~2</code> 해시 식별자가 달라져 사전 조회 점수가 바뀝니다.</li>
            </ul>
            <div style="background:#fff; border:1px solid #cbd5e1; border-radius:8px; padding:10px; text-align:center;">
              __MORGAN_2D_SVG__
              <div style="font-size:0.74rem; color:#475569; margin-top:6px; line-height:1.45;">
                <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#e63946; margin-right:3px;"></span><strong>r=0</strong> 중심(Hf) &nbsp;
                <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#f59e0b; margin-right:3px;"></span><strong>r=1</strong> 1-결합(O) &nbsp;
                <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#0284c7; margin-right:3px;"></span><strong>r=2</strong> 2-결합(C)
              </div>
            </div>
          </div>
        </div>

        <div class="grid-2" style="margin:12px 0;">
          <div class="card">
            <h4 style="margin:0 0 8px 0; color:#1e3a8a;">2.1 SAScore (2009) 모델의 수식과 정규화 규격</h4>
            <div style="font-size:0.85rem; line-height:1.68; color:#334155;">
              <p>$$\text{SAScore}_{\text{raw}} = \text{fragmentScore} - \text{complexityPenalty} + \text{score}_3 \tag{1}$$</p>
              <p>$$\text{complexityPenalty} = (n_{\text{At}}^{1.005} - n_{\text{At}}) + \log_{10}(n_{\text{Chiral}}+1) + \log_{10}(n_{\text{Bridge}}+1) + \log_{10}(n_{\text{Spiro}}+1) + \log_{10}(n_{\text{Macro}}+1) \tag{2}$$</p>
              <p>$$\text{score}_3 = \begin{cases} 0.5 \ln(n_{\text{At}} / |\text{FP}|) & \text{if } n_{\text{At}} > |\text{FP}| \\ 0 & \text{otherwise} \end{cases} \tag{3}$$</p>
              <p>$$\text{SAScore}_{\text{full}} = 11 - \frac{\text{SAScore}_{\text{raw}} - (-4.0) + 1}{2.5 - (-4.0)} \times 9 = 11 - \frac{\text{SAScore}_{\text{raw}} + 5.0}{6.5} \times 9 \tag{4}$$</p>
              <p>$$\text{SAScore}_{\text{final}} = \begin{cases} \min(10,\; 8 + \ln(\text{SAScore}_{\text{full}} - 8)) & \text{if } \text{SAScore}_{\text{full}} &gt; 8.0 \\ \max(1,\; \text{SAScore}_{\text{full}}) & \text{otherwise} \end{cases} \tag{5}$$</p>
              <p style="font-size:0.79rem; color:#64748b; margin-bottom:0;">
                * <code>fpscores.pkl.gz</code> 미등재 희귀·금속 조각은 바닥값 <code>-4.0</code>을 부여받으며, Eq. (3)의 분자 대칭 보정과 Eq. (5)의 상위 구간 로그 압축을 거쳐 $[1, 10]$ 범위로 정규화됩니다. <strong>점수가 $1.0$에 가까울수록 합성이 용이하고 $10.0$에 가까울수록 합성이 난해하며, 일반적으로 $3.5$ 이하를 합성 가능 기준으로 해석합니다.</strong>
              </p>
            </div>
          </div>

          <div class="card" style="background:#fef2f2; border-color:#fecdd3;">
            <h4 style="margin:0 0 8px 0; color:#9f1239;">2.2 정적 빈도·위상 패널티의 구조적 한계 및 역합성 경로 괴리</h4>
            <ul style="font-size:0.84rem; color:#334155; line-height:1.62; padding-left:18px; margin:0 0 10px 0;">
              <li style="margin-bottom:6px;"><strong>상용 카이랄 빌딩블록에 대한 과도한 감점:</strong> 천연물 유래 물질이나 아미노산·당 유도체는 다수의 카이랄 중심 $n_{\text{Chiral}}$과 가교 고리 $n_{\text{Bridge}}$를 포함하여 Eq. (2)의 <code>complexityPenalty</code>에 의해 점수가 <code>6~8</code>점 이상으로 급등합니다. 그러나 실제 유기합성에서는 L-프롤린이나 타타르산과 같은 상용 카이랄 출발물질을 도입하여 1~2단계 커플링 반응만으로 손쉽게 합성되는 경우가 많습니다.</li>
              <li><strong>분자 생성 모델의 보상 해킹 현상:</strong> 강화학습이나 확산 생성 모델의 목적 함수로 <code>SAScore</code>만 단독 최소화하면, 생성기는 PubChem 출현 빈도가 높은 페닐 고리나 선형 알킬 사슬만을 단순 반복 연결합니다. 그 결과 산출되는 <code>SAScore</code>는 <code>1.5~2.0</code>으로 매우 낮지만, 입체 장애와 산화 불안정성 및 결합 형성을 위한 반응 작용기 결여로 인해 실제 합성이 불가능한 분자가 양산됩니다.</li>
            </ul>
            <div style="background:#fff; border:1px solid #fda4af; border-radius:6px; padding:9px 11px; font-size:0.79rem; color:#334155; line-height:1.58;">
              <strong style="color:#9f1239;">실증 검증 &mdash; Gao &amp; Coley (<em>J. Chem. Inf. Model.</em> 2020, 60, 5714&ndash;5723):</strong>
              <ul style="margin:4px 0 0 16px; padding:0;">
                <li>생성 모델이 제안한 수만 종의 분자에 대해 ASKCOS 기반 다단계 역합성 탐색을 수행한 결과, <code>SAScore &le; 2.5</code>로 통과한 분자 가운데 실제 역합성 경로가 발견된 비율은 <strong>42%</strong>에 그쳤습니다.</li>
                <li>나머지 <strong>58%</strong>는 상용 빌딩블록까지 도달하는 정방향 반응 단계가 존재하지 않았으며, 이는 <code>SAScore</code>가 부분구조의 출현 친숙도만 측정할 뿐 실제 실행 가능한 합성 경로의 존재 여부를 담보하지 못함을 입증합니다.</li>
              </ul>
            </div>
          </div>
        </div>

        <div class="card" style="background:#fffbeb; border-color:#fde68a; margin:12px 0;">
          <h4 style="margin:0 0 8px 0; color:#92400e;">2.3 유기금속 프리커서의 SMILES 재구성 이슈</h4>
          <ol style="font-size:0.84rem; color:#334155; line-height:1.62; padding-left:18px; margin:0 0 10px 0;">
            <li><strong>킬레이트(&rarr;)·$\eta^5$-Cp 토폴로지 재구성:</strong> <code>acac</code>를 단좌 개방형으로 두는지 6원 킬레이트로 닫는지(<code>Hf(acac)4</code>: <code>4.193 &rarr; 5.485</code>), <code>Cp</code>를 $\eta^1$(<code>5.432</code>)에서 $\eta^5$ 5중 배위(<code>6.431</code>)로 연결하는지에 따라 <code>fragmentScore</code>와 고리 패널티가 급변합니다. 외부 발표 순위(<code>5.179 &lt; 5.491 &lt; 6.74</code>)도 고배위 재구성 방향과 일치합니다.</li>
            <li><strong><code>AddHs</code> &middot; 입체/전하 처리:</strong> <code>AddHs</code> 미제거 시 수소가 $n_{\text{At}}$와 지문 조각에 합산되어 점수가 급등하고(<code>+1.62</code>), 카이랄 중심 인식 및 이온/양쪽성이온 표기에 따라 조회 조각이 바뀝니다.</li>
          </ol>
          <div style="background:#fff; border:1px solid #fcd34d; border-radius:6px; padding:9px 11px; font-size:0.78rem; color:#334155; line-height:1.58;">
            <strong style="color:#92400e;">재구성 방식별 명시적 계산 비교 (RDKit 실측 항 대입):</strong>
            <ul style="margin:4px 0 0 16px; padding:0;">
              <li><strong>Hf(OMe)<sub>3</sub>(acac) [Open]:</strong> $\text{frag}=-0.206,\; \text{comp}=0.186,\; s_3=0.000\; (n_{\text{At}}\le |\text{FP}|) \implies \text{raw}=-0.392 \implies 11 - \frac{-0.392+5.0}{6.5}\times 9 = \mathbf{4.620}$</li>
              <li><strong>Hf(OMe)<sub>3</sub>(acac) [Chelate &rarr;]:</strong> 고리 닫힘으로 $\text{frag}=-0.694,\; \text{comp}=0.186,\; s_3=0.000 \implies \text{raw}=-0.880 \implies 11 - \frac{-0.880+5.0}{6.5}\times 9 = \mathbf{5.296}$</li>
              <li><strong>Hf(OMe)<sub>3</sub>(acac) [Chelate+AddHs]:</strong> H 16개 포함($n_{\text{At}}=30$) $\text{frag}=-2.616,\; \text{comp}=0.515,\; s_3=0.053 \implies \text{raw}=-3.078 \implies \mathbf{6.918}$</li>
              <li><strong>HfCp(OMe)(acac)<sub>2</sub> [&eta;<sup>5</sup>-Cp]:</strong> $\text{frag}=-1.056,\; \text{comp}=0.644\text{ (ring }0.301\text{ 포함)},\; s_3=0.000 \implies \text{raw}=-1.700 \implies \mathbf{6.431}$</li>
            </ul>
          </div>
        </div>

        <div class="card" style="margin:12px 0;">
          <h4 style="margin:0 0 10px 0; color:#0f172a;">
            2.4 Hf 프리커서 3종 SMILES 재구성별 <code>SAScore (2009)</code> 모델 실측 비교 vs 단독(미배위) 프로리간드($HL_i$) 점수
          </h4>
          <div class="table-wrap" style="margin-bottom:14px;">
            <table>
              <thead>
                <tr>
                  <th>반도체 프리커서 (Precursor)</th>
                  <th>SMILES 그래프 재구성 방식</th>
                  <th style="white-space:nowrap; min-width:96px;">원자 수 ($n_{\text{At}}$)</th>
                  <th>Size Penalty</th>
                  <th>SAScore (2009)</th>
                  <th>재구성 편차 원인 (RDKit 2023.09.6 실측)</th>
                </tr>
              </thead>
              <tbody>
                __VANILLA_COMPLEX_ROWS__
              </tbody>
            </table>
          </div>

          <h4 style="margin:10px 0 8px 0; color:#0f766e; font-size:0.94rem;">
            단독(미배위) 프로리간드($HL_i$) 기준 SAScore (2009) 모델 실측
          </h4>
          <div class="table-wrap">
            <table>
              <thead>
                <tr><th>리간드 음이온</th><th>프로리간드</th><th>SMILES</th><th style="white-space:nowrap; min-width:96px;">원자 수 ($n_{\text{At}}$)</th><th style="white-space:nowrap;">SAScore (2009)</th><th>해설</th></tr>
              </thead>
              <tbody>
                __VANILLA_LIGAND_ROWS__
              </tbody>
            </table>
          </div>
        </div>

        <div class="card" style="background:#f8fafc; margin-top:12px;">
          <h4 style="margin:0 0 8px 0; color:#0f172a;">2.5 Summary</h4>
          <p style="font-size:0.86rem; color:#334155; line-height:1.74; margin:0;">
            <code>SAScore (2009)</code>는 유기 분자에서도 상용 카이랄 빌딩블록과 실제 반응 경로 존재성을 반영하지 못해 58%의 역합성 괴리(Gao &amp; Coley, 2020)와 보상 해킹을 유발하며, 반도체 도메인에서는 <code>MeOH (2.701)</code>, <code>Hacac (2.173~2.676)</code>, <code>HCp (3.307)</code>처럼 널리 쓰이는 리간드를 <code>Hf</code> 중심 원자와 결합해 입력할 때 킬레이트 고리 닫힘, $\eta^5$ 배위 표기, 이온 해리, <code>AddHs</code> 처리 등 그래프 재구성 방식에 따라 동일 화합물의 점수가 <code>3.906</code>에서 <code>6.918</code>까지 크게 요동칩니다. 이는 단일 분자 정적 빈도 기반 모델이 실제 합성 경로보다 <code>SMILES</code> 표기 관습과 PubChem 유기 조각 결손에 좌우됨을 시사하며, 따라서 Chapter 3의 상용 시약·반응 분리형 <code>BR-SAScore (2024)</code> 모델 도입과 금속-리간드 분리 기반의 <strong>Pro-ligand + 물리 Descriptor(<code>SAPIS</code>)</strong> 결합 체계로의 전환이 필수적입니다.
          </p>
        </div>
      </section>

      <!-- SECTION 3: BR-SASCORE UNIFIED SINGLE SECTION -->
      <section id="brt-sec2">
        <div class="section-header">
          <div class="section-tag">Chapter 03 &bull; Base Model Deep Dive II</div>
          <h2 class="section-title">3. BR-SAScore (2024) 모델의 원리 및 베이스 모델 선정 근거</h2>
        </div>

        <p style="line-height:1.72; color:#334155;">
          <code>BR-SAScore (2024)</code> 모델은 역합성의 두 핵심 축인 eMolecules 상용 출발물질 $B$와 USPTO 특허 반응 형성 결합 $R$을 분리하여 평가하는 3세대 하이브리드 모델입니다. 다른 세대의 모델과 달리 미지원 원소가 입력되어도 특성 벡터 소실이나 단일 미지 노드 치환 없이 결정론적으로 동작하고, 흔한 알킬 치환기에 난해한 핵심 구조가 희석되는 현상을 음수 병목 평균으로 방지합니다. 또한 평가 대상 프리커서의 정답 조각은 배제한 채 실제 구매 가능한 출발 시약 29종의 조각만 <code>_fscores</code> 해시 사전의 음수 영역 하한인 <code>max(existing, -0.47)</code>로 상향 보정함으로써, 신경망 재학습이나 데이터 누설 없이 유기 프로리간드 및 주족 프리커서 평가 도메인을 확장할 수 있어 본 연구의 베이스 모델로 선정되었습니다.
        </p>

        <!-- 3.1 Architecture Flowchart -->
        <h3 id="brt-sub21" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:18px 0 8px 0;">
          3.1 B/R 분리 아키텍처 파이프라인
        </h3>
        <div class="mermaid-card clickable-card" id="card-brt-pipeline" onclick="openDiagramModal('brt-pipeline-container', 'BR-SAScore (2024) B/R 분리 아키텍처 파이프라인')" title="클릭하여 팝업 확대 보기">
          <div class="mermaid-toolbar">
            <span class="mermaid-title">다이어그램: BR-SAScore (2024) B/R 분리 아키텍처 파이프라인 (클릭 시 확대)</span>
            <span class="click-to-expand-badge">클릭하여 팝업 확대</span>
          </div>
          <div class="mermaid-viewport" id="viewport-brt-pipeline">
            <div class="mermaid-svg-container" id="brt-pipeline-container" style="width:100%;">
              <div class="mermaid" id="brt-pipeline-svg">
flowchart TD
  M["Target Molecule m"] --> MFP["Morgan Fingerprint (r=2, useChirality=True)"]
  M --> COMP["ComplexityPenalty (Size, Stereo, Ring, Macro)"]
  MFP --> B["BScore Dict: eMolecules 231M (244,572 keys) (Eq. 9)"]
  MFP --> R["RScore Dict: USPTO ~1M, 2^-d (175,186 keys) (Eq. 10)"]
  B --> MAX["Fusion: max(B_i, R_i) else -6 (Eq. 11)"]
  R --> MAX
  MAX --> NEG["Negative Mean: avg over Score_i &lt; 0 (Eq. 12)"]
  NEG --> RAW["Raw = BR_frag - Comp (Eq. 8)"]
  COMP --> RAW
  RAW --> NORM["Scale [1, 10] (Buffer=1.0, Div=7) (Eq. 13)"]
              </div>
            </div>
          </div>
        </div>

        <!-- 3.2 Equations (8)-(13) -->
        <h3 id="brt-sub22" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:20px 0 8px 0;">
          3.2 핵심 수식 체계
        </h3>
        <div class="grid-2" style="margin-top:8px;">
          <div class="card" style="border-color:#93c5fd; background:#f8fbff;">
            <div style="font-size:0.85rem; line-height:1.68; color:#334155;">
              <p>$$\text{BR-SAScore} = \text{BR-fragmentScore} - \text{complexityPenalty} \tag{8}$$</p>
              <p>$$\text{BScore}_i = \text{Scale}_{[-3,3]}\!\left(\log_{10} \frac{\text{count}_i^{(B)}}{0.001 \cdot N_B}\right) \tag{9}$$</p>
              <p>$$\text{RScore}_i = \text{Scale}_{[-3,3]}\!\left(\log_{10} \frac{\text{count}_i^{(R)}}{0.001 \cdot N_R}\right), \quad w(d) = 2^{-d} \tag{10}$$</p>
              <p>$$\text{Score}_i = \begin{cases} \max(\text{BScore}_i, \text{RScore}_i) & \text{if } i \in B \cup R \\ -6 & \text{otherwise (Unseen)} \end{cases} \tag{11}$$</p>
              <p>$$\text{BR-fragmentScore} = \frac{1}{|\{k : \text{Score}_k &lt; 0\}|} \sum_{\text{Score}_k &lt; 0} \text{Score}_k \tag{12}$$</p>
              <p>$$\text{BR-SAScore}_{\text{norm}} = 10 - 9 \cdot \frac{\text{BR-SAScore} - (-6 - \text{PenaltyBuffer})}{0 - (-6 - \text{PenaltyBuffer})} \tag{13}$$</p>
            </div>
          </div>

          <div class="callout-info" style="margin:0;">
            <strong>수식 설계의 3가지 핵심 논리</strong>
            <ul style="margin:8px 0 0 18px; line-height:1.68; font-size:0.84rem;">
              <li><strong>최댓값 선택 논리 ($\max(\text{BScore}_i, \text{RScore}_i)$):</strong> 특정 부분구조가 상용 시약 카탈로그에서 직접 구매 가능하거나 특허 화학 반응을 통해 손쉽게 형성될 수 있다면, 두 경로의 점수인 <code>BScore</code>와 <code>RScore</code> 가운데 더 유리한 값을 채택합니다.</li>
              <li><strong>병목 중심 음수 평균 (Eq. 12) 및 상위 0.1% 기준선의 근거:</strong> Eq. (9)~(10)에서 점수 <code>0</code>은 전체 데이터셋 출현 빈도의 <strong>상위 0.1% 기준선</strong>입니다. 분자 지문 조각의 출현 빈도는 극단적인 <strong>멱법칙 분포</strong>를 따르기 때문에 전체 사전 조각의 99% 이상이 수 회 이하로만 등장하는 희귀 조각에 해당합니다. 만약 이 기준선을 상위 10% 수준으로 완화하면 수 회밖에 등장하지 않은 드문 조각들까지 모두 합성이 쉬운 양수 점수 영역으로 분류되어 모델의 난이도 변별력이 무너집니다. 이러한 분포 특성을 고려하여 원조 <code>SAScore (2009)</code>는 전체 누적 출현 횟수의 약 80%를 차지하는 최상위 0.1% 핵심 조각군을 0점 기준선으로 설정하였으며, <code>BR-SAScore</code> 역시 이 경험적 기준을 그대로 계승하였습니다. 이에 따라 $\text{Score}_i \ge 0$인 조각은 메틸기처럼 상용 시약과 특허 반응에 흔히 존재하여 합성에 제약을 주지 않는 부분구조를 의미하고, $\text{Score}_i &lt; 0$인 조각은 출현 빈도가 낮아 전체 합성의 병목으로 작용하는 난해한 부분구조를 의미합니다. 흔한 조각들까지 모두 평균 계산에 포함하면 소수의 난해한 핵심 조각이 지닌 난이도가 희석되어 전체 분자를 합성하기 쉽다고 오판할 수 있으므로, 0점 미만인 병목 조각만을 선별하여 평균을 산출함으로써 난해한 부분구조에 대한 민감도를 극대화합니다. 분자 내 모든 조각의 점수가 0점 이상인 경우에는 합성 병목이 없는 것으로 판단하여 조각 평균 점수를 0으로 처리합니다.</li>
              <li><strong>미등재 패널티 <code>-6</code> &amp; <code>PenaltyBuffer = 1</code> (Eq. 13):</strong> 사전에 등재되지 않은 미지 조각에는 사전 내 최소 점수인 $-3.0$의 두 배에 해당하는 <code>-6.0</code>을 부여하고, 위상 복잡도 감점이 더해질 때 원점수 하한이 $-7.0$까지 확보되도록 <code>PenaltyBuffer = 1</code>을 두어 분모 $7$을 기준으로 $[1, 10]$ 구간에 정규화합니다.</li>
            </ul>
          </div>
        </div>

        <!-- 3.3 Training Datasets & Pickle Audit -->
        <h3 id="brt-sub23" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          3.3 학습&middot;구축 데이터셋 해부
        </h3>
        <div class="table-wrap" style="margin:10px 0;">
          <table>
            <thead>
              <tr>
                <th>데이터셋 / 파일</th>
                <th>원천 데이터 규모</th>
                <th>조각 추출 &amp; 가중 규칙</th>
                <th style="white-space:nowrap; min-width:125px;">최종 사전 키 수</th>
                <th>점수 범위 / 음수 비율</th>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td><strong>B-Dataset</strong><br><code>BScores_emolecules.pkl.gz</code></td>
                <td><strong>2.31억 종 (231M)</strong><br>eMolecules 무료 카탈로그</td>
                <td>분자 내 총 $k$개 R=2 조각에 $1/k$ 정규화<br>단일 출현($\le 1$) 제거</td>
                <td style="text-align:right; font-weight:800; color:#0369a1;">__B_COUNT__</td>
                <td>[__B_MIN__, __B_MAX__]<br>평균 __B_MEAN__ (음수 __B_NEG_PCT__%)</td>
              </tr>
              <tr>
                <td><strong>R-Dataset</strong><br><code>RScores_uspto.pkl.gz</code></td>
                <td><strong>약 100만 건 (~1M)</strong><br>Lowe USPTO 특허 반응</td>
                <td>반응 중심($d=0$)~2결합 거리($d \in \{0,1,2\}$)에 $w = 2^{-d}$ 가중<br>단일 출현($\le 1$) 제거</td>
                <td style="text-align:right; font-weight:800; color:#b45309;">__R_COUNT__</td>
                <td>[__R_MIN__, __R_MAX__]<br>평균 __R_MEAN__ (음수 __R_NEG_PCT__%)</td>
              </tr>
              <tr style="background:#f0fdf4;">
                <td><strong>통합 BR 사전</strong><br><code>BRScores_uspto_emolecules.pkl.gz</code></td>
                <td><strong>eMolecules + USPTO 융합</strong><br>(논문 보고값 <code>331,332</code> 일치)</td>
                <td>$\text{BR}[i] = \max(B[i], R[i])$<br>겹침 __OV_BOTH__개 중 B 우세 __OV_B_WINS__ vs R 우세 __OV_R_WINS__</td>
                <td style="text-align:right; font-weight:900; color:#0f766e;">__BR_COUNT__</td>
                <td>[__BR_MIN__, __BR_MAX__]<br>평균 __BR_MEAN__ (음수 __BR_NEG_PCT__%)</td>
              </tr>
            </tbody>
          </table>
        </div>
        <div class="card" style="background:#fff; padding:14px 16px;">
          <p style="font-size:0.84rem; color:#475569; margin:0 0 6px 0;">
            전체 <strong>__BR_COUNT__개</strong> 조각 중 <strong>__BR_NEG_CNT__개(__BR_NEG_PCT__%)가 0 미만</strong>이며, 전체 키의 <strong>__OV_R_PCT__%(__OV_R_ONLY__개)</strong>는 상용 빌딩블록($B$) 없이 화학 반응($R$)을 통해서만 공급되는 환경입니다.
          </p>
          __DIST_SVG__
        </div>

        <!-- 3.4 Code Audit & Worked Examples -->
        <h3 id="brt-sub24" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          3.4 소스코드 구현 해부 및 원자기여도 분석
        </h3>
        <div class="grid-2" style="margin:10px 0;">
          <div class="card">
            <h4 style="margin:0 0 8px 0; color:#0f172a;">BRSAScore.py 핵심 4대 디테일 (Lines 30~139)</h4>
            <ol style="font-size:0.83rem; color:#334155; line-height:1.65; padding-left:18px; margin:0;">
              <li><strong>Radius-2 엄격 필터 (<code>if vs[0][1] != 2: continue</code>):</strong> 반경 $0, 1$을 배제하므로 <code>MeOH (CO)</code>처럼 반경 2 환경이 0개인 소분자는 <code>BR-fragmentScore = 0</code>이 되어 <code>1.01</code>로 산출됩니다.</li>
              <li><strong>원자 기여도 맵 (<code>contribution_to_svg</code>):</strong> 중심 원자별 <code>fscore</code>를 기록해 <code>-contribution/6</code> 비율로 <code>OrRd</code> 컬러맵 음영을 씌워 희귀(-6.0) 병목 원자를 즉시 시각화합니다.</li>
              <li><strong><code>multicyclePenalty</code> &amp; 거대고리:</strong> 2개 이상 고리 공유 원자 <code>log10(nMultiRingAtoms + 1)</code>과 <code>len(ring) &gt; 6</code> 항이 <code>complexityPenalty</code>에 합산됩니다.</li>
              <li><strong>대칭 보정 <code>score3</code>:</strong> <code>nAtoms &gt; len(fps)</code>일 때 <code>0.5 * log(nAtoms / len(fps))</code>를 더해 대칭 분자(<code>TEOS</code>, <code>TDMAH</code>)의 반복 조각 중복 부담을 일부 완화합니다.</li>
            </ol>
          </div>
          <div class="card" style="background:#f8fafc;">
            <h4 style="margin:0 0 8px 0; color:#0f172a;">논문 Fig. 2 기준 검증: Aspirin vs AI 분자 [Gao et al.]</h4>
            <ul style="font-size:0.83rem; color:#334155; line-height:1.66; margin:0; padding-left:18px;">
              <li><strong>Aspirin:</strong> 유효 반경-2 원자 3개 $(-0.82 - 0.98 - 0.23)/3 = -0.68$, 복잡도 $0.17$ &rarr; Raw $-0.85$ &rarr; <strong>정규화 2.09</strong> (논문 = 본 저장소 실측 2.09 정확 일치).</li>
              <li><strong>AI-proposed 분자 [Gao et al., 5]:</strong> 평가 조각 전부 $B \cup R$ 미등재(<code>-6.0</code>), 복잡도 $0.10$ &rarr; Raw $-6.10$ &rarr; <strong>8.84 (Hard)</strong>.</li>
              <li><strong>케토-에놀 민감도 실측 (한계 사례):</strong> <code>Hacac (아세틸아세톤, keto)</code>는 실제 상용 시약으로 <code>SAScore (2009)</code> 모델에서 <code>2.17</code>이지만, <code>BR-SAScore</code>에선 중앙 <code>-CH2-</code> 1개가 사전 미등재(<code>-6.0</code>)에 걸려 <code>5.20</code>으로 과대 평가되는 한계를 보입니다(<code>compute_ligand_sa()</code>에서 <code>min(BR, SA)=2.17</code>로 복구).</li>
            </ul>
          </div>
        </div>
        <div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(255px,1fr)); gap:12px; margin-top:10px;">
          __CARDS_HTML__
        </div>

        <!-- 3.5 Why BR-SAScore Beats Competitors & Semiconductor Domain Bridge -->
        <h3 id="brt-sub25" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          3.5 베이스 모델로 BR-SAScore (2024) 선정 근거
        </h3>
        <p style="line-height:1.72; color:#334155;">
          완제품 금속 프리커서 분자를 기존 유기 SA 모델에 그대로 입력하면 금속&ndash;리간드 결합에 대한 물리화학적 이해 없이 사전 미등재로 기각되거나 단순 치환기 개수만 보고 통과됩니다. 예를 들어 바닐라 <code>BR-SAScore</code>가 불안정 대조군인 <code>TiEt4</code>(<code>8.77</code>)나 <code>WPh5</code>(<code>8.79</code>)에 높은 난이도 점수를 부여하는 현상 역시 $\beta$-수소 탈락이나 열역학적 불안정성을 모델이 알고 차단한 결과가 아니라, 정상 양산 프리커서인 <code>TDMAH</code>(<code>8.62</code>)나 <code>Hf(acac)4</code>(<code>9.04</code>)와 똑같이 금속&ndash;리간드 반경-2 조각이 유기 사전에 등재되어 있지 않아 일괄적으로 <code>-6.0</code> 감점을 받았기 때문입니다. 따라서 본 연구의 <code>SAPIS</code> 체계에서 금속 배위 결합의 열역학·입체 안정성 판정은 Track 2의 <code>SevenNet-0</code> 물리 지표가 전담하며, Track 1 베이스 SA 모델의 고유 역할은 <strong>분리된 유기 프로리간드 $HL_i$, 주족 공유결합 실란, 및 유기 인히비터의 실제 시약 조달 가능성과 유기 합성 난이도를 왜곡 없이 정량화하는 것</strong>입니다. 이 기준에서 18종 후보 모델을 대조했을 때 <code>BR-SAScore</code>가 최적의 Track 1 베이스 엔진으로 선정된 <strong>4가지 핵심 근거</strong>는 다음과 같습니다.
        </p>

        <div class="grid-2" style="margin:12px 0;">
          <div class="card" style="background:#f8fbff; border-color:#bfdbfe;">
            <h4 style="margin:0 0 8px 0; color:#1e40af;">BR-SAScore 베이스 선정 4대 핵심 근거</h4>
            <ol style="font-size:0.83rem; color:#334155; line-height:1.68; padding-left:18px; margin:0;">
              <li style="margin-bottom:6px;"><strong>상용 출발 시약($B$)과 합성 반응 결합($R$)의 분리 인식:</strong> 원조 <code>SAScore (2009)</code>는 의약품 위주인 PubChem 내 단일 출현 빈도만 조회하므로, 실제로 값싸게 시판되는 기초 시약이나 소분자 리간드인 <code>MeOH</code>(<code>2.70</code>), <code>Me2NH</code>(<code>2.59</code>)에도 감점을 매겨 반도체 출발 시약 29종 통과율이 <strong>31.0%(9/29)</strong>에 그칩니다. 반면 <code>BR-SAScore</code>는 eMolecules 상용 시약 사전과 USPTO 반응 중심 사전을 분리 조회하여 시판 시약·리간드 원형(<code>MeOH 1.01</code>, <code>Me2NH 1.02</code>, 시약 통과율 기본 <code>72.4%</code>에서 도메인 시약 반영 후 <strong><code>96.6%(28/29)</code></strong>)과 다단계 반응이 필요한 구조를 명확히 구분합니다. 한편 4세대 경로 생성 모델(<code>RxnFlow</code>&middot;<code>SynFormer</code>&middot;<code>ChemProjector</code>&middot;<code>ReaSyn</code>)은 내장된 72~115종 유기 반응 템플릿 밖의 규소·헤테로원자 리간드에서 경로 도출 자체가 중단되어 임의의 후보 리간드에 대한 연속적인 난이도 정량 지표로 사용할 수 없습니다.</li>
              <li style="margin-bottom:6px;"><strong>벌키 알킬 치환기에 의한 병목 희석 방지 (Eq. 12 음수 조각 평균):</strong> 반도체 프리커서용 리간드(아미디네이트, 구아니디네이트, 치환 시클로펜타디에닐, 실릴아민)는 기화성과 입체 차폐 확보를 위해 다수의 메틸·에틸·이소프로필·<em>tert</em>-뷰틸 사슬을 포함합니다. 전체 조각 산술 평균을 쓰는 <code>SAScore</code>나 지문 빈도 합산에 의존하는 <code>SCScore</code>&middot;<code>RAscore</code>는 이처럼 흔한 알킬 조각이 늘어날수록 리간드 중심 골격에 합성 불가능한 비정상 고리(<code>Silacyclobutadiene</code> 등)나 난해한 결합이 포함되어 있어도 점수가 쉽게 희석됩니다. 반면 <code>BR-SAScore</code>는 흔한 알킬 조각($\text{Score}_k \ge 0$)을 평균 분모에서 제외하고 <strong><code>0</code> 미만인 합성 병목 조각만 선별하여 평균</strong>하므로, 알킬 치환기 개수에 가려지지 않고 리간드 핵심 작용기의 실제 합성 난이도만 일관되게 평가합니다.</li>
              <li style="margin-bottom:6px;"><strong>비유기 원소 무결성 및 신경망 재학습 없는 국소 사전 확장성:</strong> 주족 원소(<code>Si</code>, <code>Ge</code>, <code>B</code>, <code>P</code>, <code>Sn</code>) 공유결합 프리커서나 실릴 치환 리간드(<code>-SiMe3</code> 포함 리간드 등)를 입력할 때, 그래프 신경망 계열인 <code>GASA</code>는 미지원 원소 8종(102종 중 53종)의 원자 특성 벡터를 영벡터로 소실시키고 <code>DFRscore</code>는 모든 비유기 원소를 단일 미지 노드(<code>ELSE</code>)로 통합하여 원소 구별력을 상실합니다. 또한 <code>SCScore</code>, <code>RAscore</code>, <code>DeepSA</code>, <code>GASA</code>에 반도체 상용 시약 29종을 반영하려면 전체 데이터셋 재구축과 가중치 재학습이 필요해 기존 유기 분자의 결정 경계까지 흔들리는 반면, <code>BR-SAScore</code>는 모든 원소에 대해 결정론적 해시 키를 유지하며 <code>_fscores</code> 사전에서 해당 시약 조각 하한만 <code>max(existing, -0.47)</code>로 국소 보정하여 <strong>기존 프리커서 점수 악화 0건</strong>으로 도메인 시약 지식을 투명하게 편입합니다.</li>
              <li><strong>원자 단위 시약($B$)·반응($R$) 병목 원인 분해 및 <code>SAScore</code>와의 상호보완성:</strong> 분자 대칭 보정 항(<code>score3</code>)은 원조 <code>SAScore (2009)</code>에도 동일하게 포함되어 있으며, 프로리간드를 분리해 Track 2(<code>SevenNet-0</code>) 물리 지표와 결합하는 구조 역시 원리상 다른 유기 SA 모델에도 적용할 수 있습니다. 그럼에도 블랙박스 스칼라 출력 모델(<code>SCScore</code>, <code>RAscore</code>, <code>DeepSA</code>)이나 단일 빈도만 보여 주는 <code>SAScore</code> 대비 <code>BR-SAScore</code>가 지니는 고유한 강점은 <strong>각 원자의 반경-2 환경마다 시약 존재 점수(<code>BScore</code>)와 반응 형성 점수(<code>RScore</code>)를 개별 역추적</strong>할 수 있다는 점입니다. 이를 통해 후보 리간드의 난이도가 높을 때 상용 출발물질이 부재한 것인지 연결 반응 전례가 드문 것인지를 즉시 구분할 수 있으며, 반경-2 엄격 필터로 인해 <code>Hacac</code> 케토형(<code>5.20</code>) 같은 소분자 호변이성질체가 과대 감점되는 단독 한계는 <code>compute_ligand_sa()</code>에서 원조 <code>SAScore</code>(<code>2.17</code>)와 <code>min(BR-SAScore, SAScore)</code>로 결합하여 보완합니다.</li>
            </ol>
          </div>

          <div class="card" style="background:#fff;">
            <h4 style="margin:0 0 8px 0; color:#0f172a;">1G~4G 후보 모델 vs BR-SAScore 실측 비교 요약</h4>
            <div class="table-wrap">
              <table>
                <thead>
                  <tr>
                    <th>모델 (세대)</th>
                    <th>102종 프리커서 / 29종 시약</th>
                    <th>Neg 10종 오탐(FPR)</th>
                    <th>Track 1 베이스 모델 탈락 사유 및 비교</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><strong>SAScore (2009)</strong> (1G)</td>
                    <td>2/102 (2.0%)<br>9/29 (31.0%)</td>
                    <td>3/10 (30.0%)</td>
                    <td>PubChem 단일 빈도로 상용 시약 69.0% 기각(<code>MeOH 2.70</code>); 전체 조각 평균으로 알킬 희석 발생(단, 호변이성질체 보조용으로 결합)</td>
                  </tr>
                  <tr>
                    <td><strong>SYBA</strong> (1G)</td>
                    <td>83/102 (81.4%)<br>27/29 (93.1%)</td>
                    <td><strong>8/10 (80.0%)</strong></td>
                    <td>ZINC 기반 로그오즈 합산이 시클로펜타디에닐 고리에 극단적 음수 감점을 부여하고 불안정 구조를 과도 통과시킴</td>
                  </tr>
                  <tr>
                    <td><strong>SCScore / RAscore</strong> (2G)</td>
                    <td>62~90/102<br>29/29 (100%)</td>
                    <td><strong>9/10 (90.0%)</strong></td>
                    <td>블랙박스 단일 스칼라 출력으로 시약·반응 병목 원인 분해 불가 및 재학습 없는 도메인 시약 보정 불가</td>
                  </tr>
                  <tr>
                    <td><strong>MPScore</strong> (2G)</td>
                    <td><strong>0/102 In-Domain</strong><br>9/29 (31.0%)</td>
                    <td>N/A (유기 전용)</td>
                    <td>유기 다공성 케이지·고분자 빌딩블록 학습 편향으로 반도체 리간드·실란 시약 커버리지 부족</td>
                  </tr>
                  <tr>
                    <td><strong>GASA / DFRscore / DeepSA</strong> (3G)</td>
                    <td>GASA 53/102 미지원<br>DFR 102/102 ELSE</td>
                    <td>실측 불능/편향</td>
                    <td>미지원 원소의 영벡터 소실(<code>GASA</code>), 비유기 원소의 단일 <code>ELSE</code> 노드 통합(<code>DFRscore</code>), 신경망 재학습 필요</td>
                  </tr>
                  <tr>
                    <td><strong>RxnFlow / SynFormer / ReaSyn</strong> (4G)</td>
                    <td><strong>0/102 M-L 규칙</strong><br>(유기 치환기만 매칭)</td>
                    <td>N/A (0% 경로)</td>
                    <td>72~115종 유기 반응 템플릿 밖의 규소·헤테로원자 리간드에서 경로 생성 중단; 리간드 합성 경로 탐색 보조로만 분업</td>
                  </tr>
                  <tr style="background:#ecfdf5; font-weight:700;">
                    <td><strong>BR-SAScore (3G 베이스 선정)</strong></td>
                    <td><strong>19/102 (18.6%)<br>Strict 28/29 (96.6%)</strong></td>
                    <td><strong>3/10 (30.0%)</strong><br><span style="font-weight:400; font-size:0.74rem; color:#475569;">(금속 7종 기각은 M-L 미등재 탓)</span></td>
                    <td><strong>시약(B)·반응(R) 원자별 분리 진단 &middot; 음수 병목 평균으로 알킬 희석 방지 &middot; 무재학습 <code>max(existing,-0.47)</code> 국소 시약 확장(악화 0건)</strong></td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
        </div>

        <div class="grid-3" style="margin:12px 0 4px 0;">
          <div class="card">
            <h4 style="margin:0 0 6px 0; color:#b91c1c;">베이스 단독 한계 1: M-L 조각 0% 중첩</h4>
            <p style="font-size:0.83rem; color:#475569; margin:0; line-height:1.6;">
              <code>BR-SAScore</code> 역시 유기 데이터셋(eMolecules+USPTO) 기반이라 <code>Hf-N</code>, <code>Hf-O</code>, <code>Zr-O</code>, <code>Mo-Cl</code> 반경-2 환경이 없어 프리커서 분자 그대로 입력하면 27종 통과율이 <strong>11.1% (3/27)</strong>에 머무릅니다(<code>TDMAH 8.62</code>, <code>Hf(acac)4 9.04</code>).
            </p>
          </div>
          <div class="card">
            <h4 style="margin:0 0 6px 0; color:#b45309;">주의사항: 0.0 덮어쓰기 및 데이터 누설 방지</h4>
            <p style="font-size:0.83rem; color:#475569; margin:0; line-height:1.6;">
              평가 대상 프리커서의 정답은 제외하고 실제 상용 시약 29종만 반영해야 데이터 누설(Data Leakage) 없는 공정 평가가 가능하며, 시약 조각을 <code>0.0</code>으로 덮어쓰면 Eq. (12) <code>&lt; 0</code> 필터에서 빠져 Cp/아미디네이트 3종이 악화되므로(<code>audit_reproduce_all.py</code> Check 5) 반드시 음수 영역 하한(<code>max(existing, -0.47)</code>)으로만 보정합니다.
            </p>
          </div>
          <div class="card">
            <h4 style="margin:0 0 6px 0; color:#0f766e;">최종 결합: Pro-ligand + SAPIS</h4>
            <p style="font-size:0.83rem; color:#475569; margin:0; line-height:1.6;">
              <code>compute_ligand_sa()</code>가 금속 중심을 떼고 리간드 $HL_i$에만 <code>min(BR-SAScore, SAScore)</code>를 적용해 케토-에놀 오판(<code>Hacac 5.20 &rarr; 2.17</code>)까지 복구하며, 배위&middot;열안정성은 <code>SevenNet-0</code> 물리 지표(<code>%V_bur</code>, <code>BDE_min</code>, <code>E_coh</code>)로 보완합니다.
            </p>
          </div>
        </div>

        <div class="callout-info" style="margin:12px 0; background:#f8fafc; border-color:#cbd5e1;">
          <strong style="font-size:0.94rem; color:#0f172a;">*Pro-ligand 분리 (Pro-ligand Decoupling)</strong>
          <div class="brt-split-grid">
            <ul style="margin:0 0 0 18px; line-height:1.68; font-size:0.84rem; color:#334155;">
              <li><strong>Step 1 &mdash; 배위 중심(금속) 제거:</strong> 분자 그래프에서 중심 금속 원자(적색 강조, 예: <code>Ti</code>, <code>Hf</code>)와 이에 직접 연결된 금속&ndash;리간드 결합(주황 강조)을 절단해 유기 리간드 골격만 분리합니다(<code>compute_ligand_sa()</code>).</li>
              <li><strong>Step 2 &mdash; 수소 포화 및 중성 프로리간드($HL_i$) 복원:</strong> 끊어진 배위 자리의 라디칼을 제거하고 암묵적 수소(<code>H</code>)로 원자가를 채워, 실제로 구매·합성 가능한 중성 프로리간드 분자(예: <code>-NMe2</code> &rarr; <code>Me2NH</code>, <code>-Cp</code> &rarr; <code>HCp</code>)로 변환합니다. (단순 무기종 <code>[H]</code>, <code>Cl</code>, <code>O</code>는 제외)</li>
              <li><strong>Step 3 &mdash; 율속 리간드 기준 점수 산출:</strong> 각 프로리간드 $HL_i$에 <code>min(BR-SAScore, SAScore)</code>를 적용해 케토-에놀 등 사전 편향을 상호 보완한 뒤, 가장 합성하기 어려운 <strong>최댓값($\max_i$)</strong>을 <code>ligand_sa</code>로 채택합니다(예: $\max(1.02, 3.31) = \mathbf{3.31}$). 이후 금속&ndash;리간드 결합·열안정성은 <code>SevenNet-0</code> 물리 지표가 전담합니다.</li>
            </ul>
            <div style="background:#fff; border:1px solid #cbd5e1; border-radius:8px; padding:10px; text-align:center;">
              __PROLIGAND_DIAGRAM_HTML__
              <div style="font-size:0.73rem; color:#475569; margin-top:6px; line-height:1.45;">
                <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#e63946; margin-right:3px;"></span><strong>제거 대상</strong> 중심 금속 &nbsp;
                <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#f59e0b; margin-right:3px;"></span><strong>절단 결합</strong> M&ndash;L &nbsp;
                <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#0f766e; margin-right:3px;"></span><strong>율속 리간드</strong> max = 3.31
              </div>
            </div>
          </div>
        </div>

        <div class="card" style="background:#f8fafc; border-color:#cbd5e1; margin-top:12px;">
          <h4 style="margin:0 0 8px 0; color:#0f172a;">
            바닐라 BR-SAScore vs 반도체 DB 확장 외삽 정합성 원리 및 명시적 수치 계산 (<code>CpTi(NMe2)3</code> 실측)
          </h4>
          <p style="font-size:0.84rem; color:#334155; line-height:1.68; margin:0 0 10px 0;">
            바닐라 <code>BR-SAScore</code>에 반도체 상용 시약 29종과 프로리간드 분리 체계를 결합할 때 미지의 신규 프리커서에 대한 순위 및 판정 일관성, 즉 <strong>외삽 정합성</strong>이 향상되는 핵심 원리는 <strong>Eq. (12)의 음수 조각 집합 $\text{Score}_k &lt; 0$의 분모 보존과 해시 사전의 국소 불변성</strong>에 있습니다. 신경망 모델은 새로운 시약 29종을 추가 학습하면 전체 결정 경계가 이동하여 기존 유기·무기 분자의 점수까지 흔들리지만, <code>BR-SAScore</code>는 미등재 금속-리간드 조각을 <code>-6.0</code>으로 유지한 채 <strong>실제 상용 시약이 공유하는 리간드 조각만 음수 영역 내부 하한인 <code>max(existing, -0.47)</code>로 안전하게 상향 보정</strong>합니다. 따라서 음수 조각 개수인 분모 $n_{\text{neg}}$가 그대로 유지되면서 점수가 단조롭게 개선되며, 여기에 프로리간드 분리를 결합하면 금속 중심의 미등재 <code>-6.0</code> 붕괴까지 근본적으로 해소되어 신규 원소 프리커서로도 일관되게 외삽됩니다.
          </p>
          <div style="background:#fff; border:1px solid #cbd5e1; border-radius:8px; padding:12px 14px; font-size:0.8rem; color:#334155; line-height:1.65;">
            <strong style="color:#0369a1;">실측 대입 비교: <code>CpTi(NMe2)3</code> ($\text{C}_{11}\text{H}_{23}\text{N}_3\text{Ti}$, $n_{\text{At}}=15,\; |\text{FP}|=21 \implies s_3=0.000,\; \text{comp}=0.2045$)</strong>
            <div class="brt-split-grid">
              <ul style="margin:0 0 0 18px; padding:0;">
                <li><strong>① 바닐라 BR-SAScore (Default):</strong> 고유 반경-2 조각 7개 중 금속 인접 5개(적색 영역)는 <code>-6.0</code>, Cp 고리 2개(청색 영역)는 <code>-1.016, -1.324</code> ($n_{\text{neg}}=7$).<br>
                  $\text{BR-frag} = \frac{5(-6.0) - 1.016 - 1.324}{7} = -4.620 \implies \text{raw} = -4.825 \implies 10 - 9\frac{-4.825+7.0}{7.0} = \mathbf{7.203}$</li>
                <li><strong>② 반도체 시약 DB 확장 (<code>Strict Unleaked</code>, <code>max(existing, -0.47)</code>):</strong> 시약 사전의 Cp 조각 2개(청색)가 <code>-0.47</code>로 보정되나 여전히 <code>&lt; 0</code>이므로 분모 $n_{\text{neg}}=7$이 보존됨.<br>
                  $\text{BR-frag} = \frac{5(-6.0) + 2(-0.47)}{7} = -4.420 \implies \text{raw} = -4.625 \implies 10 - 9\frac{-4.625+7.0}{7.0} = \mathbf{6.946}$ (<code>-0.257</code> 개선, 기존 27종 악화 0건)</li>
                <li><strong>③ (오류 대조) <code>0.0</code> 덮어쓰기 함정:</strong> Cp 2개 조각을 <code>0.0</code>으로 덮어쓰면 Eq. (12) <code>&lt; 0</code> 필터에서 탈락해 분모가 $n_{\text{neg}}=5$로 축소됨.<br>
                  $\text{BR-frag} = \frac{5(-6.0)}{5} = -6.000 \implies \text{raw} = -6.205 \implies \mathbf{8.977}$ (<code>+1.774</code> 역전 악화)</li>
                <li><strong>④ 완전 외삽 (Pro-ligand 분리 + <code>SAPIS</code>):</strong> <code>Ti</code> 중심을 분리한 <code>Me2NH</code>(<code>BR=1.02</code>)와 <code>HCp</code>(<code>BR=3.33</code>)는 모두 사전 내에 등재되어 <code>ligand_sa = 3.31</code>($\le 3.5$ 통과)을 산출하고, <code>SevenNet-0</code> 물리 지표가 금속-리간드 안정성을 결합해 미지 프리커서에서도 높은 외삽 정합성을 달성함.</li>
              </ul>
              <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; padding:10px; text-align:center;">
                __CPTI_2D_SVG__
                <div style="font-size:0.74rem; color:#475569; margin-top:6px; line-height:1.5;">
                  <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#e63946; margin-right:3px;"></span><strong>금속 인접 미등재</strong> <code>-6.0</code> (5개 조각) &nbsp;
                  <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#0284c7; margin-right:3px;"></span><strong>Cp 시약 확장</strong> <code>-0.47</code> (2개 조각)
                  <div style="font-size:0.71rem; color:#64748b; margin-top:2px;">* 비강조 말단 메틸기(<code>-CH3</code>)는 <code>&ge; 0</code>으로 Eq. (12) 분모 제외</div>
                </div>
              </div>
            </div>
            <div style="margin-top:10px; padding:9px 12px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; font-size:0.78rem; color:#334155; line-height:1.6;">
              <strong style="color:#0f172a;">BR-SAScore 수치 해석 기준 (1~10 척도, 낮을수록 합성 용이):</strong>
              <span style="display:inline-block; padding:1px 7px; border-radius:999px; background:#dcfce7; color:#166534; font-weight:700; margin:0 4px;">1.0 &le; 점수 &le; 3.5 : 합성 용이 (Easy / Pass)</span>
              <span style="display:inline-block; padding:1px 7px; border-radius:999px; background:#fef3c7; color:#92400e; font-weight:700; margin:0 4px;">3.5 &lt; 점수 &le; 6.0 : 중간 난이도 (Moderate)</span>
              <span style="display:inline-block; padding:1px 7px; border-radius:999px; background:#fee2e2; color:#991b1b; font-weight:700; margin:0 4px;">6.0 &lt; 점수 &le; 10.0 : 합성 난해 (Hard / Fail)</span>
              <div style="margin-top:4px; color:#475569;">
                * <strong>판정 기준 적용:</strong> <code>SAScore (2009)</code>와 동일하게 <strong>$3.5$ 이하</strong>를 표준 합성 가능 임계치로 사용합니다. 위 예시의 ①번 <code>7.203</code>점과 ②번 <code>6.946</code>점은 시약 사전 확장으로 일부 점수가 개선되더라도 금속-리간드 미등재 패널티 <code>-6.0</code>으로 인해 여전히 합성 난해 구간인 <code>6.0</code> 초과 영역에 머무르지만, ④번 프로리간드 분리 후 산출된 <code>ligand_sa = 3.31</code>은 <strong><code>&le; 3.5</code> 기준을 통과</strong>하여 실제 합성 가능성을 정확히 반영합니다.
              </div>
            </div>
          </div>
        </div>
      </section>
    """

    replacements = {
        "__SA_2009_DOI__": SA_2009_DOI,
        "__BR_PAPER_DOI__": BR_PAPER_DOI,
        "__VANILLA_COMPLEX_ROWS__": vanilla_complex_rows_html,
        "__VANILLA_LIGAND_ROWS__": vanilla_ligand_rows_html,
        "__MORGAN_2D_SVG__": morgan_2d_svg,
        "__CPTI_2D_SVG__": cpti_2d_svg,
        "__PROLIGAND_DIAGRAM_HTML__": proligand_diagram_html,
        "__DIST_SVG__": dist_svg,
        "__CARDS_HTML__": cards_html,
        "__B_COUNT__": f"{b_stat['count']:,}",
        "__B_MIN__": str(b_stat["min"]),
        "__B_MAX__": str(b_stat["max"]),
        "__B_MEAN__": str(b_stat["mean"]),
        "__B_NEG_PCT__": str(b_stat["neg_pct"]),
        "__R_COUNT__": f"{r_stat['count']:,}",
        "__R_MIN__": str(r_stat["min"]),
        "__R_MAX__": str(r_stat["max"]),
        "__R_MEAN__": str(r_stat["mean"]),
        "__R_NEG_PCT__": str(r_stat["neg_pct"]),
        "__BR_COUNT__": f"{br_stat['count']:,}",
        "__BR_MIN__": str(br_stat["min"]),
        "__BR_MAX__": str(br_stat["max"]),
        "__BR_MEAN__": str(br_stat["mean"]),
        "__BR_NEG_PCT__": str(br_stat["neg_pct"]),
        "__BR_NEG_CNT__": f"{br_stat['neg_count']:,}",
        "__BR_NONNEG__": str(br_stat["non_neg_count"]),
        "__OV_BOTH__": f"{ov['both']:,}",
        "__OV_B_WINS__": f"{ov['b_wins_in_overlap']:,}",
        "__OV_R_WINS__": f"{ov['r_wins_in_overlap']:,}",
        "__OV_R_ONLY__": f"{ov['r_only']:,}",
        "__OV_R_PCT__": f"{100 * ov['r_only'] / ov['union']:.1f}",
    }
    html_out = html_template
    for key, val in replacements.items():
        html_out = html_out.replace(key, val)
    return html_out
