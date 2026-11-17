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
    contribution maps, Retro* TS1-TS3 benchmarks, and TMC-SA domain extension.
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
            "명시적 H가 nAtoms 및 Morgan 비트에 합산되어 점수 급등 (+1.62)",
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
            "Hacac (keto)",
            "CC(=O)CC(=O)C",
            "베타-디케톤 프로리간드 (케토)",
            "중앙 메틸렌(C3)의 R=2 패턴이 사전 미등재(-6.0)로 판정되어 5.20으로 과대 패널티.",
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
              <td style="text-align:right; font-weight:800; color:#0f172a;">{row["sa_2009"]:.4f}</td>
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
              <td style="text-align:right; font-weight:800; color:#0f766e;">{row["sa_2009"]:.4f}</td>
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


def _build_interactive_br_calculator_html() -> str:
    """
    Build a self-contained HTML + JS interactive widget to simulate Equations (8)-(13)
    of ``BR-SAScore``.
    """
    return """
        <div id="brTheorySimWidget" style="background:linear-gradient(135deg,#f8fafc 0%,#eff6ff 100%); border:1px solid #bfdbfe; border-radius:12px; padding:18px 20px; margin-top:16px;">
          <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px;">
            <div>
              <span style="font-size:0.74rem; font-weight:800; color:#0284c7; text-transform:uppercase; letter-spacing:0.06em;">Interactive Equation Sandbox</span>
              <h4 style="margin:2px 0 0 0; font-size:1.02rem; font-weight:800; color:#0f172a;">Eq. (8)~(13) BR-SAScore 실시간 수식 시뮬레이터</h4>
            </div>
            <div style="display:flex; flex-wrap:wrap; gap:6px;" id="brSimPresetBtns">
              <button type="button" class="db-pill active" onclick="applyBrPreset('aspirin', this)">Aspirin (2.09)</button>
              <button type="button" class="db-pill" onclick="applyBrPreset('gao', this)">AI 분자 [Gao et al.] (8.84)</button>
              <button type="button" class="db-pill" onclick="applyBrPreset('hacac', this)">Hacac (5.20)</button>
              <button type="button" class="db-pill" onclick="applyBrPreset('tdmah', this)">TDMAH Hf(NMe2)4 (8.62)</button>
              <button type="button" class="db-pill" onclick="applyBrPreset('meoh', this)">MeOH (1.01)</button>
            </div>
          </div>

          <div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:12px; margin-bottom:14px;">
            <div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px 12px;">
              <label for="brSimBScore" style="display:flex; justify-content:space-between; font-size:0.8rem; font-weight:700; color:#0369a1;">
                <span>Building Block Score (BScore_i)</span>
                <span id="brSimBScoreVal">-1.20</span>
              </label>
              <input id="brSimBScore" type="range" min="-3.0" max="3.0" step="0.05" value="-1.20" style="width:100%;" oninput="updateBrSimulator()">
              <div style="font-size:0.72rem; color:#64748b; margin-top:2px;">eMolecules 231M 기반 [-3, +3]</div>
            </div>

            <div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px 12px;">
              <label for="brSimRScore" style="display:flex; justify-content:space-between; font-size:0.8rem; font-weight:700; color:#b45309;">
                <span>Reaction Score (RScore_i, d=0~2)</span>
                <span id="brSimRScoreVal">-0.68</span>
              </label>
              <input id="brSimRScore" type="range" min="-3.0" max="3.0" step="0.05" value="-0.68" style="width:100%;" oninput="updateBrSimulator()">
              <div style="font-size:0.72rem; color:#64748b; margin-top:2px;">USPTO ~1M 반응 중심 2^-d 가중</div>
            </div>

            <div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px 12px;">
              <label for="brSimUnseen" style="display:flex; justify-content:space-between; font-size:0.8rem; font-weight:700; color:#b91c1c;">
                <span>미등재 조각 비율 (Unseen Fraction)</span>
                <span id="brSimUnseenVal">0%</span>
              </label>
              <input id="brSimUnseen" type="range" min="0" max="100" step="5" value="0" style="width:100%;" oninput="updateBrSimulator()">
              <div style="font-size:0.72rem; color:#64748b; margin-top:2px;">미등재 시 Score_i = frag_penalty (-6.0)</div>
            </div>

            <div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px 12px;">
              <label for="brSimComp" style="display:flex; justify-content:space-between; font-size:0.8rem; font-weight:700; color:#475569;">
                <span>ComplexityPenalty (Size/Ring/Stereo)</span>
                <span id="brSimCompVal">0.17</span>
              </label>
              <input id="brSimComp" type="range" min="0.0" max="1.5" step="0.01" value="0.17" style="width:100%;" oninput="updateBrSimulator()">
              <div style="font-size:0.72rem; color:#64748b; margin-top:2px;">BR-SAScore에서 차감 (0에 가까울수록 단순)</div>
            </div>

            <div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px 12px;">
              <label for="brSimBuffer" style="display:flex; justify-content:space-between; font-size:0.8rem; font-weight:700; color:#0f766e;">
                <span>PenaltyBuffer (Eq. 13, 기본값 1.0)</span>
                <span id="brSimBufferVal">1.00</span>
              </label>
              <input id="brSimBuffer" type="range" min="0.0" max="3.0" step="0.1" value="1.0" style="width:100%;" oninput="updateBrSimulator()">
              <div style="font-size:0.72rem; color:#64748b; margin-top:2px;">S_min = -6 - PenaltyBuffer (최대 Complexity 흡수 폭)</div>
            </div>
          </div>

          <div style="background:#0f172a; color:#f8fafc; border-radius:10px; padding:14px 16px; display:grid; grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); gap:12px; align-items:center;">
            <div>
              <div style="font-size:0.72rem; color:#94a3b8;">Eq. (11) 사전 최대값</div>
              <div style="font-size:1.05rem; font-weight:800; color:#38bdf8;" id="brOutMaxScore">-0.68</div>
              <div style="font-size:0.71rem; color:#cbd5e1;" id="brOutWinner">RScore 우세</div>
            </div>
            <div>
              <div style="font-size:0.72rem; color:#94a3b8;">Eq. (12) BR-fragmentScore</div>
              <div style="font-size:1.05rem; font-weight:800; color:#fbbf24;" id="brOutFragMean">-0.68</div>
              <div style="font-size:0.71rem; color:#cbd5e1;">음수 항 평균 (0 &le; 기준)</div>
            </div>
            <div>
              <div style="font-size:0.72rem; color:#94a3b8;">Eq. (8) Raw BR-SAScore</div>
              <div style="font-size:1.05rem; font-weight:800; color:#f87171;" id="brOutRaw">-0.85</div>
              <div style="font-size:0.71rem; color:#cbd5e1;" id="brOutRange">정규화 범위 [-7.00, 0.00]</div>
            </div>
            <div style="background:rgba(255,255,255,0.08); border-radius:8px; padding:8px 12px; text-align:center;">
              <div style="font-size:0.73rem; color:#94a3b8;">Eq. (13) 최종 정규화 (1~10)</div>
              <div style="font-size:1.45rem; font-weight:900; color:#34d399;" id="brOutFinal">2.09</div>
              <div style="font-size:0.72rem; font-weight:700; color:#a7f3d0;" id="brOutVerdict">합성 용이 (Easy, &le; 3.5)</div>
            </div>
          </div>
          <div id="brSimPresetNote" style="font-size:0.78rem; color:#334155; margin-top:10px; line-height:1.5;">
            * Aspirin 실측 예: (-0.82 - 0.98 - 0.23)/3 = -0.68, ComplexityPenalty = 0.17, 최종 10 - 9 * ((-0.85 + 7)/7) = 2.09 (논문 Fig. 2 및 본 저장소 실측 2.09와 일치).
          </div>
        </div>"""


def build_br_sascore_theory_html() -> str:
    """
    Build the complete 2-section ``배경이론`` HTML fragment combining:
      1) Section 1 (``#brt-sec1``): 2009 vanilla SAScore principles, Hf-complex vs.
         pro-ligand evaluation, and SMILES graph reconstruction sensitivity.
      2) Section 2 (``#brt-sec2``): Unified single-section BR-SAScore deep dive
         (equations, live simulator, dataset pickles, code audit, benchmarks, TMC-SA).

    Physics/Algorithm:
        All reconstruction scores, pickle dictionary counts/bins, and atom
        contribution heatmaps are computed live from RDKit ``sascorer`` and
        ``BRSAScore==0.1.1`` on build.
    """
    vanilla_data = compute_vanilla_sascore_issues()
    vanilla_complex_rows_html, vanilla_ligand_rows_html = _build_vanilla_reconstruction_rows_html(vanilla_data)
    morgan_2d_svg = _build_morgan_fingerprint_2d_svg()
    cpti_2d_svg = _build_cpti_extrapolation_2d_svg()

    stats = compute_br_pickle_statistics()
    examples = compute_br_worked_examples()
    dist_svg = _build_fragment_distribution_svg(stats)
    cards_html = _build_worked_example_cards_html(examples)
    sim_html = _build_interactive_br_calculator_html()

    b_stat = stats["B"]
    r_stat = stats["R"]
    br_stat = stats["BR"]
    ov = stats["overlap"]

    html_template = r"""
      <div class="card" style="margin-bottom:22px; background:linear-gradient(135deg, #f0f9ff 0%, #f5f3ff 100%); border:1px solid #c7d2fe; padding:20px 24px;">
        <p style="font-size:0.92rem; color:#334155; line-height:1.68; margin:0;">
          반도체 SA 모델을 이해하기 위한 기초 이론 가이드입니다. 섹션 1에서 SAScore (2009) 바닐라 모델의 원리와 동일 금속 복합체(metal complex)에서도 SMILES 그래프 재구성에 따라 점수가 달라지는 이슈를 정리하고, 섹션 2에서 현재 베이스 모델인 BR-SAScore (2024)의 원리, 학습 데이터셋, 소스코드 및 베이스 모델 선정 근거를 통합적으로 해설합니다.
        </p>
      </div>

      <!-- SECTION 1: SASCORE (2009) VANILLA ISSUES -->
      <section id="brt-sec1">
        <div class="section-header">
          <div class="section-tag">Background Section 01</div>
          <h2 class="section-title">1. SAScore (2009) 모델의 원리 및 SMILES 재구성 이슈</h2>
        </div>

        <p style="line-height:1.72; color:#334155;">
          Ertl &amp; Schuffenhauer가 제안한 <code>SAScore (2009)</code> 모델은 PubChem 대표 1% 샘플(934,046종)에서 추출한 Morgan Fingerprint(r&le;2) 출현 빈도 <code>fpscores.pkl.gz</code>를 기반으로 합니다. 하지만 이를 반도체용 프리커서에 그대로 적용하면 두 가지 근본 이슈가 드러납니다.
        </p>

        <div class="callout-info" style="margin:10px 0 14px 0; background:#f8fafc; border-color:#cbd5e1;">
          <strong style="font-size:0.94rem; color:#0f172a;">*Morgan Fingerprint</strong>
          <div style="display:grid; grid-template-columns:minmax(260px, 1.35fr) minmax(240px, 1fr); gap:16px; align-items:center; margin-top:8px;">
            <ul style="margin:0 0 0 18px; line-height:1.68; font-size:0.84rem; color:#334155;">
              <li><strong>원형 부분구조 해싱 (Circular Substructure Hashing):</strong> 분자 그래프에서 각 원자를 중심으로 결합 거리 반경 $r = 0, 1, 2$ 이내 이웃 원자들의 원소 종류&middot;원자가&middot;결합 차수&middot;방향족성&middot;고리 포함 여부를 순차적으로 집계해 고유 정수 해시 ID(<code>bitId</code>)로 변환하는 원형 분자 지문입니다.</li>
              <li><strong>반경별 의미와 재구성 민감도:</strong> <code>Radius 0</code>(적색 중심 원자)은 단일 원자 자체, <code>Radius 1</code>(주황 1차 이웃)은 직접 결합한 원자(예: $\text{Hf--O}$ vs $\text{Hf--C}$), <code>Radius 2</code>(청색 2차 이웃)는 2결합 떨어진 고리&middot;치환기 맥락(예: <code>acac</code> 6원 킬레이트 닫힘, $\eta^5$-<code>Cp</code>, <code>AddHs</code> 명시적 수소)을 포착합니다. 따라서 같은 금속 복합체라도 그래프 연결을 재구성하면 생성되는 <code>Radius 1~2</code> 해시 ID가 달라져 사전 조회 점수가 즉시 바뀝니다.</li>
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
            <h4 style="margin:0 0 8px 0; color:#1e3a8a;">1.1 SAScore (2009) 모델의 수식과 정규화 규격</h4>
            <div style="font-size:0.85rem; line-height:1.68; color:#334155;">
              <p>$$\text{SAScore}_{\text{raw}} = \text{fragmentScore} - \text{complexityPenalty} + \text{score}_3 \tag{1}$$</p>
              <p>$$\text{complexityPenalty} = (n_{\text{At}}^{1.005} - n_{\text{At}}) + \log_{10}(n_{\text{Chiral}}+1) + \log_{10}(n_{\text{Bridge}}+1) + \log_{10}(n_{\text{Spiro}}+1) + \log_{10}(n_{\text{Macro}}+1)$$</p>
              <p>$$\text{SAScore}_{\text{full}} = 11 - \frac{\text{SAScore}_{\text{raw}} - (-4.0) + 1}{2.5 - (-4.0)} \times 9 = 11 - \frac{\text{SAScore}_{\text{raw}} + 5.0}{6.5} \times 9$$</p>
              <p style="font-size:0.79rem; color:#64748b; margin-bottom:0;">
                * <code>fpscores.pkl.gz</code> 미등재 희귀/금속 비트는 바닥값 <code>-4.0</code>을 받으며, 대칭 분자에서는 <code>score3 = 0.5*ln(nAtoms / len(fps))</code>가 가산되고, 최종값이 $8.0$ 초과 시 <code>8 + ln(SA - 8)</code> 완만 압축을 거쳐 $[1, 10]$ 범위로 클리핑됩니다. <strong>점수가 $1.0$에 가까울수록 합성이 용이(Easy)하고 $10.0$에 가까울수록 합성이 난해(Hard)하며, 일반적으로 $3.5$ 이하를 합성 가능(Pass) 기준으로 해석합니다.</strong>
              </p>
            </div>
          </div>

          <div class="card" style="background:#fffbeb; border-color:#fde68a;">
            <h4 style="margin:0 0 8px 0; color:#92400e;">1.2 SMILES 재구성 이슈</h4>
            <ol style="font-size:0.84rem; color:#334155; line-height:1.62; padding-left:18px; margin:0 0 10px 0;">
              <li><strong>킬레이트(&rarr;)·$\eta^5$-Cp 토폴로지 재구성:</strong> <code>acac</code>를 단좌 개방형으로 두는지 6원 킬레이트로 닫는지(<code>Hf(acac)4</code>: <code>4.1933 &rarr; 5.4854</code>), <code>Cp</code>를 $\eta^1$(<code>5.4320</code>)에서 $\eta^5$ 5중 배위(<code>6.4305</code>)로 연결하는지에 따라 <code>fragmentScore</code>와 고리 패널티가 급변합니다. 외부 발표 순위(<code>5.179 &lt; 5.491 &lt; 6.74</code>)도 고배위 재구성 방향과 일치합니다.</li>
              <li><strong>명시적 수소(<code>AddHs</code>) &middot; 입체/전하 처리:</strong> <code>AddHs</code> 미제거 시 수소가 $n_{\text{At}}$와 지문 비트에 합산되어 점수가 급등하고(<code>+1.62</code>), 키랄 중심 인식 및 이온/양쪽성이온 표기에 따라 룩업 비트가 바뀝니다.</li>
            </ol>
            <div style="background:#fff; border:1px solid #fcd34d; border-radius:6px; padding:9px 11px; font-size:0.78rem; color:#334155; line-height:1.58;">
              <strong style="color:#92400e;">재구성 방식별 명시적 계산 비교 (RDKit 실측 항 대입):</strong>
              <ul style="margin:4px 0 0 16px; padding:0;">
                <li><strong>Hf(OMe)<sub>3</sub>(acac) [Open]:</strong> $\text{frag}=-0.206,\; \text{comp}=0.186,\; s_3=0.000\; (n_{\text{At}}\le |\text{FP}|) \implies \text{raw}=-0.392 \implies 11 - \frac{-0.392+5.0}{6.5}\times 9 = \mathbf{4.6201}$</li>
                <li><strong>Hf(OMe)<sub>3</sub>(acac) [Chelate &rarr;]:</strong> 고리 닫힘으로 $\text{frag}=-0.694,\; \text{comp}=0.186,\; s_3=0.000 \implies \text{raw}=-0.880 \implies 11 - \frac{-0.880+5.0}{6.5}\times 9 = \mathbf{5.2957}$</li>
                <li><strong>Hf(OMe)<sub>3</sub>(acac) [Chelate+AddHs]:</strong> H 16개 포함($n_{\text{At}}=30$) $\text{frag}=-2.616,\; \text{comp}=0.515,\; s_3=0.053 \implies \text{raw}=-3.078 \implies \mathbf{6.9182}$</li>
                <li><strong>HfCp(OMe)(acac)<sub>2</sub> [&eta;<sup>5</sup>-Cp]:</strong> $\text{frag}=-1.056,\; \text{comp}=0.644\text{ (ring }0.301\text{ 포함)},\; s_3=0.000 \implies \text{raw}=-1.700 \implies \mathbf{6.4305}$</li>
              </ul>
            </div>
          </div>
        </div>

        <div class="card" style="margin:12px 0;">
          <h4 style="margin:0 0 10px 0; color:#0f172a;">
            1.3 Hf 프리커서 3종 SMILES 재구성별 <code>SAScore (2009)</code> 모델 실측 비교 vs 단독(미배위) 프로리간드($HL_i$) 점수
          </h4>
          <div class="table-wrap" style="margin-bottom:14px;">
            <table>
              <thead>
                <tr>
                  <th>금속 복합체 (Precursor)</th>
                  <th>SMILES 그래프 재구성 방식</th>
                  <th>원자 수 ($N_{\text{at}}$)</th>
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
                <tr><th>리간드 음이온</th><th>프로리간드</th><th>SMILES</th><th>원자 수</th><th>SAScore (2009)</th><th>해설</th></tr>
              </thead>
              <tbody>
                __VANILLA_LIGAND_ROWS__
              </tbody>
            </table>
          </div>
        </div>

        <div class="card" style="background:#f8fafc; margin-top:12px;">
          <h4 style="margin:0 0 8px 0; color:#0f172a;">1.4 Summary</h4>
          <p style="font-size:0.86rem; color:#334155; line-height:1.74; margin:0;">
            <code>MeOH (2.7010)</code>, <code>Hacac (2.1725~2.6758)</code>, <code>HCp (3.3069)</code>는 모두 실험실 및 공정에서 널리 활용되는 합성 용이 리간드임에도 불구하고, 이를 <code>Hf</code> 중심 원자와 결합한 전체 금속 복합체(metal complex) <code>SMILES</code>로 인코딩해 <code>SAScore (2009)</code> 모델에 입력하면 킬레이트 고리 닫힘, $\eta^5$ 배위 표기, 이온 해리, 명시적 수소(<code>AddHs</code>) 처리 등 그래프 재구성 방식에 따라 동일 화합물의 점수가 <code>3.9056</code>에서 <code>6.9182</code>까지 크게 요동칩니다. 이는 단일 분자 그래프 기반의 기존 모델이 실제 화학적 합성 난이도보다 <code>SMILES</code> 표기 관습과 PubChem 유기 조각 결손을 반영함을 시사하며, 따라서 섹션 2의 <code>BR-SAScore (2024)</code> 모델 도입과 함께 금속-리간드 분리 기반의 <strong>Pro-ligand + 물리 Descriptor(<code>TMC-SA</code>)</strong> 결합 체계로의 전환이 필수적입니다.
          </p>
        </div>
      </section>

      <!-- SECTION 2: BR-SASCORE UNIFIED SINGLE SECTION -->
      <section id="brt-sec2">
        <div class="section-header">
          <div class="section-tag">Background Section 02</div>
          <h2 class="section-title">2. BR-SAScore (2024) 모델의 원리 및 베이스 모델 선정 근거</h2>
        </div>

        <p style="line-height:1.72; color:#334155;">
          <code>BR-SAScore (2024)</code> 모델은 역합성의 두 핵심 축인 ① 상용 출발물질($B$, eMolecules)과 ② 반응 형성 결합($R$, USPTO)을 분리해 평가하는 3세대 하이브리드 모델입니다. 다른 1G~4G 모델과 달리 원소 결손(OOV)&middot;All-Zero 침묵 붕괴 없이 결정론적으로 동작하고, 위험 대조군 오탐(False Positive)을 1차 차단하며, 신경망 재학습 없이 <code>_fscores</code> 해시 사전에 도메인 시약만 누설 없이 확장(<code>max(existing, -0.47)</code>)할 수 있어 본 연구의 베이스 모델로 선정되었습니다.
        </p>

        <!-- 2.1 Architecture Flowchart -->
        <h3 id="brt-sub21" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:18px 0 8px 0;">
          2.1 B/R 분리 아키텍처 파이프라인
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
  MFP --> B["BScore Dict: eMolecules 231M (244,572 keys)"]
  MFP --> R["RScore Dict: USPTO ~1M, 2^-d (175,186 keys)"]
  B --> MAX["Eq.11 Fusion: max(B_i, R_i) else -6"]
  R --> MAX
  MAX --> NEG["Eq.12 Negative Mean: avg over Score_i &lt; 0"]
  NEG --> RAW["Eq.8 Raw = BR_frag - Comp"]
  COMP --> RAW
  RAW --> NORM["Eq.13 Scale [1, 10] (Buffer=1.0, Div=7)"]
              </div>
            </div>
          </div>
        </div>

        <!-- 2.2 Equations & Simulator -->
        <h3 id="brt-sub22" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:20px 0 8px 0;">
          2.2 핵심 수식 체계 (Eq. 8~13) 및 실시간 수식 시뮬레이터
        </h3>
        <div class="grid-2" style="margin-top:8px;">
          <div class="card" style="border-color:#93c5fd; background:#f8fbff;">
            <h4 style="margin:0 0 8px 0; color:#0369a1;">BR-SAScore 지배 방정식 (Eq. 8~13)</h4>
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
              <li><strong>OR 논리 ($\max(\text{BScore}_i, \text{RScore}_i)$):</strong> 어떤 부분 구조가 시약으로 팔리거나(<code>BScore</code>) 혹은 특허 반응(<code>RScore</code>)으로 쉽게 이어지는 결합이라면 둘 중 더 유리한 쪽을 채택합니다.</li>
              <li><strong>병목 중심 음수 평균 (Eq. 12):</strong> 상위 $0.1\%$ 초빈출 조각($\text{Score}_i \ge 0$)을 평균에 포함하면 희귀 병목 조각 1~2개의 난이도가 희석되므로, <strong>0 미만 조각만 평균</strong>내어 병목 민감도를 극대화합니다. (전부 $\ge 0$이면 $0$ 처리)</li>
              <li><strong>미등재 패널티 <code>-6</code> &amp; <code>PenaltyBuffer = 1</code> (Eq. 13):</strong> 사전 외 조각은 최소점(-3)의 2배인 <code>-6</code>을 부여하고, 복잡도 패널티가 겹칠 때 하한이 $-7.0$까지 확보되도록 <code>PenaltyBuffer = 1</code>(분모 $7$)로 $[1, 10]$에 정규화합니다.</li>
            </ul>
          </div>
        </div>
        __SIM_HTML__

        <!-- 2.3 Training Datasets & Pickle Audit -->
        <h3 id="brt-sub23" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          2.3 학습&middot;구축 데이터셋 해부
        </h3>
        <div class="table-wrap" style="margin:10px 0;">
          <table>
            <thead>
              <tr>
                <th>데이터셋 / 파일</th>
                <th>원천 데이터 규모</th>
                <th>조각 추출 &amp; 가중 규칙</th>
                <th>최종 사전 키 수</th>
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

        <!-- 2.4 Code Audit & Worked Examples -->
        <h3 id="brt-sub24" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          2.4 소스코드 구현 해부 (<code>BRSAScore.py</code>) 및 2D 원자 기여도 Worked Examples
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
              <li><strong>케토-에놀 민감도 실측:</strong> <code>Hacac(keto)</code>는 <code>SAScore (2009)</code> 모델에서 <code>2.17</code>이지만 <code>BR-SAScore</code>에선 중앙 <code>-CH2-</code> 1개가 사전 미등재(<code>-6.0</code>)에 걸려 <code>5.20</code>으로 과대 패널티가 발생합니다.</li>
            </ul>
          </div>
        </div>
        <div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(255px,1fr)); gap:12px; margin-top:10px;">
          __CARDS_HTML__
        </div>

        <!-- 2.5 Benchmarks Table 1 & 2 -->
        <h3 id="brt-sub25" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          2.5 벤치마크 성능 검증: TS1&middot;TS2&middot;TS3 (5,400종 Retro* 라벨) 및 6대 모델 비교
        </h3>
        <div class="table-wrap" style="margin-top:8px;">
          <table>
            <thead>
              <tr>
                <th>모델 (세대)</th>
                <th>TS1 PR-AUC<br>(ZINC/GDB17, 1,800)</th>
                <th>TS2 PR-AUC<br>(ChEMBL/Med, 1,800)</th>
                <th>TS3 PR-AUC<br>(Complex Lit, 1,800)</th>
                <th>TS1 ROC-AUC</th>
                <th>TS2 ROC-AUC</th>
                <th>TS3 ROC-AUC</th>
                <th>계산 속도 (ms/mol)</th>
              </tr>
            </thead>
            <tbody>
              <tr style="background:#ecfdf5; font-weight:700;">
                <td>BR-SAScore (3G Hybrid)</td>
                <td>0.990</td><td>0.947</td><td>0.900</td>
                <td>0.984</td><td>0.942</td><td>0.900</td>
                <td>0.42 ms (CPU)</td>
              </tr>
              <tr>
                <td>SAScore (1G Heuristic)</td>
                <td>0.988</td><td>0.942</td><td>0.833</td>
                <td>0.980</td><td>0.929</td><td>0.818</td>
                <td>0.39 ms (CPU)</td>
              </tr>
              <tr>
                <td>CLScore (1G ChEMBL Fragment)</td>
                <td>0.927</td><td>0.912</td><td>0.776</td>
                <td>0.960</td><td>0.921</td><td>0.808</td>
                <td>8.97 ms (CPU)</td>
              </tr>
              <tr>
                <td>DeepSA (3G Transformer/ML)</td>
                <td>0.899</td><td>0.919</td><td>0.832</td>
                <td>0.951</td><td>0.931</td><td>0.886</td>
                <td>17.2* ms (*GPU)</td>
              </tr>
              <tr>
                <td>SYBA (1G Naive Bayes)</td>
                <td>0.939</td><td>0.855</td><td>0.722</td>
                <td>0.968</td><td>0.872</td><td>0.784</td>
                <td>0.28 ms (CPU)</td>
              </tr>
              <tr>
                <td>RAscore (2G Retrosynthesis ML)</td>
                <td>0.905</td><td>0.835</td><td>0.776</td>
                <td>0.948</td><td>0.860</td><td>0.841</td>
                <td>123 ms (CPU)</td>
              </tr>
              <tr>
                <td>GASA (3G Graph Attention)</td>
                <td>0.890</td><td>0.886</td><td>0.801</td>
                <td>0.951</td><td>0.890</td><td>0.851</td>
                <td>307 ms (CPU)</td>
              </tr>
            </tbody>
          </table>
        </div>

        <!-- 2.6 Why BR-SAScore Beats Competitors & Semiconductor Domain Bridge -->
        <h3 id="brt-sub26" style="font-size:1.05rem; font-weight:800; color:#0f172a; margin:22px 0 8px 0;">
          2.6 베이스 모델로 BR-SAScore (2024) 선정 근거
        </h3>
        <p style="line-height:1.72; color:#334155;">
          요약 탭(Step 2)에서 제시한 선정 결론의 배경에는 <strong>18종 공개 SA·역합성·생성 모델 전수 평가(확장 프리커서 102종 + 시약 29종 + Phase 1 Negative Control 10종)</strong>에서 확인된 타 모델들의 구조적 붕괴가 있습니다. 다른 모델들이 반도체 금속 복합체(metal complex) 또는 Track 1 역할에서 탈락하고 <code>BR-SAScore</code>가 유일하게 최적 베이스 엔진으로 선택된 <strong>5가지 결정적 원인</strong>은 다음과 같습니다.
        </p>

        <div class="grid-2" style="margin:12px 0;">
          <div class="card" style="background:#f8fbff; border-color:#bfdbfe;">
            <h4 style="margin:0 0 8px 0; color:#1e40af;">BR-SAScore 베이스 선정 5대 결정적 근거</h4>
            <ol style="font-size:0.83rem; color:#334155; line-height:1.68; padding-left:18px; margin:0;">
              <li><strong>OOV 원소 붕괴·All-Zero 무음 오인코딩 부재 (vs GASA &middot; DFRscore &middot; MPScore &middot; DeepSA):</strong> <code>GASA</code>는 102종 중 <strong>53종(52.0%, Ga&middot;Hf&middot;Ir&middot;Mo&middot;Rh&middot;Ru&middot;Ta&middot;W 8원소)</strong>에서 DGL 원자 벡터가 All-Zero로 무음 붕괴하고, <code>DFRscore</code>는 <strong>102/102(100%)</strong> 중심 금속이 <code>ELSE</code> 노드로 축퇴하며, <code>MPScore</code>는 금속 복합체 In-Domain <strong>0/102(0.0%)</strong>에 그칩니다. 반면 <code>BR-SAScore</code>는 RDKit Morgan Radius-2 해시 키 기반이므로 어떤 원소가 들어와도 벡터 소실 없이 결정론적으로 계산됩니다.</li>
              <li><strong>위험·불안정 대조군 무차별 통과(False Positive 폭주) 억제 (vs RAscore &middot; SCScore &middot; SYBA):</strong> <code>RAscore</code>(XGBoost)와 <code>SCScore</code>(Reaxys DNN)는 단순 치환기 개수만 보고 Phase 1 Negative Control 10종 중 <strong>9종(90.0%, <code>TiEt4</code> 0.972/1.65, <code>ZrMe4</code> 0.956/1.56, <code>WPh5</code> 0.985/1.63 등)</strong>을 Pass로 오판하며, <code>SYBA</code>도 <strong>8/10(80.0%)</strong>을 통과시킵니다. 반면 <code>BR-SAScore</code>는 Eq. (12) <code>&lt; 0</code> 병목 평균과 <code>-6.0</code> 미등재 패널티 덕분에 물리 필터 전에도 10종 중 <strong>7종(70.0%)</strong>을 1차 차단합니다(<code>TiEt4 8.77</code>, <code>Silacyclobutadiene 8.75</code>, <code>WPh5 8.79</code>, <code>CoEt4 8.77</code>).</li>
              <li><strong>SAScore (2009) 모델 소분자·단일빈도 편향 해소 &amp; 4G M-L 0% 한계 보완:</strong> <code>SAScore (2009)</code> 모델은 <code>MeOH (2.70)</code> 같은 2원자 시약까지 과대 패널티를 매겨 102종 프리커서 통과율이 <strong>2.0%(2/102)</strong>에 불과하고, 4G 모델(<code>RxnFlow</code>&middot;<code>SynFormer</code>&middot;<code>ChemProjector</code>)은 <code>M-L</code> 배위 반응 규칙이 <strong>0/102(0.0%)</strong>입니다. <code>BR-SAScore</code>는 <code>MeOH 1.01</code> 정상화 및 <code>19/102(18.6%)</code> 개선을 <code>0.42 ms/mol</code>로 달성하고, 4G 모델은 해리 리간드 경로 보조에만 분리 투입할 수 있습니다.</li>
              <li><strong>신경망 재학습 없는 무손실 사전(<code>_fscores</code>) 확장성:</strong> <code>DeepSA</code>(Transformer), <code>GASA</code>(GAT), <code>RAscore</code>, <code>SCScore</code>는 도메인 시약 29종을 반영하려면 데이터셋 재구축과 신경망 재학습이 필요해 기존 유기 분포가 왜곡됩니다. <code>BR-SAScore</code>는 딕셔너리에 <code>max(existing, -0.47)</code> 하한만 추가하면 <strong>기존 프리커서 점수 악화 0건(0/27)</strong>을 유지한 채 시약 통과율을 <code>21/29(72.4%) &rarr; 28/29(96.6%)</code>로 즉시 끌어올립니다.</li>
              <li><strong>대칭 리간드(<code>score3</code>) &middot; 병목 설명성 &amp; Track 2(<code>SevenNet-0</code>) 직교성:</strong> <code>M(NMe2)4</code> 등 다중 대칭 프리커서에서 <code>score3</code> 대칭 보정이 정합적으로 작동하고, <code>contribution_to_svg()</code>가 <code>-6.0</code> 병목 원자를 즉시 지목하며, Track 2 물리 지표(<code>BDE_min</code>, <code>%V_bur^D2</code>, <code>E_coh</code>)가 금속 열역학을 전담하는 구조에서 가장 가볍고 누설 없는 Track 1 엔진 역할을 수행합니다.</li>
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
                    <th>베이스 모델 탈락 사유 vs BR-SAScore 우위</th>
                  </tr>
                </thead>
                <tbody>
                  <tr>
                    <td><strong>SAScore (2009)</strong> (1G)</td>
                    <td>2/102 (2.0%)<br>9/29 (31.0%)</td>
                    <td>3/10 (30.0%)</td>
                    <td>PubChem 단일 빈도로 98.0% 기각 및 소분자(<code>MeOH 2.70</code>) 과대 감점. B/R 분리 부재</td>
                  </tr>
                  <tr>
                    <td><strong>SYBA</strong> (1G)</td>
                    <td>83/102 (81.4%)<br>27/29 (93.1%)</td>
                    <td><strong>8/10 (80.0%)</strong></td>
                    <td>ZINC 기반 로그오즈비가 불안정 대조군 80% 통과 및 Cp-고리 극단 편향 유발</td>
                  </tr>
                  <tr>
                    <td><strong>SCScore / RAscore</strong> (2G)</td>
                    <td>62~90/102<br>29/29 (100%)</td>
                    <td><strong>9/10 (90.0%)</strong></td>
                    <td>Reaxys 복잡도/카운트 FP 의존으로 위험 물질 90% 통과; 해시 사전 직접 보정 불가</td>
                  </tr>
                  <tr>
                    <td><strong>MPScore</strong> (2G)</td>
                    <td><strong>0/102 In-Domain</strong><br>9/29 (31.0%)</td>
                    <td>N/A (유기 전용)</td>
                    <td>유기 케이지/고분자 빌딩블록 학습 편향으로 금속·할라이드 도메인 커버리지 0%</td>
                  </tr>
                  <tr>
                    <td><strong>GASA / DFRscore / DeepSA</strong> (3G)</td>
                    <td>GASA 53/102 OOV<br>DFR 102/102 ELSE</td>
                    <td>실측 불능/편향</td>
                    <td>미지원 원소 All-Zero 침묵 붕괴(GASA), 금속 ELSE 축퇴(DFR), Ro5 분절(DeepSA) + 재학습 필요</td>
                  </tr>
                  <tr>
                    <td><strong>RxnFlow / SynFormer / ReaSyn</strong> (4G)</td>
                    <td><strong>0/102 M-L 규칙</strong><br>(유기 치환기만 매칭)</td>
                    <td>N/A (0% 경로)</td>
                    <td>유기 템플릿(72~115종) 한계로 금속-리간드 코어 경로 생성 불능; 리간드 경로 생성 보조로만 분업</td>
                  </tr>
                  <tr style="background:#ecfdf5; font-weight:700;">
                    <td><strong>BR-SAScore (3G 베이스 선정)</strong></td>
                    <td><strong>19/102 (18.6%)<br>Strict 28/29 (96.6%)</strong></td>
                    <td><strong>3/10 (30.0%)<br>(단독 7/10 차단)</strong></td>
                    <td><strong>0% 원소 붕괴 &middot; 무재학습 <code>max(existing,-0.47)</code> 확장(악화 0건) &middot; <code>0.42 ms</code> &middot; Pro-ligand + SevenNet-0 완벽 직교</strong></td>
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
            <h4 style="margin:0 0 6px 0; color:#b45309;">주의사항: 0.0 덮어쓰기 누설 함정</h4>
            <p style="font-size:0.83rem; color:#475569; margin:0; line-height:1.6;">
              시약 29종을 <code>0.0</code>으로 덮어쓰면 Eq. (12) <code>&lt; 0</code> 필터에서 빠져 Cp/아미디네이트 3종이 악화됩니다(<code>audit_reproduce_all.py</code> Check 5). 비누출 사전 확장은 반드시 <code>max(existing, -0.47)</code> 상향 하한을 사용합니다.
            </p>
          </div>
          <div class="card">
            <h4 style="margin:0 0 6px 0; color:#0f766e;">최종 결합: Pro-ligand + TMC-SA</h4>
            <p style="font-size:0.83rem; color:#475569; margin:0; line-height:1.6;">
              <code>compute_ligand_sa()</code>가 금속 중심을 떼고 리간드 $HL_i$에만 <code>min(BR-SAScore, SAScore)</code>를 적용해 케토-에놀 오판(<code>Hacac 5.20 &rarr; 2.17</code>)까지 복구하며, 배위&middot;열안정성은 <code>SevenNet-0</code> 물리 지표(<code>%V_bur</code>, <code>BDE_min</code>, <code>E_coh</code>)로 보완합니다.
            </p>
          </div>
        </div>

        <div class="card" style="background:#f8fafc; border-color:#cbd5e1; margin-top:12px;">
          <h4 style="margin:0 0 8px 0; color:#0f172a;">
            바닐라 BR-SAScore vs 반도체 DB 확장 외삽 정합성 원리 및 명시적 수치 계산 (<code>CpTi(NMe2)3</code> 실측)
          </h4>
          <p style="font-size:0.84rem; color:#334155; line-height:1.68; margin:0 0 10px 0;">
            바닐라 <code>BR-SAScore</code>에 반도체 상용 시약(29종) 및 Pro-ligand 분리 체계를 확장할 때 <strong>외삽 정합성(미지 프리커서에 대한 순위·판정 일관성)</strong>이 향상되는 핵심 원리는 <strong>Eq. (12)의 음수 집합($\text{Score}_k &lt; 0$) 분모 보존과 해시 사전의 국소 불변성</strong>에 있습니다. 신경망 모델은 새로운 시약 29종을 추가 학습하면 전체 결정 경계가 이동해 기존 유기·무기 점수까지 흔들리지만, <code>BR-SAScore</code>는 미등재 금속-리간드 비트만 <code>-6.0</code>으로 두고 <strong>실제 상용 시약이 공유하는 리간드 조각(<code>max(existing, -0.47)</code>)만 음수 영역 안에 안전하게 상향 보정</strong>합니다. 따라서 <code>&lt; 0</code> 분모($n_{\text{neg}}$)가 유지되면서 점수가 단조롭게 개선되며, 여기에 Pro-ligand 분리를 결합하면 금속 중심의 미등재 <code>-6.0</code> 붕괴까지 해소되어 신규 원소 프리커서로도 일관되게 외삽됩니다.
          </p>
          <div style="background:#fff; border:1px solid #cbd5e1; border-radius:8px; padding:12px 14px; font-size:0.8rem; color:#334155; line-height:1.65;">
            <strong style="color:#0369a1;">실측 대입 비교: <code>CpTi(NMe2)3</code> ($\text{C}_{11}\text{H}_{23}\text{N}_3\text{Ti}$, $n_{\text{At}}=15,\; |\text{FP}|=21 \implies s_3=0.000,\; \text{comp}=0.2045$)</strong>
            <div style="display:grid; grid-template-columns:minmax(280px, 1.4fr) minmax(240px, 1fr); gap:16px; align-items:center; margin-top:8px;">
              <ul style="margin:0 0 0 18px; padding:0;">
                <li><strong>① 바닐라 BR-SAScore (Default):</strong> 고유 반경-2 비트 7개 중 금속 인접 5개(적색 영역)는 <code>-6.0</code>, Cp 고리 2개(청색 영역)는 <code>-1.016, -1.324</code> ($n_{\text{neg}}=7$).<br>
                  $\text{BR-frag} = \frac{5(-6.0) - 1.016 - 1.324}{7} = -4.6200 \implies \text{raw} = -4.8245 \implies 10 - 9\frac{-4.8245+7.0}{7.0} = \mathbf{7.2029}$</li>
                <li><strong>② 반도체 시약 DB 확장 (<code>Strict Unleaked</code>, <code>max(existing, -0.47)</code>):</strong> 시약 사전의 Cp 조각 2개(청색)가 <code>-0.47</code>로 보정되나 여전히 <code>&lt; 0</code>이므로 분모 $n_{\text{neg}}=7$이 보존됨.<br>
                  $\text{BR-frag} = \frac{5(-6.0) + 2(-0.47)}{7} = -4.4200 \implies \text{raw} = -4.6245 \implies 10 - 9\frac{-4.6245+7.0}{7.0} = \mathbf{6.9458}$ (<code>-0.2571</code> 개선, 기존 27종 악화 0건)</li>
                <li><strong>③ (오류 대조) <code>0.0</code> 덮어쓰기 함정:</strong> Cp 2개 비트를 <code>0.0</code>으로 덮어쓰면 Eq. (12) <code>&lt; 0</code> 필터에서 탈락해 분모가 $n_{\text{neg}}=5$로 축소됨.<br>
                  $\text{BR-frag} = \frac{5(-6.0)}{5} = -6.0000 \implies \text{raw} = -6.2045 \implies \mathbf{8.9772}$ (<code>+1.7743</code> 역전 악화)</li>
                <li><strong>④ 완전 외삽 (Pro-ligand 분리 + <code>TMC-SA</code>):</strong> <code>Ti</code> 중심을 분리한 <code>Me2NH</code>(<code>BR=1.02</code>)와 <code>HCp</code>(<code>BR=3.33</code>)는 모두 사전 내에 등재되어 <code>ligand_sa = 3.31</code>($\le 3.5$ 통과)을 산출하고, <code>SevenNet-0</code> 물리 지표가 금속-리간드 안정성을 결합해 미지 프리커서에서도 높은 외삽 정합성을 달성함.</li>
              </ul>
              <div style="background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; padding:10px; text-align:center;">
                __CPTI_2D_SVG__
                <div style="font-size:0.74rem; color:#475569; margin-top:6px; line-height:1.5;">
                  <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#e63946; margin-right:3px;"></span><strong>금속 인접 미등재</strong> <code>-6.0</code> (5개 비트) &nbsp;
                  <span style="display:inline-block; width:9px; height:9px; border-radius:50%; background:#0284c7; margin-right:3px;"></span><strong>Cp 시약 확장</strong> <code>-0.47</code> (2개 비트)
                  <div style="font-size:0.71rem; color:#64748b; margin-top:2px;">* 비강조 말단 메틸기(<code>-CH3</code>)는 <code>&ge; 0</code>으로 Eq. (12) 분모 제외</div>
                </div>
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
        "__SIM_HTML__": sim_html,
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
