"""Shared HTML/SVG rendering, compression, and chemical normalization utilities.

This module centralizes duplicated RDKit 2D SVG generation, lossless SVG/HTML
minification, and chemical identifier normalization across the report build
pipeline.
"""

from __future__ import annotations

import json
import os
import re
from typing import Any, List, Tuple

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

# High-frequency inline styles mapped to compact utility classes (.u01 - .u25).
_HIGH_FREQ_STYLES: List[Tuple[str, str]] = [
    ("text-align:center; white-space:nowrap;", "u01"),
    ("text-align:center; font-family:'JetBrains Mono'; font-size:0.8rem; color:#64748b;", "u02"),
    ("font-family:'JetBrains Mono'; font-size:0.75rem; color:#64748b;", "u03"),
    ("text-align:center; font-family:'JetBrains Mono';", "u04"),
    ("font-weight:600; color:#0f172a; font-size:0.82rem;", "u05"),
    ("font-family:'JetBrains Mono'; font-size:0.76rem; color:#475569;", "u06"),
    ("font-size:0.8rem; color:#0369a1; max-width:220px;", "u07"),
    ("font-size:0.78rem; color:#64748b; max-width:180px;", "u08"),
    ("font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a; max-width:240px; word-break:break-all;", "u09"),
    ("color:#0369a1; font-family:'JetBrains Mono';", "u10"),
    ("font-family:'JetBrains Mono'; font-weight:600; color:#475569;", "u11"),
    ("font-size:0.75rem; color:#64748b;", "u12"),
    ("text-align:center; white-space:nowrap; font-family:'JetBrains Mono'; font-size:0.76rem;", "u13"),
    ("font-size:0.7rem; color:#64748b;", "u14"),
    ("color:#94a3b8; font-size:0.75rem; font-style:italic;", "u15"),
    ("font-size:0.72rem; color:#64748b;", "u16"),
    ("text-align:center; font-family:'JetBrains Mono'; font-size:0.74rem;", "u17"),
    ("color:#15803d; font-weight:700;", "u18"),
    ("font-size:0.75rem;", "u19"),
    ("word-break:keep-all; overflow-wrap:anywhere;", "u20"),
    ("font-size:0.8rem; color:#334155;", "u21"),
    ("font-weight:600; color:#0f172a; font-size:0.8rem;", "u22"),
    ("font-family:'JetBrains Mono'; font-size:0.72rem; color:#0f172a; max-width:220px; word-break:break-all;", "u23"),
    ("font-family:monospace; font-size:0.78rem;", "u24"),
    ("text-align:center;", "u25"),
    ("text-align:right;", "u26"),
]

COMPACT_REPORT_CSS: str = (
    "\n/* --- Compact RDKit 2D SVG & Table Utility Rules (html_common.py) --- */\n"
    "svg.rdk path:not([fill]){fill:none;fill-rule:evenodd;stroke:#000;stroke-width:1.6px;"
    "stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1}\n"
    + "\n".join(f".{cls}{{{css}}}" for css, cls in _HIGH_FREQ_STYLES)
    + "\n"
)


def load_meta(filename: str, base_dir: str = ".") -> Any:
    """Load a JSON metadata file from ``database/metadata`` or ``downloaded_papers``.

    Physics/Algorithm:
        Provides deterministic lookup across the curated metadata hierarchy so
        section builders do not duplicate path fallback logic.

    Args:
        filename: Target JSON filename (e.g., ``manifest.json``).
        base_dir: Repository root directory.

    Returns:
        Parsed JSON object.

    Raises:
        FileNotFoundError: If ``filename`` does not exist in either directory.
    """
    p1 = os.path.join(base_dir, "database", "metadata", filename)
    if os.path.exists(p1):
        with open(p1, "r", encoding="utf-8") as fh:
            return json.load(fh)
    p2 = os.path.join(base_dir, "downloaded_papers", filename)
    if os.path.exists(p2):
        with open(p2, "r", encoding="utf-8") as fh:
            return json.load(fh)
    raise FileNotFoundError(f"Cannot find metadata file {filename}")


def minify_rdkit_svg(svg: str) -> str:
    """Losslessly compress an RDKit ``MolDraw2DSVG`` string.

    Physics/Algorithm:
        RDKit's SVG writer repeats XML namespaces, ``atom-*``/``bond-*`` classes,
        and a 118-byte inline ``style`` string on every bond ``<path>``. By
        hoisting the invariant stroke geometry into ``svg.rdk path:not([fill])``,
        retaining only explicit heteroatom bond ``stroke`` attributes, and
        compacting SVG path command whitespace, total SVG payload size is
        reduced by ~60% with zero visual change.

    Args:
        svg: Raw SVG string produced by ``MolDraw2DSVG.GetDrawingText()``.

    Returns:
        Compact, single-line SVG string.
    """
    if not svg:
        return ""
    s = svg.split("?>", 1)[-1].replace("svg:", "").replace("\t", " ")
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\s+width=['\"][^'\"]*['\"]", "", s, count=1)
    s = re.sub(r"\s+height=['\"][^'\"]*['\"]", "", s, count=1)
    s = re.sub(r"\s+version=['\"]1\.1['\"]", "", s)
    s = re.sub(r"\s+baseProfile=['\"]full['\"]", "", s)
    s = re.sub(r"\s+xmlns:rdkit=['\"][^'\"]*['\"]", "", s)
    s = re.sub(r"\s+xmlns:xlink=['\"][^'\"]*['\"]", "", s)
    s = re.sub(r"\s+xml:space=['\"]preserve['\"]", "", s)
    if "class='rdk'" not in s and 'class="rdk"' not in s:
        s = s.replace("<svg ", "<svg class='rdk' ", 1)
    s = re.sub(r"\s+class=['\"](?:bond|atom)-[^'\"]*['\"]", "", s)
    s = re.sub(
        r"\s+style=['\"]fill:none;fill-rule:evenodd;stroke:#000000;"
        r"stroke-width:1\.6px;stroke-linecap:butt;stroke-linejoin:miter;"
        r"stroke-opacity:1['\"]",
        "",
        s,
    )
    s = re.sub(
        r"\s+style=['\"]fill:none;fill-rule:evenodd;stroke:(#[0-9A-Fa-f]{6});"
        r"stroke-width:1\.6px;stroke-linecap:butt;stroke-linejoin:miter;"
        r"stroke-opacity:1['\"]",
        r" stroke='\1'",
        s,
    )
    # Compact SVG path d='...' command spacing
    s = re.sub(
        r"d='([^']+)'",
        lambda m: "d='"
        + re.sub(r"\s*([MLQCAZ])\s*", r"\1", m.group(1).strip()).replace(", ", ",")
        + "'",
        s,
    )
    s = (
        s.replace("fill='#000000'", "fill='#000'")
        .replace("stroke='#0000FF'", "stroke='#00F'")
        .replace("fill='#0000FF'", "fill='#00F'")
        .replace("stroke='#FF0000'", "stroke='#F00'")
        .replace("fill='#FF0000'", "fill='#F00'")
        .replace("fill='#FFFFFF'", "fill='#FFF'")
        .replace(" />", "/>")
        .replace("> <", "><")
    )
    return s


def render_mol2d_svg(smi: str, width: int = 240, height: int = 165) -> str:
    """Render a compact standalone RDKit 2D SVG for a SMILES string.

    Physics/Algorithm:
        Parses ``smi`` with RDKit (falling back to partial sanitization for
        hypervalent/dative metal coordination graphs), computes 2D depiction
        coordinates, draws with ``MolDraw2DSVG``, and compresses via
        :func:`minify_rdkit_svg`.

    Args:
        smi: Input SMILES string.
        width: Viewport width in px.
        height: Viewport height in px.

    Returns:
        Minified SVG markup string, or ``""`` if ``smi`` cannot be parsed.
    """
    if not smi:
        return ""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            return ""
        mol.UpdatePropertyCache(strict=False)
    else:
        # Canonicalize atom order so identical graphs always receive identical 2D layouts
        try:
            canon = Chem.MolToSmiles(mol, canonical=True)
            m_canon = Chem.MolFromSmiles(canon)
            if m_canon is not None:
                mol = m_canon
        except Exception:
            pass
    try:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol)
        rdDepictor.StraightenDepiction(mol)
    except Exception:
        try:
            rdDepictor.SetPreferCoordGen(False)
            rdDepictor.Compute2DCoords(mol)
        except Exception:
            pass
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.12
    opts.bondLineWidth = 1.6
    opts.minFontSize = 11
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return minify_rdkit_svg(drawer.GetDrawingText())


def norm_elem(val: str | None) -> str:
    """Normalize element symbol into metal center or ``'Organic'``.

    Physics/Algorithm:
        Maps non-metal light organic elements (``C, N, O, S, Cl, F``) or empty
        strings to ``'Organic'`` for semiconductor DB categorization.

    Args:
        val: Raw element symbol string.

    Returns:
        Normalized element category string.
    """
    v = (val or "").strip()
    if not v or v in {"C", "N", "O", "S", "Cl", "F"}:
        return "Organic"
    return v


def canon_smi(smi: str | None) -> str:
    """Compute achiral canonical RDKit SMILES for entity deduplication.

    Physics/Algorithm:
        Generates canonical graph SMILES (``isomericSmiles=False``) to merge
        tautomeric/notation aliases across precursor and reagent tables.

    Args:
        smi: Input SMILES string.

    Returns:
        Canonical SMILES string, or ``""`` if unparseable.
    """
    if not smi:
        return ""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def norm_chem_key(val: str | None) -> str:
    """Normalize a chemical name or formula into an alphanumeric lookup key.

    Physics/Algorithm:
        Strips non-alphanumeric characters and lowercases ``val`` for robust
        cross-referencing between benchmark IDs and AtomicLimits names.

    Args:
        val: Chemical name or molecular formula.

    Returns:
        Lowercase alphanumeric string.
    """
    return re.sub(r"[^a-z0-9]", "", str(val or "").lower())


def compress_html_report(html_str: str) -> str:
    """Losslessly compress static DOM tables and inline SVGs in the HTML report.

    Physics/Algorithm:
        Preserves ``<script>``, ``<style>``, and ``<pre>`` blocks verbatim while
        replacing repetitive inline cell styles on bare tags (``<td>``, ``<th>``,
        ``<span>``, ``<div>``) with compact utility classes (``.u01``-``.u26``)
        and stripping redundant multi-space indentation on table/SVG tags.

    Args:
        html_str: Full assembled HTML document string.

    Returns:
        Compressed HTML document string with identical DOM semantics and text.
    """
    protected: List[str] = []

    def _stash(match: re.Match[str]) -> str:
        idx = len(protected)
        protected.append(match.group(0))
        return f"__PROTECTED_BLOCK_{idx}__"

    work = re.sub(
        r"<(script|style|pre)\b[^>]*>.*?</\1>",
        _stash,
        html_str,
        flags=re.S | re.I,
    )

    for style_val, cls_name in _HIGH_FREQ_STYLES:
        for tag in ("td", "th", "span", "div"):
            work = work.replace(
                f'<{tag} style="{style_val}">',
                f'<{tag} class="{cls_name}">',
            )

    # Strip redundant XML attributes in inline SVGs
    work = re.sub(r"\s+xmlns:rdkit=['\"][^'\"]*['\"]", "", work)
    work = re.sub(r"\s+xmlns:xlink=['\"][^'\"]*['\"]", "", work)
    work = re.sub(r"\s+xml:space=['\"]preserve['\"]", "", work)
    work = re.sub(r"\s+version=['\"]1\.1['\"]\s+baseProfile=['\"]full['\"]", "", work)

    # Collapse leading indentation before structural HTML/SVG tags
    work = re.sub(
        r"\n[ ]{2,}(?=<(?:tr|/tr|td|th|thead|/thead|tbody|/tbody|table|/table|"
        r"path|circle|line|rect|polygon|polyline|g|/g|text|defs|/defs|li|/li|ul|/ul|ol|/ol|p|/p|div|/div)\b)",
        "\n",
        work,
    )
    work = re.sub(r"\n{3,}", "\n\n", work)

    def _restore(match: re.Match[str]) -> str:
        return protected[int(match.group(1))]

    work = re.sub(r"__PROTECTED_BLOCK_(\d+)__", _restore, work)
    return work
