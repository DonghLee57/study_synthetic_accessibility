"""Build a semiconductor precursor inventory from the AtomicLimits ALD database.

Usage
-----
    python scripts/fetch_external_precursor_dbs.py           # rebuild from local raw snapshot
    python scripts/fetch_external_precursor_dbs.py --fetch   # download a fresh snapshot first

Data policy
-----------
* The raw AtomicLimits dump (``database/external/atomiclimits_raw.json``) is NOT
  committed: it is a full copy of a third-party database (including contributor
  names) obtained from the site's web API, and its redistribution terms have not
  been confirmed. Only derived aggregates and a snapshot record (date, URL,
  SHA-256, counts) are tracked. Cite the database as
  https://www.atomiclimits.com/alddatabase (DOI 10.6100/alddatabase).
* No physical property, CAS number or SMILES is invented here. The master file
  carries identifiers/properties only for the 27 project-curated baseline
  precursors (``database/precursors.json``); for the added AtomicLimits entries
  those fields are null and flagged for a later PubChem/datasheet lookup.

Method
------
* Every reactant field (A-D) of every process is read, not only reactant A:
  in multi-component processes the metal precursor is often listed in B/C
  (e.g. AlMe3 appears 76x as A but 155x in total).
* A reactant is a *precursor* if it contains a metal/metalloid centre
  (see ``NON_CENTER``); everything else (H2O, O3, plasmas, H2S, ...) is a
  *co-reactant*. For P/S/Se sources carrying a metal(loid) ligand
  (Se(SiMe3)2, P(SnMe3)3) the chalcogen/pnictogen is reported as the centre.
* The centre element is the metal/metalloid symbol at the lowest parenthesis
  depth (first occurrence), so ``Pb(N(SiMe3)2)2`` -> Pb, ``(Et3Si)2Te`` -> Te.
* Spelling variants are merged through ``atomiclimits_aliases.json``.
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import ssl
import urllib.request
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXT = os.path.join(ROOT, "database", "external")
RAW_PATH = os.path.join(EXT, "atomiclimits_raw.json")
SNAPSHOT_PATH = os.path.join(EXT, "atomiclimits_snapshot.json")
ALIAS_PATH = os.path.join(EXT, "atomiclimits_aliases.json")
INVENTORY_PATH = os.path.join(EXT, "atomiclimits_precursors_extracted.json")
MASTER_PATH = os.path.join(EXT, "expanded_semiconductor_precursors_master.json")
CURATED_ID_PATH = os.path.join(EXT, "curated_target_element_identifiers.json")
BASELINE_PATH = os.path.join(ROOT, "database", "precursors.json")

API_URL = "https://www.atomiclimits.com/alddatabase/api/processes.php"
MASTER_SIZE = 120
FOCUS_ELEMENTS = ("Ir", "Rh", "As", "Ga", "In", "Ru")

ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu"
).split()
# Elements that never count as the precursor centre (ligand / co-reactant atoms).
NON_CENTER = {"H", "He", "C", "N", "O", "F", "Ne", "P", "S", "Cl", "Ar", "Se", "Br", "Kr", "I", "Xe", "Rn", "At"}
CENTER = set(ELEMENTS) - NON_CENTER
# Chalcogen/pnictogen sources such as Se(SiMe3)2 or P(SnMe3)3: they only qualify as
# precursors because of the metal/metalloid they carry, but the film element they
# deliver is the depth-0 P/S/Se atom, so that atom is ranked as the centre.
ANION_SOURCE = {"P", "S", "Se"}
# Element symbols that collide with ligand abbreviations (iPr/nPr, OAc, Np = neopentyl,
# Pa): only accepted as the very first token of the name, e.g. Pr(thd)3.
AMBIGUOUS = {"Pr", "Ac", "Np", "Pa"}

# Baseline precursor id -> canonical AtomicLimits name (None = not present in the DB).
BASELINE_TO_ATOMICLIMITS = {
    "PREC_SI_01": "SiH3(NiPr2)",          # DIPAS
    "PREC_SI_02": "SiH2(NHtBu)2",         # BTBAS
    "PREC_SI_03": "SiH(NMe2)3",           # 3DMAS
    "PREC_SI_04": "Si(NMe2)4",            # 4DMAS
    "PREC_SI_05": "N(SiH3)3",             # TSA
    "PREC_SI_06": "Si2Cl6",               # HCDS
    "PREC_SI_07": "Si(OEt)4",             # TEOS
    "PREC_SI_08": None,                   # SAM.24: baseline formula C7H19NSi does not match BDEAS; left unmapped
    "PREC_HK_01": "Hf(NEtMe)4",           # TEMAH
    "PREC_HK_02": "Hf(NMe2)4",            # TDMAH
    "PREC_HK_03": "Hf(OtBu)4",
    "PREC_HK_04": "Zr(NEtMe)4",           # TEMAZ
    "PREC_HK_05": "Zr(NMe2)4",            # TDMAZ
    "PREC_HK_06": "Zr(CpMe)2(OMe)Me",     # ZMA
    "PREC_HK_07": "ZrCp(NMe2)3",
    "PREC_TI_01": "Ti(NMe2)4",            # TDMAT
    "PREC_TI_02": "Ti(OiPr)4",            # TTIP
    "PREC_TI_03": None,                   # CpTi(NMe2)3: only substituted-Cp analogues in the DB
    "PREC_AL_01": "AlMe3",                # TMA
    "PREC_AL_02": "AlMe2H",               # DMAH
    "PREC_MET_01": "W(NtBu)2(NMe2)2",     # BTBMW
    "PREC_MET_02": "Mo(NtBu)2(NMe2)2",    # BTBMMo
    "PREC_MET_03": "Ru(CpEt)2",           # Ru(EtCp)2
    "PREC_MET_04": "Co(iPrAMD)2",         # Co(AMD)2 = bis(N,N'-diisopropylacetamidinato)cobalt
    "PREC_MET_05": "CCTBA",
    "PREC_MET_06": "Ta(NtBu)(NEt2)3",     # TBTDET
    "PREC_EUV_01": None,                  # tBuSn(OiPr)3: not in the DB
}


def fetch_snapshot():
    """Download the AtomicLimits process list (TLS verification ON)."""
    ctx = ssl.create_default_context()
    req = urllib.request.Request(API_URL, headers={"User-Agent": "synth_acc-precursor-db/1.0"})
    with urllib.request.urlopen(req, context=ctx, timeout=60) as resp:
        body = resp.read()
    data = json.loads(body.decode("utf-8"))
    os.makedirs(EXT, exist_ok=True)
    with open(RAW_PATH, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_aliases():
    with open(ALIAS_PATH, encoding="utf-8") as f:
        cfg = json.load(f)
    aliases = {re.sub(r"\s+", "", k): v for k, v in cfg["aliases"].items()}
    overrides = {k: v for k, v in cfg.get("center_overrides", {}).items() if not k.startswith("_")}
    for name in cfg.get("coreactant_overrides", {}).get("names", []):
        overrides[name] = None
    return aliases, overrides


def canonical(name, aliases):
    key = re.sub(r"\s+", "", name)
    # "[X]" wrapper without a trailing count is pure notation
    if key.startswith("[") and key.endswith("]"):
        key = key[1:-1]
    return aliases.get(key, key)


def normalize_coreactant(name):
    s = re.sub(r"\s+", " ", name).strip()
    s = re.sub(r"\bPlasma\b", "plasma", s)
    return "Electrons" if s.lower() == "electrons" else s


def center_elements(name, overrides):
    """Return centre symbols ordered by (paren depth, position); [] = co-reactant."""
    if name in overrides:
        return [overrides[name]] if overrides[name] else []
    found, depth = [], 0
    rankable = CENTER | ANION_SOURCE
    for m in re.finditer(r"[\(\[\{\)\]\}]|[A-Z][a-z]*", name):
        tok = m.group(0)
        if tok in "([{":
            depth += 1
        elif tok in ")]}":
            depth = max(0, depth - 1)
        elif tok in rankable and (tok not in AMBIGUOUS or m.start() == 0):
            found.append((depth, m.start(), tok))
    if not any(el in CENTER for _, _, el in found):
        return []
    seen, out = set(), []
    for _, _, el in sorted(found):
        if el not in seen:
            seen.add(el)
            out.append(el)
    return out


def build_inventory(raw, aliases, overrides):
    procs = raw["processes"]
    pids = {p["process_id"] for p in procs}
    dois_by_pid = defaultdict(set)
    orphan_refs = 0
    for ref in raw["references"]:
        pid, doi = ref.get("process_id"), (ref.get("reference_doi") or "").strip()
        if pid not in pids:
            orphan_refs += 1
            continue
        if doi:
            dois_by_pid[pid].add(doi)

    inv = defaultdict(lambda: {"aliases": set(), "pids": set(), "pids_A": set(),
                               "materials": set(), "co_reactants": set(), "co_precursors": set()})
    centers = {}
    for p in procs:
        pid = p["process_id"]
        reactants = []
        for slot in "ABCD":
            v = (p.get("process_reactant" + slot) or "").strip()
            if not v:
                continue
            c = canonical(v, aliases)
            els = center_elements(c, overrides)
            reactants.append((slot, v, c, els))
        prec = [r for r in reactants if r[3]]
        co = [normalize_coreactant(r[1]) for r in reactants if not r[3]]
        for slot, v, c, els in prec:
            e = inv[c]
            centers[c] = els
            e["aliases"].add(v)
            e["pids"].add(pid)
            if slot == "A":
                e["pids_A"].add(pid)
            if p.get("process_material", "").strip():
                e["materials"].add(p["process_material"].strip())
            e["co_reactants"].update(co)
            e["co_precursors"].update(o[2] for o in prec if o[2] != c)

    out = []
    for name, e in inv.items():
        dois = sorted(set().union(*[dois_by_pid[x] for x in e["pids"]]))
        els = centers[name]
        out.append({
            "precursor": name,
            "aliases_in_source": sorted(a for a in e["aliases"] if a != name),
            "center_element": els[0],
            "center_elements": els,
            "is_elemental": name in ELEMENTS,
            "process_count": len(e["pids"]),
            "process_count_as_reactantA": len(e["pids_A"]),
            "materials": sorted(e["materials"]),
            "co_reactants": sorted(e["co_reactants"]),
            "co_precursors": sorted(e["co_precursors"]),
            "process_ids": sorted(e["pids"], key=int),
            "unique_doi_count": len(dois),
            "dois": dois,
        })
    out.sort(key=lambda x: (-x["process_count"], x["precursor"]))
    stats = {
        "processes": len(procs),
        "reference_records": len(raw["references"]),
        "unique_dois": len({(r.get("reference_doi") or "").strip() for r in raw["references"]} - {""}),
        "orphan_reference_records": orphan_refs,
        "materials": len({p["process_material"].strip() for p in procs}),
        "unique_precursors_after_alias_merge": len(out),
        "unique_reactantA_strings_raw": len({p["process_reactantA"].strip() for p in procs}),
    }
    return out, stats


def load_curated_identifiers():
    if not os.path.exists(CURATED_ID_PATH):
        return {}
    with open(CURATED_ID_PATH, encoding="utf-8") as f:
        data = json.load(f)
    return data.get("identifiers", data) if isinstance(data, dict) else {}


def _build_external_master_record(entry_id: str, inv_item: dict, origin: str, curated_ids: dict) -> dict:
    """Construct a master record from an AtomicLimits inventory item and optional curated identifiers."""
    cur = curated_ids.get(inv_item["precursor"], {})
    has_smi = bool(cur.get("smiles"))
    return {
        "id": entry_id,
        "name": inv_item["precursor"],
        "chemical_name": cur.get("chemical_name"),
        "formula": cur.get("formula"),
        "center_element": inv_item["center_element"],
        "category": cur.get("category"),
        "smiles": cur.get("smiles"),
        "cas_number": cur.get("cas_number"),
        "molecular_weight": cur.get("molecular_weight"),
        "vapor_pressure": None,
        "thermal_decomposition_temp_c": None,
        "commercial_vendors": [],
        "origin": origin,
        "property_source": cur.get("identifier_source") if has_smi else None,
        "identifier_status": (
            "verified: CAS/literature + RDKit formula & 3D conformer validated"
            if has_smi
            else "unverified: needs PubChem/datasheet lookup"
        ),
        "atomiclimits_name": inv_item["precursor"],
        "atomiclimits_aliases": inv_item["aliases_in_source"],
        "atomiclimits_process_count": inv_item["process_count"],
        "atomiclimits_materials": inv_item["materials"],
        "atomiclimits_co_reactants": inv_item["co_reactants"],
        "atomiclimits_unique_doi_count": inv_item["unique_doi_count"],
        "sample_dois": inv_item["dois"][:5],
    }


def build_master(inventory, baseline, curated_ids=None):
    if curated_ids is None:
        curated_ids = load_curated_identifiers()
    by_name = {e["precursor"]: e for e in inventory}
    master, used = [], set()
    for b in baseline:
        if b["id"] not in BASELINE_TO_ATOMICLIMITS:
            raise KeyError(f"baseline id {b['id']} missing from BASELINE_TO_ATOMICLIMITS")
        al_name = BASELINE_TO_ATOMICLIMITS[b["id"]]
        al = by_name.get(al_name) if al_name else None
        if al_name and al is None:
            raise KeyError(f"{b['id']}: mapped name {al_name!r} not found in AtomicLimits inventory")
        if al:
            used.add(al_name)
        props = b.get("physical_properties", {})
        master.append({
            "id": b["id"],
            "name": b["name"],
            "chemical_name": b["chemical_name"],
            "formula": b["formula"],
            "center_element": b["metal_center"],
            "category": b["category"],
            "smiles": b["smiles"],
            "cas_number": b["cas_number"],
            "molecular_weight": b["molecular_weight"],
            "vapor_pressure": props.get("vapor_pressure_torr"),
            "thermal_decomposition_temp_c": props.get("thermal_decomposition_temp_c"),
            "commercial_vendors": b.get("commercial_vendors", []),
            "origin": "project baseline (database/precursors.json)",
            "property_source": "database/precursors.json (project-curated; per-value sources not recorded)",
            "atomiclimits_name": al_name,
            "atomiclimits_process_count": al["process_count"] if al else 0,
            "atomiclimits_materials": al["materials"] if al else [],
            "atomiclimits_co_reactants": al["co_reactants"] if al else [],
            "atomiclimits_unique_doi_count": al["unique_doi_count"] if al else 0,
        })

    for e in inventory:
        if len(master) >= MASTER_SIZE:
            break
        if e["precursor"] in used or e["is_elemental"]:
            continue
        master.append(
            _build_external_master_record(
                f"EXT-ALD-{len(master) + 1:03d}",
                e,
                "AtomicLimits ALD database",
                curated_ids,
            )
        )
        used.add(e["precursor"])

    # Append remaining non-elemental entries for priority elements (Ir, Rh, As, Ga, In, Ru)
    focus_set = set(FOCUS_ELEMENTS)
    for e in inventory:
        if e["precursor"] in used or e["is_elemental"]:
            continue
        if e["center_element"] not in focus_set:
            continue
        master.append(
            _build_external_master_record(
                f"EXT-ALD-{len(master) + 1:03d}",
                e,
                "AtomicLimits ALD database (focus-element expansion)",
                curated_ids,
            )
        )
        used.add(e["precursor"])
    return master


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fetch", action="store_true", help="download a fresh AtomicLimits snapshot first")
    ap.add_argument(
        "--from-inventory",
        action="store_true",
        help="rebuild master DB from tracked extracted inventory when atomiclimits_raw.json is absent",
    )
    args = ap.parse_args()

    prev = {}
    if os.path.exists(SNAPSHOT_PATH):
        with open(SNAPSHOT_PATH, encoding="utf-8") as f:
            prev = json.load(f)

    fetched_at = None
    if args.fetch:
        fetched_at = fetch_snapshot()

    use_inventory_only = args.from_inventory or (not args.fetch and not os.path.exists(RAW_PATH))
    if use_inventory_only:
        if not os.path.exists(INVENTORY_PATH):
            raise SystemExit(f"{INVENTORY_PATH} and {RAW_PATH} not found; run with --fetch")
        with open(INVENTORY_PATH, encoding="utf-8") as f:
            inventory = json.load(f)
        stats = {
            k: prev[k]
            for k in (
                "processes",
                "reference_records",
                "unique_dois",
                "orphan_reference_records",
                "materials",
                "unique_precursors_after_alias_merge",
                "unique_reactantA_strings_raw",
            )
            if k in prev
        }
        digest = prev.get("raw_sha256", "unknown")
        fetched_at = prev.get("fetched_at", "unknown")
    else:
        with open(RAW_PATH, encoding="utf-8") as f:
            raw = json.load(f)
        aliases, overrides = load_aliases()
        inventory, stats = build_inventory(raw, aliases, overrides)
        digest = sha256_of(RAW_PATH)
        if fetched_at is None:
            fetched_at = prev.get("fetched_at") if prev.get("raw_sha256") == digest else "unknown"

    with open(BASELINE_PATH, encoding="utf-8") as f:
        baseline = json.load(f)
    curated_ids = load_curated_identifiers()
    master = build_master(inventory, baseline, curated_ids=curated_ids)

    snapshot = {
        "source": "AtomicLimits ALD database, https://www.atomiclimits.com/alddatabase (DOI 10.6100/alddatabase)",
        "endpoint": API_URL,
        "fetched_at": fetched_at,
        "raw_sha256": digest,
        "raw_tracked_in_git": False,
        "redistribution_terms": "not confirmed; raw dump kept local only",
        **stats,
        "master_entries": len(master),
        "master_top_frequency_cutoff": MASTER_SIZE,
        "master_focus_elements": list(FOCUS_ELEMENTS),
        "master_baseline_entries": sum(1 for m in master if m["origin"].startswith("project")),
        "master_baseline_matched_in_atomiclimits": sum(1 for m in master if m["origin"].startswith("project") and m["atomiclimits_name"]),
        "master_curated_smiles_entries": sum(1 for m in master if m.get("smiles")),
    }

    outputs = [(MASTER_PATH, master), (SNAPSHOT_PATH, snapshot)]
    if not use_inventory_only:
        outputs.insert(0, (INVENTORY_PATH, inventory))
    for path, obj in outputs:
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            json.dump(obj, f, indent=2, ensure_ascii=False)
            f.write("\n")
    print(json.dumps(snapshot, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
