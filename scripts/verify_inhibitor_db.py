"""Verify database/external/as_ald_inhibitors.json.

Checks per entry:
  * RDKit: SMILES parses; formula, MW, Crippen logP, heavy atoms, rotatable bonds
    match the stored values.
  * PubChem (online): the CAS number resolves to the same InChIKey as the SMILES.
  * Crossref (online): every non-null DOI resolves and its first-author surname
    matches the stored ``authors`` field.

Usage:
    python scripts/verify_inhibitor_db.py            # RDKit + PubChem + Crossref
    python scripts/verify_inhibitor_db.py --offline  # RDKit only
Exit code 1 if any check fails.
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem import rdMolDescriptors as rdd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(ROOT, "database", "external", "as_ald_inhibitors.json")
UA = {"User-Agent": "synth_acc-inhibitor-verify/1.0"}


def get_json(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30) as r:
        return json.load(r)


def check_descriptors(e):
    m = Chem.MolFromSmiles(e["smiles"])
    if m is None:
        return ["SMILES does not parse"], None
    errs = []
    calc = {
        "formula": rdd.CalcMolFormula(m),
        "molecular_weight": round(Descriptors.MolWt(m), 2),
        "logp": round(Crippen.MolLogP(m), 2),
        "heavy_atom_count": m.GetNumHeavyAtoms(),
        "rotatable_bonds": rdd.CalcNumRotatableBonds(m),
    }
    for k, v in calc.items():
        stored = e.get(k)
        ok = abs(stored - v) <= 0.05 if isinstance(v, float) else stored == v
        if not ok:
            errs.append(f"{k}: stored {stored!r} != RDKit {v!r}")
    return errs, Chem.MolToInchiKey(m)


def check_cas(e, inchikey):
    cas = e.get("cas_number")
    if not cas:
        return []
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(cas) + "/property/InChIKey/JSON")
    try:
        keys = {p["InChIKey"] for p in get_json(url)["PropertyTable"]["Properties"]}
    except Exception as exc:  # noqa: BLE001
        return [f"CAS {cas}: PubChem lookup failed ({exc})"]
    return [] if inchikey in keys else [f"CAS {cas}: PubChem InChIKey {sorted(keys)} != SMILES {inchikey}"]


def check_doi(ref):
    doi = (ref or {}).get("doi")
    if not doi:
        return []
    try:
        msg = get_json("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
    except Exception as exc:  # noqa: BLE001
        return [f"DOI {doi}: not resolvable on Crossref ({exc})"]
    first = ((msg.get("author") or [{}])[0].get("family") or "").lower()
    claimed = (ref.get("authors") or "").lower()
    if first and first not in claimed:
        return [f"DOI {doi}: Crossref first author '{first}' not in stored authors '{ref.get('authors')}'"]
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="skip PubChem/Crossref checks")
    args = ap.parse_args()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    with open(DB_PATH, encoding="utf-8") as f:
        entries = json.load(f)["inhibitors"]
    failed = 0
    for e in entries:
        errs, ik = check_descriptors(e)
        if not args.offline and ik:
            errs += check_cas(e, ik)
            errs += check_doi(e.get("literature_reference"))
            time.sleep(0.3)
        status = "OK" if not errs else "FAIL"
        failed += bool(errs)
        print(f"[{status}] {e['id']} {e['abbreviation']}" + "".join(f"\n    - {x}" for x in errs))
    print(f"\n{len(entries) - failed}/{len(entries)} entries passed"
          + (" (offline: RDKit checks only)" if args.offline else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
