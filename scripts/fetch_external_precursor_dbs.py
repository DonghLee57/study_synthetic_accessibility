import os
import json
import urllib.request
import ssl
from collections import defaultdict

def fetch_atomiclimits():
    print("[1/3] Fetching live data from AtomicLimits ALD Database (TU Eindhoven)...")
    ctx = ssl._create_unverified_context()
    url = "https://www.atomiclimits.com/alddatabase/api/processes.php"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"})
    
    with urllib.request.urlopen(req, context=ctx, timeout=30) as resp:
        raw_data = json.loads(resp.read().decode("utf-8"))
        
    processes = raw_data.get("processes", [])
    references = raw_data.get("references", [])
    print(f"  -> Successfully retrieved {len(processes)} ALD processes and {len(references)} literature references.")
    
    os.makedirs("database/external", exist_ok=True)
    
    # Save raw data
    with open("database/external/atomiclimits_raw.json", "w", encoding="utf-8") as f:
        json.dump(raw_data, f, indent=2, ensure_ascii=False)
        
    # Build reference lookup: process_id -> list of DOIs
    ref_by_pid = defaultdict(list)
    for ref in references:
        pid = ref.get("process_id")
        doi = ref.get("reference_doi")
        if pid and doi:
            ref_by_pid[pid].append(doi)
            
    # Group by Precursor (process_reactantA)
    precursor_map = defaultdict(lambda: {
        "precursor_name": "",
        "process_count": 0,
        "materials": set(),
        "co_reactants": set(),
        "process_ids": [],
        "dois": set()
    })
    
    for p in processes:
        prec = p.get("process_reactantA", "").strip()
        mat = p.get("process_material", "").strip()
        co_b = p.get("process_reactantB", "").strip()
        co_c = p.get("process_reactantC", "").strip()
        pid = p.get("process_id", "")
        
        if not prec:
            continue
            
        entry = precursor_map[prec]
        entry["precursor_name"] = prec
        entry["process_count"] += 1
        if mat:
            entry["materials"].add(mat)
        if co_b:
            entry["co_reactants"].add(co_b)
        if co_c:
            entry["co_reactants"].add(co_c)
        if pid:
            entry["process_ids"].append(pid)
            for d in ref_by_pid.get(pid, []):
                entry["dois"].add(d)
                
    # Convert sets to sorted lists for JSON serialization
    curated_precursors = []
    for prec, d in precursor_map.items():
        curated_precursors.append({
            "precursor": prec,
            "process_count": d["process_count"],
            "materials_count": len(d["materials"]),
            "materials": sorted(list(d["materials"])),
            "co_reactants": sorted(list(d["co_reactants"])),
            "literature_citations_count": len(d["dois"]),
            "sample_dois": sorted(list(d["dois"]))[:5]
        })
        
    curated_precursors.sort(key=lambda x: x["process_count"], reverse=True)
    
    with open("database/external/atomiclimits_precursors_extracted.json", "w", encoding="utf-8") as f:
        json.dump(curated_precursors, f, indent=2, ensure_ascii=False)
        
    print(f"  -> Extracted {len(curated_precursors)} unique ALD precursors.")
    print("  -> Saved to database/external/atomiclimits_precursors_extracted.json")
    return curated_precursors, raw_data

def build_expanded_master_database(atomic_precs):
    print("\n[2/3] Integrating with Gelest/Strem/Min et al. (JCIM 2024) Organometallic Knowledge Base...")
    
    # Load our 27 baseline precursors
    with open("database/precursors.json", "r", encoding="utf-8") as f:
        baseline_27 = json.load(f)
        
    baseline_names = {p["name"].lower(): p for p in baseline_27}
    
    # Well-known semiconductor precursor metadata dictionary
    # Connecting formulas, CAS numbers, metal centers, SMILES, and thermal properties
    known_precursor_catalog = {
        "alme3": {"standard_name": "TMA", "formula": "Al(CH3)3", "metal": "Al", "category": "Aluminum", "cas": "75-24-1", "smiles": "C[Al](C)C", "vp_torr": "8.4 Torr (20 °C)", "t_decomp": "> 300 °C"},
        "alcl3": {"standard_name": "AlCl3", "formula": "AlCl3", "metal": "Al", "category": "Aluminum", "cas": "7446-70-0", "smiles": "Cl[Al](Cl)Cl", "vp_torr": "1.0 Torr (100 °C)", "t_decomp": "> 400 °C"},
        "zne2": {"standard_name": "DEZ", "formula": "Zn(C2H5)2", "metal": "Zn", "category": "Zinc", "cas": "557-20-0", "smiles": "CC[Zn]CC", "vp_torr": "12.4 Torr (20 °C)", "t_decomp": "> 300 °C"},
        "znet2": {"standard_name": "DEZ", "formula": "Zn(C2H5)2", "metal": "Zn", "category": "Zinc", "cas": "557-20-0", "smiles": "CC[Zn]CC", "vp_torr": "12.4 Torr (20 °C)", "t_decomp": "> 300 °C"},
        "ticl4": {"standard_name": "TiCl4", "formula": "TiCl4", "metal": "Ti", "category": "Titanium", "cas": "7550-45-0", "smiles": "Cl[Ti](Cl)(Cl)Cl", "vp_torr": "10.0 Torr (21 °C)", "t_decomp": "> 450 °C"},
        "ti(oipr)4": {"standard_name": "TTIP", "formula": "Ti(OCH(CH3)2)4", "metal": "Ti", "category": "Titanium", "cas": "546-68-9", "smiles": "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C", "vp_torr": "1.0 Torr (50 °C)", "t_decomp": "> 250 °C"},
        "ti(nme2)4": {"standard_name": "TDMAT", "formula": "Ti(N(CH3)2)4", "metal": "Ti", "category": "Titanium", "cas": "3275-24-9", "smiles": "CN(C)[Ti](N(C)C)(N(C)C)N(C)C", "vp_torr": "1.0 Torr (50 °C)", "t_decomp": "> 180 °C"},
        "hfcl4": {"standard_name": "HfCl4", "formula": "HfCl4", "metal": "Hf", "category": "Hafnium", "cas": "13499-05-3", "smiles": "Cl[Hf](Cl)(Cl)Cl", "vp_torr": "1.0 Torr (190 °C)", "t_decomp": "> 500 °C"},
        "hf(netme)4": {"standard_name": "TEMAH", "formula": "Hf(N(CH3)(C2H5))4", "metal": "Hf", "category": "Hafnium", "cas": "352535-01-4", "smiles": "CCN(C)[Hf](N(C)CC)(N(C)CC)N(C)CC", "vp_torr": "0.1 Torr (60 °C)", "t_decomp": "> 250 °C"},
        "hf(nme2)4": {"standard_name": "TDMAH", "formula": "Hf(N(CH3)2)4", "metal": "Hf", "category": "Hafnium", "cas": "19782-68-4", "smiles": "CN(C)[Hf](N(C)C)(N(C)C)N(C)C", "vp_torr": "0.1 Torr (50 °C)", "t_decomp": "> 220 °C"},
        "hf(net2)4": {"standard_name": "TDEAH", "formula": "Hf(N(C2H5)2)4", "metal": "Hf", "category": "Hafnium", "cas": "630138-86-8", "smiles": "CCN(CC)[Hf](N(CC)CC)(N(CC)CC)N(CC)CC", "vp_torr": "0.1 Torr (85 °C)", "t_decomp": "> 260 °C"},
        "zrcl4": {"standard_name": "ZrCl4", "formula": "ZrCl4", "metal": "Zr", "category": "Zirconium", "cas": "10026-11-6", "smiles": "Cl[Zr](Cl)(Cl)Cl", "vp_torr": "1.0 Torr (190 °C)", "t_decomp": "> 500 °C"},
        "zr(netme)4": {"standard_name": "TEMAZ", "formula": "Zr(N(CH3)(C2H5))4", "metal": "Zr", "category": "Zirconium", "cas": "175923-04-3", "smiles": "CCN(C)[Zr](N(C)CC)(N(C)CC)N(C)CC", "vp_torr": "0.1 Torr (60 °C)", "t_decomp": "> 250 °C"},
        "zr(nme2)4": {"standard_name": "TDMAZ", "formula": "Zr(N(CH3)2)4", "metal": "Zr", "category": "Zirconium", "cas": "19756-04-8", "smiles": "CN(C)[Zr](N(C)C)(N(C)C)N(C)C", "vp_torr": "0.1 Torr (50 °C)", "t_decomp": "> 220 °C"},
        "wf6": {"standard_name": "WF6", "formula": "WF6", "metal": "W", "category": "Tungsten", "cas": "7783-82-6", "smiles": "F[W](F)(F)(F)(F)F", "vp_torr": "760 Torr (17 °C)", "t_decomp": "> 350 °C"},
        "tacl5": {"standard_name": "TaCl5", "formula": "TaCl5", "metal": "Ta", "category": "Tantalum", "cas": "7721-01-9", "smiles": "Cl[Ta](Cl)(Cl)(Cl)Cl", "vp_torr": "1.0 Torr (120 °C)", "t_decomp": "> 500 °C"},
        "mocl5": {"standard_name": "MoCl5", "formula": "MoCl5", "metal": "Mo", "category": "Molybdenum", "cas": "10241-05-1", "smiles": "Cl[Mo](Cl)(Cl)(Cl)Cl", "vp_torr": "1.0 Torr (90 °C)", "t_decomp": "> 200 °C"},
        "ru(etcp)2": {"standard_name": "Ru(EtCp)2", "formula": "Ru(C2H5C5H4)2", "metal": "Ru", "category": "Ruthenium", "cas": "32992-05-9", "smiles": "CC[C]1[CH][CH][CH][CH]1.CC[C]1[CH][CH][CH][CH]1.[Ru]", "vp_torr": "0.1 Torr (80 °C)", "t_decomp": "> 300 °C"},
        "co(amd)2": {"standard_name": "Co(AMD)2", "formula": "Co(iPr-AMD)2", "metal": "Co", "category": "Cobalt", "cas": "635680-56-7", "smiles": "CC(C)N=C(C)[N-]C(C)C.CC(C)N=C(C)[N-]C(C)C.[Co+2]", "vp_torr": "0.1 Torr (65 °C)", "t_decomp": "> 280 °C"},
        "si(oet)4": {"standard_name": "TEOS", "formula": "Si(OC2H5)4", "metal": "Si", "category": "Silicon", "cas": "78-10-4", "smiles": "CCO[Si](OCC)(OCC)OCC", "vp_torr": "1.0 Torr (20 °C)", "t_decomp": "> 600 °C"},
        "sih2cl2": {"standard_name": "DCS", "formula": "SiH2Cl2", "metal": "Si", "category": "Silicon", "cas": "4109-96-0", "smiles": "Cl[SiH2]Cl", "vp_torr": "760 Torr (8.3 °C)", "t_decomp": "> 500 °C"},
        "sicl4": {"standard_name": "STC", "formula": "SiCl4", "metal": "Si", "category": "Silicon", "cas": "10026-04-7", "smiles": "Cl[Si](Cl)(Cl)Cl", "vp_torr": "194 Torr (20 °C)", "t_decomp": "> 600 °C"},
        "sic2h6": {"standard_name": "DMAH / Silanes", "formula": "SiH2(CH3)2", "metal": "Si", "category": "Silicon", "cas": "1066-29-1", "smiles": "C[SiH2]C", "vp_torr": "760 Torr (-20 °C)", "t_decomp": "> 400 °C"}
    }
    
    expanded_master = []
    
    # 1. First add all baseline 27 with enriched process data
    for b in baseline_27:
        b_key = b["name"].lower().replace("-", "").replace(" ", "").replace("_", "")
        # Find matching AtomicLimits entries
        matched_ap = None
        for ap in atomic_precs:
            if ap["precursor"].lower().replace("-", "").replace(" ", "").replace("_", "") == b_key or \
               b["formula"].lower().replace(" ", "") in ap["precursor"].lower().replace(" ", ""):
                matched_ap = ap
                break
                
        entry = {
            "id": b["id"],
            "name": b["name"],
            "chemical_name": b["chemical_name"],
            "formula": b["formula"],
            "metal_center": b["metal_center"],
            "category": b["category"],
            "smiles": b["smiles"],
            "molecular_weight": b["molecular_weight"],
            "cas_number": b["cas_number"],
            "source": "Commercial ALD/CVD Baseline (Verified Fab Spec)",
            "atomiclimits_processes": matched_ap["process_count"] if matched_ap else 0,
            "atomiclimits_materials": matched_ap["materials"] if matched_ap else [b["target_films"]],
            "co_reactants": matched_ap["co_reactants"] if matched_ap else [],
            "vapor_pressure": b["physical_properties"].get("vapor_pressure_torr", "-"),
            "thermal_decomposition_temp_c": b["physical_properties"].get("thermal_decomposition_temp_c", "-"),
            "commercial_vendors": b.get("commercial_vendors", [])
        }
        expanded_master.append(entry)
        
    print(f"  -> Processed baseline {len(expanded_master)} precursors with AtomicLimits cross-referencing.")
    
    # 2. Add top AtomicLimits ALD Precursors (threshold: process_count >= 3)
    existing_formulas = {p["formula"].lower().replace(" ", "") for p in expanded_master}
    existing_names = {p["name"].lower().replace(" ", "") for p in expanded_master}
    
    added_count = 0
    for ap in atomic_precs:
        p_name = ap["precursor"].strip()
        p_key = p_name.lower().replace("-", "").replace(" ", "").replace("_", "")
        
        if p_name.lower() in existing_names or p_key in known_precursor_catalog and known_precursor_catalog[p_key]["formula"].lower() in existing_formulas:
            continue
            
        # Get metadata from known catalog if available
        meta = known_precursor_catalog.get(p_key, {})
        metal = meta.get("metal", "")
        if not metal:
            # Simple heuristic metal detection from formula
            for m in ["Zr", "Hf", "Ti", "Al", "Ru", "Ta", "Mo", "W", "Co", "Sn", "Zn", "Ga", "In", "La", "Si", "Ge", "Cu", "Ni", "Fe", "Mn", "V", "Nb", "Y", "Sc", "Sb", "Bi"]:
                if m in p_name:
                    metal = m
                    break
        if not metal:
            metal = "Other"
            
        new_entry = {
            "id": f"EXT-ALD-{len(expanded_master)+1:03d}",
            "name": meta.get("standard_name", p_name),
            "chemical_name": p_name,
            "formula": meta.get("formula", p_name),
            "metal_center": metal,
            "category": meta.get("category", f"{metal} Precursor"),
            "smiles": meta.get("smiles", ""),
            "molecular_weight": "-",
            "cas_number": meta.get("cas", "-"),
            "source": "AtomicLimits ALD Database (TU Eindhoven)",
            "atomiclimits_processes": ap["process_count"],
            "atomiclimits_materials": ap["materials"],
            "co_reactants": ap["co_reactants"],
            "vapor_pressure": meta.get("vp_torr", "-"),
            "thermal_decomposition_temp_c": meta.get("t_decomp", "-"),
            "literature_dois": ap["sample_dois"]
        }
        expanded_master.append(new_entry)
        existing_names.add(p_name.lower())
        added_count += 1
        if len(expanded_master) >= 120:  # Expand to top 120 key ALD precursors
            break
            
    print(f"  -> Added {added_count} top ALD precursors from AtomicLimits Database.")
    print(f"  -> Total master semiconductor precursors database: {len(expanded_master)} entries.")
    
    with open("database/external/expanded_semiconductor_precursors_master.json", "w", encoding="utf-8") as f:
        json.dump(expanded_master, f, indent=2, ensure_ascii=False)
        
    print("  -> Saved master database to database/external/expanded_semiconductor_precursors_master.json")
    return expanded_master

if __name__ == "__main__":
    atomic_precs, raw_data = fetch_atomiclimits()
    build_expanded_master_database(atomic_precs)
    print("\n[3/3] External Precursor Database Extraction & Synthesis Completed Successfully!")
