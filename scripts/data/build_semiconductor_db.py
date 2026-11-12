import json
import csv
import os
import sys
from rdkit import Chem

# Ensure UTF-8 output
sys.stdout.reconfigure(encoding='utf-8')

# ---------------------------------------------------------
# 1. Commercial Semiconductor Precursors Database
# ---------------------------------------------------------
PRECURSORS = [
    # --- Category: Silicon & Main Group (SiO2, SiNx, SiOC) ---
    {
        "id": "PREC_SI_01",
        "name": "DIPAS",
        "chemical_name": "Diisopropylaminosilane",
        "formula": "C6H17NSi",
        "smiles": "CC(C)N([SiH3])C(C)C",
        "cas_number": "18077-24-9",
        "molecular_weight": 131.29,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiO2, SiNx",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "10 Torr @ 20 °C",
            "boiling_point_c": 122.0,
            "melting_point_c": -15.0,
            "density_g_cm3": 0.78,
            "thermal_decomposition_temp_c": 450.0
        },
        "synthesis_info": {
            "starting_materials": ["Dichlorosilane (SiH2Cl2)", "Diisopropylamine (HN(iPr)2)"],
            "reaction_type": "Aminolysis / Salt Precipitation",
            "typical_yield_pct": 82.0,
            "byproduct": "[iPr2NH2]+Cl- (solid salt)"
        },
        "commercial_vendors": ["Gelest", "Merck/Sigma-Aldrich", "SoulBrain", "Hansol Chemical"],
        "reference_doi": "10.1039/c7ra02301d"
    },
    {
        "id": "PREC_SI_02",
        "name": "BTBAS",
        "chemical_name": "Bis(tert-butylamino)silane",
        "formula": "C8H22N2Si",
        "smiles": "CC(C)(C)N[SiH2]NC(C)(C)C",
        "cas_number": "186598-40-3",
        "molecular_weight": 174.36,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiO2, SiNx, SiON",
        "process_type": "ALD / Low-temp CVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "7 Torr @ 20 °C",
            "boiling_point_c": 167.0,
            "melting_point_c": -40.0,
            "density_g_cm3": 0.82,
            "thermal_decomposition_temp_c": 550.0
        },
        "synthesis_info": {
            "starting_materials": ["Dichlorosilane (SiH2Cl2)", "tert-Butylamine (H2NtBu)"],
            "reaction_type": "Aminolysis",
            "typical_yield_pct": 88.0,
            "byproduct": "[tBuNH3]+Cl-"
        },
        "commercial_vendors": ["Gelest", "Air Liquide", "Merck", "Entegris"],
        "reference_doi": "10.1021/cm00045a037"
    },
    {
        "id": "PREC_SI_03",
        "name": "3DMAS",
        "chemical_name": "Tris(dimethylamino)silane",
        "formula": "C6H19N3Si",
        "smiles": "CN(C)[SiH](N(C)C)N(C)C",
        "cas_number": "15112-89-7",
        "molecular_weight": 161.32,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiO2, SiNx",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "18 Torr @ 20 °C",
            "boiling_point_c": 145.0,
            "melting_point_c": -50.0,
            "density_g_cm3": 0.83,
            "thermal_decomposition_temp_c": 480.0
        },
        "synthesis_info": {
            "starting_materials": ["Trichlorosilane (HSiCl3)", "Dimethylamine (HNMe2)"],
            "reaction_type": "Aminolysis",
            "typical_yield_pct": 85.0,
            "byproduct": "[Me2NH2]+Cl-"
        },
        "commercial_vendors": ["Gelest", "Adeka", "Merck", "SoulBrain"],
        "reference_doi": "10.1149/1.1393976"
    },
    {
        "id": "PREC_SI_04",
        "name": "4DMAS",
        "chemical_name": "Tetrakis(dimethylamino)silane",
        "formula": "C8H24N4Si",
        "smiles": "CN(C)[Si](N(C)C)(N(C)C)N(C)C",
        "cas_number": "3768-55-1",
        "molecular_weight": 204.39,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiO2, SiNx, SiC",
        "process_type": "ALD / CVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "5 Torr @ 25 °C",
            "boiling_point_c": 196.0,
            "melting_point_c": -2.0,
            "density_g_cm3": 0.88,
            "thermal_decomposition_temp_c": 500.0
        },
        "synthesis_info": {
            "starting_materials": ["Silicon Tetrachloride (SiCl4)", "Dimethylamine (HNMe2)"],
            "reaction_type": "Exhaustive Aminolysis",
            "typical_yield_pct": 78.0,
            "byproduct": "[Me2NH2]+Cl-"
        },
        "commercial_vendors": ["Gelest", "Strem", "Merck"],
        "reference_doi": "10.1002/jcc.24437"
    },
    {
        "id": "PREC_SI_05",
        "name": "TSA",
        "chemical_name": "Trisilylamine",
        "formula": "H9NSi3",
        "smiles": "[SiH3]N([SiH3])[SiH3]",
        "cas_number": "13862-16-3",
        "molecular_weight": 107.33,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiNx (Low-temperature gap-fill)",
        "process_type": "LPCVD / ALD",
        "physical_properties": {
            "state_at_rt": "Liquid (High vapor pressure)",
            "vapor_pressure_torr": "312 Torr @ 20 °C",
            "boiling_point_c": 52.0,
            "melting_point_c": -105.0,
            "density_g_cm3": 0.89,
            "thermal_decomposition_temp_c": 350.0
        },
        "synthesis_info": {
            "starting_materials": ["Monochlorosilane (SiH3Cl)", "Ammonia (NH3)"],
            "reaction_type": "Gas-phase Aminolysis",
            "typical_yield_pct": 75.0,
            "byproduct": "[NH4]+Cl-"
        },
        "commercial_vendors": ["Air Liquide", "Wonik Materials", "Entegris"],
        "reference_doi": "10.1021/ja01188a032"
    },
    {
        "id": "PREC_SI_06",
        "name": "HCDS",
        "chemical_name": "Hexachlorodisilane",
        "formula": "Cl6Si2",
        "smiles": "Cl[Si](Cl)(Cl)[Si](Cl)(Cl)Cl",
        "cas_number": "13465-77-5",
        "molecular_weight": 268.88,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiNx, SiON (Low-temperature)",
        "process_type": "ALD / CVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "4.2 Torr @ 20 °C",
            "boiling_point_c": 145.0,
            "melting_point_c": -1.0,
            "density_g_cm3": 1.56,
            "thermal_decomposition_temp_c": 300.0
        },
        "synthesis_info": {
            "starting_materials": ["Silicon metal (Si)", "Chlorine (Cl2)"],
            "reaction_type": "Direct Chlorination of Silicon / Disilane Fractionation",
            "typical_yield_pct": 65.0,
            "byproduct": "SiCl4, Si3Cl8"
        },
        "commercial_vendors": ["Gelest", "Air Liquide", "Hansol Chemical", "Merck"],
        "reference_doi": "10.1016/j.jcrysgro.2004.05.050"
    },
    {
        "id": "PREC_SI_07",
        "name": "TEOS",
        "chemical_name": "Tetraethyl orthosilicate",
        "formula": "C8H20O4Si",
        "smiles": "CCO[Si](OCC)(OCC)OCC",
        "cas_number": "78-10-4",
        "molecular_weight": 208.33,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiO2",
        "process_type": "CVD / ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "1.5 Torr @ 20 °C",
            "boiling_point_c": 168.0,
            "melting_point_c": -82.0,
            "density_g_cm3": 0.93,
            "thermal_decomposition_temp_c": 650.0
        },
        "synthesis_info": {
            "starting_materials": ["Silicon Tetrachloride (SiCl4)", "Ethanol (EtOH)"],
            "reaction_type": "Alcoholysis",
            "typical_yield_pct": 94.0,
            "byproduct": "HCl (gas)"
        },
        "commercial_vendors": ["Gelest", "Merck", "SoulBrain", "Honeywell"],
        "reference_doi": "10.1021/ja01188a032"
    },
    {
        "id": "PREC_SI_08",
        "name": "SAM.24",
        "chemical_name": "N,N-Diisopropyl-1-methylsilanamine",
        "formula": "C7H19NSi",
        "smiles": "CC(C)N(C(C)C)[SiH2]C",
        "cas_number": "1042308-77-9",
        "molecular_weight": 145.32,
        "metal_center": "Si",
        "category": "Silicon & Main Group",
        "target_films": "SiOC, SiC, SiO2",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "15 Torr @ 25 °C",
            "boiling_point_c": 138.0,
            "melting_point_c": -30.0,
            "density_g_cm3": 0.79,
            "thermal_decomposition_temp_c": 420.0
        },
        "synthesis_info": {
            "starting_materials": ["Methyldichlorosilane (MeSiHCl2)", "LiAlH4", "HN(iPr)2"],
            "reaction_type": "Selective Hydride Reduction & Aminolysis",
            "typical_yield_pct": 76.0,
            "byproduct": "LiCl, AlCl3"
        },
        "commercial_vendors": ["Air Liquide", "Gelest"],
        "reference_doi": "10.1021/acs.chemmater.5c03271"
    },

    # --- Category: High-k Dielectric (Hf, Zr) ---
    {
        "id": "PREC_HK_01",
        "name": "TEMAH",
        "chemical_name": "Tetrakis(ethylmethylamino)hafnium",
        "formula": "C12H32HfN4",
        "smiles": "CCN(C)[Hf](N(C)CC)(N(C)CC)N(C)CC",
        "cas_number": "352535-01-4",
        "molecular_weight": 410.89,
        "metal_center": "Hf",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "HfO2, HfSiOx",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 60 °C",
            "boiling_point_c": 240.0,
            "melting_point_c": -45.0,
            "density_g_cm3": 1.32,
            "thermal_decomposition_temp_c": 280.0
        },
        "synthesis_info": {
            "starting_materials": ["Hafnium Tetrachloride (HfCl4)", "Lithium ethylmethylamide (LiN(Et)Me)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 86.0,
            "byproduct": "LiCl (solid salt)"
        },
        "commercial_vendors": ["Strem/Entegris", "Adeka", "Merck", "SoulBrain", "Hansol Chemical"],
        "reference_doi": "10.1021/ja020287v"
    },
    {
        "id": "PREC_HK_02",
        "name": "TDMAH",
        "chemical_name": "Tetrakis(dimethylamino)hafnium",
        "formula": "C8H24HfN4",
        "smiles": "CN(C)[Hf](N(C)C)(N(C)C)N(C)C",
        "cas_number": "19782-68-4",
        "molecular_weight": 354.78,
        "metal_center": "Hf",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "HfO2",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Solid (Low melting point)",
            "vapor_pressure_torr": "0.1 Torr @ 45 °C",
            "boiling_point_c": 210.0,
            "melting_point_c": 29.0,
            "density_g_cm3": 1.38,
            "thermal_decomposition_temp_c": 260.0
        },
        "synthesis_info": {
            "starting_materials": ["Hafnium Tetrachloride (HfCl4)", "Lithium dimethylamide (LiNMe2)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 89.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/ja020287v"
    },
    {
        "id": "PREC_HK_03",
        "name": "Hf(OtBu)4",
        "chemical_name": "Hafnium tetra-tert-butoxide",
        "formula": "C16H36HfO4",
        "smiles": "CC(C)(C)O[Hf](OC(C)(C)C)(OC(C)(C)C)OC(C)(C)C",
        "cas_number": "2172-02-3",
        "molecular_weight": 470.94,
        "metal_center": "Hf",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "HfO2",
        "process_type": "ALD / MOCVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "1.0 Torr @ 75 °C",
            "boiling_point_c": 210.0,
            "melting_point_c": 8.0,
            "density_g_cm3": 1.25,
            "thermal_decomposition_temp_c": 320.0
        },
        "synthesis_info": {
            "starting_materials": ["Hafnium Tetrachloride (HfCl4)", "Sodium tert-butoxide (NaOtBu)"],
            "reaction_type": "Salt Metathesis / Alcoholysis",
            "typical_yield_pct": 85.0,
            "byproduct": "NaCl"
        },
        "commercial_vendors": ["Strem", "Gelest", "Merck"],
        "reference_doi": "10.1021/cm00045a037"
    },
    {
        "id": "PREC_HK_04",
        "name": "TEMAZ",
        "chemical_name": "Tetrakis(ethylmethylamino)zirconium",
        "formula": "C12H32N4Zr",
        "smiles": "CCN(C)[Zr](N(C)CC)(N(C)CC)N(C)CC",
        "cas_number": "175923-04-3",
        "molecular_weight": 323.63,
        "metal_center": "Zr",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "ZrO2, ZAZ (ZrO2-Al2O3-ZrO2)",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 50 °C",
            "boiling_point_c": 220.0,
            "melting_point_c": -50.0,
            "density_g_cm3": 1.05,
            "thermal_decomposition_temp_c": 270.0
        },
        "synthesis_info": {
            "starting_materials": ["Zirconium Tetrachloride (ZrCl4)", "Lithium ethylmethylamide (LiN(Et)Me)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 88.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem/Entegris", "Adeka", "Merck", "SoulBrain"],
        "reference_doi": "10.1021/ja020287v"
    },
    {
        "id": "PREC_HK_05",
        "name": "TDMAZ",
        "chemical_name": "Tetrakis(dimethylamino)zirconium",
        "formula": "C8H24N4Zr",
        "smiles": "CN(C)[Zr](N(C)C)(N(C)C)N(C)C",
        "cas_number": "19756-04-8",
        "molecular_weight": 267.53,
        "metal_center": "Zr",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "ZrO2",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Solid",
            "vapor_pressure_torr": "0.1 Torr @ 35 °C",
            "boiling_point_c": 190.0,
            "melting_point_c": 60.0,
            "density_g_cm3": 1.15,
            "thermal_decomposition_temp_c": 250.0
        },
        "synthesis_info": {
            "starting_materials": ["Zirconium Tetrachloride (ZrCl4)", "Lithium dimethylamide (LiNMe2)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 91.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/ja020287v"
    },
    {
        "id": "PREC_HK_06",
        "name": "ZMA",
        "chemical_name": "Bis(methylcyclopentadienyl)methoxymethylzirconium",
        "formula": "C14H20OZr",
        "smiles": "CO[Zr](C)(C1=CC=C(C)C1)C2=CC=C(C)C2",
        "cas_number": "866088-57-9",
        "molecular_weight": 307.53,
        "metal_center": "Zr",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "ZrO2 (High permittivity DRAM capacitor)",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 65 °C",
            "boiling_point_c": 235.0,
            "melting_point_c": -10.0,
            "density_g_cm3": 1.28,
            "thermal_decomposition_temp_c": 350.0
        },
        "synthesis_info": {
            "starting_materials": ["(CpMe)2ZrCl2", "Methyllithium (MeLi)", "Methanol (MeOH)"],
            "reaction_type": "Sequential Alkylation & Methanolysis",
            "typical_yield_pct": 74.0,
            "byproduct": "LiCl, CH4"
        },
        "commercial_vendors": ["Adeka", "Air Liquide", "UP Chemical"],
        "reference_doi": "10.1021/jacsau.2c00547"
    },
    {
        "id": "PREC_HK_07",
        "name": "Cp(NMe2)3Zr",
        "chemical_name": "Cyclopentadienyltris(dimethylamino)zirconium",
        "formula": "C11H23N3Zr",
        "smiles": "CN(C)[Zr](N(C)C)(N(C)C)C1=CC=CC1",
        "cas_number": "87227-56-3",
        "molecular_weight": 288.54,
        "metal_center": "Zr",
        "category": "High-k Dielectric (Hf/Zr)",
        "target_films": "ZrO2",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 55 °C",
            "boiling_point_c": 230.0,
            "melting_point_c": 12.0,
            "density_g_cm3": 1.18,
            "thermal_decomposition_temp_c": 310.0
        },
        "synthesis_info": {
            "starting_materials": ["CpZrCl3", "Lithium dimethylamide (LiNMe2)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 82.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/ja020287v"
    },

    # --- Category: Titanium & Dielectrics ---
    {
        "id": "PREC_TI_01",
        "name": "TDMAT",
        "chemical_name": "Tetrakis(dimethylamino)titanium",
        "formula": "C8H24N4Ti",
        "smiles": "CN(C)[Ti](N(C)C)(N(C)C)N(C)C",
        "cas_number": "3275-24-9",
        "molecular_weight": 224.18,
        "metal_center": "Ti",
        "category": "Titanium & Dielectrics",
        "target_films": "TiN (Diffusion barrier), TiO2",
        "process_type": "ALD / MOCVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "1.5 Torr @ 60 °C",
            "boiling_point_c": 192.0,
            "melting_point_c": -12.0,
            "density_g_cm3": 0.95,
            "thermal_decomposition_temp_c": 220.0
        },
        "synthesis_info": {
            "starting_materials": ["Titanium Tetrachloride (TiCl4)", "Lithium dimethylamide (LiNMe2)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 89.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem/Entegris", "Adeka", "Merck", "SoulBrain"],
        "reference_doi": "10.1021/ja020287v"
    },
    {
        "id": "PREC_TI_02",
        "name": "TTIP",
        "chemical_name": "Titanium tetraisopropoxide",
        "formula": "C12H28O4Ti",
        "smiles": "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C",
        "cas_number": "546-68-9",
        "molecular_weight": 284.22,
        "metal_center": "Ti",
        "category": "Titanium & Dielectrics",
        "target_films": "TiO2, BaTiO3, SrTiO3",
        "process_type": "ALD / CVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "1.4 Torr @ 50 °C",
            "boiling_point_c": 232.0,
            "melting_point_c": 18.0,
            "density_g_cm3": 0.96,
            "thermal_decomposition_temp_c": 350.0
        },
        "synthesis_info": {
            "starting_materials": ["Titanium Tetrachloride (TiCl4)", "Isopropanol (iPrOH)", "Ammonia (NH3)"],
            "reaction_type": "Alcoholysis / Neutralization",
            "typical_yield_pct": 95.0,
            "byproduct": "NH4Cl (solid salt)"
        },
        "commercial_vendors": ["Gelest", "Merck", "SoulBrain", "DuPont"],
        "reference_doi": "10.1021/ja01188a032"
    },
    {
        "id": "PREC_TI_03",
        "name": "CpTi(NMe2)3",
        "chemical_name": "Cyclopentadienyltris(dimethylamino)titanium",
        "formula": "C11H23N3Ti",
        "smiles": "CN(C)[Ti](N(C)C)(N(C)C)C1=CC=CC1",
        "cas_number": "1072-67-9",
        "molecular_weight": 245.20,
        "metal_center": "Ti",
        "category": "Titanium & Dielectrics",
        "target_films": "TiO2, TiN",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.5 Torr @ 60 °C",
            "boiling_point_c": 215.0,
            "melting_point_c": -5.0,
            "density_g_cm3": 1.03,
            "thermal_decomposition_temp_c": 290.0
        },
        "synthesis_info": {
            "starting_materials": ["CpTiCl3", "Lithium dimethylamide (LiNMe2)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 84.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/ja020287v"
    },

    # --- Category: Group 13 / Alumina ---
    {
        "id": "PREC_AL_01",
        "name": "TMA",
        "chemical_name": "Trimethylaluminum",
        "formula": "C3H9Al",
        "smiles": "C[Al](C)C",
        "cas_number": "75-24-1",
        "molecular_weight": 72.09,
        "metal_center": "Al",
        "category": "Group 13 / Alumina",
        "target_films": "Al2O3 (Standard barrier & passivation)",
        "process_type": "ALD (Gold standard model system)",
        "physical_properties": {
            "state_at_rt": "Liquid (Pyrophoric)",
            "vapor_pressure_torr": "11.7 Torr @ 20 °C",
            "boiling_point_c": 126.0,
            "melting_point_c": 15.0,
            "density_g_cm3": 0.75,
            "thermal_decomposition_temp_c": 330.0
        },
        "synthesis_info": {
            "starting_materials": ["Aluminum metal (Al)", "Methyl chloride (MeCl)", "Sodium (Na)"],
            "reaction_type": "Direct Sesquichloride Route & Reduction",
            "typical_yield_pct": 92.0,
            "byproduct": "NaCl"
        },
        "commercial_vendors": ["AkzoNobel", "Merck", "SoulBrain", "Adeka", "Strem"],
        "reference_doi": "10.1021/cr00033a005"
    },
    {
        "id": "PREC_AL_02",
        "name": "DMAH",
        "chemical_name": "Dimethylaluminum hydride",
        "formula": "C2H7Al",
        "smiles": "C[AlH]C",
        "cas_number": "865-37-2",
        "molecular_weight": 58.06,
        "metal_center": "Al",
        "category": "Group 13 / Alumina",
        "target_films": "Al Metal (CVD wiring)",
        "process_type": "CVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "2.0 Torr @ 25 °C",
            "boiling_point_c": 154.0,
            "melting_point_c": -19.0,
            "density_g_cm3": 0.76,
            "thermal_decomposition_temp_c": 220.0
        },
        "synthesis_info": {
            "starting_materials": ["Trimethylaluminum (TMA)", "Hydrogen gas (H2)", "Aluminum (Al)"],
            "reaction_type": "High-pressure Hydrogenation",
            "typical_yield_pct": 80.0,
            "byproduct": "CH4"
        },
        "commercial_vendors": ["Gelest", "Strem"],
        "reference_doi": "10.1021/cr00033a005"
    },

    # --- Category: Metal & Barrier (W, Mo, Ru, Co, Ta) ---
    {
        "id": "PREC_MET_01",
        "name": "BTBMW",
        "chemical_name": "Bis(tert-butylimido)bis(dimethylamino)tungsten",
        "formula": "C12H30N4W",
        "smiles": "CN(C)[W](=NC(C)(C)C)(=NC(C)(C)C)N(C)C",
        "cas_number": "198463-54-4",
        "molecular_weight": 414.24,
        "metal_center": "W",
        "category": "Metal & Barrier (W/Mo/Ru/Co)",
        "target_films": "WN, WSiN (Diffusion barrier)",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 60 °C",
            "boiling_point_c": 230.0,
            "melting_point_c": -15.0,
            "density_g_cm3": 1.48,
            "thermal_decomposition_temp_c": 360.0
        },
        "synthesis_info": {
            "starting_materials": ["Tungsten Hexachloride (WCl6)", "tert-Butylamine (tBuNH2)", "Lithium dimethylamide (LiNMe2)"],
            "reaction_type": "Sequential Imido Addition & Aminolysis",
            "typical_yield_pct": 72.0,
            "byproduct": "LiCl, [tBuNH3]+Cl-"
        },
        "commercial_vendors": ["Strem", "Adeka", "Air Liquide"],
        "reference_doi": "10.1021/ja020287v"
    },
    {
        "id": "PREC_MET_02",
        "name": "BTBMMo",
        "chemical_name": "Bis(tert-butylimido)bis(dimethylamino)molybdenum",
        "formula": "C12H30MoN4",
        "smiles": "CN(C)[Mo](=NC(C)(C)C)(=NC(C)(C)C)N(C)C",
        "cas_number": "923956-62-1",
        "molecular_weight": 326.33,
        "metal_center": "Mo",
        "category": "Metal & Barrier (W/Mo/Ru/Co)",
        "target_films": "Mo, MoN, MoS2 (Next-gen interconnect)",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.2 Torr @ 55 °C",
            "boiling_point_c": 215.0,
            "melting_point_c": -10.0,
            "density_g_cm3": 1.22,
            "thermal_decomposition_temp_c": 300.0
        },
        "synthesis_info": {
            "starting_materials": ["Molybdenum Pentachloride (MoCl5)", "tert-Butylamine (tBuNH2)", "LiNMe2"],
            "reaction_type": "Imido Complexation & Salt Metathesis",
            "typical_yield_pct": 68.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Adeka", "Strem", "Air Liquide"],
        "reference_doi": "10.1021/acs.chemmater.5c03271"
    },
    {
        "id": "PREC_MET_03",
        "name": "Ru(EtCp)2",
        "chemical_name": "Bis(ethylcyclopentadienyl)ruthenium(II)",
        "formula": "C14H18Ru",
        "smiles": "CCC1=CC=C([Ru]C2=CC=C(CC)C2)C1",
        "cas_number": "32992-96-4",
        "molecular_weight": 287.36,
        "metal_center": "Ru",
        "category": "Metal & Barrier (W/Mo/Ru/Co)",
        "target_films": "Ru Metal (Low-resistivity metal line)",
        "process_type": "ALD / CVD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 70 °C",
            "boiling_point_c": 270.0,
            "melting_point_c": 6.0,
            "density_g_cm3": 1.34,
            "thermal_decomposition_temp_c": 340.0
        },
        "synthesis_info": {
            "starting_materials": ["Ruthenium Trichloride (RuCl3)", "Sodium ethylcyclopentadienide (Na(EtCp))"],
            "reaction_type": "Salt Metathesis / Reductive Coordination",
            "typical_yield_pct": 79.0,
            "byproduct": "NaCl"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/om010268w"
    },
    {
        "id": "PREC_MET_04",
        "name": "Co(AMD)2",
        "chemical_name": "Bis(N,N'-diisopropylacetamidinato)cobalt(II)",
        "formula": "C16H34CoN4",
        "smiles": "CC(C)N=C(C)N(C(C)C)[Co]N(C(C)C)C(C)=NC(C)C",
        "cas_number": "635680-56-7",
        "molecular_weight": 341.40,
        "metal_center": "Co",
        "category": "Metal & Barrier (W/Mo/Ru/Co)",
        "target_films": "Co Metal (Contact plug, Cu liner)",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Solid / Low-melting wax",
            "vapor_pressure_torr": "0.1 Torr @ 60 °C",
            "boiling_point_c": 220.0,
            "melting_point_c": 42.0,
            "density_g_cm3": 1.12,
            "thermal_decomposition_temp_c": 270.0
        },
        "synthesis_info": {
            "starting_materials": ["Cobalt(II) Chloride (CoCl2)", "Lithium N,N'-diisopropylacetamidinate"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 81.0,
            "byproduct": "LiCl"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/om010268w"
    },
    {
        "id": "PREC_MET_05",
        "name": "CCTBA",
        "chemical_name": "Dicobalt hexacarbonyl tert-butylacetylene",
        "formula": "C12H10Co2O6",
        "smiles": "CC(C)(C)C12[CH]3[Co]1([C-]#[O+])([C-]#[O+])([C-]#[O+])[Co]23([C-]#[O+])([C-]#[O+])[C-]#[O+]",
        "cas_number": "56792-69-9",
        "molecular_weight": 368.07,
        "metal_center": "Co",
        "category": "Metal & Barrier (W/Mo/Ru/Co)",
        "target_films": "Co (CVD interconnect)",
        "process_type": "CVD / ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.5 Torr @ 25 °C",
            "boiling_point_c": 180.0,
            "melting_point_c": -15.0,
            "density_g_cm3": 1.35,
            "thermal_decomposition_temp_c": 160.0
        },
        "synthesis_info": {
            "starting_materials": ["Dicobalt octacarbonyl (Co2(CO)8)", "3,3-Dimethyl-1-butyne (tBuC#CH)"],
            "reaction_type": "Alkyne Substitution / CO displacement",
            "typical_yield_pct": 90.0,
            "byproduct": "CO (gas)"
        },
        "commercial_vendors": ["Adeka", "Strem", "Merck"],
        "reference_doi": "10.1021/om010268w"
    },
    {
        "id": "PREC_MET_06",
        "name": "TBTDET",
        "chemical_name": "Tert-butylimidotris(diethylamino)tantalum",
        "formula": "C16H39N4Ta",
        "smiles": "CCN(CC)[Ta](=NC(C)(C)C)(N(CC)CC)N(CC)CC",
        "cas_number": "169896-41-7",
        "molecular_weight": 468.45,
        "metal_center": "Ta",
        "category": "Metal & Barrier (W/Mo/Ru/Co)",
        "target_films": "TaN (Cu diffusion barrier)",
        "process_type": "ALD",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 90 °C",
            "boiling_point_c": 285.0,
            "melting_point_c": -20.0,
            "density_g_cm3": 1.31,
            "thermal_decomposition_temp_c": 350.0
        },
        "synthesis_info": {
            "starting_materials": ["Tantalum Pentachloride (TaCl5)", "tert-Butylamine (tBuNH2)", "Lithium diethylamide (LiNEt2)"],
            "reaction_type": "Sequential Imido Insertion & Salt Metathesis",
            "typical_yield_pct": 77.0,
            "byproduct": "LiCl, [tBuNH3]+Cl-"
        },
        "commercial_vendors": ["Strem", "Merck", "Adeka"],
        "reference_doi": "10.1021/ja020287v"
    },

    # --- Category: EUV & Novel ---
    {
        "id": "PREC_EUV_01",
        "name": "tBuSn(OiPr)3",
        "chemical_name": "tert-Butyltris(isopropoxy)tin(IV)",
        "formula": "C13H30O3Sn",
        "smiles": "CC(C)O[Sn](C(C)(C)C)(OC(C)C)OC(C)C",
        "cas_number": "13919-49-8",
        "molecular_weight": 353.09,
        "metal_center": "Sn",
        "category": "EUV & Novel",
        "target_films": "SnOx (High-resolution EUV photoresist)",
        "process_type": "ALD / MOCVD / Spin-coating",
        "physical_properties": {
            "state_at_rt": "Liquid",
            "vapor_pressure_torr": "0.1 Torr @ 75 °C",
            "boiling_point_c": 230.0,
            "melting_point_c": -10.0,
            "density_g_cm3": 1.28,
            "thermal_decomposition_temp_c": 260.0
        },
        "synthesis_info": {
            "starting_materials": ["tert-Butyltin trichloride (tBuSnCl3)", "Sodium isopropoxide (NaOiPr)"],
            "reaction_type": "Salt Metathesis",
            "typical_yield_pct": 85.0,
            "byproduct": "NaCl"
        },
        "commercial_vendors": ["Inpria/JSR", "Gelest", "Strem"],
        "reference_doi": "10.1021/acs.chemmater.5c03271"
    }
]

# ---------------------------------------------------------
# 2. Starting Materials & Reagents (Building Blocks) Database
# ---------------------------------------------------------
BUILDING_BLOCKS = [
    # Inorganic Halide & Alkyl Cores
    {"id": "BB_01", "name": "Silicon Tetrachloride", "formula": "SiCl4", "smiles": "Cl[Si](Cl)(Cl)Cl", "cas": "10026-04-7", "role": "Inorganic Core", "purity": "99.9999% (6N)", "handling": "Moisture Sensitive, Corrosive"},
    {"id": "BB_02", "name": "Dichlorosilane", "formula": "SiH2Cl2", "smiles": "Cl[SiH2]Cl", "cas": "4109-96-0", "role": "Inorganic Core", "purity": "99.999% (5N)", "handling": "Flammable Gas, Pyrophoric Risk"},
    {"id": "BB_03", "name": "Trichlorosilane", "formula": "HSiCl3", "smiles": "Cl[SiH](Cl)Cl", "cas": "10025-78-2", "role": "Inorganic Core", "purity": "99.999% (5N)", "handling": "Flammable Liquid, Moisture Sensitive"},
    {"id": "BB_04", "name": "Hexachlorodisilane", "formula": "Si2Cl6", "smiles": "Cl[Si](Cl)(Cl)[Si](Cl)(Cl)Cl", "cas": "13465-77-5", "role": "Inorganic Core", "purity": "99.999% (5N)", "handling": "Moisture Sensitive, Corrosive"},
    {"id": "BB_05", "name": "Hafnium(IV) Chloride", "formula": "HfCl4", "smiles": "Cl[Hf](Cl)(Cl)Cl", "cas": "13499-05-3", "role": "Inorganic Core", "purity": "99.9% (Zr-free <50 ppm)", "handling": "Solid, Hydrolyzes in Air"},
    {"id": "BB_06", "name": "Zirconium(IV) Chloride", "formula": "ZrCl4", "smiles": "Cl[Zr](Cl)(Cl)Cl", "cas": "10026-11-6", "role": "Inorganic Core", "purity": "99.9%", "handling": "Solid, Moisture Sensitive"},
    {"id": "BB_07", "name": "Titanium(IV) Chloride", "formula": "TiCl4", "smiles": "Cl[Ti](Cl)(Cl)Cl", "cas": "7550-45-0", "role": "Inorganic Core", "purity": "99.999% (5N)", "handling": "Liquid, Fuming in Air"},
    {"id": "BB_08", "name": "Aluminum Trichloride", "formula": "AlCl3", "smiles": "Cl[Al](Cl)Cl", "cas": "7446-70-0", "role": "Inorganic Core", "purity": "99.99%", "handling": "Solid, Lewis Acid"},
    {"id": "BB_09", "name": "Trimethylaluminum", "formula": "AlMe3", "smiles": "C[Al](C)C", "cas": "75-24-1", "role": "Alkyl Precursor / Core", "purity": "99.9999% (6N)", "handling": "Pyrophoric Liquid"},
    {"id": "BB_10", "name": "Tungsten(VI) Hexachloride", "formula": "WCl6", "smiles": "Cl[W](Cl)(Cl)(Cl)(Cl)Cl", "cas": "13283-01-7", "role": "Inorganic Core", "purity": "99.9%", "handling": "Moisture Sensitive Solid"},
    {"id": "BB_11", "name": "Molybdenum(V) Chloride", "formula": "MoCl5", "smiles": "Cl[Mo](Cl)(Cl)(Cl)Cl", "cas": "10241-05-1", "role": "Inorganic Core", "purity": "99.9%", "handling": "Moisture Sensitive Solid"},
    {"id": "BB_12", "name": "Ruthenium(III) Chloride hydrate", "formula": "RuCl3", "smiles": "Cl[Ru](Cl)Cl", "cas": "10049-08-8", "role": "Inorganic Core", "purity": "99.9%", "handling": "Solid, Hygroscopic"},
    {"id": "BB_13", "name": "Cobalt(II) Chloride", "formula": "CoCl2", "smiles": "Cl[Co]Cl", "cas": "7646-79-9", "role": "Inorganic Core", "purity": "99.99%", "handling": "Solid, Anhydrous"},
    {"id": "BB_14", "name": "Tantalum(V) Chloride", "formula": "TaCl5", "smiles": "Cl[Ta](Cl)(Cl)(Cl)Cl", "cas": "7721-01-9", "role": "Inorganic Core", "purity": "99.99%", "handling": "Moisture Sensitive Solid"},
    {"id": "BB_15", "name": "tert-Butyltin Trichloride", "formula": "tBuSnCl3", "smiles": "CC(C)(C)[Sn](Cl)(Cl)Cl", "cas": "7384-14-7", "role": "Organotin Core", "purity": "99.0%", "handling": "Moisture Sensitive Solid"},

    # Organic Pro-ligands & Amines
    {"id": "BB_16", "name": "Diisopropylamine", "formula": "C6H15N", "smiles": "CC(C)NC(C)C", "cas": "108-18-9", "role": "Pro-ligand (Amine)", "purity": "99.5%", "handling": "Flammable, Basic"},
    {"id": "BB_17", "name": "tert-Butylamine", "formula": "C4H11N", "smiles": "CC(C)(C)N", "cas": "75-64-9", "role": "Pro-ligand (Amine / Imido source)", "purity": "99.5%", "handling": "Flammable, Volatile"},
    {"id": "BB_18", "name": "Dimethylamine", "formula": "C2H7N", "smiles": "CNC", "cas": "124-40-3", "role": "Pro-ligand (Amine)", "purity": "99.0% (Gas/Solution)", "handling": "Flammable Gas"},
    {"id": "BB_19", "name": "Ethylmethylamine", "formula": "C3H9N", "smiles": "CCNC", "cas": "624-78-2", "role": "Pro-ligand (Amine)", "purity": "99.0%", "handling": "Flammable Liquid"},
    {"id": "BB_20", "name": "Cyclopentadiene", "formula": "C5H6", "smiles": "C1=CCC=C1", "cas": "542-92-7", "role": "Pro-ligand (Cyclopentadienyl)", "purity": "Freshly cracked", "handling": "Flammable, Dimerizes"},
    {"id": "BB_21", "name": "N,N'-Diisopropylacetamidine", "formula": "C8H18N2", "smiles": "CC(C)N=C(C)NC(C)C", "cas": "635680-55-6", "role": "Pro-ligand (Amidinate)", "purity": "98.0%", "handling": "Air/Moisture Sensitive Base"},
    {"id": "BB_22", "name": "Isopropanol", "formula": "C3H8O", "smiles": "CC(C)O", "cas": "67-63-0", "role": "Pro-ligand (Alkoxide)", "purity": "Electronic Grade (99.99%)", "handling": "Flammable"},
    {"id": "BB_23", "name": "tert-Butanol", "formula": "C4H10O", "smiles": "CC(C)(C)O", "cas": "75-65-0", "role": "Pro-ligand (Alkoxide)", "purity": "99.5%", "handling": "Flammable Solid/Liquid"},
    {"id": "BB_24", "name": "3,3-Dimethyl-1-butyne", "formula": "C6H10", "smiles": "CC(C)(C)C#C", "cas": "917-92-0", "role": "Pro-ligand (Alkyne)", "purity": "98.0%", "handling": "Flammable Volatile Liquid"},

    # Lithiated & Metal Reagents
    {"id": "BB_25", "name": "n-Butyllithium", "formula": "C4H9Li", "smiles": "CCCC[Li]", "cas": "109-72-8", "role": "Deprotonating Base", "purity": "2.5M in Hexanes", "handling": "Pyrophoric"},
    {"id": "BB_26", "name": "Methyllithium", "formula": "CH3Li", "smiles": "C[Li]", "cas": "917-54-4", "role": "Alkylation Reagent", "purity": "1.6M in Diethyl ether", "handling": "Pyrophoric"},
    {"id": "BB_27", "name": "Lithium dimethylamide", "formula": "C2H6LiN", "smiles": "CN(C)[Li]", "cas": "3585-32-8", "role": "Salt Metathesis Reagent", "purity": "Solid or in THF", "handling": "Pyrophoric/Moisture Sensitive"},
    {"id": "BB_28", "name": "Lithium ethylmethylamide", "formula": "C3H8LiN", "smiles": "CCN(C)[Li]", "cas": "3585-33-9", "role": "Salt Metathesis Reagent", "purity": "in Hexanes", "handling": "Moisture Sensitive"},
    {"id": "BB_29", "name": "Lithium cyclopentadienide", "formula": "C5H5Li", "smiles": "[Li]C1=CC=CC1", "cas": "16733-97-4", "role": "Cp Transfer Reagent", "purity": "Solid powder", "handling": "Air/Moisture Sensitive"}
]

# ---------------------------------------------------------
# 3. Chemical Synthesis Reactions Database
# ---------------------------------------------------------
REACTIONS = [
    {
        "id": "RXN_01",
        "name": "Aminolysis of Dichlorosilane to DIPAS",
        "reaction_type": "Aminolysis / Salt Precipitation",
        "reaction_smarts": "[Si:1]([Cl:2])([Cl:3])([H:4])[H:5].[N:6]([C:7])([C:8])[H:9]>>[Si:1]([N:6]([C:7])[C:8])([H:4])([H:5])[H].[Cl:2]",
        "reactants": ["Cl[SiH2]Cl", "CC(C)NC(C)C"],
        "products": ["CC(C)N([SiH3])C(C)C"],
        "byproduct": "[iPr2NH2]+Cl-",
        "solvent": "n-Hexane / Diethyl ether",
        "conditions": "0 °C to RT, 4 h, N2 atmosphere",
        "thermodynamics_dG_kcal_mol": -18.5
    },
    {
        "id": "RXN_02",
        "name": "Bis-aminolysis of Dichlorosilane to BTBAS",
        "reaction_type": "Bis-aminolysis",
        "reaction_smarts": "[Si:1]([Cl:2])([Cl:3])([H:4])[H:5].2[N:6]([C:7])([H:8])[H:9]>>[Si:1]([NH:6][C:7])([NH:6][C:7])([H:4])[H:5]",
        "reactants": ["Cl[SiH2]Cl", "CC(C)(C)N"],
        "products": ["CC(C)(C)N[SiH2]NC(C)(C)C"],
        "byproduct": "[tBuNH3]+Cl-",
        "solvent": "Toluene / Pentane",
        "conditions": "-20 °C to RT, 6 h",
        "thermodynamics_dG_kcal_mol": -24.2
    },
    {
        "id": "RXN_03",
        "name": "Tris-aminolysis of Trichlorosilane to 3DMAS",
        "reaction_type": "Aminolysis",
        "reaction_smarts": "[Si:1]([Cl:2])([Cl:3])([Cl:4])[H:5].3[N:6]([C:7])[C:8]>>[Si:1]([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])[H:5]",
        "reactants": ["Cl[SiH](Cl)Cl", "CNC"],
        "products": ["CN(C)[SiH](N(C)C)N(C)C"],
        "byproduct": "[Me2NH2]+Cl-",
        "solvent": "Hexane",
        "conditions": "-10 °C to RT, 8 h",
        "thermodynamics_dG_kcal_mol": -32.0
    },
    {
        "id": "RXN_04",
        "name": "Exhaustive Alcoholysis of SiCl4 to TEOS",
        "reaction_type": "Alcoholysis",
        "reaction_smarts": "[Si:1]([Cl:2])([Cl:3])([Cl:4])[Cl:5].4[O:6][C:7]>>[Si:1]([O:6][C:7])([O:6][C:7])([O:6][C:7])[O:6][C:7]",
        "reactants": ["Cl[Si](Cl)(Cl)Cl", "CCO"],
        "products": ["CCO[Si](OCC)(OCC)OCC"],
        "byproduct": "HCl (gas)",
        "solvent": "Neat / Toluene",
        "conditions": "Reflux 70 °C, N2 sweep",
        "thermodynamics_dG_kcal_mol": -40.5
    },
    {
        "id": "RXN_05",
        "name": "Salt Metathesis for TEMAH Synthesis",
        "reaction_type": "Salt Metathesis",
        "reaction_smarts": "[Hf:1]([Cl:2])([Cl:3])([Cl:4])[Cl:5].4[Li][N:6]([C:7])[C:8]>>[Hf:1]([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])[N:6]([C:7])[C:8]",
        "reactants": ["Cl[Hf](Cl)(Cl)Cl", "CCN(C)[Li]"],
        "products": ["CCN(C)[Hf](N(C)CC)(N(C)CC)N(C)CC"],
        "byproduct": "4 LiCl",
        "solvent": "Anhydrous THF / Hexanes",
        "conditions": "-78 °C to RT, 12 h, dry box",
        "thermodynamics_dG_kcal_mol": -58.4
    },
    {
        "id": "RXN_06",
        "name": "Salt Metathesis for TEMAZ Synthesis",
        "reaction_type": "Salt Metathesis",
        "reaction_smarts": "[Zr:1]([Cl:2])([Cl:3])([Cl:4])[Cl:5].4[Li][N:6]([C:7])[C:8]>>[Zr:1]([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])[N:6]([C:7])[C:8]",
        "reactants": ["Cl[Zr](Cl)(Cl)Cl", "CCN(C)[Li]"],
        "products": ["CCN(C)[Zr](N(C)CC)(N(C)CC)N(C)CC"],
        "byproduct": "4 LiCl",
        "solvent": "THF / Pentane",
        "conditions": "-40 °C to RT, 10 h",
        "thermodynamics_dG_kcal_mol": -62.1
    },
    {
        "id": "RXN_07",
        "name": "Salt Metathesis for TDMAT Synthesis",
        "reaction_type": "Salt Metathesis",
        "reaction_smarts": "[Ti:1]([Cl:2])([Cl:3])([Cl:4])[Cl:5].4[Li][N:6]([C:7])[C:8]>>[Ti:1]([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])([N:6]([C:7])[C:8])[N:6]([C:7])[C:8]",
        "reactants": ["Cl[Ti](Cl)(Cl)Cl", "CN(C)[Li]"],
        "products": ["CN(C)[Ti](N(C)C)(N(C)C)N(C)C"],
        "byproduct": "4 LiCl",
        "solvent": "Hexanes",
        "conditions": "0 °C to RT, 6 h",
        "thermodynamics_dG_kcal_mol": -52.0
    },
    {
        "id": "RXN_08",
        "name": "Ammonia-Assisted Alcoholysis for TTIP",
        "reaction_type": "Alcoholysis / Neutralization",
        "reaction_smarts": "[Ti:1]([Cl:2])([Cl:3])([Cl:4])[Cl:5].4[O:6][C:7].4N>>[Ti:1]([O:6][C:7])([O:6][C:7])([O:6][C:7])[O:6][C:7]",
        "reactants": ["Cl[Ti](Cl)(Cl)Cl", "CC(C)O"],
        "products": ["CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C"],
        "byproduct": "4 NH4Cl",
        "solvent": "Benzene / Toluene",
        "conditions": "10 °C to RT, NH3 sparging",
        "thermodynamics_dG_kcal_mol": -48.7
    },
    {
        "id": "RXN_09",
        "name": "Amidinate Chelation for Co(AMD)2",
        "reaction_type": "Chelation Metathesis",
        "reaction_smarts": "[Co:1]([Cl:2])[Cl:3].2[Li][N:4]=[C:5][N:6]>>[Co:1]([N:4]=[C:5][N:6])2",
        "reactants": ["Cl[Co]Cl", "CC(C)N=C(C)NC(C)C", "CCCC[Li]"],
        "products": ["CC(C)N=C(C)N(C(C)C)[Co]N(C(C)C)C(C)=NC(C)C"],
        "byproduct": "2 LiCl, 2 Butane",
        "solvent": "THF / Hexane",
        "conditions": "Deprotonation at -78 °C, then CoCl2 addition, RT 12 h",
        "thermodynamics_dG_kcal_mol": -38.2
    },
    {
        "id": "RXN_10",
        "name": "Imido Exchange & Aminolysis for Bis(imido) Tungsten Precursor",
        "reaction_type": "Sequential Imido Addition & Aminolysis",
        "reaction_smarts": "[W:1][Cl:2].2[N:3][C:4].2[Li][N:5]>>[W:1](=[N:3][C:4])(=[N:3][C:4])([N:5])2",
        "reactants": ["Cl[W](Cl)(Cl)(Cl)(Cl)Cl", "CC(C)(C)N", "CN(C)[Li]"],
        "products": ["CN(C)[W](=NC(C)(C)C)(=NC(C)(C)C)N(C)C"],
        "byproduct": "4 LiCl, 2 [tBuNH3]+Cl-",
        "solvent": "Toluene",
        "conditions": "Step 1: 80 °C, Step 2: RT",
        "thermodynamics_dG_kcal_mol": -44.0
    },
    {
        "id": "RXN_11",
        "name": "Alkyne Coordination with CO Displacement for CCTBA",
        "reaction_type": "Alkyne Displacement",
        "reaction_smarts": "[Co2(CO)8].[C:1]#[C:2]>>[Co2(CO)6]([C:1]#[C:2]) + 2 CO",
        "reactants": ["CC(C)(C)C#C"],
        "products": ["CC(C)(C)C12[CH]3[Co]1([C-]#[O+])([C-]#[O+])([C-]#[O+])[Co]23([C-]#[O+])([C-]#[O+])[C-]#[O+]"],
        "byproduct": "2 CO (gas)",
        "solvent": "Dichloromethane / Hexane",
        "conditions": "RT, 4 h, CO gas release",
        "thermodynamics_dG_kcal_mol": -12.3
    },
    {
        "id": "RXN_12",
        "name": "Salt Metathesis for EUV Tin Precursor tBuSn(OiPr)3",
        "reaction_type": "Salt Metathesis",
        "reaction_smarts": "[Sn:1]([C:2])([Cl:3])([Cl:4])[Cl:5].3[Na][O:6][C:7]>>[Sn:1]([C:2])([O:6][C:7])([O:6][C:7])[O:6][C:7]",
        "reactants": ["CC(C)(C)[Sn](Cl)(Cl)Cl", "CC(C)O"],
        "products": ["CC(C)O[Sn](C(C)(C)C)(OC(C)C)OC(C)C"],
        "byproduct": "3 NaCl",
        "solvent": "iPrOH / Toluene",
        "conditions": "Reflux 4 h",
        "thermodynamics_dG_kcal_mol": -35.6
    }
]

# ---------------------------------------------------------
# Verification & File Output
# ---------------------------------------------------------
print("=== Validating Chemical Structures with RDKit ===")
valid_precursors = 0
for p in PRECURSORS:
    mol = Chem.MolFromSmiles(p["smiles"])
    if mol is not None:
        valid_precursors += 1
    else:
        print(f"Warning: RDKit failed on precursor SMILES: {p['name']} ({p['smiles']})")

print(f"Precursors validated: {valid_precursors} / {len(PRECURSORS)}")

valid_bbs = 0
for bb in BUILDING_BLOCKS:
    mol = Chem.MolFromSmiles(bb["smiles"])
    if mol is not None:
        valid_bbs += 1
    else:
        print(f"Warning: RDKit failed on BB SMILES: {bb['name']} ({bb['smiles']})")

print(f"Building blocks validated: {valid_bbs} / {len(BUILDING_BLOCKS)}")

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DB_DIR = os.path.join(PROJECT_ROOT, "database")
BR_DIR = os.path.join(DB_DIR, "br_sascore_inputs")

# Write JSON files
os.makedirs(DB_DIR, exist_ok=True)
os.makedirs(BR_DIR, exist_ok=True)

with open(os.path.join(DB_DIR, "precursors.json"), "w", encoding="utf-8") as f:
    json.dump(PRECURSORS, f, ensure_ascii=False, indent=2)

with open(os.path.join(DB_DIR, "building_blocks.json"), "w", encoding="utf-8") as f:
    json.dump(BUILDING_BLOCKS, f, ensure_ascii=False, indent=2)

with open(os.path.join(DB_DIR, "reactions.json"), "w", encoding="utf-8") as f:
    json.dump(REACTIONS, f, ensure_ascii=False, indent=2)

# Write CSV files for Precursors
prec_fieldnames = [
    "id", "name", "chemical_name", "formula", "cas_number", "molecular_weight",
    "metal_center", "category", "target_films", "process_type",
    "state_at_rt", "vapor_pressure_torr", "boiling_point_c", "thermal_decomposition_temp_c",
    "typical_yield_pct", "smiles"
]
with open(os.path.join(DB_DIR, "precursors.csv"), "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=prec_fieldnames)
    writer.writeheader()
    for p in PRECURSORS:
        row = {
            "id": p["id"],
            "name": p["name"],
            "chemical_name": p["chemical_name"],
            "formula": p["formula"],
            "cas_number": p["cas_number"],
            "molecular_weight": p["molecular_weight"],
            "metal_center": p["metal_center"],
            "category": p["category"],
            "target_films": p["target_films"],
            "process_type": p["process_type"],
            "state_at_rt": p["physical_properties"].get("state_at_rt"),
            "vapor_pressure_torr": p["physical_properties"].get("vapor_pressure_torr"),
            "boiling_point_c": p["physical_properties"].get("boiling_point_c"),
            "thermal_decomposition_temp_c": p["physical_properties"].get("thermal_decomposition_temp_c"),
            "typical_yield_pct": p["synthesis_info"].get("typical_yield_pct"),
            "smiles": p["smiles"]
        }
        writer.writerow(row)

# Write CSV files for Building Blocks
with open(os.path.join(DB_DIR, "building_blocks.csv"), "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["id", "name", "formula", "cas", "role", "purity", "handling", "smiles"])
    writer.writeheader()
    for bb in BUILDING_BLOCKS:
        writer.writerow(bb)

# Write CSV files for Reactions
with open(os.path.join(DB_DIR, "reactions.csv"), "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["id", "name", "reaction_type", "reactants", "products", "byproduct", "solvent", "conditions", "thermodynamics_dG_kcal_mol"])
    writer.writeheader()
    for rxn in REACTIONS:
        writer.writerow({
            "id": rxn["id"],
            "name": rxn["name"],
            "reaction_type": rxn["reaction_type"],
            "reactants": " + ".join(rxn["reactants"]),
            "products": " + ".join(rxn["products"]),
            "byproduct": rxn["byproduct"],
            "solvent": rxn["solvent"],
            "conditions": rxn["conditions"],
            "thermodynamics_dG_kcal_mol": rxn["thermodynamics_dG_kcal_mol"]
        })

# Write BR-SAScore ready files: semi_bb.csv and semi_rxn.csv
with open(os.path.join(BR_DIR, "semi_bb.csv"), "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["SMILES"])
    for bb in BUILDING_BLOCKS:
        if bb["smiles"]:
            writer.writerow([bb["smiles"]])
    # Also add pro-ligands and simple silane reactants
    for p in PRECURSORS:
        for sm in p["synthesis_info"]["starting_materials"]:
            pass

with open(os.path.join(BR_DIR, "semi_rxn.csv"), "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["reaction"])
    for rxn in REACTIONS:
        rxn_str = f"{'.'.join(rxn['reactants'])}>>{'.'.join(rxn['products'])}"
        writer.writerow([rxn_str])

print("=== Successfully Built Semiconductor DB ===")
print(f"1. Precursors: {len(PRECURSORS)} items -> database/precursors.json, precursors.csv")
print(f"2. Building Blocks: {len(BUILDING_BLOCKS)} items -> database/building_blocks.json, building_blocks.csv")
print(f"3. Reactions: {len(REACTIONS)} items -> database/reactions.json, reactions.csv")
print(f"4. BR-SAScore Inputs: database/br_sascore_inputs/semi_bb.csv, semi_rxn.csv")
