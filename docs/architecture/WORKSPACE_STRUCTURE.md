# Workspace Architecture & Module Organization

This document outlines the directory layout and module responsibilities following code refactoring and topic-based documentation structuring.

```
synth_acc/
├── generate_html.py                          # Primary HTML build script (reads template.html + database/ -> index.html)
├── template.html                             # Dual-tab HTML master template (Summary + 15-section Report)
├── index.html                                # Compiled standalone HTML report
├── synthetic_accessibility_report.html       # Synchronized distribution copy of index.html
├── README.md                                 # Top-level scientific framework, equations, and quickstart
├── report.md                                 # Formal engineering & verification report
├── docs/
│   ├── architecture/
│   │   └── WORKSPACE_STRUCTURE.md            # Directory architecture and execution flow map
│   ├── audits/
│   │   ├── LOCAL_BENCHMARK_ENVIRONMENT_AUDIT.md  # Python runtime, external repos, patches, and 15-model audit
│   │   └── PHASE2_SEVENN_VALIDATION_AUDIT.md # Audit of Phase 1/2 scripts, hardcoding removal, All-sample/CV/LOO review
│   └── literature/
│       └── BIBLIOGRAPHY.md                   # Section-indexed bibliography (34 DOIs) and BibTeX archive
├── scripts/
│   ├── benchmark/                            # Model setup, SA benchmarks, SevenNet features & ML validation
│   │   ├── setup_external_models.py          # Clones pinned upstream repositories into external/
│   │   ├── patch_external_models.py          # Applies NumPy/Python compatibility patches to external/
│   │   ├── run_all_benchmarks.py             # Master 15-model x 27-precursor SA & SMARTS benchmark runner
│   │   ├── audit_reproduce_all.py            # Independent recalculation and counterfactual leakage check (9/9)
│   │   ├── run_phase1_baseline.py            # Evaluates 55 molecules across SAScore, SCScore, SYBA, RAscore, BR-SAScore
│   │   ├── run_phase2_sevenn_features.py     # Computes SevenNet-0 single-point energies, bond stretches, dimer & %V_bur
│   │   └── run_model_tuning_comparison.py    # Un-hardcoded continuous physical screening & All-Sample/5Fold/LOOCV suite
│   ├── data/                                 # Curated DBs, literature downloaders & verification scripts
│   │   ├── build_semiconductor_db.py         # Builds precursors.json, building_blocks.json, reactions.json, CSVs
│   │   ├── fetch_external_precursor_dbs.py   # Extracts AtomicLimits ALD precursor aggregates & aliases
│   │   ├── verify_inhibitor_db.py            # Offline RDKit & online PubChem/Crossref verification of AS-ALD inhibitors
│   │   ├── download_literature.py            # CrossRef/arXiv DOI verifier & open-access PDF downloader
│   │   ├── verify_repos_and_precursors.py    # Audits GitHub stars/licenses and semiconductor literature DOIs
│   │   └── verify_vasp_frameworks.py         # Audits 10 VASP workflow & automation repositories
│   └── patches/                              # Unified diff compatibility patches for external models
├── database/                                 # Ground-truth JSON/CSV datasets and reproducible benchmark artifacts
├── tests/
│   └── test_benchmark_reproducibility.py     # Automated unittest suite (9 formal tests)
├── external/                                 # Pinned local upstream SA/retrosynthesis repositories (gitignored)
└── downloaded_papers/                        # Local PDF/XML open-access literature archive (gitignored)
```

## Execution Pipeline Sequence

1. **Database & Metadata Construction**:
   - `python scripts/data/build_semiconductor_db.py`
   - `python scripts/data/verify_inhibitor_db.py --offline`
   - `python scripts/data/fetch_external_precursor_dbs.py`
2. **External Model Setup & 15-Model Empirical Benchmark**:
   - `python scripts/benchmark/setup_external_models.py`
   - `python scripts/benchmark/run_all_benchmarks.py`
   - `python scripts/benchmark/audit_reproduce_all.py`
3. **Phase 1 & Phase 2 SevenNet Evaluation**:
   - `python scripts/benchmark/run_phase1_baseline.py`
   - `python scripts/benchmark/run_phase2_sevenn_features.py` *(CPU ~21 min for full 55 molecules)*
   - `python scripts/benchmark/run_model_tuning_comparison.py`
4. **Automated Tests & HTML Report Generation**:
   - `python -m unittest discover tests`
   - `python generate_html.py`
