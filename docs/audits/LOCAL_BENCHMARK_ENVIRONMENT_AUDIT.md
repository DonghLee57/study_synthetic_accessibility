# Local Python Environment, Repository Inventory & Empirical Benchmark Audit Trail

> **Auditor Notice**: This document provides a reproducible empirical audit trail detailing the local runtime environment, repository commits, weights/resource artifacts, execution protocols, local git patches, and empirical data files used to evaluate **15 unique open-source SA/CASP architectures (16 benchmark profiles including BR-SAScore reagent-tuned counterfactuals)** across **27 industrial semiconductor thin-film precursors** (Phase 1) and the curated **55-molecule evaluation dataset** (Phase 2). For the Phase 2 SevenNet-0 physical feature and LOOCV validation audit, see [`PHASE2_SEVENN_VALIDATION_AUDIT.md`](./PHASE2_SEVENN_VALIDATION_AUDIT.md).

---

## 1. Executive Metadata & Audit Scope

- **Date of Execution & Audit**: November 6, 2026 (Phase 2 & Path Audit: November 2026)
- **Lead Role**: Scientific Intelligence & Validation Lead
- **Workspace Directory Root**: `C:\Users\user\Downloads\dev_w_antigravity\synth_acc`
- **Audit Target**: Verification of all 15 unique open-source model architectures across the **4-Stage Execution Pipeline**:
  - **Stage 1 (Syntax & Ingestion Layer)**
  - **Stage 2 (Representation & Embedding Layer)**
  - **Stage 3 (Reaction Policy & Coordination Manifold Layer)**
  - **Stage 4 (Physical Synthesizability & Route Output Layer, $N_{\text{pass}}$)**
- **Target Chemical Dataset**: 27 high-purity semiconductor thin-film precursors (Silicon, High-k, Metal, and EUV organometallics) compiled in [`database/precursors.json`](file:///C:/Users/user/Downloads/dev_w_antigravity/synth_acc/database/precursors.json).

---

## 2. Local Python Runtime Environment

All benchmarks and empirical code audits were executed in the following localized Windows 64-bit Python environment:

### 2.1 Interpreter & System Specifications

| Parameter | Local Environment Value |
| :--- | :--- |
| **Operating System** | Microsoft Windows 11 Home (build 10.0.26300, 64-bit, AMD64) |
| **Default Shell** | Windows PowerShell 5.1 / Git Bash |
| **Python Interpreter** | `C:\Users\user\AppData\Local\Microsoft\WindowsApps\PythonSoftwareFoundation.Python.3.11_qbz5n2kfra8p0\python.exe` |
| **Python Version** | `3.11.9 (tags/v3.11.9:de54cf5, Apr 2 2024, 10:12:12) [MSC v.1938 64 bit (AMD64)]` |
| **Working Directory** | `C:\Users\user\Downloads\dev_w_antigravity\synth_acc` |

### 2.2 Core Package & Library Versions (from `pip list`)

| Package Name | Installed Version | Primary Role in Benchmark | Execution Mode |
| :--- | :---: | :--- | :--- |
| **`rdkit`** | `2023.9.6` | Cheminformatics core, SMILES parsing, Morgan fingerprints, built-in SAScore | Wheel (PyPI) - Executed live |
| **`torch`** | `2.4.1+cpu` | Deep learning backend (CPU-only build; `torch.version.cuda` is `None`) | Wheel (PyTorch CPU) |
| **`torch-geometric`**| `2.7.0` | PyG core; note that GASA and DFRscore require DGL, not PyG | Wheel (PyPI) |
| **`xgboost`** | `1.7.6` | Gradient boosted tree engine for RAscore retrosynthetic classifier | Wheel (PyPI) - Executed via Booster |
| **`aizynthfinder`** | `4.4.1` | Monte Carlo Tree Search (MCTS) retrosynthetic planning engine | Wheel (PyPI) - Executed live |
| **`BRSAScore`** | `0.1.1` | Decoupled building block & reaction synthesizability framework | Wheel (PyPI, SNU MICC) - Executed live |
| **`syba`** | `1.0.2a0` | Naive Bayesian synthesizability classifier with count table | Editable (`external/syba`) - Executed live |
| **`rdchiral`** | `1.1.0` | SMARTS reaction template matching engine (Retro*, Coley group) | Wheel (PyPI) - Executed live |
| **`transformers`** | `4.48.3` | Hugging Face Transformers for ReaSyn and DeepSA tokenizers | Wheel (PyPI) - Executed live |
| **`scikit-learn`** | `1.8.0` | Classical ML evaluators, ROC metrics, and featurizers | Wheel (PyPI) |
| **`scipy`** | `1.14.0` | Scientific computing, sparse matrix math, Bayesian distributions | Wheel (PyPI) |
| **`numpy`** | `1.26.4` | Array-based tensor processing and standalone SCScore matrix math | Wheel (PyPI) |
| **`pymatgen`** | `2025.10.7` | Materials and crystal structure processing (VASP synergy) | Wheel (PyPI) |
| **`ase`** | `3.29.0` | Atomic Simulation Environment (DFT calculator interfaces) | Wheel (PyPI) |

### 2.3 Missing Runtime Dependencies & Execution Constraints

The following libraries are not installed in the local environment, which restricts end-to-end model execution for 3 specific models. Their evaluations are derived from static source code analysis, tokenizer execution, and featurizer diagnostics:

| Missing Package | Target Model | Execution Consequence & Verified Failure Mechanism |
| :--- | :--- | :--- |
| **`dgl`**, **`dgllife`** | **GASA**, **DFRscore** | GNN inference cannot run end-to-end. Static audit of `external/gasa/model/gasa_utils.py` reveals that unlisted elements ($Hf, Ta, Ru, W, Mo$) pass through `dgllife.utils.atom_type_one_hot(..., encode_unknown=False)`. Rather than triggering an unhandled exception crash, DGL-Life emits an **all-zero feature vector**, causing **silent mis-encoding and complete representation loss**. |
| **`autogluon`** | **DeepSA** | `DeepSA_predict.py` cannot execute tabular neural prediction. Evaluated via ELECTRA SMILES subword tokenizer and Lipinski Rule of 5 analysis. |
| **`tensorflow`** | **RAscore (DNN variant)** | The DNN variant (`RAscore_NN.py`) cannot run. The gradient boosted variant (`RAscore_XGB.py`) executes live using the direct XGBoost Booster API. |

---

## 3. Cloned Repositories, Hash Pinning & Local Compatibility Patches

All 12 open-source model repositories are located in `external/`. Because `external/` contains large binary weight checkpoints and count tables (>150 MB), it is excluded from Git tracking via `.gitignore`. 

To clone and configure all repositories deterministically in a single step, run:
```bash
python scripts/benchmark/setup_external_models.py
```

### 3.1 Repository Inventory & Commit Manifest

| Directory | Model Name & Venue | Git Remote URL | Branch | Commit Hash (HEAD) | Critical Weight / Resource Files |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **`external/scscore`** | **SCScore** (MIT, 2018) | `https://github.com/connorcoley/scscore.git` | `master` | `37090a6aa8220408f572de20224c33c221312d16` | `models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz` (NumPy standalone runner; evaluates boolean Morgan fingerprints) |
| **`external/syba`** | **SYBA** (UCT Prague, 2020) | `https://github.com/lich-uct/syba.git` | `master` | `30929eb2d3b9a03bef326a1f5180f26f55b19a29` | `syba/resources/syba.csv.gz` (~25 MB count table; **methodological deviation** from upstream default `syba4.csv.gz` [~115 MB] for 10x faster execution and consistency) |
| **`external/rascore`** | **RAscore** (Reymond, 2021) | `https://github.com/reymond-group/RAscore.git` | `master` | `cb77db503ee5cbf0e8bb8963df6e5b76b3a94f06` | `RAscore/models/XGB_chembl_ecfp_counts/model.pkl` |
| **`external/rxnflow`** | **RxnFlow** (KAIST, ICLR 2025) | `https://github.com/SeonghwanSeo/RxnFlow.git` | `main` | `530d098b56625381395fe7a30a2a10b7aacb8077` | `data/templates/explore/reaction.yaml` (30 rules), `data/templates/cgflow/reaction.yaml` (42 rules) |
| **`external/reasyn`** | **BioNeMo ReaSyn** (NVIDIA, ICLR 2026 Oral) | `https://github.com/NVIDIA-BioNeMo/ReaSyn.git` | `reasyn_v2` | `f28270ca3a5873b2bf7eb52821c5202b111ef9f4` | `reasyn/chem/featurize.py` (153-token SMILES vocab covering 118 elements; `NUM_REACTIONS = 115`) |
| **`external/synnet`** | **SynNet** (NeurIPS 2021) | `https://github.com/wenhao-gao/SynNet.git` | `master` | `3ccf4fcd6eb0605cb60f6376ad7b4a7accaf1549` | `data/assets/reaction-templates/hb.txt` (HB 91 rules) |
| **`external/chemprojector`** | **ChemProjector** (ICLR 2024) | `https://github.com/luost26/ChemProjector.git` | `main` | `dd3afda12d33437f9c41bc7b246903628e60429a` | `data/reaction_templates_hb.txt` (HB 91 rules), `chemprojector/chem/smiles_vocab.txt` |
| **`external/synformer`** | **SynFormer** (2023) | `https://github.com/wenhao-gao/synformer.git` | `main` | `bef02bc4481adf3ed792090fc6169331d4d0553a` | `data/rxn_templates/comprehensive.txt` (115 rules) |
| **`external/retro_star`** | **Retro\*** (ICML 2020) | `https://github.com/binghong-ml/retro_star.git` | `master` | `60e7a24213810706a7075d113704c9320158d726` | `retro_star/packages/mlp_retrosyn/`, `rdchiral` USPTO 50k rules |
| **`external/dfrscore`** | **DFRscore** (2023) | `https://github.com/Hwoo-Kim/DFRscore.git` | `main` | `5b1a4a4157105089bb3c818fd4a136fe255c80b8` | `scripts/modelScripts/model.py`, `preprocessing.py` |
| **`external/deepsa`** | **DeepSA** (2023) | `https://github.com/Shihang-Wang-58/DeepSA.git` | `main` | `eaaa6a52e7edd82fb4ec3823e45c6f6f373e8ef2` | `DeepSA_predict.py` |
| **`external/gasa`** | **GASA** (2022) | `https://github.com/cadd-synthetic/GASA.git` | `master` | `b4f46553da749c4dc474a7fe1fc6041770d44659` | `model/gasa_utils.py` (`allowable_set`) |

### 3.2 Local Compatibility Patches (Documented Diffs)

A fresh clone from GitHub requires compatibility patches for Python 3.11 and NumPy 1.26. These patches are automated via `python scripts/benchmark/patch_external_models.py` and saved in `scripts/patches/`:

#### 1. `external/scscore` Patch (`scripts/patches/scscore_numpy_compat.patch`)
- **Issue**: `sigmoid(x)` used `math.exp(-x)` which fails when `x` is a NumPy array/scalar. Furthermore, NumPy 1.24+ removed `np.bool`, and Python 3 removed `cPickle`.
- **Diff**:
```diff
--- a/scscore/standalone_model_numpy.py
+++ b/scscore/standalone_model_numpy.py
@@ -23,7 +23,7 @@ FP_len = 1024
 FP_rad = 2
 
 def sigmoid(x):
-  return 1 / (1 + math.exp(-x))
+  return 1.0 / (1.0 + math.exp(-float(np.squeeze(x))))
 
 class SCScorer():
@@ -53,3 +53,3 @@
-                    useChirality=True), dtype=np.bool)
+                    useChirality=True), dtype=bool)
@@ -98,3 +98,3 @@
-            import cPickle as pickle
+            import pickle
```

#### 2. `external/syba` Patch (`scripts/patches/syba_resources.patch`)
- **Issue**: `syba4.csv.gz` (~115 MB compressed, ~500 MB raw text) requires >40 seconds per import. Switching default table to `syba.csv.gz` (~25 MB) loads in ~4 seconds and reproduces the verified benchmark score matrix.
- **Diff**:
```diff
--- a/syba/syba.py
+++ b/syba/syba.py
@@ -50,3 +50,3 @@
-        with gzip.open(os.path.join(this_dir, "resources", "syba4.csv.gz"), mode="rt") as counts:
+        with gzip.open(os.path.join(this_dir, "resources", "syba.csv.gz"), mode="rt") as counts:
             self.fitFromCountFile(counts)
```

---

## 4. Input Datasets & Benchmark Artifacts

The benchmark executes over version-controlled datasets located in `database/`:

1. **Precursors Primary Database**: [`database/precursors.json`](file:///C:/Users/user/Downloads/dev_w_antigravity\synth_acc/database/precursors.json) (27 industrial precursors).
2. **BR-SAScore Domain Dictionaries**:
   - Building Blocks: [`database/br_sascore_inputs/semi_bb.csv`](file:///C:/Users/user/Downloads/dev_w_antigravity\synth_acc/database/br_sascore_inputs/semi_bb.csv) (Gelest chlorosilanes, metal chlorides $HfCl_4, TiCl_4, ZrCl_4, AlCl_3$, amines, alkoxides, alkyl lithiums).
   - Reaction Rules: [`database/br_sascore_inputs/semi_rxn.csv`](file:///C:/Users/user/Downloads/dev_w_antigravity\synth_acc/database/br_sascore_inputs/semi_rxn.csv) (Salt metathesis, transamination, alcoholysis).
3. **Master Benchmark Execution Script**:
   - [`scripts/benchmark/run_all_benchmarks.py`](file:///C:/Users/user/Downloads/dev_w_antigravity/synth_acc/scripts/benchmark/run_all_benchmarks.py): Regenerates all Phase 1 benchmark databases from scratch.
4. **Empirical Evaluation Outputs**:
   - [`database/precursor_real_benchmark_results.json`](file:///C:/Users/user/Downloads/dev_w_antigravity/synth_acc/database/precursor_real_benchmark_results.json): Per-precursor scores across all benchmarked models.
   - [`database/expanded_models_benchmark_results.json`](file:///C:/Users/user/Downloads/dev_w_antigravity/synth_acc/database/expanded_models_benchmark_results.json): Diagnostic profiles across all 15 unique models (16 profiles).
   - [`database/benchmark_summary_stats.json`](file:///C:/Users/user/Downloads/dev_w_antigravity/synth_acc/database/benchmark_summary_stats.json): Aggregated summary counts and pass rates.
   - [`database/benchmark_phase2_tuned_comparison.json`](file:///C:/Users/user/Downloads/dev_w_antigravity/synth_acc/database/benchmark_phase2_tuned_comparison.json): Phase 2 un-hardcoded uniform screening and ML cross-validation (All-Sample, 5-Fold CV, LOOCV) across 55 molecules.

---

## 5. Verification Harness & Execution Protocol

Two automated testing mechanisms are committed to the repository for independent verification:

### 5.1 Formal Python Unit Test Suite (`unittest`)

Run the complete test suite across all recalculations:
```bash
python -m unittest tests.test_benchmark_reproducibility
# OR
python -m unittest discover tests
```

**Output**:
```text
Ran 9 tests in ~6.0s
OK
```

### 5.2 Standalone Audit Script (`scripts/benchmark/audit_reproduce_all.py`)

Run the detailed audit script to verify live recalculations for SAScore, SCScore, SYBA, RAscore, BR-SAScore (Default & Custom), ReaSyn lossless roundtrip, GASA silent mis-encoding, and DFRscore/ChemProjector manifolds:
```bash
python scripts/benchmark/audit_reproduce_all.py
```

---

## 6. Detailed 4-Stage Execution Diagnostics Summary (16 Models)

| Model | Stage 1 (Syntax Ingestion) | Stage 2 (Representation) | Stage 3 (Reaction Policy) | Stage 4 (True $N_{\text{pass}}$) | Exact Breakdown Barrier & Physical Root Cause |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **BR-SAScore (Strict Unleaked)** | **100% (27/27)** | **`max(existing, -0.47)`** | **11.1% (3/27)** | **11.1% (3/27)** | **Stage 3 Barrier (0 Worsened under `max(existing, -0.47)`)**: Because `BRSAScore.calculateScore` evaluates only radius-2 fragments (`vs[0][1] == 2`) and averages **only negative-scoring fragments (`fscore < 0`)**—where official `BScores` span `-3.0` to `+3.0` (median `-2.54`, 65 positive entries, 0 exact `0.0` keys)—updating reagent fragments via `max(existing, -0.47)` keeps them inside the negative-only average so **0 precursors score worse than Default** (e.g., `Cp(NMe2)3Zr` `7.20 -> 6.95`, `Co(AMD)2` `8.29 -> 7.91`), while the pass rate remains **3/27** (TSA, TMA, DMAH) across `-0.47`, `-0.1`, and `-0.01` because $M\text{-N}/M\text{-C}$ center fragments remain `-6.0`. (By contrast, naive `0.0` overwriting drops organic fragments from `nf` and worsens 3 Cp/amidinate precursors.) |
| **BR-SAScore (Leaked Baseline)** | **100% (27/27)** | **Memorized ($s=0$)** | **100.0% (27/27)** | **100.0% (27/27)** | **Counterfactual Audit Baseline**: Illustrates the circular data leakage identified by Claude audit where pre-registering all 27 target products into the dictionary zeroes out even the metal-coordination fragments (`score1 = 0`), yielding artificial 27/27 passes (24/27 [88.9%] in earlier commit `fb7305a` when 24 products were listed). |
| **BR-SAScore (Default)**| **100% (27/27)** | **100% Valid** | **11.1% (3/27)** | **11.1% (3/27)** | **Stage 3 Breakdown**: eMolecules/USPTO dictionaries lack Si-N/M-N, assigning max $-6.0$ penalty to 24 precursors ($3.82 \sim 10.00$ pts). |
| **SAScore (RDKit)** | **100% (27/27)** | **100% Morgan** | N/A (Scoring) | **7.4% (2/27)** | **Stage 2·4 Breakdown**: Negative log-frequency penalties for near-zero PubChem drug counts cause severe False Negatives (TSA 6.91, DMAH 6.31, connected Ru(EtCp)2 4.70); only TEOS (2.84) and 4DMAS (3.37) pass. |
| **SCScore (Reaxys)** | **100% (27/27)** | **100% 1024-bit** | N/A (DNN) | **77.8% (21/27)** | **Stage 4 Bias**: Early-stage reactant registration in Reaxys (`full_reaxys_model_1024bool`) yields low complexity scores (`1.01~1.89` for 21/27), not reflecting pyrophoric hazards (TMA). |
| **SYBA (Bayesian)** | **100% (27/27)** | **100% Subtrees** | N/A (Bayesian) | **85.2% (23/27)** | **Stage 2·4 Misclassification**: ZINC drug library (`syba.csv.gz`, method deviation noted) lacks metal-Cp fragments; metallocene/Cp complexes receive strongly negative scores (`CpTi`/`CpZr` $-24.98$, `ZMA` $-18.21$, `Ru(EtCp)2` $-17.31$). |
| **RAscore (AiZynth)** | **100% (27/27)** | **Severe Distortion** | N/A (Classifier) | **0.0% Valid (100% FP)** | **Stage 2 Breakdown**: ECFP6 projects metal centers as generic alkyl/amine carbons. 22/27 score ≥ 0.9 (mean `0.908`, range `0.602~0.988`); exceptions: CCTBA 0.602, BTBMMo 0.629, Co(AMD)2 0.639, BTBMW 0.753, TBTDET 0.875. |
| **RxnFlow (ICLR 2025)** | **100% (27/27)** | **100% Synthons** | **33.3% Substruct / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Breakdown**: Live RDKit SMARTS matching of 72 Explore (30) / CGFlow (42) reverse templates matches 9/27 generic alkyl-amine product substructures (`reductive_amination_ketone`, `reductive_amination_aliphatic_aldehyde`, `amine_halide_coupling`) and 0 metal rules in set. |
| **BioNeMo ReaSyn** | **100% (118 elem)** | **100% Sequence** | **88.9% Substruct (85.2% ≥2-atom) / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Breakdown**: Universal 118-element tokenizer ingests all atoms (100% lossless roundtrip, 0 UNK), while SynFormer 115 patent templates contain 0 metal rules and match 24/27 generic product substructures (23/27 when excluding 1-atom `[N:1]` template `#84` that trivially matches carbon-free `TSA`). |
| **AiZynthFinder** | **100% (27/27)** | **100% `Molecule`** | **Unexecuted (`null`)** | **N/A (`null`)** | **Stage 3 Static Analysis**: Live `aizynthfinder.chem.Molecule(smiles=smi).sanitize()` succeeds for 27/27; MCTS policy/search not executed locally (weights not downloaded; `policy_proposals`, `search_depth`, `solved` set to `null`). |
| **ChemProjector** | **100% (Atoms&le;100)** | **100% Latent** | **85.2% Substruct / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Breakdown**: Live RDKit SMARTS matching of Hartenfeller-Button (HB) 91 templates matches 23/27 generic product substructures (all except TSA, HCDS, TMA, DMAH; includes 2-atom patterns `[N:1][#6:2]`, `[O:1][C:2]`), with 0 metal rules in set. |
| **SynNet** | **100% (27/27)** | **100% Morgan (`>0` on)** | **85.2% Substruct / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Breakdown**: Nonzero 2048-bit Morgan fingerprint verified live (`GetNumOnBits() > 0` for 27/27). Shares the exact HB 91 reaction template file with ChemProjector (23/27 generic product substructure matches, 0 metal rules in set). |
| **SynFormer** | **100% (27/27)** | **100% Native Vocab** | **88.9% (&ge;2: 85.2%, &ge;4: 40.7%) / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Breakdown**: Native `synformer.chem.featurize.tokenize_smiles` passes 27/27. Live RDKit `AllChem.ReactionFromSmarts` product matching across 115 `comprehensive.txt` templates matches 24/27 all-pattern (including carbon-free `TSA` via 1-atom `[N:1]` #84), 23/27 for $\ge 2$-atom patterns (12 of which match only 2–3 atom patterns like `N-C`), and 11/27 for $\ge 4$-atom patterns, with 0 metal rules in set. |
| **Retro\*** | **100% (27/27)** | **100% `rdchiral`** | **Unexecuted (`null`)** | **N/A (`null`)** | **Stage 3 Static Analysis**: Live `rdchiralReactants(smi)` initialization succeeds for 27/27; Neural A* search not executed locally (USPTO templates/weights not downloaded; `uspto_rules_matched` and `solved` set to `null`). |
| **DFRscore** | **100% (27/27)** | **Severe Collapse** | N/A (GAT Reg) | **0.0% Grounded** | **Stage 2 Breakdown**: All metals (Si, Ti, Zr, Hf, Ru, W, Mo, Al, Co, Ta, Zn, Sn) collapse into a single "ELSE" token, erasing chemical identity. |
| **DeepSA** | **100% (27/27)** | **Severe Mismatch** | N/A (Multimodal) | **0.0% Valid** | **Stage 2·4 Breakdown**: Blocked by missing `autogluon`. Theoretical analysis shows Lipinski Rule of 5 and ELECTRA subword splitting ([Zr], [Hf] fragmented) distort scores. |
| **GASA** | **25.9% Unsupported (7/27)** | **74.1% Distorted** | N/A (GAT Class) | **0.0% Valid** | **Stage 1 Unsupported Elements & Silent Mis-encoding**: Hf, Ta, Ru, W and Mo are absent from the allowable atom set (7 precursors). GASA was not executed end-to-end (DGL/DGL-Life missing). Its featurizer calls `atom_type_one_hot(..., encode_unknown=False)`, which in DGL-Life returns an **all-zero element vector** for an unknown atom instead of raising an error. The actual failure mode is **silent mis-encoding / feature zeroing**, not an immediate runtime crash. |

### 6.1 ReaSyn Tokenizer Roundtrip Analysis & Offset Resolution

In previous benchmark logs, ReaSyn entries were flagged as `roundtrip_exact: false`. An in-depth code audit of `external/reasyn/reasyn/chem/featurize.py` discovered that this was caused by an indexing offset within ReaSyn's decoding utility:
- `_smiles_token_to_id` indexes tokens starting at `1`.
- `_smiles_id_to_token` indexes tokens starting at `TokenType.MOL_END + 1` (which equals `4`).
- Calling `decode_smiles(tokenize_smiles(smi))` therefore produces an offset character sequence (e.g. `C` [id 6] is looked up at id 6 in `_smiles_id_to_token`, which maps to `Li`).

When decoded using the consistent 1-based token mapping (`_smiles_vocab`), **all 27 precursors achieve 100% lossless character-for-character reconstruction with 0 unknown tokens**. The tokenizer itself possesses full semantic ingestion capacity for all 118 elements; the breakdown occurs strictly at **Stage 3** (0 transition metal reaction templates in the 115 patent reaction space).

---

## 7. Archival & Audit Sign-Off

- **Verification Harness Status**: 100% of live recalculations and diagnostic assertions verified.
- **Formal Unit Test Suite**: `tests/test_benchmark_reproducibility.py` (9 passing tests).
- **Master Benchmark Runner**: `scripts/benchmark/run_all_benchmarks.py` (committed and reproducible).
- **Phase 2 & CV Audit Runner**: `scripts/benchmark/run_model_tuning_comparison.py` (audit documented in [`PHASE2_SEVENN_VALIDATION_AUDIT.md`](./PHASE2_SEVENN_VALIDATION_AUDIT.md)).
- **Remote Git Repository**: `https://github.com/DonghLee57/study_synthetic_accessibility.git`
- **Audit Verification State**: Validated against Python 3.11.9, RDKit 2023.9.6, and all 12 cloned sub-repositories in `external/`.
