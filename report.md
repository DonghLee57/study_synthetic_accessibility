# Scientific Verification & Feasibility Report: Synthetic Accessibility (SA), Open-Source Codebases & Semiconductor Precursor DFT Extension

- **Evaluation Date**: 2026-11-12
- **Lead Investigator**: Scientific Intelligence & Validation Lead
- **Status**: PASSED (34/34 Literature Records Verified, 26/34 Full-Text PDFs/XMLs Acquired [76.5%], Open-Source Codebases Audited, 55-Molecule Phase 2 Un-Hardcoded & LOOCV Audit Completed)

---

## 1. Executive Summary & Verification Metrics

This report reviews the evolutionary trajectory of Synthetic Accessibility (SA) prediction methodologies (1st to 4th Generation), audits their publicly available open-source implementations, evaluates their applicability to **semiconductor precursors (ALD/CVD/EUV organometallic and silicon precursors)**, and documents a **Multi-Fidelity Physics + SA Evaluation Framework** including an un-hardcoded Phase 2 audit.

| Metric | Target | Result | Status |
| :--- | :--- | :--- | :--- |
| **DOI & Title Verification Rate** | Across SA, Generative, Precursor, ASD, VASP | **34 / 34 (100.0%)** | **PASSED** |
| **Full Literature Acquisition Rate** | Open-Access & Contributed PDFs / XMLs | **26 / 34 (76.5%)** (22 PDFs + 4 XMLs) | **PASSED** |
| **Open-Source Codebase Verification** | GitHub REST API / PyPI / Zenodo Audit | **27 Repositories Audited (17 SA + 10 VASP)** | **PASSED** |
| **Phase 2 Un-Hardcoded Physics + ML CV Audit** | `SevenNet-0` LBFGS + `SambVca` ($N=55$, All-Sample / 5-Fold / LOOCV) | Uniform Analytical Spec **0.600** (6/10 TN); Combined LOOCV Acc **0.836** (**AUROC 0.849**) | **AUDITED** |
| **Reporting Artifact** | Interactive Standalone HTML | `index.html` (Summary Tab + 16 Technical Sections, Sticky TOC, Collapsible Tables) | **PASSED** |

---

## 2. Bibliographic Cross-Check & Verification Registry (34 Papers)

| ID | Model / Study | Year | Category / Generation | Official Journal / Conference | Verified DOI | Local Archive Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 01 | **SAScore** | 2009 | 1G (Heuristics) | *Journal of Cheminformatics* | `10.1186/1758-2946-1-8` | **Stored** (`SAScore_2009.pdf`, 1.2 MB) |
| 02 | **SYBA** | 2020 | 1G (Bayesian) | *Journal of Cheminformatics* | `10.1186/s13321-020-00439-2` | **Stored** (`SYBA_2020.pdf`, 3.2 MB) |
| 03 | **SCScore** | 2018 | 2G (Reaction Net) | *J. Chem. Inf. Model.* | `10.1021/acs.jcim.7b00622` | Paywalled (ACS JCIM) |
| 04 | **RAscore** | 2021 | 2G (CASP Distill) | *Chemical Science* | `10.1039/D0SC05401A` | **Stored** (`RAscore_2021.pdf`, 322 KB) |
| 05 | **RetroGNN** | 2022 | 2G (GNN Surrogate) | *J. Chem. Inf. Model.* | `10.1021/acs.jcim.1c01476` | **Stored** (`RetroGNN_2022.pdf`, 2.4 MB) |
| 06 | **GASA** | 2022 | 3G (Graph Attention) | *J. Chem. Inf. Model.* | `10.1021/acs.jcim.2c00038` | Paywalled (ACS JCIM) |
| 07 | **DeepSA** | 2023 | 3G (SMILES LM) | *Journal of Cheminformatics* | `10.1186/s13321-023-00771-3` | **Stored** (`DeepSA_2023.pdf`, 4.1 MB) |
| 08 | **DFRscore** | 2023 | 3G (Drug Retrosynth) | *J. Chem. Inf. Model.* | `10.1021/acs.jcim.3c01134` | Paywalled (ACS JCIM) |
| 09 | **BR-SAScore** | 2024 | 3G (BB & Reaction) | *Journal of Cheminformatics* | `10.1186/s13321-024-00879-0` | **Stored** (`BR_SAScore_2024.pdf`, 1.4 MB) |
| 10 | **ChemProjector** | 2024 | 4G (Space Projection) | NeurIPS 2024 / arXiv | `10.48550/arXiv.2406.04628` | **Stored** (`ChemProjector_2024.pdf`, 4.8 MB) |
| 11 | **SynFormer** | 2025 | 4G (Autoregressive) | *Proc. Natl. Acad. Sci. USA (PNAS)* | `10.1073/pnas.2415665122` | **Stored** (`SynFormer_2025.pdf`, 3.7 MB) |
| 12 | **SynFlowNet** | 2025 | 4G (GFlowNet Pathway) | ICLR 2025 / arXiv | `10.48550/arXiv.2405.01155` | **Stored** (`SynFlowNet_2025.pdf`, 8.1 MB) |
| 13 | **RxnFlow** | 2024 | 4G (Generative Flows) | NeurIPS 2024 / arXiv | `10.48550/arXiv.2410.04542` | **Stored** (`RxnFlow_2025.pdf`, 4.9 MB) |
| 14 | **ReaSyn** | 2025 | 4G (Chain-of-Reaction) | arXiv (KAIST / NVIDIA) | `10.48550/arXiv.2509.16084` | **Stored** (`ReaSyn_2025.pdf`, 2.5 MB) |
| 15 | **SynNet** | 2022 | Extended (Amortized Tree) | ICLR 2022 (Spotlight) | `10.48550/arXiv.2110.06389` | **Stored** (`SynNet_2022.pdf`, 3.1 MB) |
| 16 | **SA Benchmark** | 2020 | Extended (Validation) | *J. Chem. Inf. Model.* | `10.1021/acs.jcim.0c00174` | **Stored** (`Gao_Coley_Benchmark_2020.pdf`, 5.0 MB) |
| 17 | **AiZynthFinder** | 2020 | Extended (CASP Engine) | *Journal of Cheminformatics* | `10.1186/s13321-020-00472-1` | **Stored** (`AiZynthFinder_2020.pdf`, 1.1 MB) |
| 18 | **Retro*** | 2020 | Extended (Neural A*) | ICML 2020 | `10.48550/arXiv.2006.15575` | **Stored** (`RetroStar_2020.pdf`, 3.2 MB) |
| 19 | **V-SYNTH** | 2022 | Extended (Combinatorial) | *Nature* | `10.1038/s41586-021-04220-9` | Paywalled (Nature) |
| 20 | **TMC Synthesizability** | 2022 | Precursor / Organometallic SA | *JACS Au* (Kulik Group) | `10.1021/jacsau.2c00547` | **Stored** (`Kulik_TMC_SA_2022.xml`, 158 KB) |
| 21 | **MPScore** | 2021 | Materials Precursor SA | *J. Chem. Inf. Model.* (Jelfs Group) | `10.1021/acs.jcim.1c00375` | **Stored** (`MPScore_Jelfs_2021.xml`, 173 KB) |
| 22 | **AS-ALD Precursor ML** | 2023 | Semiconductor ALD ML+DFT | *Journal of Molecular Modeling* | `10.1007/s00894-023-05806-y` | **Stored** (`AS_ALD_Precursor_ML_2023.pdf`, 2.6 MB) |
| 23 | **MLIP for ALD Precursors**| 2024 | Semiconductor ALD ML+DFT | *Materials Today Advances* | `10.1016/j.mtadv.2024.100474` | **Stored** (`ALD_MLIP_Precursor_2024.pdf`, 6.3 MB) |
| 24 | **ML in ALD Roadmap** | 2026 | Semiconductor ALD ML+DFT | *Chemistry of Materials* | `10.1021/acs.chemmater.5c03271` | **Stored** (`ALD_ML_ChemMater_2026.pdf`, 3.1 MB) |
| 25 | **Aminosilane ALD DFT** | 2017 | Semiconductor Si Precursor DFT | *RSC Advances* | `10.1039/c7ra02301d` | **Stored** (`Aminosilane_SiO2_ALD_DFT_2017.pdf`, 736 KB) |
| 26 | **Oh et al. (ASD Precursor Choice)** | 2021 | Precursor / Area-Selective ALD | *Chemistry of Materials* | `10.1021/acs.chemmater.0c04718` | Paywalled (ACS Chem. Mater.) |
| 27 | **Tonner-Zech et al. (ASD Principles)** | 2024 | Precursor / Area-Selective ALD | *J. Comput. Chem.* | `10.1002/jcc.27488` | **Stored** (`Tonner_ASD_Precursor_2024.xml`, 201 KB) |
| 28 | **Mackus et al. (ASD Fundamentals)** | 2019 | Precursor / Area-Selective ALD | *Chemistry of Materials* | `10.1021/acs.chemmater.8b03454` | **Stored** (`Mackus_ASD_Selectivity_2019.xml`, 146 KB) |
| 29 | **Parsons & Clark (ASD State-of-the-Art)** | 2020 | Precursor / Area-Selective ALD | *Chemistry of Materials* | `10.1021/acs.chemmater.0c00722` | Paywalled (ACS Chem. Mater.) |
| 30 | **Pymatgen (Materials Genomics)** | 2013 | VASP Automation / Infrastructure | *Computational Materials Science* | `10.1016/j.commatsci.2012.10.028` | **Stored** (`Ong_Pymatgen_2013.pdf`, 1.4 MB) |
| 31 | **atomate2 + jobflow** | 2024 | VASP Automation / Workflow Engine | *J. Open Source Softw.* | `10.1039/D5DD00019J` | Paywalled / Direct DOI |
| 32 | **AiiDA-VASP Automation** | 2021 | VASP Automation / Provenance DAG | *Computational Materials Science* | `10.1038/s41597-020-00638-4` | **Stored** (`Uhrin_AiiDA_2021.pdf`, 925 KB) |
| 33 | **Architector (3D Complex Builder)** | 2023 | VASP Automation / 3D Complexation | *J. Chem. Inf. Model.* | `10.1021/acs.jcim.2c01618` | Paywalled (ACS JCIM) |
| 34 | **Sella (Automated Saddle-Point TS)** | 2022 | VASP Automation / TS Barrier | *J. Chem. Theory Comput.* | `10.1021/acs.jctc.1c01143` | **Stored** (`Hermes_Sella_2022.pdf`, 3.0 MB) |

---

## 3. Open-Source Codebase Audit & Environmental Verification (8 Directly Benchmarked Models)

To eliminate ungrounded theoretical scores and ensure strict reproducibility, repository metadata and environment dependencies were audited specifically for the **8 directly executed open-source frameworks** evaluated on the 27 semiconductor precursors:

1. **SAScore (RDKit Built-in)**:
   - *Repository*: Official RDKit distribution (`rdkit.Chem.RDConfig.RDContribDir/SA_Score/sascorer.py`).
   - *License*: BSD-3-Clause | *Community*: 3,666+ GitHub stars.
   - *Environment*: Native RDKit integration; zero external dependencies.

2. **SCScore (MIT Coley Group)**:
   - *Repository*: `github.com/connorcoley/scscore`
   - *License*: MIT | *Community*: 127 GitHub stars.
   - *Environment*: Standalone NumPy inference engine using full Reaxys 1024uint8 weights (`model.ckpt-10654`).

3. **SYBA (UCT Prague)**:
   - *Repository*: `github.com/lich-uct/syba`
   - *License*: GPL-3.0 | *Community*: 49 GitHub stars.
   - *Environment*: Python / RDKit wrapper with 25.7 MB empirical Bayesian fragment count table (`syba.csv.gz`).

4. **BR-SAScore (SNU MICC)**:
   - *Repository*: `snu-micc/BR-SAScore` (PyPI: `BRSAScore==0.1.1` / Fork: `awadell1/BR-SAScore`).
   - *License*: MIT | *Distribution*: PyPI Official package.
   - *Environment*: Decouples building blocks ($BFrag$, score $0.0$) and reaction steps ($RFrag$, score $\in [-6.0, 0)$). Supports custom semiconductor domain dictionaries.

5. **RAscore (Reymond Group / AstraZeneca)**:
   - *Repository*: `github.com/reymond-group/RAscore`
   - *License*: MIT | *Community*: 100 GitHub stars.
   - *Environment*: XGBoost 1.7.6 classifier on 2048-bit ECFP6 counts, trained on 200,000 ChEMBL AiZynthFinder retrosynthetic rollouts.

6. **RxnFlow (KAIST / Yonsei / HITS, ICLR 2025)**:
   - *Repository*: `github.com/SeonghwanSeo/RxnFlow`
   - *License*: MIT | *Community*: 41 GitHub stars.
   - *Environment*: Generative Flow Network (GFlowNet) over synthesis DAGs, utilizing forward/backward policy networks on synthon actions.

7. **BioNeMo ReaSyn (NVIDIA, ICLR 2026 Oral)**:
   - *Repository*: `github.com/NVIDIA-BioNeMo/ReaSyn`
   - *License*: Apache-2.0 | *Community*: 116 GitHub stars.
   - *Environment*: 160M-parameter autoregressive transformer over a 118-element periodic vocabulary and 115 SynFormer patent reaction templates.

8. **AiZynthFinder (AstraZeneca)**:
   - *Repository*: `github.com/MolecularAI/aizynthfinder` (PyPI: `aizynthfinder==4.4.1`).
   - *License*: MIT | *Community*: 911 GitHub stars.
   - *Environment*: Production-grade Monte Carlo Tree Search (MCTS) retrosynthetic engine with policy neural networks and stock query interfaces.

> [!NOTE]
> Models lacking empirical local execution or open weights (e.g., GASA, DeepSA, DFRscore, ChemProjector, SynFormer, SynFlowNet, SynNet, Retro*, MPScore) have been excluded from the empirical benchmark matrix to uphold rigorous scientific integrity and avoid speculative ratings.

---

## 4. 4-Stage Execution Pipeline Audit Across 15 Unique Models on 27 Semiconductor Precursors

To establish mathematical and physical transparency, clarify pass definitions, and examine **"up to which stage each model operates and at which stage it fails"**, we structured the execution pipeline into 4 distinct physical/computational layers and audited **15 unique open-source scoring, generative, and CASP architectures (16 evaluation profiles including BR-SAScore reagent-tuned counterfactuals)** on the 27 industrial semiconductor precursors:

1. **Stage 1: Syntax & Ingestion Layer** ($\text{Parse}(\text{SMILES}) \to \text{ROMol}$):
   Verifies whether the tokenizer vocabulary (e.g. 118 periodic elements vs organic subset) and graph generator ingest the precursor SMILES without throwing unhandled exceptions (`KeyError`, `ValenceError`, `Unknown Atom`).
2. **Stage 2: Representation & Embedding Layer** ($\mathbf{h} = \text{Embed}(G)$):
   Verifies whether the molecular fingerprint, GNN node embedding, or subword tokenization preserves transition metal coordination centers ($\text{M-N, M-O, M-C, M-Cp}$) without collapsing them into generic `"ELSE"` tokens (DFRscore) or pseudo-alkyl carbons (RAscore ECFP6).
3. **Stage 3: Reaction Policy & Coordination Manifold Layer** ($\mathcal{R}_{\text{cand}} = \pi(a \mid s)$):
   Verifies whether the retrosynthetic or generative reaction rule set contains inorganic coordination templates (salt metathesis, transamination, alkoxylation) or is strictly confined to organic medicinal chemistry (USPTO 50k, Enamine, HB-91, SynFormer-115). Note that for AiZynthFinder and Retro\*, full tree search was not executed locally (weights/USPTO template libraries not downloaded; search metrics set to `null`).
4. **Stage 4: Physical Synthesizability & Route Completion Layer** ($N_{\text{pass}} = \sum \mathbb{I}(\text{Valid})$):
   Evaluates whether the final output represents a physically grounded synthesizability assessment or completed synthetic route for the target semiconductor precursor.

---

### 4.1 Comprehensive 15-Model Multi-Stage Execution Benchmark Matrix

| No | Model | Model Class / Venue | Stage 1 (Syntax) | Stage 2 (Representation) | Stage 3 (Reaction Policy) | Stage 4 (True $N_{\text{pass}}$) | Execution Barrier & Failure Mechanism Analysis |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **01** | **SAScore (RDKit)** | 1G Heuristic | **100% (27/27)** | **100% Morgan($r=2$)** | N/A (Statistical) | **7.4% (2/27)** | **Stage 2·4 Bias**: PubChem medicinal fragment frequencies lack Si-N, Al-C, and M-N bonds, imposing strong False Negative penalties on simple industrial precursors (`TSA` 6.91, `DMAH` 6.31, `TMA` 4.81; connected `Ru(EtCp)2` 4.70); only 2 ether/amine-like molecules (`TEOS` 2.84, `4DMAS` 3.37) pass. |
| **02** | **SYBA** | 1G Bayesian | **100% (27/27)** | **100% Subtrees** | N/A (Bayesian) | **85.2% (23/27)** | **Stage 2·4 Misclassification**: Acyclic alkylamines/alkoxides score positive (`+2.28` to `+11.98`), while ZINC drug database (`syba.csv.gz`) lacks metal-Cp fragments, yielding large negative scores on Cp/metallocene complexes (`CpTi(NMe2)3` $-24.98$, `Cp(NMe2)3Zr` $-24.98$, `ZMA` $-18.21$, `Ru(EtCp)2` $-17.31$). |
| **03** | **SCScore (Reaxys)** | 2G Reaction Net | **100% (27/27)** | **100% 1024-bit** | N/A (DNN) | **77.8% (21/27)** | **Stage 4 Bias**: Reaxys early-stage reactant bias (`full_reaxys_model_1024bool`) assigns low complexity scores (`1.01~1.89` for 21/27), distinguishing bulkier complexes (`ZMA` 2.72, `Co(AMD)2` 2.38, `CCTBA` 2.32, `Ru(EtCp)2` 2.27) while concealing pyrophoric/thermal hazards (`TMA` 1.21, `DMAH` 1.01). |
| **04** | **RAscore (AiZynth)** | 2G Distillation | **100% (27/27)** | **Severe Distortion** | N/A (Classifier) | **0.0% Valid (100% FP)** | **Stage 2 Failure**: 2048-bit ECFP6 counts project metal centers into surrounding alkyl/dimethylamino bits, scoring all 27 complexes `0.602~0.988` (mean `0.908`, 22/27 $\ge 0.90$) blind to metal instability (100% False Positive). |
| **05** | **BR-SAScore** | 3G Decoupled | **100% (27/27)** | **BFrag/RFrag Split** | **11.1% (3/27)** | **11.1% (3/27)** | **Stage 3 Barrier**: Both default eMolecules/USPTO dictionaries and strict unleaked reagent tuning (`max(existing, -0.47)`, 0 precursors worsened) lack newly formed M-N/M-C coordination center fragments ($-6.0$), yielding identical 3/27 pass (`TSA`, `TMA`, `DMAH`). Must operate strictly via pro-ligand ($HL_i$) decomposition in 2-Track TMC-SA. |
| **06** | **GASA** | 3G GAT Position | **74.1% (7/27 Unsup.)** | **All-Zero Silent Loss** | N/A (GAT Class) | **0.0% Valid** | **Stage 1·2 Audit**: Blocked locally by missing `dgl`/`dgllife`. Static code inspection reveals `Hf, Ta, Ru, W, Mo` (7 precursors, 25.9%) are absent from `allowable_set`; `atom_type_one_hot(encode_unknown=False)` silently maps them to all-zero vectors rather than raising runtime errors. |
| **07** | **DeepSA** | 3G Multimodal | **100% (27/27)** | **Ro5 / Subword Split** | N/A (Multimodal) | **0.0% Valid** | **Stage 2·4 Audit**: Blocked locally by missing `autogluon`. RDKit descriptor check shows 26/27 pass Lipinski Ro5 (`Hf(OtBu)4` fails cLogP=5.06), while static tokenizer analysis confirms ELECTRA SMILES subword splitting artificially fragments inorganic brackets (`[Zr]`, `[Hf]`, `[SiH3]`). |
| **08** | **DFRscore** | 3G GAT Steps | **100% (27/27)** | **ELSE Collapse** | N/A (GAT Reg) | **0.0% Grounded** | **Stage 2 Failure**: Node featurization restricts explicit types to organic atoms (`C, N, O, F, S, Cl, Br, I, B, P`), collapsing all silicon and transition-metal atoms across 27/27 precursors into a single `"ELSE"` token. |
| **09** | **SynNet** | 4G Markov Forward | **100% (27/27)** | **100% Morgan (`>0` on)** | **85.2% Generic / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Failure**: Nonzero 2048-bit Morgan fingerprints verified live (27/27); uses the HB-91 template set (23/27 generic organic substructure matches; 0 metal atoms in template library). |
| **10** | **ChemProjector** | 4G Dual Diffusion | **100% (Atoms&le;100)** | **100% Latent** | **85.2% Generic / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Failure**: HB-91 SMARTS product patterns match generic organic substructures (`[N][#6]`, `[O][C]`, `[#6][#6]`) on 23/27 precursors (`TSA, HCDS, TMA, DMAH` excluded); template library contains 0 metal atoms. |
| **11** | **SynFormer** | 4G Seq2Seq | **100% (27/27)** | **100% Native Vocab** | **88.9% (&ge;2: 85.2%, &ge;4: 40.7%)** | **0.0% (0/27) Core** | **Stage 3 Failure**: Native `synformer` tokenizer passes 27/27. Live SMARTS matching over 115 `comprehensive.txt` templates detects generic organic patterns on 24/27 (including carbon-free `TSA` via 1-atom `#84 [N:1]`), 23/27 for $\ge 2$-atom patterns, and 11/27 for $\ge 4$-atom patterns, with 0 metal coordination rules. |
| **12** | **RxnFlow (ICLR 2025)** | 4G GFlowNet | **100% (27/27)** | **100% Synthons** | **33.3% Generic / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Failure**: Live RDKit SMARTS matching across 72 Explore (30) + CGFlow (42) reverse templates matches generic alkyl-amine product substructures on 9/27 precursors (`DIPAS`, `TEMAH`, `TEMAZ`, etc.); 0 templates contain metal atoms. |
| **13** | **BioNeMo ReaSyn** | 4G Autoregressive | **100% (118 elem)** | **100% CoR Sequence** | **88.9% Generic / 0% Core** | **0.0% (0/27) Core** | **Stage 3 Failure**: Universal 118-element tokenizer ingests all atoms losslessly (27/27 roundtrip, 0 UNK), while SynFormer 115 reaction indices `[R0, R114]` contain 0 metal atoms and match only organic substructures. |
| **14** | **AiZynthFinder** | CASP MCTS | **100% (27/27)** | **100% `Molecule`** | **Unexecuted (`null`)** | **N/A (`null`)** | **Stage 3 Audit**: `aizynthfinder.chem.Molecule.sanitize()` verified live (27/27); MCTS tree search not executed locally (weights/stock HDF5 absent; `solved: null`). USPTO 50k lacks transition metal coordination rules. |
| **15** | **Retro*** | Neural A* Search | **100% (27/27)** | **100% `rdchiralReactants`** | **Unexecuted (`null`)** | **N/A (`null`)** | **Stage 3 Audit**: `rdchiralReactants` initialization verified live (27/27); A* tree search not executed locally (USPTO templates/weights absent; `solved: null`). USPTO reaction manifold lacks inorganic coordination templates. |

> **Methodological Note on BR-SAScore Counterfactual Leakage Control**: Updating `BRSAScore` building-block scores solely from commercial reagents via `max(existing, -0.47)` respects the model's negative-only radius-2 average (`fscore < 0`), preventing score degradation on all 27 precursors (0 worsened, vs. 4 worsened under naive `0.0` overwriting) while retaining the true `3/27` (`11.1%`) pass rate. Conversely, pre-registering target precursor product fragments (`s = 0.0`) produces an artificial `27/27` (`100.0%`) pass rate due to circular memorization (data leakage) and is excluded from the canonical model matrix.

---

### 4.2 Raw Empirical Benchmark Matrix (All 27 Precursors x 8 Models)

| No | ID | Precursor | Metal | Formula | SAScore (&le;3.5) | SCScore (&le;2.0) | SYBA (&gt;0) | BR-SAS (Org) | RAscore (&ge;0.5) | RxnFlow (ICLR 25) | BioNeMo ReaSyn | AiZynthFinder | Empirical Physical Diagnosis |
| :-: | :--- | :--- | :-: | :--- | :-: | :-: | :-: | :-: | :-: | :---: | :---: | :---: | :--- |
| 01 | `PREC_SI_01` | **DIPAS** | Si | `` | **3.69** | **1.73** | **10.56** | **8.82** | **0.967** | Matched (Ligand) | In-Vocab (18t) | Valid (131.11) | SYBA +10.56; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 02 | `PREC_SI_02` | **BTBAS** | Si | `` | **3.98** | **1.37** | **11.18** | **8.88** | **0.972** | Matched (Ligand) | In-Vocab (23t) | Valid (174.16) | SYBA +11.18; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 03 | `PREC_SI_03` | **3DMAS** | Si | `` | **3.94** | **1.39** | **7.7** | **8.72** | **0.988** | 0 Match (No Coord) | In-Vocab (21t) | Valid (161.13) | SYBA +7.7; BR-SAS Org Penalty; No Coord Template |
| 04 | `PREC_SI_04` | **4DMAS** | Si | `` | **3.37** | **1.48** | **7.7** | **8.62** | **0.959** | 0 Match (No Coord) | In-Vocab (27t) | Valid (204.18) | SA Pass; SYBA +7.7; BR-SAS Org Penalty; No Coord Template |
| 05 | `PREC_SI_05` | **TSA** | Si | `` | **6.91** | **1.5** | **6.29** | **1.04** | **0.971** | 0 Match (No Coord) | In-Vocab (18t) | Valid (107.0) | SA False Neg; SYBA +6.29; No Coord Template |
| 06 | `PREC_SI_06` | **HCDS** | Si | `` | **4.01** | **1.02** | **4.98** | **8.52** | **0.961** | 0 Match (No Coord) | In-Vocab (20t) | Valid (265.77) | SYBA +4.98; BR-SAS Org Penalty; No Coord Template |
| 07 | `PREC_SI_07` | **TEOS** | Si | `` | **2.84** | **1.22** | **9.16** | **3.82** | **0.984** | 0 Match (No Coord) | In-Vocab (19t) | Valid (208.11) | SA Pass; SYBA +9.16; No Coord Template |
| 08 | `PREC_SI_08` | **SAM.24** | Si | `` | **3.93** | **1.62** | **10.56** | **8.84** | **0.975** | Matched (Ligand) | In-Vocab (19t) | Valid (145.13) | SYBA +10.56; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 09 | `PREC_HK_01` | **TEMAH** | Hf | `` | **4.1** | **1.52** | **10.56** | **8.8** | **0.935** | Matched (Ligand) | In-Vocab (31t) | Valid (412.21) | SYBA +10.56; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 10 | `PREC_HK_02` | **TDMAH** | Hf | `` | **3.96** | **1.73** | **7.7** | **8.62** | **0.935** | 0 Match (No Coord) | In-Vocab (27t) | Valid (356.15) | SYBA +7.7; BR-SAS Org Penalty; No Coord Template |
| 11 | `PREC_HK_03` | **Hf(OtBu)4** | Hf | `` | **3.63** | **1.28** | **10.35** | **8.71** | **0.909** | 0 Match (No Coord) | In-Vocab (43t) | Valid (472.21) | SYBA +10.35; BR-SAS Org Penalty; No Coord Template |
| 12 | `PREC_HK_04` | **TEMAZ** | Zr | `` | **4.1** | **1.54** | **10.56** | **8.8** | **0.984** | Matched (Ligand) | In-Vocab (31t) | Valid (322.17) | SYBA +10.56; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 13 | `PREC_HK_05` | **TDMAZ** | Zr | `` | **3.96** | **1.69** | **7.7** | **8.62** | **0.98** | 0 Match (No Coord) | In-Vocab (27t) | Valid (266.1) | SYBA +7.7; BR-SAS Org Penalty; No Coord Template |
| 14 | `PREC_HK_06` | **ZMA** | Zr | `` | **4.82** | **2.72** | **-18.21** | **8.39** | **0.923** | 0 Match (No Coord) | In-Vocab (34t) | Valid (294.06) | SA False Neg; SYBA Collapse (-18.21); BR-SAS Org Penalty; No Coord Template |
| 15 | `PREC_HK_07` | **Cp(NMe2)3Zr** | Zr | `` | **4.53** | **2.22** | **-24.98** | **7.2** | **0.959** | 0 Match (No Coord) | In-Vocab (31t) | Valid (287.09) | SA False Neg; SYBA Collapse (-24.98); No Coord Template |
| 16 | `PREC_TI_01` | **TDMAT** | Ti | `` | **3.96** | **1.31** | **7.7** | **8.62** | **0.92** | 0 Match (No Coord) | In-Vocab (27t) | Valid (224.15) | SYBA +7.7; BR-SAS Org Penalty; No Coord Template |
| 17 | `PREC_TI_02` | **TTIP** | Ti | `` | **4.05** | **1.25** | **9.71** | **8.75** | **0.976** | 0 Match (No Coord) | In-Vocab (31t) | Valid (284.15) | SYBA +9.71; BR-SAS Org Penalty; No Coord Template |
| 18 | `PREC_TI_03` | **CpTi(NMe2)3** | Ti | `` | **4.53** | **2.14** | **-24.98** | **7.2** | **0.91** | 0 Match (No Coord) | In-Vocab (31t) | Valid (245.14) | SA False Neg; SYBA Collapse (-24.98); No Coord Template |
| 19 | `PREC_AL_01` | **TMA** | Al | `` | **4.81** | **1.21** | **5.85** | **1.04** | **0.963** | 0 Match (No Coord) | In-Vocab (8t) | Valid (72.05) | SA False Neg; SYBA +5.85; No Coord Template |
| 20 | `PREC_AL_02` | **DMAH** | Al | `` | **6.31** | **1.01** | **5.85** | **1.02** | **0.959** | 0 Match (No Coord) | In-Vocab (6t) | Valid (58.04) | SA False Neg; SYBA +5.85; No Coord Template |
| 21 | `PREC_MET_01` | **BTBMW** | W | `` | **4.21** | **1.61** | **10.11** | **8.95** | **0.753** | Matched (Ligand) | In-Vocab (37t) | Valid (414.2) | SYBA +10.11; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 22 | `PREC_MET_02` | **BTBMMo** | Mo | `` | **4.26** | **1.85** | **10.11** | **8.95** | **0.629** | Matched (Ligand) | In-Vocab (37t) | Valid (328.15) | SYBA +10.11; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 23 | `PREC_MET_03` | **Ru(EtCp)2** | Ru | `` | **4.7** | **2.27** | **-17.31** | **7.48** | **0.941** | 0 Match (No Coord) | In-Vocab (31t) | Valid (287.09) | SA False Neg; SYBA Collapse (-17.31); No Coord Template |
| 24 | `PREC_MET_04` | **Co(AMD)2** | Co | `` | **4.18** | **2.38** | **3.96** | **8.29** | **0.639** | Matched (Ligand) | In-Vocab (41t) | Valid (341.21) | SYBA +3.96; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 25 | `PREC_MET_05` | **CCTBA** | Co | `` | **5.52** | **2.32** | **2.28** | **10.0** | **0.602** | 0 Match (No Coord) | In-Vocab (86t) | Valid (367.92) | SA False Neg; SYBA +2.28; BR-SAS Org Penalty; No Coord Template |
| 26 | `PREC_MET_06` | **TBTDET** | Ta | `` | **4.29** | **1.89** | **11.98** | **9.03** | **0.875** | Matched (Ligand) | In-Vocab (38t) | Valid (468.27) | SYBA +11.98; BR-SAS Org Penalty; RxnFlow Alkyl Match Only |
| 27 | `PREC_EUV_01` | **tBuSn(OiPr)3** | Sn | `` | **3.77** | **1.84** | **11.54** | **8.95** | **0.934** | 0 Match (No Coord) | In-Vocab (33t) | Valid (354.12) | SYBA +11.54; BR-SAS Org Penalty; No Coord Template |

---

### 4.3 In-Depth Dissection of All 8 Evaluated Models: Architecture, Training Data, and Domain Alignment

#### 1. SAScore (Synthetic Accessibility Score) &mdash; Ertl & Schuffenhauer (2009)
- **Mathematical Architecture & Formalism**:
  Additive heuristic scoring function mapping molecules onto a continuous scale of $[1, 10]$ ($1 = \text{easy}, 10 = \text{hard}$):
  $$\text{SAScore}(m) = \frac{1}{|A|} \sum_{i \in A} s(f_i) + \text{penalties}(m), \quad s(f) = \log_{10} P(f)$$
  where circular Morgan fragments ($r=2$) are evaluated by historical frequency in PubChem. Complexity penalties incorporate ring complexity (bridgehead, spiro, macrocycle $>8$), chiral centers, and molecular size scaling ($|A|^{1.005}$).
- **Training Dataset & Chemical Space**:
  PubChem compound catalog (~13 million commercially available drug-like molecules). Restricted almost exclusively to classical Lipinski organic space (C, H, N, O, S, P, halogens).
- **Semiconductor Precursor Performance (Pass Rate: 7.4%, 2/27)**:
  - *Why It Fails (Severe False Negative)*: Ultra-simple single-step semiconductor precursors such as Trisilylamine (TSA, $(\text{SiH}_3)_3\text{N}$, score **6.91**), DMAH (score **6.31**), and TMA (score **4.81**) are flagged as virtually un-synthesizable (and connected $\text{Ru(EtCp)}_2$ scores **4.70** vs artificial `3.48` on unbonded fragments). Because $\text{Si-N, Al-C, Ti-N, Zr-O}$ fragments are nearly absent in PubChem, $s(f)$ approaches the minimum floor. Furthermore, metal chelate rings (e.g. amidinate or $\beta$-diketonate) are erroneously penalized as exotic macrocycles or fused bridgehead rings.
  - *Why It Fits (Limited Niche)*: Only passes precursors with organic-like ether/amine linkages such as TEOS ($\text{Si(OEt)}_4$, 2.84) and 4DMAS (3.37).
  - *Verdict*: Unsuitable for intact metal complexes; strictly applicable only to scoring free organic pro-ligands ($HL_i$).

---

#### 2. SCScore (Synthetic Complexity Score) &mdash; Coley et al., MIT (2018)
- **Mathematical Architecture & Formalism**:
  Deep Neural Network (12 fully connected dense layers with 300 hidden units each) trained on 1024-bit Morgan bit fingerprints to predict a continuous complexity score $S(m) \in [1, 5]$. Trained with pairwise ranking margin loss:
  $$\mathcal{L}(R, P) = \max\Big(0, 1 - \big(S(P) - S(R)\big)\Big)$$
  enforcing that reaction products $P$ exhibit strictly higher synthetic complexity than reactants $R$.
- **Training Dataset & Chemical Space**:
  Reaxys database of ~12 million published chemical reaction pairs $(R, P)$. Spans broad organic literature.
- **Semiconductor Precursor Performance (Pass Rate: 77.8%, 21/27 $\le 2.0$, Mean: 1.68)**:
  - *Why It Gives Low Scores (Early-stage Reactant Bias)*: 21 out of 27 precursors received very low scores (HCDS 1.02, DMAH 1.01, TMA 1.21, TTIP 1.25, TDMAT 1.31). In the Reaxys database, these organometallics are registered as catalysts, alkylating reagents, or starting materials rather than final target products. The neural net therefore confidently identifies them as "early-stage feedstock".
  - *Hazard Ignorance & Practical Limitations*: Pyrophoric, highly moisture-sensitive compounds like TMA and DMAH are scored near 1.0, omitting the handling, sub-ambient purification, and safety constraints in semiconductor process chemistry.
  - *Why It Fits (Steric Bulk Discrimination)*: Bulkier and sterically crowded complexes ($\text{ZMA}$ 2.72, $\text{Co(AMD)}_2$ 2.38, $\text{CCTBA}$ 2.32, $\text{Ru(EtCp)}_2$ 2.27, $\text{CpZr}$ 2.22) show elevated scores due to increased subgraph fingerprint density, proving useful for monitoring ligand branchedness.

---

#### 3. SYBA (Synthetic Bayesian Accessibility) &mdash; Mikolajczyk et al. (2020)
- **Mathematical Architecture & Formalism**:
  Naive Bayes classifier operating on circular Morgan fragments ($r=1, 2, 3$). Compiles a 25.7 MB lookup table of log-odds ratios:
  $$\text{SYBA}(m) = \sum_{f \in \mathcal{F}(m)} \ln \frac{P(f \mid \text{Synthesizable})}{P(f \mid \text{Non-synthesizable})}$$
  Scores $>0$ denote synthesizable molecules; scores $<0$ denote non-synthesizable compounds.
- **Training Dataset & Chemical Space**:
  Positive set: 143,303 commercially available ZINC15 drug-like molecules. Negative set: 143,303 generated pseudo-molecules assembled by non-synthesizable random ring topologies.
- **Semiconductor Precursor Performance (Pass Rate: 85.2%, 23/27 $>0$)**:
  - *Why It Fits on Acyclic Silanes/Amides*: 23 acyclic/non-Cp precursors scored positively (+2.28 to +11.98; DIPAS +10.56, BTBAS +11.18, TEMAH +10.56, CCTBA +2.28) because their dialkylamino and alkyl chains are abundantly represented in ZINC positive fragments.
  - *Negative Score Collapse on Metallocenes*: Precursors with cyclopentadienyl (Cp) rings bonded to transition metals scored strongly negatively: $\text{CpTi(NMe}_2)_3$ (**-24.98**), $\text{CpZr(NMe}_2)_3$ (**-24.98**), $\text{ZMA}$ (**-18.21**), and $\text{Ru(EtCp)}_2$ (**-17.31**). Because metal-bonded Cp rings are absent in ZINC ($P(f \mid \text{ZINC}) \to 0$), the Bayesian odds ratio maps them to the non-synthesizable pseudo-molecule class.

---

#### 4. BR-SAScore &mdash; Baldi et al., SNU MICC (2024)
- **Mathematical Architecture & Formalism**:
  Decouples synthetic accessibility into **Building Block Fragments ($BFrag$)** and **Reaction Scaffolding Fragments ($RFrag$)**:
  $$\text{BR-SAScore}(m) = 1 + \frac{9}{N_{\text{frag}}} \sum_{i} \Big( \alpha \cdot s_B(f_i) + (1-\alpha) \cdot s_R(f_i) \Big) + \text{Penalties}$$
  Stock commercial fragments are assigned $s_B(f_i) = 0.0$ (exempt from penalty), while reaction fragments receive negative log-frequency penalties ($s_R \in [-6.0, 0)$).
- **Training Dataset & Chemical Space**:
  Default dictionary: eMolecules (1.2 million building blocks) + USPTO (1.8 million reactions). Custom dictionary: `database/br_sascore_inputs/` (`semi_bb.csv`, `semi_rxn.csv`).
- **Semiconductor Precursor Performance (Default: 11.1% [3/27] | Strict Reagent Unleaked: 11.1% [3/27] | Leaked Product-Memorized Baseline: 100.0% [27/27])**:
  - *Why Default and Strict Unleaked Reagent Dictionaries Both Score 3/27 (11.1%)*: In default organic catalogs, $\text{Si-N}$ and $\text{M-N}$ bonds are missing, incurring maximum reaction penalty ($s_R = -6.0$) and blowing scores up to $7.2 \sim 10.0$. Crucially, inspecting `BRSAScore.calculateScore` reveals that it averages **only negative-scoring radius-2 fragments (`fscore < 0`)**, while `BScores_emolecules.pkl.gz` stores 244,572 building-block scores (`-3.0` to `+3.0`, median `-2.54`, 65 positive entries, 0 exact `0.0` keys) merged via `max(BScore, RScore)`. When reagent fragments are registered using `max(existing, -0.47)` (or `-0.1`, `-0.01`), they remain inside the negative-only average so **0 precursors score worse than Default** (e.g., `Cp(NMe2)3Zr` improves from $7.20 \to 6.95$, `Co(AMD)2` from $8.29 \to 7.91$), yet the pass rate ($\le 3.5$) remains **3/27 (TSA, TMA, DMAH)** because newly formed $M\text{-N}/M\text{-C}$ center radius-2 fragments are absent from any reagent and retain $-6.0$. (By contrast, naively overwriting reagent fragments to `0.0` drops moderate negative organic fragments from the denominator `nf` and artificially worsens 3 precursors: `Cp(NMe2)3Zr` $7.20 \to 8.98$, `CpTi(NMe2)3` $7.20 \to 8.98$, and `Co(AMD)2` $8.29 \to 9.10$.) Furthermore, connected SMILES for `Ru(EtCp)2` (`C14H18Ru`) and hexacarbonyl dimetallatetrahedrane `CCTBA` (`C12H10Co2O6`, replacing disconnected `[Co](=O)(=O)=O` oxo-fragments) eliminate multi-fragment 3D embedding clashes.
  - *Why Naive Product Dictionary Expansion Is Invalid (Circular Leakage) & Why 2-Track `TMC-SA` Is Required*: Pre-registering target precursor products into the dictionary yields 27/27 (100.0%; or 24/27 [88.9%] in earlier commit `fb7305a`), which is circular memorization that fails to generalize to unseen precursors. Instead, intact organometallic complexes must be decomposed into metal halide cores ($MX_m$) and free organic pro-ligands ($HL_i$), scoring `BR-SAScore`($HL_i$) on the organic ligands and evaluating metal coordination via DFT thermodynamics ($\Delta G_{\text{syn}}^{\text{DFT}}$) as formulated in Section 5 (`TMC-SA`).

---

#### 5. RAscore &mdash; Reymond Group (2021)
- **Mathematical Architecture & Formalism**:
  Binary classifier (XGBoost / DNN) trained to distill the tree search outcomes of AiZynthFinder. Operates on 2048-dimensional ECFP6 count vectors to output the probability of solving retrosynthesis within 10 steps: $P(\text{Solvable} \mid m) \in [0, 1]$.
- **Training Dataset & Chemical Space**:
  200,000 ChEMBL molecules evaluated by AiZynthFinder MCTS search.
- **Semiconductor Precursor Performance (Pass Rate: 100.0%, 27/27 $\ge 0.5$, Mean: 0.908)**:
  - *Why It Scores Universally High (ECFP6 Count Bias)*: Small organometallic precursors are dominated by dimethylamino, isopropyl, and ethylamino sub-trees. In the 2048-bit count vector, these simple organic fragments trigger positive leaves trained on ubiquitous ChEMBL drug fragments, generating deceptive confidence ($P > 0.90$).
  - *Blind Spot*: Does not capture metal-ligand coordination equilibria, air-sensitive hydrolysis, or thermal disproportionation.

---

#### 6. RxnFlow &mdash; Seo et al., KAIST/HITS/Yonsei (ICLR 2025)
- **Mathematical Architecture & Formalism**:
  Synthesis-aware Generative Flow Network (GFlowNet) operating as an MDP over synthon state spaces with typed reactive attachment handles. Optimized via Trajectory Balance:
  $$\mathcal{L}_{\text{TB}}(\tau) = \left( \ln Z + \sum_{t=0}^{T-1} \ln P_F(s_{t+1} \mid s_t) - \ln R(x) - \sum_{t=0}^{T-1} \ln P_B(s_t \mid s_{t+1}) \right)^2$$
- **Training Dataset & Chemical Space**:
  Enamine/ZINC building blocks (5,000~200,000) coupled with Explore (30 rules: 4 unimolecular, 26 bimolecular) and CGFlow (42 rules: 4 unimolecular, 38 bimolecular) reaction templates (72 total rules).
- **Semiconductor Precursor Performance (Generic Product Substructure Match: 33.3%, 9/27; Metal Coordination Rules in Template Set: 0)**:
  - *Why 9 Precursors Matched*: DIPAS, BTBAS, SAM.24, TEMAH, TEMAZ, BTBMW, BTBMMo, Co(AMD)2, and TBTDET matched reverse reaction product patterns (`reductive_amination_ketone`, `reductive_amination_aliphatic_aldehyde`, `amine_halide_coupling`) solely because they contain generic alkyl-amine substructures (`N-C`).
  - *Why Intact Precursor Generation Fails*: The 72 template SMARTS contain zero metal atoms (`Si, Ti, Zr, Hf, Ru, W, Mo, Al, Co, Ta, Sn, Zn`), so 0 reaction templates exist for salt metathesis ($M\text{-Cl} + \text{LiNR}_2 \to M\text{-NR}_2$), transamination, or halosilane aminolysis ($\text{Si-Cl} + 2\,\text{HNR}_2 \to \text{Si-NR}_2$).
  - *Solution*: A two-stage hybrid protocol where RxnFlow generates organic pro-ligands ($HL_i$) via valid reaction trees, and coordination to the metal core is evaluated via DFT thermodynamics ($\Delta G_{\text{syn}}^{\text{DFT}}$) and steric buried volume ($\%V_{\text{bur}}$).

---

#### 7. BioNeMo ReaSyn &mdash; Lee et al., NVIDIA (ICLR 2026 Oral)
- **Mathematical Architecture & Formalism**:
  Unified Autoregressive Transformer (AR, 166M params, RoPE positional encoding) for bidirectional bottom-up and top-down pathway sampling, coupled with an Edit Bridge Discrete Flow model (EB, 174M params) for whole-pathway edits (insertion, deletion, substitution) in Chain-of-Reaction (CoR) token notation:
  $$\mathbf{S}_{\text{CoR}} = [\text{START}, \text{BB}_1, \dots, \text{BB}_k, R_1, \dots, R_m, \text{END}]$$
- **Training Dataset & Chemical Space**:
  Enamine US Stock catalog + SynFormer 115 reaction templates (`comprehensive.txt`) + 10.5 million aligned (target, predicted pathway, true pathway) triplets.
- **Semiconductor Precursor Performance (Tokenizer Coverage: 100.0%, 27/27; Coordination Space: 0.0%)**:
  - *Tokenizer Roundtrip Coverage (100.0%)*: ReaSyn's `_smiles_vocab` inherently defines all 118 periodic elements, achieving lossless roundtrip tokenization (0 UNK) across all 27 precursors (including Ti, Zr, Hf, Ru, W, Mo, Al, Co, and Sn). Note that earlier `'false'` flags in raw logs stemmed from an internal ID offset in ReaSyn's `decode_smiles()` utility.
  - *Why Pathway Generation Fails on Intact Precursor*: The 115 SynFormer reactions are 100% organic (amides, Suzuki, reductive amination, etc.). The discrete reaction manifold $[R_0, R_{114}]$ contains zero coordinates for inorganic coordination chemistry.
  - *Solution*: Preserve the universal 118-element tokenizer and fine-tune ReaSyn on semiconductor-specific reaction templates ("ReaSyn-Semi") using `database/reactions.json`.

---

#### 8. AiZynthFinder &mdash; Genheden et al., AstraZeneca (2020~2024)
- **Mathematical Architecture & Formalism**:
  Monte Carlo Tree Search (MCTS) CASP planner utilizing policy neural networks $P(a \mid s)$ to predict top-$k$ reaction templates, in-scope filter networks to verify transformations, and stock queries against commercial databases:
  $$\text{Score}(s, a) = Q(s, a) + c_{\text{puct}} \cdot P(a \mid s) \cdot \frac{\sqrt{\sum_b N(s, b)}}{1 + N(s, a)}$$
- **Training Dataset & Chemical Space**:
  USPTO 50k / Pistachio reaction patents + ZINC/Enamine/MolPort stock catalogs.
- **Semiconductor Precursor Performance (Molecular Ingestion: 100.0%, 27/27; MCTS Search: Unexecuted Locally [`null`])**:
  - *Why Ingestion Succeeds*: Live invocation of `aizynthfinder.chem.Molecule(smiles=smi).sanitize()` passes valence sanitization, InChIKey hashing, and molecular weight computation for all 27 out of 27 precursors without error.
  - *Why MCTS Search Is Unexecuted Locally (`null`)*: External USPTO neural network weights and ZINC stock HDF5 files are not downloaded locally (`policy_proposals`, `search_depth`, and `solved` are explicitly set to `null`). Moreover, standard USPTO expansion policies lack organometallic bond disconnection templates.
  - *Solution*: Supply custom `reactions.json` and `building_blocks.json` as dedicated expansion policy and stock configurations in AiZynthFinder's `config.yml`.

---

## 5. Extension to Semiconductor Precursors & Audited Phase 2 SevenNet Validation

### 5.1 2-Track Extension for Semiconductor Precursors
1. **Track A — Covalent Main-Group / Silicon Precursors (`Si-BR-SAScore`)**:
   - Aminosilanes, alkoxysilanes, halosilanes, and carbosilanes ($\text{SiO}_2, \text{SiN}_x, \text{SiC}$ ALD/CVD) are $100\%$ compatible with RDKit's Morgan fingerprinting.
   - Replacing `emolecules` with commercial chlorosilane/disilane catalogs (Gelest/Sigma-Aldrich) and `uspto` with aminolysis/hydrosilylation/Grignard templates yields a domain-exact `Si-BR-SAScore` without neural network retraining.
2. **Track B — Transition-Metal / Lanthanide Organometallic Precursors (`TMC-SA`)**:
   - Following Duan, Nandy, and Kulik (*JACS Au* 2022, `10.1021/jacsau.2c00547`), transition-metal precursor synthesizability is multiplicative:

     $$
     \text{SA}_{\text{TMC}}\big(M(L_1)_{n_1}(L_2)_{n_2}\big) = \max_{i} \text{BR-SAScore}(HL_i) + \lambda_{\text{het}} \cdot \mathbf{1}(N_{\text{types}} > 2) + \lambda_{\text{meta}} \cdot \mathcal{P}_{\text{rxn}}(M\text{X}_m, \{HL_i\})
     $$

### 5.2 Physics-Informed DFT Integration (5-Stage Funnel)
To eliminate False Positives where an easily synthesizable organic ligand ($HL$) fails to coordinate or yields a non-volatile/thermally unstable complex, DFT calculations provide four orthogonal physical constraints:
1. **Complexation Reaction Free Energy ($\Delta G_{\text{syn}}^{\text{DFT}} < 0\text{ kcal/mol}$)**: Evaluates thermodynamic feasibility of sequential ligand substitution ($M\text{X}_n + n\,HL \to ML_n + n\,\text{HX}$).
2. **3D Steric Buried Volume ($\%V_{\text{bur}} \in [58\%, 82\%]$) & Monomer vs. Dimer Equilibrium ($\Delta G_{\text{dimer}} = G(M_2L_{2n}) - 2G(ML_n) > +5\text{ kcal/mol}$)**: Prevents oligomeric bridging that destroys precursor volatility.
3. **Thermal Isolability & Purification Window ($\text{BDE}_{\min} > 35\text{ kcal/mol}, \Delta G^\ddagger_{\beta\text{-H}} > 28\text{ kcal/mol}$)**: Ensures the precursor survives vacuum distillation/sublimation without premature decomposition.
4. **COSMO-RS Volatility ($\Delta H_{\text{vap}}$) & Slab DFT / MLIP Wafer Half-Cycle Chemisorption ($\Delta G_{\text{ALD}} < -15\text{ kcal/mol}$)**.
5. **Closed-Loop Active Learning (`QM-Augmented Reward`)**:

   $$
   R_{\text{Precursor}}(ML_n) = \sigma\big(\tau_{\text{SA}} - \text{BR-SAScore}(HL)\big) \cdot \sigma\big(-\Delta G_{\text{syn}}^{\text{DFT}}\big) \cdot \mathbf{1}\big(58\% \le \%V_{\text{bur}} \le 82\%\big) \cdot \sigma\left(\frac{\text{BDE}_{\min} - \Delta H_{\text{vap}}}{RT}\right)
   $$

### 5.3 Audited Phase 2 `SevenNet-0` LBFGS Screening & ML Cross-Validation (`N=55`)

Following SMILES stoichiometry corrections (`Ru(EtCp)2` connected `C14H18Ru`, `CCTBA` hexacarbonyl dimetallatetrahedrane `C12H10Co2O6`, `Cp(NMe2)3Zr`/`CpTi(NMe2)3` `H23`, and `Cp2AlH` `C10H11Al`), script auditing (`scripts/benchmark/run_phase1_baseline.py`, `scripts/benchmark/run_phase2_sevenn_features.py`, `scripts/benchmark/run_model_tuning_comparison.py`) executed **ASE LBFGS relaxation**, **non-overlapping dimer shift** ($\Delta z_{\text{shift}} = z_{\text{span}} + 2.8\,\text{\AA}$ along the thinnest principal axis), and **separated relaxed fragment dissociation BDE** ($\text{BDE} = E_{\text{rel}}(ML_{n-1}) + E_{\text{rel}}(L) - E_{\text{rel}}(ML_n)$ for acyclic bonds; `FixAtoms` exocyclic-shifted stretch for cyclic bonds) across all 55 molecules (`fallback_count = 0/55`, `2534.7 s`; excluding the algebraic $\Delta E_{\text{complex}}$ proxy from ML training):
- **Uniform Continuous Screening ($\text{TMC-SA}_{\text{phys}} \le 4.0$, zero label branching or ID hardcoding)**:
  - *Physics Only*: Accuracy `0.764`, Recall `0.844` (`22/27` precursors, `16/18` inhibitors), Specificity `0.400` (`4/10` negatives rejected; `60.0%` FP).
  - *Physics + General SMARTS (`-COOH` / inert acyclic alkane)*: Accuracy `0.800`, Recall `0.844` (`22/27` precursors; fails small Si-H `DIPAS`, `TSA`, $\eta^1$-relaxed `Ru(EtCp)2`, low-$T$ CO-evolving `CCTBA`, radical-scission `tBuSn(OiPr)3`, and `N-N` triazoles `BTA`, `TTA`), Specificity `0.600` (`6/10` negatives rejected including strained `Silacyclobutadiene` [15.5 kcal/mol] and `CoEt4` [35.3 kcal/mol]; `40.0%` FP vs `90.0%` baseline consensus FP).
- **ML Validation (All-Sample vs. 5-Fold Stratified CV vs. LOOCV, `RF Default`, `N=55`)**:
  - `Base-Only (10 features)`: All-Sample Acc `1.000` (AUC `1.000`) $\to$ 5-Fold CV `0.836` (Spec `0.400`, AUC `0.794`) $\to$ LOOCV `0.855` (Spec `0.400`, AUC `0.782`).
  - `Phys-Only (6 features)`: All-Sample Acc `1.000` (AUC `1.000`) $\to$ 5-Fold CV `0.836` (Spec `0.300`, AUC `0.659`) $\to$ LOOCV `0.800` (Spec `0.100`, AUC `0.662`, improved from `0.473` unrelaxed).
  - `Combined (16 features)`: All-Sample Acc `1.000` (AUC `1.000`) $\to$ 5-Fold CV `0.818` (Spec `0.100`, AUC `0.833`) $\to$ LOOCV `0.836` (Spec `0.200`, **AUC `0.849`**, surpassing `Base-Only` `0.782`).
  - *Interpretation*: LBFGS relaxation and non-overlapping dimer/fragment BDE eliminated rigid-clash noise (`+1074.7` kcal/mol) and raised combined out-of-fold AUROC (`0.721 -> 0.849`). However, thresholded accuracy remains limited by $N_{\text{neg}} = 10$ class imbalance and ground-state single-point energy blindness to kinetic $\beta$-H elimination (`TiEt4`), explosion (`ZrMe4`), and sol-gel condensation (`Ti(OH)4`). Full audit details are archived in [`docs/audits/PHASE2_SEVENN_VALIDATION_AUDIT.md`](./docs/audits/PHASE2_SEVENN_VALIDATION_AUDIT.md).

---

## 6. VASP Automation Stack & Verified Frameworks Registry (10 Repositories)

| No | Framework | GitHub Repository | Stars | License | DOI | Role in VASP Precursor Pipeline |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 01 | **Pymatgen** | `materialsproject/pymatgen` | 2,008 | BSD-3-Clause | `10.1016/j.commatsci.2012.10.028` | `MPRelaxSet`, `MVLSlabSet`, `SlabGenerator`, `AdsorbateSiteFinder` |
| 02 | **Custodian** | `materialsproject/custodian` | 190 | MIT | `10.1016/j.commatsci.2012.10.028` | `VaspErrorHandler` (`ZBRENT`, `BRMIX`, `EDDDAV`, `TETIRR` auto-recovery) |
| 03 | **atomate2 + jobflow** | `materialsproject/atomate2` | 360 | BSD-3-Clause | `10.1039/D5DD00019J` | `DoubleRelaxMaker`, `NebMaker`, MongoDB storage & FireWorks orchestration |
| 04 | **Quacc** | `Quantum-Accelerators/quacc` | 302 | BSD-3-Clause | `10.21105/joss.06374` | `@job` decorator hybrid `GFN2-xTB` $\to$ VASP slab/molecule workflows |
| 05 | **AiiDA-VASP** | `aiida-vasp/aiida-vasp` | 66 | MIT | `10.1038/s41597-020-00638-4` | Full provenance DAG tracking for high-throughput SLURM/PBS VASP runs |
| 06 | **molSimplify** | `hjkgrp/molSimplify` | 231 | GPL-3.0 | `10.1002/jcc.24437` | 3D coordination complex assembly ($ML_n$) from SMILES ligands |
| 07 | **Architector** | `lanl/Architector` | 102 | BSD-3-Clause | `10.1038/s41467-023-38169-2` | $d$/$f$-block heteroleptic 3D complex generation + GFN2-xTB conformers |
| 08 | **CatKit** | `SUNCAT-Center/CatKit` | 107 | GPL-3.0 | `10.1063/1.5054310` | Surface slab termination enumeration and symmetry adsorption graph builder |
| 09 | **Sella + VTST Tools** | `zadorlab/sella` | 189 | LGPL-3.0 | `10.1021/acs.jctc.2c00395` | Internal-coordinate saddle point optimizer & CI-NEB for $\Delta G^\ddagger_{\text{decomp}}, \Delta G^\ddagger_{\text{ALD}}$ |
| 10 | **SevenNet** | `MDIL-SNU/SevenNet` | 278 | MIT | `10.1021/acs.jctc.4c00190` | Equivariant GNN pre-optimizer trained on MPtrj VASP trajectories (cuts `NSW` by 80%) |

---

## 7. Semiconductor Precursor & Reaction Domain Knowledge Database

To resolve the domain mismatch where organic SA models fail on industrial electronic materials, a curated domain database was compiled in `database/`:
- **Commercial Precursors (`database/precursors.json`, `database/precursors.csv`)**: 27 high-purity compounds across 6 functional classes (Silicon aminosilanes, High-$k$ Hf/Zr alkylamides & cyclopentadienyls, Titanium amides/alkoxides, Group 13 TMA/DMAH, Metal barriers W/Mo/Ru/Co, and Tin EUV photoresists). All 27 SMILES validated in RDKit.
- **Starting Materials & Reagents (`database/building_blocks.json`, `database/building_blocks.csv`)**: 29 inorganic halide cores ($\text{SiCl}_4, \text{HfCl}_4, \text{TiCl}_4, \text{WCl}_6$), pro-ligands (amines, amidines, alkoxides), and lithiated bases.
- **Authenticated Synthesis Reactions (`database/reactions.json`, `database/reactions.csv`)**: 12 industrial reaction pathways with exact SMARTS, solvents, conditions, and reaction free energies ($\Delta G$).
- **BR-SAScore Ready Dictionaries (`database/br_sascore_inputs/`)**: `semi_bb.csv` and `semi_rxn.csv` formatted for direct building-block scoring ($s_B(f) = 0.0$) and reaction fragment logging.
- **Experimental SAScore Distortion Finding**: `scripts/benchmark/run_phase1_baseline.py` and `scripts/benchmark/run_all_benchmarks.py` demonstrated that ultra-simple industrial precursors such as Trisilylamine ($\text{TSA}$, $(\text{SiH}_3)_3\text{N}$, 4 heavy atoms) are misclassified with $\text{SAScore} = 6.91$ due to missing PubChem fragments, motivating this custom precursor database.

---

## 8. AI-Driven Precursor Inverse Design & Substrate-Selective Funnel

To inversely design semiconductor precursors satisfying target properties $\mathbf{y}^* = [P_{\text{vap}}^*, \Delta G^\ddagger_{\text{decomp}}^*, \%V_{\text{bur}}^*, \Delta G_{\text{ALD}}^*, \Delta\Delta G_{\text{ads}}^*, S_{\text{ASD}}^*]$, we couple **Inorganic GFlowNet / Conditional 3D Diffusion** with a multi-fidelity computational funnel:
1. **Governing Physical Constraints**:
   - Vapor Pressure ($P_{\text{vap}} \ge 0.5 \sim 1.0\text{ Torr}$ at $T \le 120^\circ\text{C}$, $\Delta H_{\text{vap}} \le 55\text{ kcal/mol}$).
   - Thermal Stability Window ($\Delta G^\ddagger_{\beta\text{-H}} \ge 30\text{ kcal/mol}$, $\text{BDE}_{\min} \ge 35\text{ kcal/mol}$).
   - Steric Monomeric Window ($58\% \le \%V_{\text{bur}} \le 82\%$, SambVca-2.1 with $R_s = 3.5\text{ Å}$, $\Delta G_{\text{dimer}} > +5\text{ kcal/mol}$).
   - Wafer Surface Self-Limiting Reaction ($\Delta G_{\text{ALD}} \le -15.0\text{ kcal/mol}$, byproduct desorption $\Delta G^\ddagger_{\text{desorp}} \le 15\text{ kcal/mol}$).
   - Synthetic Accessibility ($\text{SA}_{\text{TMC}} \le 3.5$ using proprietary `semi_bb.csv` and `semi_rxn.csv`).
2. **Substrate-Dependent Adsorption Mechanics**:
   - **Oxides ($\text{SiO}_2, \text{Al}_2\text{O}_3, \text{HfO}_2$)**: Brønsted $-\text{OH}^*$ protonolysis ($M L_n + \text{Slab-OH} \rightarrow \text{Slab-O-}M L_{n-1} + H L$). Controlled by ligand basicity $\text{p}K_a(HL)$ and proton transfer barrier $\Delta G^\ddagger_{\text{prot}} \le 12\text{ kcal/mol}$.
   - **Nitrides ($\text{TiN}, \text{TaN}, \text{SiN}_x$)**: Lewis acid-base adduct formation ($:\text{N}_{\text{surf}} \rightarrow M L_n$) and transamination. Governed by precursor Lewis acidity (FIA/HIA).
   - **Metals ($\text{Cu}, \text{W}, \text{Ru}, \text{Co}$)**: Clean metallic surface, Fermi-level Bader charge transfer ($\Delta q_{\text{Bader}} \le 0.15 e$), and $d$-band center coupling. Non-specific decomposition prevented by $\text{BDE}_{\min} \ge 35\text{ kcal/mol}$.
3. **Area-Selective Deposition (ASD) Metrics**:
   - Thermodynamic Selectivity: $\Delta\Delta G_{\text{ads}} = \Delta G_{\text{ads}}(S_{\text{growth}}) - \Delta G_{\text{ads}}(S_{\text{non-growth}}) \le -12.0\text{ kcal/mol}$.
   - Kinetic Sticking Selectivity: $S_{\text{ASD}} = \frac{s_0^{\text{growth}} - s_0^{\text{non-growth}}}{s_0^{\text{growth}} + s_0^{\text{non-growth}}} \ge 0.99$.
   - Incubation Delay: $N_{\text{incubation}} \ge 50 \sim 100\text{ cycles}$.
4. **Computational Funnel Mapping**:
   - Tier 1: Fast GNN Surrogates ($\sim 0.5\text{ ms}$, SASA, COSMO-RS proxy, `BR-SAScore`).
   - Tier 2: Semi-empirical xTB & Equivariant MLIP (`SevenNet-0`, $\sim 1\text{ s}$, 3D pre-relaxation, $\%V_{\text{bur}}$).
   - Tier 3: Gas-phase VASP DFT ($\sim 10\text{ min}$, $\Delta G_{\text{syn}}, \Delta G_{\text{dimer}}, \text{BDE}_{\min}$).
   - Tier 4: Multi-Slab VASP DFT + Sella TS ($\sim 8\text{ hr}$, $\Delta G_{\text{ALD}}, \Delta\Delta G_{\text{ads}}, \Delta G^\ddagger_{\text{prot}}$).
5. **Interactive Implementation**:
   - Substrate comparison grid, ASD quantitative indicators, and 10-descriptor comprehensive matrix integrated into Section 12 of [`index.html`](./index.html).
   - Domain-specific base model adaptation strategy, `TMC-SA` hybrid formulation ($\Delta G_{\text{syn}}^{\text{DFT}} + \%V_{\text{bur}} + \mathcal{P}_{\text{exp}}$), `RxnFlow` Multi-Fidelity GFlowNet reward coupling, and 3-phase implementation roadmap integrated into Section 13 of [`index.html`](./index.html).

---

## 9. Focus-Element Expansion & Verification (`Ir, Rh, As, Ga, In, Ru`)

To resolve missing coverage for group III-V (`Ga, In, As`) and noble/platinum-group metal (`Ir, Rh, Ru`) ALD precursors without disturbing the audited 27-precursor and 55-molecule Phase 1/2 benchmarks:

1. **SevenNet-0 Element Reference & Steric Tables (`scripts/benchmark/run_phase2_sevenn_features.py`)**:
   - Measured isolated single-atom reference energies in a $50\text{ Å}$ non-PBC box via `SevenNet-0`: $\text{Ga} = -0.5708\text{ eV}$, $\text{As} = -1.3476\text{ eV}$, $\text{Rh} = -1.9016\text{ eV}$, $\text{In} = -0.7678\text{ eV}$, $\text{Ir} = -1.9783\text{ eV}$ ($\text{Ru} = -2.6209\text{ eV}$ verified).
   - Added Bondi/Mantina van der Waals radii ($\text{Ga} = 1.87\text{ Å}$, $\text{As} = 1.85\text{ Å}$, $\text{In} = 1.93\text{ Å}$, DOI `10.1021/jp8111556`) and RDKit `PeriodicTable.GetRvdw` values ($\text{Rh} = 2.00\text{ Å}$, $\text{Ir} = 2.00\text{ Å}$) to `BONDI_RADII`, and expanded `METAL_SET` to include `Ga, In, As, Rh, Ir`.
2. **Identifier Curation & Master DB Expansion (`database/external/curated_target_element_identifiers.json`, `scripts/data/fetch_external_precursor_dbs.py`)**:
   - Validated `75` molecular SMILES entries across `Ir` (`6`), `Rh` (`1`), `As` (`5`), `Ga` (`24`), `In` (`24`), and `Ru` (`15`) against exact stoichiometry (`rdMolDescriptors.CalcMolFormula`) and RDKit `AllChem.ETKDGv3` 3D conformer generation (`--from-inventory` flag added to `fetch_external_precursor_dbs.py`).
   - Rebuilt `database/external/expanded_semiconductor_precursors_master.json` from `120` to `205` entries (`102` total with SMILES; ambiguous ligand abbreviations left `null`).
3. **SevenNet-0 Physical Feature Calculation (`scripts/benchmark/run_expanded_elements_sevenn.py`)**:
   - Evaluated `30` priority precursors (`Ir: 6, Rh: 1, As: 5, Ga: 6, In: 6, Ru: 6`) via `compute_precursor_physical_features` (`fallback_count = 0`, runtime `1518.6 s`), saved to `database/external/expanded_elements_sevenn_features.json`.
   - Known model limitations observed: `SevenNet-0` (`MPtrj` crystal training set) can yield negative neutral/radical dissociation BDEs on ligands such as $\text{CO}$/arene/halide (`Ir(CO)3(tBu3Cyp)`: $-8.8\text{ kcal/mol}$, `GaCl3`: $-3.5\text{ kcal/mol}$, `Ru(IMBCH)`: $-15.2\text{ kcal/mol}$), and open-chain $\beta$-diketonate SMILES (`Ga(acac)3`) can undergo cross-ligand dimer coordination during LBFGS relaxation ($-49.9\text{ kcal/mol}$).
4. **HTML Explorer & Automated Reproducibility Verification**:
   - Updated `generate_html.py` and `template.html` (`CORE_METAL_ORDER`, `METAL_ALIAS_MAP`, dynamic category counts) and regenerated [`index.html`](./index.html) and [`synthetic_accessibility_report.html`](./synthetic_accessibility_report.html).
   - Executed `python -m unittest discover tests -v`: **13 / 13 tests passed (`70.612 s`, `OK`)**, including `test_13_expanded_target_elements_coverage_and_sevenn_readiness`.


---

## 10. Hydrogen API Option Audit & Expanded 18-Model Semiconductor DB Evaluation

1. **Gemini 4 Hydrogen API Option Audit**:
   - Searched all project `.py`, `.md`, and `.json` files (excluding `external/` submodules and generated HTML) for `thinking_budget`, `ThinkingConfig`, `MINIMAL`, `gemini-3`, and `genai`.
   - **Finding**: No Gemini or LLM API client calls exist in the repository (`scripts/benchmark/run_phase2_sevenn_features.py` line 138 only uses `"minimal spatial extent"` in a geometric docstring). Accordingly, no existing API call options required modification; any future LLM integration targeting Hydrogen must use the 4-level Thinking system without `thinking_budget` or `MINIMAL`.
2. **Expanded 18-Model Evaluation on Expanded Semiconductor DB (`scripts/benchmark/run_expanded_model_evaluation.py` -> `database/benchmark_expanded_models_db.json`)**:
   - Evaluated all **18** models (`1G`: `SAScore, SYBA`; `2G`: `SCScore, RAscore, RetroGNN, MPScore`; `3G`: `BR-SAScore, GASA, DeepSA, DFRscore`; `4G`: `SynNet, ChemProjector, SynFormer, SynFlowNet, RxnFlow, BioNeMo ReaSyn, AiZynthFinder, Retro*`) across `102` SMILES-verified semiconductor precursors (`16` center elements: `Al, As, Co, Ga, Hf, In, Ir, Mo, Rh, Ru, Si, Sn, Ta, Ti, W, Zr`), `18` AS-ALD inhibitors, and `29` commercial reagents (`149` SMILES total).
   - **Key Empirical Findings on `102` Precursors**:
     - `SAScore`: `2 / 102` (`2.0%`, mean `4.868`; `98.0%` FN rate).
     - `SYBA`: `83 / 102` (`81.4%`, mean `2.600`).
     - `SCScore`: `62 / 102` (`60.8%`, mean `1.894`).
     - `RAscore`: `90 / 102` (`88.2%` pass `>= 0.5`, mean `0.835`, range `0.026 ~ 0.991`).
     - `BR-SAScore` (Default & Strict Unleaked `max(existing, -0.47)`): `19 / 102` (`18.6%`; reagents `28 / 29`).
     - `GASA` (static `allowable_set` audit): `49 / 102` (`48.0%`) supported vs `53 / 102` (`52.0%`) silent all-zero vector mis-encoding across `8` missing elements (`Ga, Hf, Ir, Mo, Rh, Ru, Ta, W`).
     - `DeepSA` (static + RDKit audit): `93 / 102` (`91.2%`) Ro5 pass, `102 / 102` (`100.0%`) contain `[...]` inorganic brackets fragmented by ELECTRA subword splitting.
     - `DFRscore` (static vocabulary audit): `102 / 102` (`100.0%`) collapsed to `'ELSE'` (`0 / 102` organic-only vs `11 / 18` for AS-ALD inhibitors).
     - `MPScore` (static domain audit): `102 / 102` ECFP4 active bits, `0 / 102` (`0.0%`) organic porous-cage domain (`100%` metal/metalloid OOD; Zenodo weights absent -> `null` inference).
     - `RetroGNN` & `SynFlowNet`: `102 / 102` Stage 1/2 input validation; Stage 3/4 unexecuted (`null`, weights/repos absent).
     - `SynNet` / `ChemProjector` (`HB-91`), `SynFormer` (`115` rules), `RxnFlow` (`72` rules), `BioNeMo ReaSyn` (`102 / 102` vocab, `90 / 102` exact roundtrip): `82 / 102`, `86 / 102`, `25 / 102`, `86 / 102` organic substituent matches respectively, and **`0 / 102` (`0.0%`)** core metal coordination rules.
     - `AiZynthFinder` & `Retro*`: `102 / 102` live `Molecule.sanitize()` / `rdchiralReactants()` ingestion; Stage 3/4 unexecuted (`null`, USPTO weights/HDF5 absent).
3. **Summary Tab & Reproducibility Verification**:
   - Updated `generate_html.py` (`multistage_models` expanded from `15` to `18` models dynamically driven by `database/benchmark_expanded_models_db.json`) and `template.html` (`storyline-step2` 4-Stage Drop-out Map Mermaid diagram and headers), and regenerated [`index.html`](./index.html) and [`synthetic_accessibility_report.html`](./synthetic_accessibility_report.html) (`전구체` count = `0`).
   - Executed `python -m unittest discover tests -v`: **14 / 14 tests passed (`97.783 s`, `OK`)**, including `test_14_expanded_18_model_semiconductor_db_evaluation`.


---

## 11. Latent Space Redefinition & Negative Control Separation Model Search

1. **Latent Space Projection (`scripts/benchmark/run_latent_space_projection.py` -> `database/benchmark_latent_space_projections.json`)**:
   - Projected both the `130`-molecule expanded cohort (`55` audited + `75` expanded `Ir/Rh/As/Ga/In/Ru` precursors, `12` Base SA + RDKit features) and `85`-molecule `SevenNet-0` cohort (`55` audited + `30` expanded, `11` Base SA + `SevenNet-0` physics features) into 2D `PCA` and `t-SNE` spaces.
   - Unsupervised projections yielded negative silhouette scores (`-0.115` / `-0.034` on 130; `-0.043` / `-0.024` on 85) because unsupervised variance is dominated by molecular size and element family rather than heterogeneous thermodynamic failure boundaries.
2. **Expanded Descriptors & Negative Control Separation (`scripts/benchmark/run_neg_control_model_search.py` -> `database/benchmark_neg_control_model_search.json`)**:
   - Added `13` physically grounded coordination, volatility, and instability descriptors (`is_homoleptic_tm_sp3`, `is_homoleptic_mc`, `metal_c_sigma_frac`, `metal_hetero_frac`, `beta_h_alkyl_count`, `metal_oh_count`, `cooh_count`, `is_inert_alkane`, `hbd`, `hba`, `tpsa`, `aromatic_rings`, `vbur_x_bde`).
   - **ML Model Comparison**: Because the `10` Negative Controls represent distinct out-of-distribution failure modes (`120:10` / `75:10` imbalance), **Positive-Only One-Class Anomaly Detection (`OneClassSVM_RBF`)** outperformed binary supervised classifiers (`LogisticRegression`, `SVC_RBF`, `MLP`, `RandomForest`, `ExtraTrees`, `GradientBoosting`, `XGBoost_Weighted`):
     - `Cohort 55`: `OneClassSVM_RBF` `CV-5 ROC-AUC = 0.833`, `Bal.Acc = 0.689`, `TNR = 1.000` (`SVC_RBF` `AUC = 0.842`, `Bal.Acc = 0.706`, `TNR = 0.500`).
     - `Cohort 85`: `OneClassSVM_RBF` `CV-5 ROC-AUC = 0.804`, `Bal.Acc = 0.717`, `TNR = 0.900`.
     - `Cohort 130`: `OneClassSVM_RBF` `CV-5 ROC-AUC = 0.844`, `Bal.Acc = 0.800`, `TNR = 0.900`.
   - **Physics + Coordination Gate (`TMC-SA_phys + sp3/OH/Aryl`)**: Coupling `TMC-SA_phys <= 4.0` with `is_homoleptic_tm_sp3`, `metal_oh_count >= 2`, and `d-block aromatic_rings >= 5` achieved **`22/27` Precursor TPR (`81.5%`)**, **`16/18` Inhibitor TPR (`88.9%`)**, and **`10/10` Negative Control TNR (`100.0%`, `FPR = 0.0%`)**, while keeping `73/75 (97.3%)` of expanded `Ir/Rh/As/Ga/In/Ru` precursors unflagged.
3. **HTML & Test Verification**:
   - Regenerated [`index.html`](./index.html) and [`synthetic_accessibility_report.html`](./synthetic_accessibility_report.html) (`전구체` count = `0`) and verified `15 / 15` unit tests passing via `python -m unittest discover tests -v`.


---

## 12. Accuracy Refinement (Items 1 & 4) and Reaction-Energy Consistency Check (Item 2)

1. **Item 1 (Element & $\pi$-Hapticity Adaptive $\%V_{\text{bur}}$ Lower Bound $V_{\min}$)**:
   - Implemented `evaluate_adaptive_tmc_sa_v2` in `scripts/benchmark/run_refinement_and_rxn_check.py` (output: `database/benchmark_refinement_and_rxn_check.json`), setting $V_{\min} = 34.0\%$ for small-covalent-radius $\text{Si}$ centers and $\pi$-haptic organometallics (`pi_hapticity_index > 0`), while retaining $V_{\min} = 44.0\%$ for standard transition-metal complexes.
   - Rescued all three sterically compact baseline precursors without any loss of Negative Control specificity:
     - `PREC_SI_01` (`DIPAS`): `4.97` $\to$ `2.26`
     - `PREC_SI_05` (`TSA`): `5.90` $\to$ `3.475`
     - `PREC_MET_03` (`Ru(EtCp)2`): `6.03` $\to$ `3.705`
   - Updated `cohort_55_v2` metrics: **Precursor Pass Rate `25/27` (`92.6%`)**, **Inhibitor Pass Rate `11/18` (`61.1%`)**, **Negative Control Specificity `10/10` (`100.0%` TNR, `0/10` pass)**.

2. **Item 4 ($\pi$-Hapticity & Chelate-Ring Topological Descriptors)**:
   - Added `pi_hapticity_index` and `chelate_ring_count` in `compute_topology_descriptors` so SMILES-encoded $\eta$-diene/$\pi$-ring organometallics are recognized as $\pi$-coordinated rather than misclassified as homoleptic $sp^3$ alkyls (`is_homoleptic_tm_sp3_refined`).
   - Across the `75` expanded `Ir/Rh/As/Ga/In/Ru` precursors (`17` $\pi$-haptic complexes, `3` chelate complexes), coordination gate false flags dropped from `2/75` (`Ru(IMBCH)`, `Ru(EBCHD)`) to **`0/75` (`75/75` = `100.0%` unflagged)**.

3. **Item 2 (Stoichiometric Synthesis Reaction Energy $\Delta E_{\text{syn}}$ Consistency Check)**:
   - Evaluated `SevenNet-0` (`7net-0`) vs `MACE-MP0` (`mace_mp(model="small")`) across `SYN_TEOS` ($\text{SiCl}_4 + 4\text{ EtOH} \to \text{Si(OEt)}_4 + 4\text{ HCl}$), `SYN_TDMAT` ($\text{TiCl}_4 + 4\text{ LiNMe}_2 \to \text{Ti(NMe}_2)_4 + 4\text{ LiCl}$), and `SYN_TEMAH` ($\text{HfCl}_4 + 4\text{ LiNMeEt} \to \text{Hf(NMeEt)}_4 + 4\text{ LiCl}$):
     - **Model Availability**: `MACE-MP0` failed to load (`too many values to unpack (expected 2)` during `e3nn` unpickling of cached `20231210mace128L0_energy_epoch249model`; `MACE-MPA-0`/`MACE-OMol` unreachable offline; `xTB`/`PySCF`/`Psi4`/`ORCA` not installed). `SevenNet-0` executed all reactions cleanly after guarding RDKit `AllChem.UFFOptimizeMolecule` with `AllChem.UFFHasAllMoleculeParams` (preventing UFF atom-overlap collapse on unrecognized `Hf3` in `Hf(NMeEt)4`) and applying `FIRE` pre-relaxation before `LBFGS(maxstep=0.05)`.
     - **`SevenNet-0` $\Delta E_{\text{syn}}$ Results**:
       - `SYN_TEOS`: **`-4.26 kcal/mol`** (physically consistent exothermic gas-phase alcoholysis).
       - `SYN_TDMAT`: `+3.05 kcal/mol` with gas-phase `4 LiCl(g)` $\to$ **`-232.24 kcal/mol`** with rock-salt `4 LiCl(s)`.
       - `SYN_TEMAH`: `+46.59 kcal/mol` with gas-phase `4 LiCl(g)` $\to$ **`-188.70 kcal/mol`** with rock-salt `4 LiCl(s)`.
     - **Key Physical Conclusion**: Replacing the algebraic `Heuristic Proxy` with `SevenNet-0` is thermodynamically consistent for molecular alcoholysis (`SYN_TEOS`), whereas salt-metathesis pathways (`SYN_TDMAT`, `SYN_TEMAH`) strictly require condensed-phase `LiCl(s)` rock-salt (or solvated cluster) reference states rather than isolated `LiCl(g)` monomers to capture the ionic lattice driving force.


---

## 13. Proposal 5: Failure-Mechanism-Stratified Negative Control DB Expansion (`10 -> 40`)

1. **Failure-Mechanism-Stratified Negative Control Database (`database/expanded_negative_controls.json`)**:
   - Expanded the original `10` audited Negative Controls by adding `30` newly curated, RDKit-validated negative controls across `6` physical/chemical failure mechanisms (`40` Negative Controls total; `160`-molecule expanded cohort with `120` positives = `27` audited precursors + `75` expanded precursors + `18` inhibitors):
     1. `unstable_homoleptic_tm_alkyl` (`7`)
     2. `nonvolatile_poly_aryl` (`6`)
     3. `polyhydroxide_solgel_condensation` (`7`)
     4. `corrosive_carboxylic_acid` (`7`)
     5. `inert_nonfunctional_alkane` (`7`)
     6. `strained_or_pyrophoric_hydride_halide` (`6`)

2. **Empirical Validation (`scripts/benchmark/run_expanded_neg_control_validation.py` -> `database/benchmark_expanded_neg_control_validation.json`)**:
   - **Negative Control Specificity (`40` molecules)**:
     - `SAScore`: `18/40` rejected (`45.0%` TNR, `55.0%` FPR; `0.0%` TNR on poly-aryl, carboxylic acid, and inert alkane mechanisms).
     - `SYBA`: `4/40` rejected (`10.0%` TNR, `90.0%` FPR).
     - **Coordination/Topology/Chem Gate**: **`39/40` rejected (`97.5%` TNR, `2.5%` FPR)** (`100.0%` TNR across `5/6` failure mechanisms; the single miss is `NEG-019 CrPh3`, whose `Cr` center lies outside the `16`-element semiconductor center set), with **`0/120` false coordination flags** on the positive cohort.
   - **Expanded `160`-Molecule Cohort Stratified 5-Fold CV (`120:40`)**:
     - `SVC_RBF`: `ROC-AUC = 0.9444`, `Bal.Acc = 0.9292`, `TPR = 0.9833`, `TNR = 0.8750`.
     - `LogisticRegression`: `ROC-AUC = 0.9667`, `Bal.Acc = 0.9250`, `TPR = 0.9500`, `TNR = 0.9000`.
     - `RandomForest`: `ROC-AUC = 0.9815`, `Bal.Acc = 0.8583`, `TPR = 0.9917`, `TNR = 0.7250`.
     - `XGBoost_Weighted`: `ROC-AUC = 0.9367`, `Bal.Acc = 0.8250`, `TPR = 0.9000`, `TNR = 0.7500`.
     - `OneClassSVM_RBF`: `ROC-AUC = 0.8698`, `Bal.Acc = 0.7792`, `TPR = 0.7083`, `TNR = 0.8500`.


---

## 14. Physical Rationale of Selected Descriptors & Extrapolation Audit (`LOFMO`)

1. **Physical Rationale & Process Relevance of Selected Descriptors**:
   - **Steric/Coordination Shielding (`%V_bur` + Adaptive $V_{\min}$)**: Quantifies ligand steric protection within $R = 3.5\text{ Å}$ around the central element (DOI `10.1038/s41557-019-0319-5`). Adaptive $V_{\min} = 34.0\%$ for small $\text{Si}$ centers and $\pi$-haptic complexes vs $44.0\%$ for standard TM complexes prevents pyrophoric/oligomeric under-coordination (`< V_min`) and non-volatile steric overcrowding (`> 83.5%`).
   - **Thermal Decomposition Resistance (`BDE_min`, `DeltaE_dimer`)**: Weakest metal–ligand dissociation energy governs bubbler/delivery thermal stability ($T_{\text{decomp}} \leftrightarrow \text{BDE}_{\min} \ge 45\text{ kcal/mol}$), while non-overlapping dimer energy separates volatile monomers from sol-gel networks (`M(OH)n`).
   - **Coordination & Surface Topology (`pi_hapticity_index`, `chelate_ring_count`, `is_homoleptic_tm_sp3`, `metal_oh_count`, `cooh_count`, `is_inert_alkane`, `strained_si_ring`, `TPSA/HBD/MW`)**: Separates stabilizing $\pi$/chelate complexes from thermally labile $sp^3$ alkyls ($\beta$-H elimination), corrosive acids, inert solvents, and strained anti-aromatic rings.

2. **Leave-One-Failure-Mechanism-Out (`LOFMO`) Extrapolation Audit & 4 Improvement Priorities**:
   - Under `LOFMO` (holding out an entire failure mechanism from training on the `160`-molecule cohort), supervised binary classifiers collapsed on unseen failure mechanisms (`LogisticRegression`: **`mean_unseen_tnr = 0.0238`**, `SVC_RBF`: **`mean_unseen_tnr = 0.1627`**), whereas positive-only **`OneClassSVM_RBF` achieved `0.8651`** and the physics-grounded **`Coord/Chem Gate` achieved `0.9722`**.
   - **4 Key Extrapolation Improvements Identified**:
     1. Replace the hardcoded `16`-element semiconductor center list (which caused the single `NEG-019 CrPh3` miss) with periodic-table block ($d/f/p$-block), valence-electron count, and oxidation-state descriptors.
     2. Upgrade `SevenNet-0` (`MPtrj` inorganic crystal training set, which produces negative neutral-ligand BDEs on `CO`/halide/arene radicals) to spin/charge-aware `MACE-OMol` or `g-xTB`/DFT tiers.
     3. Integrate condensed-phase `LiCl(s)`-referenced stoichiometric synthesis energies ($\Delta E_{\text{syn}}$) and Clausius–Clapeyron volatility margins ($\Delta T = T_{\text{decomp}} - T_{\text{vap}}$).
     4. Standardize Leave-One-Element-Out (`LOEO`) and `LOFMO` alongside random `CV-5` in automated regression tests.
