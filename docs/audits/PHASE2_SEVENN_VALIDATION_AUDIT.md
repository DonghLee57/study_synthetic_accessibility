# Phase 2 SevenNet Physical Feature & Model Tuning Validation Audit

> **Audit Scope**: Rigorous code, data-lineage, and statistical validation review of `scripts/benchmark/run_phase1_baseline.py`, `scripts/benchmark/run_phase2_sevenn_features.py`, and `scripts/benchmark/run_model_tuning_comparison.py`, comparing All-Sample, 5-Fold Stratified CV, and Leave-One-Out Cross-Validation (LOOCV).

---

## 1. Executive Audit Summary

A complete code-level and numerical audit of the Phase 1/2 pipeline identified **four critical defects** in the prior iteration and resolved them with reproducible, un-hardcoded runs:

1. **Silent `BR-SAScore` Failure in Phase 1 Baseline (`run_phase1_baseline.py`)**:
   - Calling `BRSAScorer.calculateScore(mol)` with an RDKit `Mol` object and calling `round()` on the returned `(score, contrib)` tuple triggered a caught `TypeError` across **55/55 molecules**, assigning a constant fallback `br_sascore = 9.99` (producing `0.0000` feature importance).
   - **Resolution**: Fixed to pass canonical SMILES (`smi`), unpack `(score, _)`, and additionally compute true pro-ligand cleaved `ligand_sa` via RDKit metal/silicon bond cleavage ($\min(\text{BR-SAScore}(HL), \text{SAScore}(HL))$).
2. **Stale JSON Snapshot Mismatch (`37` Discrepancies)**:
   - `database/external/sevenn_physical_features.json` had been updated by a full `1290.8 s` recalculation run (commit `0ce6879`), whereas `database/benchmark_phase2_tuned_comparison.json` still held values from an earlier `426.5 s` run (`e8745ce`, 37 mismatched entries across `bde_min` and `delta_e_dimer`), and `run_phase2_sevenn_features.py` still carried an `if vbur >= 65.0` linear bypass.
   - **Resolution**: Synchronized `run_phase2_sevenn_features.py` with the direct SevenNet dimer calculation logic and regenerated `database/benchmark_phase2_tuned_comparison.json` from the latest `1290.8 s` physical feature dataset.
3. **SMILES Stoichiometry Fix & Upgraded Relaxation / Dimer / BDE Algorithms (`database/precursors.json`, `run_phase2_sevenn_features.py`)**:
   - **SMILES & Formula Parity**: Fixed `CCTBA (PREC_MET_05)` from disconnected `[Co](=O)(=O)=O` oxo-fragments (`C12H9Co2O6` typo missing terminal alkyne H) to coordinated hexacarbonyl dimetallatetrahedrane (`C12H10Co2O6`, CAS `56792-69-9`), connected `Ru(EtCp)2 (PREC_MET_03)` (`C14H18Ru`), corrected `Cp(NMe2)3Zr`/`CpTi(NMe2)3` formulas (`H24 -> H23`), and fixed `Cp2AlH (NEG-010)` SMILES (`C10H11Al`). **Result: `0 / 55` fallbacks (`-100.0 eV` / `4.10 eV/atom` eliminated).**
   - **LBFGS Relaxation**: Replaced unrelaxed ETKDGv3 single-point evaluation with ASE `LBFGS` geometry optimization (`fmax=0.04 eV/A`, up to 100 steps for monomers, 35 steps for dimers) on CPU (`torch.set_num_threads(4)`), followed by SambVca `%V_bur` evaluation on the relaxed coordinates.
   - **Non-Overlapping Dimer Separation**: Replaced buggy `0.5 * z_span + 2.8` shift (which caused nuclear interpenetration up to `+1074.7 kcal/mol` on `WPh5`) with principal-axis alignment + $\Delta z_{\text{shift}} = z_{\text{span}} + d_{\text{vdW}}$ ($d_{\text{vdW}} = 2.8\,\text{\AA}$) guaranteed non-overlapping gap (`dimer_init_gap_ang >= 2.8 A` verified across all 55 molecules) + dimer LBFGS relaxation.
   - **Fragment Dissociation & Subtree Ring Stretch BDE**: Replaced single-atom (`a2`) displacement (which caused intramolecular compression artifacts such as `DIPAS 37.9 kcal/mol`) with true relaxed separated fragment dissociation $\text{BDE} = E_{\text{rel}}(ML_{n-1}) + E_{\text{rel}}(L) - E_{\text{rel}}(ML_n)$ for acyclic bonds, and exocyclic subtree displacement + `FixAtoms([a1, a2])` framework relaxation for cyclic bonds.
   - `delta_e_complex_kcal_mol` remains flagged as an algebraic proxy excluded from supervised ML validation and analytical rules.
4. **Hardcoded Labels & ID Fallbacks in Analytical Rules (`run_model_tuning_comparison.py`)**:
   - Removed all label branching, constant ligand assumptions, and hardcoded ID checks. Re-evaluated uniformly across all 55 molecules (`tmc_sa_phys <= 4.0`).

---

## 2. Un-Hardcoded Analytical Screening Results ($N = 55$, Relaxed `7net-0`)

Using genuine `ligand_sa` plus continuous Fermi-Dirac smooth penalties on LBFGS-relaxed `E_coh`, `BDE_min`, `%V_bur` (restricted structurally to `center_atom` $\in$ metal/metalloid set), and non-overlapping `Delta E_dimer`, with uniform acceptance criterion $\text{TMC-SA}_{\text{phys}} \le 4.0$:

| Group / Evaluation Scope | Baseline (`SAScore <= 3.5` / Consensus) | Uniform Physics Only (`<= 4.0`) | Physics + General SMARTS Chem Screen |
| :--- | :---: | :---: | :---: |
| **Overall Accuracy ($N=55$)** | — | **0.764** (`TP=38, FN=7, FP=6, TN=4`) | **0.800** (`TP=38, FN=7, FP=4, TN=6`) |
| **Semiconductor Precursors ($N=27$)** | 2/27 (7.4%, FN 92.6%) | **22/27 (81.5%, FN 18.5%)** | **22/27 (81.5%, FN 18.5%)** |
| **AS-ALD Inhibitors ($N=18$)** | 15/18 (83.3%) | **16/18 (88.9%)** | **16/18 (88.9%)** |
| **Negative Controls Rejected ($N=10$)** | 1/10 rejected (90.0% consensus FP; 90.0% SCScore FP) | **4/10 rejected (60.0% FP)** (`NEG-003, 008, 009, 010`) | **6/10 rejected (40.0% FP)** (`+ NEG-005, 006`) |

**Key Algorithmic Shifts & Remaining Failure Modes**:
- **Resolved False Negatives/Collisions**: `DIPAS` BDE rose from collision-compressed `37.9` to `55.5 kcal/mol` (`DIPAS`/`TSA` now fail only due to universal `%V_bur < 44%` dimer penalty on small volatile hydrides); `WPh5` dimer energy normalized from `+616.0` (single-point) / `+1074.7` to `0.0 kcal/mol`.
- **Newly Captured Negative Controls**: Strained `NEG-003` (`Silacyclobutadiene`, relaxed ring BDE `15.5 kcal/mol`, `TMC-SA = 5.43`) and `NEG-008` (`CoEt4`, fragment BDE `35.3 kcal/mol`, `TMC-SA = 4.48`) are now correctly rejected (`Spec 0.500 -> 0.600`).
- **Remaining Precursor/Inhibitor FNs (`5/27` + `2/18`) & Missed Negatives (`4/10`)**: Small volatile silanes (`DIPAS`, `TSA`) and $\eta^1$-embedded `Ru(EtCp)2` (`%V_bur 36.9%`) trigger the universal metal under-coordination dimer penalty; `CCTBA` (`11.1 kcal/mol`), EUV `tBuSn(OiPr)3` (`3.8 kcal/mol`), and triazoles `BTA`/`TTA` (`34.8`/`34.6 kcal/mol`) reflect facile CO/radical/N-N bond cleavage. Meanwhile, `NEG-001` (`TiEt4`, $\beta$-H elimination), `NEG-002` (`ZrMe4`), `NEG-004` (`WPh5`, LBFGS-relaxed `%V_bur 79.1%`), and `NEG-007` (`Ti(OH)4`, `TMC-SA 3.70`) require kinetic transition-state barriers ($\Delta G^\ddagger$) and element-specific thresholds.

---

## 3. Machine Learning Validation: All-Sample vs. 5-Fold Stratified CV vs. LOOCV

Removing the algebraic complexation proxy and training unbiased classifiers on genuine Phase 1 (`fixed BR-SAScore` + `ligand_sa`) and LBFGS-relaxed Phase 2 (`7net_E_coh`, `7net_BDE_min`, `SambVca_Vbur`, `7net_dE_Dimer`, `Rg`, `V_mol`) features yields substantially improved ranking discrimination (LOOCV AUROC):

### 3.1 Full Dataset ($N = 55$: 45 Valid Positive, 10 Negative Controls; Majority Baseline Acc = `0.818`)

| Classifier & Feature Set | All-Sample (In-Sample) | Stratified 5-Fold CV | Leave-One-Out CV (LOOCV) |
| :--- | :---: | :---: | :---: |
| **RF (`default`) — Base Only (10 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.836`, Spec `0.400`, AUROC `0.794` | Acc `0.855`, Spec `0.400`, AUROC `0.782` |
| **RF (`default`) — Phys Only (6 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.836`, Spec `0.300`, AUROC `0.659` | Acc `0.800`, Spec `0.100`, AUROC `0.662` |
| **RF (`default`) — Combined (16 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.818`, Spec `0.100`, **AUROC `0.833`** | Acc `0.836`, Spec `0.200`, **AUROC `0.849`** |
| **RF (`balanced`) — Base Only (10 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.873`, Spec `0.500`, AUROC `0.773` | **Acc `0.855`, Spec `0.300`, AUROC `0.778`** |
| **RF (`balanced`) — Combined (16 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.818`, Spec `0.100`, AUROC `0.828` | Acc `0.855`, Spec `0.200`, AUROC `0.813` |
| **LogReg (`L2, balanced`) — Base Only** | Acc `0.691`, AUROC `0.787` | Acc `0.673`, Spec `0.400`, AUROC `0.449` | Acc `0.600`, Spec `0.500`, AUROC `0.422` |
| **LogReg (`L2, balanced`) — Combined** | Acc `0.782`, AUROC `0.851` | Acc `0.655`, Spec `0.400`, AUROC `0.398` | Acc `0.618`, Spec `0.400`, AUROC `0.367` |

### 3.2 Precursor-Domain Subset ($N = 35$: 27 Precursors, 8 Organometallic/Silane Negatives; Majority Acc = `0.771`)

| Classifier & Feature Set | All-Sample (In-Sample) | Stratified 5-Fold CV | Leave-One-Out CV (LOOCV) |
| :--- | :---: | :---: | :---: |
| **RF (`default`) — Base Only (10 feat)** | Acc `1.000`, AUROC `1.000` | **Acc `0.857`, Spec `0.625`, AUROC `0.847`** | **Acc `0.800`, Spec `0.375`, AUROC `0.836`** |
| **RF (`default`) — Phys Only (6 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.800`, Spec `0.375`, AUROC `0.817` | Acc `0.743`, Spec `0.250`, AUROC `0.782` |
| **RF (`default`) — Combined (16 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.771`, Spec `0.375`, AUROC `0.831` | Acc `0.743`, Spec `0.250`, AUROC `0.836` |

---

## 4. Scientific Interpretation & Conclusions

1. **Out-of-Fold AUROC Improvement vs. Thresholded Overfitting**: Fixing SMILES stoichiometry (`0/55` fallbacks) and computing LBFGS-relaxed geometries, non-overlapping dimers, and separated fragment BDE raised `Phys-Only` LOOCV AUROC (`0.473 -> 0.662`) and **elevated `Combined` LOOCV AUROC (`0.721 -> 0.849`), outperforming `Base-Only` (`0.782`)** with `40.6%` physical feature Gini importance (`#1 7net_BDE_min = 17.33%`). Nevertheless, default probability cutoff (`0.50`) still suffers minority-class recall collapse under $N_{\text{neg}} = 10$ imbalance (`LOOCV Spec 0.200`).
2. **Physical Limitations Exposed Without Noise**: Ground-state MLIP descriptors (`E_coh`, `BDE_min`, `%V_bur`, `Delta E_dimer`) alone cannot identify kinetic $\beta$-H elimination (`TiEt4`), thermal detonation (`ZrMe4`), or sol-gel hydrolysis without DFT saddle-point activation barriers ($\Delta G^\ddagger$), and applying metal under-coordination dimer penalties indiscriminately across covalently bonded Si(IV) hydrides/aminosilanes causes false negatives (`DIPAS`, `TSA`).
3. **Operational Recommendation & Test Coverage**: SevenNet + SambVca descriptors should act as **element-class-conditioned physical guardrails** (`%V_bur`, `BDE_min`, `E_coh`) paired with SMARTS reactivity filters while building a DFT transition-state failure corpus ($N_{\text{neg}} \ge 100$). Automated verification now spans **12 reproducible unit tests** in `tests/test_benchmark_reproducibility.py` covering Phase 1 (`N=55` formula/score parity), `sevenn_physical_features.json` zero-fallback & `gap >= 2.8 A` invariants, and exact `run_model_tuning_comparison.py` CV/LOOCV reproduction.
