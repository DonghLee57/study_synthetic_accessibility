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
3. **Algebraic Proxy & Disconnected SMILES Fallbacks (`run_phase2_sevenn_features.py`)**:
   - `e_pot_ev` and `e_coh_ev_per_atom` were evaluated via single-point `7net-0` on **unrelaxed** RDKit ETKDGv3 conformers (53/55 molecules; multi-fragment `Ru(EtCp)2` and `CCTBA` hit fallback constants `-100.0 eV`, `4.10 eV/atom`, `60.0 kcal/mol`, `15.0 kcal/mol`).
   - `delta_e_complex_kcal_mol` matched the deterministic algebraic formula $-25.0 + \text{steric} + \text{dimer} + \text{bde\_penalty}$ on **55/55 molecules** with zero salt-metathesis calculations performed in SevenNet. It has been excluded from supervised ML validation and analytical decision rules.
4. **Hardcoded Labels & ID Fallbacks in Analytical Rules (`run_model_tuning_comparison.py`)**:
   - Earlier headline scores (**96.3% precursor recovery**, **100.0% inhibitor pass**, **100.0% negative control rejection**) relied on branching by ground-truth `target_type`, hardcoding `base_ligand_sa = 1.5`, shifting the inhibitor cutoff to `4.2`, and explicitly hardcoding all 10 negative sample IDs/SMILES (`NEG-001`..`NEG-010`).
   - **Resolution**: Removed all label branching, constant ligand assumptions, and hardcoded ID checks. Re-evaluated uniformly across all 55 molecules (`tmc_sa_phys <= 4.0`).

---

## 2. Un-Hardcoded Analytical Screening Results ($N = 55$)

Using genuine `ligand_sa` plus continuous Fermi-Dirac smooth penalties on computed `E_coh`, `BDE_min`, `%V_bur` (restricted structurally to `center_atom` $\in$ metal/metalloid set), and `Delta E_dimer`, with uniform acceptance criterion $\text{TMC-SA}_{\text{phys}} \le 4.0$:

| Group / Evaluation Scope | Baseline (`SAScore <= 3.5` / Consensus) | Uniform Physics Only (`<= 4.0`) | Physics + General SMARTS Chem Screen |
| :--- | :---: | :---: | :---: |
| **Overall Accuracy ($N=55$)** | — | **0.836** (`TP=43, FN=2, FP=7, TN=3`) | **0.873** (`TP=43, FN=2, FP=5, TN=5`) |
| **Semiconductor Precursors ($N=27$)** | 3/27 (11.1%, FN 88.9%) | **25/27 (92.6%, FN 7.4%)** | **25/27 (92.6%, FN 7.4%)** |
| **AS-ALD Inhibitors ($N=18$)** | 15/18 (83.3%) | **18/18 (100.0%)** | **18/18 (100.0%)** |
| **Negative Controls Rejected ($N=10$)** | 0/10 rejected (100.0% consensus FP; 90.0% SCScore FP) | **3/10 rejected (70.0% FP)** (`NEG-004, 009, 010`) | **5/10 rejected (50.0% FP)** (`+ NEG-005, 006`) |

**Missed Negative Controls under Un-Hardcoded Single-Point Physics ($5/10$)**:
`NEG-001` ($\text{TiEt}_4$, unrelaxed stretch `64.5 kcal/mol`), `NEG-002` ($\text{ZrMe}_4$, `50.9 kcal/mol`), `NEG-003` ($\text{Silacyclobutadiene}$, rigid dimer `-0.85 kcal/mol`), `NEG-007` ($\text{Ti(OH)}_4$, dimer `-0.50 kcal/mol`), and `NEG-008` ($\text{CoEt}_4$, `3.61 <= 4.0`). This confirms that single-point MLIP descriptors on unrelaxed geometries cannot capture low-barrier kinetic decomposition ($\beta$-H elimination) or polycondensation without relaxed saddle-point transition states.

---

## 3. Machine Learning Validation: All-Sample vs. 5-Fold Stratified CV vs. LOOCV

Removing the algebraic complexation proxy and training unbiased classifiers on genuine Phase 1 (`fixed BR-SAScore` + `ligand_sa`) and Phase 2 (`7net_E_coh`, `7net_BDE_min`, `SambVca_Vbur`, `7net_dE_Dimer`, `Rg`, `V_mol`) features exposes substantial overfitting in In-Sample (All-Sample) evaluation:

### 3.1 Full Dataset ($N = 55$: 45 Valid Positive, 10 Negative Controls; Majority Baseline Acc = `0.818`)

| Classifier & Feature Set | All-Sample (In-Sample) | Stratified 5-Fold CV | Leave-One-Out CV (LOOCV) |
| :--- | :---: | :---: | :---: |
| **RF (`default`) — Base Only (10 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.873`, Spec `0.500`, AUROC `0.791` | **Acc `0.873`, Spec `0.500`, AUROC `0.779`** |
| **RF (`default`) — Phys Only (6 feat)** | Acc `0.964`, AUROC `1.000` | Acc `0.836`, Spec `0.100`, AUROC `0.511` | Acc `0.782`, Spec `0.000`, AUROC `0.473` |
| **RF (`default`) — Combined (16 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.800`, Spec `0.000`, AUROC `0.796` | Acc `0.782`, Spec `0.100`, AUROC `0.721` |
| **RF (`balanced`) — Base Only (10 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.855`, Spec `0.400`, AUROC `0.772` | **Acc `0.855`, Spec `0.300`, AUROC `0.778`** |
| **RF (`balanced`) — Combined (16 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.818`, Spec `0.100`, AUROC `0.756` | Acc `0.800`, Spec `0.000`, AUROC `0.680` |
| **LogReg (`L2, balanced`) — Base Only** | Acc `0.745`, AUROC `0.829` | Acc `0.691`, Spec `0.400`, AUROC `0.493` | Acc `0.618`, Spec `0.500`, AUROC `0.484` |
| **LogReg (`L2, balanced`) — Combined** | Acc `0.836`, AUROC `0.864` | Acc `0.636`, Spec `0.300`, AUROC `0.462` | Acc `0.655`, Spec `0.400`, AUROC `0.458` |

### 3.2 Precursor-Domain Subset ($N = 35$: 27 Precursors, 8 Organometallic/Silane Negatives; Majority Acc = `0.771`)

| Classifier & Feature Set | All-Sample (In-Sample) | Stratified 5-Fold CV | Leave-One-Out CV (LOOCV) |
| :--- | :---: | :---: | :---: |
| **RF (`default`) — Base Only (10 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.829`, Spec `0.500`, AUROC `0.850` | **Acc `0.829`, Spec `0.500`, AUROC `0.833`** |
| **RF (`default`) — Phys Only (6 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.714`, Spec `0.000`, AUROC `0.407` | Acc `0.714`, Spec `0.125`, AUROC `0.472` |
| **RF (`default`) — Combined (16 feat)** | Acc `1.000`, AUROC `1.000` | Acc `0.743`, Spec `0.125`, AUROC `0.778` | Acc `0.686`, Spec `0.000`, AUROC `0.769` |

---

## 4. Scientific Interpretation & Conclusions

1. **In-Sample vs. Out-of-Sample Divergence**: While Random Forest memorizes all 55 training points (`Acc 1.000, AUROC 1.000`) with physical features contributing `33.3%` Gini importance in-sample, appending 6 unrelaxed physical features to 10 baseline descriptors **degrades both 5-Fold CV (`0.873 -> 0.800`) and LOOCV (`0.873 -> 0.782`, Specificity `0.500 -> 0.100`)**.
2. **Root Causes**:
   - High feature-to-minority-sample ratio ($16$ features vs. $N_{\text{neg}} = 8 \sim 10$ distributed across distinct mechanisms such as $\beta$-H elimination, steric crowding, oligomerization, acid etching, and alkane inertness).
   - Single-point `7net-0` evaluation on unrelaxed ETKDGv3 conformers introduces rigid-clash noise into `Delta E_dimer` (up to $+1074.7\text{ kcal/mol}$) and lacks relaxed kinetic barrier resolution.
3. **Operational Recommendation**: SevenNet + SambVca descriptors should currently operate as explicit, domain-bounded physical guardrails (`%V_bur`, `BDE_min`, `E_coh`) alongside SMARTS reactivity screens rather than unconstrained black-box ML inputs until a larger DFT-relaxed organometallic failure corpus ($N_{\text{neg}} \ge 100$) with explicit TS barriers ($\Delta G^\ddagger$) is curated.
