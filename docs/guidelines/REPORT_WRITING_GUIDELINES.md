# Report Authoring, Terminology & Verification Guidelines

> **Scope**: Mandatory standard for all HTML reports ([`template.html`](../../template.html), [`generate_html.py`](../../generate_html.py), [`index.html`](../../index.html)), verification summaries ([`report.md`](../../report.md)), and audit documentation in this repository.

---

## 1. Epistemological & Tone Standard

1. **Neutral, Data-Driven Objectivity**:
   - Report empirical measurements, failure modes, and methodological boundaries concisely without promotional or hyperbolic adjectives (`완벽`, `치명적`, `재앙`, `극단적`, `혁신`, `유일한`, `전무후무`, `무비판`, `100% 보장`, etc.).
   - State limitations explicitly alongside positive metrics: unrelaxed single-point calculations vs. DFT geometry optimization/TS barriers, fallback geometries on disconnected SMILES, algebraic proxy terms vs. explicit quantum reactions, and small-sample CV/LOOCV generalization drops vs. in-sample accuracy.
2. **Evidence Grounding**:
   - All quantitative claims in prose, summary cards, and diagrams must match the latest JSON artifacts in [`database/`](../../database) dynamically rendered by [`generate_html.py`](../../generate_html.py).
   - Never inject target-label branching, molecule-ID hardcoding, or circular product memorization into scoring or comparison tables.

---

## 2. Declarative Prose Requirement (No Interrogative Phrasing)

1. **Prohibit Question-Form Headings and Sentences**:
   - Do not use interrogative sentences or rhetorical questions (`?`, `~할까요?`, `~는가?`, `~인가?`, `질문 핵심`, `Why ...?`) in section headers, table column titles, callouts, diagram labels, or body text.
2. **Convert All Analysis to Declarative Statements (`평서문`)**:
   - Replace question headers with direct thematic labels stating the finding or mechanism.
   - Example:
     - Disallowed: `왜 비환식 전구체는 통과하고, 샌드위치 착체는 큰 음수로 산출되는가?`
     - Required: `비환식 전구체와 Cp 샌드위치 착체 간 SYBA 점수 분화 원인`
     - Disallowed: `어디까지 정상 패스하는가?` / `어디서 실패·기각되는가?`
     - Required: `유효 통과 범위 (Valid Domain)` / `실패 원인 및 기각 한계 (Failure Mechanism)`

---

## 3. Terminology Policy: English-Only Statistical & Modeling Terms

All technical terminology belonging to **statistics, evaluation metrics, machine learning, and computational model architectures** must be written strictly in **English**. Do **not** translate them into Korean and do **not** attach Korean parenthetical glosses alongside the English terms.

### 3.1 Statistical & Evaluation Metrics (Use English Only)
- **Confusion Matrix & Outcomes**: `Confusion Matrix`, `True Positive (TP)`, `False Positive (FP)`, `False Negative (FN)`, `True Negative (TN)` *(never `혼동 행렬`, `진양성`, `위양성`, `위음성`, `진음성`)*.
- **Classification & Validation Metrics**: `Accuracy`, `Precision`, `Recall`, `Sensitivity`, `Specificity`, `F1-Score`, `AUROC`, `Pass Rate`, `Fail Rate`, `Rejection Rate`, `False Positive Rate (FPR)`, `False Negative Rate (FNR)`, `Majority Baseline`, `Error Rate` *(never `정확도`, `정밀도`, `재현율`, `민감도`, `특이도`, `통과율`, `기각률`, `오판율`, `차단율`, `다수 클래스 기준선`)*.
- **Validation Regimes & Diagnostics**: `In-Sample (All-Sample)`, `Cross-Validation`, `Stratified 5-Fold CV`, `Leave-One-Out CV (LOOCV)`, `Out-of-Fold`, `Overfitting`, `Underfitting`, `Data Leakage`, `Circular Memorization`, `Counterfactual Audit`, `Sensitivity Analysis`, `Sample Size (N)` *(never `교차 검증`, `과적합`, `데이터 누수`, `순환 기억`, `표본 수`)*.

### 3.2 Machine Learning & Modeling Terminology (Use English Only)
- **Model Families & Architectures**: `ML`, `Machine Learning`, `Deep Learning`, `DNN`, `GNN`, `GAT`, `XGBoost`, `Random Forest`, `Logistic Regression`, `Naive Bayes`, `MLP`, `Transformer`, `Autoregressive (AR)`, `Seq2Seq`, `Diffusion Model`, `GFlowNet`, `Generative Model`, `Surrogate Model`, `RL`, `CASP`, `MCTS`, `A* Search`, `MLIP` *(never `머신러닝`, `기계학습`, `심층 신경망`, `그래프 어텐션`, `랜덤 포레스트`, `나이브 베이지안/베이즈`, `생성 모델`, `확산 모델`, `자기회귀`, `트랜스포머`, `대리 모델`, `강화학습`, `몬테카를로 트리 탐색`, `트리 탐색`, `앙상블`, `분류기`, `회귀`)*.
- **Features, Embeddings & Optimization**: `Baseline`, `Base Model`, `Feature`, `Feature Set`, `Feature Importance`, `Gini Importance`, `Descriptor`, `Circular Fingerprint` (`Morgan Fingerprint`, `ECFP4`, `ECFP6`), `Fragment Frequency`, `One-Hot Encoding`, `Node Embedding`, `Latent Space`, `Tokenizer`, `Tokenization`, `Vocabulary`, `Out-of-Vocabulary (OOV)`, `Subword Splitting`, `Loss Function`, `Objective Function`, `Reward Function`, `Reward Hacking`, `Log-Odds Ratio`, `Regularization`, `Hyperparameter Tuning`, `Active Learning`, `Multi-Fidelity Gating`, `Fermi-Dirac Sigmoid`, `Smooth Penalty`, `Continuous Penalty`, `Proxy`, `Single-Point`, `Fallback` *(never `피처`, `피처 중요도`, `베이스라인`, `단편 빈도`, `지문`, `잠재 공간`, `토크나이저`, `토큰화`, `임베딩`, `손실 함수`, `목적함수`, `보상함수`, `승산비/로그공산비`, `능동 학습`, `다중정밀도`, `연속 페널티`, `시그모이드`, `과소평가/과대평가`)*.

### 3.3 Permitted Korean Context
- Natural Korean grammar particles, verbs, and semiconductor/chemistry domain prose remain in Korean (e.g., `반도체 전구체`, `결합 해리에너지`, `입체 차폐`, `이합체화`, `올리고머화`, `표면 흡착`, `양자화학 계산`, `유기 프로-리간드`, `판정 기준`).
- Do **not** wrap Korean chemistry nouns with redundant English synonyms when the term is purely domain-chemical rather than a statistical/model term, unless introducing a specific chemical property acronym once (e.g., `$\%V_{\text{bur}}$`, `$\text{BDE}_{\min}$`).

---

## 4. Technical Formatting & Rendering Integrity

1. **LaTeX / KaTeX Safety**:
   - In Python string literals inside [`generate_html.py`](../../generate_html.py), escape all LaTeX backslashes (`\\text`, `\\frac`, `\\le`, etc.) or use raw strings so control characters (`\t`, `\b`, `\f`, `\r`) are never written into [`index.html`](../../index.html).
2. **Mermaid Compatibility**:
   - Avoid literal `%` inside Mermaid node labels (write `percent`), avoid unbalanced brackets/HTML tags inside node strings, and verify rendering via headless Edge SVG inspection (`9/9` diagrams, `0` syntax errors).
3. **Repository & Execution Standards**:
   - Do not leave ad-hoc verification scripts in the workspace root; place durable scripts under [`scripts/benchmark/`](../../scripts/benchmark) or [`scripts/data/`](../../scripts/data) and run scratch inspections outside the repository.

---

## 5. Pre-Commit Verification Checklist

Every report update must pass all automated checks after the final edit:
1. `python generate_html.py` compiles [`index.html`](../../index.html) and [`synthetic_accessibility_report.html`](../../synthetic_accessibility_report.html) with 0 unresolved placeholders and 0 ASCII control characters in `$...$`.
2. Automated style scan confirms 0 interrogative prose lines (`?`), 0 prohibited hyperbolic phrases, and 0 Korean statistical/modeling glosses across [`template.html`](../../template.html) and [`generate_html.py`](../../generate_html.py).
3. Headless Edge Mermaid verification confirms `9/9` rendered SVGs with 0 syntax errors.
4. `python -m unittest discover tests` and `python scripts/benchmark/audit_reproduce_all.py` pass without error.
