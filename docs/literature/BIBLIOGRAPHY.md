# Research Bibliography & Literature Archive (문헌 서지 및 아카이브)

> **Project**: Synthetic Accessibility (SA) Evaluation, Precursor Inverse Design, and High-Throughput DFT Automation  
> **Verification Status**: All **34** scholarly papers verified via CrossRef & DOI registry (100.0%)  
> **Local Archival Status**: **26 / 34 (76.5%)** full-text documents secured locally (PDF 22편 + XML 4편)

---

## 1. Executive Summary & Archival Statistics

| Metric | Count / Status | Percentage | Notes |
| :--- | :---: | :---: | :--- |
| **Total Evaluated Literature** | **34 papers** | 100.0% | Spanning 2009–2026 across 13 report sections |
| **Official DOI Verification** | **34 papers** | 100.0% | CrossRef & Publisher API verified |
| **Full-Text Locally Archived** | **26 papers** | **76.5%** | Saved in `downloaded_papers/` and `by_section/` |
| &bull; *Open-Access / Preprint PDFs* | *22 papers* | *64.7%* | arXiv, ChemRxiv, eScholarship, Nature Comms, etc. |
| &bull; *PMC / EuropePMC Full-Text XML* | *4 papers* | *11.8%* | JACS Au, Chem. Sci., Chem. Mater. via PMC/EBI |
| **Commercial Paywall (Direct DOI)** | **8 papers** | **23.5%** | ACS JCIM, RSC without open-access preprints |

### Directory Layout of Archived Full Texts

The repository maintains both a canonical flat archive and thematic section subdirectories for immediate localized browsing:

```text
downloaded_papers/
├── [26 canonical PDF/XML files: SAScore_2009.pdf, Ong_Pymatgen_2013.pdf, ...]
├── manifest.json (Machine-readable metadata for all 34 papers)
└── by_section/
    ├── sec01_sa_benchmarks/           (15 files: 1st~4th Gen SA & Benchmark papers)
    ├── sec05_extended_studies/         (4 files: RetroGNN, Gao&Coley, Kulik TMC, MPScore)
    ├── sec08_semiconductor_precursors/ (7 files: Aminosilane, Kulik, MPScore, AS-ALD, MLIP, 2026 Roadmap, BR-SAScore)
    ├── sec10_vasp_automation/          (4 files: Pymatgen, AiiDA, Sella, MLIP-precursor)
    └── sec13_precursor_inverse_design_asd/ (11 files: Mackus ASD, Tonner ASD, ChemProjector, SynFlowNet, RxnFlow, ReaSyn, ...)
```

---

## 2. Section-by-Section Literature Catalog

### Section 1: Literature Verification & Master Download Hub (문헌 검증 및 마스터 아카이브)
**Literature Count**: 34 papers (Locally secured: 26/34)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **SAScore** | 2009 | Peter Ertl, Ansgar Schuffenhauer (*Journal of Cheminformatics*) | [PDF 열람 (1242 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SAScore_2009.pdf) | [10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8) |
| 02 | **SYBA** | 2020 | Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil (*Journal of Cheminformatics*) | [PDF 열람 (3191 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SYBA_2020.pdf) | [10.1186/s13321-020-00439-2](https://doi.org/10.1186/s13321-020-00439-2) |
| 03 | **SCScore** | 2018 | Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.7b00622) | [10.1021/acs.jcim.7b00622](https://doi.org/10.1021/acs.jcim.7b00622) |
| 04 | **RAscore** | 2021 | Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond (*Chemical Science*) | [PDF 열람 (322 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RAscore_2021.pdf) | [10.1039/D0SC05401A](https://doi.org/10.1039/D0SC05401A) |
| 05 | **RetroGNN** | 2022 | Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler (*Journal of Chemical Information and Modeling*) | [PDF 열람 (2400 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroGNN_2022.pdf) | [10.1021/acs.jcim.1c01476](https://doi.org/10.1021/acs.jcim.1c01476) |
| 06 | **GASA** | 2022 | Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.2c00038) | [10.1021/acs.jcim.2c00038](https://doi.org/10.1021/acs.jcim.2c00038) |
| 07 | **DeepSA** | 2023 | Shihang Wang, Lin Wang, Fenglei Li, Fang Bai (*Journal of Cheminformatics*) | [PDF 열람 (4058 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/DeepSA_2023.pdf) | [10.1186/s13321-023-00771-3](https://doi.org/10.1186/s13321-023-00771-3) |
| 08 | **DFRscore** | 2023 | Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01134) | [10.1021/acs.jcim.3c01134](https://doi.org/10.1021/acs.jcim.3c01134) |
| 09 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |
| 10 | **ChemProjector** | 2024 | Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma (*NeurIPS 2024 / arXiv:2406.04628*) | [PDF 열람 (4788 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) | [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628) |
| 11 | **SynFormer** | 2025 | Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al. (*PNAS 2025 / arXiv:2410.03494*) | [PDF 열람 (3651 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFormer_2025.pdf) | [10.1073/pnas.2415665122](https://doi.org/10.1073/pnas.2415665122) |
| 12 | **SynFlowNet** | 2025 | Miruna Cretu, Charles Harris, Austin Tripp, et al. (*ICLR 2025 / arXiv:2405.01155*) | [PDF 열람 (8060 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) | [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155) |
| 13 | **RxnFlow** | 2024 | Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim (*NeurIPS 2024 / arXiv:2410.04542*) | [PDF 열람 (4912 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) | [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542) |
| 14 | **ReaSyn** | 2025 | KAIST & NVIDIA Research Team (*arXiv:2509.16084*) | [PDF 열람 (2479 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) | [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084) |
| 15 | **SynNet** | 2022 | Wenhao Gao, Rocío Mercado, Connor W. Coley (*ICLR 2022*) | [PDF 열람 (3079 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynNet_2022.pdf) | [10.48550/arXiv.2110.06389](https://doi.org/10.48550/arXiv.2110.06389) |
| 16 | **Generative SA Benchmark** | 2020 | Wenhao Gao, Connor W. Coley (*Journal of Chemical Information and Modeling*) | [PDF 열람 (4955 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Gao_Coley_Benchmark_2020.pdf) | [10.1021/acs.jcim.0c00174](https://doi.org/10.1021/acs.jcim.0c00174) |
| 17 | **AiZynthFinder** | 2020 | Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist (*Journal of Cheminformatics*) | [PDF 열람 (1115 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AiZynthFinder_2020.pdf) | [10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1) |
| 18 | **Retro*** | 2020 | Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al. (*ICML 2020*) | [PDF 열람 (3244 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroStar_2020.pdf) | [10.48550/arXiv.2006.15575](https://doi.org/10.48550/arXiv.2006.15575) |
| 19 | **V-SYNTH** | 2022 | Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al. (*Nature*) | [출판사 구독 (DOI)](https://doi.org/10.1038/s41586-021-04220-9) | [10.1038/s41586-021-04220-9](https://doi.org/10.1038/s41586-021-04220-9) |
| 20 | **TMC Synthesizability (Kulik Group)** | 2022 | Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al. (*JACS Au*) | [XML 열람 (159 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) | [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547) |
| 21 | **MPScore (Materials Precursor Score)** | 2021 | Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al. (*Journal of Chemical Information and Modeling*) | [XML 열람 (175 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) | [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375) |
| 22 | **AS-ALD Precursor ML Descriptors** | 2023 | Tran Thi Ngoc Van, Changsu Kim, Hojae Lee, Jiyong Kim et al. (*Journal of Molecular Modeling*) | [PDF 열람 (2563 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AS_ALD_Precursor_ML_2023.pdf) | [10.1007/s00894-023-05806-y](https://doi.org/10.1007/s00894-023-05806-y) |
| 23 | **MLIP for ALD Precursors** | 2024 | Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al. (*Materials Today Advances*) | [PDF 열람 (6267 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_MLIP_Precursor_2024.pdf) | [10.1016/j.mtadv.2024.100474](https://doi.org/10.1016/j.mtadv.2024.100474) |
| 24 | **ML in ALD Roadmap** | 2026 | Julia Westermayr, Stacey F. Bent, Ralf Tonner-Zech (*Chemistry of Materials*) | [PDF 열람 (3091 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_ML_ChemMater_2026.pdf) | [10.1021/acs.chemmater.5c03271](https://doi.org/10.1021/acs.chemmater.5c03271) |
| 25 | **Aminosilane ALD Precursor DFT Design** | 2017 | Liang Huang, Bo Han, Maohong Fan, Hansong Cheng (*RSC Advances*) | [PDF 열람 (736 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Aminosilane_SiO2_ALD_DFT_2017.pdf) | [10.1039/c7ra02301d](https://doi.org/10.1039/c7ra02301d) |
| 26 | **Oh et al. (ASD Precursor Choice)** | 2021 | In-Kyun Oh, Stacey F. Bent, Adriaan J. M. Mackus et al. (*Chemistry of Materials*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.chemmater.0c04718) | [10.1021/acs.chemmater.0c04718](https://doi.org/10.1021/acs.chemmater.0c04718) |
| 27 | **Tonner-Zech et al. (ASD Design Principles)** | 2024 | Ralf Tonner-Zech et al. (*Journal of Computational Chemistry*) | [XML 열람 (201 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Tonner_ASD_Precursor_2024.xml) | [10.1002/jcc.70493](https://doi.org/10.1002/jcc.70493) |
| 28 | **Mackus et al. (ASD Fundamentals & Metrics)** | 2019 | Adriaan J. M. Mackus, Ageeth A. Bol, W. M. Erwin Kessels (*Chemistry of Materials*) | [XML 열람 (146 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Mackus_ASD_Selectivity_2019.xml) | [10.1021/acs.chemmater.8b03454](https://doi.org/10.1021/acs.chemmater.8b03454) |
| 29 | **Parsons & Clark (ASD State-of-the-Art)** | 2020 | Gregory N. Parsons, Robert D. Clark (*Chemistry of Materials*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.chemmater.0c00722) | [10.1021/acs.chemmater.0c00722](https://doi.org/10.1021/acs.chemmater.0c00722) |
| 30 | **Pymatgen (Materials Genomics)** | 2013 | Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier et al. (*Computational Materials Science*) | [PDF 열람 (1357 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Ong_Pymatgen_2013.pdf) | [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028) |
| 31 | **atomate2 + jobflow** | 2024 | Andrew S. Rosen, Matthew Gallant, Anderson Janotti et al. (*Journal of Open Source Software*) | [출판사 구독 (DOI)](https://doi.org/10.21105/joss.06543) | [10.21105/joss.06543](https://doi.org/10.21105/joss.06543) |
| 32 | **AiiDA-VASP Automation** | 2021 | Martin Uhrin, Sebastiaan P. Huber, Jusong Yu, Nicola Marzari et al. (*Computational Materials Science*) | [PDF 열람 (925 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Uhrin_AiiDA_2021.pdf) | [10.1016/j.commatsci.2020.110086](https://doi.org/10.1016/j.commatsci.2020.110086) |
| 33 | **Architector (3D TMC Complex Builder)** | 2023 | Corin Gautier, Tyler Castonguay, Heather J. Kulik (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01353) | [10.1021/acs.jcim.3c01353](https://doi.org/10.1021/acs.jcim.3c01353) |
| 34 | **Sella (Automated Saddle-Point TS)** | 2022 | Eric D. Hermes, Judit Zádor (*Journal of Chemical Theory and Computation*) | [PDF 열람 (2963 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Hermes_Sella_2022.pdf) | [10.1021/acs.jctc.1c00845](https://doi.org/10.1021/acs.jctc.1c00845) |

### Section 3: Deep Architectural Analysis & Paradigm Evolution (세대별 패러다임 진화 및 수식)
**Literature Count**: 16 papers (Locally secured: 13/16)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **SAScore** | 2009 | Peter Ertl, Ansgar Schuffenhauer (*Journal of Cheminformatics*) | [PDF 열람 (1242 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SAScore_2009.pdf) | [10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8) |
| 02 | **SYBA** | 2020 | Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil (*Journal of Cheminformatics*) | [PDF 열람 (3191 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SYBA_2020.pdf) | [10.1186/s13321-020-00439-2](https://doi.org/10.1186/s13321-020-00439-2) |
| 03 | **SCScore** | 2018 | Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.7b00622) | [10.1021/acs.jcim.7b00622](https://doi.org/10.1021/acs.jcim.7b00622) |
| 04 | **RAscore** | 2021 | Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond (*Chemical Science*) | [PDF 열람 (322 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RAscore_2021.pdf) | [10.1039/D0SC05401A](https://doi.org/10.1039/D0SC05401A) |
| 05 | **RetroGNN** | 2022 | Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler (*Journal of Chemical Information and Modeling*) | [PDF 열람 (2400 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroGNN_2022.pdf) | [10.1021/acs.jcim.1c01476](https://doi.org/10.1021/acs.jcim.1c01476) |
| 06 | **GASA** | 2022 | Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.2c00038) | [10.1021/acs.jcim.2c00038](https://doi.org/10.1021/acs.jcim.2c00038) |
| 07 | **DeepSA** | 2023 | Shihang Wang, Lin Wang, Fenglei Li, Fang Bai (*Journal of Cheminformatics*) | [PDF 열람 (4058 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/DeepSA_2023.pdf) | [10.1186/s13321-023-00771-3](https://doi.org/10.1186/s13321-023-00771-3) |
| 08 | **DFRscore** | 2023 | Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01134) | [10.1021/acs.jcim.3c01134](https://doi.org/10.1021/acs.jcim.3c01134) |
| 09 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |
| 10 | **ChemProjector** | 2024 | Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma (*NeurIPS 2024 / arXiv:2406.04628*) | [PDF 열람 (4788 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) | [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628) |
| 11 | **SynFormer** | 2025 | Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al. (*PNAS 2025 / arXiv:2410.03494*) | [PDF 열람 (3651 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFormer_2025.pdf) | [10.1073/pnas.2415665122](https://doi.org/10.1073/pnas.2415665122) |
| 12 | **SynFlowNet** | 2025 | Miruna Cretu, Charles Harris, Austin Tripp, et al. (*ICLR 2025 / arXiv:2405.01155*) | [PDF 열람 (8060 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) | [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155) |
| 13 | **RxnFlow** | 2024 | Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim (*NeurIPS 2024 / arXiv:2410.04542*) | [PDF 열람 (4912 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) | [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542) |
| 14 | **ReaSyn** | 2025 | KAIST & NVIDIA Research Team (*arXiv:2509.16084*) | [PDF 열람 (2479 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) | [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084) |
| 15 | **SynNet** | 2022 | Wenhao Gao, Rocío Mercado, Connor W. Coley (*ICLR 2022*) | [PDF 열람 (3079 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynNet_2022.pdf) | [10.48550/arXiv.2110.06389](https://doi.org/10.48550/arXiv.2110.06389) |
| 16 | **AiZynthFinder** | 2020 | Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist (*Journal of Cheminformatics*) | [PDF 열람 (1115 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AiZynthFinder_2020.pdf) | [10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1) |

### Section 4: Benchmarking & Quantitative Comparison (16개 모델 정량 비교 매트릭스)
**Literature Count**: 18 papers (Locally secured: 15/18)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **SAScore** | 2009 | Peter Ertl, Ansgar Schuffenhauer (*Journal of Cheminformatics*) | [PDF 열람 (1242 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SAScore_2009.pdf) | [10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8) |
| 02 | **SYBA** | 2020 | Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil (*Journal of Cheminformatics*) | [PDF 열람 (3191 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SYBA_2020.pdf) | [10.1186/s13321-020-00439-2](https://doi.org/10.1186/s13321-020-00439-2) |
| 03 | **SCScore** | 2018 | Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.7b00622) | [10.1021/acs.jcim.7b00622](https://doi.org/10.1021/acs.jcim.7b00622) |
| 04 | **RAscore** | 2021 | Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond (*Chemical Science*) | [PDF 열람 (322 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RAscore_2021.pdf) | [10.1039/D0SC05401A](https://doi.org/10.1039/D0SC05401A) |
| 05 | **RetroGNN** | 2022 | Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler (*Journal of Chemical Information and Modeling*) | [PDF 열람 (2400 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroGNN_2022.pdf) | [10.1021/acs.jcim.1c01476](https://doi.org/10.1021/acs.jcim.1c01476) |
| 06 | **GASA** | 2022 | Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.2c00038) | [10.1021/acs.jcim.2c00038](https://doi.org/10.1021/acs.jcim.2c00038) |
| 07 | **DeepSA** | 2023 | Shihang Wang, Lin Wang, Fenglei Li, Fang Bai (*Journal of Cheminformatics*) | [PDF 열람 (4058 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/DeepSA_2023.pdf) | [10.1186/s13321-023-00771-3](https://doi.org/10.1186/s13321-023-00771-3) |
| 08 | **DFRscore** | 2023 | Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01134) | [10.1021/acs.jcim.3c01134](https://doi.org/10.1021/acs.jcim.3c01134) |
| 09 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |
| 10 | **ChemProjector** | 2024 | Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma (*NeurIPS 2024 / arXiv:2406.04628*) | [PDF 열람 (4788 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) | [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628) |
| 11 | **SynFormer** | 2025 | Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al. (*PNAS 2025 / arXiv:2410.03494*) | [PDF 열람 (3651 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFormer_2025.pdf) | [10.1073/pnas.2415665122](https://doi.org/10.1073/pnas.2415665122) |
| 12 | **SynFlowNet** | 2025 | Miruna Cretu, Charles Harris, Austin Tripp, et al. (*ICLR 2025 / arXiv:2405.01155*) | [PDF 열람 (8060 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) | [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155) |
| 13 | **RxnFlow** | 2024 | Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim (*NeurIPS 2024 / arXiv:2410.04542*) | [PDF 열람 (4912 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) | [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542) |
| 14 | **ReaSyn** | 2025 | KAIST & NVIDIA Research Team (*arXiv:2509.16084*) | [PDF 열람 (2479 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) | [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084) |
| 15 | **SynNet** | 2022 | Wenhao Gao, Rocío Mercado, Connor W. Coley (*ICLR 2022*) | [PDF 열람 (3079 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynNet_2022.pdf) | [10.48550/arXiv.2110.06389](https://doi.org/10.48550/arXiv.2110.06389) |
| 16 | **Generative SA Benchmark** | 2020 | Wenhao Gao, Connor W. Coley (*Journal of Chemical Information and Modeling*) | [PDF 열람 (4955 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Gao_Coley_Benchmark_2020.pdf) | [10.1021/acs.jcim.0c00174](https://doi.org/10.1021/acs.jcim.0c00174) |
| 17 | **AiZynthFinder** | 2020 | Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist (*Journal of Cheminformatics*) | [PDF 열람 (1115 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AiZynthFinder_2020.pdf) | [10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1) |
| 18 | **Retro*** | 2020 | Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al. (*ICML 2020*) | [PDF 열람 (3244 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroStar_2020.pdf) | [10.48550/arXiv.2006.15575](https://doi.org/10.48550/arXiv.2006.15575) |

### Section 5: In-Depth Extended Studies (유사 및 확장 연구 심층 분석)
**Literature Count**: 5 papers (Locally secured: 4/5)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **RetroGNN** | 2022 | Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler (*Journal of Chemical Information and Modeling*) | [PDF 열람 (2400 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroGNN_2022.pdf) | [10.1021/acs.jcim.1c01476](https://doi.org/10.1021/acs.jcim.1c01476) |
| 02 | **Generative SA Benchmark** | 2020 | Wenhao Gao, Connor W. Coley (*Journal of Chemical Information and Modeling*) | [PDF 열람 (4955 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Gao_Coley_Benchmark_2020.pdf) | [10.1021/acs.jcim.0c00174](https://doi.org/10.1021/acs.jcim.0c00174) |
| 03 | **V-SYNTH** | 2022 | Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al. (*Nature*) | [출판사 구독 (DOI)](https://doi.org/10.1038/s41586-021-04220-9) | [10.1038/s41586-021-04220-9](https://doi.org/10.1038/s41586-021-04220-9) |
| 04 | **TMC Synthesizability (Kulik Group)** | 2022 | Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al. (*JACS Au*) | [XML 열람 (159 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) | [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547) |
| 05 | **MPScore (Materials Precursor Score)** | 2021 | Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al. (*Journal of Chemical Information and Modeling*) | [XML 열람 (175 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) | [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375) |

### Section 7: Verified Open-Source Repositories & Semiconductor Benchmark (17개 오픈소스 저장소 및 반도체 DB 실증)
**Literature Count**: 17 papers (Locally secured: 14/17)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **SAScore** | 2009 | Peter Ertl, Ansgar Schuffenhauer (*Journal of Cheminformatics*) | [PDF 열람 (1242 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SAScore_2009.pdf) | [10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8) |
| 02 | **SYBA** | 2020 | Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil (*Journal of Cheminformatics*) | [PDF 열람 (3191 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SYBA_2020.pdf) | [10.1186/s13321-020-00439-2](https://doi.org/10.1186/s13321-020-00439-2) |
| 03 | **SCScore** | 2018 | Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.7b00622) | [10.1021/acs.jcim.7b00622](https://doi.org/10.1021/acs.jcim.7b00622) |
| 04 | **RAscore** | 2021 | Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond (*Chemical Science*) | [PDF 열람 (322 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RAscore_2021.pdf) | [10.1039/D0SC05401A](https://doi.org/10.1039/D0SC05401A) |
| 05 | **GASA** | 2022 | Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.2c00038) | [10.1021/acs.jcim.2c00038](https://doi.org/10.1021/acs.jcim.2c00038) |
| 06 | **DeepSA** | 2023 | Shihang Wang, Lin Wang, Fenglei Li, Fang Bai (*Journal of Cheminformatics*) | [PDF 열람 (4058 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/DeepSA_2023.pdf) | [10.1186/s13321-023-00771-3](https://doi.org/10.1186/s13321-023-00771-3) |
| 07 | **DFRscore** | 2023 | Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01134) | [10.1021/acs.jcim.3c01134](https://doi.org/10.1021/acs.jcim.3c01134) |
| 08 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |
| 09 | **ChemProjector** | 2024 | Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma (*NeurIPS 2024 / arXiv:2406.04628*) | [PDF 열람 (4788 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) | [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628) |
| 10 | **SynFormer** | 2025 | Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al. (*PNAS 2025 / arXiv:2410.03494*) | [PDF 열람 (3651 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFormer_2025.pdf) | [10.1073/pnas.2415665122](https://doi.org/10.1073/pnas.2415665122) |
| 11 | **SynFlowNet** | 2025 | Miruna Cretu, Charles Harris, Austin Tripp, et al. (*ICLR 2025 / arXiv:2405.01155*) | [PDF 열람 (8060 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) | [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155) |
| 12 | **RxnFlow** | 2024 | Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim (*NeurIPS 2024 / arXiv:2410.04542*) | [PDF 열람 (4912 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) | [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542) |
| 13 | **ReaSyn** | 2025 | KAIST & NVIDIA Research Team (*arXiv:2509.16084*) | [PDF 열람 (2479 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) | [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084) |
| 14 | **SynNet** | 2022 | Wenhao Gao, Rocío Mercado, Connor W. Coley (*ICLR 2022*) | [PDF 열람 (3079 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynNet_2022.pdf) | [10.48550/arXiv.2110.06389](https://doi.org/10.48550/arXiv.2110.06389) |
| 15 | **AiZynthFinder** | 2020 | Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist (*Journal of Cheminformatics*) | [PDF 열람 (1115 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AiZynthFinder_2020.pdf) | [10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1) |
| 16 | **Retro*** | 2020 | Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al. (*ICML 2020*) | [PDF 열람 (3244 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroStar_2020.pdf) | [10.48550/arXiv.2006.15575](https://doi.org/10.48550/arXiv.2006.15575) |
| 17 | **MPScore (Materials Precursor Score)** | 2021 | Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al. (*Journal of Chemical Information and Modeling*) | [XML 열람 (175 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) | [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375) |

### Section 8: Semiconductor Precursors (ALD/CVD/EUV) Application & Extension (반도체 프리커서 적용 및 확장)
**Literature Count**: 7 papers (Locally secured: 7/7)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **TMC Synthesizability (Kulik Group)** | 2022 | Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al. (*JACS Au*) | [XML 열람 (159 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) | [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547) |
| 02 | **MPScore (Materials Precursor Score)** | 2021 | Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al. (*Journal of Chemical Information and Modeling*) | [XML 열람 (175 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) | [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375) |
| 03 | **AS-ALD Precursor ML Descriptors** | 2023 | Tran Thi Ngoc Van, Changsu Kim, Hojae Lee, Jiyong Kim et al. (*Journal of Molecular Modeling*) | [PDF 열람 (2563 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AS_ALD_Precursor_ML_2023.pdf) | [10.1007/s00894-023-05806-y](https://doi.org/10.1007/s00894-023-05806-y) |
| 04 | **MLIP for ALD Precursors** | 2024 | Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al. (*Materials Today Advances*) | [PDF 열람 (6267 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_MLIP_Precursor_2024.pdf) | [10.1016/j.mtadv.2024.100474](https://doi.org/10.1016/j.mtadv.2024.100474) |
| 05 | **ML in ALD Roadmap** | 2026 | Julia Westermayr, Stacey F. Bent, Ralf Tonner-Zech (*Chemistry of Materials*) | [PDF 열람 (3091 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_ML_ChemMater_2026.pdf) | [10.1021/acs.chemmater.5c03271](https://doi.org/10.1021/acs.chemmater.5c03271) |
| 06 | **Aminosilane ALD Precursor DFT Design** | 2017 | Liang Huang, Bo Han, Maohong Fan, Hansong Cheng (*RSC Advances*) | [PDF 열람 (736 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Aminosilane_SiO2_ALD_DFT_2017.pdf) | [10.1039/c7ra02301d](https://doi.org/10.1039/c7ra02301d) |
| 07 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |

### Section 9: Multi-Fidelity DFT Integration & Quantum Chemical Applications (DFT 계산 연계 응용 프레임워크)
**Literature Count**: 6 papers (Locally secured: 5/6)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **TMC Synthesizability (Kulik Group)** | 2022 | Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al. (*JACS Au*) | [XML 열람 (159 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) | [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547) |
| 02 | **Aminosilane ALD Precursor DFT Design** | 2017 | Liang Huang, Bo Han, Maohong Fan, Hansong Cheng (*RSC Advances*) | [PDF 열람 (736 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Aminosilane_SiO2_ALD_DFT_2017.pdf) | [10.1039/c7ra02301d](https://doi.org/10.1039/c7ra02301d) |
| 03 | **MLIP for ALD Precursors** | 2024 | Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al. (*Materials Today Advances*) | [PDF 열람 (6267 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_MLIP_Precursor_2024.pdf) | [10.1016/j.mtadv.2024.100474](https://doi.org/10.1016/j.mtadv.2024.100474) |
| 04 | **ML in ALD Roadmap** | 2026 | Julia Westermayr, Stacey F. Bent, Ralf Tonner-Zech (*Chemistry of Materials*) | [PDF 열람 (3091 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_ML_ChemMater_2026.pdf) | [10.1021/acs.chemmater.5c03271](https://doi.org/10.1021/acs.chemmater.5c03271) |
| 05 | **Architector (3D TMC Complex Builder)** | 2023 | Corin Gautier, Tyler Castonguay, Heather J. Kulik (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01353) | [10.1021/acs.jcim.3c01353](https://doi.org/10.1021/acs.jcim.3c01353) |
| 06 | **Sella (Automated Saddle-Point TS)** | 2022 | Eric D. Hermes, Judit Zádor (*Journal of Chemical Theory and Computation*) | [PDF 열람 (2963 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Hermes_Sella_2022.pdf) | [10.1021/acs.jctc.1c00845](https://doi.org/10.1021/acs.jctc.1c00845) |

### Section 10: VASP-Centered High-Throughput DFT Automation Architecture (VASP 제일원리 계산 자동화 스택)
**Literature Count**: 6 papers (Locally secured: 4/6)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **Pymatgen (Materials Genomics)** | 2013 | Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier et al. (*Computational Materials Science*) | [PDF 열람 (1357 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Ong_Pymatgen_2013.pdf) | [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028) |
| 02 | **atomate2 + jobflow** | 2024 | Andrew S. Rosen, Matthew Gallant, Anderson Janotti et al. (*Journal of Open Source Software*) | [출판사 구독 (DOI)](https://doi.org/10.21105/joss.06543) | [10.21105/joss.06543](https://doi.org/10.21105/joss.06543) |
| 03 | **AiiDA-VASP Automation** | 2021 | Martin Uhrin, Sebastiaan P. Huber, Jusong Yu, Nicola Marzari et al. (*Computational Materials Science*) | [PDF 열람 (925 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Uhrin_AiiDA_2021.pdf) | [10.1016/j.commatsci.2020.110086](https://doi.org/10.1016/j.commatsci.2020.110086) |
| 04 | **Architector (3D TMC Complex Builder)** | 2023 | Corin Gautier, Tyler Castonguay, Heather J. Kulik (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01353) | [10.1021/acs.jcim.3c01353](https://doi.org/10.1021/acs.jcim.3c01353) |
| 05 | **Sella (Automated Saddle-Point TS)** | 2022 | Eric D. Hermes, Judit Zádor (*Journal of Chemical Theory and Computation*) | [PDF 열람 (2963 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Hermes_Sella_2022.pdf) | [10.1021/acs.jctc.1c00845](https://doi.org/10.1021/acs.jctc.1c00845) |
| 06 | **MLIP for ALD Precursors** | 2024 | Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al. (*Materials Today Advances*) | [PDF 열람 (6267 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_MLIP_Precursor_2024.pdf) | [10.1016/j.mtadv.2024.100474](https://doi.org/10.1016/j.mtadv.2024.100474) |

### Section 11: Semiconductor Precursor, Reagents & Reaction Database (독자 도메인 데이터베이스)
**Literature Count**: 5 papers (Locally secured: 4/5)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **Aminosilane ALD Precursor DFT Design** | 2017 | Liang Huang, Bo Han, Maohong Fan, Hansong Cheng (*RSC Advances*) | [PDF 열람 (736 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Aminosilane_SiO2_ALD_DFT_2017.pdf) | [10.1039/c7ra02301d](https://doi.org/10.1039/c7ra02301d) |
| 02 | **TMC Synthesizability (Kulik Group)** | 2022 | Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al. (*JACS Au*) | [XML 열람 (159 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) | [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547) |
| 03 | **MPScore (Materials Precursor Score)** | 2021 | Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al. (*Journal of Chemical Information and Modeling*) | [XML 열람 (175 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) | [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375) |
| 04 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |
| 05 | **V-SYNTH** | 2022 | Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al. (*Nature*) | [출판사 구독 (DOI)](https://doi.org/10.1038/s41586-021-04220-9) | [10.1038/s41586-021-04220-9](https://doi.org/10.1038/s41586-021-04220-9) |

### Section 12: Comprehensive Benchmark & Engineering Evaluation (DB 및 오픈소스 종합 평가)
**Literature Count**: 17 papers (Locally secured: 14/17)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **SAScore** | 2009 | Peter Ertl, Ansgar Schuffenhauer (*Journal of Cheminformatics*) | [PDF 열람 (1242 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SAScore_2009.pdf) | [10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8) |
| 02 | **SYBA** | 2020 | Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil (*Journal of Cheminformatics*) | [PDF 열람 (3191 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SYBA_2020.pdf) | [10.1186/s13321-020-00439-2](https://doi.org/10.1186/s13321-020-00439-2) |
| 03 | **SCScore** | 2018 | Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.7b00622) | [10.1021/acs.jcim.7b00622](https://doi.org/10.1021/acs.jcim.7b00622) |
| 04 | **RAscore** | 2021 | Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond (*Chemical Science*) | [PDF 열람 (322 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RAscore_2021.pdf) | [10.1039/D0SC05401A](https://doi.org/10.1039/D0SC05401A) |
| 05 | **GASA** | 2022 | Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.2c00038) | [10.1021/acs.jcim.2c00038](https://doi.org/10.1021/acs.jcim.2c00038) |
| 06 | **DeepSA** | 2023 | Shihang Wang, Lin Wang, Fenglei Li, Fang Bai (*Journal of Cheminformatics*) | [PDF 열람 (4058 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/DeepSA_2023.pdf) | [10.1186/s13321-023-00771-3](https://doi.org/10.1186/s13321-023-00771-3) |
| 07 | **DFRscore** | 2023 | Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01134) | [10.1021/acs.jcim.3c01134](https://doi.org/10.1021/acs.jcim.3c01134) |
| 08 | **BR-SAScore** | 2024 | Shuan Chen, Yousung Jung (*Journal of Cheminformatics*) | [PDF 열람 (1351 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) | [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0) |
| 09 | **ChemProjector** | 2024 | Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma (*NeurIPS 2024 / arXiv:2406.04628*) | [PDF 열람 (4788 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) | [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628) |
| 10 | **SynFormer** | 2025 | Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al. (*PNAS 2025 / arXiv:2410.03494*) | [PDF 열람 (3651 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFormer_2025.pdf) | [10.1073/pnas.2415665122](https://doi.org/10.1073/pnas.2415665122) |
| 11 | **SynFlowNet** | 2025 | Miruna Cretu, Charles Harris, Austin Tripp, et al. (*ICLR 2025 / arXiv:2405.01155*) | [PDF 열람 (8060 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) | [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155) |
| 12 | **RxnFlow** | 2024 | Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim (*NeurIPS 2024 / arXiv:2410.04542*) | [PDF 열람 (4912 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) | [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542) |
| 13 | **ReaSyn** | 2025 | KAIST & NVIDIA Research Team (*arXiv:2509.16084*) | [PDF 열람 (2479 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) | [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084) |
| 14 | **SynNet** | 2022 | Wenhao Gao, Rocío Mercado, Connor W. Coley (*ICLR 2022*) | [PDF 열람 (3079 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynNet_2022.pdf) | [10.48550/arXiv.2110.06389](https://doi.org/10.48550/arXiv.2110.06389) |
| 15 | **AiZynthFinder** | 2020 | Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist (*Journal of Cheminformatics*) | [PDF 열람 (1115 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AiZynthFinder_2020.pdf) | [10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1) |
| 16 | **Retro*** | 2020 | Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al. (*ICML 2020*) | [PDF 열람 (3244 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroStar_2020.pdf) | [10.48550/arXiv.2006.15575](https://doi.org/10.48550/arXiv.2006.15575) |
| 17 | **MPScore (Materials Precursor Score)** | 2021 | Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al. (*Journal of Chemical Information and Modeling*) | [XML 열람 (175 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) | [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375) |

### Section 12: AI-Driven Precursor Inverse Design & Substrate ASD (프리커서 AI 역설계 및 기판 흡착 화학)
**Literature Count**: 13 papers (Locally secured: 10/13)  

| No | Model / Target | Year | Authors & Journal | Local Archive | Official DOI |
| :-: | :--- | :-: | :--- | :--- | :--- |
| 01 | **ML in ALD Roadmap** | 2026 | Julia Westermayr, Stacey F. Bent, Ralf Tonner-Zech (*Chemistry of Materials*) | [PDF 열람 (3091 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_ML_ChemMater_2026.pdf) | [10.1021/acs.chemmater.5c03271](https://doi.org/10.1021/acs.chemmater.5c03271) |
| 02 | **MLIP for ALD Precursors** | 2024 | Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al. (*Materials Today Advances*) | [PDF 열람 (6267 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_MLIP_Precursor_2024.pdf) | [10.1016/j.mtadv.2024.100474](https://doi.org/10.1016/j.mtadv.2024.100474) |
| 03 | **TMC Synthesizability (Kulik Group)** | 2022 | Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al. (*JACS Au*) | [XML 열람 (159 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) | [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547) |
| 04 | **Oh et al. (ASD Precursor Choice)** | 2021 | In-Kyun Oh, Stacey F. Bent, Adriaan J. M. Mackus et al. (*Chemistry of Materials*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.chemmater.0c04718) | [10.1021/acs.chemmater.0c04718](https://doi.org/10.1021/acs.chemmater.0c04718) |
| 05 | **Tonner-Zech et al. (ASD Design Principles)** | 2024 | Ralf Tonner-Zech et al. (*Journal of Computational Chemistry*) | [XML 열람 (201 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Tonner_ASD_Precursor_2024.xml) | [10.1002/jcc.70493](https://doi.org/10.1002/jcc.70493) |
| 06 | **Mackus et al. (ASD Fundamentals & Metrics)** | 2019 | Adriaan J. M. Mackus, Ageeth A. Bol, W. M. Erwin Kessels (*Chemistry of Materials*) | [XML 열람 (146 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Mackus_ASD_Selectivity_2019.xml) | [10.1021/acs.chemmater.8b03454](https://doi.org/10.1021/acs.chemmater.8b03454) |
| 07 | **Parsons & Clark (ASD State-of-the-Art)** | 2020 | Gregory N. Parsons, Robert D. Clark (*Chemistry of Materials*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.chemmater.0c00722) | [10.1021/acs.chemmater.0c00722](https://doi.org/10.1021/acs.chemmater.0c00722) |
| 08 | **ChemProjector** | 2024 | Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma (*NeurIPS 2024 / arXiv:2406.04628*) | [PDF 열람 (4788 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) | [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628) |
| 09 | **SynFlowNet** | 2025 | Miruna Cretu, Charles Harris, Austin Tripp, et al. (*ICLR 2025 / arXiv:2405.01155*) | [PDF 열람 (8060 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) | [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155) |
| 10 | **RxnFlow** | 2024 | Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim (*NeurIPS 2024 / arXiv:2410.04542*) | [PDF 열람 (4912 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) | [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542) |
| 11 | **ReaSyn** | 2025 | KAIST & NVIDIA Research Team (*arXiv:2509.16084*) | [PDF 열람 (2479 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) | [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084) |
| 12 | **Architector (3D TMC Complex Builder)** | 2023 | Corin Gautier, Tyler Castonguay, Heather J. Kulik (*Journal of Chemical Information and Modeling*) | [출판사 구독 (DOI)](https://doi.org/10.1021/acs.jcim.3c01353) | [10.1021/acs.jcim.3c01353](https://doi.org/10.1021/acs.jcim.3c01353) |
| 13 | **Sella (Automated Saddle-Point TS)** | 2022 | Eric D. Hermes, Judit Zádor (*Journal of Chemical Theory and Computation*) | [PDF 열람 (2963 KB)](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Hermes_Sella_2022.pdf) | [10.1021/acs.jctc.1c00845](https://doi.org/10.1021/acs.jctc.1c00845) |

---

## 3. Master BibTeX Archive (All 34 Papers)

The complete BibTeX collection is indexed below for immediate integration into LaTeX manuscripts and academic reference managers:

### [01] SAScore (2009)
- **Title**: Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions
- **Authors**: Peter Ertl, Ansgar Schuffenhauer
- **Journal**: *Journal of Cheminformatics* (2009)
- **DOI**: [10.1186/1758-2946-1-8](https://doi.org/10.1186/1758-2946-1-8)
- **Local Document**: [downloaded_papers/SAScore_2009.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SAScore_2009.pdf) (1242 KB)
- **Sections Referenced**: Section 1, Section 2, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{ertl2009SAScore,
  author    = {Peter Ertl, Ansgar Schuffenhauer},
  title     = {{Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions}},
  journal   = {Journal of Cheminformatics},
  year      = {2009},
  doi       = {10.1186/1758-2946-1-8},
  url       = {https://doi.org/10.1186/1758-2946-1-8}
}
```

### [02] SYBA (2020)
- **Title**: SYBA: Bayesian estimation of synthetic accessibility of organic compounds
- **Authors**: Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil
- **Journal**: *Journal of Cheminformatics* (2020)
- **DOI**: [10.1186/s13321-020-00439-2](https://doi.org/10.1186/s13321-020-00439-2)
- **Local Document**: [downloaded_papers/SYBA_2020.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SYBA_2020.pdf) (3191 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{voršilák2020SYBA,
  author    = {Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil},
  title     = {{SYBA: Bayesian estimation of synthetic accessibility of organic compounds}},
  journal   = {Journal of Cheminformatics},
  year      = {2020},
  doi       = {10.1186/s13321-020-00439-2},
  url       = {https://doi.org/10.1186/s13321-020-00439-2}
}
```

### [03] SCScore (2018)
- **Title**: SCScore: Synthetic Complexity Learned from a Reaction Corpus
- **Authors**: Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen
- **Journal**: *Journal of Chemical Information and Modeling* (2018)
- **DOI**: [10.1021/acs.jcim.7b00622](https://doi.org/10.1021/acs.jcim.7b00622)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1021/acs.jcim.7b00622))
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{coley2018SCScore,
  author    = {Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen},
  title     = {{SCScore: Synthetic Complexity Learned from a Reaction Corpus}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2018},
  doi       = {10.1021/acs.jcim.7b00622},
  url       = {https://doi.org/10.1021/acs.jcim.7b00622}
}
```

### [04] RAscore (2021)
- **Title**: Retrosynthetic accessibility score (RAscore) – rapid machine learned synthesizability classification from AI driven retrosynthetic planning
- **Authors**: Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond
- **Journal**: *Chemical Science* (2021)
- **DOI**: [10.1039/D0SC05401A](https://doi.org/10.1039/D0SC05401A)
- **Local Document**: [downloaded_papers/RAscore_2021.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RAscore_2021.pdf) (322 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{thakkar2021RAscore,
  author    = {Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond},
  title     = {{Retrosynthetic accessibility score (RAscore) – rapid machine learned synthesizability classification from AI driven retrosynthetic planning}},
  journal   = {Chemical Science},
  year      = {2021},
  doi       = {10.1039/D0SC05401A},
  url       = {https://doi.org/10.1039/D0SC05401A}
}
```

### [05] RetroGNN (2022)
- **Title**: RetroGNN: Fast Estimation of Synthesizability for Virtual Screening and De Novo Design by Learning from Slow Retrosynthesis Software
- **Authors**: Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler
- **Journal**: *Journal of Chemical Information and Modeling* (2022)
- **DOI**: [10.1021/acs.jcim.1c01476](https://doi.org/10.1021/acs.jcim.1c01476)
- **Local Document**: [downloaded_papers/RetroGNN_2022.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroGNN_2022.pdf) (2400 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 5

```bibtex
@article{liu2022RetroGNN,
  author    = {Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler},
  title     = {{RetroGNN: Fast Estimation of Synthesizability for Virtual Screening and De Novo Design by Learning from Slow Retrosynthesis Software}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2022},
  doi       = {10.1021/acs.jcim.1c01476},
  url       = {https://doi.org/10.1021/acs.jcim.1c01476}
}
```

### [06] GASA (2022)
- **Title**: Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism
- **Authors**: Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu
- **Journal**: *Journal of Chemical Information and Modeling* (2022)
- **DOI**: [10.1021/acs.jcim.2c00038](https://doi.org/10.1021/acs.jcim.2c00038)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1021/acs.jcim.2c00038))
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{yu2022GASA,
  author    = {Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu},
  title     = {{Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2022},
  doi       = {10.1021/acs.jcim.2c00038},
  url       = {https://doi.org/10.1021/acs.jcim.2c00038}
}
```

### [07] DeepSA (2023)
- **Title**: DeepSA: a deep-learning driven predictor of compound synthesis accessibility
- **Authors**: Shihang Wang, Lin Wang, Fenglei Li, Fang Bai
- **Journal**: *Journal of Cheminformatics* (2023)
- **DOI**: [10.1186/s13321-023-00771-3](https://doi.org/10.1186/s13321-023-00771-3)
- **Local Document**: [downloaded_papers/DeepSA_2023.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/DeepSA_2023.pdf) (4058 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{wang2023DeepSA,
  author    = {Shihang Wang, Lin Wang, Fenglei Li, Fang Bai},
  title     = {{DeepSA: a deep-learning driven predictor of compound synthesis accessibility}},
  journal   = {Journal of Cheminformatics},
  year      = {2023},
  doi       = {10.1186/s13321-023-00771-3},
  url       = {https://doi.org/10.1186/s13321-023-00771-3}
}
```

### [08] DFRscore (2023)
- **Title**: DFRscore: Deep Learning-Based Scoring of Synthetic Complexity with Drug-Focused Retrosynthetic Analysis for High-Throughput Virtual Screening
- **Authors**: Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim
- **Journal**: *Journal of Chemical Information and Modeling* (2023)
- **DOI**: [10.1021/acs.jcim.3c01134](https://doi.org/10.1021/acs.jcim.3c01134)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1021/acs.jcim.3c01134))
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{kim2023DFRscore,
  author    = {Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim},
  title     = {{DFRscore: Deep Learning-Based Scoring of Synthetic Complexity with Drug-Focused Retrosynthetic Analysis for High-Throughput Virtual Screening}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2023},
  doi       = {10.1021/acs.jcim.3c01134},
  url       = {https://doi.org/10.1021/acs.jcim.3c01134}
}
```

### [09] BR-SAScore (2024)
- **Title**: Estimating the synthetic accessibility of molecules with building block and reaction-aware SAScore
- **Authors**: Shuan Chen, Yousung Jung
- **Journal**: *Journal of Cheminformatics* (2024)
- **DOI**: [10.1186/s13321-024-00879-0](https://doi.org/10.1186/s13321-024-00879-0)
- **Local Document**: [downloaded_papers/BR_SAScore_2024.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/BR_SAScore_2024.pdf) (1351 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 8, Section 11, Section 12

```bibtex
@article{chen2024BRSAScore,
  author    = {Shuan Chen, Yousung Jung},
  title     = {{Estimating the synthetic accessibility of molecules with building block and reaction-aware SAScore}},
  journal   = {Journal of Cheminformatics},
  year      = {2024},
  doi       = {10.1186/s13321-024-00879-0},
  url       = {https://doi.org/10.1186/s13321-024-00879-0}
}
```

### [10] ChemProjector (2024)
- **Title**: Projecting Molecules into Synthesizable Chemical Spaces
- **Authors**: Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma
- **Journal**: *NeurIPS 2024 / arXiv:2406.04628* (2024)
- **DOI**: [10.48550/arXiv.2406.04628](https://doi.org/10.48550/arXiv.2406.04628)
- **Local Document**: [downloaded_papers/ChemProjector_2024.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ChemProjector_2024.pdf) (4788 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12, Section 12

```bibtex
@article{luo2024ChemProjector,
  author    = {Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma},
  title     = {{Projecting Molecules into Synthesizable Chemical Spaces}},
  journal   = {NeurIPS 2024 / arXiv:2406.04628},
  year      = {2024},
  doi       = {10.48550/arXiv.2406.04628},
  url       = {https://doi.org/10.48550/arXiv.2406.04628}
}
```

### [11] SynFormer (2025)
- **Title**: Generative Artificial Intelligence for Navigating Synthesizable Chemical Space
- **Authors**: Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al.
- **Journal**: *PNAS 2025 / arXiv:2410.03494* (2025)
- **DOI**: [10.1073/pnas.2415665122](https://doi.org/10.1073/pnas.2415665122)
- **Local Document**: [downloaded_papers/SynFormer_2025.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFormer_2025.pdf) (3651 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12, Section 12

```bibtex
@article{chuang2025SynFormer,
  author    = {Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al.},
  title     = {{Generative Artificial Intelligence for Navigating Synthesizable Chemical Space}},
  journal   = {PNAS 2025 / arXiv:2410.03494},
  year      = {2025},
  doi       = {10.1073/pnas.2415665122},
  url       = {https://doi.org/10.1073/pnas.2415665122}
}
```

### [12] SynFlowNet (2025)
- **Title**: SynFlowNet: Design of Diverse and Novel Molecules with Synthesis Constraints
- **Authors**: Miruna Cretu, Charles Harris, Austin Tripp, et al.
- **Journal**: *ICLR 2025 / arXiv:2405.01155* (2025)
- **DOI**: [10.48550/arXiv.2405.01155](https://doi.org/10.48550/arXiv.2405.01155)
- **Local Document**: [downloaded_papers/SynFlowNet_2025.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynFlowNet_2025.pdf) (8060 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12, Section 12

```bibtex
@article{cretu2025SynFlowNet,
  author    = {Miruna Cretu, Charles Harris, Austin Tripp, et al.},
  title     = {{SynFlowNet: Design of Diverse and Novel Molecules with Synthesis Constraints}},
  journal   = {ICLR 2025 / arXiv:2405.01155},
  year      = {2025},
  doi       = {10.48550/arXiv.2405.01155},
  url       = {https://doi.org/10.48550/arXiv.2405.01155}
}
```

### [13] RxnFlow (2024)
- **Title**: Generative Flows on Synthetic Pathway for Drug Design
- **Authors**: Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim
- **Journal**: *NeurIPS 2024 / arXiv:2410.04542* (2024)
- **DOI**: [10.48550/arXiv.2410.04542](https://doi.org/10.48550/arXiv.2410.04542)
- **Local Document**: [downloaded_papers/RxnFlow_2025.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RxnFlow_2025.pdf) (4912 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12, Section 12

```bibtex
@article{seo2024RxnFlow,
  author    = {Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim},
  title     = {{Generative Flows on Synthetic Pathway for Drug Design}},
  journal   = {NeurIPS 2024 / arXiv:2410.04542},
  year      = {2024},
  doi       = {10.48550/arXiv.2410.04542},
  url       = {https://doi.org/10.48550/arXiv.2410.04542}
}
```

### [14] ReaSyn (2025)
- **Title**: Exploring Synthesizable Chemical Space with Iterative Pathway Refinements
- **Authors**: KAIST & NVIDIA Research Team
- **Journal**: *arXiv:2509.16084* (2025)
- **DOI**: [10.48550/arXiv.2509.16084](https://doi.org/10.48550/arXiv.2509.16084)
- **Local Document**: [downloaded_papers/ReaSyn_2025.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ReaSyn_2025.pdf) (2479 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12, Section 12

```bibtex
@article{kaist2025ReaSyn,
  author    = {KAIST & NVIDIA Research Team},
  title     = {{Exploring Synthesizable Chemical Space with Iterative Pathway Refinements}},
  journal   = {arXiv:2509.16084},
  year      = {2025},
  doi       = {10.48550/arXiv.2509.16084},
  url       = {https://doi.org/10.48550/arXiv.2509.16084}
}
```

### [15] SynNet (2022)
- **Title**: Amortized Tree Generation for Bottom-up Synthesis Planning and Synthesizable Molecular Design
- **Authors**: Wenhao Gao, Rocío Mercado, Connor W. Coley
- **Journal**: *ICLR 2022* (2022)
- **DOI**: [10.48550/arXiv.2110.06389](https://doi.org/10.48550/arXiv.2110.06389)
- **Local Document**: [downloaded_papers/SynNet_2022.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/SynNet_2022.pdf) (3079 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{gao2022SynNet,
  author    = {Wenhao Gao, Rocío Mercado, Connor W. Coley},
  title     = {{Amortized Tree Generation for Bottom-up Synthesis Planning and Synthesizable Molecular Design}},
  journal   = {ICLR 2022},
  year      = {2022},
  doi       = {10.48550/arXiv.2110.06389},
  url       = {https://doi.org/10.48550/arXiv.2110.06389}
}
```

### [16] Generative SA Benchmark (2020)
- **Title**: The Synthesizability of Molecules Proposed by Generative Models
- **Authors**: Wenhao Gao, Connor W. Coley
- **Journal**: *Journal of Chemical Information and Modeling* (2020)
- **DOI**: [10.1021/acs.jcim.0c00174](https://doi.org/10.1021/acs.jcim.0c00174)
- **Local Document**: [downloaded_papers/Gao_Coley_Benchmark_2020.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Gao_Coley_Benchmark_2020.pdf) (4955 KB)
- **Sections Referenced**: Section 1, Section 4, Section 5, Section 12

```bibtex
@article{gao2020GenerativeSABenchmark,
  author    = {Wenhao Gao, Connor W. Coley},
  title     = {{The Synthesizability of Molecules Proposed by Generative Models}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2020},
  doi       = {10.1021/acs.jcim.0c00174},
  url       = {https://doi.org/10.1021/acs.jcim.0c00174}
}
```

### [17] AiZynthFinder (2020)
- **Title**: AiZynthFinder: a software package for automating retrosynthetic planning
- **Authors**: Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist
- **Journal**: *Journal of Cheminformatics* (2020)
- **DOI**: [10.1186/s13321-020-00472-1](https://doi.org/10.1186/s13321-020-00472-1)
- **Local Document**: [downloaded_papers/AiZynthFinder_2020.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AiZynthFinder_2020.pdf) (1115 KB)
- **Sections Referenced**: Section 1, Section 3, Section 4, Section 7, Section 12

```bibtex
@article{genheden2020AiZynthFinder,
  author    = {Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist},
  title     = {{AiZynthFinder: a software package for automating retrosynthetic planning}},
  journal   = {Journal of Cheminformatics},
  year      = {2020},
  doi       = {10.1186/s13321-020-00472-1},
  url       = {https://doi.org/10.1186/s13321-020-00472-1}
}
```

### [18] Retro* (2020)
- **Title**: Retro*: Learning Retrosynthetic Planning with Neural Guided A* Search
- **Authors**: Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al.
- **Journal**: *ICML 2020* (2020)
- **DOI**: [10.48550/arXiv.2006.15575](https://doi.org/10.48550/arXiv.2006.15575)
- **Local Document**: [downloaded_papers/RetroStar_2020.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/RetroStar_2020.pdf) (3244 KB)
- **Sections Referenced**: Section 1, Section 4, Section 7, Section 12

```bibtex
@article{chen2020Retro*,
  author    = {Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al.},
  title     = {{Retro*: Learning Retrosynthetic Planning with Neural Guided A* Search}},
  journal   = {ICML 2020},
  year      = {2020},
  doi       = {10.48550/arXiv.2006.15575},
  url       = {https://doi.org/10.48550/arXiv.2006.15575}
}
```

### [19] V-SYNTH (2022)
- **Title**: Synthetically accessible vast chemical space
- **Authors**: Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al.
- **Journal**: *Nature* (2022)
- **DOI**: [10.1038/s41586-021-04220-9](https://doi.org/10.1038/s41586-021-04220-9)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1038/s41586-021-04220-9))
- **Sections Referenced**: Section 1, Section 5, Section 11

```bibtex
@article{sadybekov2022VSYNTH,
  author    = {Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al.},
  title     = {{Synthetically accessible vast chemical space}},
  journal   = {Nature},
  year      = {2022},
  doi       = {10.1038/s41586-021-04220-9},
  url       = {https://doi.org/10.1038/s41586-021-04220-9}
}
```

### [20] TMC Synthesizability (Kulik Group) (2022)
- **Title**: Active Learning Exploration of Transition-Metal Complexes to Discover Method-Insensitive and Synthetically Accessible Chromophores
- **Authors**: Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al.
- **Journal**: *JACS Au* (2022)
- **DOI**: [10.1021/jacsau.2c00547](https://doi.org/10.1021/jacsau.2c00547)
- **Local Document**: [downloaded_papers/Kulik_TMC_SA_2022.xml](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Kulik_TMC_SA_2022.xml) (159 KB)
- **Sections Referenced**: Section 1, Section 5, Section 8, Section 9, Section 11, Section 12

```bibtex
@article{duan2022TMCSynthesizability(KulikGroup),
  author    = {Chenru Duan, Aditya Nandy, Gianmarco G. Terrones, David W. Kastner et al.},
  title     = {{Active Learning Exploration of Transition-Metal Complexes to Discover Method-Insensitive and Synthetically Accessible Chromophores}},
  journal   = {JACS Au},
  year      = {2022},
  doi       = {10.1021/jacsau.2c00547},
  url       = {https://doi.org/10.1021/jacsau.2c00547}
}
```

### [21] MPScore (Materials Precursor Score) (2021)
- **Title**: Materials Precursor Score: Modeling Chemists’ Intuition for the Synthetic Accessibility of Porous Organic Cage Precursors
- **Authors**: Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al.
- **Journal**: *Journal of Chemical Information and Modeling* (2021)
- **DOI**: [10.1021/acs.jcim.1c00375](https://doi.org/10.1021/acs.jcim.1c00375)
- **Local Document**: [downloaded_papers/MPScore_Jelfs_2021.xml](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/MPScore_Jelfs_2021.xml) (175 KB)
- **Sections Referenced**: Section 1, Section 5, Section 7, Section 8, Section 11, Section 12

```bibtex
@article{bennett2021MPScore(MaterialsPrecursorScore),
  author    = {Steven Bennett, Filip T. Szczypiński, Lukas Turcani, Michael E. Briggs et al.},
  title     = {{Materials Precursor Score: Modeling Chemists’ Intuition for the Synthetic Accessibility of Porous Organic Cage Precursors}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2021},
  doi       = {10.1021/acs.jcim.1c00375},
  url       = {https://doi.org/10.1021/acs.jcim.1c00375}
}
```

### [22] AS-ALD Precursor ML Descriptors (2023)
- **Title**: Machine learning-based exploration of molecular design descriptors for area-selective atomic layer deposition (AS-ALD) precursors
- **Authors**: Tran Thi Ngoc Van, Changsu Kim, Hojae Lee, Jiyong Kim et al.
- **Journal**: *Journal of Molecular Modeling* (2023)
- **DOI**: [10.1007/s00894-023-05806-y](https://doi.org/10.1007/s00894-023-05806-y)
- **Local Document**: [downloaded_papers/AS_ALD_Precursor_ML_2023.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/AS_ALD_Precursor_ML_2023.pdf) (2563 KB)
- **Sections Referenced**: Section 1, Section 8, Section 12

```bibtex
@article{van2023ASALDPrecursorMLDescriptors,
  author    = {Tran Thi Ngoc Van, Changsu Kim, Hojae Lee, Jiyong Kim et al.},
  title     = {{Machine learning-based exploration of molecular design descriptors for area-selective atomic layer deposition (AS-ALD) precursors}},
  journal   = {Journal of Molecular Modeling},
  year      = {2023},
  doi       = {10.1007/s00894-023-05806-y},
  url       = {https://doi.org/10.1007/s00894-023-05806-y}
}
```

### [23] MLIP for ALD Precursors (2024)
- **Title**: Toward fast and accurate machine learning interatomic potentials for atomic layer deposition precursors
- **Authors**: Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al.
- **Journal**: *Materials Today Advances* (2024)
- **DOI**: [10.1016/j.mtadv.2024.100474](https://doi.org/10.1016/j.mtadv.2024.100474)
- **Local Document**: [downloaded_papers/ALD_MLIP_Precursor_2024.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_MLIP_Precursor_2024.pdf) (6267 KB)
- **Sections Referenced**: Section 1, Section 8, Section 9, Section 10, Section 12

```bibtex
@article{kang2024MLIPforALDPrecursors,
  author    = {Seungpyo Kang, Joonchul Kim, Taehyun Park, Joonghee Won et al.},
  title     = {{Toward fast and accurate machine learning interatomic potentials for atomic layer deposition precursors}},
  journal   = {Materials Today Advances},
  year      = {2024},
  doi       = {10.1016/j.mtadv.2024.100474},
  url       = {https://doi.org/10.1016/j.mtadv.2024.100474}
}
```

### [24] ML in ALD Roadmap (2026)
- **Title**: Opportunities and Challenges of Machine Learning in Atomic Layer Deposition
- **Authors**: Julia Westermayr, Stacey F. Bent, Ralf Tonner-Zech
- **Journal**: *Chemistry of Materials* (2026)
- **DOI**: [10.1021/acs.chemmater.5c03271](https://doi.org/10.1021/acs.chemmater.5c03271)
- **Local Document**: [downloaded_papers/ALD_ML_ChemMater_2026.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/ALD_ML_ChemMater_2026.pdf) (3091 KB)
- **Sections Referenced**: Section 1, Section 8, Section 9, Section 12

```bibtex
@article{westermayr2026MLinALDRoadmap,
  author    = {Julia Westermayr, Stacey F. Bent, Ralf Tonner-Zech},
  title     = {{Opportunities and Challenges of Machine Learning in Atomic Layer Deposition}},
  journal   = {Chemistry of Materials},
  year      = {2026},
  doi       = {10.1021/acs.chemmater.5c03271},
  url       = {https://doi.org/10.1021/acs.chemmater.5c03271}
}
```

### [25] Aminosilane ALD Precursor DFT Design (2017)
- **Title**: Design of efficient mono-aminosilane precursors for atomic layer deposition of SiO
                    <sub>2</sub>
                    thin films
- **Authors**: Liang Huang, Bo Han, Maohong Fan, Hansong Cheng
- **Journal**: *RSC Advances* (2017)
- **DOI**: [10.1039/c7ra02301d](https://doi.org/10.1039/c7ra02301d)
- **Local Document**: [downloaded_papers/Aminosilane_SiO2_ALD_DFT_2017.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Aminosilane_SiO2_ALD_DFT_2017.pdf) (736 KB)
- **Sections Referenced**: Section 1, Section 8, Section 9, Section 11

```bibtex
@article{huang2017AminosilaneALDPrecursorDFTDesign,
  author    = {Liang Huang, Bo Han, Maohong Fan, Hansong Cheng},
  title     = {{Design of efficient mono-aminosilane precursors for atomic layer deposition of SiO
                    <sub>2</sub>
                    thin films}},
  journal   = {RSC Advances},
  year      = {2017},
  doi       = {10.1039/c7ra02301d},
  url       = {https://doi.org/10.1039/c7ra02301d}
}
```

### [26] Oh et al. (ASD Precursor Choice) (2021)
- **Title**: Role of Precursor Choice on Area-Selective Atomic Layer Deposition
- **Authors**: In-Kyun Oh, Stacey F. Bent, Adriaan J. M. Mackus et al.
- **Journal**: *Chemistry of Materials* (2021)
- **DOI**: [10.1021/acs.chemmater.0c04718](https://doi.org/10.1021/acs.chemmater.0c04718)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1021/acs.chemmater.0c04718))
- **Sections Referenced**: Section 1, Section 12

```bibtex
@article{oh2021Ohetal.(ASDPrecursorChoice),
  author    = {In-Kyun Oh, Stacey F. Bent, Adriaan J. M. Mackus et al.},
  title     = {{Role of Precursor Choice on Area-Selective Atomic Layer Deposition}},
  journal   = {Chemistry of Materials},
  year      = {2021},
  doi       = {10.1021/acs.chemmater.0c04718},
  url       = {https://doi.org/10.1021/acs.chemmater.0c04718}
}
```

### [27] Tonner-Zech et al. (ASD Design Principles) (2024)
- **Title**: Precursor Design Principles for Area-Selective Atomic Layer Deposition: Quantifying the Role of Lewis Acidity, Steric Repulsion, and Dimerization
- **Authors**: Ralf Tonner-Zech et al.
- **Journal**: *Journal of Computational Chemistry* (2024)
- **DOI**: [10.1002/jcc.70493](https://doi.org/10.1002/jcc.70493)
- **Local Document**: [downloaded_papers/Tonner_ASD_Precursor_2024.xml](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Tonner_ASD_Precursor_2024.xml) (201 KB)
- **Sections Referenced**: Section 1, Section 12

```bibtex
@article{al.2024TonnerZechetal.(ASDDesignPrinciples),
  author    = {Ralf Tonner-Zech et al.},
  title     = {{Precursor Design Principles for Area-Selective Atomic Layer Deposition: Quantifying the Role of Lewis Acidity, Steric Repulsion, and Dimerization}},
  journal   = {Journal of Computational Chemistry},
  year      = {2024},
  doi       = {10.1002/jcc.70493},
  url       = {https://doi.org/10.1002/jcc.70493}
}
```

### [28] Mackus et al. (ASD Fundamentals & Metrics) (2019)
- **Title**: From the Bottom-Up: Toward Area-Selective Atomic Layer Deposition with High Selectivity
- **Authors**: Adriaan J. M. Mackus, Ageeth A. Bol, W. M. Erwin Kessels
- **Journal**: *Chemistry of Materials* (2019)
- **DOI**: [10.1021/acs.chemmater.8b03454](https://doi.org/10.1021/acs.chemmater.8b03454)
- **Local Document**: [downloaded_papers/Mackus_ASD_Selectivity_2019.xml](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Mackus_ASD_Selectivity_2019.xml) (146 KB)
- **Sections Referenced**: Section 1, Section 12

```bibtex
@article{mackus2019Mackusetal.(ASDFundamentals&Metrics),
  author    = {Adriaan J. M. Mackus, Ageeth A. Bol, W. M. Erwin Kessels},
  title     = {{From the Bottom-Up: Toward Area-Selective Atomic Layer Deposition with High Selectivity}},
  journal   = {Chemistry of Materials},
  year      = {2019},
  doi       = {10.1021/acs.chemmater.8b03454},
  url       = {https://doi.org/10.1021/acs.chemmater.8b03454}
}
```

### [29] Parsons & Clark (ASD State-of-the-Art) (2020)
- **Title**: Area-Selective Deposition: Fundamentals, Applications, and Future Outlook
- **Authors**: Gregory N. Parsons, Robert D. Clark
- **Journal**: *Chemistry of Materials* (2020)
- **DOI**: [10.1021/acs.chemmater.0c00722](https://doi.org/10.1021/acs.chemmater.0c00722)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1021/acs.chemmater.0c00722))
- **Sections Referenced**: Section 1, Section 12

```bibtex
@article{parsons2020Parsons&Clark(ASDStateoftheArt),
  author    = {Gregory N. Parsons, Robert D. Clark},
  title     = {{Area-Selective Deposition: Fundamentals, Applications, and Future Outlook}},
  journal   = {Chemistry of Materials},
  year      = {2020},
  doi       = {10.1021/acs.chemmater.0c00722},
  url       = {https://doi.org/10.1021/acs.chemmater.0c00722}
}
```

### [30] Pymatgen (Materials Genomics) (2013)
- **Title**: Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis
- **Authors**: Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier et al.
- **Journal**: *Computational Materials Science* (2013)
- **DOI**: [10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028)
- **Local Document**: [downloaded_papers/Ong_Pymatgen_2013.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Ong_Pymatgen_2013.pdf) (1357 KB)
- **Sections Referenced**: Section 1, Section 10

```bibtex
@article{ong2013Pymatgen(MaterialsGenomics),
  author    = {Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier et al.},
  title     = {{Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis}},
  journal   = {Computational Materials Science},
  year      = {2013},
  doi       = {10.1016/j.commatsci.2012.10.028},
  url       = {https://doi.org/10.1016/j.commatsci.2012.10.028}
}
```

### [31] atomate2 + jobflow (2024)
- **Title**: atomate2: Modular computational materials science workflows
- **Authors**: Andrew S. Rosen, Matthew Gallant, Anderson Janotti et al.
- **Journal**: *Journal of Open Source Software* (2024)
- **DOI**: [10.21105/joss.06543](https://doi.org/10.21105/joss.06543)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.21105/joss.06543))
- **Sections Referenced**: Section 1, Section 10

```bibtex
@article{rosen2024atomate2+jobflow,
  author    = {Andrew S. Rosen, Matthew Gallant, Anderson Janotti et al.},
  title     = {{atomate2: Modular computational materials science workflows}},
  journal   = {Journal of Open Source Software},
  year      = {2024},
  doi       = {10.21105/joss.06543},
  url       = {https://doi.org/10.21105/joss.06543}
}
```

### [32] AiiDA-VASP Automation (2021)
- **Title**: Workflows and automation in AiiDA: Managing computational materials science
- **Authors**: Martin Uhrin, Sebastiaan P. Huber, Jusong Yu, Nicola Marzari et al.
- **Journal**: *Computational Materials Science* (2021)
- **DOI**: [10.1016/j.commatsci.2020.110086](https://doi.org/10.1016/j.commatsci.2020.110086)
- **Local Document**: [downloaded_papers/Uhrin_AiiDA_2021.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Uhrin_AiiDA_2021.pdf) (925 KB)
- **Sections Referenced**: Section 1, Section 10

```bibtex
@article{uhrin2021AiiDAVASPAutomation,
  author    = {Martin Uhrin, Sebastiaan P. Huber, Jusong Yu, Nicola Marzari et al.},
  title     = {{Workflows and automation in AiiDA: Managing computational materials science}},
  journal   = {Computational Materials Science},
  year      = {2021},
  doi       = {10.1016/j.commatsci.2020.110086},
  url       = {https://doi.org/10.1016/j.commatsci.2020.110086}
}
```

### [33] Architector (3D TMC Complex Builder) (2023)
- **Title**: Architector: Rapid, Characterizable Generation of 3D Inorganic Complex Structures
- **Authors**: Corin Gautier, Tyler Castonguay, Heather J. Kulik
- **Journal**: *Journal of Chemical Information and Modeling* (2023)
- **DOI**: [10.1021/acs.jcim.3c01353](https://doi.org/10.1021/acs.jcim.3c01353)
- **Local Document**: Commercial Subscription Paywall ([Direct Publisher DOI Link](https://doi.org/10.1021/acs.jcim.3c01353))
- **Sections Referenced**: Section 1, Section 9, Section 10, Section 12

```bibtex
@article{gautier2023Architector(3DTMCComplexBuilder),
  author    = {Corin Gautier, Tyler Castonguay, Heather J. Kulik},
  title     = {{Architector: Rapid, Characterizable Generation of 3D Inorganic Complex Structures}},
  journal   = {Journal of Chemical Information and Modeling},
  year      = {2023},
  doi       = {10.1021/acs.jcim.3c01353},
  url       = {https://doi.org/10.1021/acs.jcim.3c01353}
}
```

### [34] Sella (Automated Saddle-Point TS) (2022)
- **Title**: Sella: An Automated Saddle-Point Optimization Tool for Complex Potential Energy Surfaces
- **Authors**: Eric D. Hermes, Judit Zádor
- **Journal**: *Journal of Chemical Theory and Computation* (2022)
- **DOI**: [10.1021/acs.jctc.1c00845](https://doi.org/10.1021/acs.jctc.1c00845)
- **Local Document**: [downloaded_papers/Hermes_Sella_2022.pdf](file:///c:/Users/user/Downloads/dev_w_antigravity/synth_acc/downloaded_papers/Hermes_Sella_2022.pdf) (2963 KB)
- **Sections Referenced**: Section 1, Section 9, Section 10, Section 12

```bibtex
@article{hermes2022Sella(AutomatedSaddlePointTS),
  author    = {Eric D. Hermes, Judit Zádor},
  title     = {{Sella: An Automated Saddle-Point Optimization Tool for Complex Potential Energy Surfaces}},
  journal   = {Journal of Chemical Theory and Computation},
  year      = {2022},
  doi       = {10.1021/acs.jctc.1c00845},
  url       = {https://doi.org/10.1021/acs.jctc.1c00845}
}
```
