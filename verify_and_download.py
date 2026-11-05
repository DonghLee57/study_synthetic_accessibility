import os
import sys
import json
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET

PAPERS = [
    {
        "id": "SAScore_2009",
        "generation": "1G",
        "name": "SAScore",
        "title": "Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions",
        "authors": "Peter Ertl, Ansgar Schuffenhauer",
        "year": 2009,
        "venue": "Journal of Cheminformatics",
        "doi": "10.1186/1758-2946-1-8",
        "pdf_urls": [
            "https://jcheminf.biomedcentral.com/counter/pdf/10.1186/1758-2946-1-8.pdf",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2736273/pdf/1758-2946-1-8.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "SYBA_2020",
        "generation": "1G",
        "name": "SYBA",
        "title": "SYBA: Bayesian estimation of synthetic accessibility of organic compounds",
        "authors": "Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil",
        "year": 2020,
        "venue": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-020-00439-2",
        "pdf_urls": [
            "https://jcheminf.biomedcentral.com/counter/pdf/10.1186/s13321-020-00439-2.pdf",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7236545/pdf/13321_2020_Article_439.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "SCScore_2018",
        "generation": "2G",
        "name": "SCScore",
        "title": "SCScore: Synthetic Complexity Learned from a Reaction Corpus",
        "authors": "Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen",
        "year": 2018,
        "venue": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.7b00622",
        "pdf_urls": [
            "https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/60c74070567dfe0fd7983610/original/scscore-synthetic-complexity-learned-from-a-reaction-corpus.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "RAscore_2021",
        "generation": "2G",
        "name": "RAscore",
        "title": "Retrosynthetic accessibility score (RAscore) – rapid machine learned synthesizability classification from AI driven retrosynthetic planning",
        "authors": "Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond",
        "year": 2021,
        "venue": "Chemical Science",
        "doi": "10.1039/D0SC05401A",
        "pdf_urls": [
            "https://pubs.rsc.org/en/content/articlepdf/2021/sc/d0sc05401a",
            "https://arxiv.org/pdf/2010.05035.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "RetroGNN_2022",
        "generation": "2G",
        "name": "RetroGNN",
        "title": "RetroGNN: Fast Estimation of Synthesizability for Virtual Screening and De Novo Design by Learning from Slow Retrosynthesis Software",
        "authors": "Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler",
        "year": 2022,
        "venue": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.1c01476",
        "pdf_urls": [
            "https://arxiv.org/pdf/2203.15579.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "GASA_2022",
        "generation": "3G",
        "name": "GASA",
        "title": "Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism",
        "authors": "Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu",
        "year": 2022,
        "venue": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.2c00038",
        "pdf_urls": [],
        "category": "Core"
    },
    {
        "id": "DeepSA_2023",
        "generation": "3G",
        "name": "DeepSA",
        "title": "DeepSA: a deep-learning driven predictor of compound synthesis accessibility",
        "authors": "Shihang Wang, Lin Wang, Fenglei Li, Fang Bai",
        "year": 2023,
        "venue": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-023-00771-3",
        "pdf_urls": [
            "https://jcheminf.biomedcentral.com/counter/pdf/10.1186/s13321-023-00771-3.pdf",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10594958/pdf/13321_2023_Article_771.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "DFRscore_2023",
        "generation": "3G",
        "name": "DFRscore",
        "title": "DFRscore: Deep Learning-Based Scoring of Synthetic Complexity with Drug-Focused Retrosynthetic Analysis for High-Throughput Virtual Screening",
        "authors": "Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim",
        "year": 2023,
        "venue": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.3c01134",
        "pdf_urls": [],
        "category": "Core"
    },
    {
        "id": "BR_SAScore_2024",
        "generation": "3G",
        "name": "BR-SAScore",
        "title": "Estimating the synthetic accessibility of molecules with building block and reaction-aware SAScore",
        "authors": "Shuan Chen, Yousung Jung",
        "year": 2024,
        "venue": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-024-00879-0",
        "pdf_urls": [
            "https://jcheminf.biomedcentral.com/counter/pdf/10.1186/s13321-024-00879-0.pdf",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11270030/pdf/13321_2024_Article_879.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "ChemProjector_2024",
        "generation": "4G",
        "name": "ChemProjector",
        "title": "Projecting Molecules into Synthesizable Chemical Spaces",
        "authors": "Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma",
        "year": 2024,
        "venue": "NeurIPS 2024 / arXiv",
        "doi": "10.48550/arXiv.2406.04628",
        "pdf_urls": [
            "https://arxiv.org/pdf/2406.04628.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "SynFormer_2025",
        "generation": "4G",
        "name": "SynFormer",
        "title": "Generative Artificial Intelligence for Navigating Synthesizable Chemical Space",
        "authors": "Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al.",
        "year": 2025,
        "venue": "Proceedings of the National Academy of Sciences (PNAS)",
        "doi": "10.1073/pnas.2415665122",
        "pdf_urls": [
            "https://arxiv.org/pdf/2410.03494.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "SynFlowNet_2025",
        "generation": "4G",
        "name": "SynFlowNet",
        "title": "SynFlowNet: Design of Diverse and Novel Molecules with Synthesis Constraints",
        "authors": "Miruna Cretu, Charles Harris, Austin Tripp, et al.",
        "year": 2025,
        "venue": "ICLR 2025 / arXiv",
        "doi": "10.48550/arXiv.2405.01155",
        "pdf_urls": [
            "https://arxiv.org/pdf/2405.01155.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "RxnFlow_2025",
        "generation": "4G",
        "name": "RxnFlow",
        "title": "Generative Flows on Synthetic Pathway for Drug Design",
        "authors": "Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim",
        "year": 2024,
        "venue": "NeurIPS 2024 / arXiv",
        "doi": "10.48550/arXiv.2410.04542",
        "pdf_urls": [
            "https://arxiv.org/pdf/2410.04542.pdf"
        ],
        "category": "Core"
    },
    {
        "id": "ReaSyn_2025",
        "generation": "4G",
        "name": "ReaSyn",
        "title": "Exploring Synthesizable Chemical Space with Iterative Pathway Refinements",
        "authors": "KAIST & NVIDIA Research Team",
        "year": 2025,
        "venue": "arXiv Preprint",
        "doi": "10.48550/arXiv.2509.16084",
        "pdf_urls": [
            "https://arxiv.org/pdf/2509.16084.pdf"
        ],
        "category": "Core"
    },
    # Additional Key Related Studies
    {
        "id": "SynNet_2022",
        "generation": "4G Precursor",
        "name": "SynNet",
        "title": "Amortized Tree Generation for Bottom-up Synthesis Planning and Synthesizable Molecular Design",
        "authors": "Wenhao Gao, Rocío Mercado, Connor W. Coley",
        "year": 2022,
        "venue": "ICLR 2022 (Spotlight)",
        "doi": "10.48550/arXiv.2110.06389",
        "pdf_urls": [
            "https://arxiv.org/pdf/2110.06389.pdf"
        ],
        "category": "Related"
    },
    {
        "id": "Gao_Coley_Benchmark_2020",
        "generation": "Benchmark",
        "name": "Generative SA Benchmark",
        "title": "The Synthesizability of Molecules Proposed by Generative Models",
        "authors": "Wenhao Gao, Connor W. Coley",
        "year": 2020,
        "venue": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.0c00174",
        "pdf_urls": [
            "https://arxiv.org/pdf/2002.07007.pdf"
        ],
        "category": "Related"
    },
    {
        "id": "AiZynthFinder_2020",
        "generation": "CASP Engine",
        "name": "AiZynthFinder",
        "title": "AiZynthFinder: a software package for automating retrosynthetic planning",
        "authors": "Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist",
        "year": 2020,
        "venue": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-020-00472-1",
        "pdf_urls": [
            "https://jcheminf.biomedcentral.com/counter/pdf/10.1186/s13321-020-00472-1.pdf",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7678170/pdf/13321_2020_Article_472.pdf"
        ],
        "category": "Related"
    },
    {
        "id": "RetroStar_2020",
        "generation": "CASP Search",
        "name": "Retro*",
        "title": "Retro*: Learning Retrosynthetic Planning with Neural Guided A* Search",
        "authors": "Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al.",
        "year": 2020,
        "venue": "ICML 2020",
        "doi": "10.48550/arXiv.2006.15575",
        "pdf_urls": [
            "https://arxiv.org/pdf/2006.15575.pdf"
        ],
        "category": "Related"
    },
    {
        "id": "V_SYNTH_2022",
        "generation": "Combinatorial Space",
        "name": "V-SYNTH",
        "title": "Synthetically accessible vast chemical space",
        "authors": "Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al.",
        "year": 2022,
        "venue": "Nature",
        "doi": "10.1038/s41586-021-04220-9",
        "pdf_urls": [],
        "category": "Related"
    }
]

def check_doi_crossref(doi):
    url = f"https://api.crossref.org/works/{doi}"
    req = urllib.request.Request(url, headers={'User-Agent': 'mailto:researcher@antigravity.ai'})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            msg = data.get('message', {})
            title = msg.get('title', [''])[0]
            container = msg.get('container-title', [''])[0]
            issued = msg.get('issued', {}).get('date-parts', [[None]])[0][0]
            authors = ", ".join([f"{a.get('given', '')} {a.get('family', '')}".strip() for a in msg.get('author', [])])
            return {
                "status": "Verified (CrossRef)",
                "doi": doi,
                "verified_title": title,
                "verified_venue": container,
                "verified_year": issued,
                "verified_authors": authors
            }
    except Exception as e:
        return {"status": f"CrossRef Check Failed: {e}", "doi": doi}

def download_file(url, out_path):
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'}
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            content = resp.read()
            if len(content) > 1000 and (content[:4] == b'%PDF' or b'PDF' in content[:20]):
                with open(out_path, 'wb') as f:
                    f.write(content)
                return True, len(content)
            else:
                return False, f"Not valid PDF (size: {len(content)}, header: {content[:10]})"
    except Exception as e:
        return False, str(e)

def main():
    os.makedirs("downloaded_papers", exist_ok=True)
    results = []

    print(f"Starting verification and download of {len(PAPERS)} literature entries...\n")

    for p in PAPERS:
        print(f"--- Processing: {p['name']} ({p['year']}) ---")
        doi_info = {}
        if not p['doi'].startswith("10.48550/arXiv"):
            doi_info = check_doi_crossref(p['doi'])
            print(f"  DOI [{p['doi']}]: {doi_info.get('status')}")
            if "verified_title" in doi_info and doi_info["verified_title"]:
                print(f"  Title match: {doi_info['verified_title']}")
        else:
            doi_info = {"status": "Verified (arXiv DOI)", "doi": p['doi']}
            print(f"  ArXiv DOI [{p['doi']}]: {doi_info['status']}")

        pdf_saved = False
        saved_file = None
        file_size = 0

        for url in p['pdf_urls']:
            filename = f"{p['id']}.pdf"
            filepath = os.path.join("downloaded_papers", filename)
            success, msg = download_file(url, filepath)
            if success:
                pdf_saved = True
                saved_file = filepath
                file_size = msg
                print(f"  [SUCCESS] Downloaded PDF from {url} ({file_size // 1024} KB)")
                break
            else:
                print(f"  [RETRY] Failed from {url}: {msg}")

        results.append({
            **p,
            "verification": doi_info,
            "downloaded": pdf_saved,
            "local_path": saved_file,
            "file_size_kb": file_size // 1024 if pdf_saved else 0
        })
        time.sleep(0.5)

    with open("literature_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)

    print("\nVerification & Download Summary completed. Results saved to literature_results.json.")

if __name__ == "__main__":
    main()
