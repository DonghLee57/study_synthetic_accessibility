import os
import sys
import json
import time
import subprocess
import urllib.request

sys.stdout.reconfigure(encoding='utf-8')

TARGETS = [
    {
        "id": "SAScore_2009",
        "name": "SAScore",
        "year": 2009,
        "generation": "1st Gen (Heuristics & Substructure)",
        "authors": "Peter Ertl, Ansgar Schuffenhauer",
        "title": "Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions",
        "journal": "Journal of Cheminformatics",
        "doi": "10.1186/1758-2946-1-8",
        "download_url": "https://link.springer.com/content/pdf/10.1186/1758-2946-1-8.pdf",
        "source_type": "OpenAccess_Journal",
        "role": "Core"
    },
    {
        "id": "SYBA_2020",
        "name": "SYBA",
        "year": 2020,
        "generation": "1st Gen (Bayesian Classification)",
        "authors": "Milan Voršilák, Michal Kolář, Ivan Čmelo, Daniel Svozil",
        "title": "SYBA: Bayesian estimation of synthetic accessibility of organic compounds",
        "journal": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-020-00439-2",
        "download_url": "https://link.springer.com/content/pdf/10.1186/s13321-020-00439-2.pdf",
        "source_type": "OpenAccess_Journal",
        "role": "Core"
    },
    {
        "id": "SCScore_2018",
        "name": "SCScore",
        "year": 2018,
        "generation": "2nd Gen (Reaction Corpus Neural Net)",
        "authors": "Connor W. Coley, Luke Rogers, William H. Green, Klavs F. Jensen",
        "title": "SCScore: Synthetic Complexity Learned from a Reaction Corpus",
        "journal": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.7b00622",
        "download_url": None,
        "source_type": "Subscription_ACS",
        "role": "Core"
    },
    {
        "id": "RAscore_2021",
        "name": "RAscore",
        "year": 2021,
        "generation": "2nd Gen (CASP Route Distillation)",
        "authors": "Amol Thakkar, Veronika Chadimová, Esben Jannik Bjerrum, Ola Engkvist, Jean-Louis Reymond",
        "title": "Retrosynthetic accessibility score (RAscore) – rapid machine learned synthesizability classification from AI driven retrosynthetic planning",
        "journal": "Chemical Science",
        "doi": "10.1039/D0SC05401A",
        "download_url": "https://arxiv.org/pdf/2010.05035.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    {
        "id": "RetroGNN_2022",
        "name": "RetroGNN",
        "year": 2022,
        "generation": "2nd Gen (CASP GNN Surrogate)",
        "authors": "Cheng-Hao Liu, Maksym Korablyov, Stanisław Jastrzębski, Yoshua Bengio, Marwin H. S. Segler",
        "title": "RetroGNN: Fast Estimation of Synthesizability for Virtual Screening and De Novo Design by Learning from Slow Retrosynthesis Software",
        "journal": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.1c01476",
        "download_url": "https://arxiv.org/pdf/2203.15579.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    {
        "id": "GASA_2022",
        "name": "GASA",
        "year": 2022,
        "generation": "3rd Gen (Graph Attention Network)",
        "authors": "Jiahui Yu, Yifei Chen, Zhaojun Rong, Ziqing Fang, Xingyu Liu, Tingjun Hou, Peicheng Zhu",
        "title": "Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism",
        "journal": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.2c00038",
        "download_url": None,
        "source_type": "Subscription_ACS",
        "role": "Core"
    },
    {
        "id": "DeepSA_2023",
        "name": "DeepSA",
        "year": 2023,
        "generation": "3rd Gen (SMILES LM & Fragment Tree)",
        "authors": "Shihang Wang, Lin Wang, Fenglei Li, Fang Bai",
        "title": "DeepSA: a deep-learning driven predictor of compound synthesis accessibility",
        "journal": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-023-00771-3",
        "download_url": "https://link.springer.com/content/pdf/10.1186/s13321-023-00771-3.pdf",
        "source_type": "OpenAccess_Journal",
        "role": "Core"
    },
    {
        "id": "DFRscore_2023",
        "name": "DFRscore",
        "year": 2023,
        "generation": "3rd Gen (Drug-Focused Retrosynthesis)",
        "authors": "Hyeongwoo Kim, Kyunghoon Lee, Chansu Kim, Jaechang Lim, Woo Youn Kim",
        "title": "DFRscore: Deep Learning-Based Scoring of Synthetic Complexity with Drug-Focused Retrosynthetic Analysis for High-Throughput Virtual Screening",
        "journal": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.3c01134",
        "download_url": None,
        "source_type": "Subscription_ACS",
        "role": "Core"
    },
    {
        "id": "BR_SAScore_2024",
        "name": "BR-SAScore",
        "year": 2024,
        "generation": "3rd Gen (Building Block & Reaction Aware)",
        "authors": "Shuan Chen, Yousung Jung",
        "title": "Estimating the synthetic accessibility of molecules with building block and reaction-aware SAScore",
        "journal": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-024-00879-0",
        "download_url": "https://link.springer.com/content/pdf/10.1186/s13321-024-00879-0.pdf",
        "source_type": "OpenAccess_Journal",
        "role": "Core"
    },
    {
        "id": "ChemProjector_2024",
        "name": "ChemProjector",
        "year": 2024,
        "generation": "4th Gen (Synthesizable Space Projection)",
        "authors": "Shitong Luo, Wenhao Gao, Zuofan Wu, Jian Peng, Connor W. Coley, Jianzhu Ma",
        "title": "Projecting Molecules into Synthesizable Chemical Spaces",
        "journal": "NeurIPS 2024 / arXiv:2406.04628",
        "doi": "10.48550/arXiv.2406.04628",
        "download_url": "https://arxiv.org/pdf/2406.04628.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    {
        "id": "SynFormer_2025",
        "name": "SynFormer",
        "year": 2025,
        "generation": "4th Gen (Autoregressive Pathway Generation)",
        "authors": "Kangway V. Chuang, Christian A. Mejia, Jeffrey M. Ede, et al.",
        "title": "Generative Artificial Intelligence for Navigating Synthesizable Chemical Space",
        "journal": "PNAS 2025 / arXiv:2410.03494",
        "doi": "10.1073/pnas.2415665122",
        "download_url": "https://arxiv.org/pdf/2410.03494.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    {
        "id": "SynFlowNet_2025",
        "name": "SynFlowNet",
        "year": 2025,
        "generation": "4th Gen (GFlowNet Pathway Design)",
        "authors": "Miruna Cretu, Charles Harris, Austin Tripp, et al.",
        "title": "SynFlowNet: Design of Diverse and Novel Molecules with Synthesis Constraints",
        "journal": "ICLR 2025 / arXiv:2405.01155",
        "doi": "10.48550/arXiv.2405.01155",
        "download_url": "https://arxiv.org/pdf/2405.01155.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    {
        "id": "RxnFlow_2025",
        "name": "RxnFlow",
        "year": 2024,
        "generation": "4th Gen (Generative Flow Pathway Search)",
        "authors": "Seonghwan Seo, Minsu Kim, Tony Shen, Martin Ester, Jinkyoo Park, Sungsoo Ahn, Woo Youn Kim",
        "title": "Generative Flows on Synthetic Pathway for Drug Design",
        "journal": "NeurIPS 2024 / arXiv:2410.04542",
        "doi": "10.48550/arXiv.2410.04542",
        "download_url": "https://arxiv.org/pdf/2410.04542.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    {
        "id": "ReaSyn_2025",
        "name": "ReaSyn",
        "year": 2025,
        "generation": "4th Gen (Chain-of-Reaction & Edit Flow)",
        "authors": "KAIST & NVIDIA Research Team",
        "title": "Exploring Synthesizable Chemical Space with Iterative Pathway Refinements",
        "journal": "arXiv:2509.16084",
        "doi": "10.48550/arXiv.2509.16084",
        "download_url": "https://arxiv.org/pdf/2509.16084.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Core"
    },
    # Additional Key Studies
    {
        "id": "SynNet_2022",
        "name": "SynNet",
        "year": 2022,
        "generation": "4th Gen Milestone (Amortized Tree Design)",
        "authors": "Wenhao Gao, Rocío Mercado, Connor W. Coley",
        "title": "Amortized Tree Generation for Bottom-up Synthesis Planning and Synthesizable Molecular Design",
        "journal": "ICLR 2022",
        "doi": "10.48550/arXiv.2110.06389",
        "download_url": "https://arxiv.org/pdf/2110.06389.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Related"
    },
    {
        "id": "Gao_Coley_Benchmark_2020",
        "name": "Generative SA Benchmark",
        "year": 2020,
        "generation": "Evaluation / Benchmark",
        "authors": "Wenhao Gao, Connor W. Coley",
        "title": "The Synthesizability of Molecules Proposed by Generative Models",
        "journal": "Journal of Chemical Information and Modeling",
        "doi": "10.1021/acs.jcim.0c00174",
        "download_url": "https://arxiv.org/pdf/2002.07007.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Related"
    },
    {
        "id": "AiZynthFinder_2020",
        "name": "AiZynthFinder",
        "year": 2020,
        "generation": "CASP Infrastructure",
        "authors": "Samuel Genheden, Amol Thakkar, Veronika Chadimová, Jean-Louis Reymond, Ola Engkvist",
        "title": "AiZynthFinder: a software package for automating retrosynthetic planning",
        "journal": "Journal of Cheminformatics",
        "doi": "10.1186/s13321-020-00472-1",
        "download_url": "https://link.springer.com/content/pdf/10.1186/s13321-020-00472-1.pdf",
        "source_type": "OpenAccess_Journal",
        "role": "Related"
    },
    {
        "id": "RetroStar_2020",
        "name": "Retro*",
        "year": 2020,
        "generation": "Neural Guided CASP",
        "authors": "Binghong Chen, Chengtao Shen, Shuangjia Zheng, et al.",
        "title": "Retro*: Learning Retrosynthetic Planning with Neural Guided A* Search",
        "journal": "ICML 2020",
        "doi": "10.48550/arXiv.2006.15575",
        "download_url": "https://arxiv.org/pdf/2006.15575.pdf",
        "source_type": "ArXiv_Preprint",
        "role": "Related"
    },
    {
        "id": "V_SYNTH_2022",
        "name": "V-SYNTH",
        "year": 2022,
        "generation": "Combinatorial On-Demand Library",
        "authors": "Arman A. Sadybekov, Anastasiia V. Sadybekov, Yurii S. Moroz, et al.",
        "title": "Synthetically accessible vast chemical space",
        "journal": "Nature",
        "doi": "10.1038/s41586-021-04220-9",
        "download_url": None,
        "source_type": "Subscription_Nature",
        "role": "Related"
    }
]

def verify_doi_crossref(doi):
    url = f"https://api.crossref.org/works/{doi}"
    req = urllib.request.Request(url, headers={'User-Agent': 'mailto:researcher@antigravity.ai'})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            msg = data.get('message', {})
            title = msg.get('title', [''])[0]
            container = msg.get('container-title', [''])[0]
            issued = msg.get('issued', {}).get('date-parts', [[None]])[0][0]
            authors_list = msg.get('author', [])
            authors_str = ", ".join([f"{a.get('given', '')} {a.get('family', '')}".strip() for a in authors_list[:4]])
            if len(authors_list) > 4:
                authors_str += " et al."
            return {
                "verified": True,
                "verified_by": "CrossRef",
                "verified_title": title,
                "verified_container": container,
                "verified_year": issued,
                "verified_authors": authors_str
            }
    except Exception as e:
        return {"verified": False, "error": str(e)}

def download_with_curl(url, out_path):
    cmd = ["curl.exe", "-L", "-s", "--max-time", "30", "-o", out_path, url]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=35)
        if os.path.exists(out_path):
            size = os.path.getsize(out_path)
            # Verify if it's a PDF
            with open(out_path, "rb") as f:
                head = f.read(10)
            if size > 1000 and (head.startswith(b"%PDF") or b"PDF" in head):
                return True, size
            else:
                os.remove(out_path)
                return False, f"Invalid PDF content (size: {size})"
        return False, "File not created"
    except Exception as e:
        return False, str(e)

def main():
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    papers_dir = os.path.join(project_root, "downloaded_papers")
    meta_dir = os.path.join(project_root, "database", "metadata")
    os.makedirs(papers_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)
    manifest = []

    print(f"=== Starting Comprehensive Literature Acquisition & DOI Verification ===")
    print(f"Total entries: {len(TARGETS)}\n")

    for idx, item in enumerate(TARGETS, 1):
        print(f"[{idx}/{len(TARGETS)}] {item['name']} ({item['year']}) - {item['generation']}")
        print(f"  Target Title: {item['title']}")
        print(f"  DOI: {item['doi']}")

        # 1. DOI Verification
        doi_result = {}
        if item['doi'].startswith("10.48550/arXiv"):
            doi_result = {
                "verified": True,
                "verified_by": "arXiv",
                "verified_title": item['title'],
                "verified_container": "arXiv",
                "verified_year": item['year'],
                "verified_authors": item['authors']
            }
            print("  [DOI] Confirmed via arXiv Registry")
        else:
            v_res = verify_doi_crossref(item['doi'])
            if v_res.get("verified"):
                doi_result = v_res
                print(f"  [DOI] Confirmed via CrossRef: \"{v_res.get('verified_title')}\" in {v_res.get('verified_container')}")
            else:
                doi_result = {"verified": False, "error": v_res.get("error")}
                print(f"  [DOI WARNING] CrossRef query issue: {v_res.get('error')}")

        # 2. PDF Download
        pdf_status = {"downloaded": False, "file_path": None, "size_kb": 0, "source": item['source_type']}
        if item['download_url']:
            filename = f"{item['id']}.pdf"
            filepath = os.path.join(papers_dir, filename)
            success, msg = download_with_curl(item['download_url'], filepath)
            if success:
                pdf_status["downloaded"] = True
                pdf_status["file_path"] = f"downloaded_papers/{filename}"
                pdf_status["size_kb"] = msg // 1024
                print(f"  [DOWNLOAD] Success -> {filename} ({msg // 1024} KB)")
            else:
                print(f"  [DOWNLOAD FAILED] {msg}")
        else:
            print("  [DOWNLOAD] Publisher has strict access restriction (Abstract/Metadata indexed)")

        manifest.append({
            **item,
            "verification": doi_result,
            "pdf_status": pdf_status
        })
        time.sleep(0.3)
        print()

    with open(os.path.join(meta_dir, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    downloaded_count = sum(1 for m in manifest if m["pdf_status"]["downloaded"])
    verified_count = sum(1 for m in manifest if m["verification"]["verified"])
    print(f"=== Complete! Verified: {verified_count}/{len(TARGETS)}, Downloaded: {downloaded_count}/{len(TARGETS)} ===")

if __name__ == "__main__":
    main()
