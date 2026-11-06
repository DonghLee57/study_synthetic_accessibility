#!/usr/bin/env python3
"""
Automated Setup & Cloning of External Benchmark Repositories
-----------------------------------------------------------
Clones and pins external model repositories to exact audited commits in external/,
and applies compatibility patches for scscore and syba.
"""

import os
import subprocess
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTERNAL_DIR = os.path.join(PROJECT_ROOT, "external")

REPOSITORIES = {
    "scscore": {
        "url": "https://github.com/connorcoley/scscore.git",
        "commit": "37090a6aa8220408f572de20224c33c221312d16",
        "branch": "master"
    },
    "syba": {
        "url": "https://github.com/lich-uct/syba.git",
        "commit": "30929eb2d3b9a03bef326a1f5180f26f55b19a29",
        "branch": "master"
    },
    "rascore": {
        "url": "https://github.com/reymond-group/RAscore.git",
        "commit": "cb77db503ee5cbf0e8bb8963df6e5b76b3a94f06",
        "branch": "master"
    },
    "rxnflow": {
        "url": "https://github.com/SeonghwanSeo/RxnFlow.git",
        "commit": "530d098b56625381395fe7a30a2a10b7aacb8077",
        "branch": "main"
    },
    "reasyn": {
        "url": "https://github.com/NVIDIA-BioNeMo/ReaSyn.git",
        "commit": "f28270ca3a5873b2bf7eb52821c5202b111ef9f4",
        "branch": "reasyn_v2"
    },
    "synnet": {
        "url": "https://github.com/wenhao-gao/SynNet.git",
        "commit": "3ccf4fcd6eb0605cb60f6376ad7b4a7accaf1549",
        "branch": "master"
    },
    "chemprojector": {
        "url": "https://github.com/luost26/ChemProjector.git",
        "commit": "dd3afda12d33437f9c41bc7b246903628e60429a",
        "branch": "main"
    },
    "synformer": {
        "url": "https://github.com/wenhao-gao/synformer.git",
        "commit": "bef02bc4481adf3ed792090fc6169331d4d0553a",
        "branch": "main"
    },
    "retro_star": {
        "url": "https://github.com/binghong-ml/retro_star.git",
        "commit": "60e7a24213810706a7075d113704c9320158d726",
        "branch": "master"
    },
    "dfrscore": {
        "url": "https://github.com/Hwoo-Kim/DFRscore.git",
        "commit": "5b1a4a4157105089bb3c818fd4a136fe255c80b8",
        "branch": "main"
    },
    "deepsa": {
        "url": "https://github.com/Shihang-Wang-58/DeepSA.git",
        "commit": "eaaa6a52e7edd82fb4ec3823e45c6f6f373e8ef2",
        "branch": "main"
    },
    "gasa": {
        "url": "https://github.com/cadd-synthetic/GASA.git",
        "commit": "b4f46553da749c4dc474a7fe1fc6041770d44659",
        "branch": "master"
    }
}


def setup_repos():
    os.makedirs(EXTERNAL_DIR, exist_ok=True)
    for name, info in REPOSITORIES.items():
        repo_path = os.path.join(EXTERNAL_DIR, name)
        if not os.path.exists(repo_path):
            print(f"[CLONE] Cloning {name} from {info['url']}...")
            subprocess.run(["git", "clone", info["url"], repo_path], check=True)
            print(f"[CHECKOUT] Pinning {name} to commit {info['commit']}...")
            subprocess.run(["git", "-C", repo_path, "checkout", info["commit"]], check=True)
        else:
            print(f"[EXISTS] {name} already present at {repo_path}.")

    # Apply patches
    patch_script = os.path.join(PROJECT_ROOT, "scripts", "patch_external_models.py")
    if os.path.isfile(patch_script):
        subprocess.run([sys.executable, patch_script], check=True)


if __name__ == "__main__":
    setup_repos()
