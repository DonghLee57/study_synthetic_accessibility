#!/usr/bin/env python3
"""
Patch External Models for Modern Python 3.11 & NumPy 1.26 Compatibility
------------------------------------------------------------------------
This script applies necessary compatibility patches to cloned external repositories:
1. external/scscore:
   - Modifies `sigmoid(x)` in `scscore/standalone_model_numpy.py` to handle numpy arrays via `np.squeeze(x)`.
   - Replaces deprecated `np.bool` with `bool` for NumPy 1.24+ compatibility.
   - Replaces Python 2 `import cPickle as pickle` with `import pickle`.
2. external/syba:
   - Configures `fitDefaultScore()` in `syba/syba.py` to load `syba.csv.gz` (~25 MB, 4 s)
     instead of `syba4.csv.gz` (~115 MB, 45 s), reproducing the audited benchmark metrics.
"""

import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def patch_scscore():
    file_path = os.path.join(PROJECT_ROOT, "external", "scscore", "scscore", "standalone_model_numpy.py")
    if not os.path.isfile(file_path):
        print(f"[SKIP] scscore file not found: {file_path}")
        return

    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    modified = False

    # 1. sigmoid
    if "return 1 / (1 + math.exp(-x))" in content:
        content = content.replace(
            "return 1 / (1 + math.exp(-x))",
            "return 1.0 / (1.0 + math.exp(-float(np.squeeze(x))))"
        )
        modified = True
        print("[PATCH] scscore: patched sigmoid() to handle numpy scalar/array.")

    # 2. np.bool deprecation
    if "dtype=np.bool)" in content:
        content = content.replace("dtype=np.bool)", "dtype=bool)")
        modified = True
        print("[PATCH] scscore: replaced np.bool with builtin bool.")

    # 3. cPickle import
    if "import cPickle as pickle" in content:
        content = content.replace("import cPickle as pickle", "import pickle")
        modified = True
        print("[PATCH] scscore: replaced cPickle with standard pickle.")

    if modified:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        print("[OK] external/scscore patched successfully.")
    else:
        print("[OK] external/scscore already up-to-date.")


def patch_syba():
    file_path = os.path.join(PROJECT_ROOT, "external", "syba", "syba", "syba.py")
    if not os.path.isfile(file_path):
        print(f"[SKIP] syba file not found: {file_path}")
        return

    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    modified = False
    if 'os.path.join(this_dir, "resources", "syba4.csv.gz")' in content:
        content = content.replace(
            'os.path.join(this_dir, "resources", "syba4.csv.gz")',
            'os.path.join(this_dir, "resources", "syba.csv.gz")'
        )
        modified = True
        print("[PATCH] syba: set default score table to syba.csv.gz.")

    if modified:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        print("[OK] external/syba patched successfully.")
    else:
        print("[OK] external/syba already up-to-date.")


if __name__ == "__main__":
    print("=== Patching External Repositories for Modern Environment Compatibility ===")
    patch_scscore()
    patch_syba()
    print("=== All External Patches Verified ===")
