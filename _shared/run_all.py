#!/usr/bin/env python3
from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    "p2_2_gmm.py",
    "p2_3_pca12.py",
    "p2_4_jaccard_null.py",
    "p2_5_scc_null.py",
    "p2_7_proximity_sensitivity.py",
    "p2_8_phylogeny_pgls.py",
    "p2_9_spatial_decomposition.py",
    "p2_10_directionality.py",
    "p2_11_dominance.py",
    "p2_12_metabolism.py",
]


def main() -> None:
    for script in SCRIPTS:
        subprocess.run(["python", str(ROOT / script)], check=True)


if __name__ == "__main__":
    main()

