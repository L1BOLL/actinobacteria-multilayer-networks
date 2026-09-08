#!/usr/bin/env python3
"""Regenerate every supplementary table and figure from the raw inputs.

    pip install -r requirements.txt
    python run_all.py

Stages are independent and are reported individually. A stage whose optional
input data is absent (currently only the spatial decomposition) records a skip in
supplementary/report.md and the run continues, so one missing input never costs
you the rest of the outputs.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parent / "code"

# prep_phylogeny and build_tree produce data/phylogeny/tree.nwk, which p2_8 consumes.
STAGES = [
    "prep_phylogeny.py",
    "build_tree.py",
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
    "p2_13_taxon_sensitivity.py",
    "p2_14_reciprocity_transitivity_null.py",
    "generate_figs.py",
    "make_paper_figures.py",
    "figure_2_spatial.py",
    "figure_3_overlap.py",
    "figure_4_strategy.py",
    "tables_1_2.py",
    "verify_reported_values.py",
]


def main() -> None:
    results: list[tuple[str, int]] = []
    for stage in STAGES:
        print(f"\n=== {stage} ===", flush=True)
        rc = subprocess.run([sys.executable, str(CODE / stage)], cwd=CODE).returncode
        results.append((stage, rc))
        if rc != 0:
            print(f"!!! {stage} exited {rc}; continuing", flush=True)

    print("\n=== summary ===")
    for stage, rc in results:
        print(f"  {'ok  ' if rc == 0 else 'FAIL'}  {stage}")

    failed = [s for s, rc in results if rc != 0]
    if failed:
        print(f"\n{len(failed)} of {len(results)} stages failed: {', '.join(failed)}")
        sys.exit(1)
    print(f"\nall {len(results)} stages completed")


if __name__ == "__main__":
    main()
