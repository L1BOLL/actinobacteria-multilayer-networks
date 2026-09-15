#!/usr/bin/env python3
"""Regenerate every table and figure from the raw inputs.

    pip install -r requirements.txt
    python run_all.py

Stages are independent and are reported individually. A stage whose optional input
is absent records a skip in results/report.md and the run continues, so one
missing input never costs you the rest of the outputs.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parent / "code"

# prep_phylogeny and build_tree produce data/phylogeny/tree.nwk, which
# phylogenetic_signal.py consumes.
STAGES = [
    "prep_phylogeny.py",
    "build_tree.py",
    "mixture_models.py",
    "ordination.py",
    "layer_overlap.py",
    "connectivity.py",
    "layer_descriptors.py",
    "phylogenetic_signal.py",
    "directionality.py",
    "dominance.py",
    "metabolic_niche.py",
    "taxon_sensitivity.py",
    "reciprocity_transitivity.py",
    "spatial_configuration.py",
    "figures_architecture.py",
    "figures_phylogeny.py",
    "figure_plates.py",
    "figure_phenotype_gallery.py",
    "figure_layer_overlap.py",
    "figure_strategy_space.py",
    "figure_role_stability.py",
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
