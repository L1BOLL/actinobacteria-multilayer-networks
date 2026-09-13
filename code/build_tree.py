#!/usr/bin/env python3
"""Align the cleaned 16S and infer a maximum-likelihood tree.

Pip-only, no external binaries and no sudo, so the supplement stays standalone:

  MAFFT  (pymafft)      alignment            -> data/phylogeny/16S_aligned.fasta
  trimAl (pytrimal)     trim ragged ends     -> data/phylogeny/16S_trimmed.fasta
  VeryFastTree          GTR+GAMMA ML tree    -> data/phylogeny/tree.nwk

VeryFastTree is a validated reimplementation of FastTree 2 producing equivalent
trees, used here because RAxML has no pip distribution. The tool names in the
manuscript Methods must match what is actually run here.

If a tree produced by the coauthor's own MAFFT + RAxML pipeline is dropped in at
data/phylogeny/tree.nwk, this script leaves it alone (unless --force).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from data_io import ROOT, SEED

# Strain IDs contain spaces ("MS1 13"). Newick has no way to carry those, and
# FastTree truncates names at the first whitespace -- which silently collides
# "MS1 13" and "MS1 14" into a single tip called "MS1". Tips are therefore
# written with underscores and mapped back when the tree is consumed.
TREE_LABEL_SEP = "_"


def to_tree_label(name: str) -> str:
    return name.replace(" ", TREE_LABEL_SEP)


def from_tree_label(label: str) -> str:
    return label.replace(TREE_LABEL_SEP, " ")


PHYLO_DIR = ROOT / "data" / "phylogeny"
FASTA_IN = PHYLO_DIR / "16S_clean.fasta"
ALIGNED = PHYLO_DIR / "16S_aligned.fasta"
TRIMMED = PHYLO_DIR / "16S_trimmed.fasta"
TREE_OUT = PHYLO_DIR / "tree.nwk"
PROVENANCE = PHYLO_DIR / "tree_provenance.txt"


def align(fasta_text: str) -> str:
    """MAFFT alignment; falls back to FAMSA if the MAFFT bindings misbehave."""
    try:
        import pymafft

        # L-INS-i-style accuracy settings: 16S at within-genus divergence is
        # nearly ungapped, so refinement iterations are cheap and worth it.
        result = pymafft.align_fasta_string(fasta_text, strategy="fftnsi", maxiterate=1000)
        return result.to_fasta(), "MAFFT (pymafft, fftnsi, maxiterate=1000)"
    except Exception as exc:  # noqa: BLE001 - fall back rather than abort the pipeline
        print(f"  pymafft failed ({exc}); falling back to FAMSA", file=sys.stderr)
        from pyfamsa import Aligner, Sequence

        records = []
        name = None
        buf: list[str] = []
        for line in fasta_text.splitlines():
            if line.startswith(">"):
                if name:
                    records.append(Sequence(name.encode(), "".join(buf).encode()))
                name, buf = line[1:].strip(), []
            elif name:
                buf.append(line.strip())
        if name:
            records.append(Sequence(name.encode(), "".join(buf).encode()))

        aligned = Aligner(guide_tree="upgma").align(records)
        out = [f">{s.id.decode()}\n{s.sequence.decode()}" for s in aligned]
        return "\n".join(out) + "\n", "FAMSA (pyfamsa, UPGMA guide tree)"


def trim(aligned_path: Path, out_path: Path) -> str:
    """Strip ragged Sanger ends, which otherwise dominate the distance signal."""
    from pytrimal import Alignment, AutomaticTrimmer

    alignment = Alignment.load(str(aligned_path))
    trimmed = AutomaticTrimmer(method="gappyout").trim(alignment)
    trimmed.dump(str(out_path), format="fasta")
    return f"trimAl gappyout: {alignment.sequences[0].__len__()} -> {len(trimmed.sequences[0])} columns"


def infer_tree(trimmed_path: Path, out_path: Path) -> str:
    import veryfasttree

    # nt nucleotide, gtr GTR model, gamma GAMMA rate heterogeneity,
    # seed for reproducibility of the SPR/NNI search.
    veryfasttree.run(
        str(trimmed_path),
        nt=True,
        gtr=True,
        gamma=True,
        seed=SEED % 32767,
        out=str(out_path),
        nopr=True,
    )
    return "VeryFastTree -nt -gtr -gamma (FastTree-2 equivalent ML, GTR+GAMMA)"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="rebuild even if tree.nwk exists")
    args = parser.parse_args()

    PHYLO_DIR.mkdir(parents=True, exist_ok=True)

    if not FASTA_IN.exists():
        print(f"missing {FASTA_IN}; run prep_phylogeny.py first")
        return

    if TREE_OUT.exists() and not args.force:
        print(f"{TREE_OUT.relative_to(ROOT)} already exists; leaving it alone (--force to rebuild)")
        return

    fasta_text = FASTA_IN.read_text(encoding="utf-8")
    n_seq = fasta_text.count(">")
    fasta_text = "\n".join(
        f">{to_tree_label(line[1:].strip())}" if line.startswith(">") else line
        for line in fasta_text.splitlines()
    ) + "\n"
    print(f"aligning {n_seq} sequences...")
    aligned_text, aligner_note = align(fasta_text)
    ALIGNED.write_text(aligned_text, encoding="utf-8")
    print(f"  {aligner_note}")

    trim_note = trim(ALIGNED, TRIMMED)
    print(f"  {trim_note}")

    tree_note = infer_tree(TRIMMED, TREE_OUT)
    print(f"  {tree_note}")

    PROVENANCE.write_text(
        "Tree provenance (regenerate with: python code/build_tree.py --force)\n\n"
        f"sequences  : {n_seq} (from {FASTA_IN.name})\n"
        f"alignment  : {aligner_note}\n"
        f"trimming   : {trim_note}\n"
        f"inference  : {tree_note}\n"
        f"seed       : {SEED}\n\n"
        "These are the tool names that must appear in the manuscript Methods.\n"
        "If the coauthor supplies a MAFFT + RAxML tree, overwrite tree.nwk and\n"
        "update this file; phylogenetic_signal.py consumes it unchanged.\n",
        encoding="utf-8",
    )
    print(f"wrote {TREE_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
