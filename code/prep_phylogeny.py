#!/usr/bin/env python3
"""Sanitize the 16S FASTA and join it to the interaction-matrix strain set.

The delivered FASTA has three defects that silently corrupt downstream analysis
if left alone (all logged to data/phylogeny/label_map.csv):

1. `MS8 6` carries ~77 bp of sequence concatenated onto its header line, so those
   bases are dropped by any naive parser.
2. `S705A` appears twice; the two records differ in length (1390 vs 1468 bp).
3. Tip labels use non-breaking spaces and drop the apostrophes present in the
   .xlsx indices (`S713'V` vs `S713V`).

Emits data/phylogeny/16S_clean.fasta with tip labels in the canonical key space.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import pandas as pd

from data_io import ROOT, canonical_strain, ensure_output_dirs, load_tensor

FASTA_IN = ROOT / "data" / "annotation" / "16S_strepto_R"
PHYLO_DIR = ROOT / "data" / "phylogeny"
FASTA_OUT = PHYLO_DIR / "16S_clean.fasta"
LABEL_MAP = PHYLO_DIR / "label_map.csv"

# Genus is parsed from the FASTA description; these four are not Streptomyces and
# drive the taxonomy wording fix plus the n=56 sensitivity re-run (taxon_sensitivity).
NON_STREPTOMYCES = {"Saccharothrix", "Lentzea", "Amycolatopsis"}

VALID_DNA = set("ACGTRYSWKMBDHVN")


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    buf: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(buf)))
            header, buf = line[1:], []
        elif header is not None:
            buf.append(line.strip())
    if header is not None:
        records.append((header, "".join(buf)))
    return records


def split_runon_header(header: str) -> tuple[str, str]:
    """Recover sequence that was concatenated onto a header line (defect 1)."""
    match = re.search(r"\s([ACGTRYSWKMBDHVN]{20,})\s*$", header)
    if not match:
        return header, ""
    return header[: match.start()].rstrip(), match.group(1)


def parse_label(header: str) -> tuple[str, str]:
    """Return (genus, strain_label) from a FASTA description line."""
    # Some records are prefixed with a GenBank accession (e.g. "OP131863.1 ...").
    text = re.sub(r"^[A-Z]{2}\d+\.\d+\s+", "", header).strip()
    genus = text.split()[0] if text.split() else ""
    label = text.split(" strain ", 1)[1] if " strain " in text else text
    return genus, canonical_strain(label)


def main() -> None:
    ensure_output_dirs()
    PHYLO_DIR.mkdir(parents=True, exist_ok=True)

    _, _, node_ids = load_tensor()
    matrix_keys = {canonical_strain(n): n for n in node_ids}

    rows: list[dict[str, object]] = []
    best: dict[str, tuple[str, str, str]] = {}  # key -> (genus, label, seq)

    for header, seq in parse_fasta(FASTA_IN):
        clean_header, recovered = split_runon_header(header)
        seq = (recovered + seq).upper()
        seq = "".join(c for c in seq if c in VALID_DNA)
        genus, label = parse_label(clean_header)

        note = []
        if recovered:
            note.append(f"recovered {len(recovered)} bp from header")
        if genus in NON_STREPTOMYCES:
            note.append(f"non-Streptomyces ({genus})")
        if label not in matrix_keys:
            note.append("no interaction matrix")

        kept = True
        if label in best:
            # Defect 2: keep the longer record, log the discard.
            prev = best[label]
            if len(seq) > len(prev[2]):
                note.append(f"duplicate: replaces shorter {len(prev[2])} bp record")
                for r in rows:
                    if r["strain"] == label and r["kept"]:
                        r["kept"] = False
                        r["note"] = (str(r["note"]) + "; superseded by longer record").strip("; ")
            else:
                note.append(f"duplicate: discarded, shorter than kept {len(prev[2])} bp record")
                kept = False

        if kept:
            best[label] = (genus, label, seq)

        rows.append({
            "fasta_header": header[:80],
            "strain": label,
            "genus": genus,
            "length_bp": len(seq),
            "in_interaction_matrix": label in matrix_keys,
            "kept": kept,
            "note": "; ".join(note),
        })

    # Restrict the tree to strains that also have interaction data; a tip with no
    # traits contributes nothing to the lambda test and only adds topology noise.
    usable = {k: v for k, v in best.items() if k in matrix_keys}
    with FASTA_OUT.open("w", encoding="utf-8") as fh:
        for key in sorted(usable):
            fh.write(f">{key}\n")
            seq = usable[key][2]
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")

    pd.DataFrame(rows).to_csv(LABEL_MAP, index=False)

    unmatched_fasta = sorted(set(best) - set(matrix_keys))
    unmatched_matrix = sorted(set(matrix_keys) - set(best))
    non_strep = sorted(k for k, v in usable.items() if v[0] in NON_STREPTOMYCES)

    print(f"records read           : {len(rows)}")
    print(f"unique strains kept    : {len(best)}")
    print(f"written to tree input  : {len(usable)}  -> {FASTA_OUT.relative_to(ROOT)}")
    print(f"in FASTA, no matrix    : {unmatched_fasta}")
    print(f"in matrix, no FASTA    : {unmatched_matrix}")
    print(f"non-Streptomyces kept  : {non_strep}")
    print(f"label audit            : {LABEL_MAP.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
