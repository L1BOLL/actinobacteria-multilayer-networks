#!/usr/bin/env python3
"""Figure S13 - gallery of confrontation plates, one column per phenotype family.

The source is a single photograph montage, `data/plates/phenotype_gallery.png`,
laid out as 4 rows x 6 columns of plates. This script only adds the column
headers; it does not crop or rearrange the tiles, so replacing the source file
with a higher-resolution scan of the same montage regenerates the figure.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "plates" / "phenotype_gallery.png"
OUT = ROOT / "results" / "figures"

COLUMNS = ["IG", "IAC_RAC", "IRAC", "IC", "RP / IRP", "CCAM / CCVM / CS"]
HEADER_FONT = 7.5
N_ROWS = 4


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    im = Image.open(SRC).convert("RGB")
    w, h = im.size
    fig_w = 7.2
    fig_h = fig_w * h / w + 0.35
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
    ax = fig.add_axes([0, 0, 1, h / w * fig_w / fig_h])
    ax.imshow(im, interpolation="lanczos")
    ax.set_axis_off()
    hdr = fig.add_axes([0, h / w * fig_w / fig_h, 1, 1 - h / w * fig_w / fig_h])
    hdr.set_axis_off()
    for k, label in enumerate(COLUMNS):
        hdr.text((k + 0.5) / len(COLUMNS), 0.45, label, ha="center", va="center",
                 fontsize=HEADER_FONT, fontweight="bold", family="DejaVu Sans")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"figure_phenotype_gallery.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"source {w}x{h} px; columns {COLUMNS}")


if __name__ == "__main__":
    main()
