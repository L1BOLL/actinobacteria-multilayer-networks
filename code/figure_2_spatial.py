#!/usr/bin/env python3
"""Final Figure 2 — phenotype-specific spatial dependence of IG and CS.

    !!  INHERITED INPUTS  !!

Every number below is transcribed from `pre_final/Paper Coexistence_final_v1.0.docx`.
None of it is recomputed, because the per-configuration matrices this analysis ran
on are not in BN_paper_methods/ — `p2_9_spatial_decomposition.py` skips for want of
`data/matrices_xls_direct/`, `_indirect/` and `_distant/`. Drop those in, run p2_9,
and this script should be replaced by one that reads its output.

What this script *does* verify, since the counts fully determine them:
  - density        = edges / 3540           (3540 = 60 x 59 ordered non-self pairs)
  - odds ratio     = lost / gained
  - delta density  = density_close - density_separated
Mismatches are raised, not silently plotted. One is already known and corrected
here: the manuscript states IG close-proximity density as 0.217, which is the
*main-layer* IG density (767 edges, Table 1). The spatial arm has 785 edges, i.e.
0.222 — and 0.222 is the value its own reported delta of 0.148 is built from.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from data_io import FIGURE_FORMATS

OUT = Path(__file__).resolve().parents[1] / "supplementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86", "#ebebe8"
N_PAIRS = 60 * 59  # 3540 ordered non-self strain pairs

# --------------------------------------------------------------------------- #
# inherited constants
# --------------------------------------------------------------------------- #
SPATIAL = {
    "IG": {
        "edges_close": 785,
        "edges_separated": 261,
        "persistent": 249,
        "lost": 536,
        "gained": 12,
        "delta_density": 0.148,
        "delta_ci": (0.136, 0.160),
        "odds_ratio": 44.67,
        "or_ci": (25.21, 79.15),
        "p_perm": 0.0001,
        "p_mcnemar_text": "$P < 10^{-100}$",
    },
    "CS": {
        "edges_close": 144,
        "edges_separated": 133,
        "persistent": 59,
        "lost": 85,
        "gained": 74,
        "delta_density": 0.0031,
        "delta_ci": (-0.0040, 0.0102),
        "odds_ratio": 1.15,
        "or_ci": (0.84, 1.57),
        "p_perm": 0.432,
        "p_mcnemar_text": "$P = 0.43$",
    },
}


def check() -> None:
    """Fail loudly rather than plot a figure that disagrees with its own numbers."""
    for layer, d in SPATIAL.items():
        dens_close = d["edges_close"] / N_PAIRS
        dens_sep = d["edges_separated"] / N_PAIRS
        assert d["persistent"] + d["lost"] == d["edges_close"], f"{layer}: persistent+lost != close edges"
        assert d["persistent"] + d["gained"] == d["edges_separated"], f"{layer}: persistent+gained != separated edges"
        assert abs(d["lost"] / d["gained"] - d["odds_ratio"]) < 0.01, f"{layer}: lost/gained != reported OR"
        assert abs((dens_close - dens_sep) - d["delta_density"]) < 5e-4, (
            f"{layer}: density difference {dens_close - dens_sep:.4f} != reported delta {d['delta_density']}"
        )
        d["density_close"], d["density_separated"] = dens_close, dens_sep
        print(f"  {layer}: close {dens_close:.4f}  separated {dens_sep:.4f}  delta {dens_close - dens_sep:.4f}  OK")


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans", "font.size": 8.5,
        "axes.titlesize": 9.5, "axes.labelsize": 8.5, "axes.titlelocation": "left", "axes.titlepad": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": INK2, "axes.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 7.5,
    })


def tag(ax, letter: str) -> None:
    ax.text(-0.18, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="left", color=INK)


def main() -> None:
    print("checking inherited constants against each other")
    check()
    style()
    layers = ["IG", "CS"]
    colors = {"IG": ORANGE, "CS": BLUE}

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.5))

    # (a) density under each configuration --------------------------------- #
    ax = axes[0]
    for i, L in enumerate(layers):
        d = SPATIAL[L]
        ax.plot([0, 1], [d["density_close"], d["density_separated"]], "-o", color=colors[L],
                markersize=6, linewidth=1.8, label=L, zorder=3)
        ax.annotate(f"{d['edges_close']}", (0, d["density_close"]), textcoords="offset points",
                    xytext=(-6, 6), ha="right", fontsize=7.5, color=colors[L])
        ax.annotate(f"{d['edges_separated']}", (1, d["density_separated"]), textcoords="offset points",
                    xytext=(6, 6), ha="left", fontsize=7.5, color=colors[L])
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Close\nproximity", "Spatially\nseparated"])
    ax.set_xlim(-0.35, 1.35); ax.set_ylim(0, 0.25)
    ax.set_ylabel("Directed network density")
    ax.set_title("Density collapses for IG, not for CS")
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.legend(loc="center right")
    tag(ax, "a")

    # (b) fate of the close-proximity edge set ----------------------------- #
    ax = axes[1]
    cats = ["persistent", "lost", "gained"]
    labels = ["Persistent", "Lost after\nseparation", "Gained after\nseparation"]
    w, x = 0.36, np.arange(3)
    for i, L in enumerate(layers):
        vals = [SPATIAL[L][c] for c in cats]
        bars = ax.bar(x + (i - 0.5) * w, vals, width=w, color=colors[L], label=L, zorder=3)
        ax.bar_label(bars, fontsize=7, padding=2, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Ordered strain pairs")
    ax.set_ylim(0, 620)
    ax.set_title("IG loses edges; CS exchanges them")
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.legend()
    tag(ax, "b")

    # (c) matched-pair odds ratio ------------------------------------------ #
    ax = axes[2]
    for i, L in enumerate(layers):
        d = SPATIAL[L]
        y = 1 - i
        ax.plot([d["or_ci"][0], d["or_ci"][1]], [y, y], color=colors[L], linewidth=2.0, solid_capstyle="round")
        ax.plot(d["odds_ratio"], y, "o", color=colors[L], markersize=7, zorder=3)
        ax.annotate(f"OR = {d['odds_ratio']:.2f}  ({d['or_ci'][0]:.2f}–{d['or_ci'][1]:.2f})",
                    (d["odds_ratio"], y), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=7.5, color=INK2)
        ax.annotate(d["p_mcnemar_text"], (d["odds_ratio"], y), textcoords="offset points",
                    xytext=(0, -17), ha="center", fontsize=7.5, color=MUTED)
    ax.axvline(1.0, color=INK2, linewidth=0.9, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlim(0.5, 200); ax.set_ylim(-0.75, 1.75)
    ax.set_yticks([1, 0]); ax.set_yticklabels(layers)
    ax.set_xlabel("Matched-pair odds ratio, loss : gain (log scale)")
    ax.set_title("Only IG departs from unity")
    ax.grid(axis="x", color=GRID, linewidth=0.7, zorder=0)
    tag(ax, "c")

    fig.tight_layout()
    for ext in FIGURE_FORMATS:
        fig.savefig(OUT / f"figure2_spatial_dependence.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT/'figure2_spatial_dependence.png'} / .pdf")


if __name__ == "__main__":
    main()
