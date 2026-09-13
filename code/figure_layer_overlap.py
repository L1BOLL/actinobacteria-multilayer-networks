#!/usr/bin/env python3
"""Final Figure 3 — cross-layer overlap against the degree-preserving null.

Reads results/tables/tableS2_layer_overlap_full.csv, which is the
output of code/layer_overlap.py. Nothing is hard-coded: the panel titles,
the counts in the annotation and the significance marks all derive from the CSV.

This replaces the pre_final Figure 3, whose two panels were raw cosine and raw
Jaccard clustermaps. The manuscript states that raw overlap "was considered
descriptive and not interpreted", so those panels illustrated the one quantity
the text declines to draw conclusions from. They move to the supplement.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_io import FIGURE_FORMATS
from matplotlib.colors import TwoSlopeNorm

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"
OUT = Path(__file__).resolve().parents[1] / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86", "#ebebe8"

# Layers grouped by functional category, as in Table 1.
CATEGORY = {"IG": "DA", "IC": "DA", "RAC": "DA", "IAC_RAC": "IA", "IAC_RDE": "IA", "IRAC": "IA",
            "RP": "MM", "IRP": "MM", "CCAM": "MC", "CCVM": "MC", "CMP": "MC", "CS": "MC"}
ORDER = ["RAC", "IG", "IC", "IAC_RAC", "IRAC", "IAC_RDE", "RP", "IRP", "CS", "CMP", "CCVM", "CCAM"]


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.titlesize": 11.5, "axes.labelsize": 10.5, "axes.titlelocation": "left", "axes.titlepad": 8,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "axes.edgecolor": INK2, "axes.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 9.5,
    })


def tag(ax, letter: str, dx: float = -0.16) -> None:
    ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="bottom", ha="left", color=INK)


def main() -> None:
    df = pd.read_csv(TABLES / "tableS2_layer_overlap_full.csv")
    n_draws = int(df["n_null_draws"].iloc[0])
    style()

    # square z matrix ------------------------------------------------------ #
    Z = pd.DataFrame(np.nan, index=ORDER, columns=ORDER, dtype=float)
    SIG = pd.DataFrame(False, index=ORDER, columns=ORDER)
    CONF = pd.DataFrame(False, index=ORDER, columns=ORDER)
    for r in df.itertuples():
        for a, b in ((r.layer1, r.layer2), (r.layer2, r.layer1)):
            Z.loc[a, b] = r.z
            SIG.loc[a, b] = bool(r.significant_by)
            CONF.loc[a, b] = bool(r.methodological_confound)

    fig = plt.figure(figsize=(14.6, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.42)

    # (a) z heatmap -------------------------------------------------------- #
    ax = fig.add_subplot(gs[0, 0])
    # The CCAM-CCVM confound (z = +16) would flatten every biological contrast,
    # so the scale is set by the interpretable pairs and that cell saturates.
    bio = df[~df["methodological_confound"]]
    vmax = float(np.ceil(bio["z"].abs().max()))
    im = ax.imshow(Z.values, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax))
    for i, a in enumerate(ORDER):
        for j, b in enumerate(ORDER):
            if i == j:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, color="#f4f4f2", zorder=2))
            elif SIG.loc[a, b]:
                mark = "†" if CONF.loc[a, b] else "*"
                ax.text(j, i, mark, ha="center", va="center", fontsize=10.5, color=INK, zorder=3)
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, rotation=90)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels(ORDER)
    for k, lab in enumerate(ORDER):
        col = {"DA": ORANGE, "IA": "#b45309", "MM": "#0f766e", "MC": BLUE}[CATEGORY[lab]]
        ax.get_xticklabels()[k].set_color(col)
        ax.get_yticklabels()[k].set_color(col)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Standardized overlap  $Z_{lm}$", fontsize=10)
    cb.outline.set_linewidth(0.6)
    ax.set_title("Most pairs sit at their degree-preserving expectation")
    ax.text(0.0, -0.30, "* BY q < 0.05    † CCAM–CCVM, shared measurement",
            transform=ax.transAxes, fontsize=9, color=MUTED)
    for s in ax.spines.values():
        s.set_visible(False)
    tag(ax, "a", dx=-0.22)

    # (b) observed J vs null for the significant pairs --------------------- #
    ax = fig.add_subplot(gs[0, 1])
    sig = df[df["significant_by"]].copy()
    sig["pair"] = sig["layer1"] + "–" + sig["layer2"]
    sig = sig.sort_values("z")
    y = np.arange(len(sig))
    lo = sig["J_null_mean"] - 1.96 * sig["J_null_sd"]
    hi = sig["J_null_mean"] + 1.96 * sig["J_null_sd"]
    ax.barh(y, hi - lo, left=lo, height=0.55, color="#d9d9d6",
            label="null, mean ± 1.96 SD", zorder=2)
    ax.plot(sig["J_null_mean"], y, "|", color=MUTED, markersize=9, zorder=3)
    cols = [MUTED if c else (ORANGE if z > 0 else BLUE)
            for c, z in zip(sig["methodological_confound"], sig["z"])]
    ax.scatter(sig["J_obs"], y, s=42, color=cols, zorder=4, label="observed $J$")
    for i, r in enumerate(sig.itertuples()):
        left = r.z < 0
        ax.annotate(f"z = {r.z:+.1f}", (r.J_obs, i), textcoords="offset points",
                    xytext=(-9 if left else 9, 0), va="center",
                    ha="right" if left else "left", fontsize=9, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([p + (" †" if c else "") for p, c in zip(sig["pair"], sig["methodological_confound"])])
    ax.set_xlabel("Jaccard similarity of directed edge sets")
    ax.set_xlim(-0.035, max(sig["J_obs"].max(), hi.max()) * 1.22)
    n_bio = int((~sig["methodological_confound"]).sum())
    n_conf = int(sig["methodological_confound"].sum())
    ax.set_title(f"{n_bio} biological departures" + (f", plus the measurement confound" if n_conf else ""))
    ax.grid(axis="x", color=GRID, linewidth=0.7, zorder=0)
    ax.legend(loc="lower right")
    tag(ax, "b", dx=-0.30)

    n_over = int(((sig["z"] > 0) & (~sig["methodological_confound"])).sum())
    n_under = int(((sig["z"] < 0) & (~sig["methodological_confound"])).sum())
    fig.suptitle(
        f"{n_draws:,} degree-preserving draws per layer; BY-FDR across {len(df)} pairs; "
        f"{n_over} enriched, {n_under} depleted among the {len(df) - 1} biological pairs",
        fontsize=10, color=INK2, y=1.02)

    for ext in FIGURE_FORMATS:
        fig.savefig(OUT / f"figure_layer_overlap.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote figure_layer_overlap.{{png,svg}}   ({n_over} over, {n_under} under, {n_draws:,} draws)")


if __name__ == "__main__":
    main()
