#!/usr/bin/env python3
"""Strain interaction strategies in the twelve-layer space.

Five panels: the ordination of the twelve-layer outgoing profile, its loadings,
the variance spectrum against the parallel-analysis null, the out-degree against
participation continuum, and the held-out likelihood by mixture component count.

Retention uses the 95th percentile of the parallel-analysis null rather than its
mean, so a component whose variance clears only the mean is drawn as marginal.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from data_io import FIGURE_FORMATS, LAYER_ORDER, SEED, load_tensor

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "figures"
TAB = REPO / "results" / "tables"
BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86", "#ebebe8"
GREY = "#b9b9b4"


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans",
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
        "axes.titlelocation": "left", "axes.titlepad": 8,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.edgecolor": INK2, "axes.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 8.5,
    })


def cbar(fig, ax, mappable, label):
    """Colorbar in an axis carved out of the parent, so it cannot be relocated
    by the layout engine into a neighbouring panel."""
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.08)
    cb = fig.colorbar(mappable, cax=cax)
    cb.set_label(label, fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(0.6)
    return cb


def panel_tag(ax, letter):
    ax.text(-0.17, 1.06, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="left", color=INK)


def main() -> None:
    style()
    tensor, layer_ids, nodes = load_tensor()
    out_deg = tensor.sum(axis=1).T.astype(float)          # strains x layers
    total = out_deg.sum(axis=1)
    safe = np.where(total == 0, 1.0, total)
    participation = 1.0 - np.sum((out_deg / safe[:, None]) ** 2, axis=1)

    X = StandardScaler().fit_transform(out_deg)
    pca = PCA(random_state=SEED).fit(X)
    Z = pca.transform(X)
    evr = pca.explained_variance_ratio_

    var = pd.read_csv(TAB / "ordination_variance.csv")
    gmm = pd.read_csv(TAB / "mixture_model_selection.csv")
    gmm = gmm[gmm["embedding"].str.startswith("12D")]

    fig = plt.figure(figsize=(14.0, 8.8), )
    gs = fig.add_gridspec(2, 6, hspace=0.42, wspace=1.5,
                          left=0.06, right=0.97, top=0.93, bottom=0.08)
    axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]),
            fig.add_subplot(gs[0, 4:6]), fig.add_subplot(gs[1, 0:3]),
            fig.add_subplot(gs[1, 3:6])]

    # a - the ordination the dimensionality claim is about
    ax = axes[0]
    sc = ax.scatter(Z[:, 0], Z[:, 1], c=total, cmap="viridis", s=58,
                    edgecolor="white", linewidth=0.7, zorder=3)
    cbar(fig, ax, sc, "Total outgoing interactions")
    ax.axhline(0, color=GRID, lw=0.8, zorder=1); ax.axvline(0, color=GRID, lw=0.8, zorder=1)
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    ax.set_title("Strains fill the space continuously", color=INK)
    panel_tag(ax, "a")

    # b - what the retained axes are made of
    ax = axes[1]
    L = pca.components_[:3].T                       # layers x 3
    vmax = float(np.abs(L).max())
    im = ax.imshow(L, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["PC1", "PC2", "PC3"])
    ax.set_yticks(range(len(LAYER_ORDER))); ax.set_yticklabels(LAYER_ORDER, fontsize=8)
    for i in range(L.shape[0]):
        for j in range(3):
            if abs(L[i, j]) >= 0.30:
                ax.text(j, i, f"{L[i, j]:.2f}", ha="center", va="center", fontsize=7.5,
                        color="white" if abs(L[i, j]) > 0.45 else INK)
    cbar(fig, ax, im, "Loading")
    ax.set_title("Loadings of the leading axes", color=INK)
    for sp in ax.spines.values():
        sp.set_visible(False)
    panel_tag(ax, "b")

    # c - spectrum against both retention rules
    ax = axes[2]
    k = np.arange(1, 13)
    obs = var["variance_explained"].to_numpy() * 100
    p95 = var["parallel_null_p95"].to_numpy() * 100
    mean = var["parallel_null_mean"].to_numpy() * 100
    lo = var["var_ci_lo"].to_numpy() * 100
    hi = var["var_ci_hi"].to_numpy() * 100
    keep = obs > p95
    colors = [BLUE if kp else (ORANGE if o > m else GREY) for kp, o, m in zip(keep, obs, mean)]
    ax.bar(k, obs, color=colors, edgecolor="white", linewidth=0.7, zorder=2)
    ax.errorbar(k, obs, yerr=[obs - lo, hi - obs], fmt="none", ecolor=INK2,
                elinewidth=0.9, capsize=2.5, zorder=4)
    ax.plot(k, p95, "-o", color=INK, ms=3.5, lw=1.4, zorder=5, label="null 95th pct.")
    ax.plot(k, mean, "--", color=MUTED, lw=1.2, zorder=5, label="null mean")
    ax.set_xticks(k[::2])
    ax.set_xlabel("Principal component (12-layer profile)")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title("Two components clear the null", color=INK)
    ax.legend(loc="upper right")
    ax.annotate("PC3 marginal:\nabove the mean,\nbelow the 95th pct.", (3, obs[2]),
                textcoords="offset points", xytext=(26, 26), fontsize=8, color=ORANGE,
                linespacing=1.4, arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8))
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    panel_tag(ax, "c")

    # d - the specialist-generalist continuum the text cites
    ax = axes[3]
    sc = ax.scatter(total, participation, c=Z[:, 0], cmap="coolwarm", s=58,
                    edgecolor="white", linewidth=0.7, zorder=3)
    cbar(fig, ax, sc, "PC1 score")
    ax.set_xlabel("Total outgoing interactions")
    ax.set_ylabel("Multilayer participation coefficient")
    ax.set_title("Gradual, not partitioned", color=INK)
    ax.text(0.97, 0.05, f"out-degree {int(total.min())}–{int(total.max())}\n"
                        f"participation {participation.min():.2f}–{participation.max():.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=INK2, linespacing=1.4)
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    panel_tag(ax, "d")

    # e - cluster number is not resolved
    ax = axes[4]
    marks = {"full": ("o", BLUE, "-"), "diag": ("s", ORANGE, "--"), "spherical": ("^", MUTED, ":")}
    floor = -20.0
    for cov, (m, c, ls) in marks.items():
        sub = gmm[gmm["covariance"] == cov].sort_values("k")
        y = sub["cv_loglik"].to_numpy()
        yc = np.clip(y, floor, None)
        ax.plot(sub["k"], yc, ls, marker=m, color=c, ms=5, lw=1.4, label=f"{cov} covariance")
        below = y < floor
        if below.any():
            ax.plot(sub["k"].to_numpy()[below], np.full(below.sum(), floor), "v",
                    color=c, ms=6, clip_on=False)
        best = sub.loc[sub["cv_loglik"].idxmax()]
        ax.plot([best["k"]], [np.clip(best["cv_loglik"], floor, None)], marker=m,
                color=c, ms=12, mfc="none", mew=1.8)
    ax.set_ylim(floor, -10.5)
    ax.set_xlabel("Mixture components $k$")
    ax.set_ylabel("Held-out log-likelihood (10-fold CV)")
    ax.set_title("Selected $k$ depends on the model", color=INK)
    ax.legend(loc="lower left", ncol=1)
    ax.text(0.97, 0.97, "CV-optimal $k$:  full 1,  diag 3,  spherical 2\n"
                        "dip test on PC1 (12-layer): $P$ = 0.97\n"
                        "▼ below axis; full covariance reaches −3673 at $k$ = 8",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            color=INK2, linespacing=1.5)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    panel_tag(ax, "e")

    for ext in FIGURE_FORMATS:
        fig.savefig(OUT / f"figure_strategy_space.{ext}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote figure_strategy_space  (PC1 {evr[0]*100:.1f}%, PC2 {evr[1]*100:.1f}%, "
          f"retained by 95th pct: {int(keep.sum())})")


if __name__ == "__main__":
    main()
