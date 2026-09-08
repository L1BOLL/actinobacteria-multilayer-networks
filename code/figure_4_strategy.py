#!/usr/bin/env python3
"""Final Figure 4 — the ecological strategy space is continuous, not discrete.

Replaces the pre_final Figure 4, which coloured strains as "Cluster 0-5" while the
manuscript's conclusion was that no discrete clusters are supported. The repo's own
generate_figs.make_figure_4 has the same problem with four clusters. Here strains are
coloured by a continuous quantity (total out-degree) and the two panels that carry the
discreteness question are shown alongside the ordination, so the figure states the
same thing the text does.

All inputs are GEN outputs of BN_paper_methods/run_all.py:
    data/ecological_embeddings_4D.csv                    4-D category profile
    data/ecological_identity_with_PCs_and_clusters.csv   out-degrees, participation
    supplementary/tables/p2_3_pca12_variance.csv         12-D spectrum + Horn null
    supplementary/tables/p2_2_gmm_modelselect.csv        BIC / ICL / CV by k
    supplementary/tables/p2_2_diptest.csv                dip test
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_io import FIGURE_FORMATS
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
DATA, TABLES = REPO / "data", REPO / "supplementary" / "tables"
OUT = Path(__file__).resolve().parents[1] / "supplementary" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, INK, INK2, MUTED, NULL_BAND, GRID = (
    "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86", "#d9d9d6", "#ebebe8")
CATS = ["DA", "IA", "MM", "MC"]


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans", "font.size": 8.5,
        "axes.titlesize": 9.5, "axes.labelsize": 8.5, "axes.titlelocation": "left", "axes.titlepad": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": INK2, "axes.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 7.5,
    })


def tag(ax, letter: str, dx: float = -0.17) -> None:
    ax.text(dx, 1.05, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="left", color=INK)


def main() -> None:
    emb = pd.read_csv(DATA / "ecological_embeddings_4D.csv", index_col=0)
    ident = pd.read_csv(DATA / "ecological_identity_with_PCs_and_clusters.csv").set_index("strain")
    var12 = pd.read_csv(TABLES / "p2_3_pca12_variance.csv")
    gmm = pd.read_csv(TABLES / "p2_2_gmm_modelselect.csv")
    dip = pd.read_csv(TABLES / "p2_2_diptest.csv")
    style()

    # 4-D PCA on z-standardized profiles, as the Methods specify. Running it on the
    # raw embedding gives 54.2 / 23.3, which is where the old figure's axis labels
    # came from and why they disagreed with the text.
    X = StandardScaler().fit_transform(emb[CATS].values)
    pca = PCA(n_components=2).fit(X)
    S = pca.transform(X)
    pc1, pc2 = pca.explained_variance_ratio_[:2] * 100

    outdeg = ident.loc[emb.index, [c for c in ident.columns if c.startswith("outdeg_")]].sum(axis=1)

    fig = plt.figure(figsize=(12.6, 4.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.30, 0.92, 0.92], wspace=0.36)

    # (a) ordination, coloured continuously ------------------------------- #
    ax = fig.add_subplot(gs[0, 0])
    sc = ax.scatter(S[:, 0], S[:, 1], c=outdeg.values, cmap="viridis", s=46,
                    edgecolor="white", linewidth=0.6, zorder=3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("Total outgoing interactions", fontsize=8)
    cb.outline.set_linewidth(0.6)
    # Scale loadings to the point cloud, then widen the axes to contain the
    # labelled arrow tips - otherwise a long loading (IA) lands outside the axes.
    scale = 0.78 * np.abs(S).max()
    tips = pca.components_[:2].T * scale
    xs = np.concatenate([S[:, 0], tips[:, 0] * 1.18])
    ys = np.concatenate([S[:, 1], tips[:, 1] * 1.18])
    padx, pady = 0.10 * np.ptp(xs), 0.10 * np.ptp(ys)
    for k, cat in enumerate(CATS):
        vx, vy = tips[k, 0], tips[k, 1]
        ax.annotate("", xy=(vx, vy), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, linewidth=1.3))
        ax.text(vx * 1.10, vy * 1.10, cat, color=ORANGE, fontsize=9, fontweight="bold",
                ha="center", va="center")
    ax.axhline(0, color=GRID, linewidth=0.8, zorder=0)
    ax.axvline(0, color=GRID, linewidth=0.8, zorder=0)
    ax.set_xlim(xs.min() - padx, xs.max() + padx)
    ax.set_ylim(ys.min() - pady, ys.max() + pady)
    ax.set_xlabel(f"PC1 ({pc1:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({pc2:.1f}%)")
    ax.set_title("Strains fill the space continuously")
    tag(ax, "a", dx=-0.20)

    # (b) 12-D spectrum vs Horn null --------------------------------------- #
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(1, len(var12) + 1)
    keep = var12["significant_vs_parallel"].values.astype(bool)
    ax.bar(x, var12["variance_explained"] * 100,
           color=[BLUE if k else NULL_BAND for k in keep], zorder=3, width=0.72)
    ax.plot(x, var12["parallel_null_mean"] * 100, "o-", color=ORANGE, markersize=3.6,
            linewidth=1.3, label="Horn parallel-analysis null", zorder=4)
    n_keep = int(keep.sum())
    ax.set_xticks(x[::2])
    ax.set_xlabel("Principal component (12-D profile)")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title(f"{n_keep} components exceed the null")
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.legend(loc="upper right")
    tag(ax, "b")

    # (c) held-out likelihood by k, for every covariance parameterization --- #
    # Showing only full covariance would imply k = 1 is the settled answer. It is
    # the answer under that parameterization; the constrained ones prefer 2-3, and
    # the figure has to say so or it overstates the paper's own conclusion.
    ax = fig.add_subplot(gs[0, 2])
    styles = {"full": (BLUE, "-", "o"), "diag": (ORANGE, "--", "s"), "spherical": (MUTED, ":", "^")}
    picks = []
    for cov, (col, ls, mk) in styles.items():
        sub = gmm[(gmm["embedding"] == "12D_outdegree_profile") & (gmm["covariance"] == cov)].sort_values("k")
        if sub.empty:
            continue
        ax.plot(sub["k"], sub["cv_loglik"], marker=mk, linestyle=ls, color=col,
                markersize=3.6, linewidth=1.4, label=f"{cov} covariance")
        k_best = int(sub.loc[sub["cv_loglik"].idxmax(), "k"])
        ax.plot(k_best, sub["cv_loglik"].max(), marker=mk, color=col, markersize=9,
                markerfacecolor="white", markeredgewidth=1.6, zorder=5)
        picks.append(f"{cov} {k_best}")
    # Full covariance collapses to -3673 by k = 8, which would flatten the ~6 units
    # that separate the three optima. Clip to the region containing every optimum
    # and say in the axis what falls outside it.
    d12 = gmm[gmm["embedding"] == "12D_outdegree_profile"]["cv_loglik"]
    best = d12.max()
    lo, hi = best - 7.0, best + 1.2
    off = d12[d12 < lo]
    ax.set_ylim(lo, hi)
    for cov, (col, ls, mk) in styles.items():
        sub = gmm[(gmm["embedding"] == "12D_outdegree_profile") & (gmm["covariance"] == cov)].sort_values("k")
        below = sub[sub["cv_loglik"] < lo]
        ax.plot(below["k"], [lo] * len(below), marker="v", linestyle="none",
                color=col, markersize=4.5, clip_on=False, zorder=6)
    ax.set_xlabel("Mixture components $k$")
    ax.set_ylabel("Held-out log-likelihood (10-fold CV)")
    ax.set_title("Selected $k$ depends on the covariance model")
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
    ax.legend(loc="lower left")
    if len(off):
        ax.text(0.5, -0.22, f"▾ below axis ({len(off)} points; full covariance reaches {off.min():.0f} at $k$ = 8)",
                transform=ax.transAxes, fontsize=7, color=MUTED, ha="center")
    dip12 = dip[(dip["embedding"] == "12D_outdegree_profile") & (dip["component"] == "PC1")]["pvalue"].iloc[0]
    ax.text(0.97, 0.94, "CV-optimal $k$:  " + ",  ".join(picks) + f"\ndip test on PC1: $P$ = {dip12:.2f}",
            transform=ax.transAxes, fontsize=7.3, color=INK2, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="none", alpha=0.88))
    tag(ax, "c")

    for ext in FIGURE_FORMATS:
        fig.savefig(OUT / f"figure4_strategy_space.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote figure4_strategy_space.{{png,svg}}")
    print(f"  4-D PCA (z-standardized): PC1 {pc1:.1f}%  PC2 {pc2:.1f}%  cumulative {pc1 + pc2:.1f}%")
    print(f"  Horn retains {n_keep} components; CV-optimal k by covariance: {', '.join(picks)}; dip P={dip12:.2f}")


if __name__ == "__main__":
    main()
