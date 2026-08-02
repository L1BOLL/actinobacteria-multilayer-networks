#!/usr/bin/env python3
"""Publication figures for the annotation-data results (h3 and metabolism).

Produces three figures, each as PNG (300 dpi) and PDF (vector, for typesetting):

  fig_h3_phylogeny            main text  - ecological strategy is decoupled from 16S
  fig_metabolism              main text  - assimilation does not predict network position
  figS_permutation_calibration supplement - why the asymptotic lambda test fails here

Design notes, because these figures make a *negative* claim and are easy to get
wrong: a bare bar chart of Pagel's lambda reads as "strong signal" no matter what
the p-values say, since the point estimates run up to 0.95. Every panel therefore
plots the observed statistic **against the null it was tested on**, so the reader
sees the comparison the inference actually made rather than the raw estimate.
Non-significant effects are never coloured as though they were significant.
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
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

from data_io import CATEGORY_MAP, FIG_DIR, LAYER_ORDER, SEED, canonical_strain, ensure_phase2_dirs, load_tensor
from p2_8_phylogeny_pgls import (
    N_PERM,
    build_traits,
    load_distance_matrix,
    pagel_lambda,
    vcv_from_linkage,
    vcv_from_tree,
    TREE_PATH,
)
from p2_12_metabolism import load_profiles

# --- palette (validated categorical slots 1-2 + neutral ink) ------------------ #
BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a86"
NULL_BAND = "#d9d9d6"
GRID = "#ebebe8"

PRETTY = {
    **{f"outdeg_{l}": f"out-degree {l}" for l in LAYER_ORDER},
    **{f"cat_{c}": f"{c} category mean" for c in CATEGORY_MAP},
    "total_out_degree": "total out-degree",
    "participation": "participation coefficient",
    "PC1_12D": "PC1 (12-D profile)",
    "PC2_12D": "PC2 (12-D profile)",
    "PC3_12D": "PC3 (12-D profile)",
}


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "axes.titlelocation": "left",
        "axes.titlepad": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.edgecolor": INK2,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.labelcolor": INK,
        "legend.frameon": False,
        "legend.fontsize": 7.5,
    })


def panel_tag(ax, letter: str) -> None:
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"{stem}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {stem}.png / .pdf")


# --------------------------------------------------------------------------- #
# shared data
# --------------------------------------------------------------------------- #
def phylo_inputs():
    tensor, _, node_ids = load_tensor()
    traits = build_traits(tensor, node_ids)
    D = load_distance_matrix()
    common = sorted(set(traits.index) & set(D.index))
    Dg = D.loc[common, common].to_numpy()
    traits = traits.loc[common]

    keys = [canonical_strain(n) for n in node_ids]
    idx = [keys.index(s) for s in common]
    profile = tensor.sum(axis=1).T.astype(float)[idx]
    z = (profile - profile.mean(0)) / np.where(profile.std(0) == 0, 1.0, profile.std(0))
    Dp = squareform(pdist(z))
    return traits, Dg, Dp, common


# --------------------------------------------------------------------------- #
# Figure 1 (main) - h3
# --------------------------------------------------------------------------- #
def figure_h3(pgls: pd.DataFrame, mantel: pd.DataFrame, Dg: np.ndarray, Dp: np.ndarray) -> None:
    ml = pgls[pgls["topology"].str.startswith("ML")].copy()
    ml["label"] = ml["trait"].map(lambda t: PRETTY.get(t, t))
    ml = ml.sort_values("pagel_lambda").reset_index(drop=True)

    fig = plt.figure(figsize=(11.6, 4.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.85], wspace=0.42)

    # --- A: observed lambda inside its own permutation null -------------------
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(ml))
    ax.barh(y, ml["null_lambda_p95"], height=0.72, color=NULL_BAND, zorder=1)
    ax.scatter(ml["pagel_lambda"], y, s=26, color=BLUE, zorder=3,
               edgecolor="white", linewidth=0.8)
    ax.axvspan(0.31, 0.48, color=ORANGE, alpha=0.13, zorder=0)
    ax.axvline(0.31, color=ORANGE, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax.axvline(0.48, color=ORANGE, lw=0.9, ls=(0, (4, 2)), zorder=2)
    # Null bars fill the whole plotting area, so annotation and legend go outside it.
    ax.text(0.395, len(ml) - 0.35, "previously reported\nλ = 0.31–0.48",
            ha="center", va="bottom", fontsize=7, color=ORANGE, linespacing=1.3)
    ax.set_yticks(y)
    ax.set_yticklabels(ml["label"])
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.7, len(ml) + 0.9)
    ax.set_xlabel("Pagel's $\\lambda$")
    ax.set_title("Every trait's $\\lambda$ falls inside its own null", color=INK)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="none", color=BLUE, markersize=5.5,
               markeredgecolor="white", label="observed $\\lambda$"),
        matplotlib.patches.Patch(color=NULL_BAND, label="95th pct. of tip-label null"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, handlelength=1.3)
    panel_tag(ax, "a")

    # --- B: tree-free Mantel --------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    iu = np.triu_indices(len(Dg), 1)
    x, yv = Dg[iu], Dp[iu]
    ax.scatter(x, yv, s=5, color=BLUE, alpha=0.16, edgecolor="none", zorder=2)
    b, a = np.polyfit(x, yv, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, a + b * xs, color=INK, lw=1.6, zorder=3)
    rho = float(mantel.loc[mantel["comparison"].str.startswith("12-D"), "mantel_rho"].iloc[0])
    p = float(mantel.loc[mantel["comparison"].str.startswith("12-D"), "p_value"].iloc[0])
    ax.set_xlabel("16S genetic distance")
    ax.set_ylabel("Interaction-profile distance (12-D, z-scored)")
    ax.set_title("No relationship without a tree", color=INK)
    ax.text(0.97, 0.05, f"Mantel $\\rho$ = {rho:+.3f}\n$p$ = {p:.2f}   ({len(x):,} pairs)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.8, color=INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "b")

    # --- C: 16S resolution limit ---------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    n_geno, memb = connected_components(Dg == 0, directed=False)
    sizes = np.sort(pd.Series(memb).value_counts().to_numpy())[::-1]
    colors = [BLUE if s > 1 else MUTED for s in sizes]
    ax.bar(np.arange(len(sizes)), sizes, color=colors, width=0.86)
    ax.set_xlabel(f"16S genotype (n = {n_geno})")
    ax.set_ylabel("Strains sharing the sequence")
    ax.set_title("16S cannot resolve the collection", color=INK)
    tied = int(sizes[sizes > 1].sum())
    ax.text(0.97, 0.93,
            f"{len(Dg)} strains → {n_geno} genotypes\n"
            f"{tied} strains ({tied / len(Dg):.0%}) are\nindistinguishable",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.8, color=INK2,
            linespacing=1.4)
    ax.set_xticks([])
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "c")

    save(fig, "fig_h3_phylogeny")


# --------------------------------------------------------------------------- #
# Figure 2 (main) - metabolism
# --------------------------------------------------------------------------- #
def figure_metabolism() -> None:
    strain = pd.read_csv(FIG_DIR.parent / "tables" / "p2_12_metabolism_strain.csv")
    dyad = pd.read_csv(FIG_DIR.parent / "tables" / "p2_12_metabolism_dyad.csv")
    prof = pd.read_csv(FIG_DIR.parent / "tables" / "p2_12_metabolism_profiles.csv")

    fig = plt.figure(figsize=(11.0, 4.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.05], wspace=0.62)

    # --- A: capacity vs activity ---------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(prof["assimilation_capacity"], prof["total_out_degree"],
               s=26, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.6)
    r = spearmanr(prof["assimilation_capacity"], prof["total_out_degree"])
    ax.set_xlabel("Assimilation capacity (Σ 0–4 scores, 8 sugars)")
    ax.set_ylabel("Total out-degree")
    ax.set_title("Metabolic capacity does not predict activity", color=INK)
    ax.text(0.03, 0.05, f"Spearman $\\rho$ = {r.statistic:+.3f}\n$p$ = {r.pvalue:.2f}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, color=INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "a")

    # --- B: every strain-level association, ranked ---------------------------
    ax = fig.add_subplot(gs[0, 1])
    s = strain.sort_values("spearman_rho").reset_index(drop=True)
    ax.barh(np.arange(len(s)), s["spearman_rho"], height=0.74, color=BLUE, alpha=0.85)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(np.arange(len(s)))
    ax.set_yticklabels(s["response"], fontsize=6.6)
    ax.set_xlabel("Spearman $\\rho$ with assimilation capacity")
    ax.set_title("None of 18 associations survives FDR", color=INK)
    ax.text(0.97, 0.04, f"all $q_{{BH}}$ > {s['q_BH'].min():.2f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.8, color=INK2)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "b")

    # --- C: dyad-level Mantel, honest about significance ---------------------
    ax = fig.add_subplot(gs[0, 2])
    d = dyad.sort_values("mantel_rho").reset_index(drop=True)
    ax.barh(np.arange(len(d)), d["mantel_rho"], height=0.74, color=BLUE, alpha=0.85)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["layer"], fontsize=7)
    ax.set_xlabel("Mantel $\\rho$: niche distance vs interaction")
    ax.set_title("Dissimilar pairs antagonise slightly more", color=INK)
    # Name the two leading layers without implying significance.
    lead = d.nlargest(2, "mantel_rho")
    for _, row in lead.iterrows():
        i = int(d.index[d["layer"] == row["layer"]][0])
        ax.annotate(f"  {row['layer']}  $q$ = {row['q_BH']:.2f}",
                    xy=(row["mantel_rho"], i), fontsize=7, color=INK2,
                    va="center", ha="left")
    ax.text(0.97, 0.05, "no layer reaches $q_{BH}$ < 0.05", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.8, color=INK2)
    ax.set_xlim(d["mantel_rho"].min() * 1.35, d["mantel_rho"].max() * 2.3)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "c")

    save(fig, "fig_metabolism")


# --------------------------------------------------------------------------- #
# Supplementary figure - why the asymptotic test fails
# --------------------------------------------------------------------------- #
def figure_calibration(pgls: pd.DataFrame, traits: pd.DataFrame, Dg: np.ndarray) -> None:
    rng = np.random.default_rng(SEED)
    C, order = vcv_from_tree(TREE_PATH, list(traits.index))

    # Traits the chi-squared test would have called significant.
    ml = pgls[pgls["topology"].str.startswith("ML")]
    shown = ml.nsmallest(4, "p_chi2_ref")["trait"].tolist()

    fig = plt.figure(figsize=(11.0, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.42)

    # --- A: null distribution of lambda under tip shuffling ------------------
    ax = fig.add_subplot(gs[0, 0])
    for i, trait in enumerate(shown):
        y = traits.loc[order, trait].to_numpy(float)
        null = np.array([pagel_lambda(rng.permutation(y), C)[0] for _ in range(N_PERM)])
        obs = float(ml.loc[ml["trait"] == trait, "pagel_lambda"].iloc[0])
        pp = float(ml.loc[ml["trait"] == trait, "p_perm"].iloc[0])
        parts = ax.violinplot(null, positions=[i], vert=False, widths=0.82,
                              showextrema=False, showmedians=False)
        for b in parts["bodies"]:
            b.set_facecolor(NULL_BAND); b.set_alpha(1.0); b.set_edgecolor("none")
        ax.scatter([obs], [i], s=34, color=BLUE, zorder=3, edgecolor="white", linewidth=0.9)
        ax.text(0.015, i, f"$p_{{perm}}$ = {pp:.2f}", fontsize=7.5, color=INK2,
                va="center", ha="left", transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([PRETTY.get(t, t) for t in shown])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Pagel's $\\lambda$")
    ax.set_title("Shuffling trait values reproduces the observed $\\lambda$", color=INK)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="none", color=BLUE, markersize=5.5,
               markeredgecolor="white", label="observed"),
        matplotlib.patches.Patch(color=NULL_BAND, label=f"{N_PERM} tip-label shuffles"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, handlelength=1.3)
    ax.set_ylim(-0.75, len(shown) - 0.25)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "a")

    # --- B: chi-squared p vs permutation p -----------------------------------
    ax = fig.add_subplot(gs[0, 1])
    lo = 1e-6
    xv = pgls["p_chi2_ref"].clip(lower=lo)
    yv = pgls["p_perm"].clip(lower=lo)
    ax.scatter(xv, yv, s=22, color=BLUE, alpha=0.7, edgecolor="white", linewidth=0.5, zorder=3)
    ax.plot([lo, 1], [lo, 1], color=INK, lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.axvline(0.05, color=ORANGE, lw=0.9, zorder=1)
    ax.axhline(0.05, color=MUTED, lw=0.9, zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, 1.6); ax.set_ylim(lo, 1.6)
    ax.set_xlabel("$p$ from asymptotic $\\chi^2$ test")
    ax.set_ylabel("$p$ from tip-label permutation")
    ax.set_title("The asymptotic test is anticonservative", color=INK)
    n_bad = int(((pgls["p_chi2_ref"] < 0.05) & (pgls["p_perm"] >= 0.05)).sum())
    ax.text(0.04, 0.06,
            f"{n_bad} of {len(pgls)} tests are called\nsignificant by $\\chi^2$ and\nnull by permutation",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, color=INK2,
            linespacing=1.4)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    panel_tag(ax, "b")

    save(fig, "figS_permutation_calibration")


def main() -> None:
    ensure_phase2_dirs()
    style()
    tables = FIG_DIR.parent / "tables"
    pgls = pd.read_csv(tables / "p2_8_pgls.csv")
    mantel = pd.read_csv(tables / "p2_8_mantel.csv")
    traits, Dg, Dp, _ = phylo_inputs()

    figure_h3(pgls, mantel, Dg, Dp)
    figure_metabolism()
    figure_calibration(pgls, traits, Dg)


if __name__ == "__main__":
    main()
