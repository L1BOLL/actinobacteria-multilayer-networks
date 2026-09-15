#!/usr/bin/env python3
"""Spatial dependence of IG and CS: close-proximity versus separated inocula.

Reads the four spreadsheets in `data/matrices_spatial/`, the paired experiment in
which growth inhibition (IG) and change in sporulation (CS) were scored with the
two inocula close enough for the colonies to meet and, on the same plate, with the
inocula separated. Nothing from `data/matrices_xls/` is used: the paired
close-proximity matrices are a separate scoring of the same phenotypes and are not
the direct-contact layers of Table 1. The disagreement between the two is
reported, not hidden.

The confrontation unit is the unordered strain dyad, so every resampling scheme
here moves both directed entries of a dyad together:

  * paired randomisation - the configuration label is exchanged within a dyad,
    giving the null for the density difference (Figure 6d);
  * dyad bootstrap - the interval on the density difference;
  * McNemar and the exact binomial on discordant directed pairs are reported
    beside the randomisation test, labelled, because they answer related but
    different questions and give materially different probabilities.

Figure 6 and Figure S12.
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
from scipy import stats

from data_io import (
    FIG_DIR,
    REPORT_PATH,
    SEED,
    TABLE_DIR,
    ensure_output_dirs,
    load_layer_df,
    load_layers,
    save_figure,
)
from report_utils import dataframe_to_md, replace_section
from source_paths import PACKAGE_ROOT

SPATIAL_DIR = PACKAGE_ROOT / "data" / "matrices_spatial"
FILES = {
    ("IG", "close"): "IG Short range (plate's bottom).xlsx",
    ("IG", "separated"): "IG Long range (plate's middle).xlsx",
    ("CS", "close"): "CS Short range (plate's bottom).xlsx",
    ("CS", "separated"): "CS Long range (plate's bottom).xlsx",
}
N_PERM = 10_000
N_BOOT = 10_000

BLUE, ORANGE, INK, INK2, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86"
COLOURS = {"IG": ORANGE, "CS": BLUE}


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


def panel_tag(ax, letter: str) -> None:
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left", color=INK)


def load_pair(layer: str, order: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Both configurations of one layer, aligned to the main-matrix strain order."""
    mats = []
    for cond in ("close", "separated"):
        df = load_layer_df(SPATIAL_DIR / FILES[(layer, cond)])
        if set(df.index) != set(order) or set(df.columns) != set(order):
            raise SystemExit(f"{layer}/{cond}: strain set differs from the main matrices")
        mats.append((df.reindex(index=order, columns=order).to_numpy() > 0).astype(np.int8))
    return mats[0], mats[1]


def analyse(layer: str, close: np.ndarray, sep: np.ndarray, main: np.ndarray,
            rng: np.random.Generator) -> dict:
    n = close.shape[0]
    npairs = n * (n - 1)
    iu, ju = np.triu_indices(n, 1)                  # one index per unordered dyad
    # directed entries of each dyad, both directions side by side: shape (dyads, 2)
    c = np.stack([close[iu, ju], close[ju, iu]], axis=1)
    s = np.stack([sep[iu, ju], sep[ju, iu]], axis=1)

    persistent = int(np.sum((c == 1) & (s == 1)))
    lost = int(np.sum((c == 1) & (s == 0)))
    gained = int(np.sum((c == 0) & (s == 1)))

    chi2 = (lost - gained) ** 2 / (lost + gained)
    p_mcnemar = float(stats.chi2.sf(chi2, 1))
    p_binom = float(stats.binomtest(lost, lost + gained, 0.5).pvalue)

    d_obs = (c.sum() - s.sum()) / npairs

    # Paired randomisation: exchange the configuration label of a whole dyad.
    swap = rng.random((N_PERM, len(iu))) < 0.5
    block = (c.sum(axis=1) - s.sum(axis=1)).astype(float)   # per-dyad contribution
    null = ((np.where(swap, -block, block)).sum(axis=1)) / npairs
    p_perm = float((np.sum(np.abs(null) >= abs(d_obs) - 1e-12) + 1) / (N_PERM + 1))

    # Dyad bootstrap for the interval on the density difference.
    idx = rng.integers(0, len(iu), size=(N_BOOT, len(iu)))
    boot = block[idx].sum(axis=1) / npairs
    ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))

    log_or = np.log(lost / gained)
    se = np.sqrt(1 / lost + 1 / gained)
    off = ~np.eye(n, dtype=bool)
    return dict(
        layer=layer, n_pairs=npairs,
        edges_close=int(c.sum()), edges_separated=int(s.sum()),
        density_close=c.sum() / npairs, density_separated=s.sum() / npairs,
        density_difference=d_obs, ci_low=ci[0], ci_high=ci[1],
        persistent=persistent, lost=lost, gained=gained,
        absent_in_both=int(npairs - persistent - lost - gained),
        odds_loss_vs_gain=lost / gained,
        odds_ci_low=float(np.exp(log_or - 1.96 * se)), odds_ci_high=float(np.exp(log_or + 1.96 * se)),
        mcnemar_chi2=float(chi2), p_mcnemar=p_mcnemar, p_exact_binomial=p_binom,
        p_paired_dyad_permutation=p_perm, n_permutations=N_PERM, n_bootstrap=N_BOOT,
        entries_differing_from_main_layer=int(np.sum(close[off] != main[off])),
        seed=SEED, null=null,
    )


def fmt_p(p: float) -> str:
    if p <= 1.0 / (N_PERM + 1):
        return r"$P < 10^{-4}$"
    if p < 1e-3:
        e = int(np.floor(np.log10(p)))
        return rf"$P = {p / 10 ** e:.1f}\times10^{{{e}}}$"
    return rf"$P = {p:.3f}$"


def figure_6(res: dict[str, dict]) -> None:
    style()
    fig, axes = plt.subplots(1, 4, figsize=(15.4, 3.9), constrained_layout=True)

    ax = axes[0]
    for lab, r in res.items():
        ax.plot([0, 1], [r["density_close"], r["density_separated"]], "-o", color=COLOURS[lab],
                lw=2.0, ms=7, label=lab, zorder=3)
        ax.annotate(f"{r['edges_close']}", (0, r["density_close"]), textcoords="offset points",
                    xytext=(-6, 9), ha="right", color=COLOURS[lab], fontsize=9)
        ax.annotate(f"{r['edges_separated']}", (1, r["density_separated"]), textcoords="offset points",
                    xytext=(6, 9), ha="left", color=COLOURS[lab], fontsize=9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Close\nproximity", "Spatially\nseparated"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0, 1.2 * max(r["density_close"] for r in res.values()))
    ax.set_ylabel("Directed network density")
    ax.set_title("Density collapses for IG, not CS")
    ax.text(0, -0.28, "denominator 3,540 ordered pairs", transform=ax.transAxes, color=MUTED, fontsize=8)
    ax.legend(loc="center right"); panel_tag(ax, "a")

    ax = axes[1]
    for k, (lab, r) in enumerate(res.items()):
        y = 1 - k
        ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=COLOURS[lab], lw=3, solid_capstyle="round")
        ax.plot(r["density_difference"], y, "o", color=COLOURS[lab], ms=8)
        ax.annotate(rf"$\Delta\rho$ = {r['density_difference']:.4f}  ({r['ci_low']:.4f}–{r['ci_high']:.4f})",
                    (r["density_difference"], y), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=8.5, color=INK2)
    ax.axvline(0, color=INK2, lw=0.8, ls="--")
    ax.set_yticks([1, 0]); ax.set_yticklabels(list(res)); ax.set_ylim(-0.7, 1.7)
    ax.set_xlabel(r"$\Delta$ density (close − separated)")
    ax.set_title("With 95% dyad-bootstrap interval"); panel_tag(ax, "b")

    ax = axes[2]
    x = np.arange(3); w = 0.36
    for k, (lab, r) in enumerate(res.items()):
        vals = [r["persistent"], r["lost"], r["gained"]]
        bars = ax.bar(x + (k - 0.5) * w, vals, w, color=COLOURS[lab], label=lab)
        for b, v in zip(bars, vals):
            ax.annotate(f"{v}", (b.get_x() + b.get_width() / 2, v), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=8.5, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(["Persistent", "Lost after\nseparation", "Gained after\nseparation"])
    ax.set_ylabel("Ordered strain pairs")
    txt = "\n".join(
        rf"{lab}: OR$_{{loss:gain}}$ = {r['odds_loss_vs_gain']:.2f} ({r['odds_ci_low']:.2f}–{r['odds_ci_high']:.2f})"
        for lab, r in res.items())
    ax.text(0.98, 0.72, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK2)
    ax.set_title("IG loses edges; CS exchanges them"); ax.legend(loc="upper right"); panel_tag(ax, "c")

    ax = axes[3]
    for lab, r in res.items():
        ax.hist(r["null"], bins=60, color=COLOURS[lab], alpha=0.35, lw=0)
        ax.axvline(r["density_difference"], color=COLOURS[lab], lw=2)
    for lab, r in res.items():
        xtxt = r["density_difference"]
        ha = "right" if xtxt > 0.05 else "left"
        dx = -6 if ha == "right" else 8
        ax.annotate(f"{lab} observed\n{fmt_p(r['p_paired_dyad_permutation'])}", (xtxt, 0.72),
                    xycoords=("data", "axes fraction"), textcoords="offset points", xytext=(dx, 0),
                    ha=ha, va="center", fontsize=9, color=COLOURS[lab])
    ax.set_xlabel(r"$\Delta$ density under paired dyad randomisation")
    ax.set_ylabel(f"Draws (of {N_PERM:,})")
    ax.set_title("Only IG departs from its null"); panel_tag(ax, "d")

    save_figure(fig, "figure_spatial_dependence")
    plt.close(fig)


def figure_s12(mats: dict[tuple[str, str], np.ndarray], res: dict[str, dict]) -> None:
    """The four matrices, drawn as stored: rows = receivers, columns = senders (as Figure S1)."""
    style()
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 10.6), constrained_layout=True)
    titles = {"close": "close proximity", "separated": "spatially separated"}
    letters = iter("abcd")
    for i, layer in enumerate(("IG", "CS")):
        for j, cond in enumerate(("close", "separated")):
            ax = axes[i, j]
            A = mats[(layer, cond)]
            n = A.shape[0]
            ax.imshow(1 - A, cmap="gray", vmin=0, vmax=1, interpolation="nearest",
                      extent=(0.5, n + 0.5, n + 0.5, 0.5))
            edges = res[layer]["edges_close" if cond == "close" else "edges_separated"]
            ax.set_title(f"{layer}, {titles[cond]}  ({edges} edges)")
            ax.set_xlabel("Sender strain"); ax.set_ylabel("Receiver strain")
            ax.set_xticks([1, 20, 40, 60]); ax.set_yticks([1, 20, 40, 60])
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(True)
            panel_tag(ax, next(letters))
    save_figure(fig, "figure_spatial_matrices")
    plt.close(fig)


def main() -> None:
    ensure_output_dirs()
    layers = load_layers()
    order = layers["IG"].index.tolist()
    rng = np.random.default_rng(SEED)

    mats: dict[tuple[str, str], np.ndarray] = {}
    res: dict[str, dict] = {}
    for layer in ("IG", "CS"):
        close, sep = load_pair(layer, order)
        mats[(layer, "close")], mats[(layer, "separated")] = close, sep
        main_layer = (layers[layer].to_numpy() > 0).astype(np.int8)
        res[layer] = analyse(layer, close, sep, main_layer, rng)

    table = pd.DataFrame([{k: v for k, v in r.items() if k != "null"} for r in res.values()])
    table.to_csv(TABLE_DIR / "spatial_configuration.csv", index=False)

    figure_6(res)
    figure_s12(mats, res)

    show = table.drop(columns=["seed", "n_permutations", "n_bootstrap"]).T
    show.columns = table["layer"]
    body = f"""
Seed: {SEED}. Paired close/separated matrices from `data/matrices_spatial/`, aligned to the
main-matrix strain order. Denominator 3,540 ordered non-self pairs. The randomisation test
exchanges the configuration label of a whole unordered dyad ({N_PERM:,} draws, (b+1)/(m+1));
the interval is a {N_BOOT:,}-draw dyad bootstrap; McNemar and the exact binomial act on the
discordant directed pairs and are reported for comparison.

{dataframe_to_md(show.reset_index().rename(columns={'index': 'quantity'}), index=False)}

`entries_differing_from_main_layer` counts the directed entries in which the close-proximity
matrix disagrees with the direct-contact layer of the same phenotype in `data/matrices_xls/`:
the paired spatial dataset is a separate scoring and is never merged with the twelve layers.
"""
    replace_section(REPORT_PATH, "<!-- spatial_configuration RESULTS -->",
                    "Spatial configuration: close proximity versus separation", body)
    print(table.drop(columns=["seed"]).T.to_string())


if __name__ == "__main__":
    main()
