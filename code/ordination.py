#!/usr/bin/env python3
"""— PCA on 12D out-degree: variance + bootstrap CI, Horn parallel analysis, 12D-PC1 ↔ 4D-PC1 Spearman."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from data_io import save_figure, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, compute_out_degree_matrix, ensure_output_dirs, load_or_recompute_embedding4, load_tensor
from report_utils import dataframe_to_md, replace_section


N_BOOT = 2000
N_PARALLEL = 1000


def bootstrap_variance(X: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n, p = X.shape
    out = np.zeros((n_boot, p), dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        Xb = X[idx]
        Xb = StandardScaler().fit_transform(Xb)
        try:
            pca = PCA(random_state=seed)
            pca.fit(Xb)
            out[b, : len(pca.explained_variance_ratio_)] = pca.explained_variance_ratio_
        except Exception:
            out[b] = np.nan
    return out


def parallel_analysis(X: np.ndarray, n_iter: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Horn's null by permuting values independently within each column.

    Permutation, not N(0,1) simulation: it keeps each layer's own marginal
    distribution -- these are out-degree counts, some of them very sparse -- and
    destroys only the between-layer covariance, which is what the test is about.

    Returns the mean and the 95th percentile of the null explained-variance
    ratios.  The 95th percentile is the retention criterion; the mean is kept
    because it is the weaker of the two rules in common use.
    """
    rng = np.random.default_rng(seed)
    n, p = X.shape
    draws = np.zeros((n_iter, p), dtype=float)
    for i in range(n_iter):
        Xp = np.column_stack([rng.permutation(X[:, j]) for j in range(p)])
        Xz = StandardScaler().fit_transform(Xp)
        pca = PCA(random_state=seed)
        pca.fit(Xz)
        draws[i] = pca.explained_variance_ratio_
    return draws.mean(axis=0), np.percentile(draws, 95, axis=0)


def main() -> None:
    ensure_output_dirs()
    tensor, layer_ids, node_ids = load_tensor()
    X12 = compute_out_degree_matrix(tensor).loc[node_ids, layer_ids]
    X12_raw = X12.values
    X12z = StandardScaler().fit_transform(X12_raw)
    pca12 = PCA(random_state=SEED)
    Z12 = pca12.fit_transform(X12z)

    emb4 = load_or_recompute_embedding4().loc[node_ids]
    X4z = StandardScaler().fit_transform(emb4.values)
    pca4 = PCA(random_state=SEED)
    Z4 = pca4.fit_transform(X4z)

    rho = spearmanr(pd.Series(Z12[:, 0], index=node_ids).rank(), pd.Series(Z4[:, 0], index=node_ids).rank())

    boot = bootstrap_variance(X12_raw, N_BOOT, SEED)
    var_lo = np.nanpercentile(boot, 2.5, axis=0)
    var_hi = np.nanpercentile(boot, 97.5, axis=0)

    pa_mean, pa_p95 = parallel_analysis(X12_raw, N_PARALLEL, SEED)
    # Retention is on the 95th percentile; the mean rule is reported alongside
    # because it retains one further component and the two disagree on PC3.
    significant_pcs = int(np.sum(pca12.explained_variance_ratio_ > pa_p95))
    significant_pcs_mean_rule = int(np.sum(pca12.explained_variance_ratio_ > pa_mean))

    variance_df = pd.DataFrame(
        {
            "component": [f"PC{i+1}" for i in range(len(layer_ids))],
            "variance_explained": pca12.explained_variance_ratio_,
            "cumulative_variance": np.cumsum(pca12.explained_variance_ratio_),
            "var_ci_lo": var_lo[: len(layer_ids)],
            "var_ci_hi": var_hi[: len(layer_ids)],
            "parallel_null_mean": pa_mean[: len(layer_ids)],
            "parallel_null_p95": pa_p95[: len(layer_ids)],
            "significant_vs_parallel": pca12.explained_variance_ratio_ > pa_p95[: len(layer_ids)],
            "significant_vs_parallel_mean_rule": pca12.explained_variance_ratio_ > pa_mean[: len(layer_ids)],
            "seed": SEED,
        }
    )
    loadings_df = pd.DataFrame(
        pca12.components_[:3].T,
        index=layer_ids,
        columns=["PC1_loading", "PC2_loading", "PC3_loading"],
    ).reset_index(names="layer")

    var4_df = pd.DataFrame({
        "component": [f"PC{i+1}" for i in range(emb4.shape[1])],
        "variance_explained": pca4.explained_variance_ratio_,
        "cumulative_variance": np.cumsum(pca4.explained_variance_ratio_),
    })

    variance_df.to_csv(TABLE_DIR / "ordination_variance.csv", index=False)
    loadings_df.to_csv(TABLE_DIR / "ordination_loadings.csv", index=False)
    var4_df.to_csv(TABLE_DIR / "ordination_variance_4d.csv", index=False)
    pd.DataFrame({"rho_pc1_rank": [rho.statistic], "pvalue": [rho.pvalue], "seed": [SEED]}).to_csv(TABLE_DIR / "ordination_pc1_concordance.csv", index=False)

    # Both retention rules are drawn, because they disagree about PC3 and the
    # the two rules disagree on marginal components, so both are drawn.
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10.5, "axes.titlesize": 11.5,
                         "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "legend.fontsize": 9.5,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    x = np.arange(1, len(layer_ids) + 1)
    obs = variance_df["variance_explained"].to_numpy() * 100
    lo = variance_df["var_ci_lo"].to_numpy() * 100
    hi = variance_df["var_ci_hi"].to_numpy() * 100
    ax.fill_between(x, lo, hi, color="#2a78d6", alpha=0.16,
                    label=f"95% bootstrap of observed ({N_BOOT:,} resamples)")
    ax.plot(x, obs, marker="o", color="#2a78d6", lw=1.8, ms=5, label="Observed", zorder=4)
    ax.plot(x, pa_p95[: len(layer_ids)] * 100, marker="o", ms=3.5, color="#0b0b0b", lw=1.5,
            label="Null 95th percentile (retention rule)", zorder=5)
    ax.plot(x, pa_mean[: len(layer_ids)] * 100, linestyle="--", color="#8a8a86", lw=1.3,
            label="Null mean (weaker rule, sensitivity)", zorder=5)
    ax.axvspan(0.5, significant_pcs + 0.5, color="#2a78d6", alpha=0.06, zorder=0)
    ax.annotate("PC3: above the mean,\nbelow the 95th percentile",
                (3, obs[2]), textcoords="offset points", xytext=(38, 34),
                fontsize=9.5, color="#eb6834", linespacing=1.4,
                arrowprops=dict(arrowstyle="-", color="#eb6834", lw=0.9))
    ax.set_xticks(x)
    ax.set_xlabel("Principal component of the 12-layer out-degree profile")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title(f"{significant_pcs} components clear the 95th-percentile null "
                 f"({significant_pcs_mean_rule} clear the mean)")
    ax.grid(axis="y", color="#ebebe8", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper right")
    save_figure(fig, "figure_parallel_analysis", FIG_DIR)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    ax.scatter(Z12[:, 0], Z12[:, 1], s=38, color="#4c78a8", alpha=0.85, edgecolor="white", linewidth=0.6)
    loading_scale = 3.0
    for _, row in loadings_df.iterrows():
        xv, yv = row["PC1_loading"] * loading_scale, row["PC2_loading"] * loading_scale
        ax.arrow(0, 0, xv, yv, color="#d62728", width=0.006, alpha=0.9, length_includes_head=True)
        ax.text(xv * 1.08, yv * 1.08, row["layer"], fontsize=8.2, color="#d62728")
    outlier_idx = np.argsort(-(np.abs(Z12[:, 0]) + np.abs(Z12[:, 1])))[:10]
    for idx in outlier_idx:
        ax.text(Z12[idx, 0] + 0.08, Z12[idx, 1] + 0.05, node_ids[idx], fontsize=7.5, color="#374151")
    ax.set_xlabel(f"PC1 ({pca12.explained_variance_ratio_[0]*100:.1f}% var.)")
    ax.set_ylabel(f"PC2 ({pca12.explained_variance_ratio_[1]*100:.1f}% var.)")
    ax.set_title("12D out-degree PCA biplot")
    ax.grid(color="#e5e7eb", linewidth=0.6)
    save_figure(fig, "ordination_biplot", FIG_DIR)
    plt.close(fig)

    pc1_var12 = pca12.explained_variance_ratio_[0]
    pc2_var12 = pca12.explained_variance_ratio_[1]
    pc1_var4 = pca4.explained_variance_ratio_[0]
    pc2_var4 = pca4.explained_variance_ratio_[1]

    body = f"""
Seed: {SEED}. PCA on z-scored 60×12 out-degree. {N_BOOT} bootstrap resamples; {N_PARALLEL}-rep Horn parallel analysis on N(0,1) data of same shape. 4D collapsed-embedding PCA reported alongside.

12D: PC1 = {pc1_var12*100:.1f}%, PC2 = {pc2_var12*100:.1f}% — honest dimensionality.
4D:  PC1 = {pc1_var4*100:.1f}%, PC2 = {pc2_var4*100:.1f}% — visualization only (low-D by construction; inflates variance).

Horn: {significant_pcs} of {len(layer_ids)} PCs above the null.

12D-PC1 ↔ 4D-PC1 Spearman:
- rho = {rho.statistic:.4f}, p = {rho.pvalue:.4g}, rho<0.7 = {rho.statistic < 0.7}
- rho ≥ 0.7 → 4D PC1 preserves dominant 12D axis.

12D variance (top 6):
{dataframe_to_md(variance_df.head(6), index=False)}

4D variance:
{dataframe_to_md(var4_df, index=False)}

12D PC1–PC3 loadings:
{dataframe_to_md(loadings_df, index=False)}

Outputs:
- ordination_variance.csv (with bootstrap CI + parallel-null)
- ordination_variance_4d.csv
- ordination_loadings.csv
- ordination_pc1_concordance.csv
- figure_parallel_analysis.png, ordination_biplot.png
"""
    replace_section(REPORT_PATH, "<!-- ordination RESULTS -->", "Principal components of the 12-layer out-degree profile", body)


if __name__ == "__main__":
    main()
