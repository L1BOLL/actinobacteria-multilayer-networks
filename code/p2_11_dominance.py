#!/usr/bin/env python3
"""P2.11 — David's score + Bradley-Terry per layer; cross-layer Spearman + PCA. Direct test of h2."""
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

from data_io import save_figure, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, ensure_phase2_dirs, load_tensor
from report_utils import dataframe_to_md, replace_section


def davids_score(W: np.ndarray) -> np.ndarray:
    """De Vries 2006 NDS. W[i,j]=1 → i beat j."""
    n = W.shape[0]
    interactions = W + W.T
    interactions[interactions == 0] = 1
    P = W / interactions
    w = P.sum(axis=1)
    l = P.sum(axis=0)
    w2 = (P * w[None, :]).sum(axis=1)
    l2 = (P * l[:, None]).sum(axis=0)
    return w + w2 - l - l2


def bradley_terry_strengths(W: np.ndarray, max_iter: int = 200, tol: float = 1e-6) -> np.ndarray:
    """ML π_i by iterative scaling. P(i≻j) = π_i / (π_i + π_j)."""
    n = W.shape[0]
    pi = np.ones(n, dtype=float)
    wins = W.sum(axis=1).astype(float)
    for _ in range(max_iter):
        new_pi = np.zeros_like(pi)
        for i in range(n):
            denom = 0.0
            for j in range(n):
                if i == j:
                    continue
                games = W[i, j] + W[j, i]
                if games > 0:
                    denom += games / (pi[i] + pi[j])
            new_pi[i] = wins[i] / denom if denom > 0 else pi[i]
        new_pi = new_pi / new_pi.sum() * n
        if np.max(np.abs(new_pi - pi)) < tol:
            pi = new_pi
            break
        pi = new_pi
    return pi


def main() -> None:
    ensure_phase2_dirs()
    tensor, layer_ids, node_ids = load_tensor()

    L, N, _ = tensor.shape
    ds_matrix = np.zeros((N, L), dtype=float)
    bt_matrix = np.zeros((N, L), dtype=float)
    for li, layer in enumerate(layer_ids):
        # Sender = winner. Transpose: W[i,j]=1 → i beat j.
        W = tensor[li].T.astype(int)
        ds_matrix[:, li] = davids_score(W.astype(float))
        bt_matrix[:, li] = bradley_terry_strengths(W)

    ds_df = pd.DataFrame(ds_matrix, index=node_ids, columns=layer_ids)
    bt_df = pd.DataFrame(bt_matrix, index=node_ids, columns=layer_ids)
    ds_df.to_csv(TABLE_DIR / "p2_11_davids_per_layer.csv")
    bt_df.to_csv(TABLE_DIR / "tableS3_dominance_concordance.csv")

    L = len(layer_ids)
    corr = np.zeros((L, L), dtype=float)
    for i in range(L):
        for j in range(L):
            rho, _ = spearmanr(ds_matrix[:, i], ds_matrix[:, j])
            corr[i, j] = rho if not np.isnan(rho) else 0.0
    corr_df = pd.DataFrame(corr, index=layer_ids, columns=layer_ids)
    corr_df.to_csv(TABLE_DIR / "p2_11_cross_layer_dominance_corr.csv")

    ds_z = StandardScaler().fit_transform(ds_matrix)
    pca = PCA(random_state=SEED)
    pca.fit(ds_z)
    explained = pca.explained_variance_ratio_
    cumvar = np.cumsum(explained)
    n_pcs_80 = int(np.searchsorted(cumvar, 0.8) + 1)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    im = axes[0].imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    axes[0].set_xticks(range(L)); axes[0].set_yticks(range(L))
    axes[0].set_xticklabels(layer_ids, rotation=40, ha="right")
    axes[0].set_yticklabels(layer_ids)
    axes[0].set_title("Cross-layer Spearman rank correlation\nof David's score")
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    axes[1].plot(np.arange(1, L + 1), explained, marker="o", label="Per-PC")
    axes[1].plot(np.arange(1, L + 1), cumvar, marker="s", label="Cumulative")
    axes[1].axhline(0.8, color="#94a3b8", linestyle="--", linewidth=0.8)
    axes[1].set_xlabel("Principal component")
    axes[1].set_ylabel("Variance explained")
    axes[1].set_title("PCA of layer-specific dominance")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="y", color="#e5e7eb", linewidth=0.5)
    save_figure(fig, "p2_11_dominance", FIG_DIR)
    plt.close(fig)

    off_diag = corr[np.triu_indices_from(corr, k=1)]
    mean_cross_layer_corr = float(np.mean(off_diag))
    median_cross_layer_corr = float(np.median(off_diag))

    var_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(L)],
        "variance_explained": explained,
        "cumulative_variance": cumvar,
    })

    body = f"""
Seed: {SEED}. Per-layer DS and BT; cross-layer Spearman; PCA on z-scored 60×{L} DS.

- mean off-diag Spearman ρ = {mean_cross_layer_corr:.3f}
- median = {median_cross_layer_corr:.3f}
- PCs to 80% var: {n_pcs_80} / {L}

{dataframe_to_md(var_df.head(8), index=False)}

ρ→0 and high n_pcs_80 → multidimensional dominance (h2). ρ→1, n_pcs_80=1 → one axis.

Outputs:
- p2_11_davids_per_layer.csv
- tableS3_dominance_concordance.csv
- p2_11_cross_layer_dominance_corr.csv
- p2_11_dominance.png
"""
    replace_section(REPORT_PATH, "<!-- P2.11 RESULTS -->", "P2.11 Dominance: per-layer and cross-layer", body)


if __name__ == "__main__":
    main()
