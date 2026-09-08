#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

from data_io import save_figure, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, ensure_phase2_dirs, layer_groups, load_tensor
from null_models import generate_degree_preserving_nulls
from report_utils import dataframe_to_md, replace_section


N_DRAWS = 10000
# Methodological-confound pairs (same plate, observer, readout). Flagged, not biology.
CONFOUND_PAIRS = {
    frozenset({"CCAM", "CCVM"}),
}


def jaccard(A: np.ndarray, B: np.ndarray) -> float:
    inter = np.logical_and(A, B).sum()
    union = np.logical_or(A, B).sum()
    return float(inter / union) if union else np.nan


def main() -> None:
    ensure_phase2_dirs()
    tensor, layer_ids, _ = load_tensor()
    nulls = generate_degree_preserving_nulls(tensor, n_draws=N_DRAWS, seed=SEED)
    records = []
    for i, a in enumerate(layer_ids):
        A = tensor[i].astype(bool)
        for j in range(i + 1, len(layer_ids)):
            b = layer_ids[j]
            B = tensor[j].astype(bool)
            j_obs = jaccard(A, B)
            draws = np.array(
                [jaccard(nulls[a][k].astype(bool), nulls[b][k].astype(bool)) for k in range(nulls[a].shape[0])],
                dtype=float,
            )
            mu = float(np.nanmean(draws))
            sd = float(np.nanstd(draws, ddof=1))
            z = float((j_obs - mu) / sd) if sd > 0 else np.nan
            n_eff = np.sum(np.isfinite(draws))
            p_raw = float((np.sum(np.abs(draws - mu) >= abs(j_obs - mu)) + 1) / (n_eff + 1))
            records.append(
                {
                    "layer1": a,
                    "layer2": b,
                    "J_obs": j_obs,
                    "J_null_mean": mu,
                    "J_null_sd": sd,
                    "z": z,
                    "p_raw": p_raw,
                    "n_null_draws": n_eff,
                    "methodological_confound": frozenset({a, b}) in CONFOUND_PAIRS,
                    "seed": SEED,
                }
            )
    df = pd.DataFrame(records)
    reject_bh, p_bh = fdrcorrection(df["p_raw"].values, alpha=0.05, method="indep")
    reject_by, p_by = fdrcorrection(df["p_raw"].values, alpha=0.05, method="negcorr")
    df["p_fdr_bh"] = p_bh
    df["significant_bh"] = reject_bh
    df["p_fdr_by"] = p_by
    df["significant_by"] = reject_by
    df.to_csv(TABLE_DIR / "tableS2_layer_overlap_full.csv", index=False)

    bio_df = df[~df["methodological_confound"]]

    zmat = pd.DataFrame(np.nan, index=layer_ids, columns=layer_ids)
    pmat = pd.DataFrame(np.nan, index=layer_ids, columns=layer_ids)
    for _, row in df.iterrows():
        zmat.loc[row["layer1"], row["layer2"]] = row["z"]
        zmat.loc[row["layer2"], row["layer1"]] = row["z"]
        pmat.loc[row["layer1"], row["layer2"]] = row["p_fdr_by"]
        pmat.loc[row["layer2"], row["layer1"]] = row["p_fdr_by"]
    np.fill_diagonal(zmat.values, 0.0)
    np.fill_diagonal(pmat.values, np.nan)

    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    vmax = float(np.nanmax(np.abs(zmat.values)))
    im = ax.imshow(zmat.values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="equal")
    ax.set_xticks(np.arange(len(layer_ids)))
    ax.set_yticks(np.arange(len(layer_ids)))
    ax.set_xticklabels(layer_ids, rotation=45, ha="right")
    ax.set_yticklabels(layer_ids)
    ax.set_title("Jaccard overlap vs degree-preserving null (BY-FDR)")
    confound_layers = {layer for pair in CONFOUND_PAIRS for layer in pair}
    for i, ai in enumerate(layer_ids):
        for j, aj in enumerate(layer_ids):
            if i == j:
                continue
            val = zmat.iloc[i, j]
            p = pmat.iloc[i, j]
            is_confound = frozenset({ai, aj}) in CONFOUND_PAIRS
            if pd.notna(p) and p < 0.05 and not is_confound:
                mark = "*"
            elif is_confound and pd.notna(p) and p < 0.05:
                mark = "(c)"
            else:
                mark = ""
            ax.text(j, i, f"{val:.1f}{mark}", ha="center", va="center", fontsize=6.8, color="#111827")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard z-score")
    save_figure(fig, "p2_4_jaccard_heatmap", FIG_DIR)
    plt.close(fig)

    cc = df[((df["layer1"] == "CCAM") & (df["layer2"] == "CCVM")) | ((df["layer1"] == "CCVM") & (df["layer2"] == "CCAM"))]
    higher_bio = int(((bio_df["z"] > 0) & (bio_df["p_fdr_by"] < 0.05)).sum())
    lower_bio = int(((bio_df["z"] < 0) & (bio_df["p_fdr_by"] < 0.05)).sum())
    higher_bh = int(((bio_df["z"] > 0) & (bio_df["p_fdr_bh"] < 0.05)).sum())
    lower_bh = int(((bio_df["z"] < 0) & (bio_df["p_fdr_bh"] < 0.05)).sum())
    cc_by = bool(cc["p_fdr_by"].iloc[0] < 0.05) if not cc.empty else False
    cc_bh = bool(cc["p_fdr_bh"].iloc[0] < 0.05) if not cc.empty else False

    body = f"""
Seed: {SEED}. {N_DRAWS} degree-preserving directed-swap nulls per layer; Phipson–Smyth (M+1)/(N+1); FDR across 66 pairs by BH (p_fdr_bh) and BY (p_fdr_by, dependence-robust, default).

CCAM/CCVM (methodological confound, reported separately):
{dataframe_to_md(cc, index=False)}

Biological summary (CCAM/CCVM excluded, BY-FDR):
- higher_by = {higher_bio}
- lower_by  = {lower_bio}
- BH: higher={higher_bh}, lower={lower_bh}
- CCAM/CCVM sig BY = {cc_by} (BH: {cc_bh})

Reading: significant deviations are predominantly over-overlap, in mechanistically coherent clusters (e.g. IRAC/RAC, CCVM/CMP, IG/IRAC). The v14 "most pairs significantly lower than chance" claim is overturned.

Top positive biological z (BY-FDR):
{dataframe_to_md(bio_df.sort_values('z', ascending=False).head(10), index=False)}

Outputs:
- tableS2_layer_overlap_full.csv
- p2_4_jaccard_heatmap.png  (* = BY-FDR sig biological, (c) = confound)
"""
    replace_section(REPORT_PATH, "<!-- P2.4 RESULTS -->", "P2.4 Per-layer degree-preserving Jaccard null", body)


if __name__ == "__main__":
    main()
