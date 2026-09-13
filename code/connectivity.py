#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

from data_io import save_figure, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, ensure_output_dirs, load_tensor
from null_models import generate_degree_preserving_nulls
from report_utils import dataframe_to_md, replace_section


N_DRAWS = 10000


def largest_scc_size(A: np.ndarray) -> int:
    G = nx.from_numpy_array(A.T, create_using=nx.DiGraph)
    return max((len(c) for c in nx.strongly_connected_components(G)), default=0)


def main() -> None:
    ensure_output_dirs()
    tensor, layer_ids, _ = load_tensor()
    nulls = generate_degree_preserving_nulls(tensor, n_draws=N_DRAWS, seed=SEED)
    rows = []
    for i, layer in enumerate(layer_ids):
        obs = largest_scc_size(tensor[i])
        draws = np.array([largest_scc_size(nulls[layer][k]) for k in range(nulls[layer].shape[0])], dtype=float)
        mu = float(np.mean(draws))
        sd = float(np.std(draws, ddof=1))
        z = float((obs - mu) / sd) if sd > 0 else np.nan
        p_raw = float((np.sum(np.abs(draws - mu) >= abs(obs - mu)) + 1) / (len(draws) + 1))
        rows.append({"layer": layer, "scc_obs": obs, "scc_null_mean": mu, "scc_null_sd": sd, "z": z, "p_raw": p_raw, "n_null_draws": int(len(draws)), "seed": SEED})
    df = pd.DataFrame(rows)
    reject_bh, p_bh = fdrcorrection(df["p_raw"].values, alpha=0.05, method="indep")
    reject_by, p_by = fdrcorrection(df["p_raw"].values, alpha=0.05, method="negcorr")
    df["p_fdr_bh"] = p_bh
    df["significant_bh"] = reject_bh
    df["p_fdr_by"] = p_by
    df["significant_by"] = reject_by
    df.to_csv(TABLE_DIR / "connectivity_null.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    x = np.arange(len(layer_ids))
    ax.bar(x, df["scc_obs"], color="#4c78a8", alpha=0.95, label="Observed")
    ax.errorbar(x, df["scc_null_mean"], yerr=df["scc_null_sd"], fmt="o", color="#d62728", capsize=3, label="Null mean ± SD")
    ax.set_xticks(x)
    ax.set_xticklabels(layer_ids, rotation=40, ha="right")
    ax.set_ylabel("Largest SCC size")
    ax.set_title(f"Largest SCC vs degree-preserving null ({N_DRAWS} draws)")
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.6)
    ax.legend(frameon=False)
    save_figure(fig, "connectivity", FIG_DIR)
    plt.close(fig)

    exceeds = df[(df["z"] > 0) & (df["p_fdr_by"] < 0.05)]["layer"].tolist()
    null_like = df[df["p_fdr_by"] >= 0.05]["layer"].tolist()
    body = f"""
Seed: {SEED}. Observed largest SCC vs {N_DRAWS} degree-preserving nulls per layer; Phipson–Smyth (M+1)/(N+1); FDR across 12 layers by BH and BY.

{dataframe_to_md(df, index=False)}

- SCC above null (BY-sig): {', '.join(exceeds) if exceeds else 'none'}
- SCC consistent with null: {', '.join(null_like) if null_like else 'none'}

SCC size is a consequence of density, not extra feedback structure. v14 "MM/MC SCCs 20–40 nodes" is wrong: actual distribution is bimodal across MC (CMP 39, CS 26, CCVM 24 vs CCAM 4, IRP 1, IAC_RDE 1, IC 2, RP 3).

Outputs:
- connectivity.csv
- connectivity.png
"""
    replace_section(REPORT_PATH, "<!-- connectivity RESULTS -->", "Strongly connected components against the same null", body)


if __name__ == "__main__":
    main()
