#!/usr/bin/env python3
"""P2.14 — Reciprocity and transitivity against the degree-preserving null.

The manuscript claimed phenotype-specific deviations in reciprocity ("IG z
approaching 6") and depleted transitivity in six layers. Neither quantity had a
null model anywhere in the codebase; only the SCC did (p2_5). This stage supplies
both, reusing the cached draws so it costs nothing on top of p2_4/p2_5.

Directed edge swaps preserve in- and out-degree but not which of a dyad's two
arcs is present, so the same draws are a valid null for both statistics.
"""
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

from data_io import save_figure, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, ensure_phase2_dirs, load_tensor
from null_models import generate_degree_preserving_nulls
from p2_7_proximity_sensitivity import reciprocity_both
from report_utils import dataframe_to_md, replace_section


N_DRAWS = 10000


def transitivity(A: np.ndarray) -> float:
    """Global clustering coefficient of the undirected projection.

    ``nx.transitivity`` (closed triplets / all triplets) is the quantity the
    manuscript names. ``generate_figs.compute_layer_stats`` reports the mean
    local coefficient instead; both are tabulated so the two cannot drift.

    Computed in closed form rather than through networkx: this runs 120,000 times
    (12 layers x 10,000 draws) and building a graph object each time dominates the
    stage. For a symmetric binary U with zero diagonal,
        triangles = tr(U^3) / 6,   triads = sum_i d_i(d_i - 1) / 2,
        transitivity = 3 * triangles / triads = tr(U^3) / sum_i d_i(d_i - 1).
    ``_transitivity_reference`` checks this against networkx.
    """
    U = ((A + A.T) > 0).astype(np.float64)
    np.fill_diagonal(U, 0.0)
    d = U.sum(axis=1)
    triads = float(np.sum(d * (d - 1.0)))
    if triads == 0.0:
        return 0.0
    return float(np.trace(U @ U @ U) / triads)


def _transitivity_reference(A: np.ndarray) -> float:
    G = nx.from_numpy_array(A.T, create_using=nx.DiGraph).to_undirected()
    return float(nx.transitivity(G)) if G.number_of_nodes() else 0.0


def reciprocity_pairs(A: np.ndarray) -> float:
    """R_p = mutual / connected dyads, vectorized.

    Same quantity as ``p2_7.reciprocity_both``[0], which loops over all 1,770 dyads
    in Python. That is fine once per layer and far too slow 120,000 times.
    ``main`` asserts the two agree on every observed layer before using this one.
    """
    B = A.astype(bool)
    mutual = int(np.count_nonzero(B & B.T)) // 2
    connected = int(np.count_nonzero(B | B.T)) // 2
    return mutual / connected if connected else 0.0


def mean_local_clustering(A: np.ndarray) -> float:
    G = nx.from_numpy_array(A.T, create_using=nx.DiGraph).to_undirected()
    return float(nx.average_clustering(G)) if G.number_of_nodes() else 0.0


def null_summary(obs: float, draws: np.ndarray) -> dict[str, float]:
    mu = float(np.mean(draws))
    sd = float(np.std(draws, ddof=1))
    z = float((obs - mu) / sd) if sd > 0 else np.nan
    # Phipson-Smyth (M+1)/(N+1), two-sided, as everywhere else in the bundle.
    p_raw = float((np.sum(np.abs(draws - mu) >= abs(obs - mu)) + 1) / (len(draws) + 1))
    return {"obs": obs, "null_mean": mu, "null_sd": sd, "z": z, "p_raw": p_raw}


def main() -> None:
    ensure_phase2_dirs()
    tensor, layer_ids, _ = load_tensor()
    nulls = generate_degree_preserving_nulls(tensor, n_draws=N_DRAWS, seed=SEED)

    rows = []
    for i, layer in enumerate(layer_ids):
        A = tensor[i].astype(int)
        draws = nulls[layer]
        n_draws = draws.shape[0]

        # Guard the fast paths against the reference implementations, per layer.
        r_obs, _re, _mp, _up, _e = reciprocity_both(A)
        assert abs(reciprocity_pairs(A) - r_obs) < 1e-12, f"{layer}: reciprocity fast path disagrees"
        t_obs = transitivity(A)
        assert abs(_transitivity_reference(A) - t_obs) < 1e-9, f"{layer}: transitivity fast path disagrees"

        r_draws = np.array([reciprocity_pairs(draws[k]) for k in range(n_draws)], dtype=float)
        t_draws = np.array([transitivity(draws[k]) for k in range(n_draws)], dtype=float)

        rec = null_summary(r_obs, r_draws)
        tra = null_summary(t_obs, t_draws)
        rows.append(
            {
                "layer": layer,
                "reciprocity_obs": rec["obs"],
                "reciprocity_null_mean": rec["null_mean"],
                "reciprocity_null_sd": rec["null_sd"],
                "reciprocity_z": rec["z"],
                "reciprocity_p_raw": rec["p_raw"],
                "transitivity_obs": tra["obs"],
                "transitivity_null_mean": tra["null_mean"],
                "transitivity_null_sd": tra["null_sd"],
                "transitivity_z": tra["z"],
                "transitivity_p_raw": tra["p_raw"],
                "mean_local_clustering_obs": mean_local_clustering(A),
                "n_null_draws": int(n_draws),
                "seed": SEED,
            }
        )
    df = pd.DataFrame(rows)

    # Two separate families of 12; a degenerate null (sd = 0) cannot be tested
    # and is carried through as q = 1 rather than dropped, so the table stays 12 rows.
    for stat in ("reciprocity", "transitivity"):
        p = df[f"{stat}_p_raw"].values.astype(float)
        testable = np.isfinite(df[f"{stat}_z"].values.astype(float))
        q_bh = np.ones_like(p)
        q_by = np.ones_like(p)
        if testable.any():
            _, q_bh[testable] = fdrcorrection(p[testable], alpha=0.05, method="indep")
            _, q_by[testable] = fdrcorrection(p[testable], alpha=0.05, method="negcorr")
        df[f"{stat}_q_bh"] = q_bh
        df[f"{stat}_q_by"] = q_by
        df[f"{stat}_significant_by"] = testable & (q_by < 0.05)

    df.to_csv(TABLE_DIR / "p2_14_reciprocity_transitivity_null.csv", index=False)

    # Three panels: SCC, reciprocity, transitivity. The SCC column is read from
    # p2_5's output rather than recomputed, so the two stages cannot disagree, and
    # all three degree-preserving comparisons appear in one figure.
    scc = pd.read_csv(TABLE_DIR / "p2_5_scc_null.csv").set_index("layer")
    panels = [
        (scc.loc[layer_ids, "z"].values.astype(float),
         scc.loc[layer_ids, "significant_by"].values.astype(bool),
         "Largest SCC"),
        (df["reciprocity_z"].values.astype(float),
         df["reciprocity_significant_by"].values.astype(bool),
         "Reciprocity $R_p$"),
        (df["transitivity_z"].values.astype(float),
         df["transitivity_significant_by"].values.astype(bool),
         "Transitivity"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    x = np.arange(len(layer_ids))
    for ax, (z, sig, label) in zip(axes, panels):
        colors = ["#b91c1c" if s and v > 0 else "#1d4ed8" if s else "#cbd5e1" for s, v in zip(sig, z)]
        ax.bar(x, np.nan_to_num(z), color=colors)
        ax.axhline(0, color="#334155", linewidth=0.8)
        for thresh in (-1.96, 1.96):
            ax.axhline(thresh, color="#94a3b8", linewidth=0.7, linestyle="--")
        ax.set_xticks(x)
        ax.set_xticklabels(layer_ids, rotation=40, ha="right")
        ax.set_ylabel(f"z ({label})")
        ax.set_title(f"{label} vs degree-preserving null")
        ax.set_ylim(min(-2.6, np.nanmin(z) * 1.15), max(2.6, np.nanmax(z) * 1.15))
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.6)
    fig.suptitle(f"{N_DRAWS:,} degree-preserving draws per layer; filled bars significant after BY-FDR "
                 f"within each statistic's family of twelve", fontsize=9.5, y=1.02)
    save_figure(fig, "figureS4_network_nulls", FIG_DIR)
    plt.close(fig)

    def named(stat: str, sign: str) -> str:
        s = df[df[f"{stat}_significant_by"]]
        s = s[s[f"{stat}_z"] > 0] if sign == "over" else s[s[f"{stat}_z"] < 0]
        return ", ".join(f"{r.layer} (z = {getattr(r, f'{stat}_z'):+.2f})" for r in s.itertuples()) or "none"

    degenerate = df[~np.isfinite(df["reciprocity_z"].astype(float))]["layer"].tolist()
    body = f"""
Seed: {SEED}. Observed reciprocity ($R_p$ = mutual / connected dyads) and transitivity
(closed triplets / all triplets, undirected projection) vs {N_DRAWS:,} degree-preserving
directed-swap nulls per layer; Phipson–Smyth (M+1)/(N+1); FDR across 12 layers by BH and
BY, computed separately for each statistic.

{dataframe_to_md(df, index=False)}

- Reciprocity above null (BY-sig): {named('reciprocity', 'over')}
- Reciprocity below null (BY-sig): {named('reciprocity', 'under')}
- Transitivity above null (BY-sig): {named('transitivity', 'over')}
- Transitivity below null (BY-sig): {named('transitivity', 'under')}
- Degenerate null (no mutual dyad in any draw, z undefined): {', '.join(degenerate) if degenerate else 'none'}

These two statistics had no null model before this stage. Any manuscript sentence
about reciprocity or transitivity "relative to degree-preserving expectations"
must cite this table and no other.

Outputs:
- p2_14_reciprocity_transitivity_null.csv
- figureS4_network_nulls.png
"""
    replace_section(REPORT_PATH, "<!-- P2.14 RESULTS -->", "P2.14 Reciprocity and transitivity nulls", body)


if __name__ == "__main__":
    main()
