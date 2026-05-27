#!/usr/bin/env python3
"""P2.12 — Niche overlap (sugar) vs DA/IA/MM/MC interaction frequency. Breadth ↔ participation. Skips if metadata absent."""
from __future__ import annotations

import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from data_io import CATEGORY_MAP, FIG_DIR, REPORT_PATH, ROOT, SEED, TABLE_DIR, ensure_phase2_dirs, load_tensor
from report_utils import dataframe_to_md, replace_section


GROWTH_PATH = ROOT / "data" / "metadata" / "Y-_mu_x.xlsx"
ATB_PATH = ROOT / "data" / "metadata" / "Y-_ATB.xlsx"


def _placeholder_report() -> str:
    return f"""
Metabolism metadata missing.

Expected:
- {GROWTH_PATH.relative_to(ROOT)} (growth rate per sugar)
- {ATB_PATH.relative_to(ROOT)} (antimicrobial release per sugar)

When present: pairwise Jaccard niche overlap; per-category interaction frequency (high vs low overlap, Mann–Whitney); breadth vs participation (Spearman). Writes p2_12_metabolism_*.csv + figure.
"""


def _normalize_strain(s: object) -> str:
    s = str(s).strip()
    s = s.replace("(", " ").replace(")", " ")
    return " ".join(s.split())


def main() -> None:
    ensure_phase2_dirs()
    if not GROWTH_PATH.exists():
        replace_section(REPORT_PATH, "<!-- P2.12 RESULTS -->", "P2.12 Metabolic niche overlap and layer-specific interactions", _placeholder_report())
        return

    growth = pd.read_excel(GROWTH_PATH)
    if growth.shape[1] >= 2 and "Unnamed" in growth.columns[1]:
        growth = growth.rename(columns={growth.columns[1]: "strain"})
        growth = growth.drop(columns=[growth.columns[0]])
    else:
        growth = growth.rename(columns={growth.columns[0]: "strain"})
    growth["strain"] = growth["strain"].map(_normalize_strain)
    growth = growth.set_index("strain")
    growth = growth.apply(pd.to_numeric, errors="coerce").fillna(0)

    tensor, layer_ids, node_ids = load_tensor()

    common = sorted(set(growth.index) & set(node_ids))
    if len(common) < 20:
        body = (
            f"Insufficient overlap: sugar {len(growth)}, matrices {len(node_ids)}, common {len(common)}. "
            "Check naming (e.g. MS3(18) → MS3 18)."
        )
        replace_section(REPORT_PATH, "<!-- P2.12 RESULTS -->", "P2.12 Metabolic niche overlap and layer-specific interactions", body)
        return

    growth = growth.loc[common]
    sugar_use = (growth > 0).astype(int)
    n_sugars_per_strain = sugar_use.sum(axis=1)

    idx_map = {s: i for i, s in enumerate(common)}
    n = len(common)
    overlap = np.zeros((n, n), dtype=float)
    for i, si in enumerate(common):
        ai = sugar_use.loc[si].values
        for j, sj in enumerate(common):
            aj = sugar_use.loc[sj].values
            u = np.logical_or(ai, aj).sum()
            if u > 0:
                overlap[i, j] = np.logical_and(ai, aj).sum() / u

    layer_idx = {l: li for li, l in enumerate(layer_ids)}
    cat_pair_records = []
    for cat, layers in CATEGORY_MAP.items():
        cat_layer_idx = [layer_idx[l] for l in layers if l in layer_idx]
        agg = np.zeros((tensor.shape[1], tensor.shape[2]), dtype=int)
        for li in cat_layer_idx:
            agg = agg | (tensor[li] > 0).astype(int)
        restrict = [node_ids.index(s) for s in common]
        agg_sub = agg[np.ix_(restrict, restrict)]
        iu = np.triu_indices(n, k=1)
        overlap_pairs = overlap[iu]
        cat_pairs = agg_sub[iu] | agg_sub.T[iu]
        cat_pairs_bool = (cat_pairs > 0)
        median_o = np.median(overlap_pairs)
        high = cat_pairs_bool[overlap_pairs >= median_o]
        low = cat_pairs_bool[overlap_pairs < median_o]
        if len(high) == 0 or len(low) == 0:
            continue
        try:
            stat, p_mw = mannwhitneyu(high.astype(int), low.astype(int), alternative="greater")
        except ValueError:
            stat, p_mw = np.nan, np.nan
        cat_pair_records.append({
            "category": cat,
            "n_pairs_high_overlap": int(high.sum()),
            "n_pairs_total_high_overlap": int(len(high)),
            "frac_interact_high_overlap": float(high.mean()),
            "n_pairs_low_overlap": int(low.sum()),
            "n_pairs_total_low_overlap": int(len(low)),
            "frac_interact_low_overlap": float(low.mean()),
            "diff_frac": float(high.mean() - low.mean()),
            "mannwhitney_p": float(p_mw),
        })

    cat_df = pd.DataFrame(cat_pair_records)
    cat_df.to_csv(TABLE_DIR / "p2_12_metabolism_category.csv", index=False)

    out_deg = tensor.sum(axis=1).T
    total_out = out_deg.sum(axis=1)
    safe_total = np.where(total_out == 0, 1, total_out)
    participation = 1.0 - np.sum((out_deg / safe_total[:, None]) ** 2, axis=1)
    breadth_df = pd.DataFrame({
        "strain": node_ids,
        "metabolic_breadth": [int(n_sugars_per_strain.get(s, 0)) for s in node_ids],
        "participation": participation,
    })
    breadth_df = breadth_df[breadth_df["strain"].isin(common)]
    from scipy.stats import spearmanr
    rho_breadth = spearmanr(breadth_df["metabolic_breadth"], breadth_df["participation"])

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5), constrained_layout=True)
    if not cat_df.empty:
        x = np.arange(len(cat_df))
        w = 0.4
        axes[0].bar(x - w / 2, cat_df["frac_interact_low_overlap"], width=w, label="Low overlap", color="#94a3b8")
        axes[0].bar(x + w / 2, cat_df["frac_interact_high_overlap"], width=w, label="High overlap", color="#0f766e")
        axes[0].set_xticks(x); axes[0].set_xticklabels(cat_df["category"])
        axes[0].set_ylabel("Fraction of pairs interacting")
        axes[0].set_title("Niche overlap and interaction frequency, by category")
        axes[0].legend(frameon=False)
        axes[0].grid(axis="y", color="#e5e7eb", linewidth=0.6)
    axes[1].scatter(breadth_df["metabolic_breadth"], breadth_df["participation"], color="#4c78a8", s=40, edgecolor="white")
    axes[1].set_xlabel("Number of utilized sugars (metabolic breadth)")
    axes[1].set_ylabel("Network participation coefficient")
    axes[1].set_title(f"Spearman rho = {rho_breadth.statistic:.2f}, p = {rho_breadth.pvalue:.3g}")
    axes[1].grid(color="#e5e7eb", linewidth=0.6)
    fig.savefig(FIG_DIR / "p2_12_metabolism.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    breadth_df.to_csv(TABLE_DIR / "p2_12_metabolism_strain.csv", index=False)

    body = f"""
Seed: {SEED}. Jaccard niche overlap on utilized-sugar sets; per-category (DA/IA/MM/MC) interaction-frequency by overlap median-split (Mann–Whitney).

Strain overlap: {len(common)} / {len(node_ids)}.

Per-category:
{dataframe_to_md(cat_df, index=False)}

Breadth vs participation: Spearman ρ = {rho_breadth.statistic:.3f}, p = {rho_breadth.pvalue:.3g}.

diff_frac>0 (DA only) supports niche-overlap → constitutive antagonism (paper §6).

Outputs:
- p2_12_metabolism_category.csv
- p2_12_metabolism_strain.csv
- p2_12_metabolism.png
"""
    replace_section(REPORT_PATH, "<!-- P2.12 RESULTS -->", "P2.12 Metabolic niche overlap and layer-specific interactions", body)


if __name__ == "__main__":
    main()
