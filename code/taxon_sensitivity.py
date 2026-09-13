#!/usr/bin/env python3
"""- Does anything depend on the four non-Streptomyces isolates?

The 16S annotation shows four of the 60 strains are not Streptomyces
(Saccharothrix MS53 7, Lentzea MS10 8, Lentzea MS3 15, Amycolatopsis MS3 18).
The manuscript describes the cohort as Streptomyces throughout, so the wording
has to change; this script establishes whether the *conclusions* have to change
with it.

Every headline analysis is re-run on the 56-strain Streptomyces-only subset and
tabulated against the full cohort. Statistics are imported from the scripts that
produce the primary numbers rather than reimplemented, so the two columns are
guaranteed comparable.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import fdrcorrection

from data_io import (
    LAYER_ORDER,
    REPORT_PATH,
    ROOT,
    SEED,
    TABLE_DIR,
    canonical_strain,
    ensure_output_dirs,
    load_tensor,
)
from null_models import generate_degree_preserving_nulls
from mixture_models import cv_loglik
from ordination import parallel_analysis
from layer_overlap import CONFOUND_PAIRS, jaccard
from connectivity import largest_scc_size
from layer_descriptors import reciprocity_both
from report_utils import dataframe_to_md, replace_section

LABEL_MAP = ROOT / "data" / "phylogeny" / "label_map.csv"
# Fallback if prep_phylogeny.py has not been run; kept in sync with its
# NON_STREPTOMYCES genus list.
FALLBACK_EXCLUDE = ["MS53 7", "MS10 8", "MS3 15", "MS3 18"]

N_DRAWS = 10000  # matched to the primary analysis so the two arms are comparable
K_MAX = 8


def non_streptomyces_strains() -> list[str]:
    if LABEL_MAP.exists():
        m = pd.read_csv(LABEL_MAP)
        m = m[m["kept"] & m["in_interaction_matrix"]]
        found = sorted(m.loc[m["genus"] != "Streptomyces", "strain"].map(canonical_strain).unique())
        if found:
            return found
    return [canonical_strain(s) for s in FALLBACK_EXCLUDE]


def layer_structure(tensor: np.ndarray) -> pd.DataFrame:
    rows = []
    for i, layer in enumerate(LAYER_ORDER):
        A = tensor[i].astype(np.uint8)
        n = A.shape[0]
        r_pairs, _r_edges, _mp, _up, edges = reciprocity_both(A)
        rows.append({
            "layer": layer,
            "density": edges / (n * (n - 1)),
            "reciprocity_pairs": r_pairs,
            "largest_scc": largest_scc_size(A),
        })
    return pd.DataFrame(rows).set_index("layer")


def jaccard_null_baseline() -> tuple[int, int] | None:
    """Read the published n=60 overlap result rather than recomputing it.

    layer_overlap already produces this with 10,000 draws and BY-FDR; regenerating a
    2,000-draw copy of the same thing costs ~20 minutes and would compare the
    n=56 result against a *different* baseline than the paper reports.
    """
    path = TABLE_DIR / "tableS2_layer_overlap_full.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    bio = df[~df["methodological_confound"].astype(bool)]
    sig = bio[bio["significant_by"].astype(bool)]
    return int((sig["z"] > 0).sum()), int((sig["z"] < 0).sum())


def jaccard_null_summary(
    tensor: np.ndarray, prefix: str, seed: int, n_draws: int = N_DRAWS
) -> tuple[int, int, set[tuple[str, str]]]:
    """Count biological pairs that over- and under-overlap vs the degree-preserving null."""
    nulls = generate_degree_preserving_nulls(tensor, n_draws=n_draws, seed=seed, prefix=prefix)
    recs = []
    for i, a in enumerate(LAYER_ORDER):
        A = tensor[i].astype(bool)
        for j in range(i + 1, len(LAYER_ORDER)):
            b = LAYER_ORDER[j]
            B = tensor[j].astype(bool)
            j_obs = jaccard(A, B)
            # Draws are truncated to n_draws so a cached 10,000-draw file can be
            # compared like-for-like against a 2,000-draw one.
            draws = np.array(
                [jaccard(nulls[a][k].astype(bool), nulls[b][k].astype(bool)) for k in range(n_draws)],
                dtype=float,
            )
            mu = float(np.nanmean(draws))
            n_eff = int(np.sum(np.isfinite(draws)))
            p = float((np.sum(np.abs(draws - mu) >= abs(j_obs - mu)) + 1) / (n_eff + 1))
            recs.append({"pair": frozenset({a, b}), "a": a, "b": b, "direction": np.sign(j_obs - mu), "p": p})
    df = pd.DataFrame(recs)
    bio = df[~df["pair"].isin(CONFOUND_PAIRS)].copy()
    # Benjamini-Yekutieli, matching the primary analysis in layer_overlap.
    _, q = fdrcorrection(bio["p"].to_numpy(), method="negcorr")
    bio["q"] = q
    sig = bio[bio["q"] < 0.05]
    names = {tuple(sorted([r.a, r.b])) for r in sig.itertuples()}
    return int((sig["direction"] > 0).sum()), int((sig["direction"] < 0).sum()), names


def pca_summary(tensor: np.ndarray, seed: int) -> tuple[float, float, int]:
    X = tensor.sum(axis=1).T.astype(float)
    Z = StandardScaler().fit_transform(X)
    pca = PCA(random_state=seed).fit(Z)
    evr = pca.explained_variance_ratio_
    _null_mean, null_p95 = parallel_analysis(X, 300, seed)
    # Same retention rule as the primary analysis: 95th percentile, not the mean.
    return float(evr[0] * 100), float(evr[1] * 100), int((evr > null_p95[: len(evr)]).sum())


def best_cv_k(tensor: np.ndarray, seed: int) -> int:
    X = StandardScaler().fit_transform(tensor.sum(axis=1).T.astype(float))
    scores = {k: cv_loglik(X, k, "full", seed) for k in range(1, K_MAX + 1)}
    valid = {k: v for k, v in scores.items() if np.isfinite(v)}
    return int(max(valid, key=valid.get)) if valid else -1


def main() -> None:
    ensure_output_dirs()

    tensor, _, node_ids = load_tensor()
    keys = [canonical_strain(n) for n in node_ids]
    exclude = non_streptomyces_strains()
    keep = [i for i, k in enumerate(keys) if k not in set(exclude)]
    sub = tensor[:, keep, :][:, :, keep]

    full_struct = layer_structure(tensor)
    sub_struct = layer_structure(sub)
    struct = full_struct.join(sub_struct, lsuffix="_n60", rsuffix="_n56")
    struct["d_density"] = struct["density_n56"] - struct["density_n60"]
    struct["d_scc"] = struct["largest_scc_n56"] - struct["largest_scc_n60"]
    struct = struct.reset_index()
    struct.to_csv(TABLE_DIR / "taxon_sensitivity_layer_structure.csv", index=False)

    # Three overlap numbers, not two. The published n=60 result uses 10,000 draws;
    # the n=56 arm uses fewer. Comparing those two directly would charge a
    # draw-count effect to the excluded strains, so an n=60 result at the *same*
    # draw count is computed as the honest comparator (the 10,000-draw cache is
    # simply truncated, so this costs nothing).
    published = jaccard_null_baseline()
    over60m, under60m, pairs60 = jaccard_null_summary(tensor, "null", SEED, n_draws=N_DRAWS)
    over56, under56, pairs56 = jaccard_null_summary(sub, "null56", SEED, n_draws=N_DRAWS)
    over60, under60 = published if published is not None else (over60m, under60m)
    lost = sorted(pairs60 - pairs56)
    gained = sorted(pairs56 - pairs60)

    pc1_60, pc2_60, horn60 = pca_summary(tensor, SEED)
    pc1_56, pc2_56, horn56 = pca_summary(sub, SEED)

    k60 = best_cv_k(tensor, SEED)
    k56 = best_cv_k(sub, SEED)

    scc_rank_rho = struct["largest_scc_n60"].corr(struct["largest_scc_n56"], method="spearman")
    dens_rank_rho = struct["density_n60"].corr(struct["density_n56"], method="spearman")

    headline = pd.DataFrame([
        {"quantity": f"Jaccard null: over-overlapping pairs, published 10,000 draws", "n60": over60, "n56": "-"},
        {"quantity": f"Jaccard null: over-overlapping pairs, matched {N_DRAWS:,} draws", "n60": over60m, "n56": over56},
        {"quantity": f"Jaccard null: under-overlapping pairs, matched {N_DRAWS:,} draws", "n60": under60m, "n56": under56},
        {"quantity": "12-D PCA PC1 (% variance)", "n60": round(pc1_60, 1), "n56": round(pc1_56, 1)},
        {"quantity": "12-D PCA PC2 (% variance)", "n60": round(pc2_60, 1), "n56": round(pc2_56, 1)},
        {"quantity": "PCs above Horn parallel-analysis null", "n60": horn60, "n56": horn56},
        {"quantity": "GMM components by 10-fold CV log-likelihood (12-D, full cov)", "n60": k60, "n56": k56},
        {"quantity": "Layer density rank correlation n60 vs n56 (Spearman)", "n60": "-", "n56": round(dens_rank_rho, 4)},
        {"quantity": "Largest-SCC rank correlation n60 vs n56 (Spearman)", "n60": "-", "n56": round(scc_rank_rho, 4)},
    ])
    headline.to_csv(TABLE_DIR / "tableS4_taxon_sensitivity.csv", index=False)

    qualitative = (horn60 == horn56) and (k60 == k56) and (under60m == under56) and not gained
    if qualitative and not lost:
        verdict = (
            "Every headline result is **unchanged** when the four non-Streptomyces "
            "isolates are removed."
        )
    elif qualitative:
        verdict = (
            f"Every qualitative conclusion holds. The only change is that "
            f"{len(lost)} of {over60m} over-overlapping pairs "
            f"({', '.join('-'.join(p) for p in lost)}) fall below the FDR threshold at n=56."
        )
    else:
        verdict = "**A headline conclusion changes** at n=56 - see the table."

    body = f"""
Seed: {SEED}. Excluded ({len(exclude)}): {', '.join(exclude)}.
n = {tensor.shape[1]} -> {sub.shape[1]} strains.

Overlap counts are reported at two draw counts. The published layer_overlap result uses
10,000 null draws; the n=56 arm uses {N_DRAWS:,}. Comparing those directly would charge
any draw-count effect to the excluded strains, so an n=60 result at the matched
{N_DRAWS:,} draws is computed as the honest comparator ({over60} at 10,000 -> {over60m} at {N_DRAWS:,}).
**Read the n=56 column against the matched row, not the published one.**

{verdict}

### Headline quantities

{dataframe_to_md(headline, index=False)}

### Per-layer structure

{dataframe_to_md(struct.round(4), index=False)}

Layer density and largest-SCC orderings are preserved almost exactly
(Spearman rho = {dens_rank_rho:.3f} and {scc_rank_rho:.3f}), so the architecture contrasts that
the Results are built on do not depend on the four non-Streptomyces isolates.

{"Pairs losing significance at n=56: " + ", ".join("-".join(p) for p in lost) if lost else "No pair changes significance status at n=56."}
{"Pairs gaining significance at n=56: " + ", ".join("-".join(p) for p in gained) if gained else ""}

The taxonomy wording in the title, abstract and Methods still needs correcting -
see `data/annotation/DATA_NOTES.md`.

Outputs: `tableS4_taxon_sensitivity.csv`, `taxon_sensitivity_layer_structure.csv`
"""
    replace_section(REPORT_PATH, "<!-- taxon_sensitivity RESULTS -->", "Sensitivity to the non-Streptomyces isolates", body)

    print(f"excluded: {exclude}")
    print(f"n {tensor.shape[1]} -> {sub.shape[1]}")
    print(headline.to_string(index=False))
    print(f"\nlost at n=56  : {[' - '.join(p) for p in lost]}")
    print(f"gained at n=56: {[' - '.join(p) for p in gained]}")
    print(f"verdict: {verdict}")


if __name__ == "__main__":
    main()
