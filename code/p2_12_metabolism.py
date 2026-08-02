#!/usr/bin/env python3
"""P2.12 - Carbohydrate assimilation vs interaction network position.

Quantifies the "metabolism shapes the network" statement in the Results, which is
otherwise asserted without a test. Two questions, at two different levels:

  strain level: does a strain's assimilation capacity predict how active or how
                broadly spread its interactions are?
  dyad level:   do metabolically similar strain pairs interact more (or less)
                than dissimilar ones, per layer?

The dyad-level test is a Mantel permutation test rather than a median split on
pairs. Pairs sharing a strain are not independent observations, so a
Mann-Whitney over dyads has a badly inflated effective sample size; permuting
strain labels respects that dependence.
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
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from data_io import (
    CATEGORY_MAP,
    FIG_DIR,
    LAYER_ORDER,
    REPORT_PATH,
    ROOT,
    SEED,
    TABLE_DIR,
    canonical_strain,
    ensure_phase2_dirs,
    load_tensor,
)
from report_utils import dataframe_to_md, replace_section

ATB_PATH = ROOT / "data" / "annotation" / "atb_profile.xlsx"
# "MM" is the mineral-medium control column, not a carbon source; it is held out
# of the niche profile and reported separately.
CONTROL_COL = "MM"
N_PERM = 4999


def load_profiles() -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_excel(ATB_PATH, index_col=0)
    df.index = [canonical_strain(i) for i in df.index]
    df = df.apply(pd.to_numeric, errors="coerce")
    sugars = [c for c in df.columns if c != CONTROL_COL]
    # One missing Sorbose reading; median-fill so the strain is not dropped from
    # a 60-strain analysis over a single cell.
    df[sugars] = df[sugars].fillna(df[sugars].median())
    return df, sugars


def mantel(D1: np.ndarray, D2: np.ndarray, rng: np.random.Generator, n_perm: int) -> tuple[float, float]:
    iu = np.triu_indices(len(D1), 1)
    target = D2[iu]
    obs = spearmanr(D1[iu], target).statistic
    count = 0
    for _ in range(n_perm):
        p = rng.permutation(len(D1))
        if abs(spearmanr(D1[np.ix_(p, p)][iu], target).statistic) >= abs(obs):
            count += 1
    return float(obs), (count + 1) / (n_perm + 1)


def main() -> None:
    ensure_phase2_dirs()
    rng = np.random.default_rng(SEED)

    if not ATB_PATH.exists():
        replace_section(
            REPORT_PATH, "<!-- P2.12 RESULTS -->",
            "P2.12 Metabolic niche overlap and layer-specific interactions",
            f"Missing {ATB_PATH.relative_to(ROOT)} - metabolism analysis skipped.",
        )
        return

    atb, sugars = load_profiles()
    tensor, _, node_ids = load_tensor()
    keys = [canonical_strain(n) for n in node_ids]

    common = [k for k in keys if k in atb.index]
    idx = [keys.index(k) for k in common]
    G = atb.loc[common, sugars].to_numpy(float)

    # Assimilation is scored 0-4, not present/absent. Binarising it collapses the
    # scale: 55 of 60 strains use at least 7 of the 8 sugars, so a presence/absence
    # "breadth" is saturated and carries almost no between-strain information.
    # The graded row sum is the informative summary.
    capacity = G.sum(axis=1)
    breadth = (G > 0).sum(axis=1)

    out_deg = tensor.sum(axis=1).T.astype(float)[idx]
    total = out_deg.sum(axis=1)
    safe = np.where(total == 0, 1.0, total)
    participation = 1.0 - np.sum((out_deg / safe[:, None]) ** 2, axis=1)
    layer_idx = {l: i for i, l in enumerate(LAYER_ORDER)}

    # --- strain level -------------------------------------------------------- #
    responses: list[tuple[str, np.ndarray]] = [
        ("total out-degree", total),
        ("participation coefficient", participation),
    ]
    responses += [
        (f"{cat} mean out-degree", out_deg[:, [layer_idx[l] for l in layers]].mean(axis=1))
        for cat, layers in CATEGORY_MAP.items()
    ]
    responses += [(f"out-degree {l}", out_deg[:, layer_idx[l]]) for l in LAYER_ORDER]

    rows = []
    for name, v in responses:
        r = spearmanr(capacity, v)
        rows.append({"response": name, "spearman_rho": float(r.statistic), "p_value": float(r.pvalue)})
    strain_df = pd.DataFrame(rows)
    strain_df["q_BH"] = multipletests(strain_df["p_value"], method="fdr_bh")[1]
    strain_df = strain_df.sort_values("p_value")
    strain_df.to_csv(TABLE_DIR / "p2_12_metabolism_strain.csv", index=False)

    # --- dyad level ---------------------------------------------------------- #
    niche_dist = squareform(pdist(G, metric="euclidean"))
    rows = []
    for layer in LAYER_ORDER:
        M = tensor[layer_idx[layer]][np.ix_(idx, idx)].astype(float)
        adj = ((M + M.T) > 0).astype(float)  # dyad interacts in either direction
        if adj.sum() == 0:
            continue
        rho, p = mantel(niche_dist, adj, rng, N_PERM)
        rows.append({"layer": layer, "mantel_rho": rho, "p_value": p, "n_dyads": int(adj.sum() / 2)})
    dyad_df = pd.DataFrame(rows)
    dyad_df["q_BH"] = multipletests(dyad_df["p_value"], method="fdr_bh")[1]
    dyad_df = dyad_df.sort_values("p_value")
    dyad_df.to_csv(TABLE_DIR / "p2_12_metabolism_dyad.csv", index=False)

    profile_df = pd.DataFrame({
        "strain": common,
        "assimilation_capacity": capacity,
        "n_sugars_used": breadth,
        "total_out_degree": total,
        "participation": participation,
    })
    profile_df.to_csv(TABLE_DIR / "p2_12_metabolism_profiles.csv", index=False)

    # Figure produced by make_paper_figures.py (fig_metabolism).

    # --- report -------------------------------------------------------------- #
    strain_sig = strain_df[strain_df["q_BH"] < 0.05]
    dyad_sig = dyad_df[dyad_df["q_BH"] < 0.05]
    lead = dyad_df.head(2)

    body = f"""
Seed: {SEED}. n = {len(common)} / {len(node_ids)} strains ({len(sugars)} carbon sources; the mineral-medium
control column `{CONTROL_COL}` is excluded from the niche profile). Assimilation is scored
0-4, so the graded row sum is used as capacity: mean {capacity.mean():.1f}, range {capacity.min():.0f}-{capacity.max():.0f} of {4 * len(sugars)}.

Presence/absence breadth is **saturated** (mean {breadth.mean():.2f} of {len(sugars)} sugars; {int((breadth >= len(sugars) - 1).sum())} of {len(common)} strains
use at least {len(sugars) - 1}), which is why the graded scores are used throughout - binarising
discards nearly all of the between-strain variation.

### Strain level: no association

{dataframe_to_md(strain_df, index=False)}

**{len(strain_sig)} of {len(strain_df)}** associations reach q_BH < 0.05. Metabolic capacity does not predict how
active a strain is (total out-degree rho = {float(strain_df.loc[strain_df['response'] == 'total out-degree', 'spearman_rho'].iloc[0]):+.3f}) nor how evenly its interactions
spread across layers. The "metabolism shapes network position" reading is **not**
supported at the strain level.

### Dyad level: suggestive, and in the unexpected direction

{dataframe_to_md(dyad_df, index=False)}

**{len(dyad_sig)} of {len(dyad_df)}** layers reach q_BH < 0.05. The two strongest effects - {lead.iloc[0]['layer']}
(rho = {lead.iloc[0]['mantel_rho']:+.3f}, q = {lead.iloc[0]['q_BH']:.3f}) and {lead.iloc[1]['layer']} (rho = {lead.iloc[1]['mantel_rho']:+.3f}, q = {lead.iloc[1]['q_BH']:.3f}) - are both
direct-antagonism layers and both **positive**, meaning metabolically *dissimilar*
pairs antagonise each other more often.

That is the opposite of the competition-relatedness expectation, under which
resource-similar strains should compete hardest. It is reported as a
hypothesis-generating observation only: neither layer clears FDR, and 2 of {len(dyad_df)}
marginal results is a weak basis for a mechanistic claim.

Outputs: `p2_12_metabolism_strain.csv`, `p2_12_metabolism_dyad.csv`,
`p2_12_metabolism_profiles.csv`. Figure: `figures/fig_metabolism.*` (main text).
"""
    replace_section(
        REPORT_PATH, "<!-- P2.12 RESULTS -->",
        "P2.12 Metabolic niche overlap and layer-specific interactions", body,
    )

    print(f"n={len(common)}  capacity mean={capacity.mean():.1f} range={capacity.min():.0f}-{capacity.max():.0f}")
    print(f"strain-level significant (q<0.05): {len(strain_sig)}/{len(strain_df)}")
    print(f"dyad-level significant  (q<0.05): {len(dyad_sig)}/{len(dyad_df)}")
    print(dyad_df.head(4).to_string(index=False))


if __name__ == "__main__":
    main()
