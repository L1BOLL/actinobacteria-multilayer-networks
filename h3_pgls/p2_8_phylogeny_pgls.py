#!/usr/bin/env python3
"""P2.8 — PGLS with Pagel-λ on ecological traits. Tests h3. Skips if tree.nwk absent."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd

from data_io import REPORT_PATH, ROOT, SEED, TABLE_DIR, ensure_phase2_dirs, load_or_recompute_embedding4, load_tensor
from report_utils import dataframe_to_md, replace_section


TREE_PATH = ROOT / "data" / "phylogeny" / "tree.nwk"
TRAITS_OF_INTEREST = ["DA", "IA", "MM", "MC", "total_out_degree", "participation"]


def _placeholder_report() -> str:
    return f"""
Tree missing: {TREE_PATH.relative_to(ROOT)}

Required: Newick tree at data/phylogeny/tree.nwk, tips matching 60 strain IDs.
(MAFFT + RAxML GTR+GAMMA, 1000 bootstraps — as Methods.)

When tree present: per-trait Pagel-λ PGLS over {{{', '.join(TRAITS_OF_INTEREST)}}}; LR vs λ=0; writes p2_8_pgls.csv.

λ≈0 → independent of phylogeny; λ≈1 → heritable. Abstract λ=0.31–0.48 currently from PI's prior work.
"""


def main() -> None:
    ensure_phase2_dirs()
    if not TREE_PATH.exists():
        body = _placeholder_report()
        replace_section(REPORT_PATH, "<!-- P2.8 RESULTS -->", "P2.8 Phylogenetic signal (PGLS)", body)
        return

    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        body = (
            f"ete3 not installed. Install: pip install ete3. Tree at {TREE_PATH}."
        )
        replace_section(REPORT_PATH, "<!-- P2.8 RESULTS -->", "P2.8 Phylogenetic signal (PGLS)", body)
        return

    from scipy.optimize import minimize_scalar
    from scipy.stats import chi2

    tree = Tree(str(TREE_PATH), format=1)
    tensor, layer_ids, node_ids = load_tensor()
    emb4 = load_or_recompute_embedding4().loc[node_ids]

    out_deg = tensor.sum(axis=1).T
    total_out = out_deg.sum(axis=1)
    out_total = total_out.copy()
    out_total[out_total == 0] = 1
    participation = 1.0 - np.sum((out_deg / out_total[:, None]) ** 2, axis=1)
    traits = pd.DataFrame({
        "strain": node_ids,
        "DA": emb4["DA"].values,
        "IA": emb4["IA"].values,
        "MM": emb4["MM"].values,
        "MC": emb4["MC"].values,
        "total_out_degree": total_out,
        "participation": participation,
    }).set_index("strain")

    tree_tips = {leaf.name for leaf in tree.get_leaves()}
    common = sorted(set(traits.index) & tree_tips)
    if len(common) < 10:
        body = (
            f"Insufficient overlap: tree {len(tree_tips)}, matrices {len(traits)}, common {len(common)}. "
            "Check Newick tip labels match xlsx index."
        )
        replace_section(REPORT_PATH, "<!-- P2.8 RESULTS -->", "P2.8 Phylogenetic signal (PGLS)", body)
        return

    tree.prune(common)
    traits = traits.loc[common]

    # Brownian VCV: C[i,j] = root-to-MRCA(i,j) distance.
    n = len(common)
    name_to_idx = {name: i for i, name in enumerate(common)}
    C = np.zeros((n, n), dtype=float)
    dist_from_root = {leaf.name: tree.get_distance(leaf, topology_only=False) for leaf in tree.get_leaves()}
    for leaf_i in tree.get_leaves():
        for leaf_j in tree.get_leaves():
            mrca = tree.get_common_ancestor([leaf_i, leaf_j])
            d = tree.get_distance(mrca, topology_only=False)
            C[name_to_idx[leaf_i.name], name_to_idx[leaf_j.name]] = d

    def pagel_loglik(lam: float, y: np.ndarray) -> float:
        if lam < 0 or lam > 1:
            return 1e10
        diag = np.diag(np.diag(C))
        offd = C - diag
        Cl = diag + lam * offd
        try:
            sign, logdet = np.linalg.slogdet(Cl)
            if sign <= 0:
                return 1e10
            inv = np.linalg.inv(Cl)
        except np.linalg.LinAlgError:
            return 1e10
        ones = np.ones(n)
        beta = (ones @ inv @ y) / (ones @ inv @ ones)
        resid = y - beta
        sigma2 = (resid @ inv @ resid) / n
        if sigma2 <= 0:
            return 1e10
        ll = -0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)
        return -ll

    rows = []
    for trait in TRAITS_OF_INTEREST:
        y = traits[trait].values.astype(float)
        y = (y - y.mean()) / (y.std(ddof=0) if y.std(ddof=0) > 0 else 1.0)
        res = minimize_scalar(lambda l: pagel_loglik(l, y), bounds=(0, 1), method="bounded", options={"xatol": 1e-4})
        ll_ml = -res.fun
        lam_hat = res.x
        ll_null = -pagel_loglik(0.0, y)
        lr = 2 * (ll_ml - ll_null)
        p_val = float(chi2.sf(lr, df=1)) if lr > 0 else 1.0
        rows.append({"trait": trait, "pagel_lambda": float(lam_hat), "loglik": float(ll_ml), "lr_vs_null": float(lr), "p_value": p_val})

    df = pd.DataFrame(rows)
    df.to_csv(TABLE_DIR / "p2_8_pgls.csv", index=False)

    significant = df[df["p_value"] < 0.05]["trait"].tolist()
    body = f"""
Seed: {SEED}. PGLS, Pagel-λ ML, n={len(common)} strains.

{dataframe_to_md(df, index=False)}

- λ≈0 = independent of phylogeny; λ≈1 = heritable.
- p_value = LR vs λ=0.
- Significant: {', '.join(significant) if significant else 'none'}

Output: p2_8_pgls.csv
"""
    replace_section(REPORT_PATH, "<!-- P2.8 RESULTS -->", "P2.8 Phylogenetic signal (PGLS)", body)


if __name__ == "__main__":
    main()
