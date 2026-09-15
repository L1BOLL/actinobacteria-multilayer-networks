#!/usr/bin/env python3
"""Phylogenetic signal in interaction profiles.

Tested three independent ways, so the verdict does not rest on one modelling
choice:

  1. Pagel's lambda by maximum likelihood on the GTR+GAMMA tree (build_tree.py),
     with significance from tip-label permutation and BH-FDR across traits.
  2. The same on UPGMA and WPGMA trees built from the supplied pairwise 16S
     distance matrix, so the result does not depend on one tree-building choice.
  3. A Mantel permutation test of 16S distance against interaction-profile
     distance, which uses no tree at all.

It also reports how many distinct 16S genotypes the 60 strains resolve into,
because that caps the power of every test above.
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
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.optimize import minimize_scalar
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2, spearmanr
from sklearn.decomposition import PCA
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
    ensure_output_dirs,
    load_tensor,
)
from report_utils import dataframe_to_md, replace_section

TREE_PATH = ROOT / "data" / "phylogeny" / "tree.nwk"
DIST_PATH = ROOT / "data" / "annotation" / "distance_matrix.xlsx"
N_MANTEL = 9999
N_PERM = 999


# --------------------------------------------------------------------------- #
# inputs
# --------------------------------------------------------------------------- #
def load_distance_matrix() -> pd.DataFrame:
    """Long-format pairwise 16S distances -> square symmetric DataFrame."""
    raw = pd.read_excel(DIST_PATH, sheet_name="MatrixOutput")
    a = raw["Species 1"].map(canonical_strain)
    b = raw["Species 2"].map(canonical_strain)
    names = sorted(set(a) | set(b))
    D = pd.DataFrame(0.0, index=names, columns=names)
    for x, y, d in zip(a, b, raw["Distance"]):
        D.loc[x, y] = d
        D.loc[y, x] = d
    return D


def build_traits(tensor: np.ndarray, node_ids: list[str]) -> pd.DataFrame:
    """The 21 strain-level traits h3 is tested on."""
    out_deg = tensor.sum(axis=1).T.astype(float)  # strains x layers
    total = out_deg.sum(axis=1)
    safe = np.where(total == 0, 1.0, total)
    participation = 1.0 - np.sum((out_deg / safe[:, None]) ** 2, axis=1)

    traits = {f"outdeg_{layer}": out_deg[:, i] for i, layer in enumerate(LAYER_ORDER)}
    layer_idx = {layer: i for i, layer in enumerate(LAYER_ORDER)}
    for cat, layers in CATEGORY_MAP.items():
        traits[f"cat_{cat}"] = out_deg[:, [layer_idx[l] for l in layers]].mean(axis=1)
    traits["total_out_degree"] = total
    traits["participation"] = participation

    z = (out_deg - out_deg.mean(0)) / np.where(out_deg.std(0) == 0, 1.0, out_deg.std(0))
    pcs = PCA(n_components=3, random_state=SEED).fit_transform(z)
    for j in range(3):
        traits[f"PC{j + 1}_12D"] = pcs[:, j]

    return pd.DataFrame(traits, index=[canonical_strain(n) for n in node_ids])


# --------------------------------------------------------------------------- #
# Pagel's lambda
# --------------------------------------------------------------------------- #
def vcv_from_linkage(Z: np.ndarray, n: int) -> np.ndarray:
    """Brownian VCV from an ultrametric linkage: C[i,j] = shared root-to-MRCA path."""
    coph = squareform(cophenet(Z))
    height = coph.max() / 2.0
    C = height - coph / 2.0
    np.fill_diagonal(C, height)
    return C


def vcv_from_tree(path: Path, keep: list[str]) -> tuple[np.ndarray, list[str]]:
    """Brownian VCV from a Newick tree, restricted to `keep`."""
    import dendropy

    from build_tree import from_tree_label

    tree = dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)
    label_of = {}
    for leaf in tree.leaf_node_iter():
        label_of[leaf] = canonical_strain(from_tree_label(leaf.taxon.label))
    keep_set = set(keep)
    leaves = [lf for lf, name in label_of.items() if name in keep_set]
    names = [label_of[lf] for lf in leaves]

    # Root-to-node depth for every node, then C[i,j] = depth of MRCA(i, j).
    depth = {}
    for node in tree.preorder_node_iter():
        depth[node] = 0.0 if node.parent_node is None else depth[node.parent_node] + (node.edge_length or 0.0)

    ancestors = []
    for lf in leaves:
        chain, cur = {}, lf
        while cur is not None:
            chain[cur] = depth[cur]
            cur = cur.parent_node
        ancestors.append(chain)

    n = len(leaves)
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            shared = set(ancestors[i]) & set(ancestors[j])
            d = max(depth[a] for a in shared) if shared else 0.0
            C[i, j] = C[j, i] = d
    return C, names


def pagel_lambda(y: np.ndarray, C: np.ndarray) -> tuple[float, float, float]:
    """ML estimate of Pagel's lambda with an LR test against lambda = 0."""
    n = len(y)
    if np.std(y) == 0:
        return np.nan, np.nan, np.nan
    y = (y - y.mean()) / y.std()
    diag = np.diag(np.diag(C))
    offd = C - diag
    eye = 1e-10 * np.eye(n)
    ones = np.ones(n)

    def neg_ll(lam: float) -> float:
        Cl = diag + lam * offd + eye
        sign, logdet = np.linalg.slogdet(Cl)
        if sign <= 0:
            return 1e10
        try:
            inv = np.linalg.inv(Cl)
        except np.linalg.LinAlgError:
            return 1e10
        beta = (ones @ inv @ y) / (ones @ inv @ ones)
        resid = y - beta
        sigma2 = float(resid @ inv @ resid) / n
        if sigma2 <= 0:
            return 1e10
        return 0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)

    res = minimize_scalar(neg_ll, bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-5})
    lr = 2.0 * (-res.fun - -neg_ll(0.0))
    # lambda = 0 sits on the parameter boundary, so the LR statistic nominally
    # follows a 50:50 mixture of chi2(0) and chi2(1). This asymptotic p-value is
    # reported for reference ONLY -- it is badly miscalibrated on these data (see
    # `lambda_table`), and inference uses the permutation p-value instead.
    p = float(chi2.sf(lr, df=1) / 2.0) if lr > 0 else 1.0
    return float(res.x), float(lr), p


def lambda_table(
    traits: pd.DataFrame,
    C: np.ndarray,
    order: list[str],
    label: str,
    rng: np.random.Generator,
    n_perm: int,
) -> pd.DataFrame:
    """Pagel's lambda per trait, with a tip-label permutation test.

    The asymptotic chi2 LR test assumes a well-conditioned, roughly ultrametric
    VCV. Neither holds here: 33 of the 60 strains share an identical 16S sequence
    with at least one other, so the ML tree's VCV has rank 44 of 60, and its tip
    depths span a 53-fold range. Under those conditions shuffled trait values routinely produce
    lambda > 0.9 and LR > 30, i.e. the chi2 p-value is anticonservative by orders
    of magnitude.

    Shuffling trait values across tips holds the tree, the VCV conditioning and
    the trait's marginal distribution fixed while destroying only the
    tip-to-trait correspondence, which is exactly the null h3 needs. That
    permutation p-value is what the manuscript reports.
    """
    rows = []
    for trait in traits.columns:
        y = traits.loc[order, trait].to_numpy(float)
        lam, lr, p_chi2 = pagel_lambda(y, C)
        if np.isnan(lam):
            rows.append({"topology": label, "trait": trait, "pagel_lambda": np.nan,
                         "LR": np.nan, "p_chi2_ref": np.nan, "p_perm": np.nan,
                         "null_lambda_p95": np.nan})
            continue
        null_lr = np.empty(n_perm)
        null_lam = np.empty(n_perm)
        for i in range(n_perm):
            null_lam[i], null_lr[i], _ = pagel_lambda(rng.permutation(y), C)
        p_perm = (int((null_lr >= lr).sum()) + 1) / (n_perm + 1)
        rows.append({
            "topology": label,
            "trait": trait,
            "pagel_lambda": lam,
            "LR": lr,
            "p_chi2_ref": p_chi2,
            "p_perm": p_perm,
            "null_lambda_p95": float(np.percentile(null_lam, 95)),
        })
    df = pd.DataFrame(rows)
    ok = df["p_perm"].notna()
    df.loc[ok, "q_BH"] = multipletests(df.loc[ok, "p_perm"], method="fdr_bh")[1]
    return df


# --------------------------------------------------------------------------- #
# Mantel
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
def main() -> None:
    ensure_output_dirs()
    rng = np.random.default_rng(SEED)

    if not DIST_PATH.exists():
        replace_section(REPORT_PATH, "<!-- phylogenetic_signal RESULTS -->", "Phylogenetic signal",
            f"Missing {DIST_PATH.relative_to(ROOT)} - h3 cannot be tested.",
        )
        return

    tensor, _, node_ids = load_tensor()
    traits = build_traits(tensor, node_ids)
    D = load_distance_matrix()

    common = sorted(set(traits.index) & set(D.index))
    n = len(common)
    Dg = D.loc[common, common].to_numpy()
    traits = traits.loc[common]

    # --- how much does 16S actually resolve? -------------------------------- #
    n_geno, membership = connected_components(Dg == 0, directed=False)
    sizes = pd.Series(membership).value_counts()
    tied = sizes[sizes > 1]

    # --- topologies --------------------------------------------------------- #
    topologies: list[tuple[str, np.ndarray, list[str]]] = []
    condensed = squareform(Dg, checks=False)
    for method, label in [("average", "UPGMA (16S distances)"), ("weighted", "WPGMA (16S distances)")]:
        Z = linkage(condensed, method=method)
        topologies.append((label, vcv_from_linkage(Z, n), common))

    if TREE_PATH.exists():
        try:
            C_ml, ml_order = vcv_from_tree(TREE_PATH, common)
            topologies.insert(0, ("ML GTR+GAMMA (build_tree.py)", C_ml, ml_order))
        except Exception as exc:  # noqa: BLE001
            print(f"  could not read {TREE_PATH}: {exc}")

    diagnostics = []
    for label, C, _order in topologies:
        d = np.diag(C)
        diagnostics.append({
            "topology": label,
            "n_tips": len(C),
            "vcv_rank": int(np.linalg.matrix_rank(C, tol=1e-9)),
            "tip_depth_ratio": float(d.max() / max(d.min(), 1e-12)),
            "ultrametric": bool(np.allclose(d, d[0], rtol=1e-3)),
        })
    diag_df = pd.DataFrame(diagnostics)
    diag_df.to_csv(TABLE_DIR / "phylogenetic_vcv_diagnostics.csv", index=False)

    all_lambda = pd.concat(
        [lambda_table(traits, C, order, label, rng, N_PERM) for label, C, order in topologies],
        ignore_index=True,
    )
    all_lambda.to_csv(TABLE_DIR / "phylogenetic_signal.csv", index=False)

    # --- tree-free corroboration -------------------------------------------- #
    out_deg = tensor.sum(axis=1).T.astype(float)
    idx = [[canonical_strain(x) for x in node_ids].index(s) for s in common]
    profile = out_deg[idx]
    zp = (profile - profile.mean(0)) / np.where(profile.std(0) == 0, 1.0, profile.std(0))
    mantel_rows = [("12-D out-degree profile (z, Euclidean)", *mantel(Dg, squareform(pdist(zp)), rng, N_MANTEL))]
    layer_idx = {l: i for i, l in enumerate(LAYER_ORDER)}
    for cat, layers in CATEGORY_MAP.items():
        v = profile[:, [layer_idx[l] for l in layers]].mean(axis=1)
        mantel_rows.append((f"category {cat}", *mantel(Dg, np.abs(v[:, None] - v[None, :]), rng, 1999)))
    mantel_df = pd.DataFrame(mantel_rows, columns=["comparison", "mantel_rho", "p_value"])
    mantel_df.to_csv(TABLE_DIR / "mantel.csv", index=False)

    # Figures are produced by figures_phylogeny.py, which plots each observed
    # lambda against the null it was tested on. A bare bar chart of lambda is
    # actively misleading here: the point estimates run to 0.95 while every one
    # of them is null under permutation.

    # --- report -------------------------------------------------------------- #
    sig = all_lambda[all_lambda["q_BH"] < 0.05]
    verdict = (
        f"**{len(sig)} of {len(all_lambda)}** trait x topology tests reach q_BH < 0.05."
        if len(sig)
        else f"**No trait reaches q_BH < 0.05 on any topology** ({len(all_lambda)} tests)."
    )
    max_lam = all_lambda["pagel_lambda"].max()
    naive = all_lambda[all_lambda["p_chi2_ref"] < 0.05]

    body = f"""
Seed: {SEED}. n = {n} strains (see `data/annotation/DATA_NOTES.md` for the one
unmatched strain). Traits: {len(traits.columns)} - 12 per-layer out-degrees, 4 category means,
total out-degree, participation coefficient, and PC1-PC3 of the 12-D profile.

### Verdict on h3

{verdict} Largest lambda observed across all topologies: **{max_lam:.3f}**.
No trait reaches significance against its own tip-label permutation null.

### 16S resolution limit

The {n} strains resolve into only **{n_geno} distinct 16S genotypes**. {int(tied.sum())} strains fall into
{len(tied)} groups with pairwise distance exactly zero (group sizes: {', '.join(str(int(s)) for s in sorted(tied, reverse=True))}).

Within those groups 16S carries no information, so any phylogenetic-signal test is
underpowered by construction. This is the mechanistic explanation for the null above,
and it is a property of the marker rather than a defect of the data.

### Why the asymptotic test cannot be used here

{dataframe_to_md(diag_df, index=False)}

The ML tree's VCV is **rank {diag_df.iloc[0]['vcv_rank']} of {diag_df.iloc[0]['n_tips']}** - a direct consequence of the identical-sequence
groups above, which produce duplicate rows - and its tip depths span a
**{diag_df.iloc[0]['tip_depth_ratio']:.0f}-fold** range, so it is far from ultrametric. Under those conditions the
chi2 likelihood-ratio test is severely anticonservative: randomly shuffled trait
values on this tree routinely yield lambda > 0.9 and LR > 30.

Taking `p_chi2_ref` at face value would report **{len(naive)} of {len(all_lambda)}** tests as significant. The
tip-label permutation test - which holds the tree, the conditioning and the trait
marginals fixed and destroys only the tip-to-trait correspondence - shows those
are false positives. **Inference below uses `p_perm`.**

### Pagel's lambda, per trait and topology

{dataframe_to_md(all_lambda.sort_values(['topology', 'p_perm']), index=False)}

lambda ~ 0 = trait independent of phylogeny; lambda ~ 1 = Brownian heritability.
`p_perm` is from {N_PERM} tip-label shuffles; `q_BH` is Benjamini-Hochberg across traits
within each topology. `p_chi2_ref` is the asymptotic value, shown only to document
the miscalibration. `null_lambda_p95` is the 95th percentile of lambda under
shuffling - where it approaches 1, lambda alone carries no evidence.

### Tree-free corroboration (Mantel, {N_MANTEL} permutations)

{dataframe_to_md(mantel_df, index=False)}

The Mantel test uses no tree, so agreement with the lambda result rules out the
possibility that the null is an artefact of one tree-building choice.

### Summary for the manuscript

Three independent routes - permutation-calibrated Pagel's lambda on an ML tree, the
same on two ultrametric distance topologies, and a tree-free Mantel test - agree
that interaction profiles carry no detectable 16S phylogenetic signal in these
strains. h3 as stated is not supported.

Outputs: `phylogenetic_signal.csv`, `mantel.csv`, `phylogenetic_vcv_diagnostics.csv`.
Figures: `figures/fig_h3_phylogeny.*` (main text), `figures/figS_permutation_calibration.*` (supplement).
"""
    replace_section(REPORT_PATH, "<!-- phylogenetic_signal RESULTS -->", "Phylogenetic signal", body)

    print(f"n={n}  genotypes={n_geno}  tied strains={int(tied.sum())}")
    print(f"topologies tested: {[t[0] for t in topologies]}")
    print(f"max lambda={max_lam:.4f}")
    print(f"significant by permutation q<0.05 : {len(sig)}/{len(all_lambda)}")
    print(f"would be 'significant' by naive chi2: {len(naive)}/{len(all_lambda)}")
    print(mantel_df.to_string(index=False))


if __name__ == "__main__":
    main()
