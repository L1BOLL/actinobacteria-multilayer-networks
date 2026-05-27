#!/usr/bin/env python3

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

# Below this density, greedy-modularity labels are noise → zero community switches.
MIN_DENSITY_FOR_MODULARITY = 0.05
K_RANGE = list(range(2, 9))

from source_paths import FROZEN_ROOT, find_layer_file

ROOT = FROZEN_ROOT
OUTDIR = ROOT / "figure_outputs"
OUTDIR.mkdir(parents=True, exist_ok=True)

LAYER_FILES = {
    "CCAM": "Changes in color of aerial mycelium (CCAM).xlsx",
    "CCVM": "Changes in color of vegetative mycelium (CCVM).xlsx",
    "CMP": "Changes in mycelium production (CMP).xlsx",
    "CS": "Changes in sporulation (CS).xlsx",
    "IG": "Inhibition of growth (IG).xlsx",
    "IC": "Invasion of colony (IC).xlsx",
    "RP": "Release of a pigment (RP).xlsx",
    "IRP": "Inhibition of release of pigment (IRP).xlsx",
    "IAC_RAC": "Induction of antimicrobial compounds of the other (IAC) release of antibacterial compounds (RAC).xlsx",
    "IAC_RDE": "Induction of antimicrobial compounds of the other (IAC) release of degrading enzyme (RDE).xlsx",
    "IRAC": "Inhibition of release of antimicrobial compounds (IRAC).xlsx",
    "RAC": "Release of antimicrobial compounds (RAC).xlsx",
}
CATEGORY_MAP = {
    "CCAM": "MC", "CCVM": "MC", "CMP": "MC", "CS": "MC",
    "IG": "DA", "IC": "DA", "RP": "MM", "IRP": "MM",
    "IAC_RAC": "IA", "IAC_RDE": "IA", "IRAC": "MM", "RAC": "DA",
}
CATEGORY_COLORS = {"DA": "#7f1d1d", "IA": "#b45309", "MM": "#0f766e", "MC": "#1d4ed8"}
LAYER_ORDER = list(LAYER_FILES)
STRAIN_ALIASES = {"S138 2": "MS138 2", "S625'": "S625", "S713'v": "S713'V", "S713'w": "S713'W"}


def set_style() -> None:
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Serif", "font.size": 10.5, "axes.titlesize": 16, "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False})
    sns.set_theme(style="whitegrid")


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    df.index = [STRAIN_ALIASES.get(str(x).strip(), str(x).strip()) for x in df.index]
    df.columns = [STRAIN_ALIASES.get(str(x).strip(), str(x).strip()) for x in df.columns]
    np.fill_diagonal(df.values, 0)
    return df


def load_all_matrices() -> dict[str, pd.DataFrame]:
    matrices = {layer: load_matrix(find_layer_file(fname)) for layer, fname in LAYER_FILES.items()}
    ref_nodes = matrices[LAYER_ORDER[0]].index.tolist()
    for layer in LAYER_ORDER:
        matrices[layer] = matrices[layer].reindex(index=ref_nodes, columns=ref_nodes).fillna(0)
    return matrices


def compute_layer_stats(binary_df: pd.DataFrame) -> dict[str, float]:
    arr = binary_df.values.astype(int)
    n = arr.shape[0]
    edges = int(arr.sum())
    density = edges / (n * (n - 1))
    reciprocal_pairs = 0
    connected_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            a = arr[i, j]
            b = arr[j, i]
            if a or b:
                connected_pairs += 1
                if a and b:
                    reciprocal_pairs += 1
    # R_p = mutual / connected_pairs (main); R_e = 2·mutual / edges (cross-ref).
    reciprocity = reciprocal_pairs / connected_pairs if connected_pairs else 0.0
    reciprocity_edges = (2 * reciprocal_pairs) / edges if edges else 0.0
    graph = nx.from_pandas_adjacency(binary_df.T, create_using=nx.DiGraph)
    largest_scc = max((len(c) for c in nx.strongly_connected_components(graph)), default=0)
    clustering = nx.average_clustering(graph.to_undirected()) if graph.number_of_nodes() else 0.0
    return {
        "edges": edges,
        "density": density,
        "reciprocity": reciprocity,
        "reciprocity_edges": reciprocity_edges,
        "largest_scc": largest_scc,
        "clustering": clustering,
    }


def build_stats_df(matrices: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for layer in LAYER_ORDER:
        bdf = (matrices[layer] > 0).astype(int)
        stats = compute_layer_stats(bdf)
        stats["layer"] = layer
        stats["category"] = CATEGORY_MAP[layer]
        rows.append(stats)
    return pd.DataFrame(rows).set_index("layer").loc[LAYER_ORDER]


def similarity_matrices(matrices: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    jaccard = pd.DataFrame(index=LAYER_ORDER, columns=LAYER_ORDER, dtype=float)
    cosine = pd.DataFrame(index=LAYER_ORDER, columns=LAYER_ORDER, dtype=float)
    vectors = {layer: (matrices[layer] > 0).astype(int).values.ravel() for layer in LAYER_ORDER}
    edge_sets = {layer: set(zip(*np.where((matrices[layer] > 0).astype(int).values > 0))) for layer in LAYER_ORDER}
    for a in LAYER_ORDER:
        for b in LAYER_ORDER:
            set_a = edge_sets[a]
            set_b = edge_sets[b]
            union = len(set_a | set_b)
            jaccard.loc[a, b] = len(set_a & set_b) / union if union else 1.0
            va = vectors[a]
            vb = vectors[b]
            denom = np.linalg.norm(va) * np.linalg.norm(vb)
            cosine.loc[a, b] = float(np.dot(va, vb) / denom) if denom else 0.0
    return jaccard, cosine


def _pick_k_by_cv(X: np.ndarray, k_range: list[int], seed: int = 0, n_splits: int = 10) -> int:
    """k by 10-fold CV held-out log-lik on diag-cov GMM. (At N=60, BIC/AIC over-pick.)"""
    from sklearn.model_selection import KFold

    scores = []
    for k in k_range:
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
        fold_scores = []
        for train_idx, test_idx in kf.split(X):
            if len(train_idx) < k:
                fold_scores.append(float("-inf"))
                continue
            try:
                gmm = GaussianMixture(n_components=k, covariance_type="diag", random_state=seed, n_init=50, reg_covar=1e-4)
                gmm.fit(X[train_idx])
                fold_scores.append(float(gmm.score(X[test_idx])))
            except Exception:
                fold_scores.append(float("-inf"))
        scores.append(float(np.mean(fold_scores)))
    return int(k_range[int(np.argmax(scores))])


def build_ecological_identity(matrices: dict[str, pd.DataFrame]) -> pd.DataFrame:
    tensor = np.stack([(matrices[layer].values > 0).astype(int) for layer in LAYER_ORDER], axis=0)
    L, N, _ = tensor.shape
    out_deg_layers = np.zeros((L, N), dtype=float)
    in_deg_layers = np.zeros((L, N), dtype=float)
    pagerank_layers = np.zeros((L, N), dtype=float)
    community_labels = np.full((L, N), -1, dtype=int)
    layer_densities = np.zeros(L, dtype=float)
    nodes = matrices[LAYER_ORDER[0]].index.tolist()
    for li, layer in enumerate(LAYER_ORDER):
        M = tensor[li]
        n = M.shape[0]
        layer_densities[li] = M.sum() / (n * (n - 1))
        G = nx.from_numpy_array(M.T, create_using=nx.DiGraph)
        out_deg_dict = dict(G.out_degree(weight=None))
        in_deg_dict = dict(G.in_degree(weight=None))
        for i in range(N):
            in_deg_layers[li, i] = in_deg_dict.get(i, 0.0)
            out_deg_layers[li, i] = out_deg_dict.get(i, 0.0)
        pr = nx.pagerank(G, alpha=0.85, weight=None)
        for i in range(N):
            pagerank_layers[li, i] = pr.get(i, 0.0)
        # Skip modularity on near-empty layers.
        if layer_densities[li] < MIN_DENSITY_FOR_MODULARITY:
            continue
        G_und = G.to_undirected()
        comms = list(nx.algorithms.community.greedy_modularity_communities(G_und))
        node_to_comm = {}
        for cid, com in enumerate(comms):
            for node_id in com:
                node_to_comm[node_id] = cid
        for i in range(N):
            community_labels[li, i] = node_to_comm.get(i, -1)
    participation = np.zeros(N, dtype=float)
    sd_outdeg = np.zeros(N, dtype=float)
    sd_pagerank = np.zeros(N, dtype=float)
    comm_switch_count = np.zeros(N, dtype=int)
    # Count switches over consecutive qualifying layers only.
    qualifying_layer_indices = np.where(layer_densities >= MIN_DENSITY_FOR_MODULARITY)[0]
    for i in range(N):
        k_out_vec = out_deg_layers[:, i]
        k_sum = k_out_vec.sum()
        participation[i] = 0.0 if k_sum == 0 else 1.0 - np.sum((k_out_vec / k_sum) ** 2)
        sd_outdeg[i] = np.std(k_out_vec)
        sd_pagerank[i] = np.std(pagerank_layers[:, i])
        labels = community_labels[qualifying_layer_indices, i]
        switches = 0
        for lidx in range(1, len(labels)):
            if labels[lidx] != labels[lidx - 1] and labels[lidx] != -1 and labels[lidx - 1] != -1:
                switches += 1
        comm_switch_count[i] = switches
    cols = [f"outdeg_{layer}" for layer in LAYER_ORDER] + [f"indeg_{layer}" for layer in LAYER_ORDER]
    E_core = np.hstack([out_deg_layers.T, in_deg_layers.T])
    extra_features = np.vstack([participation, sd_outdeg, sd_pagerank, comm_switch_count]).T
    extra_cols = ["participation", "sd_outdeg", "sd_pagerank", "community_switches"]
    ecological_identity = pd.DataFrame(np.hstack([E_core, extra_features]), index=nodes, columns=cols + extra_cols)
    scaler = StandardScaler()
    E_scaled = scaler.fit_transform(ecological_identity.values)
    n_pcs = min(20, E_scaled.shape[1], E_scaled.shape[0])
    pca = PCA(n_components=n_pcs, random_state=0)
    E_pcs = pca.fit_transform(E_scaled)
    for pc_idx in range(n_pcs):
        ecological_identity[f"PC{pc_idx+1}"] = E_pcs[:, pc_idx]
    # k by 10-fold CV log-lik on diag-cov GMM.
    k_chosen = _pick_k_by_cv(E_pcs, K_RANGE, seed=0)
    cluster_labels = KMeans(n_clusters=k_chosen, random_state=0, n_init=50).fit_predict(E_pcs)
    ecological_identity["cluster"] = cluster_labels
    ecological_identity.attrs["k_chosen"] = k_chosen
    ecological_identity.attrs["layer_densities"] = dict(zip(LAYER_ORDER, layer_densities.tolist()))
    ecological_identity.attrs["modularity_excluded_layers"] = [
        layer for layer, d in zip(LAYER_ORDER, layer_densities) if d < MIN_DENSITY_FOR_MODULARITY
    ]
    ecological_identity.index.name = "strain"
    ecological_identity.to_csv(ROOT / "ecological_identity_with_PCs_and_clusters.csv")
    return ecological_identity.reset_index()


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUTDIR / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUTDIR / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.1, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom", ha="right")


def make_figure_2(matrices: dict[str, pd.DataFrame], stats_df: pd.DataFrame) -> None:
    selected = ["RAC", "IG", "CCVM", "IRP"]
    fig = plt.figure(figsize=(16, 9.9))
    gs = GridSpec(2, 8, figure=fig, height_ratios=[1.0, 0.9], hspace=0.34, wspace=0.45)
    cmap = colors.LinearSegmentedColormap.from_list("paper_heat", ["#fbfdff", "#111827"])
    for idx, layer in enumerate(selected):
        ax = fig.add_subplot(gs[0, idx * 2 : (idx + 1) * 2])
        bdf = (matrices[layer] > 0).astype(int)
        sns.heatmap(bdf, ax=ax, cmap=cmap, cbar=False, square=True, xticklabels=False, yticklabels=False, linewidths=0)
        ax.text(0.5, 1.12, f"{layer}  [{CATEGORY_MAP[layer]}]", transform=ax.transAxes, ha="center", va="bottom", fontsize=15, color=CATEGORY_COLORS[CATEGORY_MAP[layer]])
        ax.text(0.5, 1.05, f"density {stats_df.loc[layer, 'density']:.3f}   reciprocity {stats_df.loc[layer, 'reciprocity']:.3f}", transform=ax.transAxes, ha="center", va="bottom", fontsize=9, color="#475569")
    metrics = ["density", "reciprocity", "largest_scc", "clustering"]
    metric_titles = ["Edge density", "Mutual dyad fraction", "Largest strongly connected core", "Undirected clustering"]
    for midx, (metric, title) in enumerate(zip(metrics, metric_titles)):
        ax = fig.add_subplot(gs[1, midx * 2 : (midx + 1) * 2])
        palette = [CATEGORY_COLORS[CATEGORY_MAP[layer]] for layer in LAYER_ORDER]
        ax.bar(np.arange(len(LAYER_ORDER)), stats_df[metric], color=palette, edgecolor="white")
        ax.set_xticks(np.arange(len(LAYER_ORDER)))
        ax.set_xticklabels(LAYER_ORDER, rotation=45, ha="right")
        ax.set_title(title)
    save_figure(fig, "figure2_antagonistic_backbone")


def make_figure_3(matrices: dict[str, pd.DataFrame]) -> None:
    stats_df = build_stats_df(matrices)
    jaccard, cosine = similarity_matrices(matrices)
    dist = 1 - jaccard
    np.fill_diagonal(dist.values, 0)
    condensed = dist.values[np.triu_indices_from(dist.values, k=1)]
    leaves = leaves_list(linkage(condensed, method="average")) if len(condensed) else np.arange(len(LAYER_ORDER))
    ordered = [LAYER_ORDER[i] for i in leaves]
    ordered_jaccard = jaccard.loc[ordered, ordered]
    ordered_cosine = cosine.loc[ordered, ordered]
    fig = plt.figure(figsize=(15.5, 9.3))
    gs = GridSpec(2, 4, figure=fig, width_ratios=[1.05, 1.05, 1.05, 0.85], height_ratios=[0.9, 1.1], wspace=0.35, hspace=0.35)
    bar_ax = fig.add_subplot(gs[0, :3])
    palette = [CATEGORY_COLORS[CATEGORY_MAP[layer]] for layer in LAYER_ORDER]
    bar_ax.bar(np.arange(len(LAYER_ORDER)), stats_df['density'], color=palette, edgecolor='white')
    bar_ax.set_xticks(np.arange(len(LAYER_ORDER)))
    bar_ax.set_xticklabels(LAYER_ORDER, rotation=45, ha='right')
    bar_ax.set_title('Layer densities')
    jac_ax = fig.add_subplot(gs[1, 0:2])
    sns.heatmap(ordered_jaccard, ax=jac_ax, cmap=sns.color_palette('mako', as_cmap=True), vmin=0, vmax=ordered_jaccard.values.max(), cbar=False, square=True, linewidths=0.5, linecolor='white')
    cos_ax = fig.add_subplot(gs[1, 2])
    sns.heatmap(ordered_cosine, ax=cos_ax, cmap=sns.color_palette('crest', as_cmap=True), vmin=0, vmax=max(0.55, ordered_cosine.values.max()), cbar=False, square=True, linewidths=0.5, linecolor='white')
    table_ax = fig.add_subplot(gs[1, 3]); table_ax.axis('off')
    top_pairs = []
    for i, a in enumerate(LAYER_ORDER):
        for b in LAYER_ORDER[i + 1:]:
            top_pairs.append((a, b, jaccard.loc[a, b], cosine.loc[a, b]))
    top_pairs.sort(key=lambda row: row[2], reverse=True)
    y = 0.97
    for rank, (a, b, jac, cosv) in enumerate(top_pairs[:7], start=1):
        table_ax.text(0.0, y, f"{rank}. {a} - {b}", fontsize=10.2, fontweight='bold'); table_ax.text(0.0, y - 0.05, f"Jaccard {jac:.3f}", fontsize=9.1); table_ax.text(0.53, y - 0.05, f"Cosine {cosv:.3f}", fontsize=9.1); y -= 0.115
    save_figure(fig, 'figure3_layer_similarity_partitioning')


def make_figure_4(embeddings: pd.DataFrame, identity: pd.DataFrame) -> None:
    df = embeddings.merge(identity[["strain", "participation", "cluster"]], on="strain", how="left")
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(df[["DA", "IA", "MM", "MC"]])
    df["pc1"] = coords[:, 0]; df["pc2"] = coords[:, 1]
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    exemplar_names = set(); exemplar_names.update(df.nlargest(2, 'DA')['strain']); exemplar_names.update(df.nlargest(2, 'IA')['strain']); exemplar_names.update(df.nlargest(2, 'MM')['strain']); exemplar_names.update(df.nlargest(2, 'MC')['strain'])
    exemplars = df[df['strain'].isin(exemplar_names)].copy()
    fig = plt.figure(figsize=(15.9, 9.3)); gs = GridSpec(1, 3, figure=fig, width_ratios=[1.42, 0.74, 0.84], wspace=0.38)
    scatter_ax = fig.add_subplot(gs[0, 0])
    palette = sns.color_palette('Spectral', n_colors=int(df['cluster'].nunique()))
    cluster_colors = {cluster: palette[i] for i, cluster in enumerate(sorted(df['cluster'].dropna().unique()))}
    for cluster, sub in df.groupby('cluster'):
        scatter_ax.scatter(sub['pc1'], sub['pc2'], s=45 + 120 * sub['participation'], color=cluster_colors.get(cluster, '#64748b'), edgecolor='white', linewidth=0.7, alpha=0.9, label=f"Cluster {int(cluster)}")
    for i, axis_name in enumerate(['DA', 'IA', 'MM', 'MC']):
        x, y = loadings[i, 0] * 3.0, loadings[i, 1] * 3.0
        scatter_ax.add_patch(FancyArrowPatch((0, 0), (x, y), arrowstyle='-|>', mutation_scale=13, linewidth=1.8, color=CATEGORY_COLORS.get(axis_name, CATEGORY_COLORS['MM'])))
        scatter_ax.text(x * 1.08, y * 1.08, axis_name, fontsize=11, fontweight='bold')
    for _, row in exemplars.iterrows():
        scatter_ax.text(row['pc1'] + 0.08, row['pc2'] + 0.05, row['strain'], fontsize=9)
    profile_ax = fig.add_subplot(gs[0, 1])
    ordered_profiles = exemplars[['strain', 'DA', 'IA', 'MM', 'MC']].assign(total=lambda x: x[['DA', 'IA', 'MM', 'MC']].sum(axis=1)).sort_values('total', ascending=False).head(8).set_index('strain')[['DA', 'IA', 'MM', 'MC']]
    sns.heatmap(ordered_profiles, ax=profile_ax, cmap=sns.color_palette('rocket', as_cmap=True), linewidths=0.5, linecolor='white', cbar_kws={'label': 'Scaled ecological score'})
    summary_ax = fig.add_subplot(gs[0, 2]); summary_ax.axis('off')
    blocks = [('DA leaders', df.nlargest(3, 'DA')[['strain', 'DA']].values.tolist(), CATEGORY_COLORS['DA']), ('IA leaders', df.nlargest(3, 'IA')[['strain', 'IA']].values.tolist(), CATEGORY_COLORS['IA']), ('MM leaders', df.nlargest(3, 'MM')[['strain', 'MM']].values.tolist(), CATEGORY_COLORS['MM']), ('MC leaders', df.nlargest(3, 'MC')[['strain', 'MC']].values.tolist(), CATEGORY_COLORS['MC'])]
    y = 0.95
    for title, rows, color in blocks:
        summary_ax.text(0.0, y, title, fontsize=10.5, fontweight='bold', color=color); y -= 0.06
        for strain, score in rows:
            summary_ax.text(0.02, y, f"{strain}: {score:.2f}", fontsize=9.5); y -= 0.05
        y -= 0.03
    save_figure(fig, 'figure4_ecological_strategy_space')


def make_figure_5(matrices: dict[str, pd.DataFrame], identity: pd.DataFrame) -> None:
    layer_rows = []
    for layer in LAYER_ORDER:
        bdf = (matrices[layer] > 0).astype(int)
        layer_rows.append(pd.DataFrame({'strain': bdf.columns, 'layer': layer, 'outgoing': bdf.sum(axis=0).values, 'incoming': bdf.sum(axis=1).values}))
    long_df = pd.concat(layer_rows, ignore_index=True)
    long_df['asymmetry'] = long_df['outgoing'] - long_df['incoming']
    pivot = long_df.pivot(index='strain', columns='layer', values='asymmetry').loc[:, LAYER_ORDER]
    zscores = pivot.sub(pivot.mean(axis=1), axis=0).div(pivot.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    # Role-variability = sum of z-scored {participation, sd_outdeg, community_switches}.
    rv_features = identity[['strain', 'participation', 'sd_outdeg', 'community_switches']].copy()
    rv_z = rv_features[['participation', 'sd_outdeg', 'community_switches']].apply(
        lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else s * 0.0
    )
    rv_features['role_variability'] = rv_z.sum(axis=1)
    variability = rv_features.sort_values('role_variability', ascending=False)
    selected = variability.head(18)['strain'].tolist(); heat_df = zscores.loc[selected]
    fig = plt.figure(figsize=(15.5, 9.8)); gs = GridSpec(1, 3, figure=fig, width_ratios=[1.15, 0.95, 0.7], wspace=0.25)
    heat_ax = fig.add_subplot(gs[0, 0]); sns.heatmap(heat_df, ax=heat_ax, cmap=sns.diverging_palette(230, 15, as_cmap=True), center=0, linewidths=0.5, linecolor='white', cbar_kws={'label': 'Within-strain asymmetry z-score'})
    scatter_ax = fig.add_subplot(gs[0, 1]); scatter = scatter_ax.scatter(identity['sd_outdeg'], identity['community_switches'], c=identity['participation'], cmap=sns.color_palette('viridis', as_cmap=True), s=55, edgecolor='white', linewidth=0.7); fig.colorbar(scatter, ax=scatter_ax, fraction=0.046, pad=0.04)
    summary_ax = fig.add_subplot(gs[0, 2]); summary_ax.axis('off'); y = 0.96
    for _, row in variability.head(10).iterrows():
        summary_ax.text(0.0, y, f"{row['strain']}", fontsize=10, fontweight='bold'); summary_ax.text(0.0, y - 0.05, f"participation={row['participation']:.2f}  sd_out={row['sd_outdeg']:.2f}  switches={row['community_switches']:.0f}", fontsize=8.8); y -= 0.10
    save_figure(fig, 'figure5_role_switching')


def make_embedding_representation_figure(embeddings: pd.DataFrame, identity: pd.DataFrame) -> None:
    df = embeddings.merge(identity[['strain', 'participation', 'cluster']], on='strain', how='left')
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(df[['DA', 'IA', 'MM', 'MC']])
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    df['pc1'] = coords[:, 0]; df['pc2'] = coords[:, 1]
    archetypes = [('B201', 'DA-dominant'), ('MS10 14', 'IA-dominant'), ('S1705', 'MM-dominant'), ('S1350', 'MC-dominant'), ('MS53 7', 'High-total generalist'), ('S625', 'DA + IA hybrid')]
    selected = df[df['strain'].isin([name for name, _ in archetypes])].copy()
    highlight_colors = {'B201': CATEGORY_COLORS['DA'], 'MS10 14': CATEGORY_COLORS['IA'], 'S1705': CATEGORY_COLORS['MM'], 'S1350': CATEGORY_COLORS['MC'], 'MS53 7': '#111827', 'S625': '#8b5cf6'}
    fig = plt.figure(figsize=(16.2, 10.0)); gs = GridSpec(2, 4, figure=fig, width_ratios=[1.55, 1.0, 1.0, 1.0], height_ratios=[1.0, 1.0], wspace=0.5, hspace=0.38)
    scatter_ax = fig.add_subplot(gs[:, 0]); scatter_ax.scatter(df['pc1'], df['pc2'], s=50, color='#d9dde3', edgecolor='white', linewidth=0.6, alpha=0.95, zorder=1)
    for _, row in selected.iterrows():
        scatter_ax.scatter(row['pc1'], row['pc2'], s=120, color=highlight_colors[row['strain']], edgecolor='white', linewidth=1.0, zorder=3); scatter_ax.text(row['pc1'] + 0.12, row['pc2'] + 0.08, row['strain'], fontsize=9.5, color='#111827')
    for i, axis_name in enumerate(['DA', 'IA', 'MM', 'MC']):
        x, y = loadings[i, 0] * 3.2, loadings[i, 1] * 3.2
        scatter_ax.add_patch(FancyArrowPatch((0, 0), (x, y), arrowstyle='-|>', mutation_scale=14, linewidth=1.8, color=CATEGORY_COLORS.get(axis_name, CATEGORY_COLORS['MM'])))
        scatter_ax.text(x * 1.08, y * 1.08, axis_name, fontsize=11, fontweight='bold')
    categories = ['DA', 'IA', 'MM', 'MC']; angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False); angles = np.concatenate([angles, [angles[0]]])
    for idx, (strain, label) in enumerate(archetypes):
        ax = fig.add_subplot(gs[idx // 3, 1 + (idx % 3)], projection='polar')
        row = selected[selected['strain'] == strain].iloc[0]
        values = row[categories].to_numpy(dtype=float); values = np.concatenate([values, [values[0]]])
        color = highlight_colors[strain]
        ax.set_theta_offset(np.pi / 2); ax.set_theta_direction(-1); ax.plot(angles, values, color=color, linewidth=2.0); ax.fill(angles, values, color=color, alpha=0.18); ax.set_ylim(0, 10); ax.set_xticks(angles[:-1]); ax.set_xticklabels(categories, fontsize=9); ax.set_title(f"{strain}\n{label}", va='bottom', pad=14, fontsize=10.5, color=color)
    save_figure(fig, 'embedding_representation_archetypes')


def write_manifest(stats_df: pd.DataFrame) -> None:
    lines = ['Generated figures:', '  figure2_antagonistic_backbone.{png,svg}', '  figure3_layer_similarity_partitioning.{png,svg}', '  figure4_ecological_strategy_space.{png,svg}', '  figure5_role_switching.{png,svg}', '  embedding_representation_archetypes.{png,svg}', '', 'Layer statistics used in the figures:', stats_df.round(4).to_csv(index=True)]
    (OUTDIR / 'manifest.txt').write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    set_style()
    matrices = load_all_matrices()
    stats_df = build_stats_df(matrices)
    embeddings = pd.read_csv(ROOT / 'ecological_embeddings_4D.csv').rename(columns={'Unnamed: 0': 'strain'}) if (ROOT / 'ecological_embeddings_4D.csv').exists() else None
    if embeddings is None:
        raise FileNotFoundError(ROOT / 'ecological_embeddings_4D.csv')
    identity_path = ROOT / 'ecological_identity_with_PCs_and_clusters.csv'
    identity = pd.read_csv(identity_path).rename(columns={'Unnamed: 0': 'strain'}) if identity_path.exists() else build_ecological_identity(matrices)
    make_figure_2(matrices, stats_df)
    make_figure_3(matrices)
    make_figure_4(embeddings, identity)
    make_figure_5(matrices, identity)
    make_embedding_representation_figure(embeddings, identity)
    write_manifest(stats_df)
    print('Generated figures in', OUTDIR)


if __name__ == '__main__':
    main()
