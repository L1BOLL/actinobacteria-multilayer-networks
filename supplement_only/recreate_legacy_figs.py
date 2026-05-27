#!/usr/bin/env python3

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTDIR = Path(__file__).resolve().parent

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

LAYER_ORDER = list(LAYER_FILES)
GROUP_MAP = {
    "CCAM": "MC",
    "CCVM": "MC",
    "CMP": "MC",
    "CS": "MC",
    "IG": "DA",
    "IC": "DA",
    "RP": "MM",
    "IRP": "MM",
    "IAC_RAC": "IA",
    "IAC_RDE": "IA",
    "IRAC": "MM",
    "RAC": "DA",
}
GROUP_COLORS = {
    "DA": "#0072B2",
    "IA": "#D55E00",
    "MM": "#009E73",
    "MC": "#CC79A7",
}


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.titlesize": 12,
            "axes.labelsize": 9,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "grid.color": "#d9dde3",
            "grid.linewidth": 0.5,
            "grid.alpha": 1.0,
            "savefig.dpi": 400,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUTDIR / f"{stem}.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTDIR / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_layer(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    np.fill_diagonal(df.values, 0)
    return df


def load_layers() -> dict[str, pd.DataFrame]:
    layers = {name: load_layer(ROOT / fname) for name, fname in LAYER_FILES.items()}
    ref = layers[LAYER_ORDER[0]].index.tolist()
    for name in LAYER_ORDER:
        layers[name] = layers[name].reindex(index=ref, columns=ref)
    return layers


def compute_layer_summary(layers: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for layer_id in LAYER_ORDER:
        A = (layers[layer_id].values > 0).astype(int)
        n = A.shape[0]
        edge_count = int(A.sum())
        density = edge_count / (n * (n - 1))
        reciprocity_count = 0
        connected_pairs = 0
        for i in range(n):
            for j in range(i + 1, n):
                a = A[i, j]
                b = A[j, i]
                if a or b:
                    connected_pairs += 1
                    if a and b:
                        reciprocity_count += 1
        reciprocity = reciprocity_count / connected_pairs if connected_pairs else 0.0
        graph = nx.from_numpy_array(A.T, create_using=nx.DiGraph)
        largest_scc = max((len(c) for c in nx.strongly_connected_components(graph)), default=0)
        rows.append(
            {
                "layer_id": layer_id,
                "group": GROUP_MAP[layer_id],
                "density": density,
                "reciprocity": reciprocity,
                "largest_scc_size": largest_scc,
            }
        )
    return pd.DataFrame(rows).set_index("layer_id").loc[LAYER_ORDER]


def compute_jaccard(layers: dict[str, pd.DataFrame]) -> pd.DataFrame:
    L = len(LAYER_ORDER)
    J = np.zeros((L, L), dtype=float)
    for a, la in enumerate(LAYER_ORDER):
        Ea = layers[la].values != 0
        for b, lb in enumerate(LAYER_ORDER):
            Eb = layers[lb].values != 0
            inter = np.logical_and(Ea, Eb).sum()
            union = np.logical_or(Ea, Eb).sum()
            J[a, b] = inter / union if union else np.nan
    return pd.DataFrame(J, index=LAYER_ORDER, columns=LAYER_ORDER)


def compute_hubs_generalists(layers: dict[str, pd.DataFrame]) -> pd.DataFrame:
    node_ids = layers[LAYER_ORDER[0]].index.tolist()
    node_out = {node: 0 for node in node_ids}
    node_layer_out = {node: [] for node in node_ids}

    for layer_id in LAYER_ORDER:
        A = (layers[layer_id].values > 0).astype(int)
        out_degree = A.sum(axis=0)
        for idx, node in enumerate(node_ids):
            val = float(out_degree[idx])
            node_out[node] += val
            node_layer_out[node].append(val)

    records = []
    for node in node_ids:
        per_layer = np.array(node_layer_out[node], dtype=float)
        total = per_layer.sum()
        if total > 0:
            p = 1 - np.sum((per_layer / total) ** 2)
        else:
            p = 0.0
        records.append(
            {
                "node_id": node,
                "out_degree": node_out[node],
                "participation_coefficient": p,
            }
        )
    return pd.DataFrame(records)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.06, 1.03, label, transform=ax.transAxes, fontsize=16, fontweight="bold", ha="left", va="bottom")


def style_axis(ax: plt.Axes, grid_y: bool = False) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_y:
        ax.grid(axis="y", color="#d9dde3", linewidth=0.5)
        ax.set_axisbelow(True)


def figure1_layer_structure_summary(summary: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(13.8, 4.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.28)
    metric_specs = [
        ("density", "Edge density", "Density"),
        ("reciprocity", "Mutual dyad fraction", "Reciprocity"),
        ("largest_scc_size", "Largest strongly connected core", "Nodes"),
    ]
    layer_colors = [GROUP_COLORS[summary.loc[layer, "group"]] for layer in LAYER_ORDER]
    for idx, (metric, title, ylabel) in enumerate(metric_specs):
        ax = fig.add_subplot(gs[0, idx])
        values = summary.loc[LAYER_ORDER, metric].values
        ax.bar(np.arange(len(LAYER_ORDER)), values, color=layer_colors, width=0.8, edgecolor="white", linewidth=0.6)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xticks(np.arange(len(LAYER_ORDER)))
        ax.set_xticklabels([layer.replace("_", "\n") if len(layer) > 6 else layer for layer in LAYER_ORDER], rotation=40, ha="right")
        style_axis(ax, grid_y=True)
        if idx == 0:
            panel_label(ax, "a")
    handles = [Line2D([0], [0], color=GROUP_COLORS[g], lw=4, label=g) for g in ["DA", "IA", "MM", "MC"]]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=4, frameon=False, handlelength=1.2, columnspacing=1.5)
    save(fig, "fig1_layer_structure_summary_nature")


def figure2_layer_jaccard_heatmap(jaccard_df: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(8.8, 7.4))
    ax = fig.add_subplot(111)
    im = ax.imshow(jaccard_df.values, cmap="cividis", vmin=0, vmax=1, aspect="equal")
    ax.set_xticks(np.arange(len(LAYER_ORDER)))
    ax.set_yticks(np.arange(len(LAYER_ORDER)))
    ax.set_xticklabels(LAYER_ORDER, rotation=45, ha="right")
    ax.set_yticklabels(LAYER_ORDER)
    ax.set_title("Layer-layer edge-set overlap")
    for i in range(len(LAYER_ORDER)):
        for j in range(len(LAYER_ORDER)):
            value = jaccard_df.iloc[i, j]
            color = "white" if value < 0.55 else "#111827"
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7.3, color=color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Jaccard similarity")
    cbar.ax.tick_params(labelsize=7)
    ax.set_xticks(np.arange(-0.5, len(LAYER_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(LAYER_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)
    panel_label(ax, "a")
    save(fig, "fig2_layer_jaccard_heatmap_nature")


def figure4_hubs_vs_generalists(hg_df: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(8.0, 6.2))
    ax = fig.add_subplot(111)

    x = hg_df["out_degree"].to_numpy(dtype=float)
    y = hg_df["participation_coefficient"].to_numpy(dtype=float)
    x_med = float(np.median(x))
    y_med = float(np.median(y))

    xmin, xmax = x.min() - 4, x.max() + 4
    ymin, ymax = y.min() - 0.01, y.max() + 0.015
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # Quadrant washes.
    ax.axvspan(xmin, x_med, ymin=(y_med - ymin) / (ymax - ymin), ymax=1, color="#f3f4f6", zorder=0)
    ax.axvspan(x_med, xmax, ymin=(y_med - ymin) / (ymax - ymin), ymax=1, color="#ecfeff", zorder=0)
    ax.axvspan(xmin, x_med, ymin=0, ymax=(y_med - ymin) / (ymax - ymin), color="#fafaf9", zorder=0)
    ax.axvspan(x_med, xmax, ymin=0, ymax=(y_med - ymin) / (ymax - ymin), color="#fff7ed", zorder=0)

    ax.hexbin(x, y, gridsize=16, cmap="Greys", mincnt=1, linewidths=0, alpha=0.15, zorder=1)
    ax.scatter(
        x,
        y,
        s=68,
        color="#6b8fb9",
        edgecolor="white",
        linewidth=0.9,
        alpha=0.95,
        zorder=2,
    )

    ax.axvline(x_med, color="#94a3b8", linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)
    ax.axhline(y_med, color="#94a3b8", linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)
    ax.set_xlabel("Aggregate out-degree (sum over layers)")
    ax.set_ylabel("Participation coefficient")
    ax.set_title("Hubs versus generalists")
    style_axis(ax, grid_y=False)
    ax.grid(color="#e5e7eb", linewidth=0.6)
    ax.set_axisbelow(True)

    ax.text(xmin + 2.2, ymax - 0.012, "Broad generalists", fontsize=8.2, color="#4b5563", fontweight="bold")
    ax.text(x_med + 2.0, ymax - 0.012, "High-output generalists", fontsize=8.2, color="#0f766e", fontweight="bold")
    ax.text(x_med + 2.0, ymin + 0.012, "High-output specialists", fontsize=8.2, color="#c2410c", fontweight="bold")

    highlight_nodes = {
        "MS3 10": {"offset": (1.0, 0.010), "color": "#7c2d12"},
        "MS1 14": {"offset": (1.0, 0.002), "color": "#7c2d12"},
        "B201": {"offset": (1.0, -0.004), "color": "#7c2d12"},
        "S524": {"offset": (1.0, 0.004), "color": "#0f766e"},
        "S1707": {"offset": (1.0, 0.002), "color": "#4b5563"},
    }
    for _, row in hg_df.iterrows():
        if row["node_id"] in highlight_nodes:
            cfg = highlight_nodes[row["node_id"]]
            dx, dy = cfg["offset"]
            color = cfg["color"]
            ax.scatter(
                row["out_degree"],
                row["participation_coefficient"],
                s=92,
                color=color,
                edgecolor="white",
                linewidth=1.0,
                zorder=3,
            )
            ax.plot(
                [row["out_degree"], row["out_degree"] + dx * 0.78],
                [row["participation_coefficient"], row["participation_coefficient"] + dy * 0.78],
                color=color,
                linewidth=0.8,
                zorder=3,
            )
            ax.text(
                row["out_degree"] + dx,
                row["participation_coefficient"] + dy,
                row["node_id"],
                fontsize=7.8,
                color=color,
                va="center",
                ha="left",
                zorder=4,
            )
    panel_label(ax, "a")
    save(fig, "fig4_hubs_vs_generalists_nature")


def main() -> None:
    set_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    layers = load_layers()
    summary = compute_layer_summary(layers)
    jaccard = compute_jaccard(layers)
    hg_df = compute_hubs_generalists(layers)
    figure1_layer_structure_summary(summary)
    figure2_layer_jaccard_heatmap(jaccard)
    figure4_hubs_vs_generalists(hg_df)
    print(f"Wrote recreated figures to {OUTDIR}")


if __name__ == "__main__":
    main()
