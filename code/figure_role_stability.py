#!/usr/bin/env python3
"""Directional influence across layers: is a strain's role fixed or phenotype dependent?

Within a layer a strain is a net sender when its out-degree exceeds its in-degree
and a net receiver when the reverse holds; layers where the two are equal are left
unclassified.  Role instability is the number of classified layers whose class
differs from that strain's own majority class, so 0 means one role throughout and
6 means an even split.  Defined this way the measure does not depend on the order
in which the layers are taken, which a count of transitions between consecutive
layers would.
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
from mpl_toolkits.axes_grid1 import make_axes_locatable

from data_io import FIGURE_FORMATS, LAYER_ORDER, load_tensor

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "figures"
TAB = REPO / "results" / "tables"
BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8a86", "#ebebe8"


def style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans",
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
        "axes.titlelocation": "left", "axes.titlepad": 8,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.edgecolor": INK2, "axes.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "axes.labelcolor": INK, "legend.frameon": False, "legend.fontsize": 8.5,
    })


def panel_tag(ax, letter):
    ax.text(-0.09, 1.04, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="left", color=INK)


def role_table(tensor, nodes):
    out_deg = tensor.sum(axis=1).T.astype(float)       # strains x layers
    in_deg = tensor.sum(axis=2).T.astype(float)
    total = out_deg + in_deg
    with np.errstate(invalid="ignore", divide="ignore"):
        asym = np.where(total > 0, (out_deg - in_deg) / np.where(total == 0, 1, total), np.nan)

    sign = np.sign(asym)                                # +1 sender, -1 receiver, 0/nan unclassified
    sign[np.isnan(sign)] = 0
    n_sender = (sign > 0).sum(axis=1)
    n_receiver = (sign < 0).sum(axis=1)
    classified = n_sender + n_receiver
    majority = np.where(n_sender >= n_receiver, 1, -1)
    minority = np.where(majority == 1, n_receiver, n_sender)

    sd_out = out_deg.std(axis=1)
    tot = out_deg.sum(axis=1)
    safe = np.where(tot == 0, 1.0, tot)
    participation = 1.0 - np.sum((out_deg / safe[:, None]) ** 2, axis=1)

    z = (asym - np.nanmean(asym, axis=0)) / np.nanstd(asym, axis=0)
    return pd.DataFrame({
        "strain": nodes, "role_instability": minority, "classified_layers": classified,
        "sd_out_degree": sd_out, "participation": participation,
        "total_out_degree": tot,
        "majority_role": np.where(majority == 1, "sender", "receiver"),
    }), z


def main() -> None:
    style()
    tensor, layer_ids, nodes = load_tensor()
    df, z = role_table(tensor, nodes)
    df.to_csv(TAB / "role_instability.csv", index=False)

    order = np.argsort(-df["role_instability"].to_numpy())
    zz = z[order]
    labels = [nodes[i] for i in order]

    fig = plt.figure(figsize=(13.8, 8.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.30,
                          left=0.09, right=0.96, top=0.92, bottom=0.10)

    # a - every strain, every layer
    ax = fig.add_subplot(gs[0, 0])
    vmax = float(np.nanmax(np.abs(zz)))
    im = ax.imshow(zz, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(LAYER_ORDER)))
    ax.set_xticklabels(LAYER_ORDER, rotation=90, fontsize=8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.set_ylabel("Strain (ordered by role instability)")
    ax.set_title("Sender–receiver asymmetry across layers", color=INK)
    cax = make_axes_locatable(ax).append_axes("right", size="3.5%", pad=0.08)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Within-layer asymmetry (z)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8); cb.outline.set_linewidth(0.6)
    for sp in ax.spines.values():
        sp.set_visible(False)
    panel_tag(ax, "a")

    # b - the role-instability landscape, with axes that say what they are
    ax = fig.add_subplot(gs[0, 1])
    jitter = np.random.default_rng(0).normal(0, 0.055, len(df))
    sc = ax.scatter(df["sd_out_degree"], df["role_instability"] + jitter,
                    c=df["participation"], cmap="viridis", s=62,
                    edgecolor="white", linewidth=0.7, zorder=3)
    cax = make_axes_locatable(ax).append_axes("right", size="3.5%", pad=0.08)
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("Participation coefficient", fontsize=8.5)
    cb.ax.tick_params(labelsize=8); cb.outline.set_linewidth(0.6)
    ax.set_xlabel("SD of out-degree across the twelve layers")
    ax.set_ylabel("Role instability (layers disagreeing with the strain's majority role)")
    ax.set_yticks(range(int(df["role_instability"].max()) + 1))
    ax.set_title("Most strains change role between layers", color=INK)
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    # label only the extremes of the spread axis, so the annotations cannot collide
    top = pd.concat([df.nlargest(3, "sd_out_degree"), df.nsmallest(2, "sd_out_degree")])
    for _, r in top.iterrows():
        ax.annotate(r["strain"], (r["sd_out_degree"], r["role_instability"]),
                    textcoords="offset points", xytext=(8, 6), fontsize=8, color=INK2)
    n_stable = int((df["role_instability"] == 0).sum())
    ax.text(0.02, 0.02,
            f"median instability {df['role_instability'].median():.0f} of "
            f"{int(df['classified_layers'].median())} classified layers; "
            f"{n_stable} of {len(df)} strains keep one role",
            transform=ax.transAxes, fontsize=8.5, color=MUTED)
    panel_tag(ax, "b")

    for ext in FIGURE_FORMATS:
        fig.savefig(OUT / f"figure_role_stability.{ext}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  wrote figure_role_stability  "
          f"(instability {df['role_instability'].min()}–{df['role_instability'].max()}, "
          f"median {df['role_instability'].median():.0f})")


if __name__ == "__main__":
    main()
