#!/usr/bin/env python3
"""— Per-layer descriptors. Pairs-based (main) and edges-based reciprocity. Data is binary; layer constants cancel."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_io import save_figure, FIG_DIR, REPORT_PATH, SEED, TABLE_DIR, ensure_output_dirs, load_tensor
from report_utils import dataframe_to_md, replace_section


def reciprocity_both(A: np.ndarray) -> tuple[float, float, int, int, int]:
    n = A.shape[0]
    mutual_pairs = 0
    unidirectional_pairs = 0
    for a in range(n):
        for b in range(a + 1, n):
            x = A[a, b]
            y = A[b, a]
            if x and y:
                mutual_pairs += 1
            elif x or y:
                unidirectional_pairs += 1
    connected_pairs = mutual_pairs + unidirectional_pairs
    total_edges = int(A.sum())
    r_pairs = mutual_pairs / connected_pairs if connected_pairs else 0.0
    r_edges = (2 * mutual_pairs) / total_edges if total_edges else 0.0
    return r_pairs, r_edges, mutual_pairs, unidirectional_pairs, total_edges


def main() -> None:
    ensure_output_dirs()
    tensor, layer_ids, _ = load_tensor()
    rows = []
    for i, layer in enumerate(layer_ids):
        A = tensor[i].astype(int)
        n = A.shape[0]
        edges = int(A.sum())
        density = edges / (n * (n - 1))
        r_pairs, r_edges, m_pairs, uni_pairs, _ = reciprocity_both(A)
        rows.append(
            {
                "layer": layer,
                "edge_density": density,
                "reciprocity_pairs": r_pairs,
                "reciprocity_edges": r_edges,
                "mutual_pairs": m_pairs,
                "unidirectional_pairs": uni_pairs,
                "asymmetry_index_pairs": (uni_pairs / (uni_pairs + m_pairs)) if (uni_pairs + m_pairs) else 0.0,
                "mean_out_degree": float(A.sum(axis=0).mean()),
                "mean_in_degree": float(A.sum(axis=1).mean()),
                "seed": SEED,
            }
        )
    df = pd.DataFrame(rows).sort_values("edge_density", ascending=False)
    df.to_csv(TABLE_DIR / "layer_descriptors.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.6), constrained_layout=True)
    axes[0].bar(df["layer"], df["edge_density"], color="#4c78a8", edgecolor="white", linewidth=0.7)
    axes[0].set_ylabel("Edge density")
    axes[0].set_xlabel("Layer")
    axes[0].set_title("Per-layer edge density")
    axes[0].tick_params(axis="x", rotation=40)
    axes[0].grid(axis="y", color="#e5e7eb", linewidth=0.6)
    axes[0].set_axisbelow(True)

    x = np.arange(len(df))
    w = 0.4
    axes[1].bar(x - w / 2, df["reciprocity_pairs"], width=w, color="#0f766e", label="pairs-based (paper)")
    axes[1].bar(x + w / 2, df["reciprocity_edges"], width=w, color="#7f1d1d", label="edges-based")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(df["layer"], rotation=40, ha="right")
    axes[1].set_ylabel("Reciprocity")
    axes[1].set_title("Reciprocity, two conventions")
    axes[1].grid(axis="y", color="#e5e7eb", linewidth=0.6)
    axes[1].set_axisbelow(True)
    axes[1].legend(frameon=False, fontsize=9)
    save_figure(fig, "figure_layer_descriptors", FIG_DIR)
    plt.close(fig)

    D = 1e-6
    t = 518400.0
    Ldiff_cm = (D * t) ** 0.5
    Ldiff_mm = Ldiff_cm * 10.0

    body = f"""
Seed: {SEED}. Binary projection A_l = (M_l > 0). Data is intrinsically binary; each layer takes 0 or a layer constant c_l (CCAM=1; CCVM/CS/CMP/IRP/RP=2; IAC_RAC/IAC_RDE/IRAC/RAC=3; IG/IC=5). c_l is a phenotype-importance tag, not intensity.

Reciprocity:
  R_p = mutual / connected_pairs   (main, [0,1])
  R_e = 2 · mutual / edges          (= 2 R_p / (1 + R_p))
Same layer ordering under both.

Per-layer (sorted by density):
{dataframe_to_md(df, index=False)}

Diffusion order-of-magnitude:
- D ≈ 1e-6 cm²/s (small molecule); t = 6 d = 518400 s
- L_diff = √(D·t) = {Ldiff_cm:.3f} cm ≈ {Ldiff_mm:.2f} mm
- Far (~60 mm) ≫ L_diff: usable as within-plate baseline.

Outputs:
- layer_descriptors.csv  (both reciprocity, mutual / unidirectional counts)
- figure_layer_descriptors.png
"""
    replace_section(REPORT_PATH, "<!-- layer_descriptors RESULTS -->", "Per-layer connectivity descriptors", body)


if __name__ == "__main__":
    main()
