#!/usr/bin/env python3
"""— Per-layer directional asymmetry + row/column convention anchor."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd

from data_io import REPORT_PATH, TABLE_DIR, ensure_output_dirs, load_tensor


PHENOTYPE_DOC_TEMPLATE = (
    "Row i = receiver; col j = sender. M[i,j]=1 → j elicited {phenotype} in i."
)


def asymmetry_stats(A: np.ndarray) -> dict[str, float]:
    n = A.shape[0]
    out_deg = A.sum(axis=0)
    in_deg = A.sum(axis=1)
    rho_out = float(np.std(out_deg))
    rho_in = float(np.std(in_deg))
    per_node_diff = (out_deg - in_deg)
    return {
        "n": int(n),
        "edges": int(A.sum()),
        "mean_out_degree": float(out_deg.mean()),
        "mean_in_degree": float(in_deg.mean()),
        "std_out_degree": rho_out,
        "std_in_degree": rho_in,
        "max_out_in_imbalance": float(np.max(np.abs(per_node_diff))),
        "mean_out_in_imbalance": float(np.mean(np.abs(per_node_diff))),
    }


def main() -> None:
    ensure_output_dirs()
    tensor, layer_ids, _ = load_tensor()
    rows = []
    for i, layer in enumerate(layer_ids):
        A = tensor[i].astype(int)
        stats = asymmetry_stats(A)
        rows.append({
            "phenotype": layer,
            "convention_anchor": "rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent)",
            "phenotype_scoring_rule": PHENOTYPE_DOC_TEMPLATE.format(phenotype=layer),
            **stats,
        })
    df = pd.DataFrame(rows)
    df.to_csv(TABLE_DIR / "directionality_summary.csv", index=False)

    body = f"""
Convention: row i = receiver, col j = sender. M[i,j]=1 → j elicited response l in i. NetworkX is fed M.T so A[u,v] = u→v.

Asymmetry stats below. mean_out vs mean_in tells dominant scored direction; large max_out_in_imbalance = a few hub emitters/receivers.

{(__import__('report_utils').dataframe_to_md(df.drop(columns=['phenotype_scoring_rule']), index=False))}

Outputs: directionality_summary.csv. The per-layer scoring rules are given in Supplementary Methods S2.
"""
    from report_utils import replace_section
    replace_section(REPORT_PATH, "<!-- directionality RESULTS -->", "Edge-direction convention", body)


if __name__ == "__main__":
    main()
