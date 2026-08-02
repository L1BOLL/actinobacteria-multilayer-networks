#!/usr/bin/env python3
"""P2.9 — Per-(layer, config) decomposition. Skips if matrices_xls_{direct,indirect,distant}/ absent."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_io import LAYER_ORDER, REPORT_PATH, ROOT, SEED, TABLE_DIR, FIG_DIR, ensure_phase2_dirs
from report_utils import dataframe_to_md, replace_section
from source_paths import LAYER_FILES


CONFIG_DIRS = {
    "direct": ROOT / "data" / "matrices_xls_direct",
    "indirect": ROOT / "data" / "matrices_xls_indirect",
    "distant": ROOT / "data" / "matrices_xls_distant",
}


def _placeholder_report(missing: list[str]) -> str:
    missing_paths = "\n".join(f"- {CONFIG_DIRS[c].relative_to(ROOT)}/{{layer}}.xlsx" for c in missing)
    return f"""
Per-configuration matrices missing for: {', '.join(missing)}

Expected:
{missing_paths}

12 xlsx per directory, same filenames as data/matrices_xls/, scored per spatial configuration.

When present: per-(layer, config) density, reciprocity (pairs), SCC; J(direct,indirect), J(indirect,distant), J(direct,distant); spatial-sensitivity = (1 − J(direct,distant)) × mean density. Outputs p2_9_spatial_decomposition.csv + figure.
"""


def _load(path: Path) -> np.ndarray:
    df = pd.read_excel(path, index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    np.fill_diagonal(df.values, 0)
    return (df.values > 0).astype(int)


def main() -> None:
    ensure_phase2_dirs()
    missing = [c for c, d in CONFIG_DIRS.items() if not d.exists()]
    if missing:
        body = _placeholder_report(missing)
        replace_section(REPORT_PATH, "<!-- P2.9 RESULTS -->", "P2.9 Spatial decomposition (3 configurations)", body)
        return

    rows = []
    for layer, fname in LAYER_FILES.items():
        per_config = {}
        for cfg, cdir in CONFIG_DIRS.items():
            path = cdir / fname
            if not path.exists():
                per_config[cfg] = None
                continue
            per_config[cfg] = _load(path)
        if any(A is None for A in per_config.values()):
            continue
        d_direct = per_config["direct"]
        d_indir = per_config["indirect"]
        d_dist = per_config["distant"]
        n = d_direct.shape[0]
        densities = {c: int(A.sum()) / (n * (n - 1)) for c, A in per_config.items()}
        def jacc(X, Y):
            return float(np.logical_and(X, Y).sum() / np.logical_or(X, Y).sum()) if np.logical_or(X, Y).sum() else np.nan
        j_di = jacc(d_direct, d_indir)
        j_id = jacc(d_indir, d_dist)
        j_dd = jacc(d_direct, d_dist)
        spatial_sensitivity = (1.0 - j_dd) * np.mean(list(densities.values()))
        rows.append({
            "layer": layer,
            "density_direct": densities["direct"],
            "density_indirect": densities["indirect"],
            "density_distant": densities["distant"],
            "J_direct_indirect": j_di,
            "J_indirect_distant": j_id,
            "J_direct_distant": j_dd,
            "spatial_sensitivity": spatial_sensitivity,
        })

    df = pd.DataFrame(rows).sort_values("spatial_sensitivity", ascending=False)
    df.to_csv(TABLE_DIR / "p2_9_spatial_decomposition.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.barh(df["layer"], df["spatial_sensitivity"], color="#0f766e")
    ax.set_xlabel("Spatial sensitivity score = (1 - J(direct, distant)) * mean density")
    ax.set_title("Phenotypes ordered by spatial sensitivity")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p2_9_spatial_decomposition.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    body = f"""
Per-(layer, config) decomposition over direct/indirect/distant.
spatial_sensitivity = (1 − J(direct,distant)) × mean(density). High = contact-dependent.

{dataframe_to_md(df, index=False)}

Outputs:
- p2_9_spatial_decomposition.csv
- p2_9_spatial_decomposition.png
"""
    replace_section(REPORT_PATH, "<!-- P2.9 RESULTS -->", "P2.9 Spatial decomposition (3 configurations)", body)


if __name__ == "__main__":
    main()
