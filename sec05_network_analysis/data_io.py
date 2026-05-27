from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd

from source_paths import FROZEN_ROOT, find_layer_file

SEED = 20260522
ROOT = FROZEN_ROOT
PHASE2 = ROOT / "phase2_outputs"
FIG_DIR = PHASE2 / "figures"
TABLE_DIR = PHASE2 / "tables"
DATA_DIR = PHASE2 / "data"
REPORT_PATH = PHASE2 / "phase2_report.md"
CONTRADICTIONS_PATH = PHASE2 / "contradictions.md"

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
CATEGORY_MAP = {
    "DA": ["IG", "IC", "RAC"],
    "IA": ["IAC_RAC", "IAC_RDE", "IRAC"],
    "MM": ["RP", "IRP"],
    "MC": ["CCAM", "CCVM", "CMP", "CS"],
}
STRAIN_ALIASES = {
    "S138 2": "MS138 2",
    "S625'": "S625",
    "S713'v": "S713'V",
    "S713'w": "S713'W",
}


def ensure_phase2_dirs() -> None:
    for path in [PHASE2, FIG_DIR, TABLE_DIR, DATA_DIR]:
        path.mkdir(parents=True, exist_ok=True)
    # Seed all section sentinels (rewrite, not skip).
    REPORT_PATH.write_text(
        "# Phase 2 Report\n\n"
        "<!-- P2.2 RESULTS -->\n\n## P2.2 Cluster vs continuum\n\n"
        "<!-- P2.3 RESULTS -->\n\n## P2.3 PCA on 12D out-degree space\n\n"
        "<!-- P2.4 RESULTS -->\n\n## P2.4 Per-layer degree-preserving Jaccard null\n\n"
        "<!-- P2.5 RESULTS -->\n\n## P2.5 SCC null test\n\n"
        "<!-- P2.7 RESULTS -->\n\n## P2.7 Per-layer connectivity descriptors\n\n"
        "<!-- P2.8 RESULTS -->\n\n## P2.8 Phylogenetic signal (PGLS)\n\n"
        "<!-- P2.9 RESULTS -->\n\n## P2.9 Spatial decomposition (3 configurations)\n\n"
        "<!-- P2.10 RESULTS -->\n\n## P2.10 Direction-convention audit\n\n"
        "<!-- P2.11 RESULTS -->\n\n## P2.11 Dominance: per-layer and cross-layer\n\n"
        "<!-- P2.12 RESULTS -->\n\n## P2.12 Metabolic niche overlap and layer-specific interactions\n\n",
        encoding="utf-8",
    )
    if not CONTRADICTIONS_PATH.exists():
        CONTRADICTIONS_PATH.write_text("# Contradictions Log\n", encoding="utf-8")


def load_layer_df(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    df.index = [STRAIN_ALIASES.get(str(x).strip(), str(x).strip()) for x in df.index]
    df.columns = [STRAIN_ALIASES.get(str(x).strip(), str(x).strip()) for x in df.columns]
    np.fill_diagonal(df.values, 0)
    return df


def load_layers() -> dict[str, pd.DataFrame]:
    layers = {layer: load_layer_df(find_layer_file(fname)) for layer, fname in LAYER_FILES.items()}
    ref = layers[LAYER_ORDER[0]].index.tolist()
    for layer in LAYER_ORDER:
        layers[layer] = layers[layer].reindex(index=ref, columns=ref).fillna(0)
    return layers


def load_tensor() -> tuple[np.ndarray, list[str], list[str]]:
    layers = load_layers()
    node_ids = layers[LAYER_ORDER[0]].index.tolist()
    tensor = np.stack([(layers[layer].values > 0).astype(np.uint8) for layer in LAYER_ORDER], axis=0)
    return tensor, LAYER_ORDER.copy(), node_ids


def compute_out_degree_matrix(tensor: np.ndarray) -> pd.DataFrame:
    _, layer_ids, node_ids = load_tensor()
    out_mat = tensor.sum(axis=1).T.astype(float)
    return pd.DataFrame(out_mat, index=node_ids, columns=layer_ids)


def load_embedding4() -> pd.DataFrame:
    path = ROOT / 'ecological_embeddings_4D.csv'
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path).rename(columns={"Unnamed: 0": "strain"}).set_index("strain")
    return df.loc[:, ["DA", "IA", "MM", "MC"]]


def compute_embedding4_from_tensor(tensor: np.ndarray, node_ids: list[str]) -> pd.DataFrame:
    layer_index = {layer: i for i, layer in enumerate(LAYER_ORDER)}
    rows = []
    for node_idx, strain in enumerate(node_ids):
        row = {}
        for category, layers in CATEGORY_MAP.items():
            idx = [layer_index[layer] for layer in layers]
            values = []
            for li in idx:
                values.append(tensor[li, :, node_idx].sum() / (tensor.shape[1] - 1))
            row[category] = float(np.mean(values))
        rows.append(row)
    df = pd.DataFrame(rows, index=node_ids)
    scaled = (df - df.min(axis=0)) / (df.max(axis=0) - df.min(axis=0))
    scaled = scaled.fillna(0.0) * 10.0
    return scaled


def load_or_recompute_embedding4() -> pd.DataFrame:
    try:
        saved = load_embedding4()
        tensor, _, node_ids = load_tensor()
        if saved.shape == (60, 4) and list(saved.columns) == ["DA", "IA", "MM", "MC"]:
            return saved.loc[node_ids]
    except FileNotFoundError:
        pass
    tensor, _, node_ids = load_tensor()
    out = compute_embedding4_from_tensor(tensor, node_ids)
    out.to_csv(ROOT / 'ecological_embeddings_4D.csv')
    return out


def layer_groups() -> dict[str, str]:
    out = {}
    for group, layers in CATEGORY_MAP.items():
        for layer in layers:
            out[layer] = group
    return out
