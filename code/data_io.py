from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd

from source_paths import CACHE_DIR, EMBEDDING_CSV, FROZEN_ROOT, RESULTS_DIR, find_layer_file

SEED = 20260522
ROOT = FROZEN_ROOT

# Everything a reader of the paper needs lands in results/ and is committed.
OUTPUT_DIR = RESULTS_DIR
FIG_DIR = OUTPUT_DIR / "figures"
TABLE_DIR = OUTPUT_DIR / "tables"
REPORT_PATH = OUTPUT_DIR / "report.md"

# Null-model draws total ~570 MB and are fully regenerable from SEED, so they are
# cached outside results/ and excluded from version control.
DATA_DIR = CACHE_DIR

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
    # S1430 is the sequencing-record label for the isolate the interaction
    # matrices call S1226. Confirmed by the authors against the submission record
    # (youssef_prop/bn_tables.xlsx, sheet SuppTableS1_strains, caption). Before
    # this was confirmed, every phylogenetic analysis ran at n=59 with both labels
    # dropped; it now runs at n=60. See data/annotation/DATA_NOTES.md.
    "S1430": "S1226",
}

# Canonical form used to join the interaction matrices to the annotation data
# (16S FASTA, distance matrix, carbohydrate profile). The annotation sources use
# non-breaking spaces and drop the apostrophes that the .xlsx indices carry, so
# both sides are folded to a single key space before matching. See
# annotation_data/DATA_NOTES.md for the residual unmatched strain.
def canonical_strain(name: object) -> str:
    """Fold a strain label to the key space shared by matrices and annotations."""
    s = str(name).replace("\xa0", " ").replace("'", "").replace("’", "")
    s = " ".join(s.split())
    return STRAIN_ALIASES.get(s, s)


REPORT_SECTIONS = [
    ("mixture_models", "Discrete groups versus a continuum"),
    ("ordination", "Principal components of the 12-layer out-degree profile"),
    ("layer_overlap", "Cross-layer overlap against a degree-preserving null"),
    ("connectivity", "Strongly connected components against the same null"),
    ("layer_descriptors", "Per-layer connectivity descriptors"),
    ("phylogenetic_signal", "Phylogenetic signal"),
    ("directionality", "Edge-direction convention"),
    ("dominance", "Directional influence, per layer and across layers"),
    ("metabolic_niche", "Metabolic niche and network position"),
    ("taxon_sensitivity", "Sensitivity to the non-Streptomyces isolates"),
    ("reciprocity_transitivity", "Reciprocity and transitivity against the null"),
]


def ensure_output_dirs() -> None:
    """Create output dirs and make sure every report sentinel exists.

    Every stage calls this at start-up, so it must not rewrite the report
    wholesale: that would let whichever stage ran last erase the sections the
    earlier stages had just written, which is why the report
    only ever contained one populated section. Missing sentinels are appended;
    existing content is left alone, and `report_utils.replace_section` overwrites
    only the section it owns.
    """
    for path in [OUTPUT_DIR, FIG_DIR, TABLE_DIR, DATA_DIR]:
        path.mkdir(parents=True, exist_ok=True)

    text = REPORT_PATH.read_text(encoding="utf-8") if REPORT_PATH.exists() else "# Analysis report\n"
    if not text.strip():
        text = "# Analysis report\n"
    for code, title in REPORT_SECTIONS:
        sentinel = f"<!-- {code} RESULTS -->"
        if sentinel not in text:
            if not text.endswith("\n"):
                text += "\n"
            text += f"\n{sentinel}\n\n## {code} {title}\n\n"
    REPORT_PATH.write_text(text, encoding="utf-8")



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
    path = EMBEDDING_CSV
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
    out.to_csv(EMBEDDING_CSV)
    return out


def layer_groups() -> dict[str, str]:
    out = {}
    for group, layers in CATEGORY_MAP.items():
        for layer in layers:
            out[layer] = group
    return out


# Figure output formats. PNG for review and for pasting into slide decks, SVG for
# typesetting. No PDF: it duplicates SVG, and journals that want vector take SVG or
# EPS. Every figure-producing script goes through this, so the set is consistent.
FIGURE_FORMATS = ("png", "svg", "pdf")  # pdf for journals that will not take svg


def save_figure(fig, stem: str, outdir=None) -> None:
    """Write a figure to FIG_DIR (or outdir) once per format in FIGURE_FORMATS."""
    target = FIG_DIR if outdir is None else outdir
    target.mkdir(parents=True, exist_ok=True)
    for ext in FIGURE_FORMATS:
        fig.savefig(target / f"{stem}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
