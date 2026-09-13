# Which file is which paper figure

Every figure here is written by a script in `code/`, as PNG, SVG and PDF. File
names describe what the figure shows rather than where it happens to sit in the
manuscript, so that renumbering the paper does not require renaming outputs. The
mapping to paper numbers is below.

## Main text

| Paper | File stem | Written by |
|---|---|---|
| Figure 1 | `figure_plates` | `figure_plates.py` — plate photographs, from `data/plates/` |
| Figure 2 | `figure_layer_overlap` | `figure_layer_overlap.py` |
| Figure 3 | `figure_strategy_space` | `figure_strategy_space.py` |
| Figure 4 | `figure_role_stability` | `figure_role_stability.py` |
| Figure 5 | `figure_phylogeny` | `figures_phylogeny.py` |

## Supplement

| Paper | File stem | Written by |
|---|---|---|
| Figure S1 | `figure_layer_architecture` | `figures_architecture.py` |
| Figure S2 | `figure_raw_similarity` | `figures_architecture.py` |
| Figure S3 | `figure_network_nulls` | `reciprocity_transitivity.py` — SCC, reciprocity, transitivity |
| Figure S4 | `figure_parallel_analysis` | `ordination.py` |
| Figure S5 | `figure_mixture_selection` | `mixture_models.py` |
| Figure S6 | — | composed for the manuscript from per-category output |
| Figure S7 | `figure_embedding_archetypes` | `figures_architecture.py` |
| Figure S8 | — | consensus network, composed for the manuscript |
| Figure S9 | `figure_lambda_calibration` | `figures_phylogeny.py` |
| Figure S10 | — | layer-specific Mantel panels, composed for the manuscript |
| Figure S11 | `figure_metabolic_niche` | `figures_phylogeny.py` |

## Diagnostics, not used in the paper

| File stem | Shows |
|---|---|
| `connectivity` | largest SCC against its degree-preserving null, per layer |
| `dominance` | cross-layer David's-score correlation and the dominance PCA scree |
| `layer_overlap_heatmap` | Jaccard overlap alone; Figure S2 shows the same data beside cosine |
| `ordination_biplot` | 12-layer PCA biplot with layer loadings |
| `figure_layer_descriptors` | per-layer density and reciprocity bars |
