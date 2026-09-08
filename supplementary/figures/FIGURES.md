# Which file is which paper figure

Every figure here is written by a script in `code/`, as PNG and SVG. **Filenames are
the paper's own figure numbers** — `figureN_*` for the main text, `figureSN_*` for the
supplement — so this table records which script writes each one, not a translation.
A file named `p2_*` is a per-stage diagnostic that the paper does not use.

## Main text

| Paper | File stem | Written by |
|---|---|---|
| Figure 1 | — | **not generated**: plate photographs, no raw files in this deposit |
| Figure 2 | `figure2_spatial_dependence` | `figure_2_spatial.py` |
| Figure 3 | `figure3_layer_overlap_null` | `figure_3_overlap.py` |
| Figure 4 | `figure4_strategy_space` | `figure_4_strategy.py` |
| Figure 5 | `figure5_role_switching` | `generate_figs.py` |
| Figure 6 | `figure6_phylogeny` | `make_paper_figures.py` |

## Supplementary

| Paper | File stem | Written by |
|---|---|---|
| Figure S1 | `figureS1_layer_architecture` | `generate_figs.py` |
| Figure S2 | — | **not generated**: needs the per-configuration spatial matrices |
| Figure S3 | `figureS3_raw_similarity` | `generate_figs.py` — both panels: Jaccard (A) and cosine (B) |
| Figure S4 | `figureS4_network_nulls` | `p2_14_reciprocity_transitivity_null.py` — three panels: SCC, reciprocity, transitivity |
| Figure S5 | `figureS5_parallel_analysis` | `p2_3_pca12.py` |
| Figure S6 | `figureS6_mixture_selection` | `p2_2_gmm.py` |
| Figure S7 | — | **not generated**: no script produces it |
| Figure S8 | `figureS8_embedding_archetypes` | `generate_figs.py` |
| Figure S9 | — | **not generated**: consensus network, no code and no output table |
| Figure S10 | `figureS10_lambda_calibration` | `make_paper_figures.py` |
| Figure S11 | — | **not generated**: no script produces it |
| Figure S12 | `figureS12_metabolism` | `make_paper_figures.py` |

## Produced but not used in the paper

Diagnostics, kept because they are cheap and useful when checking a stage:

| File stem | What it shows |
|---|---|
| `p2_11_dominance` | cross-layer David's-score correlation + dominance PCA scree |
| `p2_4_jaccard_heatmap` | Jaccard overlap alone; Figure S3 panel A shows the same data beside cosine |
| `p2_3_pca12_biplot` | 12-D PCA biplot |
| `p2_7_proximity_density` | per-layer density and reciprocity bars |
| `p2_5_scc_null` | SCC vs null on its own; the same data is panel 1 of Figure S4 |

> **No figure colours strains by mixture assignment.** The reported conclusion is that
> no discrete clusters are supported, so a cluster-coloured ordination would contradict
> the text it illustrates. Figure 4 (`code/figure_4_strategy.py`) colours by total
> out-degree, a continuous quantity.

## Formats

PNG and SVG only, set once in `code/data_io.py`:

```python
FIGURE_FORMATS = ("png", "svg")
```

Every figure-writing script goes through `data_io.save_figure()`, so the set is
consistent. PDF was dropped: it duplicates SVG, and journals wanting vector accept SVG
or EPS.

## Why five paper figures have no file

Figure 1 is photographs. Figures S2, S7, S9 and S11 depend on inputs or code that are
not in this deposit — the per-configuration spatial matrices for S2, and no script has
ever existed for S7, S9 or S11. `p2_9_spatial_decomposition.py` is written and wired
into `run_all.py`; it skips and records the skip in `supplementary/report.md`.
