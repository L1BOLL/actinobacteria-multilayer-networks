# What each table file is

One CSV per analysis output. A table the paper prints is named for its **paper
number** (`table1_*`, `tableS2_*`); everything else is named for the
stage in `code/` that wrote it, and is an intermediate the paper quotes from rather
than reproduces.

## Tables in the paper

| Paper | File | Written by |
|---|---|---|
| Table 1 | `table1_layer_architecture.csv` | `tables_1_2.py` |
| Table 2 | `table2_layer_overlap.csv` | `tables_1_2.py` |
| Table S1 | — | strain metadata, not generated here |
| Table S2 | `tableS2_layer_overlap_full.csv` | `layer_overlap.py` |
| Table S3 | `tableS3_dominance_concordance.csv` | `dominance.py` |
| Table S4 | `tableS4_taxon_sensitivity.csv` | `taxon_sensitivity.py` |

## Everything else — per-stage outputs

Intermediate results. Each is the full output of one analysis; the paper quotes
selected values from them, and `verify_reported_values.py` checks those quotes.

| File | Contents |
|---|---|
| `mixture_summary.csv` | mixture-model summary |
| `mixture_model_selection.csv` | BIC / ICL / CV log-likelihood / silhouette for k = 1–8 × 3 covariance types |
| `mixture_cluster_centroids.csv`, `mixture_cluster_membership.csv` | the k > 1 solutions, retained for inspection only |
| `mixture_diptest.csv` | Hartigan dip test on leading components |
| `ordination_variance.csv` | 12-D variance explained + bootstrap CIs + Horn null |
| `ordination_loadings.csv` | PC1–PC3 loadings per layer |
| `ordination_variance_4d.csv` | 4-D embedding variance |
| `ordination_pc1_concordance.csv` | 12-D PC1 vs 4-D PC1 rank correlation |
| `tableS2_layer_overlap_full.csv` | all 66 layer pairs: J, null mean/sd, z, p, BH and BY q |
| `connectivity_null.csv` | largest SCC vs null, per layer |
| `layer_descriptors.csv` | per-layer density, both reciprocity conventions, degree |
| `phylogenetic_signal.csv` | Pagel's λ, 21 traits × 3 topologies, permutation and χ² p-values |
| `mantel.csv` | 16S distance vs interaction distance, overall and per category |
| `phylogenetic_vcv_diagnostics.csv` | tree tips, VCV rank, ultrametricity |
| `directionality_summary.csv` | row/column convention check across layers |
| `dominance_davids_scores.csv` | David's score per strain per layer |
| `tableS3_dominance_concordance.csv` | Bradley–Terry strengths + concordance with David's |
| `dominance_cross_layer_correlation.csv` | 12 × 12 cross-layer Spearman matrix |
| `metabolic_niche_profiles.csv` | assimilation capacity per strain |
| `metabolic_niche_strain.csv` | assimilation vs network position, 18 tests |
| `metabolic_niche_dyad.csv` | metabolic dissimilarity vs interaction, per layer |
| `tableS4_taxon_sensitivity.csv` | n = 60 vs n = 56 headline comparison |
| `taxon_sensitivity_layer_structure.csv` | per-layer density / reciprocity / SCC in both cohorts |
| `reciprocity_transitivity.csv` | reciprocity and transitivity vs null, per layer |

## Scope

Every CSV here is written by a script in `code/`; if a file is not in the tables
above, no stage produces it and it does not belong in this directory. The CSV is the
machine-readable record — the formatted table belongs in the manuscript, so Tables 1
and 2 have no Markdown copy here.
| `role_instability.csv` | per-strain role instability, out-degree spread and participation |
