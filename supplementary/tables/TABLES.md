# What each table file is

One CSV per analysis output. A table the paper prints is named for its **paper
number** (`table1_*`, `tableS2_*`); everything else keeps the `p2_*` name of the
stage in `code/` that wrote it, and is an intermediate the paper quotes from rather
than reproduces.

## Tables in the paper

| Paper | File | Written by |
|---|---|---|
| Table 1 | `table1_layer_architecture.csv` | `tables_1_2.py` |
| Table 2 | `table2_layer_overlap.csv` | `tables_1_2.py` |
| Table S1 | — | strain metadata, not generated here |
| Table S2 | `tableS2_layer_overlap_full.csv` | `p2_4_jaccard_null.py` |
| Table S3 | `tableS3_dominance_concordance.csv` | `p2_11_dominance.py` |
| Table S4 | `tableS4_taxon_sensitivity.csv` | `p2_13_taxon_sensitivity.py` |

## Everything else — per-stage outputs

Intermediate results. Each is the full output of one analysis; the paper quotes
selected values from them, and `verify_reported_values.py` checks those quotes.

| File | Contents |
|---|---|
| `p2_2_summary.csv` | mixture-model summary |
| `p2_2_gmm_modelselect.csv` | BIC / ICL / CV log-likelihood / silhouette for k = 1–8 × 3 covariance types |
| `p2_2_cluster_centroids.csv`, `p2_2_cluster_membership.csv` | the k > 1 solutions, retained for inspection only |
| `p2_2_diptest.csv` | Hartigan dip test on leading components |
| `p2_3_pca12_variance.csv` | 12-D variance explained + bootstrap CIs + Horn null |
| `p2_3_pca12_loadings.csv` | PC1–PC3 loadings per layer |
| `p2_3_pca4_variance.csv` | 4-D embedding variance |
| `p2_3_pc1_spearman.csv` | 12-D PC1 vs 4-D PC1 rank correlation |
| `tableS2_layer_overlap_full.csv` | all 66 layer pairs: J, null mean/sd, z, p, BH and BY q |
| `p2_5_scc_null.csv` | largest SCC vs null, per layer |
| `p2_7_proximity_sensitivity.csv` | per-layer density, both reciprocity conventions, degree |
| `p2_8_pgls.csv` | Pagel's λ, 21 traits × 3 topologies, permutation and χ² p-values |
| `p2_8_mantel.csv` | 16S distance vs interaction distance, overall and per category |
| `p2_8_vcv_diagnostics.csv` | tree tips, VCV rank, ultrametricity |
| `p2_10_directionality_audit.csv` | row/column convention check across layers |
| `p2_11_davids_per_layer.csv` | David's score per strain per layer |
| `tableS3_dominance_concordance.csv` | Bradley–Terry strengths + concordance with David's |
| `p2_11_cross_layer_dominance_corr.csv` | 12 × 12 cross-layer Spearman matrix |
| `p2_12_metabolism_profiles.csv` | assimilation capacity per strain |
| `p2_12_metabolism_strain.csv` | assimilation vs network position, 18 tests |
| `p2_12_metabolism_dyad.csv` | metabolic dissimilarity vs interaction, per layer |
| `tableS4_taxon_sensitivity.csv` | n = 60 vs n = 56 headline comparison |
| `p2_13_layer_structure_n56.csv` | per-layer density / reciprocity / SCC in both cohorts |
| `p2_14_reciprocity_transitivity_null.csv` | reciprocity and transitivity vs null, per layer |

## Scope

Every CSV here is written by a script in `code/`; if a file is not in the tables
above, no stage produces it and it does not belong in this directory. The CSV is the
machine-readable record — the formatted table belongs in the manuscript, so Tables 1
and 2 have no Markdown copy here.
