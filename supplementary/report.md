# Phase 2 Report

<!-- P2.2 RESULTS -->

## P2.2 Cluster vs continuum

Categories: DA={IG, IC, RAC}, IA={IAC_RAC, IAC_RDE, IRAC}, MM={RP, IRP}, MC={CCAM, CCVM, CMP, CS}. Seed: 20260522.

GMM k=1..8 fit on z-scored 4D and 12D embeddings, three covariances (full, diag, spherical).
Criteria: BIC, ICL (BIC + 2 × posterior entropy), 10-fold CV log-likelihood, silhouette on GMM hard assignments. Hartigan dip on PC1/PC2.

Diagonal cov (regularized) headline:

- `4D_codebase_embedding` (diagonal cov): BIC k=8, ICL k=8, CV-LL k=2, GMM-silhouette k=4 (s=0.478)
- `12D_outdegree_profile` (diagonal cov): BIC k=7, ICL k=7, CV-LL k=3, GMM-silhouette k=2 (s=0.239)

Best-model summary:

| embedding | covariance | k_bic | k_icl | k_cv | k_gmm_silhouette | best_gmm_silhouette | dip_pc1_p | dip_pc2_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4D_codebase_embedding | full | 4 | 4 | 1 | 3 | 0.495639 | 0.53096 | 0.987987 |
| 4D_codebase_embedding | diag | 8 | 8 | 2 | 4 | 0.478128 | 0.53096 | 0.987987 |
| 4D_codebase_embedding | spherical | 7 | 7 | 2 | 2 | 0.436014 | 0.53096 | 0.987987 |
| 12D_outdegree_profile | full | 5 | 5 | 1 | 2 | 0.350735 | 0.966953 | 0.543804 |
| 12D_outdegree_profile | diag | 7 | 7 | 3 | 2 | 0.238954 | 0.966953 | 0.543804 |
| 12D_outdegree_profile | spherical | 7 | 7 | 2 | 2 | 0.343525 | 0.966953 | 0.543804 |

Rules:
- k=1 rejected only if BIC, ICL, CV agree.
- Full-cov BIC unreliable at N=60. Check `p2_2_cluster_centroids.csv` for singletons (≤3 members) → overfit.
- Dip p>0.05 = consistent with unimodality on top axes (necessary, not sufficient).

Dip test:

| embedding | component | dip | pvalue | seed |
| --- | --- | --- | --- | --- |
| 4D_codebase_embedding | PC1 | 0.0441054 | 0.53096 | 20260522 |
| 4D_codebase_embedding | PC2 | 0.0291256 | 0.987987 | 20260522 |
| 12D_outdegree_profile | PC1 | 0.031157 | 0.966953 | 20260522 |
| 12D_outdegree_profile | PC2 | 0.0438039 | 0.543804 | 20260522 |

Outputs:
- p2_2_gmm_modelselect.csv (k × cov: BIC, ICL, CV, silhouette)
- p2_2_diptest.csv
- p2_2_summary.csv
- p2_2_cluster_centroids.csv (ICL-preferred k only)
- p2_2_cluster_membership.csv
- p2_2_bic_curves.png
<!-- P2.3 RESULTS -->

## P2.3 PCA on 12D out-degree space

Seed: 20260522. PCA on z-scored 60×12 out-degree. 2000 bootstrap resamples; 1000-rep Horn parallel analysis on N(0,1) data of same shape. 4D collapsed-embedding PCA reported alongside.

12D: PC1 = 32.5%, PC2 = 17.3% — honest dimensionality.
4D:  PC1 = 53.0%, PC2 = 23.9% — visualization only (low-D by construction; inflates variance).

Horn: 3 of 12 PCs above the null.

12D-PC1 ↔ 4D-PC1 Spearman:
- rho = 0.9337, p = 1.498e-27, rho<0.7 = False
- rho ≥ 0.7 → 4D PC1 preserves dominant 12D axis.

12D variance (top 6):
| component | variance_explained | cumulative_variance | var_ci_lo | var_ci_hi | parallel_null_mean | significant_vs_parallel | seed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PC1 | 0.324542 | 0.324542 | 0.241516 | 0.436983 | 0.151046 | True | 20260522 |
| PC2 | 0.173178 | 0.497721 | 0.133702 | 0.220967 | 0.13084 | True | 20260522 |
| PC3 | 0.118223 | 0.615944 | 0.0980192 | 0.149678 | 0.116489 | True | 20260522 |
| PC4 | 0.0822774 | 0.698221 | 0.069079 | 0.114163 | 0.104288 | False | 20260522 |
| PC5 | 0.0728622 | 0.771084 | 0.0529025 | 0.089101 | 0.0935464 | False | 20260522 |
| PC6 | 0.0547414 | 0.825825 | 0.0406466 | 0.0729012 | 0.0836163 | False | 20260522 |

4D variance:
| component | variance_explained | cumulative_variance |
| --- | --- | --- |
| PC1 | 0.529667 | 0.529667 |
| PC2 | 0.238548 | 0.768216 |
| PC3 | 0.160916 | 0.929132 |
| PC4 | 0.0708682 | 1 |

12D PC1–PC3 loadings:
| layer | PC1_loading | PC2_loading | PC3_loading |
| --- | --- | --- | --- |
| CCAM | 0.366662 | 0.254317 | -0.0328633 |
| CCVM | 0.383257 | 0.200065 | -0.158021 |
| CMP | 0.411612 | 0.152662 | 0.183155 |
| CS | 0.294942 | 0.115493 | -0.124635 |
| IG | 0.331888 | -0.377094 | -0.0674517 |
| IC | 0.251299 | -0.381998 | -0.289147 |
| RP | 0.265176 | 0.325654 | 0.365758 |
| IRP | 0.398654 | -0.061841 | -0.204521 |
| IAC_RAC | -0.18548 | 0.290428 | -0.3889 |
| IAC_RDE | -0.00232741 | -0.14586 | 0.650355 |
| IRAC | -0.156766 | 0.250362 | -0.286191 |
| RAC | -0.0270066 | 0.540695 | 0.0706832 |

Outputs:
- p2_3_pca12_variance.csv (with bootstrap CI + parallel-null)
- p2_3_pca4_variance.csv
- p2_3_pca12_loadings.csv
- p2_3_pc1_spearman.csv
- p2_3_pca12_scree.png, p2_3_pca12_biplot.png
<!-- P2.4 RESULTS -->

## P2.4 Per-layer degree-preserving Jaccard null

Seed: 20260522. 10000 degree-preserving directed-swap nulls per layer; Phipson–Smyth (M+1)/(N+1); FDR across 66 pairs by BH (p_fdr_bh) and BY (p_fdr_by, dependence-robust, default).

CCAM/CCVM (methodological confound, reported separately):
| layer1 | layer2 | J_obs | J_null_mean | J_null_sd | z | p_raw | n_null_draws | methodological_confound | seed | p_fdr_bh | significant_bh | p_fdr_by | significant_by |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CCAM | CCVM | 0.335849 | 0.129225 | 0.0127592 | 16.1941 | 9.999e-05 | 10000 | True | 20260522 | 0.000942763 | True | 0.00450115 | True |

Biological summary (CCAM/CCVM excluded, BY-FDR):
- higher_by = 8
- lower_by  = 1
- BH: higher=8, lower=1
- CCAM/CCVM sig BY = True (BH: True)

Reading: significant deviations are predominantly over-overlap, in mechanistically coherent clusters (e.g. IRAC/RAC, CCVM/CMP, IG/IRAC). The v14 "most pairs significantly lower than chance" claim is overturned.

Top positive biological z (BY-FDR):
| layer1 | layer2 | J_obs | J_null_mean | J_null_sd | z | p_raw | n_null_draws | methodological_confound | seed | p_fdr_bh | significant_bh | p_fdr_by | significant_by |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| IRAC | RAC | 0.166375 | 0.133056 | 0.00360596 | 9.23991 | 9.999e-05 | 10000 | False | 20260522 | 0.000942763 | True | 0.00450115 | True |
| CCVM | CMP | 0.108225 | 0.0473328 | 0.00739584 | 8.23332 | 9.999e-05 | 10000 | False | 20260522 | 0.000942763 | True | 0.00450115 | True |
| CCAM | CMP | 0.120743 | 0.0530328 | 0.0088099 | 7.68569 | 9.999e-05 | 10000 | False | 20260522 | 0.000942763 | True | 0.00450115 | True |
| CCVM | CS | 0.0616016 | 0.025375 | 0.00565303 | 6.40836 | 9.999e-05 | 10000 | False | 20260522 | 0.000942763 | True | 0.00450115 | True |
| IG | IRAC | 0.0923695 | 0.0647798 | 0.00571994 | 4.82342 | 9.999e-05 | 10000 | False | 20260522 | 0.000942763 | True | 0.00450115 | True |
| CCAM | CS | 0.045584 | 0.0226244 | 0.00595178 | 3.85761 | 0.00059994 | 10000 | False | 20260522 | 0.0039596 | True | 0.0189048 | True |
| CMP | CS | 0.11465 | 0.0810458 | 0.00912262 | 3.68358 | 0.00059994 | 10000 | False | 20260522 | 0.0039596 | True | 0.0189048 | True |
| CS | RAC | 0.0934685 | 0.0854961 | 0.00217159 | 3.67121 | 0.00039996 | 10000 | False | 20260522 | 0.00329967 | True | 0.015754 | True |
| IG | IC | 0.0900123 | 0.077525 | 0.00482186 | 2.58972 | 0.0107989 | 10000 | False | 20260522 | 0.0647935 | False | 0.309352 | False |
| CS | RP | 0.0406091 | 0.0297469 | 0.0043175 | 2.51587 | 0.0168983 | 10000 | False | 20260522 | 0.0929407 | False | 0.443739 | False |

Outputs:
- p2_4_jaccard_null.csv
- p2_4_jaccard_heatmap.png  (* = BY-FDR sig biological, (c) = confound)
<!-- P2.5 RESULTS -->

## P2.5 SCC null test

Seed: 20260522. Observed largest SCC vs 10000 degree-preserving nulls per layer; Phipson–Smyth (M+1)/(N+1); FDR across 12 layers by BH and BY.

| layer | scc_obs | scc_null_mean | scc_null_sd | z | p_raw | n_null_draws | seed | p_fdr_bh | significant_bh | p_fdr_by | significant_by |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CCAM | 4 | 3.8079 | 2.48081 | 0.0774344 | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| CCVM | 24 | 22.4109 | 3.41181 | 0.465765 | 0.674433 | 10000 | 20260522 | 1 | False | 1 | False |
| CMP | 39 | 38.6006 | 1.55516 | 0.256823 | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| CS | 26 | 26.36 | 2.82634 | -0.127373 | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| IG | 59 | 58.9902 | 0.0995236 | 0.0984691 | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| IC | 2 | 5.3919 | 2.62918 | -1.2901 | 0.231777 | 10000 | 20260522 | 1 | False | 1 | False |
| RP | 3 | 4.255 | 0.857233 | -1.46401 | 0.248475 | 10000 | 20260522 | 1 | False | 1 | False |
| IRP | 1 | 1 | 0 | nan | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| IAC_RAC | 46 | 46.2622 | 1.19247 | -0.219879 | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| IAC_RDE | 1 | 1 | 0 | nan | 1 | 10000 | 20260522 | 1 | False | 1 | False |
| IRAC | 41 | 40.2639 | 1.39644 | 0.527125 | 0.755824 | 10000 | 20260522 | 1 | False | 1 | False |
| RAC | 51 | 50.5771 | 0.740345 | 0.57122 | 1 | 10000 | 20260522 | 1 | False | 1 | False |

- SCC above null (BY-sig): none
- SCC consistent with null: CCAM, CCVM, CMP, CS, IG, IC, RP, IRP, IAC_RAC, IAC_RDE, IRAC, RAC

SCC size is a consequence of density, not extra feedback structure. v14 "MM/MC SCCs 20–40 nodes" is wrong: actual distribution is bimodal across MC (CMP 39, CS 26, CCVM 24 vs CCAM 4, IRP 1, IAC_RDE 1, IC 2, RP 3).

Outputs:
- p2_5_scc_null.csv
- p2_5_scc_null.png
<!-- P2.7 RESULTS -->

## P2.7 Per-layer connectivity descriptors

Seed: 20260522. Binary projection A_l = (M_l > 0). Data is intrinsically binary; each layer takes 0 or a layer constant c_l (CCAM=1; CCVM/CS/CMP/IRP/RP=2; IAC_RAC/IAC_RDE/IRAC/RAC=3; IG/IC=5). c_l is a phenotype-importance tag, not intensity.

Reciprocity:
  R_p = mutual / connected_pairs   (main, [0,1])
  R_e = 2 · mutual / edges          (= 2 R_p / (1 + R_p))
Same layer ordering under both.

Per-layer (sorted by density):
| layer | edge_density | reciprocity_pairs | reciprocity_edges | mutual_pairs | unidirectional_pairs | asymmetry_index_pairs | mean_out_degree | mean_in_degree | seed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RAC | 0.473729 | 0.296984 | 0.457961 | 384 | 909 | 0.703016 | 27.95 | 27.95 | 20260522 |
| IG | 0.216667 | 0.149925 | 0.260756 | 100 | 567 | 0.850075 | 12.7833 | 12.7833 | 20260522 |
| IAC_RAC | 0.108192 | 0.0974212 | 0.177546 | 34 | 315 | 0.902579 | 6.38333 | 6.38333 | 20260522 |
| IRAC | 0.090678 | 0.0993151 | 0.180685 | 29 | 263 | 0.900685 | 5.35 | 5.35 | 20260522 |
| CS | 0.0748588 | 0.06 | 0.113208 | 15 | 235 | 0.94 | 4.41667 | 4.41667 | 20260522 |
| CMP | 0.0734463 | 0.0441767 | 0.0846154 | 11 | 238 | 0.955823 | 4.33333 | 4.33333 | 20260522 |
| CCVM | 0.0711864 | 0.0285714 | 0.0555556 | 7 | 238 | 0.971429 | 4.2 | 4.2 | 20260522 |
| RP | 0.0409605 | 0.0211268 | 0.0413793 | 3 | 139 | 0.978873 | 2.41667 | 2.41667 | 20260522 |
| IC | 0.0330508 | 0.00862069 | 0.017094 | 1 | 115 | 0.991379 | 1.95 | 1.95 | 20260522 |
| CCAM | 0.0288136 | 0.00990099 | 0.0196078 | 1 | 100 | 0.990099 | 1.7 | 1.7 | 20260522 |
| IAC_RDE | 0.0189266 | 0 | 0 | 0 | 67 | 1 | 1.11667 | 1.11667 | 20260522 |
| IRP | 0.00536723 | 0 | 0 | 0 | 19 | 1 | 0.316667 | 0.316667 | 20260522 |

Diffusion order-of-magnitude:
- D ≈ 1e-6 cm²/s (small molecule); t = 6 d = 518400 s
- L_diff = √(D·t) = 0.720 cm ≈ 7.20 mm
- Far (~60 mm) ≫ L_diff: usable as within-plate baseline.

Outputs:
- p2_7_proximity_sensitivity.csv  (both reciprocity, mutual / unidirectional counts)
- p2_7_proximity_density.png
<!-- P2.8 RESULTS -->

## P2.8 Phylogenetic signal (PGLS)

Seed: 20260522. n = 59 strains (see `data/annotation/DATA_NOTES.md` for the one
unmatched strain). Traits: 21 - 12 per-layer out-degrees, 4 category means,
total out-degree, participation coefficient, and PC1-PC3 of the 12-D profile.

### Verdict on h3

**No trait reaches q_BH < 0.05 on any topology** (63 tests). Largest lambda observed across all topologies: **0.952**.
The manuscript's claimed range (0.31-0.48, "weak but significant") is **not reproduced**.

### 16S resolution limit

The 59 strains resolve into only **38 distinct 16S genotypes**. 31 strains fall into
10 groups with pairwise distance exactly zero (group sizes: 6, 6, 3, 3, 3, 2, 2, 2, 2, 2).

Within those groups 16S carries no information, so any phylogenetic-signal test is
underpowered by construction. This is the mechanistic explanation for the null above,
and it is a property of the marker rather than a defect of the data.

### Why the asymptotic test cannot be used here

| topology | n_tips | vcv_rank | tip_depth_ratio | ultrametric |
| --- | --- | --- | --- | --- |
| ML GTR+GAMMA (build_tree.py) | 59 | 44 | 53.0345 | False |
| UPGMA (16S distances) | 59 | 38 | 1 | True |
| WPGMA (16S distances) | 59 | 38 | 1 | True |

The ML tree's VCV is **rank 44 of 59** - a direct consequence of the identical-sequence
groups above, which produce duplicate rows - and its tip depths span a
**53-fold** range, so it is far from ultrametric. Under those conditions the
chi2 likelihood-ratio test is severely anticonservative: randomly shuffled trait
values on this tree routinely yield lambda > 0.9 and LR > 30.

Taking `p_chi2_ref` at face value would report **7 of 63** tests as significant. The
tip-label permutation test - which holds the tree, the conditioning and the trait
marginals fixed and destroys only the tip-to-trait correspondence - shows those
are false positives. **Inference below uses `p_perm`.**

### Pagel's lambda, per trait and topology

| topology | trait | pagel_lambda | LR | p_chi2_ref | p_perm | null_lambda_p95 | q_BH |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ML GTR+GAMMA (build_tree.py) | outdeg_IRP | 0.951552 | 19.2637 | 5.69258e-06 | 0.11 | 0.940388 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | total_out_degree | 0.740413 | 6.75504 | 0.00467416 | 0.175 | 0.853128 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_IAC_RDE | 0.741748 | 9.54072 | 0.00100481 | 0.191 | 0.924958 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | PC1_12D | 0.875898 | 7.21568 | 0.00361347 | 0.2 | 0.907797 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | cat_DA | 0.666428 | 3.10378 | 0.0390554 | 0.295 | 0.886254 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | PC3_12D | 0.650736 | 3.65354 | 0.0279751 | 0.295 | 0.841285 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_RAC | 0.521781 | 2.75297 | 0.0485369 | 0.337 | 0.820176 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | cat_MM | 0.796865 | 1.99878 | 0.0787131 | 0.35 | 0.912665 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_IC | 0.116574 | 1.33512 | 0.123948 | 0.392 | 0.940181 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_IG | 0.648174 | 2.14473 | 0.07153 | 0.4 | 0.883038 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_IRAC | 0.344448 | 0.941257 | 0.165977 | 0.501 | 0.792401 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_CCVM | 0.0680635 | 0.609754 | 0.21744 | 0.52 | 0.825887 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_CS | 0.394964 | 0.552942 | 0.228559 | 0.526 | 0.830933 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_CMP | 0.0877596 | 0.431276 | 0.255682 | 0.539 | 0.863936 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | outdeg_RP | 0.161374 | 0.14229 | 0.353008 | 0.547 | 0.895754 | 0.7658 |
| ML GTR+GAMMA (build_tree.py) | cat_MC | 0.0139571 | 0.0312853 | 0.429803 | 0.618 | 0.879691 | 0.811125 |
| ML GTR+GAMMA (build_tree.py) | participation | 5.62151e-06 | -7.43915e-05 | 1 | 0.708 | 0.843103 | 0.872667 |
| ML GTR+GAMMA (build_tree.py) | PC2_12D | 5.96086e-06 | -0.000113079 | 1 | 0.748 | 0.856951 | 0.872667 |
| ML GTR+GAMMA (build_tree.py) | outdeg_IAC_RAC | 5.96086e-06 | -0.000160622 | 1 | 0.85 | 0.849058 | 0.939474 |
| ML GTR+GAMMA (build_tree.py) | outdeg_CCAM | 5.96086e-06 | -0.000182855 | 1 | 0.903 | 0.866617 | 0.94815 |
| ML GTR+GAMMA (build_tree.py) | cat_IA | 0.485189 | -0.376783 | 1 | 0.974 | 0.834137 | 0.974 |
| UPGMA (16S distances) | cat_DA | 0.057413 | 0.121041 | 0.363954 | 0.061 | 0.15707 | 0.4725 |
| UPGMA (16S distances) | PC2_12D | 0.0978578 | 0.140283 | 0.354 | 0.066 | 0.158203 | 0.4725 |
| UPGMA (16S distances) | outdeg_IC | 0.0664734 | 0.215632 | 0.321194 | 0.079 | 0.421607 | 0.4725 |
| UPGMA (16S distances) | outdeg_IG | 5.96086e-06 | -3.83904e-05 | 1 | 0.09 | 0.182122 | 0.4725 |
| UPGMA (16S distances) | outdeg_CS | 4.42856e-06 | -7.65224e-05 | 1 | 0.142 | 0.135694 | 0.5964 |
| UPGMA (16S distances) | outdeg_RAC | 5.96086e-06 | -0.000113335 | 1 | 0.19 | 0.125811 | 0.665 |
| UPGMA (16S distances) | total_out_degree | 5.96086e-06 | -0.000133752 | 1 | 0.231 | 0.190465 | 0.693 |
| UPGMA (16S distances) | outdeg_CCVM | 5.96086e-06 | -0.000176001 | 1 | 0.355 | 0.145409 | 0.733385 |
| UPGMA (16S distances) | outdeg_IAC_RDE | 5.96086e-06 | -0.000185338 | 1 | 0.37 | 0.352713 | 0.733385 |
| UPGMA (16S distances) | cat_IA | 5.96086e-06 | -0.000184634 | 1 | 0.395 | 0.153808 | 0.733385 |
| UPGMA (16S distances) | PC3_12D | 5.96086e-06 | -0.000191671 | 1 | 0.432 | 0.155134 | 0.733385 |
| UPGMA (16S distances) | PC1_12D | 5.96086e-06 | -0.000193008 | 1 | 0.435 | 0.171433 | 0.733385 |
| UPGMA (16S distances) | outdeg_RP | 5.96086e-06 | -0.000191162 | 1 | 0.454 | 0.200973 | 0.733385 |
| UPGMA (16S distances) | cat_MC | 5.96086e-06 | -0.000209378 | 1 | 0.503 | 0.205993 | 0.7545 |
| UPGMA (16S distances) | outdeg_IRAC | 5.96086e-06 | -0.000219405 | 1 | 0.589 | 0.151491 | 0.784 |
| UPGMA (16S distances) | cat_MM | 5.96086e-06 | -0.00023051 | 1 | 0.659 | 0.413197 | 0.784 |
| UPGMA (16S distances) | outdeg_CMP | 5.96086e-06 | -0.000236182 | 1 | 0.67 | 0.141713 | 0.784 |
| UPGMA (16S distances) | outdeg_IRP | 5.96086e-06 | -0.000238749 | 1 | 0.672 | 0.612938 | 0.784 |
| UPGMA (16S distances) | participation | 5.96086e-06 | -0.000259061 | 1 | 0.829 | 0.155847 | 0.916263 |
| UPGMA (16S distances) | outdeg_IAC_RAC | 5.96086e-06 | -0.000277193 | 1 | 0.913 | 0.167124 | 0.926 |
| UPGMA (16S distances) | outdeg_CCAM | 5.96086e-06 | -0.00027634 | 1 | 0.926 | 0.150776 | 0.926 |
| WPGMA (16S distances) | cat_DA | 0.0584723 | 0.202106 | 0.326513 | 0.065 | 0.185643 | 0.4746 |
| WPGMA (16S distances) | outdeg_IC | 0.0609816 | 0.236684 | 0.313306 | 0.086 | 0.274619 | 0.4746 |
| WPGMA (16S distances) | PC2_12D | 0.0742029 | 0.0298384 | 0.431429 | 0.086 | 0.158785 | 0.4746 |
| WPGMA (16S distances) | outdeg_IG | 5.96086e-06 | -2.93683e-05 | 1 | 0.095 | 0.147975 | 0.4746 |
| WPGMA (16S distances) | outdeg_CS | 3.37415e-06 | -5.74524e-05 | 1 | 0.113 | 0.0940212 | 0.4746 |
| WPGMA (16S distances) | outdeg_RAC | 5.96086e-06 | -0.000112623 | 1 | 0.196 | 0.167687 | 0.654 |
| WPGMA (16S distances) | total_out_degree | 5.96086e-06 | -0.000114189 | 1 | 0.218 | 0.169435 | 0.654 |
| WPGMA (16S distances) | outdeg_CCVM | 5.96086e-06 | -0.000146141 | 1 | 0.299 | 0.191987 | 0.6545 |
| WPGMA (16S distances) | outdeg_IAC_RDE | 5.96086e-06 | -0.000154653 | 1 | 0.322 | 0.211957 | 0.6545 |
| WPGMA (16S distances) | cat_IA | 5.96086e-06 | -0.000165465 | 1 | 0.323 | 0.129446 | 0.6545 |
| WPGMA (16S distances) | PC3_12D | 5.96086e-06 | -0.00016771 | 1 | 0.363 | 0.191235 | 0.6545 |
| WPGMA (16S distances) | PC1_12D | 5.96086e-06 | -0.000170434 | 1 | 0.374 | 0.249972 | 0.6545 |
| WPGMA (16S distances) | cat_MC | 5.96086e-06 | -0.000196924 | 1 | 0.47 | 0.137143 | 0.7125 |
| WPGMA (16S distances) | outdeg_RP | 5.96086e-06 | -0.000189617 | 1 | 0.475 | 0.266174 | 0.7125 |
| WPGMA (16S distances) | outdeg_IRP | 5.96086e-06 | -0.000213707 | 1 | 0.562 | 0.57783 | 0.7665 |
| WPGMA (16S distances) | outdeg_IRAC | 5.96086e-06 | -0.000218486 | 1 | 0.584 | 0.122302 | 0.7665 |
| WPGMA (16S distances) | cat_MM | 5.96086e-06 | -0.000223502 | 1 | 0.648 | 0.292174 | 0.7805 |
| WPGMA (16S distances) | outdeg_CMP | 5.96086e-06 | -0.000231821 | 1 | 0.669 | 0.131996 | 0.7805 |
| WPGMA (16S distances) | participation | 5.96086e-06 | -0.000260752 | 1 | 0.848 | 0.140985 | 0.914 |
| WPGMA (16S distances) | outdeg_IAC_RAC | 5.96086e-06 | -0.000267495 | 1 | 0.882 | 0.131166 | 0.914 |
| WPGMA (16S distances) | outdeg_CCAM | 5.96086e-06 | -0.000273759 | 1 | 0.914 | 0.112705 | 0.914 |

lambda ~ 0 = trait independent of phylogeny; lambda ~ 1 = Brownian heritability.
`p_perm` is from 999 tip-label shuffles; `q_BH` is Benjamini-Hochberg across traits
within each topology. `p_chi2_ref` is the asymptotic value, shown only to document
the miscalibration. `null_lambda_p95` is the 95th percentile of lambda under
shuffling - where it approaches 1, lambda alone carries no evidence.

### Tree-free corroboration (Mantel, 9999 permutations)

| comparison | mantel_rho | p_value |
| --- | --- | --- |
| 12-D out-degree profile (z, Euclidean) | -0.022322 | 0.7773 |
| category DA | 0.0066869 | 0.921 |
| category IA | -0.0550371 | 0.3455 |
| category MM | 0.0366265 | 0.585 |
| category MC | 0.00618726 | 0.9245 |

The Mantel test uses no tree, so agreement with the lambda result rules out the
possibility that the null is an artefact of one tree-building choice.

### Summary for the manuscript

Three independent routes - permutation-calibrated Pagel's lambda on an ML tree, the
same on two ultrametric distance topologies, and a tree-free Mantel test - agree
that interaction profiles carry no detectable 16S phylogenetic signal in these
strains. h3 as stated is not supported.

Outputs: `p2_8_pgls.csv`, `p2_8_mantel.csv`, `p2_8_vcv_diagnostics.csv`.
Figures: `figures/fig_h3_phylogeny.*` (main text), `figures/figS_permutation_calibration.*` (supplement).
<!-- P2.9 RESULTS -->

## P2.9 Spatial decomposition (3 configurations)

Per-configuration matrices missing for: direct, indirect, distant

Expected:
- data/matrices_xls_direct/{layer}.xlsx
- data/matrices_xls_indirect/{layer}.xlsx
- data/matrices_xls_distant/{layer}.xlsx

12 xlsx per directory, same filenames as data/matrices_xls/, scored per spatial configuration.

When present: per-(layer, config) density, reciprocity (pairs), SCC; J(direct,indirect), J(indirect,distant), J(direct,distant); spatial-sensitivity = (1 − J(direct,distant)) × mean density. Outputs p2_9_spatial_decomposition.csv + figure.
<!-- P2.10 RESULTS -->

## P2.10 Direction-convention audit

Convention: row i = receiver, col j = sender. M[i,j]=1 → j elicited response l in i. NetworkX is fed M.T so A[u,v] = u→v.

Asymmetry stats below. mean_out vs mean_in tells dominant scored direction; large max_out_in_imbalance = a few hub emitters/receivers.

| phenotype | convention_anchor | n | edges | mean_out_degree | mean_in_degree | std_out_degree | std_in_degree | max_out_in_imbalance | mean_out_in_imbalance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CCAM | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 102 | 1.7 | 1.7 | 1.5308 | 4.47325 | 29 | 2.46667 |
| CCVM | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 252 | 4.2 | 4.2 | 2.06398 | 7.55601 | 27 | 5.5 |
| CMP | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 260 | 4.33333 | 4.33333 | 3.35989 | 4.94188 | 13 | 3.73333 |
| CS | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 265 | 4.41667 | 4.41667 | 1.84654 | 8.2144 | 41 | 5.53333 |
| IG | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 767 | 12.7833 | 12.7833 | 10.9089 | 7.89956 | 40 | 10.3333 |
| IC | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 117 | 1.95 | 1.95 | 4.51451 | 1.55376 | 18 | 3.06667 |
| RP | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 145 | 2.41667 | 2.41667 | 1.48651 | 8.71262 | 50 | 4.16667 |
| IRP | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 19 | 0.316667 | 0.316667 | 0.718602 | 1.2582 | 9 | 0.633333 |
| IAC_RAC | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 383 | 6.38333 | 6.38333 | 3.31709 | 7.04294 | 21 | 5.56667 |
| IAC_RDE | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 67 | 1.11667 | 1.11667 | 2.19159 | 2.21428 | 13 | 2 |
| IRAC | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 321 | 5.35 | 5.35 | 2.95423 | 6.06582 | 17 | 5.26667 |
| RAC | rows=receiver, cols=sender (data_io.load_tensor and root/_notebook.layer_stats consistent) | 60 | 1677 | 27.95 | 27.95 | 3.02999 | 23.4061 | 37 | 23 |

Outputs: p2_10_directionality_audit.csv (author-fill phenotype_scoring_rule)

Action: document bench direction per layer before submission.
<!-- P2.11 RESULTS -->

## P2.11 Dominance: per-layer and cross-layer

Seed: 20260522. Per-layer DS and BT; cross-layer Spearman; PCA on z-scored 60×12 DS.

- mean off-diag Spearman ρ = 0.063
- median = 0.043
- PCs to 80% var: 6 / 12

| PC | variance_explained | cumulative_variance |
| --- | --- | --- |
| PC1 | 0.242691 | 0.242691 |
| PC2 | 0.169049 | 0.411739 |
| PC3 | 0.141489 | 0.553229 |
| PC4 | 0.122937 | 0.676166 |
| PC5 | 0.10033 | 0.776495 |
| PC6 | 0.0682137 | 0.844709 |
| PC7 | 0.0450018 | 0.889711 |
| PC8 | 0.0373663 | 0.927077 |

ρ→0 and high n_pcs_80 → multidimensional dominance (h2). ρ→1, n_pcs_80=1 → one axis.

Outputs:
- p2_11_davids_per_layer.csv
- p2_11_bt_strengths_per_layer.csv
- p2_11_cross_layer_dominance_corr.csv
- p2_11_dominance.png
<!-- P2.12 RESULTS -->

## P2.12 Metabolic niche overlap and layer-specific interactions

Seed: 20260522. n = 60 / 60 strains (8 carbon sources; the mineral-medium
control column `MM` is excluded from the niche profile). Assimilation is scored
0-4, so the graded row sum is used as capacity: mean 12.1, range 2-23 of 32.

Presence/absence breadth is **saturated** (mean 7.33 of 8 sugars; 54 of 60 strains
use at least 7), which is why the graded scores are used throughout - binarising
discards nearly all of the between-strain variation.

### Strain level: no association

| response | spearman_rho | p_value | q_BH |
| --- | --- | --- | --- |
| MM mean out-degree | -0.198701 | 0.128016 | 0.730611 |
| participation coefficient | 0.174082 | 0.183435 | 0.730611 |
| IA mean out-degree | 0.166291 | 0.204137 | 0.730611 |
| out-degree IRP | -0.159361 | 0.223898 | 0.730611 |
| out-degree RP | -0.155098 | 0.236695 | 0.730611 |
| out-degree IRAC | 0.152884 | 0.243537 | 0.730611 |
| out-degree IAC_RAC | 0.10156 | 0.440036 | 0.960118 |
| out-degree IC | 0.0790062 | 0.548476 | 0.960118 |
| out-degree CS | 0.069993 | 0.595134 | 0.960118 |
| out-degree CMP | -0.0610679 | 0.642997 | 0.960118 |
| out-degree CCAM | -0.0517362 | 0.69463 | 0.960118 |
| out-degree IG | 0.0284965 | 0.828884 | 0.960118 |
| out-degree CCVM | 0.0283397 | 0.829811 | 0.960118 |
| out-degree IAC_RDE | 0.0233723 | 0.859307 | 0.960118 |
| MC mean out-degree | 0.0191389 | 0.884597 | 0.960118 |
| out-degree RAC | -0.0139473 | 0.915767 | 0.960118 |
| total out-degree | 0.00727089 | 0.95603 | 0.960118 |
| DA mean out-degree | 0.00659429 | 0.960118 | 0.960118 |

**0 of 18** associations reach q_BH < 0.05. Metabolic capacity does not predict how
active a strain is (total out-degree rho = +0.007) nor how evenly its interactions
spread across layers. The "metabolism shapes network position" reading is **not**
supported at the strain level.

### Dyad level: suggestive, and in the unexpected direction

| layer | mantel_rho | p_value | n_dyads | q_BH |
| --- | --- | --- | --- | --- |
| IC | 0.12148 | 0.01 | 116 | 0.0852 |
| IG | 0.131461 | 0.0142 | 667 | 0.0852 |
| IRP | 0.0589023 | 0.0998 | 19 | 0.3024 |
| RAC | 0.106127 | 0.1008 | 1293 | 0.3024 |
| IRAC | -0.0614447 | 0.138 | 292 | 0.3224 |
| CCAM | 0.0679218 | 0.1612 | 101 | 0.3224 |
| CS | -0.0633817 | 0.2394 | 250 | 0.4104 |
| IAC_RDE | -0.0336771 | 0.4116 | 67 | 0.546 |
| IAC_RAC | -0.0366641 | 0.4252 | 349 | 0.546 |
| RP | 0.0610252 | 0.455 | 142 | 0.546 |
| CMP | 0.0232456 | 0.6216 | 249 | 0.678109 |
| CCVM | 0.000648724 | 0.9936 | 245 | 0.9936 |

**0 of 12** layers reach q_BH < 0.05. The two strongest effects - IC
(rho = +0.121, q = 0.085) and IG (rho = +0.131, q = 0.085) - are both
direct-antagonism layers and both **positive**, meaning metabolically *dissimilar*
pairs antagonise each other more often.

That is the opposite of the competition-relatedness expectation, under which
resource-similar strains should compete hardest. It is reported as a
hypothesis-generating observation only: neither layer clears FDR, and 2 of 12
marginal results is a weak basis for a mechanistic claim.

Outputs: `p2_12_metabolism_strain.csv`, `p2_12_metabolism_dyad.csv`,
`p2_12_metabolism_profiles.csv`. Figure: `figures/fig_metabolism.*` (main text).
<!-- P2.13 RESULTS -->

## P2.13 Taxon sensitivity (n=56, Streptomyces only)

Seed: 20260522. Excluded (4): MS10 8, MS3 15, MS3 18, MS53 7.
n = 60 -> 56 strains.

Overlap counts are reported at two draw counts. The published P2.4 result uses
10,000 null draws; the n=56 arm uses 2,000. Comparing those directly would charge
any draw-count effect to the excluded strains, so an n=60 result at the matched
2,000 draws is computed as the honest comparator (8 at 10,000 -> 8 at 2,000).
**Read the n=56 column against the matched row, not the published one.**

Every qualitative conclusion holds. The only change is that 2 of 8 over-overlapping pairs (CCAM-CS, CMP-CS) fall below the FDR threshold at n=56.

### Headline quantities

| quantity | n60 | n56 |
| --- | --- | --- |
| Jaccard null: over-overlapping pairs, published 10,000 draws | 8 | - |
| Jaccard null: over-overlapping pairs, matched 2,000 draws | 8 | 6 |
| Jaccard null: under-overlapping pairs, matched 2,000 draws | 1 | 1 |
| 12-D PCA PC1 (% variance) | 32.5 | 31.2 |
| 12-D PCA PC2 (% variance) | 17.3 | 15.2 |
| PCs above Horn parallel-analysis null | 3 | 3 |
| GMM components by 10-fold CV log-likelihood (12-D, full cov) | 1 | 1 |
| Layer density rank correlation n60 vs n56 (Spearman) | - | 1 |
| Largest-SCC rank correlation n60 vs n56 (Spearman) | - | 0.9982 |

### Per-layer structure

| layer | density_n60 | reciprocity_pairs_n60 | largest_scc_n60 | density_n56 | reciprocity_pairs_n56 | largest_scc_n56 | d_density | d_scc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CCAM | 0.0288 | 0.0099 | 4 | 0.0305 | 0.0108 | 4 | 0.0017 | 0 |
| CCVM | 0.0712 | 0.0286 | 24 | 0.0718 | 0.0279 | 23 | 0.0006 | -1 |
| CMP | 0.0734 | 0.0442 | 39 | 0.0724 | 0.0519 | 36 | -0.001 | -3 |
| CS | 0.0749 | 0.06 | 26 | 0.0766 | 0.0583 | 24 | 0.0018 | -2 |
| IG | 0.2167 | 0.1499 | 59 | 0.2214 | 0.1348 | 55 | 0.0048 | -4 |
| IC | 0.0331 | 0.0086 | 2 | 0.0334 | 0.0098 | 2 | 0.0004 | 0 |
| RP | 0.041 | 0.0211 | 3 | 0.0373 | 0.0177 | 2 | -0.0036 | -1 |
| IRP | 0.0054 | 0 | 1 | 0.0049 | 0 | 1 | -0.0005 | 0 |
| IAC_RAC | 0.1082 | 0.0974 | 46 | 0.1058 | 0.0976 | 42 | -0.0023 | -4 |
| IAC_RDE | 0.0189 | 0 | 1 | 0.0172 | 0 | 1 | -0.0017 | 0 |
| IRAC | 0.0907 | 0.0993 | 41 | 0.0896 | 0.0952 | 38 | -0.0011 | -3 |
| RAC | 0.4737 | 0.297 | 51 | 0.4841 | 0.3067 | 46 | 0.0104 | -5 |

Layer density and largest-SCC orderings are preserved almost exactly
(Spearman rho = 1.000 and 0.998), so the architecture contrasts that
the Results are built on do not depend on the four non-Streptomyces isolates.

Pairs losing significance at n=56: CCAM-CS, CMP-CS


The taxonomy wording in the title, abstract and Methods still needs correcting -
see `data/annotation/DATA_NOTES.md`.

Outputs: `p2_13_headline_n56.csv`, `p2_13_layer_structure_n56.csv`
