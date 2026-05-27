§11 — Ecological embedding analysis, PCA, GMM

4-D PCA (visualization only):
  PC1 = 54.0%, PC1+PC2 = 77.7%

12-D out-degree PCA (honest dimensionality):
  PC1 = 32.5%  [0.242, 0.437]
  PC2 = 17.3%  [0.134, 0.221]
  PC3 = 11.8%  [0.098, 0.150]
  PC1+PC2 = 49.8%
  Horn parallel analysis: 3 components above null

12-D PC1 ↔ 4-D PC1: Spearman ρ = 0.934, p < 10⁻²⁶

12-D loadings:
  PC1+ CCAM, CCVM, CMP, CS, IG, IC, RP, IRP
  PC1− IAC_RAC, IRAC, IAC_RDE
  PC2+ RAC, RP, IAC_RAC, IRAC
  PC2− IG, IC
  PC3+ IAC_RDE (0.65)

GMM (k=1..8, cov ∈ {full, diag, spherical}, both embeddings):
  10-fold CV log-likelihood → k=1 full, k=2–3 diag/sph (every embedding)
  BIC, ICL → larger k (small-N over-selection)
  Dip test PC1, PC2 both embeddings → p > 0.5

Conclusion: continuous variation in ~3-D subspace; at most weak local 2–3 component
structure; no discrete ecotypes. (Replaces v14 k=8 full-cov BIC overfit.)

Scripts:
  p2_3_pca12.py     12-D PCA, 2000 boot CIs, 1000-draw Horn
  p2_2_gmm.py       3 cov × 4 criteria × k=1..8
  generate_figs.py  Fig. 4 panels

Outputs:
  p2_3_pca12_variance.csv, p2_3_pca4_variance.csv
  p2_3_pca12_loadings.csv, p2_3_pc1_spearman.csv
  p2_2_gmm_modelselect.csv, p2_2_diptest.csv, p2_2_summary.csv, p2_2_cluster_*.csv
  p2_3_pca12_scree.png, p2_3_pca12_biplot.png, p2_2_bic_curves.png
  figure_outputs/figure4_ecological_strategy_space.{png,svg}

Caveat: figure colours strains at k=4; CV says k=1–3. Re-render at k=3 or
flag k=4 as visualization-only.
