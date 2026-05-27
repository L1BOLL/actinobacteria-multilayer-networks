§12 — Cross-layer dominance (NEW subsection)

mean off-diagonal Spearman of DS = 0.063  (median 0.043, range [−0.659, +0.618])
7 of 66 pairs |ρ| > 0.5

Strongest anti-correlations:
  IG / RAC          −0.66
  IAC_RDE / RAC     −0.63
  IAC_RDE / IRAC    −0.58

Dominance-matrix PCA:
  PC1 = 24.3%, PC1+PC2 = 41.2%
  3 PCs for 50% variance, 6 PCs for 80%

Bradley–Terry: qualitatively identical.

Conclusion: no universal dominance axis; dominance is layer-specific.
Independent confirmation of h2 without embedding/clustering.

Script: p2_11_dominance.py
Outputs:
  phase2_outputs/tables/p2_11_davids_per_layer.csv         (60×12)
  phase2_outputs/tables/p2_11_bt_strengths_per_layer.csv   (60×12)
  phase2_outputs/tables/p2_11_cross_layer_dominance_corr.csv  (66 rows)
  phase2_outputs/figures/p2_11_dominance.png
  PCA variance/scree numbers: phase2_report.md
