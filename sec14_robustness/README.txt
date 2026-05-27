§14 — Robustness

1. 4-D PC1 ↔ 12-D PC1: Spearman ρ = 0.934, p < 10⁻²⁶
2. GMM k=1–3 across all 6 (embedding × cov) configs (CV-LL); dip test p > 0.5
3. Jaccard/SCC: 10,000-draw null; BH and BY agree qualitatively (BY default)
4. Supra-adjacency ω ∈ [0.1, 1.0]: participation variation < 10%
   (no main-text result depends on supra-adjacency)

v14 "weight thresholding" check deleted — data is {0, c_l}.

Scripts (covered earlier, re-linked here):
  p2_3_pca12.py
  p2_2_gmm.py
  p2_4_jaccard_null.py
  p2_7_proximity_sensitivity.py

Outputs: extra columns in §8 / §11 CSVs.
