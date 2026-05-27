§5 — Methods — networks analysis

Per-layer on A_l:
  ρ_l = |E_l| / (N(N−1))
  mutual pairs, R_p = mutual_pairs / connected_pairs   (main)
  R_e = 2 R_p / (1 + R_p)                              (secondary)
  SCC: Kosaraju
  out-degree, P_i over 12 layers
  betweenness B_v = Σ σ_ij(v)/σ_ij
  Jaccard J_ab = |E_a ∩ E_b| / |E_a ∪ E_b|

Null: 10,000 directed degree-preserving swaps/layer; Phipson–Smyth (M+1)/(N+1);
BY-FDR (main) across 66 pairs, BH side column.
{CCAM, CCVM} flagged as methodological confound.

Scripts (utilities for §8, §9):
  data_io.py
  null_models.py

Output: none directly.
