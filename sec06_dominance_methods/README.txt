§6 — Methods — cross-layer dominance

Win matrix: W_l[i,j] = A_l[j,i]  (sender = winner).

David's score (De Vries 2006):
  P_ij = W[i,j] / (W[i,j] + W[j,i])
  w_i  = Σ_j P_ij,         l_i  = Σ_j P_ji
  w2_i = Σ_j P_ij · w_j,   l2_i = Σ_j P_ji · l_j
  DS_i = w_i + w2_i − l_i − l2_i

Bradley–Terry: π_i by iterative ML under P(i ≻ j) = π_i / (π_i + π_j).

Output: 60×12 DS and BT matrices.
Then: 66 pairwise Spearman of per-strain DS across layers; PCA on z-scored 60×12 DS.

Script: p2_11_dominance.py
Outputs:
  phase2_outputs/tables/p2_11_davids_per_layer.csv
  phase2_outputs/tables/p2_11_bt_strengths_per_layer.csv
  phase2_outputs/tables/p2_11_cross_layer_dominance_corr.csv
  phase2_outputs/figures/p2_11_dominance.png
