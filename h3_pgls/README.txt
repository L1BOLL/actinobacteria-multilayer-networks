h3 — Phylogenetic generalized least squares

v15 abstract: Pagel's λ = 0.31–0.48 per trait. Number from PI's prior work;
this bundle does not reproduce it (tree not committed).

Script: p2_8_phylogeny_pgls.py
  per-trait Pagel's λ, LR test vs λ=0
  auto-skips if data/phylogeny/tree.nwk absent

Required input: Newick tree at data/phylogeny/tree.nwk; leaves match 60 strain IDs.

Output (when input present):
  phase2_outputs/tables/p2_8_pgls.csv
  (trait, lambda_ml, lambda_lo, lambda_hi, ll_lambda, ll_zero, lrt_p)

Before submission: rerun with tree (replace 0.31–0.48), or drop the abstract sentence.
