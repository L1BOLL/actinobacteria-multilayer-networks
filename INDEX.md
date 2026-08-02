# Claim → script → output

Every number the manuscript states, and where it comes from. Paths are relative
to the repository root. Regenerate everything with `python run_all.py`.

Section numbers refer to the manuscript; scripts live once in `code/` and several
serve more than one section.

## Section map

| § | Claim area | Script(s) | Output |
|---|---|---|---|
| 4 | Network representation | `data_io.py` | loader, no output |
| 5 | Network analysis utilities | `data_io.py`, `null_models.py` | — |
| 6 | Dominance methods | `p2_11_dominance.py` | `p2_11_davids_*`, `p2_11_bt_*` |
| 7 | Density and reciprocity | `p2_7_proximity_sensitivity.py` | `p2_7_proximity_sensitivity.csv` |
| 8 | Layer overlap and SCC nulls | `p2_4_jaccard_null.py`, `p2_5_scc_null.py` | `p2_4_*`, `p2_5_*` |
| 9 | Topology, direction convention | `p2_10_directionality.py`, `p2_7_*` | `p2_10_*.csv`, `p2_7_proximity_density.png` |
| 10 | Strain strategy | `generate_figs.py`, `p2_7_*` | `figure3_*` |
| 11 | PCA and mixture models | `p2_3_pca12.py`, `p2_2_gmm.py`, `generate_figs.py` | `p2_3_*`, `p2_2_*`, `figure4_*` |
| 12 | Cross-layer dominance | `p2_11_dominance.py` | `p2_11_*.csv`, `p2_11_dominance.png` |
| 13 | Spatial decomposition | `p2_9_spatial_decomposition.py` | *skips without per-config data* |
| 14 | Robustness | `p2_2`, `p2_3`, `p2_4`, `p2_7` | columns inside §8 / §11 CSVs |
| h3 | Phylogenetic signal | `prep_phylogeny.py`, `build_tree.py`, `p2_8_phylogeny_pgls.py` | `p2_8_*`, `p2_8_phylo_signal.png` |
| — | Metabolism | `p2_12_metabolism.py` | `p2_12_*`, `p2_12_metabolism.png` |
| — | Taxon sensitivity | `p2_13_taxon_sensitivity.py` | `p2_13_*` |

---

## Claim ↔ number ↔ output

### §8 layer overlap (10,000 null draws, BY-FDR)

| Claim | Source |
|---|---|
| mean J = 0.047 across 66 pairs | `p2_4_jaccard_null.csv`, `mean(J_obs)` |
| CCAM/CCVM confound: J = 0.336, z = +16.2, q_BY = 0.005 | same, that row |
| 8 over-overlapping biological pairs: IRAC-RAC (+9.2), CCVM-CMP (+8.2), CCAM-CMP (+7.7), CCVM-CS (+6.4), IG-IRAC (+4.8), CCAM-CS (+3.9), CMP-CS (+3.7), CS-RAC (+3.7) | same, `significant_by & z>0 & ~confound` |
| 1 under-overlapping pair: IAC_RAC-RAC (z = −5.3) | same, `significant_by & z<0` |

### §8 SCC null

| Claim | Source |
|---|---|
| \|z\| < 1.5 in every layer; no layer significant after BY-FDR | `p2_5_scc_null.csv` |

### §9 density and reciprocity

| Claim | Source |
|---|---|
| ρ: RAC 0.474, IG 0.217, IAC_RAC 0.108, IRAC 0.091, min IRP 0.005 | `p2_7_proximity_sensitivity.csv` (`density`) |
| R_p: RAC 0.297, IG 0.150, all others < 0.10 | same (`reciprocity_pairs`) |
| R_e = 2R_p/(1+R_p), identical ordering | same (`reciprocity_edges`) |
| mutual pairs: RAC 384, IG 100 | same (`mutual_pairs`) |
| largest SCC: IG 59, RAC 51, IAC_RAC 46, IRAC 41, CMP 39, CS 26, CCVM 24, CCAM 4, RP 3, IC 2, IRP 1, IAC_RDE 1 | same (`largest_scc`) |

### §10 strain strategy

| Claim | Source |
|---|---|
| P_i ∈ [0.60, 0.84] for all 60 strains — the generalist/specialist dichotomy does not separate | `data/ecological_identity_with_PCs_and_clusters.csv` |

### §11 PCA and mixture models

| Claim | Source |
|---|---|
| 12-D PCA: PC1 32.5% [24.2, 43.7], PC2 17.3% [13.4, 22.1], PC3 11.8% [9.8, 15.0]; PC1+PC2 = 49.8% | `p2_3_pca12_variance.csv` |
| Horn parallel analysis: 3 components above the null (PC4 falls below) | same (`significant_vs_parallel`) |
| 4-D PCA: PC1 53.0%, PC2 23.9%, PC1+PC2 = 76.8% | `p2_3_pca4_variance.csv` |
| 12-D PC1 ↔ 4-D PC1 Spearman ρ = 0.934, p < 10⁻²⁶ | `p2_3_pc1_spearman.csv` |
| PC1 loadings: CMP 0.41, IRP 0.40, CCVM 0.38, CCAM 0.37, IG 0.33 | `p2_3_pca12_loadings.csv` |
| PC2: RAC +0.54, IC −0.38, IG −0.38. PC3: IAC_RDE +0.65 | same |
| GMM CV log-likelihood: k = 1 (full cov); k = 2–3 (diag, spherical) | `p2_2_gmm_modelselect.csv` (`best_k_cv`) |
| BIC/ICL prefer larger k — small-N over-selection | same |
| Dip test on PC1/PC2, both embeddings: p > 0.5 | `p2_2_diptest.csv` |

> **The 4-D PCA is computed on z-standardized profiles**, as the Methods specify.
> Running it on the raw embedding gives 54.2% / 23.3%, which is where the earlier
> "54.0% / 23.7%" figure axes came from. `generate_figs.py` standardizes at both
> call sites through `embedding_pca()` so the figures and the text cannot diverge.

### §12 cross-layer dominance

| Claim | Source |
|---|---|
| mean Spearman of David's score = 0.063 (median 0.043, range [−0.66, +0.62]) | `p2_11_cross_layer_dominance_corr.csv` |
| 7 of 66 pairs \|ρ\| > 0.5 | same |
| strongest anti-correlations: IG/RAC −0.66, IAC_RDE/RAC −0.63, IAC_RDE/IRAC −0.58 | same |
| dominance PCA: PC1 24.3%, PC1+PC2 41.2%; 3 PCs for 50%, 6 for 80% | `supplementary/report.md` §P2.11 |
| Bradley–Terry and David's score concordant, ρ = 0.59–0.99 (mean 0.87) | `p2_11_bt_strengths_per_layer.csv` |

### h3 phylogenetic signal

| Claim | Source |
|---|---|
| 0 of 21 traits significant on the ML tree (999 tip-label permutations, BH) | `p2_8_pgls.csv` |
| 0 of 42 on UPGMA and WPGMA topologies | same |
| smallest q across all 63 trait × topology tests = 0.47 | same |
| Mantel, 16S distance vs 12-D profile: ρ = −0.022, p = 0.78 | `p2_8_mantel.csv` |
| per-category Mantel: \|ρ\| ≤ 0.055, all p > 0.34 | same |
| ML-tree VCV rank 44 of 59; tip depths span 53× | `p2_8_vcv_diagnostics.csv` |
| 60 isolates → 38 distinct 16S genotypes; 31 strains in 11 zero-distance groups | `supplementary/report.md` §P2.8 |

> **The asymptotic χ² test is not usable on these data.** It returns 7 of 63 tests
> as significant with λ up to 0.95, but shuffled trait values on the same tree
> reproduce λ > 0.9 and LR > 30 — the signal comes from a rank-deficient,
> strongly non-ultrametric covariance structure, not from heritability.
> Inference uses the permutation p-value (`p_perm`); `p_chi2_ref` is retained in
> the CSV only to document the discrepancy.

### Metabolism

| Claim | Source |
|---|---|
| assimilation capacity mean 12.1, range 2–23 of 32 | `p2_12_metabolism_profiles.csv` |
| 0 of 18 strain-level associations at q < 0.05; total out-degree ρ = 0.007 | `p2_12_metabolism_strain.csv` |
| presence/absence breadth saturated: 54 of 60 strains use ≥ 7 of 8 sugars | same |
| dyad level: IC ρ = +0.121 (p = 0.010), IG ρ = +0.131 (p = 0.014), both q = 0.085 | `p2_12_metabolism_dyad.csv` |

> Metabolic *dissimilarity* predicting antagonism runs opposite to the
> competition–relatedness expectation, but neither layer clears FDR. Reported as
> hypothesis-generating only.

### Taxon sensitivity (n = 56, *Streptomyces* only)

| Claim | Source |
|---|---|
| layer density ordering preserved, Spearman ρ = 1.000 | `p2_13_headline_n56.csv` |
| largest-SCC ordering preserved, ρ = 0.998 | same |
| components above parallel-analysis null: 3 → 3; CV-optimal k: 1 → 1 | same |
| 2 of 8 over-overlapping pairs lost at n = 56 (CCAM-CS, CMP-CS) | same, and `p2_13_layer_structure_n56.csv` |

> The two lost pairs were the weakest in the full cohort (both q = 0.019, the
> joint-largest of the nine significant pairs); all six pairs at q = 0.005,
> including every antagonistic-axis pair, are unaffected. The n = 60 baseline is
> reported at both 10,000 and 2,000 draws, because the n = 56 arm uses 2,000 and
> the comparison must be draw-matched.

---

## Conventions

- **Binary by design.** All per-layer analyses use `A_l = 1[M_l > 0]`. The layer
  constants `c_l` (1, 2, 3, 5) are severity tags, not per-edge weights, and enter
  only the aggregate visualisation matrix.
- **Reciprocity.** Pairs-based `R_p = mutual_pairs / connected_pairs` throughout;
  the edges-based variant is tabulated alongside for cross-reference.
- **Null model.** 10,000 directed degree-preserving edge-swap draws per layer
  (`networkx.algorithms.swap`), Phipson–Smyth (M+1)/(N+1) p-values, BY-FDR by
  default with BH in a side column.
- **Confound pair.** `{CCAM, CCVM}` reported separately from the biological summary.
- **Model selection.** k = 1..8 × {full, diag, spherical} × {BIC, ICL, 10-fold CV,
  silhouette}; the cross-validated log-likelihood is the criterion of record.
- **Phylogenetic inference.** Permutation-based, never asymptotic — see the h3
  note above.
- **Seed.** `20260522`, in `code/data_io.py`.

---

## Outstanding

1. **`S1226` / `S1430`** — one strain is in the interaction matrices, the other in
   both 16S sources; all 59 others match. Phylogenetic analyses run at n = 59
   pending confirmation. See `data/annotation/DATA_NOTES.md`.
2. **Duplicate `S705A`** in the FASTA (1390 bp and 1468 bp records); the longer is
   kept. If these are two isolates, one needs relabelling before deposition.
3. **Per-configuration matrices** for §13 — supply `data/matrices_xls_direct/`,
   `_indirect/`, `_distant/` and `p2_9_spatial_decomposition.py` runs automatically.
4. **Original RAxML tree**, if it exists — drop at `data/phylogeny/tree.nwk` and
   re-run; the Methods tool names would revert accordingly.
5. **Four non-*Streptomyces* isolates** (`MS53 7`, `MS10 8`, `MS3 15`, `MS3 18`) —
   the manuscript's cohort wording needs to match; `p2_13` shows no conclusion
   depends on them.
