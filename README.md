# Multilayer interaction networks in a soil *Streptomyces* community

Analysis code and data for the manuscript. Every number, table and figure the paper
reports is produced here from the raw inputs, which are also here.

```bash
pip install -r requirements.txt
python run_all.py
```

Everything installs from PyPI — no external binaries, no `sudo`. A cold run
regenerates the degree-preserving null models, which dominates the runtime (~2 h);
afterwards they are cached in `.cache/` (git-ignored) and reruns take minutes.

Individual stages are runnable on their own:

```bash
cd code && python p2_8_phylogeny_pgls.py
```

## Reproducibility

Seed `20260522`, set in `code/data_io.py`. Every randomised procedure — null draws,
tip-label permutations, bootstrap resamples, parallel-analysis permutations, mixture
initialisation — derives from it.

`code/verify_reported_values.py` re-derives all **121** quantities the manuscript
states from the regenerated CSVs and exits non-zero on any mismatch. It runs as the
last stage of `run_all.py`. If a number in the paper and a number in the analysis ever
disagree, this is what says so.

## Layout

| Path | Contents |
|---|---|
| `code/` | every script, exactly once |
| `data/matrices_xls/` | the 12 source phenotype matrices (60 × 60, directed) |
| `data/annotation/` | 16S sequences, pairwise distances, carbohydrate profiles, GenBank submission record, `DATA_NOTES.md` |
| `data/phylogeny/` | derived: cleaned FASTA, alignment, ML tree, provenance |
| `supplementary/tables/` | generated CSVs |
| `supplementary/figures/` | generated figures |
| `run_all.py` | single entry point |
| `INDEX.md` | claim → script → output map |

**Read `data/annotation/DATA_NOTES.md` before interpreting any phylogenetic result.**

## Pipeline

| Stage | Produces |
|---|---|
| `prep_phylogeny.py` | cleaned 16S FASTA + per-record label audit |
| `build_tree.py` | MAFFT → trimAl → VeryFastTree, GTR+GAMMA ML tree |
| `p2_2_gmm.py` | cluster-vs-continuum model selection |
| `p2_3_pca12.py` | 12-D PCA, bootstrap CIs, Horn parallel analysis |
| `p2_4_jaccard_null.py` | layer-overlap null, BY-FDR |
| `p2_5_scc_null.py` | strongly-connected-component null |
| `p2_7_proximity_sensitivity.py` | per-layer density, both reciprocity conventions |
| `p2_8_phylogeny_pgls.py` | phylogenetic signal, permutation-based |
| `p2_9_spatial_decomposition.py` | *needs per-configuration matrices; skips otherwise* |
| `p2_10_directionality.py` | direction-convention audit |
| `p2_11_dominance.py` | David's score, Bradley–Terry, cross-layer dominance |
| `p2_12_metabolism.py` | carbohydrate assimilation vs network position |
| `p2_13_taxon_sensitivity.py` | n = 56 *Streptomyces*-only robustness check |
| `p2_14_reciprocity_transitivity_null.py` | reciprocity and transitivity nulls |
| `generate_figs.py`, `make_paper_figures.py` | figure sets |
| `figure_2_spatial.py`, `figure_3_overlap.py`, `figure_4_strategy.py` | main-text figures 2–4 |
| `tables_1_2.py` | main-text tables 1 and 2 |
| `verify_reported_values.py` | re-derives every reported number |

A stage whose input data is absent records a skip in `supplementary/report.md` and the
run continues; the summary at the end lists what ran and what did not.

## Conventions

- **Binary by design.** All per-layer analyses use the binary projection
  `A_l = 1[M_l > 0]`. The layer constants `c_l` (1, 2, 3, 5) are severity tags, not
  per-edge weights, and enter only the aggregate visualisation.
- **Reciprocity.** Pairs-based `R_p = mutual / connected dyads` throughout; the
  edges-based variant is tabulated alongside for cross-reference.
- **Null model.** 10,000 directed degree-preserving edge-swap draws per layer
  (`networkx.algorithms.swap`), Phipson–Smyth (M+1)/(N+1) p-values, Benjamini–Yekutieli
  FDR by default with BH in a side column. The same draws serve as the null for layer
  overlap, SCC size, reciprocity and transitivity.
- **Confound pair.** `{CCAM, CCVM}` is corrected within its family but reported
  separately from the biological summary: both colours are scored from the same
  colonies.
- **Model selection.** k = 1..8 × {full, diag, spherical} × {BIC, ICL, 10-fold CV,
  silhouette}. Cross-validated log-likelihood is the criterion of record, and the
  selected order is reported per covariance parameterization — it is not stable across
  them.
- **Phylogenetic inference is permutation-based, never asymptotic.** On these data the
  χ² likelihood-ratio test is strongly anticonservative: it returns 7 of 63 tests as
  significant with λ up to 0.95, while shuffled traits on the same tree reproduce
  λ > 0.9. The signal comes from a rank-deficient, strongly non-ultrametric covariance
  structure, not from heritability. `p_chi2_ref` is retained in the CSV only to
  document the discrepancy.
- **Tree building.** RAxML has no pip distribution, so the tree is built with MAFFT +
  trimAl + VeryFastTree (a validated FastTree-2 equivalent).
  `data/phylogeny/tree_provenance.txt` records exactly what ran, and those are the tool
  names that belong in the Methods. Placing a RAxML tree at `data/phylogeny/tree.nwk`
  makes `build_tree.py` leave it alone and `p2_8` use it unchanged.

## Known limitations

1. **`S705A` appears twice in the 16S FASTA** (1390 bp and 1468 bp) under one label.
   `prep_phylogeny.py` keeps the longer record.
2. **Four isolates are not *Streptomyces*** (`MS53 7`, `MS10 8`, `MS3 15`, `MS3 18`).
   `p2_13` shows no conclusion depends on them.
