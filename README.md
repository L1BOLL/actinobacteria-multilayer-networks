# BN_paper_methods

Analysis code and supplementary materials for the multilayer *Streptomyces*
interaction-network manuscript. This repository is the single source of truth:
every number, table and figure in the paper is produced here, from the raw inputs
that are also here.

## Run

```bash
pip install -r requirements.txt
python run_all.py
```

Everything installs from PyPI — no external binaries, no `sudo`. A cold run
regenerates the degree-preserving null models, which dominates the runtime;
afterwards they are cached in `.cache/` (git-ignored) and reruns are fast.

Individual stages are runnable on their own:

```bash
cd code && python p2_8_phylogeny_pgls.py
```

## Layout

```
code/           every script, exactly once
data/           raw inputs + derived per-strain tables
supplementary/  everything the manuscript cites: report.md, tables/, figures/
run_all.py      single entry point
INDEX.md        claim -> script -> output map
```

`data/`:

| Path | Contents |
|---|---|
| `matrices_xls/` | the 12 source phenotype matrices (60 × 60, directed) |
| `annotation/` | 16S sequences, pairwise 16S distances, carbohydrate profiles, and `DATA_NOTES.md` |
| `phylogeny/` | derived: cleaned FASTA, alignment, ML tree, provenance |

**Read `data/annotation/DATA_NOTES.md` before interpreting any phylogenetic
result.** One strain is unmatched between the bench and sequencing datasets, and
four of the 60 isolates are not *Streptomyces*.

## Pipeline

| Stage | Produces |
|---|---|
| `prep_phylogeny.py` | cleaned 16S FASTA + per-record label audit |
| `build_tree.py` | MAFFT alignment → trimAl → GTR+GAMMA ML tree |
| `p2_2_gmm.py` | cluster-vs-continuum model selection |
| `p2_3_pca12.py` | 12-D PCA, bootstrap CIs, Horn parallel analysis |
| `p2_4_jaccard_null.py` | layer-overlap null, BY-FDR |
| `p2_5_scc_null.py` | strongly-connected-component null |
| `p2_7_proximity_sensitivity.py` | per-layer density, both reciprocity conventions |
| `p2_8_phylogeny_pgls.py` | phylogenetic signal (h3) |
| `p2_9_spatial_decomposition.py` | *needs per-configuration matrices; skips otherwise* |
| `p2_10_directionality.py` | direction-convention audit |
| `p2_11_dominance.py` | David's score, Bradley–Terry, cross-layer dominance |
| `p2_12_metabolism.py` | carbohydrate assimilation vs network position |
| `p2_13_taxon_sensitivity.py` | n = 56 *Streptomyces*-only robustness check |
| `generate_figs.py` | main-text figures |

A stage whose input data is absent records a skip in `supplementary/report.md`
and the run continues; the summary at the end lists what ran and what did not.

## Notes

- **Reproducibility.** Seed `20260522`, set in `code/data_io.py`. All randomised
  procedures (null draws, permutations, bootstraps, GMM initialisation) derive
  from it.
- **Tree building.** RAxML has no pip distribution, so the tree is built with
  MAFFT + trimAl + VeryFastTree (a validated FastTree-2 equivalent).
  `data/phylogeny/tree_provenance.txt` records exactly what ran — those are the
  tool names that must appear in the manuscript Methods. If a tree from an
  original MAFFT + RAxML pipeline is placed at `data/phylogeny/tree.nwk`,
  `build_tree.py` leaves it alone and `p2_8_phylogeny_pgls.py` uses it unchanged.
- **Binary by design.** All per-layer analyses run on the binary projection
  `A_l = 1[M_l > 0]`. The layer constants `c_l` are severity tags, not per-edge
  weights, and enter only the aggregate visualisation matrix.
