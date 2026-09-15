# Multilayer interaction networks in a soil actinobacterial collection

Analysis code and data for the manuscript. Every number, table and figure the paper
reports is produced here from the raw inputs, which are included.

```bash
pip install -r requirements.txt
python run_all.py
```

Everything installs from PyPI; there are no external binaries. A cold run
regenerates the degree-preserving null models, which dominates the runtime
(about two hours on six cores). They are cached under `.cache/` afterwards and
are fully determined by the seed in `code/data_io.py`, so a second run is quick.

## Layout

```
code/          analysis and figure modules
data/          raw inputs: interaction matrices, 16S records, plate photographs
results/       generated tables, figures and report.md
run_all.py     runs every stage in order
```

## What each module does

| module | produces |
|---|---|
| `prep_phylogeny.py`, `build_tree.py` | the 16S alignment and maximum-likelihood tree |
| `layer_descriptors.py` | per-layer density, reciprocity, out-degree |
| `layer_overlap.py` | cross-layer Jaccard overlap against a degree-preserving null |
| `connectivity.py` | largest strongly connected component against the same null |
| `reciprocity_transitivity.py` | reciprocity and transitivity against the same null |
| `ordination.py` | PCA of the 12-layer profile, bootstrap intervals, parallel analysis |
| `mixture_models.py` | Gaussian mixtures, held-out likelihood, dip tests |
| `dominance.py` | David's scores and Bradley-Terry strengths per layer |
| `directionality.py` | per-layer directional asymmetry and the edge-direction convention |
| `phylogenetic_signal.py` | Pagel's lambda by tip-label permutation, and Mantel tests |
| `metabolic_niche.py` | carbohydrate assimilation against network position |
| `taxon_sensitivity.py` | every headline result without the four non-Streptomyces isolates |
| `spatial_configuration.py` | IG and CS with the inocula close versus separated, from the paired spatial matrices |
| `tables_1_2.py` | Tables 1 and 2 |
| `figure_*.py`, `figures_*.py` | the manuscript figures |
| `verify_reported_values.py` | recomputes every number the manuscript states and compares |

## Edge direction

The interaction matrices are stored with **rows as receivers and columns as
senders**: `M[i, j] = 1` means strain *j* elicited the phenotype in strain *i*.
This is the transpose of the i to j convention the manuscript uses in prose, and
`code/directionality.py` records it. Anything reading the spreadsheets directly
must transpose, or every directed result reverses.

## Data

`data/matrices_xls/` holds one spreadsheet per phenotype, 60 x 60, one row and
column per isolate. `data/matrices_spatial/` holds the four matrices of the paired
spatial experiment (IG and CS, each scored with the inocula close and separated);
they are a separate scoring and are never merged with the twelve layers.
`data/annotation/` holds the 16S records, the pairwise distance matrix and the API assimilation profiles. `data/plates/` holds the plate
photographs behind Figure 1 and the plate montage behind Figure S13. Strain labels differ slightly between sources; they
are folded to a single key space by `canonical_strain()` in `code/data_io.py`.
