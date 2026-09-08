# Annotation data — integrity notes

Three annotation files join to the 60-strain interaction matrices:

| File | Contents | Join status |
|---|---|---|
| `16S_strepto_R` | 61 FASTA records, near-full-length 16S | 59 / 60 after cleaning |
| `distance_matrix.xlsx` | 1770 = C(60,2) pairwise 16S distances, complete | 59 / 60 |
| `atb_profile.xlsx` | 60 × 9 carbohydrate assimilation, graded 0–4 | **60 / 60** |
| `list_accession_submission.ods` | 60 GenBank submission records, PZ882263–PZ882322 | 59 / 60 |

Strain labels are folded to a shared key space by `data_io.canonical_strain()`
(strips non-breaking spaces and apostrophes, then applies `STRAIN_ALIASES`).
`phase2_code/prep_phylogeny.py` writes a per-record audit to
`data/phylogeny/label_map.csv`.

---

## OPEN — needs an author decision

### 1. `S1226` vs `S1430` — RESOLVED, applied 2026-09-07

**The two labels are one isolate.** `S1430` is the sequencing-record label for the
strain the interaction matrices call `S1226`; the alias `"S1430": "S1226"` is in
`code/data_io.py :: STRAIN_ALIASES` and every phylogeny-linked analysis now runs at
**n = 60**.

Three independent sources of 60 records each agreed, differing only on this label:
the 16S FASTA, `distance_matrix.xlsx`, and the GenBank submission record
`list_accession_submission.ods` (Seq01, accession PZ882263). Supplementary Table S1
presents them as one strain, listing Seq01 under the main-text label `S1226` and
noting in its caption that the submission record writes `S1430`. The authors
confirmed the rename.

What changed when the alias was applied:

| | n = 59 (before) | n = 60 (now) |
|---|---|---|
| tree tips | 59 | 60 |
| Mantel pairs | 1,711 | 1,770 |
| Mantel rho / P | −0.022 / 0.78 | −0.041 / 0.59 |
| minimum q over 63 tests | 0.47 | 0.39 |
| genotype collapse | 59 → 38, 31 shared (53%) | 60 → 38, 33 shared (55%) |
| ML VCV rank | 44 of 59 | 44 of 60 |

Conclusions did not change: still 0 of 21 traits significant on the ML tree and
0 of 42 on UPGMA/WPGMA. `prep_phylogeny.py` now reports no unmatched label on
either side. To reverse, remove the alias and re-run `run_all.py`.

### 2. Four strains are not *Streptomyces*

Parsed from the FASTA descriptions, all four are in the interaction set:

| Strain | Genus |
|---|---|
| `MS53 7` | *Saccharothrix* |
| `MS10 8` | *Lentzea* |
| `MS3 15` | *Lentzea* |
| `MS3 18` | *Amycolatopsis* |

The manuscript title and abstract say "soil *Streptomyces* communities" and
"60 soil-derived *Streptomyces* strains". Both need correcting to
"60 soil actinobacteria (56 *Streptomyces*, plus *Saccharothrix*, *Lentzea* ×2,
*Amycolatopsis*)". `phase2_code/p2_13_taxon_sensitivity.py` re-runs the headline
analyses without these four so the paper can state that no conclusion depends on
them.

Anyone who BLASTs the deposited sequences will find this, so it must be fixed
before submission rather than after review.

---

## RESOLVED — handled automatically by `prep_phylogeny.py`

### 3. Malformed header (`MS8 6`)

The header line `>Streptomyces sp. strain MS8 6 ATGCGGTGCTACC...` has 77 bp of
sequence concatenated onto it. A naive parser drops those bases. The sanitizer
detects a trailing DNA run on any header, splits it off, and prepends it to the
record. Logged as `recovered 77 bp from header`.

### 4. Duplicate record (`S705A`)

Two records share the label `S705A`, at 1390 bp and 1468 bp. The sanitizer keeps
the longer and logs the discard. If these are genuinely two different isolates
rather than a duplicated submission, the authors should relabel one of them —
as delivered they are indistinguishable.

---

## Resolution limit of the marker

The 60 strains collapse to **38 distinct 16S genotypes**. 33 strains fall into 11
groups whose pairwise 16S distance is exactly zero, including two groups of six
and three of three.

This is not a data defect — it is a property of 16S at within-genus resolution —
but it caps the statistical power of any phylogenetic-signal test on these
strains, and it means a "maximum-likelihood tree of 60 isolates" is largely
unresolved polytomy. It is reported alongside the h3 result rather than buried,
because it is the mechanistic explanation for that result.
