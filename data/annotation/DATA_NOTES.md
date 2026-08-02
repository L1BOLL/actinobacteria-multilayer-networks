# Annotation data — integrity notes

Three annotation files join to the 60-strain interaction matrices:

| File | Contents | Join status |
|---|---|---|
| `16S_strepto_R` | 61 FASTA records, near-full-length 16S | 59 / 60 after cleaning |
| `distance_matrix.xlsx` | 1770 = C(60,2) pairwise 16S distances, complete | 59 / 60 |
| `atb_profile.xlsx` | 60 × 9 carbohydrate assimilation, graded 0–4 | **60 / 60** |

Strain labels are folded to a shared key space by `data_io.canonical_strain()`
(strips non-breaking spaces and apostrophes, then applies `STRAIN_ALIASES`).
`phase2_code/prep_phylogeny.py` writes a per-record audit to
`data/phylogeny/label_map.csv`.

---

## OPEN — needs an author decision

### 1. `S1226` vs `S1430`

`S1226` is present in the interaction matrices but in neither 16S source.
`S1430` is present in both 16S sources but in no interaction matrix.
Every other one of the 59 remaining strains matches exactly after normalization.

This looks like one isolate renamed between the bench work and the sequencing,
but it has **not been assumed**. All phylogeny-linked analyses run at **n = 59**
with both labels excluded. The carbohydrate analysis is unaffected and runs at
n = 60.

**If the authors confirm S1226 = S1430**, add `"S1430": "S1226"` to
`STRAIN_ALIASES` in `phase2_code/data_io.py` and re-run; every downstream script
picks it up and n becomes 60. Nothing else changes.

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
