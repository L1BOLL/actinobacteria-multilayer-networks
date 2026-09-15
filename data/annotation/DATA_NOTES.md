# Annotation data — integrity notes

Three annotation files join to the 60-strain interaction matrices:

| File | Contents | Join |
|---|---|---|
| `16S_strepto_R` | 61 FASTA records, near-full-length 16S | 60 / 60 after cleaning |
| `distance_matrix.xlsx` | 1770 = C(60,2) pairwise 16S distances, complete | 60 / 60 |
| `atb_profile.xlsx` | 60 × 9 carbohydrate assimilation, graded 0–4 | 60 / 60 |
| `list_accession_submission.ods` | 60 GenBank submission records, PZ882263–PZ882322 | 60 / 60 |

Strain labels are folded to a shared key space by `data_io.canonical_strain()`
(strips non-breaking spaces and apostrophes, then applies `STRAIN_ALIASES`).
`code/prep_phylogeny.py` writes a per-record log to
`data/phylogeny/label_map.csv`.

## Label alias: `S1430` = `S1226`

The sequencing record, the distance matrix and the GenBank submission (Seq01,
PZ882263) label one isolate `S1430`; the interaction matrices call the same
isolate `S1226`. The alias `"S1430": "S1226"` in `code/data_io.py ::
STRAIN_ALIASES` joins them, so every phylogeny-linked analysis runs at n = 60.
Supplementary Table S1 lists the strain as `S1226` and notes the second label.

## Four strains are not *Streptomyces*

Parsed from the FASTA descriptions, all four are in the interaction set:

| Strain | Genus |
|---|---|
| `MS53 7` | *Saccharothrix* |
| `MS10 8` | *Lentzea* |
| `MS3 15` | *Lentzea* |
| `MS3 18` | *Amycolatopsis* |

The collection is therefore described as actinobacterial, and
`code/taxon_sensitivity.py` re-runs the headline analyses without these four
(Supplementary Table S4).

## Handled by `prep_phylogeny.py`

### Malformed header (`MS8 6`)

The header line `>Streptomyces sp. strain MS8 6 ATGCGGTGCTACC...` has 77 bp of
sequence concatenated onto it. The sanitizer detects a trailing DNA run on any
header, splits it off, and prepends it to the record. Logged as
`recovered 77 bp from header`.

### Second record labelled `S705A`

The FASTA holds two records labelled `S705A`, of 1,390 bp and 1,468 bp; they are
about 92% identical, so they are not two reads of one sequence. Only the 1,468-bp
sequence was deposited (Seq40, PZ882302) and only it is used; the sanitizer keeps
the longer record and logs the discard.

## Resolution limit of the marker

The 60 strains collapse to **38 distinct 16S genotypes**. 33 strains fall into 11
groups whose pairwise 16S distance is exactly zero, including two groups of six
and three of three.

This is a property of 16S at within-genus resolution rather than a data defect,
but it caps the statistical power of any phylogenetic-signal test on these
strains, and it means a maximum-likelihood tree of the 60 isolates is largely
unresolved polytomy. It is reported alongside the phylogenetic-signal result
because it is the mechanistic explanation for it.
