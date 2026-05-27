§10 — Strain strategy

Total out-degree spans ~1 order of magnitude across 60 strains.
High-output: MS53_7, S1350, S1705, MS10_14, S625, B201.

P_i ∈ [0.60, 0.84] for every strain. No narrow-spectrum specialists.
v14 generalist/specialist dichotomy gone.
Meaningful axis = total output, not breadth.

Scripts:
  p2_7_proximity_sensitivity.py   per-strain out-degree, 12-layer P_i
  generate_figs.py                Fig. 3 panels

Outputs:
  phase2_outputs/tables/p2_7_proximity_sensitivity.csv
  figure_outputs/figure3_layer_similarity_partitioning.{png,svg}

Caveat: Fig. 5 participation encoder ramp is wide; rescale to [0.60, 0.84]
or swap for summed out-degree.
