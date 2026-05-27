§9 — Per-layer density and reciprocity

Reciprocity (pairs-based):
  RAC    R_p = 0.297
  IG     R_p = 0.150
  others R_p < 0.1

Mutual pairs: RAC 384, IG 100; 6 layers ≤10; IC, IAC_RDE, IRP at 0–1.
Pairs-based and edges-based give identical layer orderings.

Density gradient: RAC 0.474 → IRP 0.005.

SCC sizes (bimodal across MC, no clean per-category ranking):
  CMP 39, CS 26, CCVM 24 — but CCAM 4, IRP 1, IAC_RDE 1, RP 3, IC 2
v14 "20–40 in MC/MM" removed.

Scripts:
  p2_7_proximity_sensitivity.py   density, reciprocity (both), mutual pairs, largest SCC
  p2_10_directionality.py         per-phenotype asymmetry audit

Outputs:
  phase2_outputs/tables/p2_7_proximity_sensitivity.csv
  phase2_outputs/tables/p2_10_directionality_audit.csv
  phase2_outputs/figures/p2_7_proximity_density.png
