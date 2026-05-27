§8 — Layer overlap structured by shared mechanism

mean J across 66 pairs = 0.047

Confound (excluded from biology summary):
  CCAM/CCVM   J = 0.336, z = +16.2, p_BY = 0.005

8 over (p_BY < 0.05) of 65 biological pairs:
  IRAC/RAC   +9.2
  CCVM/CMP   +8.2
  CCAM/CMP   +7.7
  CCVM/CS    +6.4
  CMP/CS     +3.7
  CCAM/CS    +3.9
  CS/RAC     +3.7
  IG/IRAC    +4.8

1 under: IAC_RAC/RAC z = −5.3

SCC null: no layer significant after BY-FDR (|z| < 1.5 every layer).
SCC size driven by density, not feedback.

(Reverses v14 "most pairs significantly lower than chance".)

Scripts:
  p2_4_jaccard_null.py   N_DRAWS=10000, CONFOUND_PAIRS={frozenset({"CCAM","CCVM"})}, BY+BH
  p2_5_scc_null.py       N_DRAWS=10000

Outputs:
  phase2_outputs/tables/p2_4_jaccard_null.csv  (J_obs, J_null_mean, J_null_sd, z, p_emp, p_bh, p_by, is_confound)
  phase2_outputs/tables/p2_5_scc_null.csv
  phase2_outputs/figures/p2_4_jaccard_heatmap.png
  phase2_outputs/figures/p2_5_scc_null.png
