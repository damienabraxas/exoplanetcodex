PREPARE solar Si  -> NOT READY -- blocking: nlte_missing 235
  F literature: litscan yes; set bibliography keys ['amarsi2017_si', 'bergemann2013_si_jband', 'deshmukh2022_si', 'elgueta2026']
  A published sets: 59 lines; adopt authors' gf on 0; missing from canonical_gf: 1
      MISSING SI_ELGUETA2026 Si I 10025.743: MISSING in canonical_gf; --apply ADDS it: 1 components from vald_solar_ir_9500_17000_hfson_raw.txt (centroid 10025.7433, VALD total -1.79) into data/linelists/ispec_ir_9200_13000/atomic_lines.tsv, canonical row at the authors' -1.79 +/- None
  B gf source PR2024_calc: 391 graded lines, classified
  B gf source NIST_ASD: 27 graded lines, classified
  B gf source PR2024_exp: 14 graded lines, classified
  B gf source AGSS21_Si_Amarsi2017: 10 graded lines, classified
  B gf source DenHartog2023: 2 graded lines, NOT in gf_error_model (UNREVIEWED)
  F2 literature need: 21 Si I lines in the near-UV band, none with a published gf uncertainty: no Codex/Deep pool until a lab or error-stated gf source is found (RYA-1237 register input) [solar_kpno_kurucz2005_corrected, solar_kpno_molecfit_corrected]
  F2 literature need: 3 Si II lines in the near-UV band, none with a published gf uncertainty: no Codex/Deep pool until a lab or error-stated gf source is found (RYA-1237 register input) [solar_kpno_kurucz2005_corrected, solar_kpno_molecfit_corrected]
  C cull candidates (saturated, not a solar-literature line): 0
  D NLTE: 235 graded line(s) not in data/nlte_grids/Si_Amarsi2020_PySME.csv -> on Sirius: scripts/extend_nlte_grid.py --element Si --write
  E holdings: 18 windows; 0 not serving / frame off rest
