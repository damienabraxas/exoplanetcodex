# Civis et al. 2013, A&A 554, A24 -- Ca I (Table 5) and Mg I (Table 6) f-values

Tables 5 and 6 are CDS-only; fetched from the VizieR mirror at vizier.cfa.harvard.edu
(asu-tsv, J/A+A/554/A24/table5, table6) on 2026-11-07. `lambda` is in nm; log gf =
log10((2 J_low + 1) f_ik). f-values are quantum defect theory (LS coupling), not measured.

- `feature_gf_vs_canonical.csv`: Civis vs our canonical gf, summed over each multiplet's
  components on both sides. Mg I vs our LAB / NIST / published-set rows: median -0.025 dex,
  MAD 0.04 (28 features). Ca I H 5p-5d: Civis 0.15-0.18 below Kurucz (16157: -0.41). Ca I K:
  within ~0.1.
- `adoption.csv`: what `scripts/rya1311_adopt_civis2013_ca.py` did to canonical_gf (33 H/K
  Ca I rows adopted, 7 ambiguous multiplets left on Kurucz).
- `CaI_*_SYNTH_1D-LTE_lines.csv`: diagnostic CRIRES+ H and K Ca I runs on the Kurucz gf
  before adoption (scratch, never published). H median 5.98 on Kurucz, 6.23 rescaled to Civis;
  K 6.29. Ca II: no H/K line reaches 0.15 depth, so none is measurable.
