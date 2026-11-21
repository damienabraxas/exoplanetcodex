# Li et al. 2023 O I transition data (RYA-1354)

Li, W., Jönsson, P., Amarsi, A. M., Li, M. C. & Grumer, J. 2023, A&A 674, A54, "Extended atomic data for oxygen
abundance analyses", doi:10.1051/0004-6361/202245645 (arXiv:2304.02310, sha256 7387421…ebba0f).

| File | Content | Source |
|---|---|---|
| `raw/tablea1.tsv` | 81 states: MCDHF and NIST energies, lifetimes (B, C gauge) | CDS J/A+A/674/A54/tablea1, VizieR CfA mirror (vizier.cfa.harvard.edu), 2026-11-18; cdsarc.cds.unistra.fr timed out |
| `raw/tablea2.tsv` | 989 E1 transitions: S, log gf, A (B and C gauge), dT, CF, accuracy class | CDS J/A+A/674/A54/tablea2, same mirror |
| `table3_experimental_gf.csv` | 12 UV lines (98.9–130.6 nm), every experimental gf Li Table 3 lists | transcribed from the arXiv PDF, p. 6 |
| `table2_experimental_lifetimes.csv` | 13 levels, every experimental lifetime Li Table 2 lists | transcribed from the arXiv PDF, p. 5 |

The experimental values are **as quoted by Li 2023**. The primary papers have not been re-read. Transcription checks:
- Table 3 against Table A2: `tests/test_rya1354_theory_calibration_oi.py`.
- Table 2 against Table A1 lifetimes: `scripts/rya1354_theory_calibration_oi.py`, which raises on any mismatch.

The NIST O I lines are in `data/reference/nist/oi_lines_asd_rya1354.csv` (sha256 f4e652eb…). They are NIST ASD, 50–6000 nm vacuum, retrieved 2026-11-18. NIST's O I transition probabilities come from Wiese, Fuhr & Deters 1996 and are themselves calculated, mostly Hibbert et al. 1991 (TP ref T5803).
