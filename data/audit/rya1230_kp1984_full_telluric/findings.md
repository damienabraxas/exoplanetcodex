# RYA-1230 -- full-coverage telluric correction of the 1984 Kitt Peak atlas

Ryan, 2026-10-01: "every band gets telluric except in the case of a space telescope."

## What was wrong
`solar_kpno_molecfit_corrected` was RAW (byte-identical to `solar_kpno`) over 8,200 of
9,800 A. RYA-940 molecfit-corrected only six registered windows, and the real bands are
wider (the 1.13 um H2O band spans ~10900-11900 A; RYA-940 fitted 11120-11560). Every
KP-molecfit product outside those windows -- the CN A-X windows 10990-11110 / 11576-11800,
C I 9061-9078 in the uncorrected 9000-9280 H2O band -- was fitted on uncorrected flux.
RYA-1191's evidence table only ever audited the six windows, so the gap was invisible.

## What was done
`scripts/rya940_kp1984_correct.py` (ESO molecfit, local Homebrew esorex 3.13.10) per
molecular complex, Kurucz 2005 `irradthu.dat.txt` as solar mask + referee below 10000 A:

| complex (A) | molecules | LSF px | columns | vs Kurucz rms |
|---|---|---|---|---|
| 3000-5800 (from joint 3000-6860) | H2O, O2 | 5.31 | H2O 2.491, O2 1.477 | 0.098 -> 0.098 |
| 5800-6860 | H2O, O2 | 4.28 | H2O 2.435, O2 1.362 | 0.023 -> 0.016 |
| 6860-7550 | H2O, O2 | 3.42 | H2O 1.989, O2 1.300 | 0.109 -> 0.048 |
| 7550-8000 | H2O, O2 | 6.09 | H2O 2.831, O2 1.356 | 0.231 -> 0.068 |
| 8000-8900 | H2O | 5.47 | H2O 2.609 | 0.095 -> 0.047 |
| 8900-10000 | H2O | 5.14 | H2O 2.525 | 0.303 -> 0.076 |
| 10000-11860 | H2O | 3.92 | H2O 2.181 | no referee |
| 11860-13000 | H2O, O2, CO2 | 2.58 | H2O 2.27, O2 1.374, CO2 1.863 | no referee |

- O2 columns (airmass proxy, well mixed) agree 1.30-1.37 across the red complexes and match
  RYA-940's 1.30-1.35. The joint blue fit gives 1.48, which is why 5800-6860 is served by its
  dedicated fit and the joint fit only for [3000, 5800).
- Blue transmission (measured, joint fit): 3000-4200 A min T >= 0.992 (no line absorption
  to speak of); 4200-5000 min 0.987; 5000-5800 min 0.927. The blue band alone cannot be fitted
  (LSF collapses, columns run away: O2 1.91, H2O 5.5/100) -- the red lines constrain it.
- Narrow/weak windows collapse the LSF (10850-11120 alone, 10000-10850 alone): fit with a
  strong neighbour instead. Script gained `--fix-lsf-ratio` and `RYA1230_FIXED_COLS`
  (unused in the installed product; recorded for the blue diagnosis).
- Each file is trimmed to its own band [lo, hi) so no pixel is served twice. The six RYA-940
  window files are in `superseded_rya940/` (not globbed by the loader).

## Measured coverage audit
`solar_kpno` vs `solar_kpno_molecfit_corrected`, 50 A chunks, 3000-13000 A:
**10000 A corrected, 0 A raw, 0 unreadable, 0 duplicate seams.**
`telluric_policy.VERIFIED_HOLDING_STATE` for the holding flipped raw -> corrected on this.

Fit manifests: `c_*/fit_manifest.json` beside this file.
