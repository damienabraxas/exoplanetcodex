# RYA-1232 — HARPS direct-Sun full-range telluric correction

`solar_harps_molecfit_corrected` was corrected for **O2 B only** (RYA-931). On sky pixels
it was byte-identical to raw at H2O 5700-6550 A and O2 gamma 6270-6300 A.

**Method** (`scripts/rya1232_harps_full_telluric.py`): RYA-931's O2 B molecfit_model fit
(unchanged ladder and gates, H2O held at the night's GDAS column) -> molecfit_calctrans over
the full exposure with the identity wavelength correction -> H2O column scale fitted
against the IAG telluric-free atlas over three H2O windows. Inside the O2 B fit window the
fitted model transmission is kept (it matches RYA-931 exactly).

Failed first attempt, recorded: one molecfit_model fit over all windows railed the LSF to
100 px on every start (solar residuals swamp the weak H2O lines). The O2-B chip's fitted
wavelength polynomial, carried to the full chip, shifted the model ~0.95 A (H2O scale 0.005)
-- hence the identity correction.

**Result, ten exposures** — mean (obs / atlas - 1) on sky pixels, raw -> RYA-931 -> RYA-1232:

| window | raw | RYA-931 | RYA-1232 |
|---|---|---|---|
| H2O 5695-5790 | -0.034 | -0.034 | -0.015 to -0.016 |
| H2O 5900-5985 | -0.070 to -0.072 | same | -0.012 |
| O2 gamma 6270-6300 | -0.145 to -0.147 | same | -0.007 to -0.009 |
| H2O 6465-6535 | -0.073 to -0.076 | same | -0.007 |
| O2 B 6867-6884 | -0.25 to -0.26 | -0.014 | -0.014 |

O2 column 1.002-1.003, reduced chi2 0.77-0.80, LSF 6.76-6.77 px (all match RYA-931);
H2O scale 1.003-1.037 x GDAS. Exposure .745 needed RYA-931's O2-only fit with the H2O
profile shape taken from .735 (same night, same GDAS file); its amount is still fitted.
Installed holding vs raw on KP sky pixels: +1-4% (H2O), +6-9% (O2 gamma), +29% (O2 B);
the telluric-free 6150 A window is unchanged.
