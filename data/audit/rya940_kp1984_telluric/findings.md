# RYA-940 — telluric-correcting the 1984 Kitt Peak solar flux atlas

## What made this different from HARPS

The atlas carries **no observation metadata at all** — no airmass, no date, no site.
Verified by grepping all 251 segments: the only non-numeric text in the whole
distribution was a saved HTTP 500 page. So the ratified observation-night-GDAS rule
is unsatisfiable, and Ryan ratified the alternative on 2026-08-21: fit the molecular
columns freely against a standard atmosphere, and referee against a telluric-free
reference instead of against provenance.

That costs the referee RYA-931 relied on. For HARPS, O2 is well mixed, so a correct
fit HAD to return a column of 1.0 — and did, ten times independently. Here the
airmass is unknown and is absorbed into the column, so that test is gone.

**What replaces it:** the airmass is unknown, but it must be the SAME unknown in every
O2 band. O2 is well mixed, so `rel_mol_col` is an airmass proxy for O2 specifically.
Three O2 bands, fitted independently and 1300 A apart, agree:

| band | implied airmass |
|---|---:|
| o2a 7594–7685 A | **1.352** |
| o2b 6867–6884 A | **1.303** |
| o2gamma 6270–6300 A | **1.306** |

Agreement to 4%. That is not something a wrong correction produces.

⚠️ For **H2O this is NOT an airmass** — the column scales airmass × (PWV / reference
PWV), and water vapour is not well mixed. The H2O bands return 1.9–2.6 against O2's
1.30–1.35; that is water vapour, not a contradiction. An earlier version of this
harness labelled the H2O column `implied_airmass` and would have published an
agreement that does not exist.

## Results

| band | range (Å) | col scale | fit rms | T_min | median before → after | %<0.5 before → after | RMS vs Kurucz 2005 |
|---|---|---:|---:|---:|---|---|---|
| h2o11120 | 11120–11560 | 2.163 | 0.0468 | 0.2158 | 0.6988 → **1.0038** | 33.86 → **0.09** | **no reference exists** |
| h2o7160 | 7160–7340 | — | — | — | **NO ADMISSIBLE FIT** | — | — |
| h2o8100 | 8100–8400 | 2.591 | 0.0269 | 0.1422 | 0.9725 → **0.9936** | 3.58 → **0.52** | 0.1546 → **0.0699** |
| h2o9280 | 9280–9600 | 2.490 | 0.0623 | 0.4255 | 0.7498 → **0.9976** | 29.12 → **0.01** | 0.4730 → **0.0604** |
| o2a | 7594–7685 | 1.352 | 0.0497 | 0.2561 | 0.7440 → **0.9908** | 32.83 → **0.55** | 0.5003 → **0.0915** |
| o2b | 6867–6884 | 1.303 | 0.0355 | 0.2213 | 0.9132 → **0.9904** | 15.80 → **1.01** | 0.3593 → **0.1210** |
| o2gamma | 6270–6300 | 1.306 | 0.0237 | 0.0873 | 0.9854 → **0.9991** | 1.12 → **0.26** | 0.0802 → **0.0495** |

**11120–11560 Å has no external referee** — Kurucz 2005 stops at 10000 Å and IAG at
10000 Å. It is corrected and marked `externally_validated: false`. Its only check is
that its H2O column scale (2.163) sits inside the family the other H2O bands return
(1.95–2.59). That is weaker than an external comparison and is not presented as more.

## The null result

**H2O 7160–7340 Å is NOT corrected.** All six starting points — spanning 2.5 to 8.1
pixels, bracketing the 5.05 the measured resolving power implies — converged to the
*identical* solution with a **zero-width LSF**. Identical from every start means this
is not a basin problem: the band's absorption is too weak to constrain the kernel at
all, and the optimiser drives it to zero because that makes the model lines sharpest.
The fit rms (0.0209) is among the best of any band, which is exactly why fit quality
alone must never be the gate. Those pixels stay under the existing per-line
clean-line selection.

## Four defects, each of which looked like a fit problem and was not

1. **A continuum was being fitted that the atlas already ships.** Column 2 is a
   *pseudo-residual flux*: unity IS the continuum. A free degree-3 polynomial across a
   window one-third saturated followed the band down (median 0.832 against an atlas
   median of 0.93), so the band was explained by continuum rather than absorption and
   the O2 column collapsed to **0.24**. This is the RYA-911 double-continuum defect in
   a new place.

2. **Solar lines must be masked out of an atmospheric fit.** Unmasked, all six starts
   railed the LSF at its 100-pixel bound — the fitter smoothing its model into mush
   trying to cover solar structure.

3. **The wavelength fit ran away.** `chip 1, coef 0` pinned at exactly −0.050000 — its
   bound, quoted to ±1e-6 — imposing a spurious **−3.3 Å** displacement. A wide
   cross-correlation found it as a sharp peak at −3.316 Å (correlation 0.735 against
   0.40 at zero lag). An earlier ±0.1 Å scan had missed it entirely and read as 'no
   shift'. An FTS atlas's wavelength scale is its single greatest strength and RYA-938
   verified this one against IAG at −0.000 Å, so it should never have been free.
   Fixing it: reduced χ² **1892 → 24.7**, residual **43% → 5.0%**.

4. **The LSF start ladder was in PIXELS.** Resolving power is an instrument property;
   FWHM in pixels is not, because the atlas sampling runs 0.0067–0.0183 Å across
   2960–13000 Å. A fixed pixel ladder is wrong everywhere except where it was tuned.
   The ladder is now derived from the MEASURED R = 148,000 (from the O2 A fit:
   6.040 px × 0.00850 Å at 7640 Å).

## On the acceptance criterion

The χ² gate was replaced. This atlas ships no uncertainties, so reduced χ² is computed
against molecfit's *assumed* `DEFAULT_ERROR` — gating on it gates on our own
assumption, and the same fit passes or fails purely by changing that constant. The
criteria are now error-independent and physical, declared before any fit ran:
a line-spread function cannot have zero width; the model must describe the data to
better than 10% fractional residual; and the column must correspond to a possible
airmass. χ² is recorded, never gated on.

## Caveat on the airmass agreement

The 1984 atlas is a **composite of FTS scans** taken at different times. Different
segments may genuinely carry different airmass and water vapour, so the three O2 bands
agreeing is evidence that those particular segments came from similar conditions — it
is not proof that a single atmosphere describes the whole atlas. Correction is done
per band for exactly that reason.

## Products

- `kp1984_corrected_11120_11560.txt` — sha256 `28508c3cf8decf46…`
- `kp1984_corrected_6270_6300.txt` — sha256 `8aa8ed41f23536cb…`
- `kp1984_corrected_6867_6884.txt` — sha256 `e66a1ffbbea7a6a0…`
- `kp1984_corrected_7594_7685.txt` — sha256 `b2e8cb03995ba472…`
- `kp1984_corrected_8100_8400.txt` — sha256 `13d77031a222dbed…`
- `kp1984_corrected_9280_9600.txt` — sha256 `7a3a37f52d56b7d9…`

Columns: air wavelength (Å), corrected pseudo-residual flux, transmission.
**NaN = quarantined** (transmission below the per-band derived `T_min`), never divided.
Products keep the atlas's own uneven sampling; the even grid was a fitting device only.
