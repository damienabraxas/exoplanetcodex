# RYA-931 — the Molecfit runtime blocker, and what was behind it

## 1. The blocker was never MIPAS

RYA-927 stopped with ESO Molecfit 4.4.4 reporting, immediately before the fit:

```
mf_config_chk_ref_atmos: Unable to open rat data file
    /opt/homebrew/Cellar/telluriccorr/4.3.3_4/share/molecfit/data/profiles/mipas/equ.fits!
molecfit_model: Bad Specifications of Molecules and Reference Atmosphere
```

That message names the last file the recipe touched, not the cause.

**Measured root cause.** The HARPS Phase-3 `ERR` column is entirely NaN — **0 of
313,028 pixels finite, in all ten direct-solar exposures**. Selecting pixels on
finite errors left a **zero-row** SCIENCE table. Molecfit then failed inside
`molecfit_data_extract_spectrum_from_cpl_table`:

```
[1/3] 'Access beyond boundaries' (11) at cpl_column_get_double
[2/3] 'Access beyond boundaries' (11) at cpl_column_get_double
[3/3] 'Input data do not match: Parse 1D spectrum from BINTABLE failed!' (13)
```

and only *afterwards* printed the reference-atmosphere line off an already-dirty
CPL error state.

**Positive control.** "MIPAS is fine" only means something if the same file
demonstrably works. The passing run loads the **identical** `equ.fits` — same
path, same md5 `cf98a77b2e84…` — and reports no failure. The only changed input
is the number of SCIENCE rows: 0 → 7,201. Recorded in `runtime_rca.json`.

An unusable *error* column must never be allowed to empty the *flux* table. The
runner now refuses to write an empty SCIENCE table and names the flux when the
band is genuinely empty.

## 2. Three further defects the first fix exposed

Each was found by measurement, and each would have produced a plausible-looking
but wrong product.

**a. esorex splits SOF lines on whitespace.** Absolute paths under the store root
`~/Documents/Exoplanet Codex/` were torn in two (`Could not open the input file
'…/Exoplanet' with tag 'Codex/…'`, rc=9). Frames are now named relatively and
the recipe runs with `cwd` set to the input directory.

**b. The default line-spread function ran away.** HARPS is fibre-fed with no
slit, but molecfit's boxcar and Lorentzian components are free by default with
nothing to constrain them. The first run drove the Lorentzian FWHM to its
100-pixel bound and paid for the over-smoothed model by collapsing the O₂ column
to **3.6% of atmospheric** — a converged fit, with products, that was physically
impossible. Disabling both and fitting a single Gaussian gives reduced χ²
**881.9 → 0.806** and the O₂ column **0.036 → 1.0009**.

`FIT_RES_GAUSS=FALSE` does **not** hold the kernel at a given width — it disables
convolution. The same 6.767 px gives reduced χ² 0.81 fitted and 894 fixed. The
LSF must stay free.

**c. `BEST_FIT_MODEL.lambda` is in vacuum.** The column is the declared air grid
times the refractive index (ratio 1.000277 here). Mapping the transmission back
by wavelength mis-registers the correction by ~1.9 Å — about 190 HARPS pixels.
The output table is row-for-row the SCIENCE table, which the flux column proves
exactly, so the row mapping is used instead and asserted before use.

Also: `mtrans` is **0** outside the fitted region, not 1. Dividing by it
wholesale would have zeroed everything blueward of 6857 Å.

## 3. The fit is bistable — refereed by physics, not by a target

Alongside the physical solution sits a degenerate one that converges just as
cleanly with a **zero-width** LSF, reduced χ² ≈ 11.5, and the O₂ column 3.7% off.
Which basin an exposure lands in is chaotic in the starting LSF width: **6.767 px
and 6.773 px send the same exposure to different answers.**

Acceptance is therefore stated in advance and physical, not tuned:

| criterion | basis |
|---|---|
| Gaussian FWHM ≥ 1 px | a line-spread function cannot have zero width |
| reduced χ² ≤ 2.0 | the model must actually describe the data |
| \|O₂ column − 1\| ≤ 0.10 | O₂ is well mixed; its column is known far better than 10% |

Only the **starting point** is varied on retry. No fit parameter, threshold or
result is adjusted. All ten exposures reach an admissible solution within five
starts.

**Corroboration, not tuning.** The fitted O₂ column lands at **1.0003–1.0018**
across ten exposures, and the fitted LSF agrees to **6.7683 ± 0.0059 px**
(R = 101,710) — an instrument property recovered independently ten times.

## 4. Lineage is proven, not assumed

Re-running `pipeline.spectra_normalize` — the same functions that built the held
product — on the ten originals reproduces the committed `solar_normalized.csv`
**byte for byte** (SHA-256 `6100681d530ffe59…`). The corrected holding therefore
differs from `solar_harps` by the correction and by nothing else.

## 5. Results

O₂ B 6867–6884 Å, normalised flux:

| metric | before | after |
|---|---:|---:|
| minimum | 0.0060 | **0.7725** |
| median | 0.8505 | 0.9648 |
| % below 0.5 | 22.06 | **0.00** |
| % below 0.95 | 86.35 | 24.39 |
| RMS vs IAG telluric-free atlas | 0.4021 | **0.0351** |

The IAG comparison is an independent discriminator: a different telescope and a
different correction method, so agreement is not self-confirmation.

Solar structure survives: max clean-line depth shift **0.0052** (tolerance 0.01),
and the raw co-add outside the fitted region is **bit-identical**.

**256 saturated-core pixels are quarantined as NaN, not divided.** The threshold
is derived per exposure — `T_min = measured δT / 0.05`, i.e. a corrected pixel may
carry at most 5% relative flux error from the transmission model — not chosen.

## 6. The RYA-904 guard caught the holding, and it was right to

Registering the corrected holding as `verified` made CI fail:
`loader_coverage.reconcile_loader_coverage` reported it VERIFIED, passing the
telluric gate, and **UNADDRESSABLE by the band harness** — data we hold, in a
state we may measure, that no loader can name. Unreachability leaves no trace;
it reads exactly like having no data.

It is wired rather than exempted. `measure_band_ew.harps_spectrum` was keyed to
one module-level path with a single-slot cache — the RYA-904 defect shape, one
instrument silently standing for one holding — so it is now keyed by **product**,
and HARPS carries two `HoldingSpec`s sharing one reader.

The corrected holding is listed **second on purpose**. It is reachable by name,
but selection order is unchanged, so no existing measurement silently switches
product. Choosing it for the affected red-edge windows is **RYA-936's decision to
make and record**, not a side effect of this wiring.

## 7. Honest caveats

- Removing the O₂ band raises the 95th-percentile continuum anchor of the knot
  window containing it, so the continuum solution shifts by up to **0.389%**
  outside the fitted region (normalised flux by ≤1.6e-3). This is a real and
  wanted consequence — the old anchor was depressed by telluric absorption — and
  it is measured and reported rather than asserted away.
- `telluric_applied=applied` is a **product-level conditioning fact, not a
  line-level science verdict**. Every candidate line in 6857–6911 Å still needs a
  local telluric disposition under the RYA-927 ratified two-level policy.
- The `ERR` column is unusable, so molecfit was weighted with its documented
  `DEFAULT_ERROR`. That is fit weighting only and **is not an EW uncertainty**.
- `pipeline/spectra_normalize.fit_continuum` now seeds its mask from the finite
  pixels. On an all-finite spectrum this is `np.ones(...)` exactly, and the
  byte-identical reproduction control proves the frozen products are unchanged.
