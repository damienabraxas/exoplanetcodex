# RYA-424 — Telluric correction as a standing data-input stage

Promotes the per-dataset one-offs (RYA-373 molecfit driver on Vesta CRIRES+,
RYA-380 per-night-GDAS recipe) into a **standing, instrument-aware stage** that every
red-optical / IR dataset passes through on data input, with a hard verification gate
before a spectrum can be marked analysis-ready.

Code: `pipeline/telluric_stage.py` (+ registries in `config/constants.py`).
Single-source molecfit correction is delegated to the RYA-373 driver
(`pipeline/crires_telluric.py`) — not re-implemented.

This stage is **upstream** of the solar full-spectrum run's Phase B (the 2.3 µm
12C/13C CO measurement), of Procyon (RYA-425), and of alpha Cen (RYA-423/432) — it must
be production-ready to *finish* the solar run, not parallel to it.

## The rule is wavelength-gated, not instrument-gated (RYA-380, codified here)

| span | regime | telluric correction |
|------|--------|---------------------|
| λ ≳ 6800 Å (IR) | `ir` | **MANDATORY** (sharp H2O / O2-A / CH4 / CO2 forest) |
| crosses 6800 Å (e.g. 6910–9500: O I 7772, K I 7699, N I 8216) | `red_optical` | **MANDATORY** |
| 3800–6800 Å | `mid_optical` | not auto-required (largely clean; case-by-case) |
| λ ≲ 3800 Å ground-based | `blue` | none (broad ozone/Rayleigh handled by normalization) |
| space-UV (HST/STIS/COS) | `space_uv` | none (above the atmosphere) |

Thresholds: `TELLURIC_LAMBDA_MIN_A = 6800`, `TELLURIC_BLUE_MAX_A = 3800`.
Functions: `telluric_regime(lo_A, hi_A)`, `requires_telluric(lo_A, hi_A, space_uv=)`.

## Engine selection — single source of truth, no silent default

`constants.TELLURIC_ENGINES` keys the engine by **instrument**. An unknown instrument
**loud-fails** (`EngineNotSelectedError`) — the engine is never guessed.

| instrument | engine | site (per-night GDAS) |
|------------|--------|------------------------|
| CRIRES, CRIRES+ | molecfit | Paranal |
| UVES-red, ESPRESSO-red | molecfit | Paranal |
| NIRPS | molecfit | La Silla |
| FEROS-red | molecfit | La Silla |
| CHIRON | molecfit | CTIO |
| **SPIRou** | **APERO + Wapiti** (NOT molecfit — permanent rule) | CFHT (Maunakea) |

`select_engine(instrument)`, `site_for_instrument(instrument)`.

## Per-night GDAS — loud-fail on silent fallback (the whole point)

molecfit's `GDAS_PROFILE=auto` silently falls back to a generic standard atmosphere
when the real 3-hourly observation-night profile is missing (the RYA-373 failure
mode). This stage:

1. `resolve_gdas_profile(mjd, instrument, work_dir)` — resolves the real nearest-3-hourly
   GDAS profile for the obs MJD at the instrument's site (per-site tarball), as a
   molecfit-ready FITS. **Raises `GDASUnavailableError`** if no real profile within ±6 h,
   if the site tarball is absent, or on a NaN MJD. No standard-atmosphere fallback.
2. `assert_real_gdas(provenance, instrument)` — post-run, the GDAS provenance recorded by
   the engine must name a real profile, not the standard-atmosphere sentinel. Catches a
   silent fall-back that happened *after* resolution.

## Verification metric + tolerance

`telluric_residual_metric(...)` (generalizes the RYA-373 D1 gate): at pixels that are
**telluric-dominated** in the engine's model (`0.02 < transmission < 0.90` — absorbing
but not saturated) **and not science-coincident** (a `science_mask` excludes the wanted
stellar lines so a line-rich target is not scored as telluric misfit), the corrected
flux must return to the local continuum. The metric is the **median |1 − corrected/continuum|**
there; `passed` requires ≥10 clean telluric pixels AND residual ≤ tolerance.

- Tolerance: `TELLURIC_RESIDUAL_TOL = 0.05` (5%), inherited from the RYA-373 D1
  telluric-specific gate. A frame that does not pass is **flagged, not silently passed**.

## analysis_ready / telluric_verified — the conditioning manifest

`TelluricManifest` (written as `<dataset>_<wlen>.telluric.json`) carries the
analysis-ready contract: `telluric_required`, `telluric_verified`, `residual`,
`tolerance`, `engine`, `site`, `gdas_profile`, `regime`, `wave_lo/hi_A`, `mjd`, plus a
derived `analysis_ready` (= not-required OR verified).

`require_telluric_verified(manifest, [wave_lo_A, wave_hi_A])` is the **abundance-path
guard**: a red-optical/IR frame that is not `telluric_verified` raises
`TelluricNotVerifiedError`. Accepts a `TelluricManifest`, a dict, or a path; if a span
is supplied it re-runs the wavelength gate as a cross-check against a stale flag.

## Vesta CRIRES+ end-to-end (the acceptance smoke)

`python -m pipeline.telluric_stage --vesta-crires` — both on-chip 12CO(2-0) frames,
extracted → routed → real per-night GDAS → molecfit telluric (topocentric) → residual
verified → manifest:

| frame | GDAS profile (per-night) | residual | tol | telluric_verified |
|-------|--------------------------|----------|-----|-------------------|
| K2192 (MJD 59906.01, 2022-11-23) | `C-70.4-24.6D2022-11-23T00.gdas.fits` | 0.061 | 0.05 | **false (flagged)** |
| K2217 (MJD 59908.01, 2022-11-25) | `C-70.4-24.6D2022-11-25T00.gdas.fits` | 0.065 | 0.05 | **false (flagged)** |

The GDAS provenance matches the observation nights (real per-night profiles, not
standard atmosphere). The residuals (~6%) sit just above the 5% tolerance, so the gate
correctly **flags both frames as NOT analysis-ready** rather than silently passing —
consistent with the RYA-373 finding that the gate is PROVISIONAL pending the FTS solar
IR atlas (RYA-162, absent) and the RYA-387 0.001 re-extraction. The stage does what it
must: it does not let a not-yet-clean IR frame into the abundance path.

## Tool availability / install notes (this box, 2026-06-23)

| tool | role | status |
|------|------|--------|
| **molecfit** (ESO, via `esorex`) | telluric engine (CRIRES+/UVES-red/ESPRESSO-red/FEROS-red/CHIRON/NIRPS) | **installed** — esorex 3.13.10, telluriccorr 4.3.3_4 (Homebrew ESO tap); recipes `molecfit_model` / `molecfit_calctrans` / `molecfit_correct` |
| per-night **GDAS** data | observation-night atmospheric profile | **installed for Paranal** — `gdas_profiles_C-70.4-24.6.tar.gz` under `…/telluriccorr/4.3.3_4/share/molecfit/data/profiles/gdas/`. La Silla / CTIO tarballs **not present** → those instruments loud-fail at `resolve_gdas_profile` until fetched (correct behaviour — no standard-atm fallback). |
| **cr2res** (CRIRES+ DRS) | raw → 1D extraction prerequisite | **not installed**. Only needed when starting from raw; the Vesta set is already cr2res `EXTRACTC` 1D IDPs, so molecfit can run directly. Install: ESO `cr2re` pipeline (esorex tap) when a raw-extraction case arrives. |
| **APERO + Wapiti** (SPIRou) | SPIRou extraction + telluric | **not installed** (`import apero` / `import wapiti` fail). SPIRou is reduced by its own DRS; the registry routes SPIRou to `apero_wapiti` and `engine_available` reports NOT INSTALLED. Install per the APERO docs + Wapiti when a SPIRou dataset is conditioned. |

`python -m pipeline.telluric_stage --availability` prints the live routing + install
table.

## Backfill — red-optical / IR datasets still to pass the stage

- **Solar full-spectrum (RYA-371) Phase B** — the Vesta CRIRES+ K-band CO arm: conditioned
  here, but still flagged (residual > tol + RYA-373 PROVISIONAL gaps). Re-validate after
  the FTS solar IR atlas (RYA-162) and the RYA-387 0.001 re-extraction land.
- **alpha Cen A/B** CRIRES + NIRPS IR (RYA-423/432) — NIRPS routes to molecfit @ La Silla
  (GDAS tarball needed); CRIRES @ Paranal ready.
- **Procyon** IR holdings (RYA-425) — the clean single-star proving ground for the stage.
- Every **red arm** red-ward of 6800 Å (UVES-red, ESPRESSO-red, FEROS-red, CHIRON) as those
  datasets are ingested.

## Engines / steps that fail loud (never faked)

`EngineNotSelectedError` (unknown instrument) · `TelluricEngineNotInstalledError`
(molecfit/APERO absent) · `GDASUnavailableError` (no real per-night profile / silent
fallback) · `TelluricNotVerifiedError` (abundance path on an unverified red-optical/IR
frame).
