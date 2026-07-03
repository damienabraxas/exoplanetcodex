# RYA-509 — Solar anchor clean-from-raw reproducibility (RYA-506 blast-radius check)

Branch `ryandamienschmitt/rya-509-solar-reanchor`. Verify-only; no tuning, no filter
loosening, no code changes, no merge. Runner: `scratch/rya509_solar_reanchor.py`.
Reference = the BANKED gold verdict `data/reference/solar/solar_abundances_v1.csv`
(via `ns.read_solar_reference('CURRENT')`), not a memorized value.

## Verdict — the solar anchor REPRODUCES FROM RAW, EXACTLY
Full from-raw chain re-ran on the Mac through the RYA-506-fixed pipeline (fork-forced
theo-EW + all-zero guard), single-thread BLAS (`OMP/OPENBLAS/MKL_NUM_THREADS=1`):

    raw FITS → spectra_normalize → lines_fit (fresh EW measure) → abundances_derive

Reproduces the banked **v1** solar verdict on **all 13 species with Δ = 0.000 dex** (LTE
*and* NLTE), identical n_lines (full table: `mac_vs_banked_v1.csv`). Headline:
- **Fe I:** A(X)=7.510 (1D-LTE) / **7.516 (NLTE)**, n=62 — Δ 0.000 vs banked.
- **Fe II:** 7.657 / 7.657, n=3 — Δ 0.000.
- Solar Fe VERDICT: **PASS** (slope · ionization · scatter). max |ΔA_NLTE| = 0.000 across
  all species.
- The RYA-506 all-zero-batch guard did **NOT** trip → the solar synthesis genuinely
  survives, it does not merely pass the gate.

**Contrast with Procyon (RYA-506):** Procyon rode persistent state and MOVED on a forced
rebuild (7.593 → 7.571); the Sun did not move at all. The solar verdict was always sound —
the anchor is clean.

## Confirm-item — the EW stage genuinely re-measures from raw (not a cached fallback)
`data/processed/solar_ew.csv` (8953 rows, all 27 elements) md5 = `230c27c4a9da5770dc6db884709b448b`.
Deleting it and re-running `lines_fit --star solar` from the FITS regenerates the **identical**
md5 `230c27c4…` — deterministic from-raw regeneration, so the measurement stage genuinely
ran (Δ=0.000 on abundances is necessary; this hash proves the *measurement* was from raw).
`ad.run('solar')` consumes `solar_ew.csv` via `_load_solar_ews` (the committed
`solar_ew_reference.csv` is a 2.7 KB curation snippet, not the measured table) → no silent
fallback path.

Stage-artifact hashes (Mac; recorded for the future Mac-vs-Sirius stage compare, RYA-511):
- `solar_normalized.csv` (312,744 rows) md5 = `df4a49cf5d69876628ebca5f2ef8cdfc`
- `solar_ew.csv` (8953 rows) md5 = `230c27c4a9da5770dc6db884709b448b`

## Honest scope of what 509 verifies
This is a **reproducibility gold-standard**: from raw, deterministic, no persistent-state
dependence. It is **NOT YET a cross-environment-robustness** gold-standard — that badge is
RYA-511 (once Sirius can run and the two environments agree). The Solar website page must
not imply cross-confirmed until 511 lands.

## Restart findings (verify-only respected — no code changed here)
1. **fork coverage is partial.** `lines_fit.py` and `spectra_normalize.py` use no
   multiprocessing and call none of the iSpec child-spawning functions (those live in
   `ispec/…/{abundances,sme,spectrum}.py`, invoked only in the derive stage), so the
   EW/normalize stages are not exposed to the RYA-506 spawn/fork bug. BUT the fork
   start-method is set only at `abundances_derive` import — standalone entries
   (`lines_fit`, `spectra_normalize`) run under the platform default. Latent fragility →
   **RYA-514** (centralize the fork setting).
2. **Sirius is not self-sufficient.** Reachable, 233 GB free, but has **no repo code, no
   iSpec, and no py3.14 numpy/pandas/scipy/astropy** — only staged data
   (grids/linelists/solar_reference/spectra). The authoritative leg cannot run without a
   build → **RYA-511 Phase 0**.
3. **Sirius py3.14 defaults to `forkserver`, not `fork`.** The "Sirius = native fork"
   premise is dead — forkserver re-imports like spawn, so the RYA-506 class applies on
   Sirius too and fork must be forced there as well → RYA-511 scope correction + RYA-514.
