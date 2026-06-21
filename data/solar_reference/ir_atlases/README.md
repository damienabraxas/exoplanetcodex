# IR reference atlases — K-band CO arm (RYA-390, under RYA-162)

Three IR reference atlases for validating RYA-373's telluric-conditioned CRIRES+ CO
output, each in a distinct role. Intaked by `scripts/intake_ir_atlases_rya390.py`;
full provenance (citation, SHA-256, coverage, conventions) in
`ir_atlases_provenance_rya390.json`.

**CO segment extracted:** 4255–4367 cm⁻¹ ≈ **22892.8–23495.4 Å air** (2.289–2.350 µm),
the CO first-overtone (Δv=2) region for the ¹³C/¹⁸O arm.

| Source | Role | Telluric | Disk | CO segment file | Band covered |
|--------|------|----------|------|-----------------|--------------|
| **ACE-FTS** (Hase, Wallace+ 2010, JQSRT 111 521) | solar truth — PRIMARY | **free** (space occultation) | ~integrated | `ace_fts_solar_co_4255_4367.csv` | full 4255–4367 |
| **NSO photatl** (Livingston & Wallace 1991, NSO TR 91-001) | terrestrial solar cross-ref | residual (cols: solar/atmospheric/total) | disk-**center** | `nso_photatl_co_4255_4367.csv` | full 4255–4367 |
| **Wallace telluric** (Wallace, Hinkle & Livingston, NOAO/NSO) | pure-telluric — validates molecfit model | **pure telluric** transmission | n/a | `wallace_telluric_co_ratio.csv` | 4299.8–4338.6 (band middle) |

## Conventions
- **Native = FTS vacuum wavenumber (cm⁻¹).** Each segment carries `wavenumber_cm-1`
  plus `wavelength_vac_A` (= 1e8/wavenumber) and `wavelength_air_A` (Birch & Downs
  1994 / Edlén — the VALD3/iSpec convention; vac−air ≈ 6.3 Å here).
- Vesta (reflected solar, the RYA-373 target) is **integrated-disk** → ACE matches it
  best; photatl is disk-center (a documented caveat).

## RYA-373 three-way validation (Part B, consumed downstream)
1. **Telluric-removal:** conditioned CRIRES CO vs **ACE-FTS** (telluric-free) → did we
   recover the true solar CO?
2. **Telluric-model:** molecfit telluric vs **Wallace** ratio → did the model match reality?
3. **Cross-instrument:** vs **photatl** (independent ground-based reduction).

## Caveats
- **Wallace ASCII telluric** (`ratio04300.txt`) covers only **4299.8–4338.6 cm⁻¹** — the
  middle of the CO band. For full-band telluric use the **photatl `atmospheric`** column
  (4248–4377) or the Wallace `atl04240`–`atl04360` eps plots; telluric line IDs in the
  Wallace `linelist_TOTAL/H2O/CH4_ext.txt`.
- **photatl `solar`** column is linearly interpolated across strong-telluric gaps
  (README); intensities are not continuous file-to-file. The `atmospheric`/`total`
  columns are the unaltered observed/telluric.

## Part B — three-way validation (`pipeline/co_validation_rya390.py`)
The validation harness consuming these atlases. Run: `python -m pipeline.co_validation_rya390`
→ `data/audit/crires_co_conditioned/rya390_co_validation.json`. **All comparisons are in
VACUUM** (CRIRES IDP / molecfit `WAVELENGTH_FRAME=VAC`); the conditioned CO is topocentric
+ RV-insufficient, so the solar checks cross-correlate to measure the reflected-solar
velocity first.

**Finding on RYA-373's provisional conditioned CO** (snapshot `data/audit/crires_co_conditioned/
vesta_crires_K_CO_K21{92,17}_topocent_PROVISIONAL.fits`): both settings are
**TELLURIC-DOMINATED** — the "corrected" product correlates **0.92/0.86 with the telluric
atlas at v≈0** but only ~0.15–0.21 with the solar atlases → telluric removal incomplete,
re-run RYA-373. The telluric-reference cross-check (Wallace vs photatl atmospheric) agrees
at **0.99**, validating the atlases + machinery. Check 2 (molecfit telluric vs Wallace) is
**BLOCKED** until RYA-373 persists the molecfit transmission (`mtrans`) in the product.

## Raw atlases (large, not in git)
`<repo-parent>/data/spectra/exoplanetcodex-data/Solar Calibration/IR Reference Atlases/`
— `ACE-FTS/` (ace-solar-spectrum.txt 16 MB + line lists), `NSO_photatl/` (wn4250–wn4350
+ README), `Wallace_telluric/` (ratio04300.txt + line lists + README.pdf). Re-fetch with
`scripts/intake_ir_atlases_rya390.py` (URLs in the provenance JSON).
