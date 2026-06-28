# Kitt Peak Solar Flux Atlas — diagnostic segments (RYA-459, under RYA-162)

**Provenance: MEASURED.** The keystone solar reference. FTS at the McMath/Pierce
Solar Telescope, NSO/Kitt Peak — Kurucz, Furenlid, Brault & Testerman 1984,
*"Solar Flux Atlas from 296 to 1300 nm"* (National Solar Observatory Atlas No. 1).

- Source: https://nispdata.nso.edu/ftp/pub/atlas/fluxatl/ (251 segments `lm0296`–`lm1296`)
- **Raw atlas lives OUTSIDE the repo** (43 MB): `<repo-parent>/data/spectra/exoplanetcodex-data/Solar Calibration/Kitt Peak Flux Atlas/`.
  Re-fetch + re-extract with `scripts/intake_solar_atlases_rya459.py`.
- Full provenance (citation, per-segment inventory, md5, coverage): `kpno_provenance_rya459.json`.

## Columns (per `kpno_*.csv`)

| column | meaning |
|--------|---------|
| `wavelength_air_A` | air wavelength (Å), native uneven FTS grid (raw is nm; ×10 here) |
| `wavelength_vac_A` | vacuum wavelength (Å), Birch & Downs 1994 / Edlén |
| `residual_flux` | pseudo-residual = **continuum-normalized** flux (0–1) |
| `irradiance_uW_cm2_nm` | observed absolute irradiance (µW/cm²/nm) |

**Normalization:** `residual_flux` is normalized (use for line work); `irradiance` is
absolute flux. Acknowledgement required on use: *"NSO/Kitt Peak FTS data used here were
produced by NSF/NOAO."*

## Extracted diagnostic segments

These are the lines HARPS-VIS (380–690 nm) cannot reach — the RYA-369 solar-N
channels + the CNO arms + the P/K/Co/Sc DATA-GAP probes. Each is a measured,
high-resolution (Δλ ≈ 0.004–0.013 Å) window:

- **N (the unblock):** `kpno_NH_3360.csv`, `kpno_CN_violet_3883.csv`,
  `kpno_NI_7442_7468.csv`, `kpno_NI_8216_8223.csv`, `kpno_NI_8680_8718.csv`
- **CNO:** `kpno_OI_6300.csv` (forbidden cross-check), `kpno_OI_777_triplet.csv` (primary O)
- **DATA-GAP probes:** `kpno_KI_7665_7699.csv`, `kpno_CoI_3845.csv`,
  `kpno_ScII_4246.csv`, `kpno_PI_10581_10596.csv` (P near-IR multiplet — an
  alternative to the FUV P lines that need HST/STIS, RYA-119)
