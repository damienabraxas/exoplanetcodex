# Measured solar EW results — committed, versioned line source

## `sol_ew_results_v1.csv` (806 resolved measured solar lines)

The **canonical, citable line source** for downstream tasks that must query a measured
line set without depending on the transient pipeline output (`data/processed/solar_ew.csv`
is gitignored). Columns:

```
element, ion, wavelength_air_A, ew_mA, ew_err_mA, profile_type, chi2, blend_flag, notes
```

### Provenance
- **Origin:** the published versioned measured-results set, `exoplanetcodex-site/assets/data/sol_ew_results_v1.csv` (v1).
- **Vendored into the code repo (RYA-396)** so scrapers/grids cite an in-repo committed file rather than the transient `solar_ew.csv`.
- **Currency verified (RYA-396):** matched against the latest pipeline run
  `data/processed/solar_ew.csv` (8951 fit rows, 2026-06-15). 802 of 806 v1 lines join
  on (element, ion, wavelength@0.01 Å) with **median |ΔEW| = 0.000 mÅ, max 0.000** — v1
  is a *resolved subset of the same run*, not stale. The 4 unmatched are wavelength-rounding
  edge cases.
- 22 elements (Al, Ba, C, Ca, Cr, Cu, Eu, Fe, Li, Mg, Mn, Na, Ni, O, P, S, Si, Sr, Ti, V, Y, Zr).
  Note: **no Co** measured lines (so Co NLTE grids cannot be queried from this set).

### Use
- **MPIA NLTE scrapers** (`build_nlte_grids_mpia.py`, param-driven): resolved
  `wavelength_air_A` for the element's lines.
- **INSPECT NLTE scrapers** (`fetch_inspect_nlte.py`, EW-mode): `ew_mA` at those lines.

Regenerate / re-version only from a committed pipeline run; never overwrite in place
(bump the version suffix), so the source stays citable.
