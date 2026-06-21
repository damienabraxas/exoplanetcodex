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

## Canonical solar-EW sources — the three roles (RYA-397)

There is no `sol_ew_results_v1_1.csv`; **`sol_ew_results_v1.csv` is the current
canonical committed set.** The three solar-EW artifacts are distinct and all *live* —
none is deprecated (RYA-396 vendored the committed one and verified currency, so there
was nothing to quarantine here):

| File | Role | Tracked? | Referenced via |
|------|------|----------|----------------|
| `data/processed/solar_ew.csv` | **Transient runtime** EW table — regenerated every pipeline run by `lines_fit.py`; the working file the abundance/diagnostic path reads. | gitignored (`data/processed/`) | `PATHS['solar_ew']` |
| `data/measured/sol_ew_results_v1.csv` | **Canonical committed/citable** measured-line set (806 lines, 22 elements) — for tasks needing a stable in-repo source (NLTE scrapers/grids, non-Fe pool curation). | committed | `data/measured/sol_ew_results_v1.csv` (literal) |
| `data/processed/solar_ew_ges_reference.csv` | **GES Fe I EW reference** — pre-stored GES EWs used for the solar Fe I leg only (avoids the NLTE EW bias the measured Fe I carries; RYA-330). | committed (force-added) | `…/solar_ew_ges_reference.csv` |

Rule: query the **committed** set (`sol_ew_results_v1.csv`) from any script that must
not depend on a transient run; the pipeline's own EW path uses `PATHS['solar_ew']`. A
dated `results/Solar/<date>/solar_ew.csv` snapshot is NOT a source — repoint any such
reference to the committed set.
