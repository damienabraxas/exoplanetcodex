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

## Canonical solar-EW sources — the three roles (RYA-397, updated RYA-408)

There is no `sol_ew_results_v1_1.csv`; **`sol_ew_results_v1.csv` is the single canonical
committed set** read by the solar Fe gate, abundance derivation, and stewardship. The
three solar-EW artifacts are distinct:

| File | Role | Tracked? | Referenced via |
|------|------|----------|----------------|
| `data/measured/sol_ew_results_v1.csv` | **CANONICAL** committed measured-line set (806 lines, 22 elements) — the single source of truth for the Fe gate, `abundances_derive`, the NLTE scrapers/grids, and non-Fe pool curation. | committed | `PATHS['solar_ew_canonical']` |
| `data/processed/solar_ew.csv` | **STAGING ONLY** — the per-run `lines_fit.py` output. Regenerable, gitignored, **never a gate/abundance EW input** (RYA-408: this was the RYA-406 incident). Promoted to the canonical only via the reviewed `scripts/promote_solar_ew.py`. | gitignored (`data/processed/`) | `PATHS['solar_ew']` |
| `data/processed/solar_ew_ges_reference.csv` | **GES Fe I EW reference** — pre-stored GES EWs used for the solar Fe I leg only (avoids the NLTE EW bias the measured Fe I carries; RYA-330). | committed (force-added) | `…/solar_ew_ges_reference.csv` |

**RYA-408 contract:** no gate or abundance derivation may take a gitignored, regenerable
file as its EW input. `abundances_derive._load_solar_ews` reads the canonical
(`PATHS['solar_ew_canonical']`) for the Fe II + non-Fe-I pool; Fe I stays on the committed
GES reference. The stewardship invariant `check_solar_ew_canonical` enforces this
(canonical present + well-formed; the loader actually reads it; and a present staging file
must not drift from the canonical on a measured EW). A dated
`results/Solar/<date>/solar_ew.csv` snapshot is NOT a source.

### Promotion (staging → canonical)
Re-version only from a committed, reviewed pipeline run via `scripts/promote_solar_ew.py`
(dry-run by default; `--apply` writes). Promotion **STOPS** — writing nothing — on any
`blend_flag` conflict between staging and the canonical: the canonical's 11 vetted
blends (incl. **O I 6300.304**, Ni-blended per RYA-104/208) are curation that the raw
`lines_fit` staging does not own and must never silently overwrite. EW measurements are
promoted; canonical coverage is preserved; new staging-only lines are reported, not
auto-added. A genuine `lines_fit` flag bug is a follow-up (fix `VETTED_BLENDS`, the single
source of `blend_flag`), not a promotion-time overwrite.
