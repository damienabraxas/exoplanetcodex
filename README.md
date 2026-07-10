# The Exoplanet Codex

> Open-science pipeline for measuring stellar elemental abundances from high-resolution spectra

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Science Questions

- What is the detailed elemental chemistry of 55 Cancri A (host to a super-Earth at 0.015 AU)?
- How do host-star abundances inform rocky planet composition and habitability?
- Can we measure [Fe/H], [Mg/Fe], [Si/Fe], [Ca/Fe], [O/Fe] and 22 other elements to ±0.05 dex precision?

## Quick Start

```bash
git clone https://github.com/damienabraxas/exoplanetcodex.git
cd exoplanetcodex
pip install -r requirements.txt

# List runnable targets (resolve via config/stars.yaml)
python run_pipeline.py --list-stars

# Run the pipeline for the solar calibrator (real spectra; no synthetic data)
python run_pipeline.py --star solar

# Stop after the stellar-parameter stage
python run_pipeline.py --star solar --validate-only

# Run tests
pytest tests/ -v
```

`run_pipeline.py` is a **thin driver**: it orchestrates the real stage `run()`s and
surfaces their own guards. It never generates data and never holds star parameters —
those come only from `config/stars.yaml` via `get_star_params()`. The real sequence is

```
spectra_normalize → lines_fit → [params_stellar*] → abundances_derive
                  → uncertainty_stack → ratios_interpret
```

`*params_stellar` runs **only** for stars that must *solve* a spectroscopic-equilibrium
parameter (Teff/log g/ξ — e.g. 55 Cnc ξ). Pinned calibrators (solar, Procyon, α Cen)
take their parameters straight from `stars.yaml` and skip it. Reading FITS happens inside
`spectra_normalize` (there is no separate acquire step), and linelist loading is internal
to `abundances_derive` (via `data/linelists/loader.py`). Stages that are not yet
implemented stop the run with a clear message rather than fabricating output.

## Pipeline

Scripts use **subject_action** naming: the noun tells you what data it operates on,
the verb tells you what it does. Related scripts sort together alphabetically.

| Script | Operates on | Action | Status |
|--------|-------------|--------|--------|
| `pipeline/spectra_normalize.py` | Raw spectra | Read FITS, co-add, continuum-normalize | **Real** (solar, Procyon) |
| `pipeline/lines_fit.py` | Normalized spectrum | Fit profiles, measure equivalent widths | **Real** (solar path) |
| `pipeline/abundances_derive.py` | EW + params | LTE (+NLTE) EW → A(X) abundance analysis | **Real** |
| `pipeline/params_stellar.py` | EW table | Solve Teff/log g/ξ via excitation+ionization equilibrium | Stub — RYA-537 |
| `pipeline/uncertainty_stack.py` | Abundances | Type A + Type B formal uncertainty budget | Stub |
| `pipeline/ratios_interpret.py` | Abundances | Solar-normalized ratios, Mg/Si, C/O, science plots | Stub |
| `pipeline/spectra_acquire.py` | Raw spectra | ESO/HARPS query + FITS inspection **library** | Not a pipeline stage — `run()` raises |
| `pipeline/lines_load.py` | Line list CSV | VALD3 + NIST loader wrapper | Not a pipeline stage — loading is internal to `abundances_derive` |

`run_pipeline.py` orchestrates only the real stages, in the sequence above. The two
non-stage modules keep their library functions but their `run()` raises (they are not
wired into the driver): `spectra_acquire` used to return synthetic demo data, which is
exactly what this driver removes.

## Repo Structure

```
exoplanetcodex/
├── run_pipeline.py           # Single entry point — runs all steps for a given star
├── config/
│   ├── constants.py           # Physics, stellar params, pipeline settings (single source of truth)
│   └── README.md
├── pipeline/
│   ├── spectra_acquire.py
│   ├── spectra_normalize.py
│   ├── params_stellar.py
│   ├── lines_load.py
│   ├── lines_fit.py
│   ├── abundances_derive.py
│   ├── uncertainty_stack.py
│   └── ratios_interpret.py
├── data/
│   ├── raw/                   # Downloaded FITS files (gitignored)
│   ├── processed/             # Intermediate pipeline outputs (gitignored)
│   ├── linelists/             # VALD3 + NIST master CSV
│   └── model_atmospheres/     # ATLAS9 model grid
├── results/
│   ├── plots/                 # Generated figures (gitignored)
│   └── tables/                # Abundance output tables (gitignored)
└── tests/
```

## Getting Real HARPS Data

1. Create a free account at the [ESO User Portal](https://www.eso.org/UserPortal)
2. Visit the [ESO Science Archive](http://archive.eso.org/wdb/wdb/adp/phase3_main/form)
3. Search: Target = `HD 75732`, Instrument = `HARPS`, Data Product = `SCIENCE`
4. Download the `ADP.*S1D*.fits` files (merged 1D spectra)
5. Place the files where the star's loader expects them (see `pipeline/spectra_normalize.py`)
6. Run `python run_pipeline.py --star <id>` — `spectra_normalize` reads the FITS directly

(`pipeline/spectra_acquire.py` provides the ESO query / FITS-inspection helpers as a
library, but it is not part of the pipeline run.)

## Contributing

This is an open-science project. Issues and PRs welcome.

If you use this pipeline, please cite the underlying data sources:
- HARPS spectra: ESO archive
- Line list: VALD3 (Ryabchikova et al. 2015) + NIST ASD
- Solar abundances: Asplund et al. (2021), A&A 653, A141
- Stellar parameters: von Braun et al. (2011), ApJ 740, 49

## License

MIT — see [LICENSE](LICENSE)
