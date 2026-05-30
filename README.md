# The Exoplanet Codex

> Open-science pipeline for measuring stellar elemental abundances from high-resolution spectra

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Science Questions

- What is the detailed elemental chemistry of 55 Cancri A (host to a super-Earth at 0.015 AU)?
- How do host-star abundances inform rocky planet composition and habitability?
- Can we measure [Fe/H], [Mg/Fe], [Si/Fe], [Ca/Fe], [O/Fe] to ±0.05 dex precision?

## Quick Start

```bash
git clone https://github.com/damienabraxas/exoplanetcodex.git
cd exoplanetcodex
pip install -r requirements.txt

# Run the full pipeline for 55 Cancri (demo mode, no FITS data needed)
python run_pipeline.py --star 55cancri

# List available targets
python run_pipeline.py --list-stars

# Run tests
pytest tests/ -v
```

## Pipeline

Scripts use **subject_action** naming: the noun tells you what data it operates on,
the verb tells you what it does. Related scripts sort together alphabetically.

| Script | Operates on | Action | Status |
|--------|-------------|--------|--------|
| `pipeline/spectra_acquire.py` | Raw spectra | Download & inspect HARPS S1D from ESO archive | Implemented |
| `pipeline/spectra_normalize.py` | Raw spectra | Continuum normalization via iterative polynomial fit | Stub |
| `pipeline/params_stellar.py` | EW table | Derive Teff/logg/vturb via excitation+ionization equilibrium | Stub |
| `pipeline/lines_load.py` | Line list CSV | Load & filter VALD3 + NIST master line list | Stub |
| `pipeline/lines_fit.py` | Normalized spectrum | Fit Voigt/Gaussian profiles, measure equivalent widths | Stub |
| `pipeline/abundances_derive.py` | EW + params | MOOG/ATLAS9 LTE abundance analysis | Stub |
| `pipeline/uncertainty_stack.py` | Abundances | Type A + Type B formal uncertainty budget | Stub |
| `pipeline/ratios_interpret.py` | Abundances | Solar-normalized ratios, Mg/Si, C/O, science plots | Stub |

All steps are wired into `run_pipeline.py` at the repo root.

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
5. Place files in `data/raw/`
6. In `pipeline/spectra_acquire.py`, switch from demo mode to Option A (real FITS)

## Contributing

This is an open-science project. Issues and PRs welcome.

If you use this pipeline, please cite the underlying data sources:
- HARPS spectra: ESO archive
- Line list: VALD3 (Ryabchikova et al. 2015) + NIST ASD
- Solar abundances: Asplund et al. (2021), A&A 653, A141
- Stellar parameters: von Braun et al. (2011), ApJ 740, 49

## License

MIT — see [LICENSE](LICENSE)
