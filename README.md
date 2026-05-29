# The Exoplanet Codex

> Open-science pipeline for measuring stellar elemental abundances from high-resolution spectra

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Science Questions

- What is the detailed elemental chemistry of 55 Cancri A (host to a super-Earth at 0.015 AU)?
- How do host-star abundances inform rocky planet composition and habitability?
- Can we measure [Fe/H], [Mg/Fe], [Si/Fe], [Ca/Fe], [O/Fe] to ±0.05 dex precision?

## Pipeline Steps

| Script | Step | Description |
|--------|------|-------------|
| `pipeline/01_acquire.py` | Acquire | Download & inspect HARPS S1D spectrum from ESO archive |
| `pipeline/02_normalize.py` | Normalize | Continuum normalization via polynomial fitting |
| `pipeline/03_linelist.py` | Line list | Load and filter VALD3 + NIST line list |
| `pipeline/04_ew_measure.py` | EW measurement | Fit Voigt/Gaussian profiles, measure equivalent widths |
| `pipeline/05_abundances.py` | Abundances | MOOG/ATLAS9 stellar atmosphere modeling, LTE abundance analysis |
| `pipeline/06_interpret.py` | Interpret | Solar-normalized abundances, plots, comparison with literature |

## Repo Structure

```
exoplanetcodex/
├── config/
│   ├── constants.py          # Single source of truth: physics, stellar params, pipeline settings
│   └── README.md
├── pipeline/
│   ├── 01_acquire.py
│   ├── 02_normalize.py
│   ├── 03_linelist.py
│   ├── 04_ew_measure.py
│   ├── 05_abundances.py
│   └── 06_interpret.py
├── data/
│   ├── raw/                  # Downloaded FITS files (gitignored)
│   ├── processed/            # Intermediate pipeline outputs (gitignored)
│   ├── linelists/            # VALD3 + NIST master CSV
│   └── model_atmospheres/    # ATLAS9 model grid
├── results/
│   ├── plots/                # Generated figures (gitignored)
│   └── tables/               # Abundance output tables (gitignored)
└── tests/
```

## Quick Start

```bash
git clone https://github.com/damienabraxas/exoplanetcodex.git
cd exoplanetcodex
pip install -r requirements.txt

# Run step 1 in demo mode (no FITS data needed)
python pipeline/01_acquire.py

# Run tests
pytest tests/ -v
```

## Getting Real HARPS Data

1. Create a free account at the [ESO User Portal](https://www.eso.org/UserPortal)
2. Visit the [ESO Science Archive](http://archive.eso.org/wdb/wdb/adp/phase3_main/form)
3. Search: Target = `HD 75732`, Instrument = `HARPS`, Data Product = `SCIENCE`
4. Download the `ADP.*S1D*.fits` files (merged 1D spectra)
5. Place files in `data/raw/`
6. In `pipeline/01_acquire.py`, switch from demo mode to Option A (real FITS)

## Contributing

This is an open-science project. Issues and PRs welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

If you use this pipeline, please cite the underlying data sources:
- HARPS spectra: ESO archive
- Line list: VALD3 (Ryabchikova et al. 2015) + NIST ASD
- Solar abundances: Asplund et al. (2021), A&A 653, A141
- Stellar parameters: von Braun et al. (2011), ApJ 740, 49

## License

MIT — see [LICENSE](LICENSE)
