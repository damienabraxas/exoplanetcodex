# Line Lists

This directory contains atomic line data for spectral analysis.

## Master Line List

File: `master_linelist.csv`

This file is not included in the repo (large, generated from external sources).
See the acquisition notes below to regenerate it.

### Column Schema

| Column | Type | Description |
|--------|------|-------------|
| `wavelength_A` | float | Air wavelength in Angstroms |
| `element` | str | Element symbol (e.g., `Fe`, `Mg`) |
| `ion` | str | Ionization stage: `I` = neutral, `II` = singly ionized |
| `excit_eV` | float | Lower-level excitation potential (eV) |
| `log_gf` | float | Log of oscillator strength × statistical weight |
| `nist_grade` | str | NIST ASD accuracy grade: A+, A, B, C, D, or E |
| `is_blend` | bool | True if flagged as potential blend in VALD3 |
| `source` | str | `VALD3`, `NIST`, or `VALD3+NIST` |

### Grade Definitions (NIST ASD)

| Grade | Accuracy | Use in pipeline |
|-------|----------|-----------------|
| A+ | < 0.3% | Yes |
| A  | < 1% | Yes |
| B  | < 3% | Yes (default minimum) |
| C  | < 10% | Optional |
| D  | < 25% | Not recommended |
| E  | > 25% | Excluded |

Default minimum grade: `B` (set in `config/constants.py` → `PIPELINE['min_nist_grade']`)

## Acquiring the Line List

### VALD3
1. Register at [VALD3](https://vald.astro.uu.se)
2. Submit an "Extract Stellar" request for 55 Cancri parameters:
   - Teff = 5196 K, log g = 4.41, [Fe/H] = +0.32
   - Wavelength range: 3780–6910 Å
   - Detection threshold: 0.01 (log(depth))
3. Download the resulting line list

### NIST ASD
1. Visit [NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database)
2. Export lines for elements of interest in the same wavelength range
3. Use `A_ki` (transition probability) to compute `log gf`

## Loader

See `data/linelists/loader.py` for the `load_linelist()` function.
See `tests/test_linelist_loader.py` for usage examples.
