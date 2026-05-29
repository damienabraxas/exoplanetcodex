# config/constants.py

Single source of truth for all pipeline constants, parameters, and file paths.
No hardcoded values anywhere else in the pipeline — import from here.

## Usage

```python
from config.constants import PHYSICS, ASTRO, SOLAR_ASPLUND2021, STAR_55CNC, PIPELINE, PATHS, MODEL

# Physical constants
c = PHYSICS['c_kms']                  # 299792.458 km/s

# Stellar parameters
teff = STAR_55CNC['teff_K']           # 5196.0 K

# Solar reference abundances
fe_solar = SOLAR_ASPLUND2021['Fe']    # 7.46 (Asplund 2021)

# File paths (Path objects, auto-created on import)
spectra_dir = PATHS['raw_spectra']
plots_dir   = PATHS['plots']

# Validate everything is intact
from config.constants import validate_constants
validate_constants()
```

## Sections

| Dict | Contents | Source |
|------|----------|--------|
| `PHYSICS` | c, h, k, G, σ, amu | NIST CODATA 2022 |
| `ASTRO` | Solar mass/radius/luminosity, AU, parsec | IAU 2015 nominal values |
| `SOLAR_ASPLUND2021` | A(X) for 25 elements | Asplund et al. 2021, A&A 653, A141 |
| `STAR_55CNC` | Teff, log g, [Fe/H], RV, distance, etc. | von Braun et al. 2011, ApJ 729, 63 |
| `PIPELINE` | Wavelength range, S/N limits, EW thresholds | Internal — see methodology.md |
| `PATHS` | All data/results directory and file paths | Auto-resolved from repo root |
| `MODEL` | ATLAS9 configuration | Castelli & Kurucz 2004 |

## Notes

**Solar Fe abundance:** `SOLAR_ASPLUND2021['Fe'] = 7.46` vs. 7.50 in Lodders (2003),
which was used in the 2010 thesis. This 0.04 dex offset is documented in the
uncertainty budget and propagates into all [Fe/H] comparisons with the thesis.

**Path creation:** All directory paths under `PATHS` are created automatically
when `config.constants` is imported. File paths (linelist_master, solar_reference,
atlas9_grid) are not created — they must be populated by the data acquisition steps.

**Updating constants:** Bump `__version__` whenever any value changes. All changes
should be traceable to a cited source in the comment on that line.
