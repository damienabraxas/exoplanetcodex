# Spectra, line data, and directory layout

## Spectral inputs

The implemented normalization paths expect reduced, merged one-dimensional
HARPS products in FITS binary tables with `WAVE`, `FLUX`, and `ERR` arrays.
Wavelengths are Angstroms. Procyon inputs must declare `SPECSYS=BARYCENT`;
the loader refuses to apply a second barycentric correction. The solar path
records and applies BERV. Resolution used for synthesis is `HARPS_R=115000`.

Raw detector frames and order extraction are outside this repository's scope.
Do not substitute an echelle order file, vacuum wavelength grid, normalized
flux, or different resolution without adding an explicit loader/instrument
profile and validation.

| Target | Source/status | Runtime location |
|---|---|---|
| Solar reflected-light calibration | HARPS direct solar feed; observation provenance in code | `$REPO_PARENT/data/spectra/exoplanetcodex-data/Solar Calibration/archive/*.fits` |
| Procyon | ESO HARPS Phase 3 ADPs; independently retrieve | `$REPO_PARENT/data/spectra/exoplanetcodex-data/Procyon/Procyon Harps/*.fits` |
| 55 Cancri A | ESO archive search by HD 75732; no approved redistributable quick-start spectrum | `data/raw/` for acquisition experiments |
| SPIRou | Loader scaffold exists | No validated reproduction |
| STIS/CRIRES/reference atlases | Audit/indicator-specific paths | See `data/audit/` and `data/solar_reference/` |

Retrieve ESO products through the
[ESO Science Archive](https://archive.eso.org/scienceportal/home) under its data
policy. Never commit restricted or very large spectra. Create a manifest:

```bash
find /absolute/path/to/spectra -type f -name '*.fits' -print0 \
  | sort -z | xargs -0 shasum -a 256 > spectra.sha256
```

Record archive observation identifiers, product category, reduction-pipeline
version, wavelength frame (air/vacuum and topocentric/barycentric), flux
convention, resolving power, S/N, retrieval date, URL, and checksum. FITS files
are gitignored.

There is no approved spectrum bundled for a publishable end-to-end quick start.
The repository-native quick start deliberately tests a vendored correction grid
instead. This is a known reproducibility gap, not hidden demo data.

## Atomic and molecular data

Target CSV line lists live in `data/linelists/`; naming follows
`linelist_<target>.csv`. `linelist_master.csv` and
`canonical_gf.csv` are central stores. VALD raw deliveries are retained where
licensing permits, with NIST cross-check files and JSON provenance for curated
products. Read [the line-list pipeline](../linelist_pipeline.md) before updates.

VALD3 is the primary extraction source and requires user registration at
[VALD](https://vald.astro.uu.se/). NIST ASD supplies independent atomic-data
checks via the [NIST Atomic Spectra Database](https://physics.nist.gov/asd).
The code preserves species/ion, wavelength, excitation potential, oscillator
strength, damping, and provenance fields. Solar/astrophysical `log gf` changes
must be made through `canonical_gf.csv` and the migration/audit tooling; never
edit only one generated target list.

Turbospectrum molecular lists are consumed from
`$ISPEC_DIR/input/linelists/turbospectrum/molecules/`. The repository also
contains target-specific ExoMol downloads under `data/linelists/molecular/`;
their presence does not make them the live Turbospectrum bundle. Record database,
isotopologue, line-list release, URL, license, checksum, conversion code, and
output format for every update.

Changing atomic or molecular data changes results. Make the update in its own
commit, regenerate derivatives, run `scripts/check_stewardship.py`, run the
tests, and compare benchmark outputs before approval.

## Layout and local overrides

```text
config/stars.yaml              canonical target parameters
data/method_policy.yaml        EW/synthesis selection
data/linelists/                versioned atomic/selected molecular data
data/nlte_grids/               versioned correction grids
data/raw/                      local FITS (ignored)
data/processed/                generated spectra/EWs (ignored)
data/model_atmospheres/        documentation only; live packs are in iSpec
results/tables/                generated tables/manifests
results/plots/                 diagnostics
```

`ISPEC_DIR` is the only supported environment-variable path override. Spectrum
locations are presently encoded in `config/constants.py` relative to the
repository parent; there is no supported local YAML override. If your layout
differs, symlink the expected target directories or propose a configuration
change—do not add a workstation-specific absolute path to version control.

There is no live `data/catalog/system_catalog.csv`. `config/stars.yaml` owns
stellar parameters, while `data/method_policy.yaml` owns method selection.
