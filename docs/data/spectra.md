# Spectra, line data, and directory layout

## Spectral inputs

The observational program spans UV/VIS/NIR/IR. Instrument capability and
preprocessing rules are canonical in
[`instrument_catalog.csv`](../../data/catalog/instrument_catalog.csv); the
human-readable strategy is [instruments.md](instruments.md). Per-target holdings
must come from a manifest joined to the System Register—not from a static table
on this page.

The most mature normalization path expects reduced, merged HARPS
one-dimensional FITS tables with `WAVE`, `FLUX`, and `ERR`. Other live modules
handle UVES, HST/STIS/COS, NIRPS, SPIRou, CRIRES+, and reference atlases with
instrument-specific contracts. Never substitute an echelle order, vacuum grid,
polarimetric product, normalized survey product, or different resolution
without the matching loader/profile and validation.

Retrieve products from the archive named in the instrument register under its
data policy. Never commit restricted or very large spectra. Create a manifest:

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
data/catalog/                  system and instrument registries
data/audit/element_status_tracker.csv  canonical element-status artifact
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

`data/catalog/system_catalog.csv` owns project target inventory;
`config/stars.yaml` owns runnable stellar parameters; `data/method_policy.yaml`
owns method selection. They are related but not interchangeable.
