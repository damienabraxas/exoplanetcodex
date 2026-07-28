# Troubleshooting

Start with:

```bash
python scripts/validate_installation.py
python scripts/validate_installation.py --full
python -m pytest -q
```

## Common failures

`ModuleNotFoundError: config` or `pipeline`: run commands from the repository
root and activate the intended environment.

`ModuleNotFoundError: ispec`: set `ISPEC_DIR` to the source checkout and follow
the upstream iSpec bootstrap. The similarly named PyPI package is not the tested
installation.

ATLAS9.Castelli or MARCS.GES missing: unpack the input bundle matching the
recorded iSpec commit. Raw upstream model files do not automatically satisfy
iSpec's pack format. Re-run the full validator and atmosphere interpolation
tests.

`bsyn_lu` missing/not executable or Turbospectrum exits: compile the engine using
iSpec's instructions; verify architecture, Fortran runtime, execute permission,
and the matching opacity/molecular data. Preserve the full engine log.

MOOG unavailable: it affects the EW baseline/comparison, not the repository
correction-grid quick start. The expected iSpec location is
`$ISPEC_DIR/synthesizer/moog/MOOGSILENT`.

No FITS files found: place independently retrieved products at the exact path
listed in [the spectra guide](data/spectra.md), then save a checksum manifest.
Check product format and filename suffix.

Procyon `SPECSYS` failure: obtain barycentric Phase 3 products. Do not remove the
guard or apply BERV twice.

55 Cancri `NotImplementedError` in normalization: this is a known unsupported
full path, despite the available acquisition demo. Use the quick-start smoke
example or implement/test target routing in a separate change.

STAR_PARAMS/canonical configuration failure: edit `config/stars.yaml`; do not
add a second parameter copy. Run:

```bash
python -c "from config.constants import validate_constants; validate_constants()"
python scripts/check_stewardship.py
```

Missing line or oscillator-strength divergence: follow
`docs/linelist_pipeline.md`, update the canonical store, regenerate target lists,
and run stewardship/GF consistency tests. Never patch only a derived CSV.

NLTE correction is NaN/unavailable: the star may be outside the grid convex hull
or the transition may have no wavelength match. Preserve the unavailable flag;
do not clamp, extrapolate, or silently label the LTE value NLTE.

C/O correction sign or scale guard failure: confirm Angstrom versus nanometre
units, absolute `A(X)` versus `[X/H]`, table/leg selection, and the correction
convention `A(non-LTE) - A(1D-LTE)`.

Multiprocessing crash/hang: reproduce serially with one process. Give each job a
private temp directory and avoid sharing mutable Fortran scratch files.

LibreSSL/urllib3 warning on macOS: use a Python build linked to current OpenSSL
for network/archive work. It does not affect the offline quick-start calculation.

## What success means

A zero exit code from installation validation establishes that expected files
and imports exist. Passing tests establishes coded invariants. Neither replaces
inspection of line-level diagnostics, method applicability, correction coverage,
provenance, and benchmark tolerances for a science result.
