# Reproduction workflows and status

Run every workflow from the repository root after saving the output of:

```bash
python scripts/validate_installation.py --full \
  --json results/tables/install-validation.json
python -m pip freeze > results/tables/environment-freeze.txt
git rev-parse HEAD
```

Commands below describe live entry points, but the data-dependent paths were not
made portable by this documentation change. Confirm their current `--help` and
preserve the full terminal log.

## Solar benchmark

Status: most developed benchmark; requires external reflected-solar HARPS FITS,
ATLAS9.Castelli, iSpec/Turbospectrum, molecular/atomic inputs, and optionally
MOOG for the EW engine comparison.

```bash
python pipeline/spectra_normalize.py solar
python pipeline/lines_fit.py solar
python pipeline/abundances_derive.py solar ATLAS9.Castelli --skip-convergence
```

Expected outputs are under `data/processed/` and `results/`. Runtime ranges from
minutes for normalization/EW work to hours for synthesis depending on line
selection and CPU. Consult `docs/pipeline_review_solar.md`: EW and synthesis
paths have method-specific limitations, and a completed command is not by itself
a science acceptance verdict.

## Procyon FGK benchmark

Status: active research workflow; requires barycentric HARPS Phase 3 products,
ATLAS9.Castelli, and the same engine assets.

```bash
python pipeline/spectra_normalize.py procyon
python pipeline/lines_fit.py procyon
python pipeline/abundances_derive.py procyon ATLAS9.Castelli --skip-convergence
```

The loader asserts `SPECSYS=BARYCENT`. Procyon lies near correction-grid/model
edges in some methods; out-of-grid results must remain flagged. Expect
minutes-to-hours and inspect the per-line tables, Fe I/Fe II balance, reduced-EW
slope, synthesis fit quality, and NLTE coverage.

## Exoplanet host: 55 Cancri A

Status: **not yet end-to-end reproducible**. It is the primary approved host
target in `config/stars.yaml`, with line data and an acquisition demo, but
`run_pipeline.py --star 55cancri` currently fails at normalization because that
module supports only solar and Procyon. Do not publish the synthetic acquisition
demo as an abundance result. Required future work includes an approved spectrum,
manifest, 55 Cancri normalization routing, and a validated method-policy run.

## M dwarf

Status: atmosphere interpolation and loader scaffolding only; no approved
benchmark reproduction. Use MARCS.GES rather than ATLAS9.Castelli, molecular
lists become essential, and continuum placement/line blending differ
substantially from FGK stars. Before a science run, add an approved target record,
instrument loader/profile, spectrum manifest, MARCS in-grid test, molecular-data
manifest, method policy, and benchmark tolerances. Never extrapolate an FGK
correction grid into an M dwarf.

## Reduced versus full runs

The repository-native reduced smoke path is:

```bash
python scripts/validate_installation.py
python scripts/quickstart_example.py
```

It validates canonical registries and performs one solar O I correction-grid
evaluation in seconds. It does not read a spectrum or atmosphere.

A full run includes spectrum acquisition/checksums, normalization, line
measurement or flux fitting, atmosphere interpolation, LTE radiative transfer,
applicable NLTE/3D corrections, diagnostics, uncertainties, and an output
manifest. Full runs are not represented by a single stable top-level command
for all targets today.
