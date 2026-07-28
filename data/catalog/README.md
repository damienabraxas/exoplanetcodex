# Canonical scientific registries

These machine-readable artifacts are the public metadata layer:

- `system_catalog.csv` — one row per stellar system; target identity and lifecycle.
- `instrument_catalog.csv` — one row per instrument/atlas capability; archives,
  wavelength coverage, preprocessing, and Codex support state.
- `../audit/element_status_tracker.csv` — canonical element/status artifact until
  a normalized `data/catalog/elements.csv` is ratified.

Instrument capability is not evidence that a target has data. Per-system
holdings manifests join to `system_catalog.csv` by `star_params_key`/stable
system identifier and to `instrument_catalog.csv` by `instrument_id`.

## Instrument schema

Wavelengths are numeric nanometers. Resolving powers are dimensionless bounds.
Pipe-delimited values are controlled multi-value fields. `codex_status` is one
of `active`, `experimental`, `candidate`, `rejected`, or `reference`.
`pipeline_ingest_status` is one of `loader_ready`, `loader_partial`,
`audit_only`, `candidate`, `reference_only`, or `rejected`.

The catalog describes:

- capability: facility, archive, bands, wavelength range, and resolution;
- support: loader/reduction state and scientific role;
- caveats: tellurics, air/vacuum conversion, calibration/iodine concerns;
- provenance: source issues, documentation, and last verification date.

Validate with:

```bash
python -m data.catalog.instruments
```

Do not add per-star file counts or claim holdings in the instrument master.
Those belong in target-specific manifests and must reference `instrument_id`.
