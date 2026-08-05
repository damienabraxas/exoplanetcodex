# Canonical scientific registries

> **Start at [`LEDGERS.md`](../../LEDGERS.md)** (repo root) — it names the full
> canonical read-set, what each ledger owns, and when to update it. The catalogs
> below are members #3–#5 of that read-set.

These machine-readable artifacts are the public metadata layer:

- `system_catalog.csv` — one row per stellar system; target identity and lifecycle.
- `instrument_catalog.csv` — one row per instrument/atlas capability; archives,
  wavelength coverage, preprocessing, and Codex support state.
- `instrument_modes.csv` — mode/grating rows where a broad instrument envelope
  would hide scientifically material coverage or resolving-power differences.
- `holdings_manifest_registry.csv` — normalized, portable joins from systems and
  instruments to repository evidence manifests; it contains no local file paths.
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

## Companion schemas

`instrument_modes.csv` uses `mode_id` as its primary key and `instrument_id` as
a foreign key. Coverage remains numeric nanometers. A mode row must stay within
the documented capability of its parent source.

`holdings_manifest_registry.csv` uses `holding_id` as its primary key and joins
`system_id` to `system_catalog.star_params_key` and `instrument_id` to the
instrument register. `manifest_path` is repository-relative; `evidence_state`
is `verified`, `audited`, `candidate`, or `rejected`. This registry asserts only
that the named manifest contains the cited evidence, never a live file count.
