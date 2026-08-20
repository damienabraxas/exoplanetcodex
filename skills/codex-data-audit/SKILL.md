---
name: codex-data-audit
description: Determine whether an observational dataset is scientifically safe and appropriate to enter the pipeline without altering evidence or performing abundance science. Use for every acquisition, delivery, staging, first measurement, DATA AUDIT, or "can we measure X on this data" request. Produces provenance, frame, telluric, coverage, SNR, and GO/CAVEAT/NO-GO dispositions.
---

# Observational data intake and audit

## Boundary and authorities

This skill establishes whether held data may enter a requested pipeline path. It does
not normalize or rewrite source evidence, fit lines, derive abundances, or choose a
scientific answer. Start at `LEDGERS.md`. Read the holdings registry before any download,
resolve instrument/mode rules from current catalogs and instrument docs, and resolve
numeric gates such as `PIPELINE['snr_min_science']` from live config at runtime. Report
the resolved value and source; never copy a mutable threshold into this skill.

## Route before touching data

1. Check `data/catalog/holdings_manifest_registry.csv` by system, instrument, product,
   and content identity.
2. Known holding: audit it. New delivery: intake and audit, then register it. Missing
   required product: report an acquisition task; do not substitute or re-download.
3. Preserve originals read-only. Derived inventories carry source checksums.

## Required identity and provenance inventory

For every distinct product record, where available: source/archive and archive/product
ID; original/local filename, byte size, checksum; observation timestamp; target;
instrument/mode; product class; reduction software/version; wavelength frame and
units/medium; and flux-column semantics. Determine these from content and authoritative
metadata, not paths or sibling tickets. Detect duplicates by checksum/content identity.
An absence needs a positive control showing the lookup would detect the field.

## Audit axes

### Product and identity

Inspect authoritative headers/metadata (`OBJECT`, instrument/mode, observation time,
product category, pipeline/version). Verify content signatures and arrays, not merely
file presence. Record ambiguity; never default object identity.

### Wavelength frame and units

Determine the delivered frame before deciding whether velocity correction is owed.
Inspect `SPECSYS` and credible velocity metadata, then validate against the applicable
loader/frame contract. Missing BERV is not permission to guess or compute one. Reflected
sources may require a two-leg source-body-observer model. Ambiguous frame metadata is a
STOP until validated. Treat units and air/vacuum medium per file.

### Telluric state when applicable

First resolve whether tellurics apply and what basis is allowed from the instrument
catalog. Determine `telluric_applied` as `applied`, `not-applied`, or `unknown` from the
full recipe chain, transmission extensions, and flux columns via
`pipeline.telluric_intake`. Record evidence and the required corrected column. A
transmission array of ones is not an applied correction; a file may carry corrected and
uncorrected columns together. `unknown` is real and never defaulted. Keep instrument
`telluric_basis` separate from per-holding `telluric_applied`.

### Science coverage

Resolve required indicators/regions from the active task, element protocol, line
accounting, and model-domain requirements. Test valid real pixels (quality mask, finite
nonzero flux), including gaps; header endpoints do not prove coverage. Healthy data
without required coverage yield `NOT-APPLICABLE`/`DATA-GAP`, never a fabricated result.

### Quality and contextual metadata

Measure SNR/normalization using current canonical gates and report the resolved
threshold/source. Verify the selected flux column's normalization state. Airmass,
transparency, and similar metadata are contextual flags unless a current authority
explicitly makes one a gate.

## Verdict and writeback

Produce a table covering identity/product, provenance, duplicate status, tellurics,
frame/units, coverage/indicators, SNR/normalization, and an overall `GO`, `CAVEAT`,
`NO-GO`, or `NOT-APPLICABLE/DATA-GAP`. Every caveat names its consequence.

For a new or newly determined holding, update the holdings registry and evidence. This
is a state change: follow `skills/codex-state-register/SKILL.md`. Route the report to the
active issue/current ledger, not a hard-coded historical ticket.

## Non-negotiable checks

- Verify tests discriminate using positive controls; enumerate filtered false positives.
- Assert on content, not existence, and preserve source checksums.
- Never infer provenance, frame, tellurics, or units from memory.
- Never perform abundance science during intake.
- Stop on unresolved scientifically material ambiguity.

Run applicable instrument/loader tests selected from current docs and code. Confirm an
IR holding marked `not-applied` or `unknown` is refused by the live telluric gate.
