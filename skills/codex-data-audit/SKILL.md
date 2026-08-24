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

### Normalisation state — MANDATORY, from the flux

Determine `normalization_state` as `normalised`, `un-normalised`, or `unknown` from the
FLUX via `pipeline.normalization_intake`, and record it per holding beside
`telluric_applied`. This is not optional and it is not satisfied by reading a label: a
filename, a README, or a `pre_normalised` flag are all the same claim written down, and a
declared flag plus a mis-routed file agree with each other perfectly while both being
wrong. **Scan what the READER RETURNS**, not what a directory listing shows — the failure
this catches is a product declared normalised whose loader opens the raw file.

Cross-check the detected state against the declared flag. **Disagreement is a LOUD STOP,
never an auto-fix and never an auto-skip**: a genuinely raw spectrum with a coincidentally
flat continuum and a normalised one with a bad blaze both need a human, and the detector
informs the declared state rather than overriding the science. A product carrying no
declaration is UNDECLARED and must be scanned before any continuum stage runs.

`unknown` is real and never defaulted — defaulting to normalised applies unity as a
continuum and inflates every EW silently. The envelope says nothing below the
continuum-limited blue edge (near-UV blanketing leaves no true continuum in any window)
or inside a registered telluric band, and the module reports `unknown` there rather than
accusing a good product. Fill values are not flux.

Keep this axis separate from `telluric_applied` (RYA-806) and `observed_conditioning`
(RYA-1006). Three conditioning axes, three columns, never collapsed.

### Science coverage

Resolve required indicators/regions from the active task, element protocol, line
accounting, and model-domain requirements. Test valid real pixels (quality mask, finite
nonzero flux), including gaps; header endpoints do not prove coverage. Healthy data
without required coverage yield `NOT-APPLICABLE`/`DATA-GAP`, never a fabricated result.

### Quality and contextual metadata

Measure SNR using current canonical gates and report the resolved threshold/source.
Normalisation has its own mandatory axis above — do not re-answer it here from a label. Airmass,
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
- Never infer provenance, frame, tellurics, units, or normalisation from memory.
- Never conclude a product is normalised from a filename, a README, or a flag. Scan it.
- Never perform abundance science during intake.
- Stop on unresolved scientifically material ambiguity.

Run applicable instrument/loader tests selected from current docs and code. Confirm an
IR holding marked `not-applied` or `unknown` is refused by the live telluric gate, and
that the recorded `normalization_state` agrees with the holding's declared
`pre_normalised` flag — `scripts/preflight_check.py` check 7 reports both.
