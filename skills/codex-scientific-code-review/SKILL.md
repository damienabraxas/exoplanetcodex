---
name: codex-scientific-code-review
description: Review scientific pipeline code for wrong-data risk, provenance loss, cross-run contamination, silent fallbacks, unit/frame/model mistakes, and superficial tests. Use for scientific code review, shared multi-instrument paths, or when a run can execute while using the wrong holding. Does not ratify new science rules or modify code unless explicitly requested.
---

# Scientific-software QA review

## Authority and scope

Start with `LEDGERS.md` and applicable science, protocol, convention, instrument, and
model sources. Review current committed behavior, not legacy prose. This skill finds and
classifies defects; it does not invent scientific values, ratify constraints, or treat a
passing execution as scientific validity.

## Review priority

Review in this order: scientific correctness/wrong-data risk; provenance and cross-run
contamination; silent failures and invalid fallbacks; units, frames, and applicability;
test coverage/reproducibility; performance/complexity; maintainability/style. Classify
findings `CRITICAL`, `WARNING`, or `SUGGESTION`, cite evidence, and route real unresolved
defects to durable Linear work rather than leaving them only in prose.

## Configuration and provenance

- Shared logic does not hardcode instrument, star, wavelength, path, model, calibration,
  or mutable threshold values when a canonical config/catalog/registry owns them.
- A hardcoded value is acceptable only for an immutable domain constant or explicit
  scientific identifier with clear ownership/citation.
- Instrument behavior is keyed by explicit instrument **and product/holding** identity,
  never by the first dataset run, filename, stale module state, or inferred label.
- Every output carries enough source identity to trace input holding, line/accounting
  artifact, model/atom/grid, units/frame, and generator commit.
- Any cache key contains every scientifically relevant identity axis (at least star,
  instrument/holding, band/window/line set, route/model/atmosphere/gf, and config vintage
  where those affect the value). Prefer immutable local state; flag module-global run state.

## Cross-instrument anti-contamination proof

For each shared multi-instrument path, trace requested identity → loader selection →
measured arrays/lines → per-line record → aggregate/product metadata. Require a canary
that fails when those disagree. Where practical, run two deliberately distinguishable
fixtures (for example Kitt Peak and HARPS) and assert that holding IDs, input arrays or
line artifacts, and emitted provenance do not collapse to one stale product. A run for
instrument B must fail if any output still names or consumes instrument A's holding,
wavelength/line-accounting artifact, cache entry, or measurement.

Static bans on hardcoded loaders are useful but insufficient: also require a behavioral
two-fixture test and end-to-end output identity. “The script ran” and “a number exists”
are not evidence.

## Failure behavior and code quality

- No broad `except:` or swallowed exceptions. Explicit recoverable scientific absence
  gets a disposition; unexpected/invalid states fail loudly with context.
- Units and wavelength/velocity frames are explicit at boundaries and in names/types
  where useful.
- Functions are single-purpose with explicit data flow and limited side effects.
- Prefer vectorized/indexed NumPy or scientific-library operations for wavelength arrays;
  small bounded control-flow loops remain acceptable when clearer and correct.
- Avoid copy/paste instrument branches when a tested abstraction/config route fits, but
  do not erase genuine instrument differences behind a false universal abstraction.
- Comments explain scientific intent and invariants. Public/complex functions document
  inputs, outputs, units, frames, failures, and provenance expectations.

## Report

For every finding include severity, affected path/line, failure mechanism, scientific
impact, existing or missing test, and the smallest safe remediation. Explicitly state
which controls were checked and whether they discriminate. A review passes only when no
scientifically material wrong-data or provenance ambiguity remains.
