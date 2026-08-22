# tau Ceti CRIRES+ — RESOLVED by RYA-957 (this file records the sequence)

**STATUS: UNBLOCKED AND REGISTERED.** `tau_cet_crires_plus` is now a row in
`data/catalog/holdings_manifest_registry.csv` with `system_id = tau_ceti`, and
`instruments.validate_all()` returns no errors.

This file is kept rather than deleted because the sequence is the useful part.

## What was blocked (RYA-952)

`instruments.validate_holdings` resolves a holding's `system_id` against the set of
**non-blank `star_params_key`** values in `data/catalog/system_catalog.csv`, and
`pipeline/system_catalog.py` requires any non-blank key to **resolve in `STAR_PARAMS`**.
tau Ceti sat in the catalogue as a `future_target` with a blank key, so registering the
holding meant first adopting real stellar parameters — a science decision that RYA-952
(data prep only) had no business making. **The row was written, the validator refused it,
and it was backed out rather than satisfied by inventing parameters.**

## What unblocked it (RYA-957)

Ryan ratified applying the project's existing GBS standard rather than inventing anything.
tau Ceti (HD 10700) and eps Eri (HD 22049) are both Gaia FGK Benchmark Stars, so their
parameters are **published on the scale `config/stars.yaml` already declares**:
Heiter+2015 (Teff/log g, Paper I) and Jofré+2014 ([Fe/H], xi, Paper III). Adopting them
APPLIES the standard; it does not invent a value. Both records were added, both
`star_params_key` pointers set, and this holding re-registered.

🔴 **One value did NOT come across as a GBS fundamental.** Heiter+2015 Table 10 prints tau
Ceti's log g in **square brackets**, which that table's own caption defines as "uncertain
and should not be used as a reference for calibration or validation purposes"; Sect. 5.4.1
gives the reason — the mass comes from "problematic evolutionary tracks predicting
unreasonable ages". So tau Ceti carries log g as a starting value and **SOLVES** it, while
eps Eri (unbracketed) pins it. Pinning a value the source paper tells you not to use as a
reference would have been the same defect in a new place.

## Still open — the CLASS, not this instance

The registry still **cannot express "we hold data for a star we have no parameters for."**
RYA-957 fixed the two stars that had a published standard to apply; it did not change the
schema. **tau Boo remains blocked and is NOT a Gaia FGK Benchmark Star** — verified against
both papers: it appears in Heiter+2015 only in the *candidate* discussion (Sect. 7) and a
candidate appendix table, never in Table 10, and it is absent from Jofré+2014 entirely. It
needs either best-literature parameters from a separate cited source, or the FK-relaxation
path. Both are separate decisions.
