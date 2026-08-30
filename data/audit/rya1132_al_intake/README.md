# RYA-1132 - Al atomic-data intake closure

This is an inventory and provenance freeze, not an abundance result. It preserves
**505** physical candidates (466 Al I,
39 Al II) from UV through IR. No canonical gf
or Solar abundance was changed.

## Verdict

- UV: `BLOCKED_ATOMIC_DATA` - the canonical UV Al inventory has no primary-lab
  gf, and FUV/NUV have no declared measurement policy.
- VIS: `FROZEN_WITH_DOCUMENTED_FALLBACKS` - physical identities and evidence
  ceilings are explicit; Burheim laboratory lines are not blurred with fallback gf.
- IR: `BLOCKED_PIPELINE_COVERAGE` - RYA-1003 telluric verification and RYA-1004
  synthesis red-edge/context coverage remain open.
- Overall: `BLOCKED_PIPELINE_COVERAGE`; no new Solar Al measurement is unblocked.

The 6696.015 Burheim transition remains physically distinct from 6696.185. The
11254.9 conflict ledger explicitly distinguishes Burheim's strong-component
`log gf=+0.327` from the observed blended-feature total near `+0.354`. IGRINS and
CRIRES+ wavelength-only evidence remains HOLD unless wavelength plus EP/levels
establish a unique physical transition.

## Reproduce

`python3 scripts/build_al_intake_rya1132.py`
