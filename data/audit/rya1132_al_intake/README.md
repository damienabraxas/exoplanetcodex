# RYA-1132 - Al atomic-data intake closure

This is an inventory and provenance freeze, not an abundance result. It preserves
**505** physical candidates (466 Al I,
39 Al II) from UV through IR. No canonical gf
or Solar abundance was changed.

## Verdict

- UV: `PARTIAL_GF_LAB_INGESTED_POLICY_BLOCKED` - all 106 Vujnovic CDS rows are
  preserved and normalized. Six finite, uncertainty-bearing Al I transitions plus
  Johnson et al. 1986 Al II 2669 are physically crossmatched and ingested as GF-LAB.
  Limits, ratio-only rows, and 3092.839 without an independent uncertainty remain
  non-promoted; FUV/NUV still lack a declared measurement policy.
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

The ingestion changes seven manifest rows from fallback to GF-LAB without mutating
`canonical_gf` or deriving an abundance. `vujnovic2002_normalized.csv` and
`vujnovic2002_crossmatch.csv` retain the complete source accounting.

The web follow-up found no new independent IR source beyond the older Buurman,
Davidson, and Buurman-Donszelmann measurements already propagated through
Burheim/NIST. In particular it did not resolve 11254.9, 7835/7836, 8772/8773,
or the 21208 A IGRINS candidate. See `web_source_followup.csv`.

## Reproduce

`python3 scripts/build_al_intake_rya1132.py`
