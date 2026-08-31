# RYA-1141 - independent QA of the RYA-1132 Al intake

Findings only. Nothing here changes a gf, a grade, a manifest row or
`canonical_gf.csv` - RYA-161's validate-don't-tune firewall. Every defect below is
REPORTED and filed as a child ticket; none is fixed in this directory.

**Overall verdict: `FAIL`** - see `verdict.md` for the per-check table and
`findings.csv` for every defect with its offending rows named.

## What refereed what

Agreement is only evidence when the referee is independent of the thing it judges,
so each check re-derives rather than re-reads:

| claim | independent referee |
| --- | --- |
| 106 Vujnovic rows, fully transcribed | the CDS `ReadMe`'s own declared record counts |
| the fixed-width column extraction | branching closure, `A_i / sum(A)` vs the separately printed `BranR` |
| no air/vacuum or cm-1/Angstrom conflation | the Al I 3p ^2^P^o^ splitting, 112.061 cm-1, recovered from vacuum wavelengths |
| the promotion join is EP-aware | an AST test for a physical-identity COMPARISON, with `nearest` as its positive control |
| every DOI names the paper it claims | the registered Crossref AUTHOR LIST, not the volume |
| `canonical_gf` was not mutated | the file list of PR #478's merge diff |

## Files

`verdict.md` / `verdict.json` - the per-check PASS/FAIL/FLAG table.
`findings.csv` - every finding, severity, subject, evidence.
`check_results.csv` - the full reasoning behind each check.
`a1_*`, `a2_*`, `a3_*`, `a4_*`, `a5_*`, `a6_*`, `b1_*`, `b3_*`, `c_*` - the
per-check evidence ledgers named in the verdict.

## Reproduce

`python3 scripts/qa_al_intake_rya1141.py --check`

DOI resolution runs from a committed Crossref cache so the battery is deterministic
offline; `--online` re-resolves every identifier live and reports any drift.
