# RYA-1217 Solar Al restart: Gate 0 CLOSED

## RYA-1176 implementation checkpoint

The frozen 505-row Al census now has a regenerated schema at
`data/audit/rya1176_al_manifest/al_line_manifest_v2.csv`. The source manifest is
preserved unchanged. Every row carries a canonical `line_set`, separate
`gf_provenance` and `selection_state`, and three independent conditioning axes.
Because conditioning is holding-dependent, those three fields are explicitly
`unknown` in the line census and are populated only by a registered holding at
product time. `pipeline.al_manifest.require_product_manifest()` refuses missing
holdings and unresolved conditioning; it never defaults unknown to safe.

Regenerate with `python scripts/rya1176_regenerate_al_manifest.py`. The output is
493 `our-all`, 5 `our-graded`, and 7 `our-deep-graded` rows, with provenance in
`data/audit/rya1176_al_manifest/provenance.json`. This fixes the manifest schema
gap while keeping Gate 0 closed until the corrected fields are propagated through
the live Al product path and the remaining upstream scientific gates pass.

Checked 2026-09-13 against fetched main `c64eccfe3f6a41462131ac0e66eb61c9323d8c11` on branch `codex/rya-1217-al-restart`.

No measurement was launched and no authoritative abundance was generated. The ticket explicitly requires the campaign to stay gated when prerequisite defects cannot be resolved by a tightly scoped run fix.

## Established blockers

- The current 505-row `data/audit/rya1132_al_intake/al_line_manifest.csv` has no `line_set`, `telluric_applied`, `normalization_state`, or `observed_conditioning` columns (RYA-1176). A risk label cannot establish conditioning state.
- RYA-1134 remains In Progress. Its branch contains `1e4dd8db` (reference-lineage setup), which is not an ancestor of main. Main contains no RYA-1134 grade-verification artifacts. Its comments explicitly leave the full 505-row adjudication unfinished. There is no final verified pool to consume.
- The current executable RYA-1141 audit returns FAIL / measurement gate CLOSED: 25 PASS, 16 FAIL, 12 FLAG; 85 findings (20 CRITICAL, 47 HIGH, 18 MEDIUM). Identity matching, HFS/component evidence, source flags, DOI provenance, evaluated-source semantics, and holding coverage remain unresolved. These require upstream scientific adjudication, not merely adding schema columns.
- RYA-1173 is merged and its AGSS21 census gate passes. The old claim that the Al reference census is entirely absent is superseded. This does not supply the missing RYA-1134 verified grades.

## Evidence and validation

`gate0_qa_control/` is the definitive fresh audit, including per-check and per-line finding CSVs. Reproduce from this main revision with:

```sh
python scripts/qa_al_intake_rya1141.py --check --out /tmp/rya1217-gate0-qa
python -m pytest tests/test_qa_al_intake_rya1141.py tests/test_al_intake_rya1132.py -q
```

The QA command exits zero despite its scientific FAIL verdict; inspect the JSON. The focused tests pass: **35 passed**. They validate the audit/intake software, not scientific measurement readiness.

The initial `gate0_qa/` run wrote into a new untracked repository directory, which the audit's whole-tree mutation check flagged. The control rerun wrote outside the repository and reports `artifacts_mutated: []`, with the same 85 scientific findings. Both runs are retained. The current executable battery has 53 checks; the historical ticket's 59-check total is not presented as a fresh result.

## Campaign disposition

- Verified pools used: none; upstream final adjudication unavailable.
- Holdings used: none. Registry inspection finds both corrected and uncorrected Solar holdings, with several unknown/blank normalization states; none is promoted into an exact-product conditioning assertion.
- Products / uncertainty audit / aggregate matrix: none, because Gate 0 fails before fitting.
- Per-line evidence: `gate0_qa_control/findings.csv` and its supporting CSVs are intake rejection evidence, not abundance measurements.
- Legacy products superseded: none operationally; all remain historical/non-authoritative for this restart until reproduced. No old values were copied or feed entries changed.
- Literature comparison: not performed; there is no newly frozen measurement to validate.
- Appendix recommendation: no new Al headline or forest product from this campaign yet.
- Next prerequisite: finish RYA-1134's verified-grade dispositions and resolve the applicable RYA-1141 child defects, including RYA-1176, before rerunning Gate 0.
- No source spectra, atomic data, feed, or current-state ledgers changed. Existing closed-gate state is confirmed, not newly signed off. No commit, PR, or merge was made.
