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

## Post-merge refresh (2026-09-13)

RYA-1176 is now merged in `main` and RYA-1134's verified pool is available there.
The refreshed eligibility matrix consumes the corrected manifest: zero engine or
holding cells carry `RYA1176_MISSING`. Gate 0 remains **CLOSED** because the
independent RYA-1141 scientific findings, exact holding/pixel validation, observed
conditioning records, and model/atom applicability are still unresolved. No
abundance run was authorized. The machine-readable checkpoint is
`gate0_post1176/status.json` with engine reason counts beside it.

## RYA-1155 coverage reconciliation (2026-09-13)

The next Gate 0 blocker is partially discharged. `pipeline.coverage` now recognizes
registered normalized spectra whose CSV is itself the manifest path, so the three
reachable CRIRES+ solar products resolve with explicit `csv_normalized` spans. Raw
archive inventories and the upstream VizieR delivery still remain intentionally
unaddressable. The Al intake builder fills legacy blank reach only from this registry;
the two census NIR intervals (13000–13195.23 A and 17493.69–19510.4 A) no longer fall
through to `OUTSIDE_CURRENT_INSTRUMENT_REACH`. Fresh RYA-1141 QA now reports `C-lines`
and `C-bands` PASS, while the overall measurement gate remains CLOSED on the other
scientific findings.

## RYA-1156 source-flag reconciliation (2026-09-13)

The Vujnovic CDS parser now preserves the documented limit and note flags (`l_e_Aki`,
`n_Aki`, `n_Lambda`, and related fields). Lower-limit uncertainties remain without a
determinate sigma, and the Al manifest carries an explicit `sigma_basis` for every
finite uncertainty. Fresh QA now reports `A1-flags` and `A5-sigma` PASS; the overall
measurement gate remains CLOSED.

The competing Vujnovic values are also retained in `competing_gf_summary` and the
conflict ledger without promotion. Fresh QA reports `A6` PASS.

The HFS/component reconciliation now carries source component counts and total log-gf
values in the manifest, and the canonical 3944.006/3961.520 rows carry counts 4 and 6.
Fresh QA reports `A3`, `A3-meta`, and `A3-rya1001` PASS.

The three misquoted DOI entries are corrected in the bibliography and follow-up ledger:
Griesmann & Kling (`10.1086/312741`), Nandakumar et al. (`10.3847/1538-4357/ad22dc`),
and Murphy & Berengut (`10.1093/mnras/stt2204`). Fresh QA reports `A5-doi` PASS.

## Established blockers

- The preserved RYA-1132 source manifest remains the audit baseline; its corrected
  RYA-1176 successor is now consumed by the downstream eligibility matrix. A risk
  label still cannot establish conditioning state, so product-time conditioning
  remains a hard gate.
- RYA-1134 is merged and its `verified_v2` pool is consumed by the refreshed
  matrix. It supplies atomic dispositions; it does not clear holding or product
  conditioning gates.
- The current executable RYA-1141 audit returns FAIL / measurement gate CLOSED: `C-lines`, `C-bands`, `A1-flags`, `A5-sigma`, `A6`, the HFS checks, and `A5-doi` now PASS; 50 findings remain (10 CRITICAL, 28 HIGH, 12 MEDIUM). Identity matching, evaluated-source semantics, raw holding policy, and model applicability remain unresolved. These require upstream scientific adjudication, not merely adding schema columns.
- RYA-1173 is merged and its AGSS21 census gate passes. The old claim that the Al reference census is entirely absent is superseded. This does not supply the missing RYA-1134 verified grades.

## Evidence and validation

`gate0_qa_control/` is the definitive fresh audit, including per-check and per-line finding CSVs. Reproduce from this main revision with:

```sh
python scripts/qa_al_intake_rya1141.py --check --out /tmp/rya1217-gate0-qa
python -m pytest tests/test_qa_al_intake_rya1141.py tests/test_al_intake_rya1132.py -q
```

The QA command exits zero despite its scientific FAIL verdict; inspect the JSON. The focused coverage/intake tests pass: **35 passed**. They validate the audit/intake software, not scientific measurement readiness.

The initial `gate0_qa/` run wrote into a new untracked repository directory, which the audit's whole-tree mutation check flagged. The control rerun wrote outside the repository and reports `artifacts_mutated: []`, with the same 85 scientific findings. Both runs are retained. The current executable battery has 53 checks; the historical ticket's 59-check total is not presented as a fresh result.

## Campaign disposition

- Verified pools used: none; upstream final adjudication unavailable.
- Holdings used: none. Registry inspection finds both corrected and uncorrected Solar holdings, with several unknown/blank normalization states; none is promoted into an exact-product conditioning assertion.
- Products / uncertainty audit / aggregate matrix: none, because Gate 0 fails before fitting.
- Per-line evidence: `gate0_qa_control/findings.csv` and its supporting CSVs are intake rejection evidence, not abundance measurements.
- Legacy products superseded: none operationally; all remain historical/non-authoritative for this restart until reproduced. No old values were copied or feed entries changed.
- Literature comparison: not performed; there is no newly frozen measurement to validate.
- Appendix recommendation: no new Al headline or forest product from this campaign yet.
- Next prerequisite: resolve the remaining applicable RYA-1141 child defects (identity, HFS, source provenance, and conditioning) before rerunning Gate 0.
- No source spectra, atomic data, feed, or current-state ledgers changed. Existing closed-gate state is confirmed, not newly signed off. No commit, PR, or merge was made.
