# Completed Reference ξ results integrated — RYA-587

The RYA-1213 shell campaign already performed the work. This update consumes its
existing outputs; it launches no synthesis and changes no abundance measurement.

- Canonical campaign: 64 applicable pools, 60 measured, four two-line pools below
  the existing `min_paired=3` rule. The latter have completed runs; they are not
  jobs to repeat merely to obtain the same two lines.
- Reference feed: 60 measured, four small-N held, four full-3D N/A. No Reference
  record remains `NOT_IN_CAMPAIGN`.
- Exactly 13 VIS product stamps change. The other 147 records are byte-identical.
  All 160 abundances, line counts, original statistical/gf-systematic values and
  measurement provenance remain unchanged.
- Corrected totals include statistics, original systematics and ξ exactly once.
  Stale missing-ξ caveats are cleared. The plot grid is rebuilt from the feed.
- `xi_term_owed.json` now names the current result and preserves its original
  all-missing assertion under `historical_snapshot`.

Reproduce from the existing copied shell outputs:

```bash
python3 scripts/rya1213_xi_campaign.py --runs /path/to/rya1213_xi
python3 scripts/rya1178_emit_fe_schema.py --xi-only --check
python3 scripts/rya1178_emit_fe_schema.py --xi-only
python3 scripts/audit_fe_existing_evidence_rya587.py
```

The new ingestion regression failed before restamping and passes after it. It
requires every measured canonical Reference pool to reach the matching product
with the correct derivative, ξ term and total. A second test checks idempotence
and preservation of measurement fields.

## Existing Fe evidence, distinct from the new schema

All 160 live Fe products have their original measurement artifact locally.
137 also have a directly adjacent budget text file; the other 23 need their
particular artifact layout resolved, not a new abundance measurement. Existing
budget terms and separate ξ sources are indexed in `fe_existing_evidence_rya587.json`.
A legacy budget text's UNMEASURED stellar label may predate its later ξ stamp.

Outside Reference, older Codex/Deep records retain 16 NOT_IN_CAMPAIGN, 17 UNMEASURED
(including four experimental mean-3D records), and six ALIASED ξ labels. They are
not changed by this 13-result integration. The old campaign's six ALIASED records
have completed derivatives but different recorded pool sizes (9 vs 8, or 3 vs 2).
Reconcile exact pool evidence before changing those labels. This audit does not
turn missing schema fields or old labels into a request for another Fe run.

The new RYA-587 JSON contract is an integration task over these existing records.
It is not evidence that the Fe measurements, statistical errors, gf terms or ξ
campaign were never performed. Scientific method questions (such as small-N
admission or MLP accuracy) are separate from missing jobs.
