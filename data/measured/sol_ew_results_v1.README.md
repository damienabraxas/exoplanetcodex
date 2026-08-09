# sol_ew_results_v1.csv — the committed solar EW pool

**This file carries no comment header, deliberately.** It never had one, and twelve
consumers read it with a bare `pd.read_csv`. Adding a `#` block broke six of them
instantly (RYA-710). The documentation lives here instead.

## Instrument provenance — added RYA-710

Ryan, 2026-08-09: *"we need to filter which instrument produced which results."*

Every measurement now names the instrument that produced it. Before this column existed
**all 808 rows were HARPS and nothing said so** — so the file could not distinguish
*"we only ever used one instrument"* from *"this is what the sky offers"*. That ambiguity
is how a single-instrument pool stayed invisible while two atlases sat unused on disk,
and how a scope limit came to be reported as a data limit.

| column | meaning |
|---|---|
| `instrument` | an `instrument_id` from `data/catalog/instrument_catalog.csv`. **Required.** |
| `instrument_provenance` | how we know |

`instrument` is enforced by `pipeline.abundances_derive._assert_instrument_stamped`,
which raises on a missing column, a blank value, or an id absent from the catalog.
**Never defaulted** — defaulting a blank to `harps` would rebuild the same blind spot
with a friendlier face.

## What the 808 rows are

All **`harps`**, and verified rather than assumed: every row falls inside HARPS
3782.6–6910.0 Å, checked with `pipeline.coverage` against the catalog, zero exceptions.

- **804** trace directly to the June HARPS `lines_fit` output (RYA-396 vendoring).
- **4** — Na I 5688.205, Eu II 6645.127, Li I 6707.840, Na I 6154.225 — entered via the
  RYA-465 recovery and the RYA-102/103 charter force-includes. They are HARPS **by
  lineage**, not by an individual re-check, and are stamped `inherited` for that reason.
  A stamp that overstates its own evidence is worse than no stamp.

## Adding rows from another instrument

Register the instrument in the catalog first, then promote through
`scripts/promote_solar_ew.py`. Per-(instrument × band) reporting is ratified
(`docs/SCIENCE_STANDARDS.md`), so a mixed-instrument pool is expected — what is not
acceptable is a row that cannot say which arm produced it.
