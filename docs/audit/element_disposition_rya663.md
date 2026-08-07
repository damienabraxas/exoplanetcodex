# Per-element disposition — pre-527 straggler sweep (RYA-663)

**GENERATED — do not hand-edit.** Regenerate with `python scripts/gen_element_disposition.py`.

Live channel: PASS **6** · CURATION-OWED **20** over 26 elements, phase_c generated 2026-07-27.

Gates (ratified RYA-561, applied via `engine_selection.evaluate_floor_promotion`): tolerance **0.1** dex, cross-engine **0.1** dex. Gate 3 is STRICT — a missing delta fails.

## Can flip to PASS now

**Ca**  — and this rests on a stale input; see below.

## Inputs read

| artifact | commit | committed |
|---|---|---|
| `data/audit/cno_synthesis/solar_phase_c_verdict.json` | `5c1e92d` | 2026-07-27T00:13:37-06:00 |
| `data/audit/rya527_two_engine/solar_two_engine_records.json` | `ae518e8` | 2026-07-27T00:50:37-06:00 |
| `data/reference/solar/solar_abundances_v1.csv` | `583cb46` | 2026-06-28T22:43:33-06:00 |
| `data/reference/solar/solar_abundances_v2.csv` | `a9dd958` | 2026-07-05T01:38:42-06:00 |

## ⚠ Stale-input evidence

The two-engine record disagrees with the live verdict channel on these elements, so it demonstrably predates them. **Every gate-3 number in this report is read from that artifact.**

- Cr: two-engine reports 5.654 but the live verdict channel measures 6.022 (delta -0.368) — the two-engine artifact predates that measurement
- Fe: two-engine reports 7.580 but the live verdict channel measures 7.466 (delta +0.114) — the two-engine artifact predates that measurement
- Li: two-engine reports 1.409 but the live verdict channel measures 0.727 (delta +0.682) — the two-engine artifact predates that measurement
- O: two-engine reports 8.730 but the live verdict channel measures 8.735 (delta -0.005) — the two-engine artifact predates that measurement
- S: two-engine reports 7.369 but the live verdict channel measures 7.486 (delta -0.117) — the two-engine artifact predates that measurement
- Si: two-engine reports 7.639 but the live verdict channel measures 7.888 (delta -0.249) — the two-engine artifact predates that measurement

## Value disagreements across artifacts (the cleanup list)

Invisible to the RYA-632 ledger guard, which compares verdicts and counts but never values. Each must be explained or corrected before the v3 freeze.

| element | spread | values |
|---|---|---|
| Cr | 0.368 | phase_c (live) **6.022**; gold (frozen) **6.022**; two-engine (review) **5.654** |
| Fe | 0.114 | phase_c (live) **7.466**; gold (frozen) **7.516**; two-engine (review) **7.580** |
| Li | 0.682 | phase_c (live) **0.727**; gold (frozen) **0.727**; two-engine (review) **1.409** |
| Mn | 0.004 | phase_c (live) **5.466**; gold (frozen) **5.470**; two-engine (review) **5.466** |
| O | 0.005 | phase_c (live) **8.735**; gold (frozen) **8.735**; two-engine (review) **8.730** |
| S | 0.117 | phase_c (live) **7.486**; two-engine (review) **7.369** |
| Si | 0.249 | phase_c (live) **7.888**; gold (frozen) **7.888**; two-engine (review) **7.639** |

## Per-element

| El | verdict | bucket | A(X) | ref | dCE | g1 | g2 | g3 | blocker |
|---|---|---|---|---|---|---|---|---|---|
| Al | CURATION-OWED | owed-HELD | 7.406 | 6.430 | — | ✗ | ✗ | UNEVALUABLE | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for Al); gate 2: |A - ref| = 0.976 > 0.1; gate 3: UNEVALUABLE — no two-engine record (the RYA-527 re-run is what produces one) |
| Ba | CURATION-OWED | measured-awaiting-freeze | 2.410 | 2.270 | — | ✓ | ✗ | UNEVALUABLE | gate 2: |A - ref| = 0.140 > 0.1; gate 3: UNEVALUABLE — no two-engine record (the RYA-527 re-run is what produces one) |
| C | PASS | PASS | 8.491 | 8.460 | — | — | — | n/a | already PASS — nothing owed |
| Ca | CURATION-OWED | owed-HELD | 6.324 | 6.300 | 0.016 | ✓ | ✓ | OK | none — promotes under the ratified three gates |
| Co | PASS | PASS | 4.965 | 4.940 | — | — | — | n/a | already PASS — nothing owed |
| Cr | CURATION-OWED | EW-pool | 6.022 | 5.620 | -1.058 | ✗ | ✗ | FAILED | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for Cr); gate 2: |A - ref| = 0.402 > 0.1; gate 3: |dCE| = 1.058 > 0.1 |
| Cu | CURATION-OWED | measured-awaiting-freeze | 4.345 | 4.180 | — | ✗ | ✗ | FAILED | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for Cu); gate 2: |A - ref| = 0.165 > 0.1; gate 3: single-engine record — zero independent confirmation |
| Eu | CURATION-OWED | owed-BLANK | — | 0.520 | — | ✗ | ✗ | UNEVALUABLE | no value exists (zero graded survivors) — needs a real line, not a decision |
| Fe | PASS | PASS | 7.466 | 7.460 | 0.084 | — | — | n/a | already PASS — nothing owed |
| K | PASS | PASS | 5.099 | 5.070 | — | — | — | n/a | already PASS — nothing owed |
| Li | CURATION-OWED | upper-limit | 0.727 | 1.050 | 0.682 | — | — | n/a | ratified UPPER_LIMIT disposition (RYA-563) — structurally never a PASS point value. Excluded from the flip denominator. |
| Mg | CURATION-OWED | owed-BLANK | 7.614 | 7.550 | — | ✓ | ✓ | FAILED | EW pool is blank (0 graded survivors) but a value exists off another channel [two-engine reported (REVIEW artifact)] — reconcile the channels before promoting; gate 3: single-engine record — zero independent confirmation |
| Mn | PASS | PASS | 5.466 | 5.420 | — | — | — | n/a | already PASS — nothing owed |
| N | CURATION-OWED | measured-awaiting-freeze | 8.188 | 7.830 | — | ✗ | ✗ | UNEVALUABLE | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for N); gate 2: |A - ref| = 0.358 > 0.1; gate 3: UNEVALUABLE — no two-engine record (the RYA-527 re-run is what produces one) |
| Na | CURATION-OWED | owed-HELD | 6.264 | 6.240 | -0.121 | ✗ | ✓ | FAILED | gate 1: Engine-B atom not RYA-534-validated (Na_gerber2023.prov.json: no verdict recorded); gate 3: |dCE| = 0.121 > 0.1 |
| Ni | CURATION-OWED | owed-HELD | 6.946 | 6.200 | -1.253 | ✓ | ✗ | FAILED | gate 2: |A - ref| = 0.746 > 0.1; gate 3: |dCE| = 1.253 > 0.1 |
| O | PASS | PASS | 8.735 | 8.690 | — | — | — | n/a | already PASS — nothing owed |
| P | CURATION-OWED | measured-awaiting-freeze | 6.610 | 5.410 | — | ✗ | ✗ | UNEVALUABLE | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for P); gate 2: |A - ref| = 1.200 > 0.1; gate 3: UNEVALUABLE — no two-engine record (the RYA-527 re-run is what produces one) |
| S | CURATION-OWED | measured-awaiting-freeze | 7.486 | 7.120 | 0.041 | ✗ | ✗ | OK | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for S); gate 2: |A - ref| = 0.366 > 0.1 |
| Sc | CURATION-OWED | measured-awaiting-freeze | 3.203 | 3.140 | — | ✗ | ✓ | UNEVALUABLE | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for Sc); gate 3: UNEVALUABLE — no two-engine record (the RYA-527 re-run is what produces one) |
| Si | CURATION-OWED | EW-pool | 7.888 | 7.510 | -0.756 | ✓ | ✗ | FAILED | gate 2: |A - ref| = 0.378 > 0.1; gate 3: |dCE| = 0.756 > 0.1 |
| Sr | CURATION-OWED | owed-HELD | — | 2.830 | — | ✓ | ✗ | FAILED | gate 2: no value to test; gate 3: single-engine record — zero independent confirmation |
| Ti | CURATION-OWED | owed-HELD | 5.471 | 4.970 | -1.326 | ✗ | ✗ | FAILED | gate 1: Engine-B atom not RYA-534-validated (Ti_gerber2023.prov.json: CHECK (not PASS): |median-anchor|=0.114 > tol 0.06. FINDING — NOT a deck/line failure: departures engaged, 3 directly-MPIA-comparable well-measured lines all ~2x the Berg… *(full text in element_disposition_rya663.json)* |
| V | CURATION-OWED | measured-awaiting-freeze | 3.917 | 3.900 | — | ✗ | ✓ | FAILED | gate 1: Engine-B atom not RYA-534-validated (no RYA-534 Engine-B grid provenance on record for V); gate 3: single-engine record — zero independent confirmation |
| Y | CURATION-OWED | owed-BLANK | — | 2.210 | — | ✗ | ✗ | UNEVALUABLE | no value exists (zero graded survivors) — needs a real line, not a decision |
| Zr | CURATION-OWED | owed-BLANK | — | 2.590 | — | ✗ | ✗ | UNEVALUABLE | no value exists (zero graded survivors) — needs a real line, not a decision |

## What each bucket means

- **PASS** (C, Co, Fe, K, Mn, O) — already PASS — nothing owed
- **owed-HELD** (Al, Ca, Na, Ni, Sr, Ti) — pool survived and a value exists in gold v1, held unfrozen by the ratified RYA-522 tier. This is a PROMOTION DECISION: run the three gates.
- **owed-BLANK** (Eu, Mg, Y, Zr) — the graded cull left ZERO survivors — there is no value to promote. Needs line-pool / gf work (a real second line), not a ratification.
- **measured-awaiting-freeze** (Ba, Cu, N, P, S, Sc, V) — measured on a dedicated synthesis / Kitt Peak channel; the value exists but is not frozen. Clears at the RYA-527 v3 re-freeze.
- **EW-pool** (Cr, Si) — a plain EW pool sitting at the gf floor — a TIER question (owed vs gf_floor), adjudicated against Ti/Cr/Si, not a measurement gap.
- **upper-limit** (Li) — ratified UPPER_LIMIT disposition (RYA-563) — structurally never a PASS point value. Excluded from the flip denominator.
