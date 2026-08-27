# RYA-1084 — what rewrote the ten RYA-1080-blocked feed products

RYA-1080 refused to commit ten Fe feed products whose source files no longer matched the
sha256 the feed recorded, and held them rather than re-ingesting blind. This is what did it.

## Verdict: a legitimate, deterministic generator. The FEED was the stale side.

| question | answer |
|---|---|
| what writes them | `scripts/derive_band_products.py` |
| is it deterministic | **yes** — see below |
| is it the RYA-1070 suite-rewrite hazard | **no** |
| is the `sigma_stat` move a `round()`/`np.round()` tie | **no** — they agree at this value |
| resolution | writer cleared → **all ten re-ingested**, 75 of 75 reconciled |

## 1. The writer

The ten rewrites fall into **five base + `ENGINE-A` pairs, 20–40 s apart** — the shape of a
derivation run, not a stray touch or an editor save:

| cluster (UTC) | product |
|---|---|
| 2026-08-26 16:38 / 16:39 | red-optical, molecfit, GRADED |
| 2026-08-26 16:43 / 16:44 | near-UV, molecfit, DEEPGRADED |
| 2026-08-26 17:04 ×2 | red-optical, kurucz2005, GRADED |
| 2026-08-26 17:08 / 17:09 | near-UV, kurucz2005, DEEPGRADED |
| 2026-08-27 13:34 ×2 | VIS, molecfit, GRADED |

**The RYA-1070 hypothesis is refuted directly.** On clean `main` (`6c1529f`), snapshotting
all 80 `band_products` files, running the full suite, and re-snapshotting: **0 files
changed**, and `git status` was empty. The suite does not rewrite band products.

## 2. It is deterministic — and that is what clears it

🔴 **The 2026-08-27 13:34 run reproduced, byte for byte, the file RYA-1051 had already
committed** (`eb35543`, "pin the current solar Fe I values AND the process"). Both the repo
copy and the regenerated Mac copy hash to `1eb939416be7`. A nondeterministic writer does not
land on the same bytes.

So the 0.0217 → 0.0218 move is attributable to the code changing between the feed's ingest
(2026-08-24 01:11) and RYA-1051 — the same window in which **RYA-1044/1045 added the `deck`
column**, which is exactly the diff six of the ten carry:

```
committed: …,route_basis
mac      : …,route_basis,deck        ← one new column, every value identical
```

## 3. The one value that moved, classified

Fe I VIS 1D-LTE GRADED, `sigma_stat` **0.0217 → 0.0218**. Not a rule disagreement:

| value | `repr` | `round(x,4)` | `np.round(x,4)` | `round_dex(x)` |
|---|---|---|---|---|
| `0.02175` | `0.02175` | **0.0217** | **0.0217** | 0.0218 |
| `nextafter(0.02175, 1)` | `0.021750000000000002` | **0.0218** | **0.0218** | 0.0218 |

**`round()` and `np.round()` agree at both points.** The split is a **one-ULP difference in
the upstream RMS**, amplified by a decimal boundary: `0.02175` as a float is
`0.0217499999999999985`, one ULP *below* the tie, so it rounds down while its immediate
neighbour rounds up.

It looked like a tie because `stat_basis` is formatted `{m:.5f}` and prints **`0.02175` for
both** — the artifact carries no evidence of which float produced it.

⚠️ `round()` and `np.round()` *do* disagree elsewhere (`0.12345` → 0.1235 vs 0.1234;
`5e-05` → 0.0001 vs 0.0), which is why the rule needed pinning at all. They just are not the
explanation here, and recording them as the cause would have sent the next person after the
wrong defect.

## 4. The determinism fix

`pipeline.error_budget.round_dex` — **one rule, single-sourced**: half-up on the shortest
round-tripping decimal (`repr`). It is stable across exactly the ULP pair that `round()`
splits, and measured against `round()` on 200,000 random values it differs on **zero** of
them: it changes nothing but the tie behaviour, which was the only ambiguous part.

The five scattered `round(x, 4)` calls in `derive_band_products.py` now route through it. A
rule declared five times is a rule nobody owns — the RYA-845 lesson.

⚠️ **What it does not do.** It does not make a genuinely different number the same. Values
straddling a boundary by more than the tie neighbourhood still round apart, correctly.
Upstream 1-ULP nondeterminism would be a separate exposure; nothing here suggests one, but a
rounding rule could not stand in for fixing it if there were.

## 5. Resolution

All ten re-ingested (`scripts/rya1084_reingest_blocked.py`), each row carrying
`reingested_by: RYA-1084` and the reason. The reconciliation now reports **75 of 75
RECONCILED, zero blocked**, and RYA-1080's guard exits 0 with 67 recorded regenerability
gaps. RYA-1080's strict xfail is removed — the ten are resolved, not tolerated.

`data/results/rya1080/rya1080_blocked.csv` is **kept**: a now-empty guard output no longer
shows what was wrong or by how much, and that evidence is the point.
