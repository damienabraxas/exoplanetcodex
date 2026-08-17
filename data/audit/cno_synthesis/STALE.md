# ⚠️ THE CNO PRODUCT SET IN THIS DIRECTORY IS STALE

**Flagged 2026-08-17 (RYA-848 item 3). Do not quote these numbers as current.**

The artifacts here were last generated **2026-06-27** (`1ddde48`, RYA-371 Phase A).
Since then `data/linelists/canonical_gf.csv` alone has taken **11 commits** — RYA-354,
RYA-492, RYA-592, RYA-822, RYA-834, RYA-837 among them — and the CNO synthesis has not
been re-run against any of it.

## This was measured, not inferred

Running the **unmodified** CNO code on today's inputs reproduces neither the values nor
the χ² of the banked run:

| | banked (2026-06-27) | unmodified code, today |
|---|---|---|
| A(N), CN_red | 7.558 | **7.456** (−0.102 dex) |
| A(O), OI_6300 | 8.800 | 8.806 |
| A(C), C2_Swan | 8.437 | 8.453 |
| **C/O** | **0.491** | **0.484** |
| red_χ², CH_Gband | 18.412 | 18.546 |

The last row is the cleanest tell: `A_X` is identical at 8.491 while `red_χ²` moves, so
the χ² surface itself changed. That is input drift, not a code defect.

## What this is NOT

It is **not** a defect introduced by RYA-848. That ticket changed only the reported
uncertainty, and proved it: a same-inputs control run holds every `A_X`, every `red_χ²`,
every `σ_sys` and C/O identical between the old and new code. The drift above is present
with the σ fix reverted.

⚠️ **The general lesson, which cost a wrong first read here:** comparing a fresh run
against a *banked artifact* measures code change **and** input drift together. Only a
control run on identical inputs isolates the code.

## What is owed

A **full re-measurement** of the CNO product set — a re-run, not a σ patch. Deliberately
deferred: the roadmap is Fe → Al → … → CNO, and re-measuring now would burn the compute
before the machinery that should govern it exists.

When it happens it must go through **RYA-847's synthesis constraint gate**, which is also
what retires the last known defect in this path: `curvature_sigma`'s result is still
passed through `np.clip(..., 0.0, 1.0)`, so a fit whose abundance the data does not
constrain reports `σ = 1.000` — a sentinel shaped like a measured 1 dex. With RYA-848's
rescale in place `CN_red` lands **on** that clip (true upper-side value ≥ 1.675 dex), so
solar N here would read `σ_stat = 1.000`. That number is left standing in these stale
artifacts precisely because they are not being re-published; it goes when CNO is genuinely
re-measured through the gate.
