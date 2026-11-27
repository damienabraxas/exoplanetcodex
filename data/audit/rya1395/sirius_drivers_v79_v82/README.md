# RYA-1395 v79–v82: Sirius-only drivers, committed verbatim

Codex ran these on Sirius between 11-25 and 11-26. They were never committed. The files here are byte copies. Each path is relative to `/mnt/codex-data/codex/rya1354` on Sirius.

| Run | What | Result (Linear RYA-1395) |
|---|---|---|
| v79 | n6 E1 only (`rtransition`, 24 h budget) from the v78 biotransform inputs | Completed in 6,592 s. The n6 strengths are anomalous (777 nm Babushkin +1.8 to +2.3 dex vs n5). HOLD |
| v81 | Private `rbiotransform` build passing the full `EVEC` instead of `EVEC(1:500)` | n5 is byte-identical. n6 is byte-identical to v79, so the slice hypothesis is refuted |
| v82 | Matrix only: f-block overlap order permutations | Joint reversal of 6f,5f,4f cuts the max LU multiplier from 2.8e5 to 1.1 |

The outputs, raw tables and receipts stay on Sirius. The successor is v83: `scripts/rya1395_order_control_v83.py`.
