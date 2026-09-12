# `solar_iag` NIR Fe I, re-measured on the guarded pool — RYA-1213

**MEASURED AND HELD. Not in the feed, and deliberately so.**

The live Codex Grade row for this cell (`A(Fe I) = 7.546`, `n = 25`, artifact
2026-08-25) predates RYA-1191's `pipeline.fit_validity` guard, and its aggregate
contains four non-convergent fits:

| λ (Å) | A(Fe) fitted | dex from solar | guard |
| --- | --- | --- | --- |
| 9350.415 | 5.708 | 1.75 | FIT-NOT-PHYSICAL |
| 9410.159 | 5.530 | 1.93 | FIT-NOT-PHYSICAL |
| 9437.793 | 4.713 | 2.75 | NON-MINIMUM |
| 9454.194 | 4.876 | 2.58 | FIT-NOT-PHYSICAL |

9437.793 is the very line RYA-1191 named as the bad fit behind the refuted
1.315 dex Kitt Peak dispersion. It is still inside this aggregate.

## How it surfaced

NIR holds no primary-lab Fe line above the 0.60 depth gate, so RYA-1213's
Reference pool and the Codex pool are **the same lines**. The Reference run
therefore re-measured this cell on current code, and the two disagreed on `n`
— a pool difference, not drift, visible only because the pools are provably
identical.

## What is in this directory

The GRADED selector re-run on commit `7ca66630`, differing from the live row by
the code vintage and nothing else:

| treatment | A | n | live row |
| --- | --- | --- | --- |
| 1D-LTE | 7.554 | 22 | 7.546, n=25 |
| ENGINE-A | 7.599 | 6 | 7.599, n=6 |
| synth-1D-LTE-gerber | 7.549 | 22 | 7.549, n=22 |

The same run is also the paired selector control: it reproduces the Reference
product to every published digit on all three treatments, so the selector
contributes **zero** and the whole 0.008 dex is the guard.

## Why it is not published

RYA-1213 owns the Reference tier, not the Codex one, and its directive routes
findings to the EOS. Replacing a published value is Ryan's call — the RYA-853
precedent, where 37 corrections with complete evidence were still held for him.

Committed rather than left on Sirius because a result that outlives its artifact
is the RYA-1034 failure exactly. Ingesting it is one command:

    python3 scripts/publish_product.py \
      --from data/results/rya1213_iag_nir_regraded/*_products.csv \
      --holding solar_iag --tier GRADED --selector GRADED --route SYNTH \
      --host Sirius --origin-path /home/damienabraxas/scratch/rya1213/data/results/band_products/ \
      --reason "RYA-1213: the live row predates RYA-1191's fit_validity guard and its aggregate carries four non-convergent fits; re-measured on the guarded pool."
