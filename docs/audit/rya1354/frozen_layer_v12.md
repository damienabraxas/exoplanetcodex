# Frozen-prior-layer sensitivity and terminal full-release failure

The n3-adapted n4 full-release campaign remains **FAILED_RETAINED**. All three
outer-layer optimizations and their fresh restarts passed independently, but
released ground failed at cycle 379: both final 2s methods reached 24/24 radial
iterations with three nodes instead of one. Last complete cycle 378 had maximum
residual 29.9. No released state or fully relaxed lifetime is admitted.

The complete terminal archive preserves 281 native paths in one 2,333,792-byte
part. Receipt SHA256:
`4ebdf6362049d42c3a2bc1f2a0fabb0bcb7d4b02e8b8e49922936fab80341c0e`.
The independently replayed six restricted stages have these maximum residuals:

| State | Optimization cycles / residual | Fresh cycles / residual |
|---|---|---|
| Ground | 32 / 5.11e-8 | 1 / 1.77e-8 |
| Quintet | 33 / 5.45e-8 | 1 / 2.21e-8 |
| Triplet | 35 / 4.30e-8 | 1 / 1.79e-8 |

All pass original ACCY 5.960464477539063e-8, H=0.03125, N=439, physical
labels, exact target configuration/nucleus, source/restart identity and all nine
frozen n<=3 records against same-input native loading baselines. Only n4 varied.
This does not establish full core relaxation or remove the historical n3
spectroscopic-component constraint limitation.

## Separately declared model

Before native rate execution, Linear comment
`12ea7ea9-0416-4d5d-a53e-5002fcd711a9` declared a separate frozen-prior-layer
model sensitivity. Plan `config/atomic/rya1354/lifetime_plan.v12.json` records the
scope. It uses only the three verified outer-layer states, never failed released
wavefunctions. The original campaign's verdict and no-intermediate-lifetime rule
remain unchanged. Native rate execution ran only on Sirius, single thread,
nice19, in `runs/lifetime_frozen_prior_n4_v1`; executed drivers are preserved.

Native E1/M2/E3 to ground and M1/E2/M3 to quintet completed. Independent replay
checks exact selected states, transition wavefunctions/mixing/CSFs/labels/nucleus,
raw native controls and complete 13-channel inventory. All original outputs are
retained in `frozen_complete_v12`: 121 paths, one 4,452,586-byte archive part.

| Gauge | Total one-photon rate (s^-1) | Model lifetime (ns) |
|---|---|---|
| B | 416251568.56400895 | 2.402393349410827 |
| C | 595114956.842919 | 1.6803476177190935 |

These are exploratory model outputs, not an accepted lifetime. Both gauges are
retained; neither is selected by proximity to experiment. No observed energies,
fitting, empirical uncertainty or accuracy conclusion enters the calculation.
The substantial gauge dependence persists. EXPOSED_NONBLIND; scientific HOLD;
adoption false; accepted lifetime, uncertainty, covariance and calibrated sigma
null. No grade, canonical line/product, reference, star, publication or merge.

## Verification and continuation

`scripts/rya1354_frozen_layer.py` shares the source audit between native execution
and offline replay; `scripts/rya1354_frozen_summary.py` checks the complete rates.
Machine records are `source_curation_v2/frozen_layer_source_v12.json` and
`source_curation_v2/frozen_lifetime_v12.json`. The existing focused suite plus
four source checks passed 109 tests; the final six-test frozen-layer suite also
passed, adding complete-decay replay and a substituted-wave rejection. Negative
checks reject altered odd-state optimization controls even with updated receipt
hashes, and reject relabelling the failed parent as successful.

Production/reference/catalog comparisons against task base are unchanged and
`git diff --check` passes. Native GRASP source remains clean. Bibliography audit
still fails on the previously recorded civis2013/mallinson2024 verification
vocabulary errors before library traversal; no library-audit success is claimed.

Next work should test model convergence with a separately declared larger
balanced active layer or independently implemented relativistic method, preserving
this frozen-layer definition and the failed full release as distinct evidence.
Do not extend a failed wave or turn gauge/model spread into calibrated error bars.
Known-grade external validation and a genuinely blind test remain owed.
