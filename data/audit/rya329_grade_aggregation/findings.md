# RYA-329 — line_score in the EW abundance aggregation (gate + weighted-median)

Branch `ryandamienschmitt/rya-329-grade-weighted-aggregation`. Implement on branch, NO merge,
re-bank HELD for the coordinated bundle. Decision (Ryan): Option 3 — line_score-weighted
**median** as primary, **but gated on line_score actually discriminating**.

## Gate result — line_score does NOT discriminate → primary STAYS plain median
The discrimination gate (run FIRST, per Ryan) fails on the solar Fe pools:

| species | n(A/B/C/D) | u.mean | p.median | w.mean | w.median | A+B | wmed−med | ls_spread | ρ_noncirc | discriminates |
|---|---|---|---|---|---|---|---|---|---|---|
| Fe I | 0/46/16/0 | 7.551 | **7.510** | 7.539 | 7.508 | 7.497 | −0.002 | 0.102 | −0.036 | **False** |
| Fe II | 0/3/0/0 | 7.714 | **7.657** | 7.715 | 7.657 | 7.714 | +0.000 | 0.007 | — | **False** |

Two independent reasons the signal is untrustworthy:
1. **3 of 5 sub-scores are inert** (never populated): `fit_chi2_score` (flat 0.50 default),
   `saturation_score` (flat 0.00), `nlte_correction_score` (flat 0.00). Only `ew_snr_score`
   and `abundance_outlier_score` vary.
2. **The apparent discrimination is circular.** `line_score` vs |a−median| looks strong
   (Spearman −0.44) but that is driven by `abundance_outlier_score`, which is itself a
   function of |a−median|. **De-circularised** (line_score minus the outlier component) the
   correlation is **−0.036 (p=0.78)** — no discrimination.

Weighting the median by a mostly-inert, partly-circular score would self-reinforce toward the
median (artificially tighten), not correct quality. So per the gate: **the primary A_X/A_X_nlte
stays the plain (grade-blind but outlier-robust) median — unchanged, no anchor move, no
re-bank.** The scoring repair (populate the 3 inert sub-scores; de-circularise) is owed to
**RYA-220** before the primary can be weighted.

## What shipped (regardless of the gate)
- `pipeline/abundances_derive.py`: `_weighted_median`, `_line_score_discriminates` (the gate),
  `_aggregation_diagnostics` (the 5 estimators); wired into `_element_grade_summary`, which now
  emits per element the loud diagnostics table + gate verdict and writes `agg_*` / `ls_*`
  columns to `solar_abundances.csv`. The primary is weighted **only if the gate passes** —
  it does not, so A_X/A_X_nlte are untouched (verified: solar Fe I A_X 7.510 / A_X_nlte 7.516,
  Fe gate still PASS).
- `tests/test_grade_weighted_aggregation_rya329.py` (4 pass): weighted-median correctness;
  gate fails on flat/circular, passes on a genuinely discriminating signal; inert detection.

## Provenance reconciliation (Ryan's note)
- The banked solar **Fe I = 7.510 is the plain MEDIAN = A_X (1D-LTE)**; **7.516 is A_X_nlte**
  (the NLTE-scale value we've been quoting from 509/336). Both stand under this fix (median
  unchanged).
- **Fe II 7.657 = the EW-path plain median** (the mean is 7.714). It is blend-biased HIGH
  (RYA-352); the pipeline already demotes it — the solar Fe gate uses the **synthesis arbiter
  7.500** for ionization (RYA-405/406), not the EW 7.657. So RYA-509 scored the EW-path Fe II,
  i.e. it ran the EW-path value RYA-406 demoted — consistent with 509 predating the 406 posture
  in the reported per-element EW number.

## Re-bank — HELD for the bundle
No re-bank here: the gate kept the primary as the median, so nothing moved. Even once RYA-220
repairs the scoring and the gate passes, per Ryan the re-bank lands as ONE coordinated pass
(confirm 406/407/336 merged + RYA-446 Fe I scatter threshold + RYA-515 provenance stamp →
one clean full-solar re-run → v1→v2 once → RYA-251 Phase-1 sign-off).
