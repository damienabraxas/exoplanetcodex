# RYA-1281: frozen IAG reference-resolution sensitivity

Predeclared2026-11-04 before evaluating the new controls. Starting HEAD468e543f,
mainf8c723a8 unchanged, stacked base4635c742/draftPR598. Mac only; no RT.
Follow the next-control plan in continuum_support_review.md and the new manifest
config/experiments/rya1281_resolution_audit.json. Both canonical and lower-g
spectra remain frozen; observations/source arrays/canonical parameters untouched.

Four labelled IAG legs: shipped Gaussian at targetR500000; Gaussian quadrature
with assumed nativeR700000 (catalog) or1000000 (Baker approximate); and no extra
smoothing. Native finite widths are assumptions, not exact measured FTS kernels.
Do not adopt, fit R, deconvolve, choose a preferred leg or infer a sinc/apodization
from the instrument name. Source/native/reference/quality fingerprints retained.
Same AIR comparison, internal RT vacuum, macro/rotation/gravitational RV, q90,
linear/3MAD/min20/max3percent/five-pass rule and original margins. No extra
telluric correction; two previous RV residual legs stay separate and unchanged.

Require exact shipped Gaussian array parity against the actual helper, every
variant trace against native fit_model_guided, and exact synthesis/shipped-IAG
record parity with previous controls. Cross baseline/variant selectors and
denominators; accounting includes induced clipping. Freeze actual nominal-IAG
and synthesis coefficient-support masks separately with no further clipping.
Keep every refusal; do not count source/geometry sensitivity as independent
uncertainty or prove physical causes. Archive derived reports/summary/plots only.

Initial exact-array test detected1.11e-16 differences from reordering equivalent
width arithmetic. Shipped width now uses the actual helper operation order;
regression requires array identity rather than a relaxed tolerance. No real
diagnostic was run with the discrepant adapter.
