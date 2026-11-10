# Kitt Peak residual and local-continuum controls

RYA-1281, 2026-11-02. Scientific HOLD. Fixed solar5772/g4.438/xi1/O8.70,
existing converged full-atom oxygen spectrum; no new SE, opacity or abundance fit.
All work on Mac. Shipped observation, prior spectra and comparisons unchanged.

## Important correction to earlier wording

The baseline retained the provider continuum without alteration. That remains a
valid frozen diagnostic, but is not the complete current ORCH continuum treatment.
`pipeline/prenormalised_guard.py` explicitly allows the bounded RYA-1230 local
window operation. `pipeline/local_continuum.py` defines the RYA-1232 nominal
robust synthesis-guided reference, with an observed-envelope estimator and an
IAG-reference sensitivity leg. A blanket prohibition on every local adjustment
would misstate that ratified exception. No production rule changed here.

The new audit calls existing estimators and holding/telluric guards, with the
same corrected2005 observation hash, AIR reader and declared model broadening.
Observed data loaded over7765-7782A supplies full +/-5A envelope neighborhoods.
No new telluric correction, baseline replacement, free continuum fit or tuning.

|AIR centre|Shipped observed mA|p95 local level|Envelope observed mA|Fixed model mA|Remaining deficit|
|---|---:|---:|---:|---:|---:|
|7771.944|77.55651|0.9962762|74.10889|71.97199|2.88%|
|7774.166|67.20039|0.9962957|63.73224|62.15361|2.48%|
|7775.388|53.79430|0.9959505|49.94709|48.77347|2.35%|

These envelope values are a continuum sensitivity diagnostic, not adopted oxygen
EWs or proof the supplied continuum is wrong. Real weak absorption also affects
local maxima. They explain why the earlier7-9% window gaps cannot be assigned
entirely to NLTE physics. No observational uncertainty/covariance supplied.

## Nominal model-guided control and coverage

Saved model spans7768-7780A. First and third centres lack complete +/-5A model
neighborhoods, so model-guided estimates are explicitly WITHHELD there. No
truncated neighborhood is silently substituted. The middle line has coverage,
but its model-guided level1.0344707 exceeds the existing3% bound and is refused.
All90 pre-clipping selected anchors lie7778.26949-7779.15925A, on the red side;
this narrow one-sided selection extrapolates to the centre. Final robust-fit
record retained. This characterizes a diagnostic limitation, not an asserted
upstream implementation defect. No refused continuum is divided into data.
Oxygen-only anchors omit weak absorption; wider full-opacity spectra are needed.

The gravitational-shift residual decomposes into blue wing[-.5,-.15],
core[-.15,+.15],red wing[+.15,+.5]A. Missing areas after envelope control:
1.2066/0.6616/0.2733mA,0.5454/0.9069/0.1333mA,0.4151/0.6762/0.0709mA.
Both zero-RV and declared gravitational-shift decompositions are archived;
partitions are diagnostic, not line identifications or independent fits.

## Reproduction and next work

`scripts/audit_bride_kitt_residuals.py --comparison data/results/rya1281/solar_bound_v1/kitt_peak/report.json --output <new-directory>`
uses registered Mac paths and project Python. Immutable external v1 retained;
v2 adds anchor geometry. Final compact report byte-identical at
`data/results/rya1281/kitt_residuals_v1/report.json`, actual generator registered.
Initial guard call used the wrong keyword and failed before any output; corrected
to actual API, no failed result accepted. No new verified science run defect.

Next widen the fixed-population output with adequate broadening padding and full
continuum neighborhoods; qualify weak bound-bound/molecular blend opacity before
using synthesis-guided anchors as a physical result. Compare the IAG reference
leg under existing holding/frame policy, then declared velocity/FTS-LSF controls
and matched microturbulence. Keep raw baseline and each control distinct. No
abundance fitting, adoption, merge/publication or all-band/full3D claim.
