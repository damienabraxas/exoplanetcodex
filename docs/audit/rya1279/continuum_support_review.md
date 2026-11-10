# RYA-1281: explicit continuum support and IAG quality provenance

2026-11-04 local; Mac only. Parent RYA-1278, draft stacked PR598. Starting
HEAD7a39500d65da3ad0668a9c59350f9df8901b6c25, refreshed main f8c723a8 unchanged,
stacked base4635c7421a7dee8671ffefe05d771c68da1cd24e. No RT or adopted result.
Apply the existing frozen-spectrum plan in continuum_anchor_review.md; all prior
sources, AIR/vacuum convention, native rule, physical inputs and refusals retained.

## Additive provenance, no estimator change

Only the audit adapter changes. It emits coefficient_support_count,
postfit_survivor_count, clipping_candidate_count, clipping_pass_count,
clipping_stop_reason and clipping_limit_reached. The native fit_model_guided,
LocalContinuum records and returned coefficients/continuum remain unchanged.
Postfit survivors mean the retained indices after accepted clips; a rejected
under-minimum candidate is counted separately, not mislabelled as the fit support.
Insufficient initial anchors emit support0/passes0/candidate null and their
retained original count. Disabling clipping emits CLIPPING_DISABLED; all-retained
termination on the fifth pass does not count as hitting a clipping limit.

Two frozen spectra/three O777 VIS lines/both references replay successfully.
Every old field in all74fit records is exactly unchanged (only adapter fingerprint
and the six additive fields differ). Twelve direct native/prior replays and
continuum array equality still pass. Counts:38ALL_RETAINED,26CLIPPING_DISABLED
(including2postclip diagnostics),6INSUFFICIENT_ORIGINAL_ANCHORS,
4ITERATION_LIMIT_AFTER_CLIP. The latter are synthesis7774 and its identical
blue-only side control in each physical case, not four independent discoveries.
Each uses73coefficient anchors and70survivors. The original source defect stays
characterized; this audit exposes it without silently refitting or fixing native
scientific behavior. Its diagnostic refit remains~-0.927mA and unadopted.

The derived pixel CSV and both figures are byte-identical to continuum_anchor_v1;
keep the existing committed copies rather than duplicate them. New report records
the updated adapter hash and source hashes. Full immutable outputs remain in
experiments.bride_continuum_support; report copied byte-for-byte to
[data/results](../../../data/results/rya1281/continuum_support_v1/report.json).

## Held-paper and IAG FITS inspection

Inspected the held PDFs Baker2020 (baker2020, ApJS247:24) and Reiners2016
(reiners2016, A&A587:A65), resolved from refs.library and fingerprinted in the
quality report. No download, new citation, raw spectrum copy or extra correction.
Baker p3 describes approximately R=1e6; p7 describes quality0 as robust,1/2 as
>10%/>25% telluric absorption and3 as saturated, with nonzero categories marking
spline-replaced regions. That is a source description, not proof that every
flag0 pixel is free from residuals or that the uncertainty array is covariance.

The held FITS has v/s/err/flags columns; no inspected header card supplies a
calibrated LSF/apodization. Read-only quality audit converts native vacuum
wavenumber to vacuum A then AIR once, matching the existing reader. All nine
scopes have exact zero flags:

| AIR centre A | +/-0.5A aperture | +/-5A continuum | +/-6A full helper input |
|---|---:|---:|---:|
|7771.944|110|1099|1319|
|7774.166|110|1098|1317|
|7775.388|110|1098|1317|

Counts overlap; they are not independent sample counts or a statistical bound.
The +/-6A scope covers all native samples used by solar_atlas_reference before
Gaussian filtering for the +/-5A window, including the helper's extra1A margin.
This rules out nonzero *delivered flags* in these helper inputs, not unidentified
correction errors, a physical opacity deficit, or a particular preferred atlas.

Unexpected source-contract issue: an initial whole-file categorical assertion
refused before creating any output. The held flags contain167615fractional values,
424distinct values spanning0-3 globally. No nonfinite flags. The paper/header
inspection does not establish why; aggregation is a hypothesis, not a finding.
The final audit preserves floating flag values and inventories the discrepancy;
never round to integers or infer categorical quality globally. All local O777
scopes remain exactly0 despite this global issue. Concrete characterized bug-ledger
row rya1281-iag-quality-v1 records the source mismatch and initial refusal. Original
file untouched; no downstream production claim affected.

## Broadening intake and next bounded control

Current solar_atlas_reference convolves the already finite-resolution IAG spectrum
with a Gaussian FWHM=lambda/catalog_R (500000 here); it does not remove native
resolution in quadrature or establish an FTS kernel. Catalog IAG R700000 and
Baker's approximate1e6 are different evidence; neither is an exact measured
kernel at777nm. Reiners Fig1 uses Gaussian degradation to HARPS110000 for a
comparison; that is not a calibration of this Kitt/IAG FTS kernel. No sinc kernel,
apodization or resolution is inferred merely from the instrument name.

Next diagnostic plan: freeze the same spectra, compare the shipped Gaussian leg
with separately labelled finite-native-resolution Gaussian controls at both the
catalog and published approximate source R, preserving zero-extra-smoothing as
an additional reference sensitivity. Require broadening-only changes, preserve
all refusals and native nominal leg, and track both continuum support and fixed-mask
denominator changes. These approximations measure sensitivity, not a calibrated
FTS correction; source-specific LSF/apodization provenance remains required before
adoption. Do not fit R, select a preferred case, or change the nominal continuum.

Broader xi/T/g/composition/interactions, independent O atom/EOS/reference, additional
O diagnostics, C/N then Fe, reference star, full3D and all bands remain owed.
Scientific HOLD; no abundance, canonical parameters, native estimator, source
spectra, production, merge or publication change. RYA-1281 In Progress and
RYA-1278 In Review.

## Commands and verification

```sh
source config/hosts/mac.env
/Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_continuum_anchors.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-continuum-support-v1
/Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_iag_quality.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-iag-quality-v1
```

Both final commands exit0. Reproduction requires fresh output directories; both
original immutable outputs remain. Actual writer invocations registered in
GENERATORS.yaml. No figures changed; both were inspected in the prior phase and
verified byte-identical here. Project-environment focused suite12passed in0.16s.
Expanded science environment suite66passed in7.01s,292existing printUpdate
warnings. It includes the previous63tests plus minimum-anchor candidate refusal,
explicit stop modes and recursive preservation of every archived scientific field.

Self-review: additive metadata only; loop arithmetic/clip decisions unchanged.
Rejected candidate count is separate from retained support. Fifth stable pass is
not called capped. Source quality flags preserve raw numeric values, no rounding;
local zero flags cannot clear broader telluric/reference gates. Native code and
physical hashes unchanged; no production behavior affected.

Final checks: MIRROR targets, v212 pointer and staged whitespace pass. Generator guard retains exactly2436existing diagnostics, zero added/removed versus the previous anchor audit. Fresh library audit194documents/127accounted/67unindexed/12off-Mac registered paths. Ledger suite still1failed/1passed on the historical RYA-587 empty science_impact row; previous current-main reproduction remains applicable (main unchanged). No unrelated ledger/gold/tracker debt repaired or hidden.

Exact expanded test command (science environment, after sourcing Mac profile):

```sh
/Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q tests/test_bride_continuum_audit_rya1281.py tests/test_bride_compute_budget_rya1281.py tests/test_bride_iteration_trace_rya1281.py tests/test_bride_ng_analysis_rya1281.py tests/test_bride_profile_geometry_rya1281.py tests/test_bride_thermal_initial_rya1281.py tests/test_bride_thermal_refinement_rya1281.py tests/test_bride_thermal_sweep_rya1281.py tests/test_bride_xi_curve_rya1281.py tests/test_bride_xi_sweep_rya1281.py tests/test_bride_solar_control_rya1281.py tests/test_bride_fe_stark_rya1281.py tests/test_bride_fe_pressure_rya1281.py tests/test_bride_trace_fe_lines_rya1281.py tests/test_bride_top_cap_rya1281.py tests/test_bride_observables_rya1281.py tests/test_bride_reference_rya1279.py tests/test_rya1230_local_continuum.py
```
