# RYA-1281: fixed-O natural FeII transfer differential

Predeclared2026-11-05 before evaluating new transfer. Startd19af347, refreshed
main9a74d843 unchanged, stacked base4635c742/draftPR598. Mac only, one thread.
User authorizes continued bounded oxygen controls; no product/adoption/merge.

Manifest config/experiments/rya1281_feii_transfer.json selects gf_079010 after
source and LTE strength qualification. Reconstruct exact source row at runtime,
require both frozen strength screens to pass, compare current main target row,
refuse ambiguity/isotopes/source changes. Add only its natural-Voigt LTE FeII
bound-bound opacity/emission, same transfer with existing41FeI pressure/Cowley
and bound background unchanged. No pressure/Stark guessed for FeII; no full-opacity
claim. No source/width/gf/abundance fit, no O/charge/SE feedback.

Both canonical/lower-g columns use saved O populations, same physical inputs,
rays/qscale/flux quadrature and dense source grid. Empty wrapper must reproduce
existing background chi/eta/sca exactly. Each baseline transfer must reproduce
saved normalized flux within1e-6, the original outer criterion, before interpreting
an extra-line result. No relaxed thresholds or best-case selection. Paired fixed
radiation/cap convergence1e-8, max1000,180s per leg; all old failures retained.

Reuse the exact hash-verified saved standalone continuum Fnu for both legs; no
refit or re-normalization. Both continuum references and both fixed RV controls
retained in post-analysis. Historical3percent continuum source remains frozen;
main's15percent model-guided rule is not exercised or represented by these fits.
All profiles in AIR; transfer internally vacuum; canonical gravitational shift
only in observed comparison. No additional telluric correction.

Two initial real-backend tests pass: empty-wrapper identity, FeII-stage-specific
coefficients/thermal emission, unchanged scattering and outside-support arrays,
deepcopy, cached-field refusal and invalid/unsupported line data. Local-opacity
ratios from the previous screen do not predict flux changes in this test.
