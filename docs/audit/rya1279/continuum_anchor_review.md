# RYA-1281: frozen-spectrum continuum-reference decomposition

Predeclared2026-11-04 before evaluating new diagnostics. Mac analysis only, no RT,
abundance/stellar-parameter/RV/LSF fit or production change. Keep the canonical
5772/g4.438 refined spectrum and the separately converged g4.338 spectrum frozen.
Use existing Kitt corrected holding, AIR comparison/internal vacuum convention,
nominal q90/3-MAD/linear/+/-5A rule, line exclusion+/-0.5A, same broadening/RV
and IAG reader/source as the previous controls. Preserve both references and
every refusal. No replacement continuum or favorable case/reference selection.

Trace the existing estimator exactly, including each clipping pass, the support
actually used for its returned coefficients and any subsequent survivors. Require
parity with native fit_model_guided and stored prior reports. Do not modify that
shared estimator. Record centre extrapolation and linear-prediction weights as
geometry only, with no noise/covariance/physical uncertainty inferred.

Use both selectors and both ratio denominators in a declared2x2 diagnostic.
Selection effects include any induced clipping. Separately freeze each nominal
coefficient-support mask and compare both denominators with no further clipping;
these are counterfactual diagnostics, not an alternative adopted estimator.
Split original nominal selections by side using the same minimum-pixel and
level bounds; insufficient sides refuse rather than acquiring invented anchors.
If the fifth clipping pass changes survivors after its coefficient fit, record
the actual support distinction and measure an explicit post-clip-refit diagnostic
without altering emitted products or silently changing the rule.

Reference source arrays stay external/read-only; archive derived selections,
relative ratios/residuals/influence and summaries, not raw external spectra.
Source/native code hashes and full identities recorded. No line/blend is
identified from wavelength proximity or a reference-near answer.

Initial HEAD9ce9be3e6aa3a28839b5ed881c61e95300e2b7cf, main f8c723a8 unchanged,
stacked base4635c742/draftPR598. Existing physical/reference/continuum/3D HOLDs
remain; no new adopted abundance or default continuum choice.
