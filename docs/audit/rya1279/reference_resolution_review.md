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

## Results

Plan/adapter/script committed before evaluation at5af54f634f550e17ef2caa92ea774ecdeb4636f1.
Both cases completed on Mac, no RT.24direct native-IAG variant fits match traced
coefficients/counts/application and continuum arrays exactly. Six synthesis and
six shipped-IAG baseline records match their prior records exactly. The shipped
Gaussian native arrays also match the actual atlas helper exactly before any
interpolation. All96crossed and96frozen-support fits pass the unchanged numerical
continuum gate; no new refusal or adopted estimator.

Canonical changes from the shipped IAG continuum in mA below. Guided and
frozen-IAG values are identical in the lower-g context because the same observed
IAG leg is used; these repeats are not independent observational evidence. The
synthesis-mask result changes slightly when its actual mask changes at7775.

| AIR centre A | IAG control | Guided | Frozen IAG mask | Frozen synthesis mask | Initial/fit membership symmetric difference |
|---|---|---:|---:|---:|---:|
|7771.944|shipped_gaussian|+0.0000000|+0.0000000|+0.0000000|0/0|
|7771.944|finite_native_catalog|+0.0139235|-0.0185783|+0.0020851|6/6|
|7771.944|finite_native_baker_approx|+0.0473515|-0.0082314|+0.0004065|4/4|
|7771.944|no_extra_smoothing|-0.0903381|-0.0342089|+0.0055877|10/10|
|7774.166|shipped_gaussian|+0.0000000|+0.0000000|+0.0000000|0/0|
|7774.166|finite_native_catalog|-0.0247719|-0.0079627|-0.0058608|2/2|
|7774.166|finite_native_baker_approx|-0.0185031|-0.0035294|-0.0019935|2/2|
|7774.166|no_extra_smoothing|-0.0864697|-0.0146498|-0.0127842|6/6|
|7775.388|shipped_gaussian|+0.0000000|+0.0000000|+0.0000000|0/0|
|7775.388|finite_native_catalog|-0.0047994|-0.0051305|-0.0006448|4/4|
|7775.388|finite_native_baker_approx|-0.0022483|-0.0022483|-0.0003120|0/0|
|7775.388|no_extra_smoothing|-0.0085923|-0.0095245|-0.0010984|8/8|

Across both cases, max absolute guided response0.090338126mA; max fixed-IAG
support0.034208901mA and fixed-synthesis support0.012784211mA. Existing reference
spread remains4.882779-10.927325mA. These declared IAG Gaussian alternatives are
therefore too small to explain that spread. This is not a bound on the full
Kitt/model instrumental profile, all possible FTS kernels or physical model error.
No native FTS kernel/apodization has been calibrated or adopted.

Selection effects are nonmonotonic even with fixed spectra. At7771 the assumed
native1e6 control changes guided area+0.047351mA but changes fixed-IAG area
-0.008231mA, with4membership changes. The2x2 finite path accounting attributes
+0.056537selector and-0.009186denominator-plus-induced-clipping, closing exactly;
not independent noise/error terms or a causal physical parameter fit. With
assumed700000 at the same line,6membership changes give+0.013924guided versus
-0.018578fixed-IAG. This is why no preferred R or monotonic inverse interpretation
is selected. The largest max-reference-flux change in +/-5A is0.000851110 under
no-extra smoothing. Native shipped sigma~0.745pixels; finite700000~0.521 and
finite1e6~0.645; these are local AIR-grid Gaussian numerical widths, not measured
FTS profiles.

## Reproduction and artifacts

```sh
source config/hosts/mac.env
/Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_reference_resolution.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-reference-resolution-v1
```

Exit0. Original output immutable in experiments.bride_resolution_audit. Reproduce
only to a fresh registered output directory. Two derived results and two figures
are byte-checked copies under data/results/rya1281/reference_resolution_v1 and
data/output/rya1281/reference_resolution_v1. Report records generator/adapters,
native rule, plan, prior report, observation/reference and frozen model identities.
No raw external spectral arrays archived. Summary records all24case/line/variant
contexts, including the zero baseline and any refusals (none here); plot/source
counts are not independent experiments. Both plots inspected visually.

[Canonical sensitivity](../../../data/output/rya1281/reference_resolution_v1/canonical_refined_resolution.png),
[lower-g sensitivity](../../../data/output/rya1281/reference_resolution_v1/lower_gravity_refined_resolution.png).

Focused project-env preflight15passes. Expanded science-env suite68passed in8.13s,
292existing Lightweaver printUpdate deprecation warnings. Exact expanded command:
use the full18module command in continuum_support_review.md plus
 tests/test_bride_resolution_audit_rya1281.py. New tests cover exact shipped
filter parity, input immutability, no-extra copy, Gaussian variance relation and
explicit refusal of deconvolution/nonfinite or invalid native R. Initial parity
failure and operation-order fix are recorded as fixed in the global bug ledger,
run_id rya1281-resolution-shipped-parity; fix SHA5af54f63 precedes analysis.
Previous native clipping and fractional-IAG-flag discrepancies stay characterized.

Self-review: observed/model/reference identities fixed, default native rule
untouched; shipped control replays exactly. Finite-width square subtraction only
for nativeR>targetR; no negative-width clipping/deconvolution. Each variant uses
its native samples for direct fit parity; no double interpolation. Frozen masks
are actual coefficient support, not reported survivors. Missing/invalid source
or identity mismatch refuses; no fallback, refit, fitted source scaling or
favorable reference chosen. Crossed and fixed fits retain all application flags.

Next bounded work: inspect independent weak-opacity/reference evidence at the
actual influential continuum anchors (not candidate wavelength proximity alone),
retaining the one-sided7771/7774 extrapolation and full opacity/source qualification.
Source-specific FTS/apodization provenance and independent O atom/EOS/reference
still required. Broader xi/T/g/composition interactions, other O diagnostics,
C/N then Fe then reference star, full3D/all-band goals remain owed. Scientific
HOLD, RYA-1281 In Progress / RYA-1278 In Review; no adopted abundance, production,
shared continuum/default, merge, publication, source or canonical parameter change.
All compute on Mac, no Sirius jobs, new RT, package installs/downloads or CI run.

Final repository checks: MIRROR targets/v213 pointer/whitespace pass; generator guard retains exactly2436existing diagnostics, zero added/removed versus the prior support audit. Ledger suite remains1failed/1passed on the historical RYA-587 empty impact row, previously reproduced on unchanged current main. Fresh library audit194documents/127accounted/67unindexed/12off-Mac registered paths. Gold/tracker and library debts remain visible, no green CI claim. Largest archived file3472978bytes, below5MiB; no binary renaming or guard bypass.
