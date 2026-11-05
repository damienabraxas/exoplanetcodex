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

## Results and interpretation

Completed on Mac only, no RT. Plan/code committed before evaluation at
563929f3d3e647d7afb1fbcd8f58964483026bcb. All12 native/stored-report parity
checks pass exactly (record values and continuum arrays), across two frozen
spectra, three O777 VIS lines, and both synthesis/IAG references.24crossed
fits,24fixed-support denominator fits and24side-control attempts retained;
refusals are explicit. This is neither other-band nor other-element validation.
The Kitt delivered vacuum-nm grid is converted once to AIR for comparison;
internal RT remains vacuum. Existing broadening/gravitational RV is fixed;
previous two-RV residual controls remain separate, not refitted here.

Area accounting in mA (II minus SS). Selector and denominator are symmetric
finite path averages; denominator includes induced robust clipping. Their sum
is the total, while interaction is an additional diagnostic, not a third term.
Neither these terms nor the geometry below is a noise/covariance/error budget.

| Frozen case | AIR centre A | Total | Selector | Denominator + induced clipping | Interaction |
|---|---:|---:|---:|---:|---:|
| canonical_refined | 7771.944 | 7.147190 | 0.283535 | 6.863655 | -1.865723 |
| canonical_refined | 7774.166 | 4.882779 | -1.746441 | 6.629219 | -0.969706 |
| canonical_refined | 7775.388 | 10.927325 | 2.080936 | 8.846389 | -4.837056 |
| lower_gravity_refined | 7771.944 | 7.155995 | 0.277618 | 6.878376 | -1.853890 |
| lower_gravity_refined | 7774.166 | 4.925416 | -1.749652 | 6.675068 | -0.963282 |
| lower_gravity_refined | 7775.388 | 10.879344 | 2.009291 | 8.870053 | -4.747715 |

With nominal coefficient support frozen (no additional clipping), changing only
the denominator gives canonical synthesis-mask deltas8.312982/7.001586/11.264917mA
and IAG-mask deltas5.930793/6.155887/6.427861mA. Lower-g values are
8.321787/7.044223/11.243910 and5.951431/6.204941/6.496195mA. The large
spread persists without anchor reselection. This does not identify which
reference is physically correct, missing blend identities, or a preferred Sun.

Both cases' synthesis7771/7774 coefficient supports are entirely blueward,
69/73points, respectively4.10169/4.10429A from the nearest anchor to centre.
The centre-prediction absolute-weight sums are14.24245/14.90902, versus
IAG2.75484/1.0; geometry shows substantial extrapolation sensitivity, not
physical uncertainty. At7775 both references bracket centre (L1=1), yet still
show large denominator sensitivity, so one-sided support alone cannot explain
all reference differences. Both anchor figures inspected visually.

Six of24side attempts refuse for zero original red anchors (synthesis7771,
IAG7771, synthesis7774 in each case). Two synthesis7775blue fits refuse the
unchanged3% level bound (levels0.933923/0.939393). All16other side attempts
pass that numerical gate, not physical qualification. IAG7774blue/red areas
72.034332/69.887831mA differ2.146501mA; no side is selected as preferred.

## Characterized native defect; no shared-estimator fix

At7774 in both spectra the fifth clipping pass fits73anchors, then drops3;
returned coefficients still use73 but n_bins reports70 survivors without a
cap flag. All actual returned values remain preserved. Explicit no-more-clipping
refits on70survivors change centre levels by-0.000990122/-0.000990267 and areas
by-0.926805/-0.927026mA. This counterfactual is not an adopted correction and
must not silently replace the nominal rule. Concrete bug recorded as characterized
in data/audit/run_bug_ledger.csv under rya1281-continuum-anchor-v1, owned by
RYA-1281. Source is pipeline/local_continuum.py, the shared ORCH estimator;
no production scientific result is asserted affected. Regression reproduces the
support mismatch from archived derived ratios and checks the returned coefficients
use the actual support, rather than treating survivor counts as the fit support.

## Reproduction, checks and artifact identity

Environment: source config/hosts/mac.env. Actual analysis command:

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_continuum_anchors.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-continuum-anchor-v1
```

Exit0; original immutable registered output is experiments.bride_continuum_anchor.
Reproduction requires a fresh registered output directory, not overwriting v1.
Two results and two figures are byte-checked archive copies. JSON records source,
adapter, native code, plan, observation, IAG, prior-report and model fingerprints;
CSV contains derived ratios/residuals/support/prediction weights, no raw external
spectral flux arrays. Both original references and all refusals remain.

[Canonical anchors](../../../data/output/rya1281/continuum_anchor_v1/canonical_refined_anchors.png),
[lower-g anchors](../../../data/output/rya1281/continuum_anchor_v1/lower_gravity_refined_anchors.png).

Project-environment initial native/adapter tests:9passed. Expanded science-env
suite:63passed in6.84s,292existing Lightweaver printUpdate deprecation warnings.
Exact test invocation is the17previous endpoint/Ng/profile/physics test modules
listed in ng_extended_review.md plus tests/test_bride_continuum_audit_rya1281.py;
full command and test output /private/tmp/rya1281-anchor-tests.log (the log contains
output; module names are recoverable from the previous review and test tree).
No installs/downloads, heavy jobs, Sirius, production or CI run initiated.

Self-review: trace mirrors the actual five-pass native loop, including its cap
behavior; independent centre-weight reconstruction passes, and stored/native
parity covers every nominal fit. Crossed effects explicitly include clipping;
frozen masks isolate denominator effects without inventing anchors. Input arrays,
physics, abundance, native estimator and old refusals remain unchanged. All
source/physical/independent-O/EOS/continuum/LSF/full3D/all-band gates stay OPEN.

Next bounded work: expose coefficient-support vs postfit-survivor counts and cap
status as separate audit provenance, proving native coefficient/continuum identity.
Then source-backed FTS/reference/opacity controls before interpreting residuals
as stellar parameters; independent matched O reference/atom remains required.
Broader xi/Teff/g/composition/interactions, additional O diagnostics, C/N then Fe,
and a reference star remain owed. No abundance fit, estimator adoption, production,
merge or publication. RYA-1281 In Progress; RYA-1278 In Review; scientific HOLD.

Final repository checks: MIRROR targets and v211 pointer pass; staged whitespace
passes. Generator guard retains exactly2436existing diagnostics with zero added
or removed against the prior final log. Bug-ledger suite remains1failed/1passed
on the old RYA-587 empty science_impact row; identical failure reproduced against
the current-main fixture. New characterized row has a concrete science_impact.
Library audit this session finds194documents,127accounted,67unindexed and12
registered paths absent on Mac; no new literature inference relies on those gaps.
Refinement-debt audit remains10not-fired/5TBD/2discharged. These repository-wide
debts remain visible; they are not reported as passing checks.
