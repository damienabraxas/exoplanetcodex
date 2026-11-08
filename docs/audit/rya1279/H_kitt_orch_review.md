# RYA-1281 paired-H oxygen controls against Kitt Peak

Scientific HOLD. The small H-normalization response survives the established
broadening/ORCH diagnostic path, but does not resolve the much larger residual
against Kitt Peak. Strong continuum-reference/support sensitivity remains.
This is a replay of the pinned native local-continuum estimator, not a full ORCH
pipeline run, abundance fit or replacement Fe/CN-blanketed product.

## Fixed scope, provenance and reproduction

Plan `config/experiments/rya1281_H_kitt_orch.json` selects six existing bare-background
O spectra: both completed v2 coupled H cases and both H cases under each of the
two fixed-population anchors. Their input/report/population/flux receipts and
matched named physical states are checked. No new RT, EOS or SE solve occurs.
The source column remains MARCS5750/logg4.5/xi1, source O8.70 and41-level atom;
it is not a fully matched5772K solar model or an independently qualified atom.

Actual native ORCH source is checked against Git blob at
de96ae51fd8c52f79e5457fc0205e4baaa4596fd, the previously isolated RYA1315
checkout. Its model bound0.15/envelope0.03, MODEL_Q90 and envelope percentile95
are checked directly; constants are not edited. Existing iSpec effects source
hash matches the prior solar comparison. Existing broadening is retained:
vmac3.8km/s, vsini1.8km/s, limb darkening0.6, R500000. Both prior RV controls
(zero and canonical solar gravitational redshift) remain, with no optimization.

The exact corrected/pre-normalized Kitt Peak holding is
`solar_kpno_kurucz2005_corrected`, file `irradrelwl.dat`, SHA256
75dd93479b5d3fe63c2e7828f5dbb3ce29f4c672f4f1147e0d30ce41bc8f5adb.
`gate_holding` and the pre-normalized guard pass. The established reader converts
delivered vacuum nm to canonical AIR Angstrom once. Model outputs are already
AIR; no second model conversion is applied here. The source observation is
unchanged; no additional telluric correction or observation wavelength shift.
The established IAG holding and actual pinned helper/reader are fingerprinted.

Model coverage7763–7783A contains complete ±5A continuum support for all three
windows. The observation is clipped to that model coverage, and full continuum
support is explicitly checked. Native iSpec broadening/Doppler interpolation
retains zeros at unused outer edges; every pixel inside any requested ±5A window
must remain positive, all outputs finite and nonnegative everywhere. Padding
indices/wavelengths are recorded, never filled or normalized away.

The synthesis and IAG reference legs each use the actual pinned estimator.
For each RV/line/reference, the baseline coupled-model-H continuum is also held
common across ALL six models. This exposes model profile changes separately
from rerunning each model's own synthesis-reference continuum. The shipped
observation path is retained too. IAG records/arrays must be exactly identical
across model labels. Refusals would be recorded without division; actual result
is72applied/0refused, with72metric rows and36paired contrasts.

Load `config/hosts/mac.env`, then use the existing project analysis environment:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/analyze_bride_H_kitt_orch.py --output <fresh-H-kitt-orch>
```

This differs from the LW virtualenv used for RT: project NumPy2.2.6,
pandas3.0.3, SciPy1.18.0 on Python3.12.13. No dependency install, benchmark,
download, Sirius or CI. Matched broadening/reference/source guards pass.
Statistics are unweighted sampled residuals and finite ±0.5A areas; no noise
model, chi-square, significance, isolated-line EW or abundance error is inferred.

## Solar diagnostic results

Predeclared plot is plus-gravity/IAG, with the common baseline continuum. All
other RV/reference/shipped/common/per-model paths are retained in summary.csv.
Model-minus-observed signed area deficits are shown below (mA):

|AIR centre(A)|Coupled model-H area|Coupled native-H area|IAG-conditioned Kitt area|Model-H deficit|Native-H deficit|
|---|---:|---:|---:|---:|---:|
|7771.944|68.74211|68.93640|81.73581|-12.99370|-12.79942|
|7774.166|59.20498|59.38202|70.81716|-11.61219|-11.43514|
|7775.388|46.28401|46.43328|57.70338|-11.41937|-11.27010|

The paired-H processed area contrasts are0.19428240/0.17704363/0.14927031mA.
They are small relative to these11–13mA conditional residuals; they do not
explain the mismatch in this test. Plus-gravity/IAG RMS residuals change from
0.01403671/0.01269352/0.01220809 to0.01377827/0.01244298/0.01200283.
The plotted line cores remain shallower than the conditioned observation, with
structured residuals. No preferred H model, optimized RV or abundance follows.

IAG continuum levels are1.00455155/1.00389236/1.00414847, exactly common to
all model labels and independent of the two model-RV choices. Rerunning each
model's synthesis-reference continuum changes the paired observed area by at
most0.0002054576422096943mA across all36contrasts. Thus the paired H response
is not being generated by a large model-induced change in continuum placement.
This does not validate the absolute placement or remove reference sensitivity.

## Large absolute continuum-reference/support sensitivity

Plus-gravity baseline synthesis-reference levels are
0.99676926/1.03447271/1.13185422. The corresponding conditioned observed areas
are74.56667/98.28582/163.97663mA, versus IAG81.73581/70.81716/57.70338mA.
The pinned15percent bound accepts them all; acceptance does not establish a
physically trustworthy continuum or support geometry. No continuum is adopted.

Actual selected/coefficient-support indices and reference/continuum arrays are
saved for all72fits. The existing trace helper reconstructs the same coefficients
as the pinned estimator under a prospective floating-point screen. Its legacy
application bound is explicitly ignored for this arithmetic trace; only the
actual pinned estimator decides application. No bound is silently changed.

|AIR centre(A), plus-gravity baseline|Synthesis blue/red coefficient support|Centre outside support hull(A)|Centre prediction weight L1|Clipping stop|
|---|---:|---:|---:|---|
|7771.944|69/0|4.10169|14.24245|ALL_RETAINED|
|7774.166|0/90|4.10349|15.18194|ALL_RETAINED|
|7775.388|0/81|4.10116|16.69177|ITERATION_LIMIT_AFTER_CLIP|

At7775.388,81pixels set the final coefficients but only80survive the last
clipping proposal. The trace preserves that distinction rather than refitting
the survivors. IAG support at7774/7775 brackets the centre(68/22 and33/57,
weight L1approximately1); IAG7771 is still blue-only88/0, extrapolating1.35244A,
L1approximately2.75484. Neither reference is declared universally validated.

The synthesis anchors therefore extrapolate from one side, roughly4.1–5A away.
This supports continued continuum/support and missing-opacity controls before
interpreting absolute solar residuals. Incomplete Fe/CN opacity, source atom,
atmosphere and collision/cap hypotheses remain confounders; no single physical
cause is established. This extends the previously characterized support/clipping
diagnostic to these bare-background H controls, not a newly adjudicated ORCH
production defect or reason to change the ratified estimator here.

## Attempts, tests, artifacts and next task

v1 stopped before input ingestion: wrong RT-only interpreter lacked pandas.
v2 used the correct environment but an overbroad positivity guard rejected
iSpec zero edge padding outside all analysis support. Source inspection confirms
`np.interp(...left=0,right=0)` and `zero_edges=True` behavior. No science pixels
or spectra were altered; the corrected domain guard requires positivity in
every requested continuum window and rejects nonfinite/negative outputs anywhere.
Two explicit regression tests cover allowed padding and forbidden in-window zeros.
Both harness issues are recorded as fixed in the bug ledger. Logs v1–v4 are
preserved at `/private/tmp/rya1281-H-kitt-orch-vN.log`.

v3 exits0. v4adds actual coefficient support/reference-array preservation and
also exits0. summary.csv, contrasts.csv and the plot are byte-identical between
v3/v4. Both original external result directories are registered; no overwrite.
Six accepted artifacts byte-match v4 originals; max2807662bytes, below5MiB.
They include diagnostic observation segments, not a copied raw atlas holding.
Actual writer/invocation is registered; no split, rename or guard bypass.

Tests:12project tests pass0.20s; the five held-LW continuum-component tests pass
1.66s separately; actual pinned native continuum tests8pass0.09s. The earlier
project invocation's LW-only skip is resolved by that separate held-LW run.
Main remains7f023fd7, already inspected PR615Ca notes; no merge/rebase.
Fresh library212documents/127accounted/85unindexed/12off-Mac,exit1; no new
paper used. Historical global ledger/tracker/gold debts remain.

Next bounded task: cross common synthesis/IAG denominator with fixed support
that brackets the line, retaining the original and exact coefficient supports
as controls. Predeclare support geometry/coverage/refusal criteria; compare
fixed-support and native selected-support outcomes with no adopted continuum,
threshold/RV/abundance tuning. Keep real missing-opacity/provider qualification
separate: primary/EOS/native-state CN gates must pass before insertion, and no
synthetic purpose guard is bypassed. Additional oxygen indicators, C/N then Fe,
reference star and adaptive3D/all-band goals remain open. Scientific HOLD;
no source/default/product/production/publication/merge or abundance admission.

Final MIRROR/version-pointer/sequence/staged-whitespace checks pass. Generator
invocation audit retains exactly2436historical diagnostics, zero added/removed.
Global bug-ledger tests remain1failed/1passed on the first old RYA587blank
science_impact row, verified byte-for-field against main; new rows are valid.
This is not a globally green claim. Postcommit stacked-base freshness and normal
push hooks are required; no bypass.
