# CN profile normalization and SI unit diagnostic

2026-11-08, RYA-1281/1278, scientific HOLD. Starting2291a074; branch
ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, worktree codex/rya1279,
stacked base4635c742. Plan/kernel/probe precommittedb40d2581.

## What was checked

Compile an independent tiny caller against unchanged native newvoigt.f, verified
byte for byte against pinned iSpec21269948750e2177fe32c5df38e1d9bb81363374.
No full native spectrum program, source formula or installed binary changes.
The prior accepted CN report and state CSVs match their committed archives; CN/all
state CSVs remain byte-identical.168native states are three previously selected
isotope lines times56depths, not all651molecular rows or all CN lines.

The compiled caller evaluates8472pairs:3864signed-offset samples (168states x23)
and4608area-quadrature points. Native REAL inputs/outputs and original double
intermediates remain; Python checks exact returned input values, count/order,
finite nonnegative results and exact signed symmetry. Damping comes from native
state, not a newly selected width. Offsets are dimensionless Doppler units, not an
assumed AIR/vacuum wavelength grid.

Across sampled offsets[-100,100], native versus SciPy Faddeeva H differences:

| Measurement | Result |
|---|---:|
|Maximum relative difference|8.21270493e-5|
|Maximum absolute H difference|3.03469376e-5|
|Worst relative state/offset|109 / -1Doppler unit|
|Native signed symmetry|exact|

No profile-accuracy admission threshold was declared or invented after the result;
this is a measured approximation difference, not a pass against a physical budget.
The diagnostic kernel uses unmodified SciPy Faddeeva, with no fitted correction.

## Normalization and source-comment discrepancy

Native newvoigt.f:20comments claim area1/sqrt(pi). Measured H area is instead
approximately sqrt(pi), consistent with native bsyn prefactor and an H/(sqrt(pi)*
Dnu) frequency-normalized profile. The comment must not supply adapter normalization.
This is a characterized documentation inconsistency; native code/source are preserved.

Gauss-Legendre integration on v=tan(pi*u/2), u in[-1,1], samples min/median/max a
and two orders. Reported values divide H area by sqrt(pi):

| a | native512 | native1024 | Faddeeva1024 |
|---|---:|---:|---:|
|1.08679572e-7|1.00000484176|1.00000465106|0.99999999933|
|2.10850139e-5|1.00000484260|1.00000464821|0.99999999933|
|0.0707858950|0.99999686839|0.99999684197|0.99999999924|

Quadrature nodes are passed through native default-REAL input; its effect is retained
in both native and reference sums. This is finite-order numerical evidence for
sampled a, not a rigorous far-wing/error bound or physical normalization calibration.
Analytic Gaussian-limit and frequency-area tests separately check the Python kernel.

## Explicit external-state SI kernel

`pipeline/bride_cn_diagnostic.py:external_state_extinction` requires supplied lower
population per m3, gf, stimulated-emission factor, damping a, Doppler width in Hz
and frequency-integrated strength in m2Hz. It returns extinction in m^-1 for
explicit dimensionless offsets. It infers no chemistry, isotope factor, partition,
frame, abundance or physical constant. No background-provider class or LTE source
function is introduced. Existing atomic/Ni kernels and default providers are untouched.

Two routes checked on the same native state:

1. Native H times pre-Voigt amplitude/standard opacity times standard cm2/g opacity
   times density g/cm3, then100 converts cm^-1 to m^-1.
2. Kernel population conversion per gram times g/cm3 times1e6 -> per m3;
   native integrated strength pi*e_esu^2/(m*c) times1e-4 -> m2Hz;
   frequency profile H/(sqrt(pi)*Dnu), gf and stimulated emission.

Native physical constants are retained for dimensional bookkeeping; replacing them
with modern LW constants is not silently included. To isolate units, the comparison
accounts for the separately measured native/Faddeeva profile ratio. That ratio is
**not applied to kernel outputs**, spectra or an adopted model. Raw native and raw
diagnostic extinctions remain separate columns in profile_samples.csv.

Maximum unit-accounting relative difference1.81866841e-7 passes predeclared
128float32epsilon=1.52587890625e-5 screening bound. The bound is arithmetic screening,
not a physical uncertainty or error budget. No threshold or coefficient adjusted.

## Verification and reproduction

Acceptedv1exit0. Compile0.55296s, native caller0.38921s, excluding setup/quadrature.
Python3.12.13, NumPy2.2.6, SciPy1.18.0, project Mac environment. No new peak-memory
measurement or end-to-end performance claim. No RT/SE/spectrum/abundance calculation.
No Sirius/CI/installation/download. All original states/sources unchanged.

Preflight8tests passed0.28s. Final12passed/1skipped0.49s, exit0; unchanged LTE
blend suite explicitly skips because project environment lacks LW. New external-state
kernel's8checks pass, covering Gaussian limit, frequency integral, population scaling,
zero population, negative widths/damping, nonfinite populations, invalid stimulated
factor and shape errors. Existing CN source/reporting tests also pass. No claim of
LW backend integration; the standalone kernel is not yet a provider.

After source config/hosts/mac.env, use OMP_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1:

```
/Users/ryanschmitt/codex/venv/bin/python scripts/check_bride_cn_profiles.py --output <fresh-output>
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q tests/test_bride_cn_diagnostic_rya1281.py tests/test_bride_cn_coefficients_rya1281.py tests/test_bride_cn_interface_rya1281.py tests/test_bride_lte_blend_rya1281.py
```

Registered external experiment holds private source/probe, raw logs and report.
Five committed outputs (report/profile samples/probe wrapper/input/stdout) byte-match
external originals, largest711648bytes. Original native source and probe binary
stay external, without renaming or guard bypass. Actual generator registered.
No execution failure in this phase; source-comment inconsistency characterized in
run bug ledger without editing native source or declaring an upstream physics defect.

## Next check

Match actual native AIR-grid linear frequency offsets versus exact vacuum-frequency
coordinates, with explicit one-time frame conversion and depth-dependent Doppler
widths. Then exercise an experimental external-state CN provider on an unchanged
Bride column with matched continuum/state inputs. Native CN pressures remain an
explicit external input; this does not solve or validate Bride chemistry. Primary
release/frame/statistical conventions, EOS/composition consistency, independent
physical reference and solar/all-band/adaptive3D gates remain open. No molecular
model admission, abundance/product/default/production/merge/publication change.

Final repository guards: MIRROR/version/sequence/whitespace pass. Generator set
exactly2436historical diagnostics, zero added/removed. Existing global bug-ledger
RYA587blank science_impact and tracker/gold debts remain; no global-green claim.
Last library audit in this continuing session remains196documents/127accounted/
69unindexed/12off-Mac; incomplete, no new literature used.
