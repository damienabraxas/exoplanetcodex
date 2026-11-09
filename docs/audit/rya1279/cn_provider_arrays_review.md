# Actual LW synthetic CN provider array checks

2026-11-08, RYA-1281/1278, scientific HOLD. Start0c5d1843; branch
ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, worktree codex/rya1279,
stacked base4635c742. Plan/implementation2c781002; finite-positive diagnostic
comparison hardening1d342708. No formal solution, RT/SE or abundance calculation.

## Integration actually exercised

New opt-in `pipeline/bride_cn_provider.py` wraps installed LW0.16.1FastBackground
on a **disposable synthetic** H-only three-depth atmosphere. It is deliberately
restricted to purpose=synthetic_array_test and an H-only atomic set. Non-synthetic packet purposes and any atom set beyond H are refused; no native
56depthpacket/solar context was used. This is not their adapter or chemistry.
No existing provider/default/model/production route is replaced.

The immutable packet binds depth/height,T,ne,nHTot,xi,vlos and the92element
abundance vector. It supplies synthetic lower populations,gf,width,a,stimulated
emission and integrated SI strength. Coordinates are explicit exact vacuum nm;
canonical AIR conversion supplies the synthetic line centre once. The kernel
adds extinction plus LTE B_nu emission on declared support, leaving scattering.
No partition, isotope factor or external chemistry is inferred.

Plan has T5000/5500/6000K,ne1/2/3e15m^-3,nH1/2/3e20m^-3,xi1000m/s and
synthetic lower populations1/2/3e12m^-3; gf.01,mass26amu,gamma1e7s^-1. These
are disposable numerical inputs, not measured or source-qualified CN populations.
Actual LW arrays are692wavelengths x3depths, with41samples inside777–779vacuum nm.
This coarse support sampling tests array wiring, not line-resolution convergence.

## Passed checks

| Check | Evidence |
|---|---|
|CN disabled|chi/eta/sca exactly equal bare installed FastBackground|
|Zero synthetic population|chi/eta/sca exactly equal bare background|
|Enabled addition|Both extinction/emission agree with requested additions within16eps floating addition bounds|
|Outside support and scattering|Exact equality to bare background|
|LTE Kirchhoff|Added eta/chi agrees with independent LW planck function to1.33226763e-15relative|
|Repeated background update|Exact arrays, no cumulative addition|
|Deepcopy|Exact values and independent array storage|
|Stale state/grid/composition|Nine negative bindings refused before any output array write|
|Duplicate inventory/mapping/purpose|Explicitly present/unknown CN baseline, wrong mapping/support and non-synthetic purpose refused|

Nine negative bindings: height,T,ne,nHTot,xi,vlos,wavelength,runtime abundance
change and an incorrectly bound packet abundance vector. Snapshot validation runs
before baseline computation, and tests verify arrays stay unchanged after refusal.
Isolated AtomicAbundance instances prevent mutation of global default composition.

Additive residuals consume at most0.0150617(extinction)/0.0128208(emission) of the
16eps floating-addition bound; these are rounding-accounting results, not physical
uncertainties. The enabled sample is explicitly positive and Planck/reference
comparisons must be finite, preventing an all-zero/NaN comparison from passing.

## Baseline molecular inventory and binary/source boundary

Actual installed FastBackground requests CH,H,H-,H2,OH in these tests and never CN.
An audit proxy refuses any CN lookup instead of assuming it can be added twice.
The factory also requires an explicit absent-CN inventory assertion. This is scoped
to this H-only synthetic baseline; it does not audit an arbitrary future background.

Inspected local LW checkout d33058f8857acf28d9187e778b6f41e49d3aad95 is clean.
LwMiddleLayer.pyx:1386ff FastBackground has CH/OH/H2 fields/hooks, and builds
atomic background continua/resonance scattering; inspected background source has
CH/OH bound-free functions. No CN hook was found in that inspected path.

| Inspected source | SHA256 |
|---|---|
|Source/LwMiddleLayer.pyx|13428eab7aa0e0e2374015895cef73f60482aa30200cf87c09fae3942c5ab57a|
|Source/FastBackground.cpp|0a9ca1862f6d1d0b4609772d3c6fd6f78da84277e297dfae4c426d17592a6b8d|
|Source/Background.cpp|d4d0f92e28ebf0dc5f82f4ffa7e6798acee9363023e2727689ca6e875031f048|

Installed0.16.1module fingerprint is separately recorded:4127408bytes,
SHA256b86652aecb2b1979e62bcfa4c3fa9095eb3ef6616b77ced9e3c9d22e4af2a952,
unchanged before/after run. Distribution has no direct_url receipt, so source
checkout equivalence to that binary is **not independently rebuild-attested**.
Runtime lookups and source inspection support this diagnostic path without
manufacturing a stronger lineage/physical validation claim. No engine rebuilt.

## Attempts, environment, checks and outputs

Preflight2actual-backend tests pass1.88s; after expanded negative controls2pass2.31s.
Predeclared CLIv1exit0. Source review then required explicit finite/nonzero Kirchhoff
comparison; v1real outputs were already finite and positive. Freshv2exit0 with
packet.json and arrays.npz byte-identical to v1. Both preserved, no physics/bound
change. Ledger records the diagnostic finite-check gap as fixed, no affected spectrum.
An exploratory introspection smoke constructed a context but its print failed
because a Cython object has no __dict__; no accepted result used from that print.

Held LW environment Python3.12.13/NumPy2.5.2/LW0.16.1; OMP/OPENBLAS1 and Context
Nthreads1. v2main-block time1.80967s includes LW import/array checks/output packing,
excludes process startup and top-level NumPy import. Darwin process peak RSS
210518016bytes covers the whole diagnostic Python process, not per-context/native
program memory. Existing LW default-config warning retained; no benchmark invoked.

Final held-LW suite13passed1.77s, exit0. Project environment10passed/1explicit
LW module skip0.31s. No claim of skipped integration passing in project environment.
Tests use actual installed Context/background machinery and independent LW Planck
function; Context initialization computes LTE equilibrium/backgrounds, no formal
solution or population iteration was called. No Sirius/CI/install/download.

After source config/hosts/mac.env:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/check_bride_cn_provider.py --output <fresh-output>
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q tests/test_bride_cn_provider_rya1281.py tests/test_bride_cn_diagnostic_rya1281.py tests/test_bride_cn_coordinates_rya1281.py tests/test_bride_lte_blend_rya1281.py
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q tests/test_bride_cn_provider_rya1281.py tests/test_bride_cn_diagnostic_rya1281.py tests/test_bride_cn_coordinates_rya1281.py
```

Registered first/accepted output directories retain original synthetic outputs.
Three compact artifacts (packet,all baseline/disabled/enabled/zero arrays and report)
are archived byte for byte fromv2, largest208138bytes. No backend binary or copied
native source committed; actual generator and invocation registered. All real source
spectra, previous native results, solar products and observed frames remain unchanged.

## Next matched-state task

Establish an explicit native56depth binding packet, including full-network hydrogen
normalization/electron density/mass density/composition and microturbulence receipts,
before permitting native CN populations in a Bride column. The current sampled CN
state export alone is not an independent H/EOS/composition match. Extend packet
purpose/atom-set support only with those receipts and matched-array controls; never
silently put5750/4.5native populations into canonical5772/4.438solar state.

Then test opacity/background arrays at that matched native state before a new O
NLTE diagnostic. Physical primary frame/isotope/partition/line-list provenance,
chemistry/reference parity, solar abundance, C/N then Fe/reference star and actual
adaptive3D/all-band goals remain OPEN. No scientific model/default/product,
production/publication/merge admission. Scientific HOLD remains.

Integration audit: shared origin/main advanced tod70d932a during this phase,
merging Ca close-out PR612and RYA1352grade audit/design PR613. No provider/CN/
wavelength code, canonical molecular parser or ratified science policy changed
in the inspected delta. Sequence now records those landings. Theoretical evidence
design is not treated as molecular-primary admission. No merge/rebase of this branch.

RYA1352also preserves an unresolved historical Codex/Deep grading authority
conflict; it is not adjudicated here. This synthetic array experiment selects no
graded pool or theoretical transition. RYA1353/1355own related O/Fe census/stress
work; no duplicate atomic-grade investigation is started.

Final guards: MIRROR/version/sequence/whitespace pass after the new main integration
record. Generator diagnostic set exactly2436historical, zero added/removed. Existing
global ledger RYA587blank science_impact and tracker/gold debts remain; no global
green claim. Last continuing-session library audit196documents/127accounted/
69unindexed/12off-Mac remains incomplete; no new literature used.
