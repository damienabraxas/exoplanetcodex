# Independent Fe I background diagnosis near277nm

RYA-1281 / parent RYA-1278, 2026-10-29. **HOLD**. Mac diagnostics; no O NLTE
rerun, opacity substitution, abundance fitting, production or native producer claim.

## Controlled build and input checks

Held iSpec source21269948750e2177fe32c5df38e1d9bb81363374 is clean for all used
Turbospectrum sources. Files are checked against Git and copied before building.
Baseline and reporting-only derivatives use the same existing GNU Fortran16.1.0
and native flags, plus the declared modern-compiler argument compatibility flag.
All original sources, binaries and DATA remain untouched. Isolated DATA copies
and the exact MARCS5750/4.5 model are hash-receipted; all92 model abundances are
explicit, global rescaling zero and xi1km/s. The iSpec cleanup wrapper is not used.

The derivative enables NTO/NTPO reporting, enlarges diagnostic fields and writes
existing Fe population/EOS and post-DETABS values. No opacity or rate formula is
changed. At both coordinate hypotheses, baseline and reporting produce
**byte-identical serialized continuum.opa files**. This proves preservation at
native output precision, not bit identity of unexported internal variables.

All56 temperatures and electron pressures equal the exact prediction from native
default-REAL input followed by its E14.7 output formatting. Maximum temperature
quantization is0.0004K; Pe relative difference from decimal input is8.06e-8.
No Pe-convergence recovery message appears. Recomputed density differs by up to
0.06660%; this EOS difference is retained and reported. It is not overwritten
with MARCS density. Each requested three-point window becomes103points because
native Babsma deliberately adds100 wavelength samples. Actual domains are saved.

## Component result

The original MARCS OPA frame is still UNVERIFIED. Both conditional interpretations
are retained. TS wavelengths are air by the inspected source; the canonical
converter maps them to backend vacuum coordinates. Local LTE coefficients and Fe
populations are recomputed at the exact exported TS sample. Native OPA absorption
is linearly sampled at the corresponding conditional coordinate, without fitting.
At the nominal2769.76A source sample, the vacuum-hypothesis TS air coordinate is
2768.942A; the mapped source coordinate is2769.759783A. Kitt Peak remains AIR.

| logtau500 | TS/MARCS absorption | TS Fe I fraction of final absorption | TS/local Fe I opacity | TS/local neutral Fe per kg | TS/local net opacity per neutral Fe |
|---|---:|---:|---:|---:|---:|
|-4.9177|1.001189|0.738630|84.3667|1.038419|81.2454|
|-3.0030|1.002852|0.603011|82.8322|1.034290|80.0861|
|-0.9737|1.005092|0.387336|81.3477|1.034009|78.6722|
|-0.0102|1.003924|0.116304|82.7752|1.053635|78.5616|

The air-hypothesis top result is similarly TS/MARCS1.000798 and Fe opacity ratio
84.2226. Thus the independent TS absorption closely reproduces MARCS at this
sample, while the local Fe contribution is much smaller. Neutral Fe differences
of only3.4-5.4% at these depths do not explain the roughly80-fold component ratio.
Net opacity per neutral atom includes **level weighting and cross-sections**;
it is not a directly measured single-level photoionization cross-section.

The held local Fe model has32levels (15Fe I,16Fe II,1Fe III), with15 explicit
continua connecting the Fe I levels to Fe II. TS's actual opacity input labels
the Fe I table as OP data with cut resonances; DETABS identifies it as Fe I bf.
These are different continuum representations. The earlier local-only Si/Mg
budgets correctly identified what dominates **our** remaining absorption, but
did not identify what MARCS supplied. The new independent component evidence
instead prioritizes Fe I background completeness. This is oxygen-background
validation, not the later Fe abundance milestone.

## Bookkeeping and failures preserved

Native DETABS reports an intermediate opacity sum. ABSKO then adds improved H
bound-free and NH opacity. The final analyzer keeps those additions separate
and checks their closure against serialized total absorption. At the top they
are2.30307e-7 and4.89117e-7m2/kg, respectively. Fe I is73.9213% of the intermediate
sum and73.8630% of the final total. Applying intermediate fractions to the final
total would misassign those additions; no such estimate is admitted here.

Attemptv1 was rejected by the new reader's incorrect positive-geometry assumption:
TS stores a signed geometric coordinate. Attemptv2 was rejected by an exact
decimal-temperature check; native single-precision conversion was then checked
explicitly rather than loosening a tolerance. Attemptv3 exposed a genuine dormant
native reporting defect:44components overflow the original36field format, causing
format reversion to a REAL descriptor for a CHARACTER field. The derivative
expands reporting fields; a compiled minimal regression reproduces the old failure
and validates the corrected capacity. Original source is preserved. Later audit
versions add Fe-state and post-DETABS receipts. Final runv8 / analysisv3 are
immutable; earlier versions and their failed outputs remain external.

## Validation and next control

128 targeted tests pass with666 known printUpdate warnings, including actual
backend Kirchhoff/closure, native format failure, parser completeness, invalid
components and post-DETABS handling. Final four native cases total0.83s on the
Mac, one thread. No owned native job remains running. Review artifacts are under
data/results/rya1281/ts_budget_v1; full logs, DATA copies, source and compiled
oracles remain external with hashes. Reproduce using registered build/run/analyze
commands and **new immutable directories**.

Use the Mac host profile and the Bride venv. Build with
`scripts/build_bride_ts_budget.py --output <new-build>`, point
CODEX_BRIDE_TS_BUDGET_BUILD at that verified build, then run
`scripts/run_bride_ts_budget.py --output <new-run>` and
`scripts/analyze_bride_ts_budget.py --run <new-run>/report.json --output <new-analysis>`.

Next obtain a source-derived Fe I opacity packet across the O lower-triplet
photoionization support207.832-277.248nm (and upper support221.784-302.645nm).
Keep active O separate, subtract the existing passive Fe I contribution before
adding an alternative, preserve LTE emissivity/Kirchhoff and existing scattering.
Compare raw specific-opacity and per-neutral-Fe normalization as explicit EOS
controls; never fit an opacity multiplier to departures or abundance. A bounded
NUV-only intervention is not complete native background replication. Then run
fresh UV J/photoionization, population/rate and visible O diagnostics on the Mac.
Actual MULTI producer EOS/opacity/rates and full3D/all-band admission remain open.
