# RYA-1281: native Ng instrumentation and bounded lower-gravity controls

Predeclared2026-11-04 before RT. Scientific HOLD. Two serial Mac attempts at
the same thermal_g4338 physical point and refined sampling: record native Ng
flags with original Ng enabled, then disable Ng. Same verified baseline starting
populations, exact certificate, physical operators, full56depths and original
1500iteration/600second, joint1e-6/conservation1e-10/inner1e-8 controls. No altered
atoms, collisions, abundances, continuum, source-domain extension or parameter fit.
No further retries with relaxed limits if either attempt fails.

Do not modify the inherited coupled solver, positivity retry, context builder,
pressure/Stark/provider kernels or source data. A forwarding Context records
the actual returned IterationUpdate.ngAccelerated flags per active species and
per native call, including retries. Record accepted/rejected transient population
validity as booleans, never clip values. Native on/off configuration and recorder
source identity are explicit in manifests. Missing flags refuse, not default false.

Test recorder return/state transparency and retry accounting on disposable fixtures
and the actual native backend before O science. Compare the recorded-on history
against the retained earlier lower-gravity timeout, requiring exact shared numeric
prefix as a no-op check. Host timing may differ; convergence count and physical
rate checks determine the outcome, not whether a wall-time deadline happened first.
Compare on/off histories, actual acceleration flags and independently rebuilt
saved R+C/particle totals. Only accepted complete spectra get aperture/profile
comparisons. Preserve every failure and all histories. No physical/production
admission from a passing replay or a faster numerical method.

Initial HEAD1ed25f7c16a12089c51e0a281427ac499a704c25; stacked base4635c742;
refreshed main f8c723a8 unchanged. Same branch/worktree and draft PR598.
New immutable root rya1281-ng-control-v1; no prior output overwritten.

## Implementation and execution

Recorder and experiment commit467d11cd precedes RT. Initial wrapper preflight
incorrectly required the not-yet-created output path to exist, exit1 before
creating a directory or starting RT. Fixed in f522be1e: codex_path resolves the
registered output identity, required inputs still use require_codex_path, and
existing/mismatched outputs refuse. Regression tests cover all three cases;
concrete setup defect appended to run_bug_ledger.csv. Original exception log
retained. This is one setup error, not a science attempt or a source fallback.

Native-backend fixtures compare plain and recorded Context trajectories with
Ng on and off. Population/conservation arrays and complete convergence histories
are array/dictionary-identical; actual Ng flags appear only in the enabled
fixture. Synthetic invalid-then-valid updates prove per-outer retry accounting
and unchanged return objects. Missing species flags refuse rather than default.
The original coupled solver, positivity retry and physical kernels stay unchanged.

57 selected tests pass in6.55s on held Mac Python, NumPy2.5.2/LW0.16.1; expected
native deprecation warnings retained. Same50 earlier regressions plus5recorder/
output tests and2analysis accounting/prefix tests. No extra stellar runs counted
as unit tests. Analysis requires an explicit Ng comparison to permit the single
ng control change; an undeclared altered policy still refuses.

Actual standalone batch invocation, after sourcing config/hosts/mac.env:

```sh
/Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/run_bride_ng_controls.py \
  --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-ng-control-v1
```

Its attempts.jsonl records the full child argv, code commit, return code, wall
time and stdout/stderr hash for both declared modes, continuing to the independent
second attempt after a first failure. Each child retains the same physical case,
certificate and starting arrays. Native controls.ng is explicitly true/false;
recorded flags describe actual calls, not requested configuration labels.
The worker/runner/recorder source bytes stay fixed throughout the batch. Completed
run telemetry includes recorder overhead; no overhead measurement is inferred
for failed runs without a completed report. Timings remain host-load dependent.

Resources: experiments.bride_ng_controls and existing thermal coefficients;
all original input bindings remain immutable. Reproduction needs a new output
resource override, never reuse of these output directories. Both runs on Mac,
one context thread, serial; no Sirius science/CI, installs or downloads.

Analysis command, with a fresh immutable output directory:

```sh
source config/hosts/mac.env
/Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/analyze_bride_ng_controls.py \
  --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-ng-review-v1
```

Controls check native-call/retry accounting, source/case/column/certificate,
actual saved starting arrays, complete joint histories, unchanged acceptance
limits and independent unscaled R+C/target-nuclei conservation. Only completed
accepted spectra receive area comparisons. Exact shared prefix checks must pass
without numerical tolerance relaxation. Failures remain explicit in all outputs.
