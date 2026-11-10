# RYA-1278 validation receipt

2026-10-25 workstation, base `0ade2170ebeda715f7b23ee8b0f57efdae5d0f3a`.
No RT/SE or abundance solve; snapshot inventory is a read-only ingestion audit.

## Commands and results

The two exact interpreters used were
`/Users/ryanschmitt/codex/venv/bin/python` (project) and
`/Users/ryanschmitt/bride_rya1118/.venv/bin/python` (held Bride evaluation).
These are execution receipts, not portable path configuration. Runtime discovery uses
`engines.bride_eval`; scientific inputs and machine-independent hashes are in
[`held_assets.json`](held_assets.json).

From the checkout, with `PYTHONPATH=.` and each interpreter:

```sh
python -m pytest -q -p no:cacheprovider tests/test_bride_foundation_rya1278.py tests/test_multi3d_read_rya1118.py tests/test_multi_atom_rya1119.py
```

- Final held Bride: **43 passed**, 2.42 s, exit 0. Includes distinct-star identity controls and actual pinned LW
  atmosphere construction with synthetic inputs; no Context/RT/SE invocation.
- Project: **37 passed, 6 skipped**, final project run, exit 0. LW is absent there; five
  existing converter checks and one new constructor check explicitly skip.
- Before adding the actual-constructor integration check: 41 passed in Bride;
  project 36 passed/5 skipped. Baseline reader/converter tests were 16 passed/
  5 skipped in project. Counts are distinct stages, not independent science proofs.

Additional project checks:

```sh
python -m pytest -q -p no:cacheprovider tests/test_run_bug_ledger_rya926.py tests/test_model_registry_rya1101.py tests/test_reference_discovery_rya1274.py
python scripts/verify_model_registry.py
python scripts/check_register_freshness.py
python scripts/generate_sources_page.py
python scripts/generate_sources_page.py --check
git diff --check
```

- Broader selected tests: **21 passed, 1 failed**, pytest exit 1. The failure is
  a **pre-existing** RYA-587 ledger row missing `science_impact`, symptom
  `Completed Reference VIS xi outputs absent from committed stamps`.
- Baseline reproduced against `git show HEAD:data/audit/run_bug_ledger.csv`, with
  the same committed test copied to `/tmp/test_rya1278_baseline_ledger.py` and
  only its LEDGER path redirected to `/tmp/rya1278-baseline-ledger.csv`:
  `python -m pytest -q -p no:cacheprovider --confcutdir=/tmp --rootdir=/tmp /tmp/test_rya1278_baseline_ledger.py`
  => **1 failed, 1 passed**, exit 1, identical assertion and row. Initial /tmp
  collection without confcutdir/rootdir hit a sandbox PermissionError; isolated
  rerun succeeded in reaching the test and confirmed the baseline failure.
- Model-registry/reference-discovery tests separately: **20 passed**, exit 0.
- Registry exit 0: nine models, seven live tokens, mandatory pair intact, no new
  token. It reports known legacy aliases as warnings, not silently omitted models.
- Register history check exit 0; known stale generated tracker notice persists
  (phase_c/gold-v3 provenance debt), not fixed or suppressed here.
- Local generated source page and its `--check`: exit 0, 208 references.
  An initial optional `--site-root /tmp/rya1278-sources-preview` attempt correctly
  refused the target because it was not a website checkout; no external site was
  written. Normal local generation was then rerun successfully. No deployment.
- Whitespace check exit 0. PR-mode freshness is checked after the implementation
  commit; its result is recorded in the Linear/PR handoff.

Read-only inventory, held interpreter, final command:

```sh
PYTHONPATH=. python scripts/audit_bride_foundation.py --snapshot-stem t5777g44m0005_00091_80x80x240 --output /tmp/rya1278-held-assets-v2.json
```

Exit 0. Captured output copied into the committed audit after inspecting it;
source hashes identify the exact generator and validator despite base HEAD at
capture. LW reported that its optional performance config was absent and chose
its default supported SIMD path. No optional benchmark or optimization ran.
The CLI's exclusive-create output contract prevents overwriting an existing file.

Standing library audit was run through the resolved `refs.library` holding with
`scripts/generate_sources_page.py --audit-library`: exit 1, 40 unindexed documents,
12 unresolved local_file declarations with candidates present. Existing RYA-1277
owns this library-wide limitation; held papers actually used were readable.
Refinement report (`python -m pipeline.refinement_debt_join --report`) exit 0,
informational: ten not-fired rows and five TBD rows remain. No CNO debt cleared.

Read-only Sirius inspection used `ssh sirius` to inspect uptime and the two
historical STAGGER directory aliases. Existing /srv holding found; /mnt alias
absent; missing queried directories made the combined command exit 2. No claim
that every remote environment was exhaustively inventoried. Load and users were
recorded in the roadmap; no heavy job launched or process killed.

## Review disposition

Reviewed units/frame, column identity and all-depth preservation, canonical star
ownership, positive/finite input rejection, immutable source copies, explicit boundary
and microturbulence arguments, absence of RT/SE calls, pinned API seam and source
provenance. Negative fixtures reject swapped shapes, reverse axes, invalid densities,
NaN/Inf, out-of-range indices and missing inputs. Actual backend construction
preserves supplied thermodynamics. These controls do not validate atom physics,
EOS equivalence, rate balance, solver numerics, observations or full 3D.

Two new ledger rows characterize the historical gravity-unit hazard and reconcile
stale engine metadata. Bundled multiplet collapse is expected library-model
structure, not a claimed LW defect; it is a mismatch to the proposed benchmark.
The unvalidated science paths remain HOLD through RYA-1279–1282. No new defect was
hidden behind the known baseline ledger failure. No protected inputs/results changed.
