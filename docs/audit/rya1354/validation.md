# RYA-1354 validation checkpoint

Base: `adefcaf073ff70035954cf45a9a2283774ac7d06`.
Branch: `ryandamienschmitt/rya-1354-oxygen-pilot-blind-mcdhf`.
Mac worktree: `codex/rya1354`; all solver builds/runs: Sirius task work root.
Final SHA/PR and actual CI state are recorded in the Linear handoff.

| Command/check | Result |
| --- | --- |
| `git fetch origin` | 0; confirmed RYA-1352/1353 merged before branching. |
| Sirius compiler/resource inventory | gfortran 15.2.0, CMake 4.2.3, 4 logical CPUs, 15 GiB RAM. Initial load 1.93/2.06/3.33; GRASP build limited to two jobs. No other worker stopped. |
| GRASP configure/build | Initial BLAS lookup failed; exact installed BLAS/LAPACK sonames resolved it. Build exit 0; upstream source clean, warnings retained. |
| GRASP quadrature and preserved serial workflow | Both exit 0; quadrature absolute discrepancy 1.93e-14, workflow `TESTS SUCCEEDED`. Not an oxygen benchmark/full suite. |
| ATSP book build | Exit 0; legacy warnings retained. Excited-state MCHF and HF run, but old angular transition result is scientifically refused. |
| Modern ATSP2K src90 build | Top-level make exits 0 despite unavailable MPI utility sub-build errors; not counted as a successful full build. Explicit serial libraries and `biotr` subsequently build. |
| Native modern Si E2 control | Process exit 0, no semantic refusal, output is byte-identical to pinned upstream `examples/tr_ls/1De1.1Se1.14_4.ls` (SHA256 4818b57d009243227117535f6550326f3a004b199c06ae50eee962f8830720fd). Reproduce with `scripts/rya1354_atsp_control.py`. |
| Final O diagnostic `scripts/rya1354_solver_smoke.py`, exact CLI in README | Exit 1 intentionally: rejected book angular path is retained. GRASP E1 and native modern biorthogonal MCHF produce finite exploratory results. No blind benchmark/convergence/calibration success claimed. |
| `python scripts/rya1354_summarize.py --run /tmp/rya1354-final-capsule/solver_smoke_v10 --out <new-file>` | Exit 0, three candidate rows, admission HOLD; exact replay covered by tests. Candidate v3 adds calculation/native-table hashes and explicit unpublished DOI status; older v1/v2 retained. |
| `python -m pytest tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q` | 48 passed in the existing Mac Python 3.12.13 environment. Tests are metadata/arithmetic/fixture tests, not Mac solver runs. |
| `git diff --exit-code adefcaf0 -- pipeline data/products data/linelists data/reference data/catalog config/stars.yaml` | Exit 0: existing production code, canonical/reference data, catalogs, stellar parameters and all product bytes unchanged. New research scripts/config are opt-in and have no production consumers. |
| `python3.12 scripts/check_register_freshness.py` | Exit 0; existing stale generated tracker notice retained. Post-commit PR mode is run separately. |
| `python -m pipeline.refinement_debt_join --report` | Exit 0; informational outstanding debt unchanged (21 registry rows, 10 not-fired ticket rows, five no-ticket rows). |
| `python3.12 scripts/generate_sources_page.py --audit-library <resolved-library>` | Exit 1 before library traversal on existing civis2013/mallinson2024 verified-state vocabulary failures. No clean-library verdict. |
| `python scripts/audit_line_keys_rya1037.py --check` | Exit 1: existing 31 findings across 12 files and excess waiver counts in nlte_corrections; none in new scripts. RYA-1353 already documented the same baseline set. |
| `git diff --check` | Exit 0. |
| Post-commit `python3.12 scripts/check_register_freshness.py --since-main` | Exit 0; zero state surfaces changed. |
| Ordinary `git push` | Exit 1: existing SEQUENCE drift omits merged RYA-1352/1353 on origin/main 7f023fd7. Version pointer/freshness pass. Used the hook's explicit `--no-verify` draft option; no merge/publication gate bypass. |

Main advanced to `7f023fd7` during the session with separate Ca feed/results and
SEQUENCE changes. This checkpoint remains based on the recorded `adefcaf0` input
state; those other workers' scientific changes are neither copied nor reverted.
The zero-production-diff assertion is against this task's base, not a claim that
the live repository stopped changing during the pilot.

Retained development failures: initial BLAS configure; optional 25.6 GB unique-label
allocation; old angular-domain refusal with false zero exit; raw book/native
coupling/header/coefficient incompatibilities; extra zero-rank operator request
from an unnecessary prompt reply. Each owned run is in a separate Sirius v1-v10
directory. No previous output was deleted or silently replaced. Final v10 requests
only E1; the exporter preserves scientific values. All owned solver processes
finished. Initial sandbox SSH/DNS failures used authorized retries; one read-only
shell quoting error was corrected without source changes.

Self-review: no experimental uncertainty manufactured, no gauge spread promoted to
sigma, no lifetime from an incomplete decay list, no NIST experimental pedigree
assigned to theory, no observer wavelength/gf/Solar abundance fitted. The old
invalid independent results remain rejected; modern limited multiplet agreement
is explicitly distinct from independent experimental validation. Two characterized
run-bug-ledger rows preserve upstream failure behaviors. Three bibliography rows
were added with primary metadata and limits; no website regeneration/deployment.

**Still owed:** curator and independently verified held-out experiments; approved
preregistration and thresholds; active-space/correlation/grid/relativistic
convergence; independently converged comparison; calibrated uncertainty and
multiplet covariance; UV/IR and separate M1/E2 calculations; scientific review.
No real experimental residuals/coverage statistics exist here. RYA-1354 must stay
unfinished; the draft PR is a reproducible checkpoint, not the full acceptance.
