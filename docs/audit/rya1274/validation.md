# RYA-1274 validation evidence

Base: `183310d8851629f43549e93269989032aa422bae`.
Host: workstation Mac; existing `codex/venv/bin/python`. No package installation.
In commands below, `PY` denotes that interpreter and `BASE` the isolated clean
worktree checked out at the full base SHA. The committed host profile was sourced
for the ORCH tests. No science engine dispatch occurred.

| Command | Exit | Observed result |
|---|---:|---|
| `git fetch origin` | 0 | Base SHA recorded; independent worktree created. |
| `PY run_pipeline.py --help` | 0 | Standalone CLI loads and documents matrix/sweep/prepare/publish flags. |
| `PY -m pytest -q -p no:cacheprovider tests/test_rya1222_run_matrix.py tests/test_rya767_run_descriptor.py tests/test_rya1233_sweep.py tests/test_rya1233_coverage_report.py tests/test_rya1234_element_verdict.py tests/test_rya1235_post_run_report.py` | 1 | 135 passed, 1 failed, 3 warnings, 66.97 s. Failure: `test_every_canonical_target_routes_to_exactly_one_dossier`; canonical Zn lacks route. Warnings: divide/overflow/invalid covariance matmul from an existing control. |
| In BASE: `PY -m pytest -q -p no:cacheprovider tests/test_rya1235_post_run_report.py::test_every_canonical_target_routes_to_exactly_one_dossier` | 1 | Identical assertion (Zn absent), 1 failed. |
| `PY -m pytest -q -p no:cacheprovider tests/test_reference_discovery_rya1274.py tests/test_ion_nlte_rya1250.py tests/test_solar_cited_nlte_rya1250.py tests/test_run_bug_ledger_rya926.py` | 1 | 30 passed, 1 failed: inherited ledger row's empty science_impact. |
| In BASE: `PY -m pytest -q -p no:cacheprovider tests/test_run_bug_ledger_rya926.py` | 1 | 1 passed, same malformed historical row failed. |
| `PY -m pytest -q -p no:cacheprovider tests/test_reference_discovery_rya1274.py tests/test_ion_nlte_rya1250.py tests/test_solar_cited_nlte_rya1250.py` | 0 | 29 passed, 3.92 s. Includes new override/missing/unindexed synthetic controls and existing ion, applicability, recipe/provenance tests. |
| `PY scripts/generate_sources_page.py --check` | 0 | Generated page matches bibliography, 206 rows. No page written. |
| `PY scripts/generate_sources_page.py --audit-library <resolved refs.library>` | 1 | 165 documents, 125 accounted, 40 unindexed; 12 local_file paths absent. Actual output diagnostics retained below; no claims of whole-library readiness. |
| `PY -m pipeline.refinement_debt_join --report` | 0 | Informational: 10 not-fired ticket rows, 5 TBD; snapshot state caveat recorded in audit. |
| `PY scripts/check_register_version_pointer.py` | 0 | One v172 pointer; 170 unique changelog rows. |
| `PY scripts/gen_state_register_targets.py --check` | 0 | Generated targets mirror unchanged/in sync. |
| `PY scripts/check_register_freshness.py` | 0 | Register at/ahead of state surfaces; inherited tracker generation notice remains nonfatal. |
| `gh api repos/damienabraxas/exoplanetcodex/actions/workflows/ci.yml --jq '{state,path,updated_at}'` | 0 | Workflow active. Initial sandbox network attempt failed; approved retry succeeded. No workflow triggered or runner status assumed. |
| `git diff --check` | 0 | Final whitespace check. |

Final post-commit PR freshness, link and protected-surface checks are recorded in the
Linear handoff. The suite is **not all green**. No full-suite/Sirius/production run
was claimed; these changes add a path-register entry and docs, not engine behavior.

## Synthetic ORCH observations

A temporary-directory fixture patched `run_matrix.load_feed/expand/graded_ions`,
`cell_process_hold/_manifest_paths/ledger_key/_band_of/LEDGER/REPORT_DIR`, the test
cell's deck emission, and `run_descriptor.resolve`. Resolver raised
`RuntimeError("synthetic unavailable dependency")`; actual `coverage(...)` returned:

```json
{"counts":{"PUBLISHED":1,"HELD":0,"NOT_RUN":0,"STALE":0},"done":true,"reason":""}
```

A second `unittest.mock.patch.object(run_matrix, "expand", return_value=[])` probe
called `element_publish.units_for("solar", "Si", report={"cells":[]},
interpreter="synthetic")`; its recorded `pools` argument was
`['set', 'codex', 'deep']`. These observations support RYA-1275/1276. They did not
execute scientific stages, write a real ledger or update a feed.

AST inspection of imports in the standalone driver, descriptor, matrix, sweep,
prepare, publisher, uncertainty contract, PySME adapter and Gerber adapter found no
OpenAI/Anthropic/LangChain/agent-service imports. This bounded check supports
independence of the inspected call path; it is not a proof about all external binaries.

## Exact unindexed-file diagnostics from the library audit

These names are diagnostics, not bibliographic identifications. Do not assign citation
metadata from a filename. RYA-1277 owns source inspection/alias classification.

- `1212.2649v1 (1).pdf`
- `1609.07283v1 (1).pdf`
- `1609.07283v1.pdf`
- `1808.09478v1 (1).pdf`
- `2301.11391v1 (1).pdf`
- `2301.11391v1 (2).pdf`
- `2301.11391v1 (3).pdf`
- `2602.14294v1 (1).pdf`
- `Codex_Process_and_Code_Map.drawio`
- `HARPS-telluric-STOP-2026-08-19.md`
- `Pehlivan_Rhodin_2017_MgI_IR.html`
- `Pehlivan_Rhodin_2017_MgI_IR_full.pdf`
- `ReadMe CH.txt`
- `ReadMe CN trasnsitions.txt`
- `ReadMe PArtion functions for molecules and atoms.txt`
- `ReadMeGaia FGK.txt`
- `TICKET-DRAFT-kpno-telluric-line-selection.md`
- `Table A1 Open in new tab Wavefunction... _ Oxford Academic.html`
- `Table A1. Open in new tab Wavefunction... _ Oxford Academic.html`
- `Table A10. Open in new tab Comparison of... _ Oxford Academic.html`
- `Table A2. Open in new tab Electric... _ Oxford Academic.html`
- `Table A3. Open in new tab Electric... _ Oxford Academic.html`
- `Table A4. Open in new tab Electric... _ Oxford Academic.html`
- `Table A5. Open in new tab Electric... _ Oxford Academic.html`
- `Table A6. Open in new tab Comparison of... _ Oxford Academic.html`
- `Table A7. Open in new tab Comparison of... _ Oxford Academic.html`
- `Table A8. Open in new tab Comparison of... _ Oxford Academic.html`
- `Table A9. Open in new tab Comparison of... _ Oxford Academic.html`
- `aa36265-19.pdf`
- `aa40445-21.pdf`
- `aa45394-22 (1).pdf`
- `aa45686-22 (1).pdf`
- `apj485713t1_ascii.txt`
- `apj485713t2_ascii.txt`
- `apj485713t3_ascii.txt`
- `apj485713t4_ascii.txt`
- `apj485713t5_ascii.txt`
- `apj485713t6_ascii.txt`
- `dataTable6.dat`
- `table7.dat.txt`
