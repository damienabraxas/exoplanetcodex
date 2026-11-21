# RYA-1354 — full-eigenvector (mixing-cut) control, v45 failure and v46

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. This is a numerical control on existing models. No new SCF, orbitals, configurations, Hamiltonian, expansion, lifetime or grade.

Continue from local `552a1aa5` (v44). The v44 push awaits Ryan, because the pre-push hook refuses on unrelated SEQUENCE/register debt. Worktree `/Users/ryanschmitt/codex/rya1354`. Plans: `config/atomic/rya1354/lifetime_plan.v45.json`, `lifetime_plan.v46.json`. Prospective declarations: Linear `12f62504-38ba-4175-b729-e0d73452fb16` (v45) and `15910208-9aad-417b-a415-b16556c7b3f5` (v46).

## Why

v44 identified one truncation it could not quantify after the fact. FAC stores each level's eigenvector only for `|mix| ≥ 1e-5·|principal|` and then renormalises (`AddToLevels`, `structure.c:3039`). In the v41 models that drops 577/370 (ground) and 612/395 (equal) components of the lower 3P2 and upper 3S1 eigenvectors.

`mix_cut` is read only at that line. `SetMixCut(0)` therefore stores every Hamiltonian component without renormalising. Its second argument is omitted, so `mix_cut2` keeps its default of 1.0 and acts only with group filters anyway. ANGZCUT=1e-8 becomes the only truncation, and the v44 component records quantify every term it drops.

## Design

The v44 driver is parameterised by its frozen index (`mix_cut`, case subset, retained failure); with default arguments it reproduces the v44 index exactly. Scope:
- **Models:** ground_coherent100 and equal_coherent100, in both gauges (4 native processes). coherent100 is the tighter v41 control, and B and C share one CI.
- **Input:** the v44 input plus a single `SetMixCut(0)` immediately before `Structure`.
- **Declared consequences:** `o.basis` and `o_export_000.cm` (printed and exported mixing) and `raw.tsv` (rates) change.
- **Must be unchanged:** every other stable v41 file byte-for-byte; level energies and labels in `o.lev` (apart from the session stamp); radial preflight, post-CI orbitals, potential and configuration state.
- **Recorded, not compared:** the coefficient transport and the independent rate audit.
- **Still required:** every v44 exact component identity.

## v45 — FAILED_RETAINED

Receipt `1a5785e0c0e445bc2db969a130d07939ac874b40b4f8edbdf7b73c4ba9f628d3`; 1:03.79, 1067592 KiB, 0 swaps.
- The canary passed: all 56 stable v41 files were byte-identical.
- ground_coherent100_B ran to completion natively, storing 916/787 components.
- The unchanged v41 CI review then raised `Declared radial CI input differs`. It requires its exact CI prefix, and the declared `SetMixCut` line sits inside that prefix. No later case ran.

This is a review-mechanics defect in the v45 driver, not a native or physics failure. The run is preserved in full in `data/audit/rya1354/fac_mixcut_control_failed_v45/` and recorded in `data/audit/run_bug_ledger.csv` (`fac_mixcut_control_v45_review_prefix`). It is not relabelled as passed.

## v46 — revision and result

Only the review mechanics change. The driver verifies that the input is exactly the v44 input plus the one declared line immediately before `Structure`. It then runs the unchanged v41 CI review on a read-only symlink view whose `input.sf` is that v44 input. The retained v45 case must reproduce byte-for-byte, including `o.basis`, `.cm`, `raw.tsv` and `components.tsv`.

A read-only dry check of the corrected review on the retained v45 outputs passed before declaration. It is not counted as the run.

```sh
# Sirius, /mnt/codex-data/codex/rya1354
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 \
  /usr/bin/time -v -o evidence/fac_mixcut_control_driver_v46/batch.time \
  python3 evidence/fac_mixcut_control_driver_v46/rya1354_fac_e1_components.py \
  --root . --source toolchains/fac-9ab6def \
  --index evidence/fac_mixcut_control_driver_v46/index.json \
  --index-sha256 ef884d92d757e85598248a294dac1c66ffe0bcd426b04f4e9c100ac4c07c59e8 \
  --probe toolchains/fac-e1-components-probe-v46 --out runs/fac_mixcut_control_v46
```

Index preparation took 33.27 s / 410124 KiB. The batch exited 0 in 2:38.70 at 1067952 KiB, with no swaps.
- The canary passed again.
- All 4 cases pass every guard: 53 stable files byte-identical to v41, level table unchanged, CI state unchanged apart from transport, and every component identity holds.
- The retained v45 ground_coherent100_B outputs reproduce byte-for-byte.
- Full eigenvectors are stored: lower/upper n_basis equals the Hamiltonian dimension (916/787 for 3P2), with 0 dropped.

## Finding

Relative change from storing the full eigenvectors (v46 vs v44, coherent100):

| Channel (lower) | ground B | ground C | equal B | equal C |
| --- | ---: | ---: | ---: | ---: |
| 3P2 amplitude | +1.6e-7 | +1.1e-6 | −1.1e-7 | +6.4e-7 |
| 3P1 amplitude | +7.1e-7 | +2.2e-6 | −6.5e-7 | +5.1e-6 |
| 3P0 amplitude | −2.7e-7 | +7.0e-8 | −1.8e-6 | −4.2e-6 |
| 1D2 rate | +1.9e-3 | −4.2e-4 | +1.1e-3 | −1.9e-3 |
| 1S0 rate | −1.0e-4 | −8.5e-5 | −1.6e-3 | −1.4e-2 |
| total rate from 3S1 | +6.0e-7 | +2.7e-6 | −9.4e-7 | +3.2e-6 |

The stored-mixing truncation moves the dominant 3P E1 amplitudes by at most about 5e-6 relative. The total decay rate from 3S1 moves by at most 3.2e-6. The ~11% B/C amplitude split (C/B 0.892 ground, 0.702 equal) and the 2p3 3p ← 2p3 3s mechanism are therefore **not** truncation artifacts. The largest configuration-pair shifts are ~1e-6 and come from minor correlation pairs.

The angular-cut-dropped amplitude stays at or below 3.9e-6 of the 3P2 amplitude. Only the cancellation-dominated intercombination lines (v44 cancellation indices 190–1200) move at the 1e-3–1e-2 level. Their absolute rates are negligible. This also shows they are numerically fragile and must never be used as calibration.

## Evidence

Lossless v3 archives (`rya1354_native_parts.load_parts(..., repo_root=)`). Files byte-identical to an earlier committed archive (v41, v44, or the retained v45 archive) are stored as manifest-pinned references. Duplicates within a single run (gauge C's CI files equal gauge B's) are stored once and self-referenced (`archive: "."`). Reference plans come from `scripts/rya1354_reference_plan.py`.

| Archive | Stored files | Parts | Compressed bytes | References | Pack + independent decode |
| --- | ---: | ---: | ---: | ---: | --- |
| `fac_mixcut_control_failed_v45` | 66 | 102 | 144390248 | 372 external, 5 self | 19.18 s, 4.8 GB RSS |
| `fac_mixcut_control_complete_v46` | 84 | 112 | 155005142 | 540 external, 8 self | 22.48 s, 4.5 GB RSS |

All parts are ≤5 MiB. Most stored bytes are the full-eigenvector `o.basis`/`.cm` for the ground model (in v45) and the equal model (in v46). These are the scientific content of this control. The machine report `source_curation_v2/fac_mixcut_control_v46.json` is regenerated from the archives by `scripts/rya1354_fac_e1_components_report.py --version v46`, and a test requires exact equality.

## Validation

```sh
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q \
  tests/test_rya1354*.py tests/test_atomic_benchmark_rya1354.py \
  tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py \
  tests/test_line_match_rya1033.py
```

**264 passed in 247.91 s** (5.9 GB max RSS, Mac Python 3.12.13, source/metadata/report checks only). New tests:
- the single-`SetMixCut` insertion and declared-change set;
- a reproduction of the v45 review failure mode, plus the strict review view that fixes it;
- the v3 self-reference round-trip and its rejection of a bad reference;
- the actual v45+v46 archive replay with machine-report equality and the mixing-cut bound.

The new ledger row passes the ledger schema checks scoped to that row. The whole-ledger test still fails on rows that already existed at `552a1aa5`; they are not rewritten. Original FAC source diff is empty, original library/executable hashes are unchanged, and no owned native process remains. Throwaway check directories `~/scratch/rya1354_v44_compilecheck` and `~/scratch/rya1354_v46_reviewcheck` on Sirius are outside the repository and run trees. Canonical scientific surfaces are unchanged. No publication, merge or Done.

## Next bounded investigation (proposal — needs its own declaration)

Both numerical side-channels that v44 left open are now closed for the dominant lines: the photon grid (v44) and the stored mixing (v46). The remaining gauge disagreement is physical: incomplete correlation in the 2p→n′p channel of the lower 3P2 level. The next step is a declared correlation expansion with p (and companion s, d) correlation orbitals beyond the spectroscopic 3p/4p, optimised for the lower level. It should track the 2p3 n′p ← 2p3 3s configuration-pair amplitude and the B/C ratio with both potentials kept, using the v44 component instrument at each layer. Primary experimental provenance/covariance review and independent physical-convergence review remain owed before any grade or adoption.
