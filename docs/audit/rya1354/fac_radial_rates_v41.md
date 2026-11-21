# RYA-1354 — expanded radial decay-rate milestone, v41

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false; accepted lifetime, uncertainty, covariance and calibrated sigma remain null. These are exploratory rates for the existing four 6,600-level models, not experimental validation, physical convergence or theoretical line grades.

## Declared calculation

Continue from v40 commit `19a32407` on `ryandamienschmitt/rya-1354-oxygen-pilot-blind-mcdhf`, worktree `/Users/ryanschmitt/codex/rya1354`. The surviving session owns the v40 handoff only; this session owns v41. The prospective declaration is Linear comment `584391ed-8ea4-4c0e-a51a-755b76c8cf6a` and `config/atomic/rya1354/lifetime_plan.v41.json`.

The full LS index has SHA256 `8ddec6ba6a2939402fad68f0f43fc39b77da3adaccd2052e1f1aea8f875105f0`. Its preparation verifies original sealed receipts, complete case-specific native labels/compositions, source C/CM and the original CI/radial replay. It retains the original full reports and links, without repeating the LS matrix/vector calculation. Preparation on Sirius took 221.47 seconds and 2485600 KiB peak RSS, no swaps. All eight actual rate input builders independently derive six lower states and 13 downward multipole requests per model.

Each fresh process restores the original paired saved potential, executes the exact qualified 114-group radial CI protocol with disabled corrections, then appends the existing complete operator enumeration for one gauge. No new Hamiltonian, optimization policy, observed-energy substitution or physical model is introduced. The CI reviewer defaults to its original exact input. Its new optional full-input argument also requires the entire original protocol as an unchanged prefix and exact equality to the declared input; a mutated Breit control is rejected in regression coverage.

Before rate totals, require all previous spinor/grid/protected-potential/configuration/signed-coefficient checks, equality to the original qualified report and byte-identical regenerated C/CM/source.pot/ci.pot. The independent rate audit checks every raw request and status, mapped transition energy, pinned formula, formatted float32 table row and source writer threshold. Tiny omitted table entries remain in the total; no missing/refused rate is filled with zero. Both gauges and both potential policies remain visible.

## Reproduction

All native work stays on Sirius under `/mnt/codex-data/codex/rya1354`, sequential, single thread and nice19. Each native process has an 1800-second limit and an owned process group; a failure stops the batch and retains every output. The fresh immutable snapshot is `evidence/fac_radial_rate_driver_v41`; fresh native output is `runs/fac_radial_rates_v41`.

```sh
# Sirius; the v41 directory already exists after execution. Never overwrite it.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 \
  /usr/bin/time -v -o evidence/fac_radial_rate_driver_v41/batch.time \
  python3 evidence/fac_radial_rate_driver_v41/rya1354_fac_radial_rates.py \
  --root . --source toolchains/fac-9ab6def \
  --index evidence/fac_radial_rate_driver_v41/index.json \
  --index-sha256 8ddec6ba6a2939402fad68f0f43fc39b77da3adaccd2052e1f1aea8f875105f0 \
  --out runs/fac_radial_rates_v41
```

Initial driver tests: five passed in 0.59 seconds, one existing archive test deselected. This checks the original protocol, explicit gauge/mode/cutoffs, dynamic additional lower-state coverage, source-case routing and rejection of altered physics. Prior v40 recovery independently passed all 236 tests in 224.49 seconds.

## Native result

All eight cases completed with exit zero, without failure, retry, omission or changed criterion. Each reproduces the original qualified CI/radial report and C/CM/source.pot/ci.pot bytes exactly. All 104 raw channels pass independent accounting: 90 printed rows plus 14 explicit source-threshold omissions. The tiny omitted E2/M3 contributions remain in the totals.

Exploratory model lifetimes in ns, retaining both gauges:

| Potential / numerical scale | Babushkin | Coulomb |
| --- | ---: | ---: |
| Ground / 10 | 1.5000537250850101 | 1.8839005393305581 |
| Ground / 100 | 1.5000543436504505 | 1.8839000115269429 |
| Equal reference / 10 | 1.5090455139790722 | 3.061008528480786 |
| Equal reference / 100 | 1.50904681885607 | 3.061008630358369 |

At scale100 the previous n3 full-valence B/C lifetimes were ground 1.4925550494452906 / 1.8666878448972763 and equal reference 1.478511960433984 / 2.9029565769155723. Adding 4s/4p changes the model rates, while tightening numerical controls has much smaller effects. Significant gauge and potential sensitivity persists. This is a limited orbital-space sensitivity experiment, not full n4/core/radial-grid convergence, an experimental lifetime determination or a calibrated error bar. No gauge or potential is selected by closeness to experiment.

The complete driver, including repeated index qualification and all native/replay steps, takes 541.23 seconds and 2499972 KiB peak RSS, no swaps. Full science/governance regression: 239 passed in 207.77 seconds. The additional actual-v41 archive test is performed after transfer; it was not counted in that run.

## Evidence and remaining work

The lossless archive preserves 781 paths in 747 parts totaling 1069988168 compressed bytes; largest part 3984482 bytes. Packing and the existing independent decoder finish before transfer, exit zero, in 183.20 seconds and 5623844 KiB peak RSS, no swaps. Complete basis/CM/C, potentials, spinors, overlaps, raw/formatted rates, original full index, scripts, source documents, logs and receipts are retained. The separate resource record is `data/audit/rya1354/source_curation_v2/fac_radial_rates_packing_v41.time.txt`.

The portable reader `scripts/rya1354_fac_radial_rates_audit.py` checks the original archive/fragment/receipt chain, frozen-index and embedded full-report hashes, case-specific original source-manifest identities, original constants, zero-correction inventories and fresh independent raw/formatted-rate replay. The machine twin is `data/audit/rya1354/source_curation_v2/fac_radial_rates_v41.json`. Native coefficient arithmetic remains on Sirius; portable validation checks evidence and rate accounting.

The first portable replay fails before results: the byte-prefix optimization cuts immediately before the words `Mixing Coefficients`, leaving the native heading's leading `============` as an invalid basis row. The original reader remains inside the immutable native archive. The failure record is `source_curation_v2/fac_radial_rates_portable_failure_v41.json`; correction declaration is Linear comment `647ac90b-3f7d-4fc7-9dbd-426a04eb2319`. The corrected slice includes the recognized marker words so the unchanged parser's original section-boundary rule executes. Decorated-heading and missing-boundary regressions pass. No native rerun, archive mutation or numerical/physical criterion change. One characterized portable-reader row is appended to the run-bug ledger.

Corrected portable replay succeeds and reproduces all eight native totals. Final focused command `python -m pytest -q tests/test_rya1354_fac_radial_rates.py` passes five tests in 14.26 seconds, including fresh complete-archive reconstruction and exact machine-twin equality. Combined with the prior 239-test regression run, this validates 241 unique science/governance tests; three focused tests overlap the earlier run. The whole-ledger test separately gives one pass and one existing failure on a historical missing `science_impact` field. The same unchanged schema test fails on current-main ledger bytes; the new v41 row alone passes it. No unrelated historical ledger rows are repaired or relabelled here.

The fetched main at milestone review is `4cdc449a5643987a3459352f99be1926392ed61f`; the task remains tied to its original scientific input base. No state surface changes, grade/product/reference/catalog/star/source-spectrum edits or scientific publication. Register freshness and version-pointer checks pass. The existing bibliography vocabulary failure and unrelated main-branch SEQUENCE drift remain disclosed; no hosted-CI or clean-library success claim.

Resume with a separately declared physical-convergence investigation of the residual potential/gauge sensitivity. Preserve both model policies and gauges, qualify independent radial/grid/core/Hamiltonian axes and state identities, and retain all failures. Existing primary experimental provenance/dependency review, calibrated uncertainties/covariances and independent scientific review remain prerequisites for theoretical grading or adoption. Do not overwrite or repeat v41 merely to select a favorable lifetime. No merge or Done.
