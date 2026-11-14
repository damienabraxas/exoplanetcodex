# Expanded FAC angular mapping and complete rates

The n3 configuration experiment now reaches complete exploratory lifetimes in both gauges. It exposes substantial correlation-space and gauge sensitivity; it does not establish physical convergence or calibrated accuracy. EXPOSED_NONBLIND/scientific HOLD remain, with accepted lifetime, uncertainty, covariance and calibrated sigma null. No grade or atomic-product adoption is justified.

## Full expanded LS replay

Declaration bffa7994-2472-41e0-8bc3-c56c96d9d50a preceded unmodified pinned GRASP jj2lsj on exact copies of all four v23 C/CM sets. CI input, all levels/full composition/C-J-T/NK convention, unique labeling OFF. Z8 support metadata does not replace the FAC nucleus or Hamiltonian. Original C/CM remain byte-identical; converter source/executable receipts and prompts retained. The angular implementation is shared with GRASP, not independent.

A separate Sirius-only auditor reconstructs every full LS vector from original FAC coefficients and saved transformation matrices. It handles wrapped full-vector lines and independently encodes the pinned PACKLS occupation/shell-parent/intermediate/final-term rules for this declared s/p/d space. All3896(singles)/63316(SD) coefficients per numerical case match exactly at native eight-decimal print precision. Every original ID/J/parity/native energy, dominant full configuration/parent/term and printed component weight is cross-checked. Transformation orthonormality residuals9.99e-16/1.55e-15 are representation diagnostics, not uncertainties. The formatter also reproduces the earlier minimal analytic control's complete LS names.

The targets are uniquely identified as 1s2 2s2 2p3(4S)3s3S1 and5S2, original IDs6/5 in all four cases. Reference-component weights are approximately.96808/.97104(singles) and.96545/.96882(SD); no purity cutoff, energy-nearness choice or admixture removal. All full compositions remain in native evidence. Sirius canaries change a full-vector sign, transformation-matrix sign, original ID and CSF occupation; all four are rejected by replay. No native file is altered by a canary.

## Rates and completeness

Declaration c406127a-7077-4d11-89fe-981125098d88 preceded eight cold native runs: two spaces x two numerical controls x separate Babushkin/Coulomb processes, explicit fully relativistic mode0 and zero table cutoffs. The unchanged v21 diagnostic frontend is reused only after binary, static library and original SFAC hashes match. No new frontend/library/source patch. Full LS replay/source qualification runs before rates; each fresh rate process passes numerical stopping and reproduces mapped C/CM bytes exactly.

All lower states and all parity/total-J-allowed ranks are enumerated from the newly verified binary-energy mapping. These cases independently yield lower0..5 and13 channels again; neither the old ID list nor channel count is assumed. A synthetic added lowerJ4 produces E3/M4/E5, verifying that the request builder expands coverage. No unsupported missing-channel zero filling is allowed.

All13 native channels now return finite rates. In particular E3 to lower0/3 is NONZERO in both expanded spaces, confirming the old minimal-basis zero proof is inapplicable. The rates are tiny (about1e-13 to8e-11/s), but explicitly retained. The only printed omission is tiny nonzero M3, whose group rate is below the pinned trd<1e-30 atomic-unit writer gate; it is retained from raw output. All12 printed rows per case match independent source-formula replay exactly at %.6E precision after native float32 amplitude quantization. Raw gf/A formulas, original state J/energies, gauge/mode/electron count, complete request uniqueness, table block counts and source receipt links are checked. Transition-energy comparison permits only the ulps introduced by reconstructing the binary common energy offset; no physical-energy tolerance is fitted.

| Space / numerical control | Babushkin lifetime (ns) | Coulomb lifetime (ns) |
|---|---:|---:|
| Minimal / coherent100 (prior) | 2.4121030293305274 | 2.316978563137304 |
| Singles / coherent10 | 4.0102352905597325 | 2.452229290977717 |
| Singles / coherent100 | 4.010241495392219 | 2.452229080616681 |
| Singles+doubles / coherent10 | 3.637387506162713 | 2.5013081094642646 |
| Singles+doubles / coherent100 | 3.6373937105892775 | 2.501307849682168 |

Tighter numerical controls are stable within each model, while adding correlation materially changes lifetimes and the gauges remain separated. These differences are not error bars. No gauge is selected against the exposed lifetime experiment, no experimental energy is substituted and no closer-to-reference choice is made. Ground-only central-potential optimization, fixed2s occupation, restricted n3 valence orbitals and nuclear/optimization differences from correlated GRASP remain explicit limitations. This is not a matched-Hamiltonian/correlation/radial convergence test against GRASP.

## Reproduction and continuation

Sirius root /mnt/codex-data/codex/rya1354. Completed runs/fac_expanded_recoupling_v24 and runs/fac_expanded_rates_v25; exact drivers evidence/fac_expanded_angular_driver_v24, fac_expanded_replay_v24, fac_expanded_rate_driver_v25 and fac_expanded_canaries_v24. No solver jobs remain. Do not overwrite/relaunch completed cases. Angular archive126paths/4810958bytes; rate+canary archive181paths/4223262bytes; each one lossless part below5MiB. Source/runtime pins and original executable remain unchanged; primary label formatter source retained separately in angular archive. Machine reports source_curation_v2/fac_expanded_ls_v24.json (compact index of full native report) and fac_expanded_rates_v25.json.

Metadata/rate replay: `python scripts/rya1354_fac_expanded_rate_audit.py`. Full numerical angular replay must run on Sirius with the snapshotted rya1354_fac_expanded_ls_audit.py and explicit --source runs/fac_space_v23 --converted runs/fac_expanded_recoupling_v24 --out <new-report>. Heavy conversion/matrix replay does not run on the Mac.

Next audit FAC's potential-optimization weighting/defaults and design a separately declared fixed-SD-space control of ground-only versus balanced-reference optimization. This can distinguish an orbital/potential sensitivity from adding configurations; do not change both silently or tune to an experimental lifetime. Requalify numerical stopping, full LS identities and all lower channels for every new case. Further core/radial/active-space convergence, independently reviewed experimental provenance/covariance and scientific review remain owed before grading/adoption or known-grade scoring. No merge/publication/Done.

Exact native full-replay command (use a new output filename):

```sh
cd /mnt/codex-data/codex/rya1354
PYTHONPATH=evidence/fac_expanded_angular_driver_v24 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 python3 evidence/fac_expanded_replay_v24/rya1354_fac_expanded_ls_audit.py --source runs/fac_space_v23 --converted runs/fac_expanded_recoupling_v24 --out /tmp/rya1354-expanded-ls-repeat.json
```

No new upstream defect was discovered; the wrapped output and richer PACKLS format were audited before scientific use. Native sign/identity/occupation canaries and raw/table corruption tests reject invalid artifacts. No new bug-ledger row is warranted; all earlier failures and characterized defects remain retained. Initial driver transfer hit the network sandbox, was retried with escalation, and succeeded; the intervening missing-script invocation created no native result directory. No model or criteria were changed in response.

Validation: `python -m pytest -q tests/test_rya1354*.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py` passed180 tests in105.56s. Final12 new LS/rate metadata tests pass again after explicit cross-parent receipt/frontend provenance links were added to the rate auditor. Sirius four full-matrix canaries all reject. Whitespace, protected pipeline/products/linelists/reference/catalog/stars comparison against task baseadefcaf073ff70035954cf45a9a2283774ac7d06 and manual new-file<=5MiB guards pass. Existing bibliography civis2013/mallinson2024 verification-state vocabulary still stops the library audit before traversal; no library-success or hosted-CI-pass claim. Starting commitd8f3725910a82dde921937ffd98c4da9d489e4a0; fetched main45901e7678ab6c8624fe659e552fba3f2681457d. All owned native jobs complete, original source checks clean; existing draft push exception for unrelated merged1352/1353 SEQUENCE debt retained with manual size guard.
