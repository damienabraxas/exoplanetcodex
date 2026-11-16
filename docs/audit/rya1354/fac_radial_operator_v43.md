# RYA-1354 — isolated radial E1 sensitivity, v42–v43

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. This investigation diagnoses existing model sensitivity; it introduces no new SCF, configuration interaction, Hamiltonian, correlated rates, lifetime or theoretical grade.

Continue from `43b5edbb65328c416f6eb99148a19a0aeea0f3c1` on the existing RYA-1354 research branch and draft PR616. Main refreshed at `4cdc449a5643987a3459352f99be1926392ed61f`; worktree `/Users/ryanschmitt/codex/rya1354`. No newer mission notes changed the scope. The v41 evidence and native calculations are not overwritten or repeated.

## What is isolated

Use the four existing orbital models and both gauges in separate cold processes. Each process queries every one of the 42 ordered valence E1 pairs from the ten active relativistic orbitals; closed 1s/2s are excluded by the frozen configuration occupations. Compare three photon-energy contexts: the model's native 3S1-to-lowest-3P2 transition energy; its paired-potential counterpart at the same numerical control; and each pair's own native bound-orbital energy gap. The reference levels are re-derived from full source mappings, rather than assumed by ID. All energies are theoretical and exposed. No observed energy is substituted into any accepted gf, A or lifetime.

The isolated command calls the unmodified pinned FAC `MultipoleRadialFR`. Photon energy is in Hartree and the native photon coordinate is `aw=alpha*E`. A supported nonrelativistic **length** integral is a labelled comparator. Pinned `MultipoleRadialNR` explicitly does not implement velocity; its Coulomb branch still uses length. It cannot be used to claim gauge agreement.

Primary `radial.c` cache keys omit gauge and photon context. The diagnostic therefore enforces one gauge/request context per cold process, each ordered pair exactly once, and an exact one-point native photon grid. `FR0` is not mixed with the production cache representation. All signed amplitudes are retained; no purity, strength or confidence cutoff is introduced.

The new SFAC command is observability only, built in a fresh task-owned frontend directory against the byte-unchanged original library. Original/new executable controls require every physical radial file to be byte-identical. Every probe must reproduce the original radial preflight, original bound-spinor headers/r/rVT0/P/Q, and all protected source/header/spatial-grid/central-potential bytes under the already-qualified W/dW/dW2 scratch rule. Full post-probe wavefunctions and saved potentials are retained.

## Retained v42 failure and numerical-domain correction

The v42 index SHA256 is `11b79832cc358e2fa68185531501acd4ddade9f42e81c33022d3b8530d38a1d6`; its parent receipt/index hashes match the committed v41 manifest. Metadata preparation takes14.40s/360180KiB. Three local control tests pass0.06s. Prospective declaration: Linear `7366df25-82d9-49d7-a992-9b8969d773cb`; plan `lifetime_plan.v42.json`.

The new frontend builds and original/new canaries pass. Ground/coherent10/B own and paired photon contexts both pass all42 pairs and physical invariants. The first orbital-gap context exits1 at inputline71, ordered pair 3d(kappa2) to 4p(kappa1). The exact-grid guard correctly rejects a native clamp. The batch stops after two successful contexts; no later case runs. Whole v42 status stays **FAILED_RETAINED**, elapsed18.28s/802396KiB. Full532 native paths, source, controls, input and logs are sealed in eight lossless parts totaling8126295 compressed bytes, with independent decoder exit0.

Pinned `radial.c` initializes `_awmin=1e-4` and `SetAWGrid` raises a requested minimum below it. The exposed `SetOption('radial:awmin',dp)` changes only that photon scalar. Original printed spinor energies give minimum requested `aw` about6.2972e-5 to6.9952e-5 across the four models. The revised v43 protocol prospectively chooses1e-8, after the unchanged orbital preflight and before queries. All actual requested points are above the new minimum; no query evaluates at the chosen minimum. This is an input-domain control, not a spatial-grid or physical-model change, fitted result, clamped-value acceptance or relaxed tolerance.

The revised index SHA256 is `bed1bf01e3750618884893818d8c4afc018cd1a1f70bdda6684b2a2a5c2b7315`; it links failed receipt `3e156ee352606723b488b070ebbdec640b91b4ab1518bbc0e1eeb93fcf6c58b5`. Preparation takes12.49s/360216KiB. Plan `lifetime_plan.v43.json`, prospective declaration Linear `ee78fb43-4904-44e5-ab30-0703bc435b5c`. New snapshot `evidence/fac_radial_operator_driver_v43`, frontend `toolchains/fac-radial-operator-probe-v43`, output `runs/fac_radial_operator_v43`. C observability and mathematical-library code are unchanged. All exact coverage/frequency/finite/sign/wave/potential criteria remain. Both earlier successful v42 raw radial files and full reports must reproduce exactly under the new control; old success is not counted as independent new physics.

## Reproduction and interpretation

```sh
# Sirius under /mnt/codex-data/codex/rya1354; fresh declared directories only.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 \
  /usr/bin/time -v -o evidence/fac_radial_operator_driver_v43/batch.time \
  python3 evidence/fac_radial_operator_driver_v43/rya1354_fac_radial_operator.py \
  --root . --source toolchains/fac-9ab6def \
  --index evidence/fac_radial_operator_driver_v43/index.json \
  --index-sha256 bed1bf01e3750618884893818d8c4afc018cd1a1f70bdda6684b2a2a5c2b7315 \
  --probe toolchains/fac-radial-operator-probe-v43 --out runs/fac_radial_operator_v43
```

Sequential, single-thread and nice19; 180s build and300s per owned native process group. Stop/preserve the first failure. Do not overwrite these directories. Whole new native files, source snapshots, build diff/logs, controls, input index and receipts must survive lossless packing and independent decoding before transfer.

The audit separates the pinned rate formula into photon-energy power `E**(2|m|+1)` and squared signed-amplitude factors. It joins nominal physical channels by configuration/parent/term, J, parity and multipole, not wavelength or ID alone, retaining all original mixtures in linked evidence. This is not a cross-potential adiabatic overlap proof. Common-energy radial comparisons isolate orbital-model dependence without fabricating a hybrid adopted rate. Normalized radial gauge ratios are diagnostics, not sigma/error bars or a validated many-body uncertainty model.

Four local control/factorization tests pass0.05s; the earlier broad run passes244tests in221.32s before the revised protocol. Final actual-evidence verification is recorded after native completion and transfer. Existing bibliography vocabulary and unrelated SEQUENCE/whole-ledger debt remain disclosed. No clean-library or hosted-CI success claim. Canonical grades/products/reference/catalog/star/source-spectrum surfaces remain unchanged; no publication, merge or Done.

## Completed diagnostic and sensitivity result

All24 revised cold contexts pass, yielding1008 finite signed radial samples with exact requested photon points. Original/new frontend physical canaries pass; both successful v42 raw files and full metadata reports reproduce exactly. All original/post-probe physical spinors and protected potential fields remain unchanged. Native batch exit0,31.57s/802380KiB, no swaps. Independent portable replay reproduces every case report and checks the full archive, source/parent/index/failure chain, unchanged earlier controls and complete42-pair contexts. No v42 failure is relabelled passed.

For the reference 3s-to-2p pairs at coherent100, the **magnitudes** show:

| Potential | Final orbital | Own-transition `abs(C/B)` | Orbital gap / photon energy | At orbital-gap photon energy `abs(C/B)` |
| --- | --- | ---: | ---: | ---: |
| Ground | 2p1/2 | 0.8493067685 | 0.8493042919 | 1.0000023829 |
| Ground | 2p3/2 | 0.8455783927 | 0.8455495460 | 1.0000262161 |
| Equal reference | 2p1/2 | 1.4738849960 | 1.4738790891 | 1.0000064642 |
| Equal reference | 2p3/2 | 1.4685955077 | 1.4685632094 | 1.0000391396 |

Signed native values remain in the archive and machine report. These internal radial results are consistent with photon-versus-orbital-energy scaling in these models. Near-unity orbital-gap ratios do not calibrate an atomic uncertainty or validate the correlated transition. The full many-electron reference-line amplitude magnitudes `abs(C/B)` are0.8921336072(ground) and0.7019146236(equal reference), so the isolated reference pair alone does not explain the many-body result. A universal energy rescaling or selection of a favorable gauge would be unjustified.

For the full reference 3P2 E1 line, changing ground to equal-reference potential gives the same `E**3` factor0.4803766705 in both gauges. Its squared-amplitude factor is2.0688038492 in B and1.2806430503 in C. Products0.9938051049 and0.6151910445 reproduce the native rate ratios within floating arithmetic. Thus nearly equal B rates involve compensating energy and amplitude changes, and cannot establish model convergence. All13 channels and tiny contributions remain accounted; the three reference 3P E1 branches carry more than0.999997 of each model total in these cases. No other channel is dropped from the audit.

The successful archive retains1660 paths in41 parts totaling42519565 compressed bytes. Full packing plus independent decoder takes4.82s/188124KiB, exit0, no swaps. All parts are<=5MiB. The separate log is `source_curation_v2/fac_radial_operator_packing_v43.time.txt`; the machine report is `source_curation_v2/fac_radial_operator_v43.json`; reader `scripts/rya1354_fac_radial_operator_audit.py`.

Final validation: **246 tests pass in215.15s**, using the existing Mac Python3.12.13 for source/metadata/report checks only:

```sh
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q \
  tests/test_rya1354*.py tests/test_atomic_benchmark_rya1354.py \
  tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py \
  tests/test_line_match_rya1033.py
```

This includes actual failed/successful archive and machine-twin replay, original successful-control equality, complete query/context checks and source-formula sensitivity factorization. Both new ledger rows pass the unchanged schema test when scoped to those rows; historical whole-ledger missing-field/vocabulary failures are not rewritten. Final original FAC source diff is empty, original library/executable hashes are unchanged, no owned native calculation/packer remains, and protected canonical science surfaces compare unchanged to original task base `adefcaf073ff70035954cf45a9a2283774ac7d06`.

The research-branch push retains the previously documented `--no-verify` draft exception for unrelated merged SEQUENCE debt, after explicit blob-size, whitespace, protected-surface and register-freshness checks. This does not authorize a merge or scientific adoption.

## Additional pinned-source concern

Inspection of pinned FAC `faclib/transition.c:TRMultipole` finds `wm` declared without initialization on the explicit nonzero-multipole path, then read by `m == -1 && wa > wm` before its first assignment. The separate m-zero branch initializes it. This characterizes an unqualified optional dominant-pair selector (`pk0/pk1`), not a quantified physical-rate defect. C undefined-behavior safety for the whole old binary is not established by repeatability. Existing original raw rates remain unadopted/HOLD and do not use returned pair IDs; this radial diagnostic never calls `TRMultipole` and enumerates all quantum pairs directly. The original library/source is intentionally unchanged. A separate characterized source finding is recorded in the run-bug ledger, with no invented runtime fix/test or physical-error claim.

## Next bounded investigation

Freeze a full signed transition-density/component audit for the existing CI wavefunctions: retain every angular coefficient and ordered radial-pair contribution, reconstruct the complete native E1 matrix element in both gauges and its actual photon-grid context, and quantify cancellation/mixing before deciding which correlation/core/radial axis to expand. Do not infer full matrix-element dominance from a reference configuration's probability weight or optional native selector. Resolve the characterized selector concern in a separately qualified observational or numerical control, preserving the pinned original. Continue independent physical-convergence and primary-experimental provenance/covariance review before any theoretical grade or adoption.
