# Oxygen emission and diffusion-boundary profile axis

RYA-1371 continuation from c603edc55ee74df23bfab8bd55b40e7afb7ae063.
Worktree `codex/rya1371`, branch `ryandamienschmitt/rya-1371-oxygen-reference-provider`.
Refreshed origin/main remains 0523151e. Draft PR620, scientific HOLD.

## Predeclared bounded contract

The next fixed-population diagnostic selects either held-LW or source-native
Planck emission for BOTH background thermal absorption/emission and bottom
thermal diffusion inputs. Both alternatives must be executed and retained;
neither is selected by smaller residual. Background chi/scattering, colliders,
active BB/BF coefficients and rates, imported n/nStar, geometry, wavelengths,
upper ghost, resource limits and tolerances remain unchanged. The existing
native-Hermite alpha control is retained in every comparison.

`EmissionProfileBoundary` uses a Python callable with the exact linear-solver
bottom condition inspected in held `FormalScalar.cpp`: B_bottom minus the
upstream-to-bottom B difference divided by trapezoidal optical depth along the
ray. Total bottom opacity is captured from actual depthData after the first
iteration and required invariant thereafter. First capture starts at B_bottom;
convergence includes the boundary refresh and may not stop at that initialization.
Wavelength/state extrapolation and negative/invalid coefficients/intensities
are refused, never clipped. Backend/default files and constants are unchanged.
The source inspection is distinct from a source-to-installed-binary attestation;
a disposable continuum context checks the callable against the actual compiled
ThermalisedRadiation boundary, including gradient and angular dependence.

Before the held source case: test invalid contracts/arrays and compiled-boundary
reproduction. Then fresh default reproduction against prior native-alpha replay,
followed by BOTH profiles at the exact inherited 19,963-node state. Preserve all
failed attempts, input/helper snapshots, source/build receipts and raw outputs.
Use existing one-thread 600s/20,000-node constraints, Mac only, no SE update.
Independently reconstruct all source radiative rates from each new saved J.
Report total-opacity/scattering/population/collision identity separately from
changed background eta, lower intensity, J and rates. No observation, abundance,
canonical-solar, UV/fullband/3D, coupled closure or adoption claim follows.

## Intake and verification

Repository task-intake, implementation-brief, scientific-code-review and
state-register procedures apply. Existing Solar Beta refinement report has
10 ticketed unresolved entries and five without resolving tickets; none is
cleared by this diagnostic. O model admission remains open. No new library
acquisition or physics source needed: held source/backend/atom receipts own
this control. Full source and frozen products remain protected.

Initial disposable test attempt failed before Context construction because a
passive-only RadiativeSet cannot generate its wavelength configuration. The
fixture now generates the configuration while H is detailed-static, then marks
H passive before construction (the existing frozen-total-provider pattern).
No native/solar run occurred on that fixture failure. No pipeline defect or
source-data change discovered at this point.

## Executed controls and independently checked results

All three fresh RT runs completed under the original resources, each in nine
radiation iterations; no population solve, native rebuild or atomic-rate edit.
The new default exactly reproduces all eight `rates.npz` arrays from the prior
`rya1371-native-alpha-joint-v1` replay.

| Case | RT seconds (excludes preflight) | Maximum physical balance residual | Max relative J change from default |
|---|---:|---:|---:|
| Existing default | 10.1852 | 2.00412354504897e-4 | 0 |
| Both emission/bottom inputs held-LW | 7.13533 | 5.59305217018223e-5 | 4.694458664e-4 |
| Both emission/bottom inputs source-native | 6.14404 | 2.00412354504897e-4 | 9.279669581e-5 |

The source-native case changes bottom diffusion input while retaining native
background emission; its global maximum residual remains unchanged even though
J and some rates change. Maximum positive-rate relative changes are
4.653734424e-4 (held-LW) and 9.211262425e-5 (source-native). Background chi/sca
and captured TOTAL bottom chi are exactly identical between both controls;
n/nStar/C, native reference R/C and wavelengths equal default exactly.
This is a controlled response to declared emission conventions, not a reason
to select a preferred profile or an independent physical-error estimate.

Each fresh saved J independently reconstructs all 256 BB / 22 BF transitions
and 31,136 positive directional-depth rates. Largest relative reconstruction
errors are 2.97540e-14 / 2.93099e-14, below the unchanged positive-summation bound
3.56799e-10. Reconstructed residuals are 5.59305217030359e-5 and
2.00412354506986e-4. Imported native LTE, active LW Einstein/inverse factors,
partition/Saha/EOS and collision conventions remain mixed and unqualified;
there is no complete thermal closure or coupled SE admission.

Evidence: `data/results/rya1371/emission_profile_v1` receipts 705 raw files,
38 compact exact copies under the unchanged 5 MiB cap. Raw arrays/oracle tables
remain external with receipts. Comparison v1 is retained; v2 adds explicitly
after-run Python component snapshots, checked against actual producing input
receipts. Comparison v1's source was not separately snapshotted before that
packaging expansion; its scalar comparison is reproduced by v2. Actual RT and
independent audit source snapshots are retained unchanged. No failed real run.
Normal held-LW configuration warning remains visible; no benchmark or install.

Reproduction (from repository root; fresh outputs required):

```sh
source config/hosts/mac.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Use the held LW interpreter resolved from engines.bride_eval/.venv/bin/python.
# RAW is the resolved experiments output root used by these immutable receipts.
LW="$CODEX_BRIDE_EVAL/.venv/bin/python"
RAW="$CODEX_WORK_ROOT/bride_mac_checks"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_grid_refined2_native_alpha.json --output "$RAW/rya1371-emission-default-v1"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_emission_lw.json --output "$RAW/rya1371-emission-lw-v1"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_emission_native.json --output "$RAW/rya1371-emission-native-v1"
"$LW" scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_emission_lw.json --replay "$RAW/rya1371-emission-lw-v1" --output "$RAW/rya1371-emission-lw-audit-v1" --absolute-replay --profile-contract
"$LW" scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_emission_native.json --replay "$RAW/rya1371-emission-native-v1" --output "$RAW/rya1371-emission-native-audit-v1" --absolute-replay --profile-contract
"$LW" scripts/compare_oxygen_emission_profiles.py --prior "$RAW/rya1371-native-alpha-joint-v1" --default "$RAW/rya1371-emission-default-v1" --lw "$RAW/rya1371-emission-lw-v1" --native "$RAW/rya1371-emission-native-v1" --lw-audit "$RAW/rya1371-emission-lw-audit-v1" --native-audit "$RAW/rya1371-emission-native-audit-v1" --output "$RAW/rya1371-emission-comparison-v2"
"$LW" scripts/archive_oxygen_reference_checks.py --raw-root "$RAW" --output data/results/rya1371/emission_profile_v1 --emission-profile-supplement
"$LW" -m pytest -q tests/test_*rya1371.py
```

Validation: 127 issue-focused tests pass in 3.03s. Disposable initial test-fixture
failure is described above; all actual RT/audit/comparison/archive commands exit0.
Library check exits1 with the same 83 unindexed documents (215 total,132 indexed,
12 bibliography paths absent on Mac). No bibliography facts invented or new
materials acquired. Existing Solar Beta debt is unchanged. Run-bug ledger:
none discovered; characterized thermal mismatch remains scientific open work.

Self-review: scoped diagnostic only, default reproduction exact, complete graph
retained, no fallback/clipping, immutable source receipts and independent rate
reconstruction checked. Coupled runner refuses the new plan key; production
adapters/configuration/products untouched. Lower-boundary implementation is
restricted to stationary fixed populations and the linear solver; extension to
moving atmospheres or coupled SE is not qualified.

Next: qualify active BB/BF inverse and emission relations together with imported
LTE/energy/partition/Saha and reverse-collision ownership under an explicit
profile. Do not fold only modified rates into a coupled solver. Numerical
convergence, joint operator/restart, canonical-solar observations, Kitt Peak AIR,
intensity/flux geometry, continuum/tellurics and instrument/fullband qualification
remain open. Fe-next resource/FeII/atom gates and C/N/RYA-1374 remain unchanged.

Verification close: all 17,420 raw receipts / 1,071 compact copies across RYA-1371
verify by length/hash/byte identity. Staged generator error set exactly matches
the existing 5,120 diagnostics; the initial unstaged check saw the 39 new files
as untracked, and the staged rerun resolves those without changing any guard.
Mac-profile path audit: 244 present / 18 pre-existing broken, unchanged. An
initial path check without sourcing mac.env reported root absence rather than
host availability; the sourced check is the applicable evidence. Conflicting
archive supplement flags correctly refuse before creating output. Register
v294 pointer, params mirror, staged whitespace, reconciliation and freshness
pass; pre-existing tracker/gold regeneration notice remains. No global clean
repository claim; immutable prior native-log whitespace remains preserved.
