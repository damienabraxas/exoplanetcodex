# Joint fixed-population oxygen LTE/inverse operator diagnostic

RYA-1371, continuation fromec57c804; draftPR620, scientific HOLD.
Task-intake/implementation-brief/scientific-review/state-register procedures
continue from this session. Latest Linear handoff unchanged; main4cdc449a.

## Predeclared contract

Retain prescribed native actual populations and total, T/ne/rho/colliders,
source atom/forward coefficients/widths, all19963wavelengths/256BB/22BF, native
Hermite alpha, one-thread600s/20000point controls and source upper ghost.
Baseline: existing held-LW background/bottom emission profile with source-native
LTE/C. Compare an explicitly experimental held-LW finite-level LTE shape
(absolute total retained) with channel-specific source inverse collisions.
Actual held library and independent log formula must reproduce at1e-12,
following the v295 arithmetic-only screen. No EOS/partition or model adoption.

**Source directions matter.** CE and CH preserve downward source coefficients;
CI and CH0 preserve upward source coefficients. Rebuild inverses per original
source record before summing, with exact per-kind forward identity. The v295
all-upward algebraic canary is not used as solver input. No rates invented for
level41, no source forward-rate/energy/constant adjustment, no residual tuning.

The compiled backend consumes one declared nStar for active BF stimulated
opacity/emissivity AND radiative integrals. BB Einstein/CRD formulas remain
unchanged; finite-profile detailed-balance limits remain explicit. Background
emission and lower diffusion inputs both retain held-LW. Context actual n is
fixed throughout; modified reference nStar does not update actual populations.
The coupled runner refuses the control until separately qualified.

Test mixed CE/CI/CH/CH0 overlapping edges, preserving per-channel source direction
and distinguishing unchanged-C versus refreshed-inverse reciprocity. Refuse
partial/mismatched profile controls and changed source state. First fresh
baseline must reproduce every earlier held-LW emission replay rate/state/field
array exactly. Then run joint fixed-population RT under identical controls.
Independently reconstruct all radiative rates from fresh saved J AND total
stationary chi/thermal eta across all nodes/depths. Coefficient arithmetic screen
256eps*(1+278) is predeclared from a conservative sum-operation count; it is not
physical uncertainty. Scattering is retained and counted once. Preserve all
attempts, producing source snapshots and immutable receipts under5MiBcap.

This milestone can qualify a fixed-state coefficient/rate diagnostic only.
Joint SE/preconditioning/restart, numerical refinement, absolute EOS/partition,
observed canonical-solar/intensity-flux/AIR/continuum/telluric/fullband/3D remain
open; Fe-next atom/resource/FeII and C/N/1374 gates unchanged.

## Executed fixed-state milestone (v296, 2026-11-16)

Fresh held-LW emission baseline reproduces all8prior `emission-lw-v1/rates.npz`
arrays exactly. Joint control keeps actual native populations fixed while the
compiled active-opacity/emissivity AND radiative-rate machinery consume the
new finite-level LTE reference together. Source record directions preserved
exactly for741CE/39CI/302CH/24CH0; each inverse regenerated before aggregation,
never by copying v295's all-upward algebraic projection. Colliders, source q,
energies and source CI exponential remain original. Thus this is coordinated
inverse-input sensitivity, not a newly derived microscopic collision model.

| Actual case | Wavelengths | Reported seconds excluding preflight | Radiation iterations | Largest physical balance residual |
|---|---:|---:|---:|---:|
| Original LTE/C, fresh emission baseline | 19963 | 6.79505 | 9 | 5.59305217018223e-5 (level41,bottom55) |
| Joint held-LW LTE / source-direction inverses | 19963 | 6.61454 | 9 | 2.72286197895160e-4 (level41,top0) |
| Same joint control, inherited midpoint grid | 17263 | 5.89668 | 9 | 2.81435582864789e-4 (level41,depth1) |

The joint convention WORSENS agreement with the native stationary populations;
it was not selected by smaller residual and is not adopted. Differences from
baseline: nStar6.12918e-4, positiveR5.45321e-4, positiveC4.98148e-4,
J9.45491e-6, totalchi2.23744e-9, totalthermaleta9.75790e-6 maximum relative.
Background chi/eta/scattering are exactly identical; actual n, native reference
R/C and fine wavelengths remain exact. Reference total oxygen is conserved to
the existing1e-12 population arithmetic screen; absolute EOS normalization is
not established by that comparison.

All43680collision pair-depths balance under declared nStar to2.45635e-16,
no one-way edges. Level41 remains collision-disconnected. Whole-graph net-flow
cancellation over41levels/56depths is below3.85e-16 of gross flow across all3
cases. This is arithmetic particle accounting, not SE stationarity or a coupled
population-conservation experiment.

### Independent opacity/emission and rate evidence

For EACH actual case, independent reconstruction of total chi and thermal eta
on every wavelength/depth cell agrees4.44089e-16 relative, below predeclared
1.58593e-11 arithmetic screen. Background scattering is retained and backend
adds sca*J exactly once; reconstructed eta is thermal only. All256BB/22BF source
transitions included. Positive directional/depth rates31136 reconstruct to
2.97540e-14 baseline /3.07532e-14 joint /3.41949e-14 coarse, below unchanged
positive-summation screens3.56799e-10 fine /3.11751e-10 coarse.

All17263coarse nodes occur exactly in19963fine grid; their J arrays are
bit-identical. n/nStar/C/sourceR/sourceC are identical between joint grids.
Maximum positiveR refinement response3.75128e-4; numerical convergence remains
OPEN. Finite-profile CRD thermal defect under held-LW Planck/actual declared
nStar is1.27800e-7 fine versus3.87732e-7 coarse, while BF consistent-Planck
canaries are below2.11e-16. BB line-centre closure does not make finite CRD
frequency integration exact. Native-Planck counter-canaries remain separate;
neither profile or mesh selected by closeness to reference.

### Reproduction and provenance

```sh
source config/hosts/mac.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LW="$CODEX_BRIDE_EVAL/.venv/bin/python"
RAW="$CODEX_WORK_ROOT/bride_mac_checks"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_emission_lw.json --output "$RAW/rya1371-joint-lte-baseline-v1"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_joint_lte.json --output "$RAW/rya1371-joint-lte-replay-v1"
"$LW" scripts/replay_hop_oxygen_radiative_rates.py --plan config/experiments/rya1371_hop_joint_lte_midpoint.json --output "$RAW/rya1371-joint-lte-midpoint-v1"
"$LW" scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_emission_lw.json --replay "$RAW/rya1371-joint-lte-baseline-v1" --output "$RAW/rya1371-joint-lte-baseline-audit-v1" --absolute-replay --profile-contract --coefficient-replay
"$LW" scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_joint_lte.json --replay "$RAW/rya1371-joint-lte-replay-v1" --output "$RAW/rya1371-joint-lte-audit-v1" --absolute-replay --profile-contract --coefficient-replay
"$LW" scripts/audit_oxygen_integrated_thermal_rates.py --plan config/experiments/rya1371_hop_joint_lte_midpoint.json --replay "$RAW/rya1371-joint-lte-midpoint-v1" --output "$RAW/rya1371-joint-lte-midpoint-audit-v1" --absolute-replay --profile-contract --coefficient-replay
"$LW" scripts/compare_oxygen_joint_lte.py --prior "$RAW/rya1371-emission-lw-v1" --baseline "$RAW/rya1371-joint-lte-baseline-v1" --joint "$RAW/rya1371-joint-lte-replay-v1" --coarse "$RAW/rya1371-joint-lte-midpoint-v1" --baseline-audit "$RAW/rya1371-joint-lte-baseline-audit-v1" --joint-audit "$RAW/rya1371-joint-lte-audit-v1" --coarse-audit "$RAW/rya1371-joint-lte-midpoint-audit-v1" --output "$RAW/rya1371-joint-lte-comparison-v1"
"$LW" scripts/archive_oxygen_reference_checks.py --raw-root "$RAW" --output data/results/rya1371/joint_lte_v1 --joint-lte-supplement
"$LW" -m pytest -q tests/test_*rya1371.py
```

Fresh outputs required. Evidence `data/results/rya1371/joint_lte_v1`:853raw
receipts/47compact exact copies under unchanged5MiBcap. Large state/coefficients
and private source-alpha oracle arrays remain external receipt-only. Actual
producer Python source snapshots are explicitly after-run, checked against
producing input fingerprints; primary comparator also preserved. No failed
actual run/fallback. Held-LW default config warning remains visible; no new
install/benchmark/backend build, Sirius/hostedCI or peak-memory measurement.
Original600s/one-thread/20000nodes and native2GBcontrols untouched.

133issue-focused tests pass2.66s, including2new mixed-channel/partial-contract
controls. Initial tests before execution also133pass3.22s. Independent audits
and whole-graph accounting pass; no physical data/abundance tuning or missing
collision repair. Run-bug ledger: none newly discovered this milestone.
Library check:215documents/132indexed/83unindexed/12absent bibliography paths,
unchanged; no fabricated bibliography or new acquisition. Existing path/generator
and tracker/gold debt reported separately. Refinement debt unchanged.

Self-review: compiled shared LTE path, whole source graph and colliders retained;
per-channel original forward directions checked separately from aggregateC;
actual n invariant; background exact; independent total coefficients and rates
verified on both grids; CRD numerical limitation retained. Scientific HOLD.

Next: resolve joint radiative-rate quadrature convergence under unchanged
packet/resource bounds and independently qualify absolute oxygen/EOS/partition
normalization before interpreting this as a physical model. A future coupled
operator/preconditioner needs its own full action, conservation, convergence and
restart checks; current runner explicitly refuses joint_lte_control. Canonical
solar/observed AIR/intensity-flux/continuum/telluric/fullband/3D and Fe-next/CN/1374
gates persist. No default/production/observational product, merge or publication.

Verification close: all18285raw receipts/1128compact copies acrossRYA1371
verify by hash/length/exact bytes. Staged generator diagnostic SET unchanged5120,
zero added/removed; Mac-profile paths244present/18existing broken. Registerv296
pointer/params mirror/reconciler/stagedwhitespace pass. Explicit CLI checks
refuse joint_lte_control in coupled runner and conflicting archive flags BEFORE
creating output. Existing tracker/gold regeneration notice retained; no global
green claim. Worktree codex/rya1371, branch
ryandamienschmitt/rya-1371-oxygen-reference-provider, baseec57c804,
predeclaration3f339371; existing draftPR620. No other worker changes overwritten.
