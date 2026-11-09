# RYA-1354 — numerical stopping and physical-model error investigation

**Exploratory, exposed nonblind work; scientific HOLD.** Start commit
`eacbdad43f8520cbc67f13b68a475570891fa4c0`, original task base
`adefcaf073ff70035954cf45a9a2283774ac7d06`. Same branch/worktree and draft PR 616.
No adopted atomic data or stellar result changes. The frozen choices and subsequent
explicit extensions are recorded in `config/atomic/rya1354/numerical_plan.v3.json`;
source snapshots in the capsules preserve each executed driver.

## Separating orbital stopping from accuracy

The pinned GRASP `scf.f90` deliberately allows either native orbital-change
convergence or relative energy stationarity. This is not an upstream failure to
implement its stated criterion. It is insufficient evidence for a claim that all
varied orbitals satisfy the function criterion. A separate **local diagnostic
variant** removes exactly the one energy-only stop statement. It retains every
other update, object, Hamiltonian, integral and native ACCY value. The pinned
source/build is untouched. The patch, original and modified source, compile/link
commands, executable hashes and clean upstream revision are retained in
`toolchains/grasp-strict-scf-v1` inside the capsule. This variant must not be
misrepresented as the unmodified upstream solver.

A one-cycle run from identical n=3 input gives **byte-identical `rwfn.out` and
`rmix.out`** in the original and strict executables. Both exhaust their deliberately
one-cycle limits, and both are retained as failures. The strict variant therefore
preserves the numerical update before stopping decisions diverge.

| Same-layer continuation | Iterations | Maximum final native SCNSTY | Native ACCY |
|---|---:|---:|---:|
| Ground n=3, damping 0.9 | 39 | 5.54e-8 | 5.96046448e-8 |
| Ground n=4, damping 0.9 | 113 | 5.93e-8 | 5.96046448e-8 |
| Ground n=5, damping 0.5, native NSIC=20 | 22 | 5.82e-8 | 5.96046448e-8 |
| Fresh n=5 verification, damping 0.5, NSIC=0 | 2 | 5.03e-8 | 5.96046448e-8 |
| n=3 half-spacing grid, explicit tighter ACCY | 40 | 8.68e-10 | 9.31322575e-10 |

All completed states above meet their declared native orbital-change criterion. `CONSIS` defines SCNSTY as maximum
radial-function change times sqrt(occupation), **not** an equation residual or a
physical error bar. Damping is bypassed natively for an orbital already below ACCY.
These continuations retain earlier layers frozen; they do not reconstruct a fully
qualified ladder or establish grid/core/physical convergence.

Two slower n=5 damping-0.9 controls were intentionally interrupted after being
superseded, not classified as convergence exhaustion. Their partial artifacts and
original FAILED_RETAINED receipts remain intact, with a separate explicit
interruption disposition. Signals were guarded by exact owned process working
directory. No other worker was stopped.

The first half-spacing attempt wrongly assumed native ACCY was automatically
recomputed after changing H. `GETSCD` initializes it before grid input; the run
retained the coarse criterion while its receipt declared the finer one. It is
rejected for the intended test, and its postprocessor fails closed. The corrected
v2 explicitly revises ACCY and checks it against the native prompt. Native
`rmcdhf.sum` confirms H=0.03125, N=439 and the same outer radius, 1892.317 Bohr.
The four diagnostic M1/E2 rates are unchanged at native printed precision after
this grid/tolerance control; frozen lower orbitals are interpolated, not reoptimized.
The bug ledger and a regression test preserve the criterion-metadata defect.

The resulting four n=3 diagnostic M1/E2 reference rates change by at most about
0.0022%, and n=4 by at most about 0.0097%, from the previous energy-stopped runs.
At n=5 the two M1 changes are below 0.0007%, while the E2 changes are about
-0.0941%/-0.0735%. Thus the large evaluated-reference discrepancies are not removed by these numerical
stopping improvements. All 16 angularly allowed one-photon channels and identical
f1/f2 directional duplicates are checked; duplicates are never summed.

## Forbidden-rate physical controls

All cases below use fixed strict-qualified n=3 orbitals. An identical-space RCI DC
control precedes added physics, isolating RCI interpolation from physical changes.
Native RCI uses its own grid, as documented in the prior sensitivity campaign.
Both recoil contributions are disabled. Stationary-mass divide-by-zero and
underflow/denormal flags remain recorded with the earlier source diagnosis and
synthetic positive-mass/no-recoil control; no new isotope interpretation is made.

| n=3 case | 1D2–3P2 energy (cm^-1) | M1 A (s^-1) | E2 A B (s^-1) |
|---|---:|---:|---:|
| Same-space DC | 17274.78 | 0.00960437 | 0.0000417565 |
| Static-limit transverse Breit | 17264.38 | 0.00706767 | 0.0000304444 |
| One 2s hole, DC | 17151.12 | 0.00945852 | 0.0000404847 |
| Two 2s holes, DC | 17090.12 | 0.00949394 | 0.0000400610 |
| Two 2s holes plus static-limit Breit | 17079.61 | 0.00688663 | 0.0000287854 |

All use valence/core SD through 3d, 1s fixed. Native LS/parent/J/parity assignments
remain checked. The large Breit sensitivity here contrasts with its small effect
on the prior permitted quintet E1 sum; that earlier sensitivity cannot be applied
to a different operator/transition. The combined case is explicitly calculated,
not the sum of independent errors. It still differs from evaluated NIST
M1 0.00563 and E2 0.0000211 s^-1. These are evaluated/theoretical references, not
new laboratory gf measurements. No observed energies were substituted or rates fitted.

The same-space Hamiltonian check is also calculated directly at n=5, rather than
transferring a percentage from n=3:

| n=5 case | 1D2–3P2 energy (cm^-1) | M1 A (s^-1) | E2 A B (s^-1) |
|---|---:|---:|---:|
| Same-space DC | 16404.72 | 0.00915918 | 0.0000374035 |
| Static-limit transverse Breit | 16394.39 | 0.00672885 | 0.0000272208 |

The 6300 M1 discrepancy falls from about +62.7% to +19.5% against evaluated
NIST; its E2 discrepancy is still about +29.0%. At 6364 the M1/E2 discrepancies
remain about +22.1%/+51.6%. The same-space RCI DC control reproduces the strict
SCF table at printed precision apart from its last rate digit. Breit produces a
26.5% change in the 6300 M1 rate, far larger than the numerical stopping correction.
This diagnoses missing physical interaction as a substantial error source, but
neither removes all error nor licenses transferring a core correction from n=3.
Larger-space core/Breit combinations, lower-layer relaxation and independent
relativistic implementations remain owed.

## Primary-lifetime correlation test and rejected level inversion

The first n=3 triplet VV-SD attempt combined a correlated triplet with the old
uncorrelated quintet destination. Their native energies reversed the known
triplet/quintet ordering: the native `to_5S` table correctly describes the opposite
decay. The raw wrapper completion receipt is **invalid for a triplet lifetime**.
Offline replay refused it; the whole attempt is retained. A leading-state identity
and downward-order guard now rejects this before lifetime transitions. The bug
ledger records the defect and its regression test.

The corrected v2 uses a matching strict-qualified n=3 VV-SD quintet destination.
It retains the independent triplet seed's n<3 orbitals, optimizes n=3 only, and
uses no observed-energy or gf input. The triplet requires 156 iterations and ends
at maximum SCNSTY 5.87e-8, below native ACCY. All 13 one-photon channels, including
weak and computed-zero channels, are accounted for.

| Conditional model | Lifetime B / C (ns) |
|---|---:|
| Previous reference-configuration calculation | 5.51246 / 2.59095 |
| Balanced n=3 valence SD | 3.01070 / 2.02292 |
| n=4 valence extension | 2.50007 / 2.18079 |
| Smith1971 direct beam-foil experiment | 1.79, reported uncertainty 0.17 |

The residuals shrink but remain +1.22070/+0.23292 ns. Gauge disagreement remains,
and the correlated atom is not fully converged. This is a diagnosed physical-model
improvement, not a calibrated uncertainty or validation pass. Smith's reported
combination/confidence and experimental dependence limitations remain in the prior
primary-source receipt; Brooks1977 reuses Smith and is not another experiment.
The comparison is exposed/nonblind. Higher-order multiphoton contributions are not
certified. Native transition energies and both gauges are retained separately. The n=4
triplet and matching quintet reach the orbital criterion after 252/92 iterations.
Its lifetime residuals remain +0.71007/+0.39079 ns; the C-gauge comparison worsens
from n=3. An increasing active space or shrinking gauge disagreement is not proof
of absolute accuracy. Inherited frozen earlier-layer differences remain explicit.

## Artifacts, reproduction and outstanding work

Capsules are immutable, fully receipt-hash checked except omitted `mcp.*` angular
caches retained on Sirius:

- `strict_controls_n34_v1.tar.gz`: strict build/software controls, n=3/4 continuations,
  multipoles, rejected lifetime v1 and separate ground core/Hamiltonian cases.
  SHA256 `b9e614ecf72bee57d425ffef163926a1a9da652520ad32fa79f8e1702ddcf915`.
- `strict_lifetime_balance_v2.tar.gz`: matching quintet, corrected n=3 lifetime,
  combined ground core/Breit case and the then-current numerical plan.
  SHA256 `2652da1d6e69a8c8094ad09c70fb95f0b5d8ccf7330274740ef136dac14ed113`.
- `strict_extensions_v3.tar.gz`: completed/verified n=5 orbitals and multipoles,
  deliberately interrupted slower runs, n=4 lifetime extension, and rejected
  first finer-grid attempt. SHA256
  `afbb007343bb160e0b7a3b138ea4c284e43a4024993c22a405685d27b6e33f16`.
- `strict_grid_correction_v4.tar.gz`: explicit tighter ACCY, verified finer grid,
  completed orbital continuation and multipoles. SHA256
  `c7f444d8132299d6a85d2c47e680932f28317c7942a816e91824ebc479ea6b52`.
- `ground_n5_hamiltonian_v1_parts.{DC,Breit_state,Breit_multipoles}.tar.gz`:
  direct n=5 DC/static Breit evidence. The original 10.7 MiB capsule was refused
  by the commit size hook. Lossless standard-tar hardlinks and three parts preserve
  **all 97 native paths byte-for-byte**, including the original receipt and warnings.
  Every part is below 5 MiB; join before receipt verification. Full archive, part and
  member hashes are in `source_curation_v2/packaging_manifest_v3.json`; the original
  archive remains on Sirius and in local `/private/tmp`. No size-hook bypass used.

Sirius uses Python 3.14.4, GNU Fortran 15.2.0 and the existing pinned toolchains;
all owned jobs use OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and nice 19. Other
workers are untouched. Mac Python 3.12.13 performs native-text replay/tests only.
The existing 21-row refinement-debt report remains informational; none is discharged.
The library audit still stops on pre-existing civis2013/mallinson2024 vocabulary
errors before traversal. GitHub-hosted CI status is infrastructure metadata only.

Representative Sirius commands (task-root relative; choose fresh output directories):

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_strict_scf.py --build-from toolchains/grasp --out toolchains/grasp-strict-scf-v1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_strict_scf.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --seed runs/forbidden_ladder_v1/n3_ground --n 3 --cycles 500 --out runs/strict_n3_v1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_ground_core.py --grasp-bin toolchains/grasp/build/bin --seed runs/strict_n3_v1 --out runs/ground_core_v1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_triplet_correlation.py --grasp-bin toolchains/grasp/build/bin --strict-bin toolchains/grasp-strict-scf-v1 --ground runs/strict_n3_v1 --triplet-seed runs/coverage_v8/triplet_3s --quintet-seed runs/strict_quintet_n3_v1 --out runs/triplet_correlation_v2
python scripts/rya1354_error_resolution.py
```

Replay output `source_curation_v2/error_resolution_v3.json` includes complete channel
inventories, numerical comparisons, rejected-attempt disposition, source-summary
hashes and explicit null uncertainty/covariance. Remaining physical errors require
larger/core/radial and independent-method checks; no precision, grade or adoption
claim follows from native orbital convergence. Distinct source review, reviewed
criteria, calibrated uncertainty and independent human scientific review remain owed.


Focused validation command (Mac native-text replay and synthetic refusal checks):

```sh
/Users/ryanschmitt/codex/venv/bin/python -m pytest tests/test_rya1354_strict_scf.py tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

The replay checks immutable capsule hashes, exact declared solver patch and one-cycle
bitwise equivalence, all physical channel inventories, leading identities, native
orbital criteria and printed ACCY, rejected reverse-direction lifetime and custom-grid
metadata. Two characterized defects were fixed in this turn's research drivers;
original failed/misqualified receipts are never rewritten. Remaining model errors are
reported rather than turned into Gaussian error bars.


Final local validation: **68 focused tests passed**; register freshness normal and
PR-mode checks pass, zero state surfaces changed. Authored whitespace checks pass.
Production pipeline/products, canonical/reference line lists, catalogs and stellar
parameters remain unchanged against the original task base. All owned solver jobs
in this campaign have finished or have an explicit retained interruption disposition.
The pre-existing SEQUENCE ledger gap still requires the hook's disclosed draft-push
exception; no hosted CI success is inferred.

Exact next work: reconstruct the lower-layer sequence with strict orbital stopping
from its first layer, expand core plus transverse interaction together at larger n,
reoptimize rather than only add fixed-orbital core configurations, extend the
balanced triplet lifetime sequence with independent-method checks, and quantify
radial/tolerance/physical convergence before uncertainty calibration. Preserve
physical identity and ordering, source independence/covariance review and human
scientific gates. RYA-1354 stays In Progress and PR 616 stays draft.
