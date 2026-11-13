# FAC numerical stopping and basis export controls

The minimal FAC oxygen model now has explicitly verified numerical stopping
controls and exported basis/mixing records. This does not admit physical LS
identities, rates, accuracy or a grade. No transition calculation ran.

## Runtime source versus documentation

At pin9ab6def, consts.h defines OPTTOL3.0, OPTNITER2048 and OPTSTABLE.5.
The same checkout's manual describes defaults1e-6/100. Runtime source governs
executable defaults; both claims are preserved. In OptimizeILoop, multiplier t
sets tol0=t*1e-4, tol1=t*1e-8 and atol0=t*1e-5, scaled further by native dynamic
stabilizers. These measure relative/absolute orbital-energy changes and potential
change, not the GRASP orbital-function ACCY.

The source evaluates the preceding cycle at loop entry and later prints a
positive native miter. Independent replay requires that marker, guard history,
muconv, and the previous and current printed numerical conditions. No additional
tolerance/epsilon waives failed conditions. Correction to v17's shorthand: the
native log has18 printed rows indexed0..17, not17 rows.

## Native negative and positive controls

Before execution, comment`e7c5cf3c-f8c8-4311-a69d-f494decf1b09` declared separate
cold-start cases with unchanged minimal O configuration/nucleus/physics and
level/basis-only scope. The capped negative control returned exit0 and wrote
levels while failing its stopping guard. Native output existence or exit status
therefore cannot establish numerical qualification.

| Case | Printed rows | Outcome |
|---|---:|---|
| multiplier3, cap1 | 2 | Exit0/output exists; rejected |
| explicit multiplier3, cap256 | 18 | Printed/native-marker conditions pass |
| multiplier.3, unchanged inner settings | 257 | Cap reached, miter-1; rejected |

The first campaign stopped and remained FAILED_RETAINED. The planned .03 case
did not run. It is archived in `fac_controls_failed_v18`:38 paths,905629bytes.
No failed wave/potential or output was used as a later initializer.

Source orbital.c EneTol exposes separate inner enerelerr/enerelerr1/eneabserr
settings, initially1e-4/1e-8/1e-5. The outer-only .3 case stalls with absolute
energy change1.70206e-6 above its scaled bound, despite smaller potential changes.
That motivates testing consistent inner/outer settings; no physical root cause or
calibrated accuracy is inferred.

Comment`a01d829e-d5c1-498e-9e59-d5ba32323f2d` separately declared a new cold-start
campaign. It tightens all three supported inner settings by10 and100 while
tightening outer multipliers by the same factors, retaining native dynamic
stabilization, cap256 and identical model inputs. No source patch or stopping
criterion relaxation occurred.

| Coherent case | Outer multiplier | Inner scale | Printed rows | Final absolute energy change |
|---|---:|---:|---:|---:|
| coherent10 | .3 | 10 | 19 | 1.06379e-7 |
| coherent100 | .03 | 100 | 20 | 2.87588e-8 |

Both pass the pinned printed/native-marker conditions. Maximum differences of
printed level energies by stable row index are6.75e-6eV (default to coherent10)
and5.0e-7eV (coherent10 to coherent100). These are numerical sensitivity
diagnostics, not physical LS mapping or uncertainty estimates. New archive
`fac_coherent_complete_v18` preserves28paths in590159bytes.

## Basis integrity and remaining mapping work

Mode0 exports the native jj basis/mixing table; mode10 exports C and G92MIX CM
records. Independent replay verifies8electrons,15CSFs,5orbitals,15levels and
7symmetry blocks, exact record dimensions, finite entries and complete unique
level indices. Maximum algebraic orthonormality residuals are8.88e-16,5.74e-16
and1.11e-15 for default/coherent10/coherent100. These diagnose exported matrix
integrity only; no physical state is assigned from them.

The manual warns of jj2lsj compatibility bugs and links doc/misc/jj2lsj.md, absent
at this pin. Do not assume that a converter run or normalized matrix establishes
correct angular phases, coupling conventions, global0-based FAC level indices,
LS parents or term assignments. This compatibility audit is the next step before
using FAC oxygen rates or matching states to GRASP.

`rya1354_fac_controls.py` owns the first declared campaign and stopping replay;
`rya1354_fac_coherent.py` owns the separate coherent campaign;
`rya1354_fac_control_audit.py` independently replays both, exact inputs, preserved
source/default discrepancy and binary exports. Machine record:
`source_curation_v2/fac_controls_v18.json`. Both source trees/executable pin remain
unchanged; source snapshots and executed drivers are retained. Two characterized
bug-ledger rows record the zero-exit/capped-output hazard and manual/runtime
default mismatch; local guards prevent result admission, upstream behavior stays.

Focused scientific/governance suite143passed100.69s, including complete native
failure/coherent replay, actual zero-exit capped rejection, wrong-threshold and
inner-input rejection, and corrupt exported-record refusal. New archive files
remain below5MiB; lossless/native hashes, whitespace, unchanged production/
reference/catalog/star comparison against task base and register checks pass.
Standing bibliography audit still fails on existing civis2013/mallinson2024
verification vocabulary before traversal; no library-success or remote-CI claim.

## Resume

All native campaigns completed; there is no pending solver from this checkpoint.
Sirius root remains /mnt/codex-data/codex/rya1354. Native outputs are
runs/fac_numerical_controls_v18 and runs/fac_coherent_controls_v18, with immutable
drivers under evidence/fac_control_driver_v18 and fac_coherent_driver_v18.
Do not overwrite/relaunch these directories. Reproduce the offline audit with:

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1354_fac_control_audit.py --out <new-audit-json>
```

Next inspect the preserved FAC export writer, compare native basis/mixing against
exported C/CM and the actual GRASP converter reader, and establish an independently
checked recoupling/physical-state mapping. If the missing compatibility note is
needed, record that evidence gap rather than applying an inferred patch. Only
after verified identities and numerical/operator completeness should a separately
declared rate calculation proceed. No observed-energy fitting or closer-gauge
selection.

EXPOSED_NONBLIND; scientific HOLD; adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. Prior GRASP/core/joint failures and n3
component limitation remain. No canonical gf/grade/product/reference/star, source
spectrum, frozen gold, publication, merge or Done change.
