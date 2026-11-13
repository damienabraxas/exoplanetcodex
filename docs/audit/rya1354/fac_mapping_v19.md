# FAC export transport and minimal-model LS mapping

Both tightly controlled FAC cases now have verified native/export/LS
representation consistency for the minimal oxygen configurations. No FAC
transition, gf or lifetime was calculated; no accuracy or grade is admitted.

Before conversion, `rya1354_fac_mapping.py` reconstructs coupled-shell identities
from FAC's exported fixed-column CSFs and matches them to the native mode0 basis.
It handles only the pinned writer's closed-shell omission, last-full-shell
fallback and zero-J intermediate-coupling carry. No physical occupation changes.
The actual basis permutation is resolved, all43 signed coefficients match binary
doubles exactly at native %.8E print precision, and15 global0-based FAC IDs remain.
Both coherent10 and coherent100 pass; a substituted coefficient sign is rejected.

## Unmodified angular conversion

Declaration`cf87be4c-0c79-44ca-8963-521000e7dc5c` preceded running the pinned
unmodified GRASP jj2lsj converter on copies of both exact FAC C/CM inputs. CI
input, all levels, zero omitted composition/print cutoff, full LS C/J/T output,
unique labeling disabled. Original FAC C/CM remain byte-identical. Converter
source snapshots, executable hash and every prompt/output are preserved in
`fac_recoupling_complete_v19`:94 paths, one140877-byte lossless part.

The converter's Z8 isodata is angular-program support metadata only. It does not
replace or match FAC nuclear physics, recompute a Hamiltonian, optimize a radial
function, or substitute an energy. FAC's central-potential/nuclear/model limits
remain. The Hamiltonian implementation is independent; this angular converter is
shared with GRASP and is not claimed independently implemented.

The actual reader stores input IVEC values as label IDs, not array offsets, and
includes the parity-offset correction already in the pinned source. Missing
manual-linked jj2lsj compatibility note remains recorded. No inferred patch was
applied. In particular, the first number in `.lsj.j` is dominant-CSF index
`max_comp`, not the original FAC ID. The audit uses documented block/vector order
and explicit `.lsj.lbl` IVEC identities instead of relabeling levels from that
display field or from energy proximity.

## Independent replay

`rya1354_fac_recoupling_audit.py` checks exact original input/controls/provenance,
native-to-export transport, seven transformation block headers, full43-entry
angular matrices and all15 original IDs/parity/J/native binary energies. It
reconstructs every full LS eigenvector by multiplying the original FAC vector by
the saved angular matrix, preserving signs/order, and requires exact native
eight-decimal coefficient formatting equality. Printed LS label coefficients and
weights are cross-checked too. Transformation orthonormality residual is8.88e-16
in each case; no value here is an accuracy uncertainty.

All15 dominant minimal-configuration LS terms appear once with the allowed
parity/J distribution. Three analytically unambiguous term blocks pass:
even J1 is3P, odd J0 is3P and odd J3 is3D, each with printed unit composition.
This validates limited term controls and transport, not every angular table in
the package or physical completeness beyond these configurations.

| Target | Original FAC ID, both cases | Parent and term | Weight from printed LS amplitude |
|---|---:|---|---:|
| Triplet upper | 6 | 2p3(4S)3s 3S1 | .99960781846 |
| Lower quintet | 5 | 2p3(4S)3s 5S2 | .99974125674 |

All smaller admixtures remain recorded. No purity cutoff, experimental energy,
observed gf/lifetime or closest-model choice selected these states. Mapping is
specific to the minimal configuration and these recorded numerical cases.
Machine report:`source_curation_v2/fac_recoupling_v19.json`.

Initial regression exposed JSON integer-key/string-key roundtrip mismatch in the
new report. Keys were canonicalized to strings before serialization; original
integer FAC IDs remain explicit. No coefficient/ID or scientific criterion changed.
All7 focused mapping checks then passed, including actual full native replay and
LS-sign/global-ID/input substitution refusals. No new upstream defect was
discovered; the missing note/display-index semantics are guarded evidence limits.

The focused scientific/governance suite150passed97.98s. After adding explicit
numerical-source requalification to the final angular auditor, all7 mapping tests
passed again. Whitespace, native archive/hash replay, <=5MiB size checks, unchanged
production/reference/catalog/star comparison and register checks pass. Standing
bibliography audit remains blocked by existing civis2013/mallinson2024 verification
vocabulary before traversal; no library-success or remote-CI-pass claim.

## Resume

All native jobs completed. Sirius root remains /mnt/codex-data/codex/rya1354;
completed run:`runs/fac_recoupling_v19`; executed drivers:
`evidence/fac_recouple_driver_v19`. Do not relaunch or overwrite. Reproduce with:

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1354_fac_recoupling_audit.py --out <new-audit-json>
```

Next declare a minimal FAC rate calculation with both gauges, explicit multipole
requests and complete lower-state/channel accounting. Recheck exact numerical
model, basis/mixing/level identities against this mapping in any new native
process; do not assume regenerated level order. Keep FAC model/correlation/
nuclear differences from GRASP explicit. Only a complete, physically identified,
numerically checked decay inventory may produce an exploratory model lifetime;
it still does not establish accuracy, calibrated uncertainty or a grade. Known-
grade external and genuinely blind validation remain owed.

Pinned FAC consts.h sets DMODE=M_NR(1), so future relativistic-operator calculations
must explicitly request mode0 using SetTransitionOptions(gauge,0); gauge1 is C
and gauge2 is B. Do not rely on the non-relativistic default or the default
all-multipole sum to claim relativistic/operator completeness. No rate calculation
has used these settings yet.

EXPOSED_NONBLIND; scientific HOLD; adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. Prior GRASP/core/joint/FAC failures and n3
component limitation retained. No canonical gf/grade/product/reference/star,
source spectrum, frozen gold, publication, merge or Done change.
