# Matched CI result and independent-code feasibility

The qualified-orbital n5 DC/Breit calculation completed. Independent replay
verifies source provenance, identical input CSFs/orbitals/nucleus, physical
states/order, native Hamiltonian/root/mixing controls, transition inputs/operators
and all thirteen channels in each case.

| Model | Gauge B lifetime (ns) | Gauge C lifetime (ns) |
|---|---|---|
| Prior n5 SCF | 2.2844837477524753 | 1.7132068154876408 |
| Matched n5 RCI DC | 2.2844837477524806 | 1.7132068154584075 |
| Matched n5 RCI static Breit | 2.284670790492177 | 1.7125401648458216 |

DC reproduces the earlier rates to reported precision. Static Breit changes the
lifetimes by about+0.00819% and-0.03891%; it does not explain the substantial gauge
spread. No gauge is selected against experiment and no model spread is converted
to calibrated uncertainty. This is the same implementation, not the independent
relativistic comparison still owed. The six new native RCI stdout/stderr files
contain no lines matching the warning filter; raw outputs remain authoritative.

Native receipt SHA256:
`bf625c0260c2f151cc5f868eccdf815ba82450be5556c819d1742f976ea24cfe`.
Manifest SHA256:
`2d9e0bf2ec2138195d8b9eb4668418b906c3afc9d2594bd664c475406945e60a`.
`qualified_ci_complete_v17` preserves280paths across19 lossless parts totaling
33135213compressed bytes. `source_curation_v2/qualified_ci_v17.json` reproduces
the audit, both gauges and native E1 energy/strength decompositions. Native RCI
internal loading/interpolation is retained. Hashed restart caches remain Sirius.

## Independent implementation setup

After finding no existing FAC holding/task installation, primary upstream
[README](https://github.com/flexible-atomic-code/fac) and
[manual](https://github.com/flexible-atomic-code/fac/blob/9ab6defb6af62e1bc8251a35ff995809abef717c/doc/manual.tex)
were inspected. FAC uses relativistic CI with basis functions from a local
central potential. That differs from the GRASP orbital model; independent
implementation is not the same as equal approximations or demonstrated accuracy.

Source pin:`9ab6defb6af62e1bc8251a35ff995809abef717c`, detached task-local clone
`toolchains/fac-9ab6def` on Sirius. No source patch, global install or GRASP change.
Single-thread/nice19 build uses gcc/gfortran15, GNU17 C, Fortran argument-mismatch
compatibility, and full retained compiler/configure/make logs including warnings.
Tracked source diff is empty. Executable SHA256:
`300b9559ef50e92ae3d1209df5576975e5b9afbf5dcbc0e96e9768f478249d51`.
The executable remains native; archived receipts retain its hash/path/bytes.

Declaration comment`2722f8ef-8e97-48e4-a2f3-bbc878635262` preceded compilation.
The supplied Fe example completed and its37 parsed level identity/energy rows
match the pinned reference exactly (maximum printed-energy difference0eV).
This is not the full PFAC-dependent upstream checker or a physical accuracy test.
`fac_intake_v17` preserves14 primary-source/configuration files in84950bytes;
`fac_complete_v17` preserves25 build/example paths in103284bytes. Plans and
executed drivers are retained. `gu2008_fac` bibliography metadata was fetched
directly from Crossref DOI10.1139/p07-197; full paper was not acquired/read and
no local reference PDF is claimed. Raw Crossref response is preserved.

## Neutral oxygen: level identities only

Declaration comment`e6a8e6eb-25db-4c01-9801-577dd3c306da` preceded a bounded
neutral-O test in `runs/fac_oxygen_levels_v17`. Input configurations are
1s2 2s2 2p4 and1s2 2s2 2p3 3s1 with closed1s. Ground-only native-default
local-potential optimization prints its diagnostics; theory-only ConfigEnergy
workflow follows the upstream example. Breit, self-energy, vacuum polarization
and normal/specific mass shifts are explicitly disabled. FAC standard nuclear
parameters remain, without a point-nucleus matching claim.

The calculation produced15 levels,8 electrons and the exact declared even5/odd10
parity/J multiplicities and configuration identities. Closed2s label omission
is accounted for from pinned ConstructLevelName/IsClosedShellNR source semantics;
it does not change the actual input occupation. Native optimization diagnostics
contain18 printed rows indexed0..17; their stopping semantics are not equated with GRASP ACCY
or admitted as a converged oxygen model here.

No FAC O transition, gf or lifetime was calculated. No LS/parent/term assignment
to GRASP or experimental energy was inferred from J/parity or proximity alone.
`fac_oxygen_complete_v17` preserves9 paths in6787bytes. Independent
`rya1354_fac_summary.py` replays archive hashes, build/pin/executable provenance,
literal upstream example rows, exact oxygen input and level inventory, and
explicitly sets numerical convergence and physical LS-term admission false.
Machine record:`source_curation_v2/fac_feasibility_v17.json`.

## Checks and next work

Focused scientific/governance suite131passed99.91s, including successful native
matched CI replay and substituted Hamiltonian/source-wave/transition-mixing
refusals. The final5-test FAC suite passed0.05s, covering actual upstream and O
inventories, nonfinite/duplicate-level rejection and substituted physics/electron
count rejection. Archive files remain below5MiB; native hashes/lossless replay,
whitespace, unchanged production/reference/catalog/star comparison and register
checks pass. Bibliography audit still exits1 on pre-existing civis2013/mallinson2024
verification vocabulary before library traversal; no library-success/remote-CI claim.

All native campaigns from this checkpoint have finished. Do not relaunch or
overwrite them. Reproduce completed audits with:

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1354_qualified_ci_audit.py --bundle data/audit/rya1354/qualified_ci_complete_v17 --out <new-ci-report>
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1354_fac_summary.py --out <new-fac-report>
```

Next: audit pinned FAC OptimizeRadial stopping/diagnostics and explicit numerical
controls, determine supported basis/mixing export and unambiguous physical-state
mapping, then declare a controlled oxygen rate calculation. Only after identity,
operator completeness and convergence evidence should any independent-method
comparison be interpreted. Preserve differences in central potential, nuclear
model, correlation spaces and state optimization; do not fit observations or
admit a closer gauge/model as accurate.

EXPOSED_NONBLIND; scientific HOLD; adoption false; accepted lifetime, uncertainty,
covariance and calibrated sigma null. Earlier joint/full-core failures and n3
component limitation remain. No canonical gf/grade/product/reference/star, source
spectrum, frozen gold, publication, merge or Done change.
