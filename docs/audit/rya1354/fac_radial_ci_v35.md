# RYA-1354 — first 4s/4p CI and exact identity replay, v33–v35

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. No LS mapping, rates or lifetime were calculated in this checkpoint.

## Native calculation and retained failures

The first ground/coherent10 case extends the fixed-potential four-valence model from 2p/3s/3p/3d to 2p/3s/3p/3d/4s/4p, keeping 1s2/2s2 fixed. Independent magnetic determinant counting predicts 114 NR groups, 6,600 levels, 16 symmetry blocks, maximum block 916, and 4,354,454 coefficient slots. This is a restricted radial extension, not full n4 or atomic convergence.

The native v33 process exited zero in 26.61 seconds, peak RSS 1,017,180 KiB (about 0.97 GiB), zero swaps. It restored the paired v28 potential, repeated the v32 radial preflight, used the same disabled Breit/SE/VP/MS/configuration corrections and inner controls, and performed one Structure. No fresh SCF was run. Full mode0 native text and mode10 C/CM exports remain intact. Original FAC/GRASP source and libraries/executables were not patched; the previously declared read-only configuration diagnostic frontend was reused.

The original v33 integration receipt remains FAILED_RETAINED: suppressed display labels were not unique. For symmetry0, basis2 (group2/cfg0) and basis135 (group50/cfg0) both display sname4p2 and coupled4p-2(0)0, while their native principal-complex nc fields are respectively 2*2.4*4 and 2*2.3*2.4*2. Source ConstructLevelName suppresses filled shells in name/sname, but nc records occupation summed by principal n, excluding closed complexes.

The separately declared metadata-only v34 replay added nc. It remains FAILED_RETAINED because the all-closed 1s2/2s2/3s2/4s2 exported CSF needs ConstructLevelName's final-shell fallback: when sname is empty, retain the outermost filled NR shell (4s2). That exact source rule was added in a separately declared v35 metadata-only replay. Neither replay reran the atomic calculation or changed its inputs or physical tolerances.

## Qualified v35 result

Exact identity uses symmetry + principal-complex nc + NR sname + coupled jj name, reconstructed from the complete exported occupations and requiring a bijection. No energy-nearness or occurrence-order matching. All 6,600 exported CSFs and original level IDs, electron/J/parity inventory, 114 configuration groups, and before/after actual zero corrections pass. All 4,118,356 printed signed coefficients match the binary export at native print precision. Every unprinted binary slot must be exactly zero; all binary values must be finite.

The old nine and new three orbital header metadata and physical r/rVT0/P/Q samples are unchanged after CI. Protected potential/header/grid bytes remain exact under the source-mapped W/dW/dW2 work-array rule. All pre-CI bound-orbital and printed-overlap guards pass. This is not radial grid convergence or an uncertainty estimate.

The streaming auditor first reproduced all four prior v28 receipts, using 30,036 KiB peak RSS and 4.49 seconds on Sirius. The final v35 replay, including original-byte hashing and all four canaries, used 396,268 KiB peak RSS and 15.54 seconds with zero swaps. Native signed text and complete binary data are retained; bounded audit memory is not evidence reduction.

Inspected primary GRASP jj2lsj source has a 999-entry weight-sorting loop and I4 explicit IDs. This proposed model's maximum block916 and maximum ID6599 fit those bounds, but no new native LS conversion or LS-label admission has occurred.

## Evidence, validation and resume

- Native failed v33 archive: 170 paths, 37 lossless parts, 116,502,838 compressed bytes.
- Failed v34 and successful v35 replay archive: 215 paths, one 165,889-byte part.
- Every part is at most 5 MiB; strict membership, part/member/fragment hashes, reconstructed receipts and paired baseline/preflight links are checked.
- Independent local metadata replay: scripts/rya1354_fac_radial_ci_audit.py; machine twin data/audit/rya1354/source_curation_v2/fac_radial_ci_v35.json. It verifies both failed statuses remain retained and the successful v35 report exactly reproduces.
- Regression controls reject signed coefficient, unknown ID, missing nonzero, duplicate coefficient, record-marker and principal-complex corruption. The complete actual extended capsule covers both the collision and fully closed-shell fallback.

Initial sandbox SCP/SSH denials were retried with escalation. A new auditor was initially copied into the old v32 snapshot directory, then moved into fresh v33 before use; no existing file was overwritten and prior archived snapshots remain authoritative. No failed native output or criterion was silently replaced.

The bibliography audit still stops before library traversal on existing civis2013/mallinson2024 verified-state vocabulary errors. Existing whole-ledger historical schema debt is unchanged. No clean-library or hosted-CI success claim. Draft push retains the documented exception for unrelated merged SEQUENCE debt, with manual file-size checks; no merge/publication/Done.

**Exact resume:** run the remaining ground/coherent100 and equal/coherent10/100 fixed-potential CI cases using the qualified v35 identity replay, preserving complete evidence and comparing orbital/potential/configuration invariants. Then extend the full signed LS audit from four to six valence shells, qualify native LS identities/compositions/formatter limits, and only then recompute both-gauge rates and lifetimes. Existing gauge/potential sensitivity remains unresolved. Experimental/known-grade validation and calibrated theoretical line grades are still owed.
