# RYA-1354 — first 6,600-level signed LS replay, v37–v39

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. This is an angular-representation diagnostic of existing FAC wavefunctions, not an independent atomic-physics calculation. No new Hamiltonian/orbitals, rates, lifetime or grading result.

## Efficient replay and corrected control qualification

The prior audit evaluates Gram and vector products with cubic Python loops. The new audit proposes NumPy/BLAS products, then falls back to the exact original Python sum for ANY mismatch in a native eight-decimal coefficient or component weight. The original strict native comparison must pass; no numerical acceptance tolerance was added. The Gram residual uses BLAS only as a reported diagnostic, without a physical-accuracy threshold. Python3.14.4/NumPy2.3.5 and single-thread environment are recorded.

Four prior n3 cases reproduce ALL nondiagnostic qualified-report fields exactly, including 329,198 signed LS coefficients per case, full IDs/J/parity/energies/parent terms/dominant weights/component counts/targets. No natural rounding fallback was needed. A deliberately bad proposed product recovers the exact old result through 16 original-sum fallbacks.

The initial v37 report incorrectly labelled its vector-sign control: the selector checked only the first number on each line, skipped tiny coefficients, and crossed into an energy identity row. Its rejection was Missing LS vector identity, not coefficient-sign rejection. That raw report/status is preserved as historical execution; its vector-control qualification is superseded, not silently accepted.

The separately declared v39 control scans only coefficient-only rows of the first vector, all tokens, before the next header; flips a nonzero coefficient; records line/column/original/mutated data; and requires the exact expected rejection mechanism. All four old positive cases, forced recovery, and seniority/vector-sign/matrix-sign/stdout-ID controls pass. The actual vector-sign rejection is Expanded LS coefficient/sign replay differs. Future launch and metadata gates require corrected v39 schema, actual coefficient location and exact rejection reasons. An actual archived-J regression verifies that every energy identity line stays unchanged.

Initial v37 qualification took36.42s/291,432KiB; corrected v39 took27.11s/285,436KiB, no swaps. The numerical LS audit module is byte-identical across v37, v38 and v39. No native rerun or physical change was needed to correct the canary.

## First expanded native conversion

The first ground/coherent10 6,600-level case uses the original full FAC C/CM byte-for-byte and the pinned unmodified GRASP jj2lsj binary/source. Same full C/J/LBL/T prompts, all levels/compositions, unique labelling disabled. The old Z8 isodata file is converter support only; it does not replace or match the FAC physical nucleus. All qualified original CI receipts, orbital/potential/configuration checks and runtime/module provenance verify before conversion.

The native converter completes in181.60seconds,400,712KiB peakRSS. Full driver including original-CI verification and complete signed replay takes318.18seconds,2,576,108KiB peakRSS. Both exit0. All 4,354,454 LS coefficients are reconstructed, including signs; all 6,600 original IDs/J/parity/energies/full parent-term compositions pass. All 5,600 I3 display overflows are independently matched to binary-CM identities and explicit I4 stdout IDs; no energy-nearness or occurrence-only guess. The Gram diagnostic is1.7763568394002505e-15; coefficient/weight fallback counts are zero.

Targets are uniquely ID6 (1s2 2s2 2p3[4S] 3s, 3S, J1, odd) and ID5 (same parent, 5S, J2, odd). Their dominant printed LS weights are approximately0.97058 and0.96446: basis composition weights, not accuracy/confidence estimates or grades.

The first launch used --old-converted instead of the parser's --old_converted and exited2 BEFORE directory creation/native execution. The timing log is preserved as driver_argparse_failure.time; correcting only that CLI token launched the same declared inputs. This was a launch typo, not a numerical failure. Native v38 receipt/report remain immutable; final metadata qualification explicitly links the corrected v39 canary to the same unchanged audit module.

## Portable metadata and evidence

Independent Mac checks verify all archives/parts/member/fragment/native-receipt hashes, source/module/control-correction chain, original C/CM/support/binary identity, actual 6,600-level native label/CM/report metadata and full composition membership/counts. Full matrix/vector arithmetic stays on Sirius. The machine twin is data/audit/rya1354/source_curation_v2/fac_radial_angular_v39.json; reader scripts/rya1354_fac_radial_angular_metadata.py.

The portable reader initially demanded raw-float equality when recomputing a Sirius-derived printed-composition weight. Twelve values differ in the last binary bit under Mac versus Sirius libm, e.g ID54:0.9747127549361765 versus0.9747127549361764. All raw values/IDs are retained. The portable derived-field check now requires exact equality at the SAME native eight-decimal print precision; this does not change the strict signed matrix/vector/component-weight replay on Sirius. No physical uncertainty/tolerance was inferred or adopted.

Archives preserve: initial v37 control history114paths/1part/198,487compressedbytes; complete native v38 conversion132paths/123parts/123,164,746bytes; corrected v39 qualification117paths/1part/204,584bytes. Every part is <=5MiB. The single-pass lossless packer and existing independent decoder preserve every selected original file/byte; packing finishes before transfer. Full native C/J/LBL/T, original inputs, all scripts/reports/control and launcher history remain available.

The full focused science/governance suite passed: **234 tests in170.63seconds**, including the complete actual expanded archive metadata check, corrected-control coverage, actual coefficient-only mutation and Sirius-only math guard.

No original FAC/GRASP source/library/executable changed. Scientific canonical surfaces remain unchanged. Bibliography audit still stops before traversal on existing civis2013/mallinson2024 verified-state vocabulary errors; historical whole-ledger missing-field debt is unchanged. No clean-library or hosted-CI success claim. Draft push retains the documented exception for unrelated merged SEQUENCE debt, after manual blob-size checks. No merge/publication/Done.

**Resume:** three remaining 6,600-level LS conversions are NOT RUN. Use the corrected v39 canary and current strict gate, unchanged math module, same original C/CM/full-output controls, and exact CLI --old_converted. Generalize metadata checks per actual case; do not assume IDs6/5 without each case's full audit. Source/driver supports selecting ground100/equal10/equal100; use fresh declared outputs. Resource profile now supports sequential runs. Preserve full compositions and lossless evidence using the new stream packer. After all four LS cases pass, recompute both-gauge rates and lifetimes. Existing gauge/potential sensitivity, physical convergence, experimental/known-grade calibration and theoretical line grades remain owed.
