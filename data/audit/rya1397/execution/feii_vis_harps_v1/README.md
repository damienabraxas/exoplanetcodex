# Solar Fe II HARPS VIS reference benchmark v1

**Executed on Sirius. Validation refit complete. Review-only products.** The ORCH dry run forecast one eligible synthesis cell. The derive stage was resolved through the ORCH descriptor and run with the exact stored 4200–6908 Å fit window and strict reuse checks. The generic matrix descriptor started from the wider holding bounds. The exact-stage runner guards the output stem and refuses any fallback fit of the untouched line.

- Seven lines refitted, one prior fit reused per fitting treatment, zero unexpected refits, zero new transitions. No listed line fell outside the executed pool.
- 1D-LTE and Engine B remained at 7.476 dex with eight accepted lines. Statistical and systematic stage terms remained 0.0842 and 0.0649 dex respectively.
- Engine A remained at 7.277 dex with two served lines and six unserved. Statistical and systematic stage terms remained 0.208 and 0.0649 dex respectively.
- All reported per-line abundances were identical at the stored precision. The source-stage outputs include the stage uncertainty terms. The assembled publication uncertainty field is blank for these stage outputs. No new publication budget was assembled.
- The committed Solar products and prior Fe stage files passed a post-run Git diff check on Sirius. The output was copied to a versioned review directory. Rollback consists of continuing to use the stored baseline.
- Storage after execution was 73 GiB free on the data mount and 239 GiB free on the home volume.

See the [refit ledger and raw products](products/), [product comparison](product_comparison.csv), [line comparison](line_comparison.csv), [ORCH dry run](orch-dry-run.json), and [generation log](exact-stage.log).

## Insertion validation

The RYA-1388 binder recognized all three V I candidates in canonical atomic data and bound zero against the stored 27-line HARPS set. A refit CSV does not admit a new line into a pool. The repository has a separate published-line admission path that preserves HFS components. The combined path still needs an end-to-end insertion test using a scientifically eligible transition, and the three V I candidates still need component reconciliation. No new transitions were fitted or adopted in this benchmark.
