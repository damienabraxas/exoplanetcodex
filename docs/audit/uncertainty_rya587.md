# RYA-587 audit — 2026-09-15

Baseline main: `40f0e12f`. Reference branch: `dbe88e26`.
Machine evidence: `data/results/rya587/audit_v1/`.

## Existing machinery and actual gaps

| Source | Reused evidence / remaining gap |
| --- | --- |
| RYA-158 | Type-A/Type-B stack and parameter lookup exist. Historical zero Type-B values with absent derivatives do not prove measured zero. Preserve original diagnostic artifact, including its NaN history. |
| RYA-282 | Fe SE/parameter layer exists; its narrower total does not satisfy all new component families. |
| RYA-283 | Complete Sun/Procyon worked budgets depend on real measurements still absent. No invented table. |
| RYA-1088/1089 | Sourced star uncertainties and solar xi sigma=0.2912 km/s exist. Do not transfer one Fe-I derivative to other pools. |
| RYA-1120 | Existing Term states, stellar join and per-pool xi machinery reused. Legacy solar Fe-I Teff bound is not the new exact-product perturbation evidence. |
| RYA-1213 | Shell outputs, stamped pairs, Reference importer and stamper exist. Later completed VIS runs are not all ingested; branch remains unmerged. |
| RYA-1214 | Atomic/molecular conversion uses common publisher and declares curvature sigma. Full base stack plus likelihood correlations, blends and CNO coupling still required. C/N/O live lists are empty. |
| Al/Si | Merged intake/measurement work routes through common publisher. No Al.json or Si.json on baseline main; no canonical complete budgets supplied. |
| RYA-1220 | N method adjudication in progress. CN/NH/C/O covariance and indicator/model decisions require scientific work, not merely new fields. |
| Existing differential modules | Matched abundance differences exist. Legacy nearest-wavelength pairing is not a complete physical-identity/exact-nominal-pool covariance proof. |
| RYA-586 | `pipeline/reconciliation_band.py` is absent on audited main. Its proposed smoke command and historical before={Fe} cannot be reproduced. |

## Completed Reference xi shell outputs

58 completed directories were copied read-only from Sirius, including controls.
The existing importer validates run stamps and recovers **64 pools, 60 MEASURED,
four UNMEASURED**, no incomplete leg pairs. The prior committed artifact contains
**51 pools, 47 MEASURED**. All previous 47 measured product stamps still match.

Thirteen VIS products remain NOT_IN_CAMPAIGN despite completed measured outputs:

- HARPS: 1D-LTE, ENGINE-A, ENGINE-B-NLTE, synth-1D-LTE-gerber.
- IAG: the same four treatments.
- KPNO molecfit: the same four treatments.
- KPNO Kurucz2005 corrected: ENGINE-B-NLTE.

`reference_products_review.json` stages recovered stamps, preserves abundances and
legacy subtotals and marks the full reported uncertainty HOLD/null. It is not a
live feed. Four full-3D products retain xi applicability history; mean-3D is not
required. The four unmeasured legacy pool entries remain held for small paired N.

Shell legs use ±0.1 km/s and multiply the derivative by adopted sigma=0.2912 km/s.
These are completed measurements under the earlier method, not actual ±1σ runs.
Nonlinearity/asymmetry and exact nominal physical-pool verification remain owed.

## Completeness and unfinished work

Inventory covers all **26 canonical element symbols**, with Fe stages represented
as separate products. Historical raw scatter, N, SE and Type-B rows are preserved
in `legacy_diagnostic_rows`. One element-wide uncertainty cannot stand in for its
separate holdings/models/indicator pools.

Under the new full contract, complete set **{} → {}**. This does not erase the
completed Fe xi work; no product yet supplies the entire component/provenance/
covariance package. The older Fe-only definition is not directly comparable.
All 92 live Fe records fail the new evidence gate. Main's original feed and Fe
historical diagnostic row remain unchanged.

Before RYA-587 can be Done:

1. Integrate recovered Reference xi evidence with the unmerged RYA-1213 work.
2. Establish physical indicator IDs and product-specific ±1σ responses for every
   required parameter, including nonlinear-response assessment.
3. Supply line-specific transition covariance and observation/model evidence;
   resolve material missing terms or retain explicit HOLD.
4. Populate Al/Si/CNO/N budgets after their independent scientific gates permit it.
5. Supply joint evidence for target multi-element ratios. Covariance arithmetic is
   tested; no target ratio is publication-ready.
6. Migrate live products through the strict gate and restore eligibility tests
   before merging. No grandfather exception is implemented.

Final targeted tests: **86 passed** for new contract and legacy Fe uncertainty,
stellar join, sourced xi and differential modules. Existing eligibility/store run:
**36 passed, four failed** because historical records lack the newly required
proof. These failures expose migration debt and are not suppressed.

The baseline eligibility/store run passes all 40 tests, confirming the four failures are the stricter contract's migration debt. The register structure check independently fails on duplicate historical v135/v136 entries on both baseline and branch. No legacy Fe numerical artifact changed.
