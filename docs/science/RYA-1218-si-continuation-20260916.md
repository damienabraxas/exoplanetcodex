# RYA-1218 Si continuation — grading, abundance, and uncertainty

## Current evidence

The completed VIS synthesis artifacts were recovered from the Sirius checkout and
copied to `data/results/rya1218/recovered_20260916/`. HARPS, IAG, and corrected Kitt
Peak each have seven accepted Si I lines on the same 1D-LTE synthesis route:

| holding | A(Si I) | n | archived statistic |
|---|---:|---:|---:|
| HARPS corrected | 7.586 | 7 | 0.0526 dex |
| IAG | 7.625 | 7 | 0.0281 dex |
| Kitt Peak corrected | 7.657 | 7 | 0.0284 dex |

These are diagnostic values, not a live Si abundance. The archived estimator is a
median while the reported statistical field is scatter divided by √n; that pairing
needs an explicit estimator/independence treatment before RYA-587 admission.

The current canonical snapshot shows five of the seven used Si I lines at a
0.100-dex offset from the published AGSS21/Scott reference values. The two
Den Hartog 2023 lines (3905.523 and 4102.936 Å) are primary-laboratory matches,
but their recovered profile fits were quarantined by the physical-width gate.
No gf was changed to reduce the offset.

Si II 6371.371 Å was measured by profile EW on HARPS, IAG, and corrected Kitt Peak.
The controlled 1D-LTE synthesis returns A(Si II)=7.562 (HARPS), 7.618 (IAG), and
7.578 (corrected Kitt Peak), one line each. The quantiser floor is explicitly
unmeasured, so these are diagnostic only. Engine-A has no served solar Si correction
for this line.

Holding-to-holding matched Si I differences are small in the accepted seven-line
pool (HARPS−IAG mean −0.0006 dex; HARPS−Kitt Peak −0.0044 dex; IAG−Kitt Peak
−0.0039 dex). Engine-B and 1D-LTE labels are byte-identical on these recovered
products, so they are not independent engine corroboration. The two-line
Engine-A subset differs from 1D-LTE by −0.0158 dex on the matched lines, but its
pool is incomplete and cannot establish a full route correction.

## Production Si I abundance synthesis

The follow-up native-Linux production run used the same eight-line EW attempt ledger
per corrected holding and the controlled synthesis route. Si II 6371.371 A is reported
as NOT-IN-SYNTH-LINELIST in this Si I route, as expected. Seven Si I lines were
accepted in each holding:

| Holding | A(Si I), 1D-LTE | line scatter | Engine-A |
|---|---:|---:|---:|
| HARPS corrected | 7.586 | 0.0526 dex | 7.613 (n=2) |
| IAG | 7.625 | 0.0281 dex | 7.658 (n=2) |
| Kitt Peak corrected | 7.657 | 0.0284 dex | 7.663 (n=2) |

These are production synthesis measurements, but they remain diagnostic under
RYA-587: the source gf ladder is still ungraded for Si I, Engine-A serves only two
lines, and the complete covariance/uncertainty budget has not been admitted.

## Uncertainty disposition

`pipeline/si_evidence.py` maps each recovered product into the canonical RYA-587
component vocabulary. Every component is HOLD: the old 0.17-dex gf field is a
placeholder, stellar responses were not run on the exact pools, fit curvature is
not a validated likelihood, and corrected-telluric/continuum states do not by
themselves quantify residual uncertainty. The report therefore admits zero Si
products and preserves every value as diagnostic evidence.

Machine output: `data/results/rya1218/evidence_20260916/` contains the product
inventory, per-line physical IDs, gf-scale deltas, estimator diagnostics, and
matched-pool comparisons. The recovery manifest records source checkout and
artifact hashes. `scripts/audit_si_evidence_rya1218.py --out <new-dir>` refuses
existing output directories.

The follow-up uncertainty audit is in
`data/results/rya1218/uncertainty_audit_20260916_v2/`. It covers 31 product rows
and 157 accepted line measurements. Diagnostic scatter/SE values are available,
but all 16 RYA-587 components remain HOLD: the same physical line pool is not
yet certified independent, transition-data covariance is absent, and no exact
pool stellar-parameter, continuum, telluric, blend, or model response has been
validated. Therefore the audit admits zero products and reports no total sigma.

## Next work

1. Finish the CRIRES Y/H/J/K uncertainty admission and attach the all-band run hashes.
2. Build exact-pool profile-likelihood or independent-line uncertainty evidence,
   including Teff/logg/ξ/metallicity responses and covariance.
3. Reconcile the 0.100-dex Si I source-scale offset against the primary lab and
   AGSS21 tables without tuning.
4. Revisit 3905/4102 through the Deep synthesis path and preserve physical-width
   refusals where they remain.
5. Only then decide whether any Si I/Si II product can enter the live feed.

## Si uncertainty budget execution (2026-09-16)

The RYA-587 budget was rebuilt from the current VIS and CRIRES+ holdings in
`data/results/rya1218/uncertainty_budget_20260916/`. It inventories 31 stored product
rows and 157 accepted line measurements. The line-scatter/√N values are retained as
measurement diagnostics only: the production estimator is a median and independence,
fit likelihood, and covariance have not been demonstrated. Consequently all 16 RYA-587
components remain `HOLD`, `admitted_products=0`, and no `sigma_reported` is emitted.

A direct exact-pool Teff response attempt was also made through the canonical abundance
engine. It could not run on this checkout because the bundled MOOG executable is not
runnable on this host (`Exec format error`). No parameter response is fabricated from
that failure; the stellar response terms remain `HOLD` until the same pool is rerun on
the compatible synthesis host.

The initial generic EW per-line response export was withdrawn after identity review: its
8-line curated EW pool did not match the seven instrument-specific Si I synthesis lines
in the production products. It must not enter the uncertainty budget. The next run is an
instrument-level synthesis perturbation on those exact seven IDs on Sirius.

## Exact production-pool stellar responses (Sirius)

The instrument synthesis harness was rerun on Sirius for the exact seven-line HARPS
production pool. Teff perturbations (5672/5872 K) give a mean response of
`+1.2429e-4 dex/K`, contributing `0.00012 dex` at the adopted 1 K solar uncertainty.
Microturbulence perturbations (0.9/1.1 km/s) give `-0.05714 dex per km/s`, contributing
`0.01664 dex` at the adopted 0.2912 km/s allowance. Per-line derivatives and pool IDs are
in `data/results/rya1218/si_response_exact/si_exact_pool_response_summary.json`.
These two terms are now measured diagnostics; the product budget remains HOLD until
transition covariance and the remaining RYA-587 components are evidenced.
