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

## Next work

1. Finish the Si II controlled synthesis outputs and attach their hashes.
2. Build exact-pool profile-likelihood or independent-line uncertainty evidence,
   including Teff/logg/ξ/metallicity responses and covariance.
3. Reconcile the 0.100-dex Si I source-scale offset against the primary lab and
   AGSS21 tables without tuning.
4. Revisit 3905/4102 through the Deep synthesis path and preserve physical-width
   refusals where they remain.
5. Only then decide whether any Si I/Si II product can enter the live feed.
