# RYA-1281: fresh solar microturbulence slice

Scientific HOLD. Ryan requested a broader sweep on2026-11-03; this is the first
local slice, not a completed parameter program. The [broader inventory](parameter_sweep_plan.md)
covers wider xi and temperature, gravity/composition, physics/reference choices,
resolution/boundaries and eventual interactions. Current1D O777 VIS experiments
are not full3D, full-spectrum, molecular CNO or adopted abundance products.

## Source preparation and fixed fields

Plan: `config/experiments/rya1281_mac_xi_slice.json`. Three explicit named cases
at xi0.9/1.0/1.1km/s, solar5772/g4.438/O8.70, same56-depth thermodynamic structure,
composition, density-integral/log-native-top column/cap, O atom and source MN/
TAUTSP collision conventions. Same41 exact source FeI pressure lines and Cowley
Stark in each solve. H2 remains omitted to retain the current atomic EOS,
continuum and O collision baseline; prior collider-only H2 test is separate.

Pinned original TS input.f reads XIFIX and its following XIC value; BABSMA's
ASCII-MARCS branch assigns XI(K)=XIC at every depth. Existing qualified isolated
point build/reference receipts match earlier ts_points_v1. Original source,
data, binaries and earlier coefficients remain unchanged.

Fresh native identity plus three solar producer cases each run baseline and
reporting builds. Producer output xi is checked against native default REAL/
E14.7 serialization:900m/s request writes899.99998,1000 and1100 write exactly
those values. T/Pe requests and output quantization, all56 rows and1601 wavelengths
are checked. Baseline/reporting serialized coefficients agree; native identity
matches its original reference exactly. FeI per-neutral bound-free coefficients
are array-identical across these xi values, as a result of actual evaluation,
not a relabelled xi1 packet. New certificates bind complete local fields and
requested xi separately;4 exact identities/62 input/schema/neighbor refusals.
Provider/interpolated-cap admission flags remain false.

New sweep helper constructs the xi-controlled base **before** binding and adding
the experimental column. Original native-only guards remain unchanged. This is
necessary because the earlier bound runner supplied its unmodified source base
to the certificate loader; that runner cannot consume a changed-xi certificate.
It would refuse, not silently accept an incorrect spectrum. New tests check
fresh900m/s context construction and a one-ULP live velocity refusal.

Source-preparation commitc9584280 precedes producer/binding work. Runner2abe4095
precedes all three O solves; helper bytes stay unchanged. All inputs and local
bindings except declared xi are checked equal to the original solar manifest.
No active abundance, EOS density or collision coefficients fitted to Kitt.

## Fresh coupled solutions and results

All three contexts start from the same checked earlier nominal O populations,
then perform new coupled O SE/radiation solves. Warm starts are not accepted
as solutions. Original1500iteration/600s, joint1e-6, conservation1e-10 and inner
radiation1e-8 controls are unchanged. Saved separate radiative and collisional
rates independently reconstruct the final physical rate matrix/residual after
fresh radiation, distinct from tiny in-loop algebraic residuals.

| xi km/s | Iterations | Total replay seconds | Fresh independent rate residual |
|---|---:|---:|---:|
| 0.9 | 71 | 91.083 | 3.026419e-7 |
| 1.0 | 16 | 28.811 | 4.137014e-7 |
| 1.1 | 76 | 93.363 | 3.021024e-7 |

Nominal1.0 normalized flux is array-identical to the previous coupled Cowley
result. The total new run time is213.26s (3.55min), serial Mac science. Process
peaks711163904/684834816/752074752bytes include imports/setup;1context thread,
loads in each report. These are not full-grid/3D scaling measurements.

Fixed +/-0.5A intrinsic diagnostic aperture areas in mA:

| AIR centre A | xi0.9 | xi1.0 | xi1.1 | Low-minus-nominal | High-minus-nominal |
|---|---:|---:|---:|---:|---:|
| 7771.944 | 71.69879 | 72.12826 | 72.59360 | -0.42947 | +0.46535 |
| 7774.166 | 62.05519 | 62.41175 | 62.79793 | -0.35656 | +0.38617 |
| 7775.388 | 48.67105 | 48.93311 | 49.21656 | -0.26206 | +0.28345 |

These responses are much larger than the prior H2 collider-only effect. The
positive/negative0.1km/s responses are slightly asymmetric: the local response
is not exactly linear. Three points cannot establish the wider curve or a
stellar-parameter uncertainty. This intervention varies transfer/SE xi at fixed
thermal structure, not a self-consistent microturbulent atmosphere response.

Canonical broadening+gravitational shift gives:

| xi km/s | 7771 area | 7774 area | 7775 area |
|---|---:|---:|---:|
| 0.9 | 71.67655 | 62.02699 | 48.67207 |
| 1.0 | 72.10550 | 62.38350 | 48.93412 |
| 1.1 | 72.57032 | 62.76962 | 49.21756 |
| Shipped Kitt | 77.55651 | 67.20039 | 53.79430 |

Kitt is AIR; internal RT is vacuum. Original observation and both canonical/zero
RV variants retained; no optimized shift, LSF, continuum or chosen xi. The
comparison does not measure an O abundance or isolated unblended EWs.

ORCH synthesis- and IAG-reference continuum legs rerun for **every case**.
Guided observed areas (xi0.9/1/1.1) are approximately
74.58881/74.58862/74.58840 for7771,
65.93373/65.93437/65.93511 for7774,
46.73590/46.77606/46.77639 for7775.
IAG-reference81.73581/70.81716/57.70338 and p95-envelope
74.10889/63.73224/49.94709 are unchanged. This much larger reference spread and
one-sided anchor qualification remain open; no favorable continuum is selected
and no second telluric correction applied. Old10K responses used an earlier
opacity subset; do not combine them into a new matched uncertainty budget.

[Full comparison](../../../data/output/rya1281/xi_slice_v1/analysis/xi_control.png)
and [local response](../../../data/output/rya1281/xi_slice_v1/analysis/xi_response.png).

## Validation, artifacts and next work

38 tests pass:4 new sweep/input/context/refusal tests plus34 bound-provider,
Stark/pressure/trace/cap/observable/reference/local-continuum regressions.
Plan extension to xi0.5 is tested without a runner rewrite; altered Teff,
abundance or depth in this xi-only phase refuses. Independent saved-R+C checks
pass for every case. No new verified accepted science/production bug; bug ledger
none discovered. Wider-axis execution requires its own validated input contract.

42 compact artifacts and8 figures are byte-identical copies at
`data/results/rya1281/xi_slice_v1` and `data/output/rya1281/xi_slice_v1`.
Full producer coefficient arrays stay in registered external holdings, not in
metadata-only certificate copies. Actual generators/helper writers and invocations
are registered. Final external roots:
`rya1281-xi-source-v1`, `rya1281-xi-coefficients-v1`,
`rya1281-xi-coupled-v1`, `rya1281-xi-review-v2`, `rya1281-xi-continuum-v1`.
Reviewv1 retained; v2 adds the response figure. Bride Python for RT/tests,
existing project Python for iSpec/pandas analysis; no installs/Sirius/CI load.

Next: widen xi to the previously declared0.5/1.5 controls with useful intermediate
sampling and representative resolution checks; source-matched wider Teff/gravity/
composition slices, followed by predeclared interaction points. Preserve failed
or uncovered cases and measure actual cost before scheduling a larger lattice.
Independent matched1D O/reference/atom lineage, full EOS/chemistry/opacity and
continuum/LSF qualification remain indispensable. C/N before Fe abundance and
another star. Scientific HOLD; no adoption/production/merge/publication change.
