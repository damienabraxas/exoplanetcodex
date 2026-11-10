# RYA-1281: paired H2 pressure-collider experiment

**Scientific HOLD.** This is an explicitly partial FeI pressure intervention on
fixed oxygen populations, not a molecular EOS/continuum/O-collision calculation.
The H2 source/domain qualification remains in [the chemistry review](h2_table_control_review.md).
Its small unassigned Eq7 reconstruction difference is unchanged; published Kp
values are used without fitting constants or claiming exact equation parity.

## Inputs and intervention

Same reconstructed full56-depth solar5772/g4.438/xi1/O8.70 bound experiment,
experimental density-integral/log-native-top column and finite cap, full O atom,
1601-point FeI continuum support,41 exact source-matched FeI lines and unchanged
native pressure recipes. Cowley Stark is fixed in both cases. Saved nominal O
populations are fixed, rather than reused as an unevaluated coupled solution.
Air-to-vacuum conversion occurs once for internal RT; output/comparison is AIR.

Two new serial Mac replays use the same4001-point7763-7783A radiation/output
sampling, radiation/cap1e-8,1000iterations/180s bounds, five flux rays and one
context thread. All source/atmosphere/atom/coefficient bindings and original
native pressure/Stark receipts are revalidated before use. Implementation
commit544c835e precedes both replays. No original source or old output changed.

- **off:** native HI/HeI, H2 explicitly omitted, exactly reproduces the previous
  fixed-population Cowley normalized flux array.
- **molecular:** published-table H2 chemistry at fixed ne/totalH supplies depleted
  HI and H2; HeI unchanged. HII and Hminus ratios participate in hydrogen-nuclei
  closure. The previous atomic reservoir did not subtract Hminus from totalH;
  this intervention includes that bookkeeping along with H2 formation.

Only the pressure term uses these replacement colliders. Native H/He/O/Fe LTE
populations, ne, continuum and O collision rates are not replaced. Tests and
runtime identity assertions discriminate those boundaries. Added extinction and
LTE emission enter the same transfer equation; continuum arrays in the paired
spectra are exactly equal. No spectrum multiplication or fitted opacity/damping.

The source collision density is HI+0.41336HeI+0.85H2. Across the actual41x56
line/depth widths, relative to native H/He pressure:

| Contribution | Minimum | Maximum |
|---|---:|---:|
| H2 addition | +0.000015057 | +0.000643371 |
| HI depletion | -0.001513813 | -0.000035387 |
| Net pressure change | -0.000870442 | -0.000020330 |

Thus adding H2 without depleting HI would get the sign wrong for this intervention.
Actual backend colliders are captured and checked against every callback backend;
width CSV includes separate additive/depletion terms and actual total widths.

## Spectral results

These are fixed +/-0.5A diagnostic aperture areas, not isolated oxygen EWs or an
abundance. All numbers below are mA; small differences are numerical experiment
outputs, not claims of astrophysical precision.

| AIR centre A | Atomic control | Molecular colliders | Difference |
|---|---:|---:|---:|
| 7771.944 | 72.13049279 | 72.13046189 | -0.00003090 |
| 7774.166 | 62.41461284 | 62.41455491 | -0.00005793 |
| 7775.388 | 48.93480417 | 48.93471074 | -0.00009343 |

Maximum normalized-flux differences inside those apertures are
5.63e-8/9.89e-8/1.43e-7. Maximum across the whole wider region is7.786e-5;
it belongs to neighboring Fe opacity and is not an oxygen-window percentage.
This bounded correction is far too small to supply the several-mA gap. It does
not bound omitted molecular continuum, O-H collisions or full-EOS effects.

Canonical broadening plus gravitational shift gives72.10771/62.38630/48.93573
versus shipped Kitt77.55651/67.20039/53.79430. Original observation unchanged.
Zero-RV and canonical-shift variants both retained, no optimized velocity or LSF.

Existing ratified ORCH continuum controls rerun with this model: guided observed
areas74.58861/65.93433/46.77588 versus IAG-reference81.73581/70.81716/57.70338;
p95-envelope74.10889/63.73224/49.94709. Guided levels approximately
0.99679289/0.99864493/0.99263676 numerically pass; initial anchor geometry remains
one-sided for first two lines. Full reference spread retained, no favorable
selection or new telluric correction. Continuum/reference uncertainty dwarfs
this pressure correction. Numerical estimator pass is not physical qualification.

[Comparison plot](../../../data/output/rya1281/h2_colliders_v1/analysis/h2_collider_control.png)
and [aperture differences](../../../data/output/rya1281/h2_colliders_v1/analysis/h2_collider_difference.png).
[Synthesis continuum](../../../data/output/rya1281/h2_colliders_v1/continuum/synthesis.png)
and [IAG continuum](../../../data/output/rya1281/h2_colliders_v1/continuum/iag.png).

## Validation, provenance and resume

Mac held Bride Python (LW0.16.1) for RT/tests; project Python for iSpec/pandas
analysis. No environment installs or Sirius jobs. Off/molecular runtimes
12.0085/11.7339s; process-lifetime peak RSS744161280/742326272bytes includes
imports/setup, not full3D/grid memory. Host loads and context threads in reports.
All owned jobs finished.

36 tests pass:14 H2/backend/source/pressure/Stark/trace controls plus22
cap/observable/reference/local-continuum regressions. Includes exact off identity,
unchanged native populations/scattering/ne, deepcopy, one-ULP Hminus cache refusal,
source tamper/domain/units and independent hydrogen root tests. Earlier analysis
invocation in Bride Python refused missing pandas; rerun in existing project
Python, no dependency installed. First regression command named a nonexistent
residual-test file and collected no tests; corrected actual suite passes.
Initial difference plot used a global y-scale obscuring aperture changes;
final reviewv3 uses local limits. Previous reviewv1/v2 remain preserved.
No verified new accepted science/production bug; bug ledger: none discovered.

External registered Mac roots:
`rya1281-h2-collider-off-v1`, `rya1281-h2-collider-molecular-v1`,
`rya1281-h2-collider-review-v3`, `rya1281-h2-collider-synthesis-v1`,
`rya1281-h2-collider-iag-v1`. Sixteen compact artifacts and four figures are
byte-identical copies; actual generator invocations are in GENERATORS.yaml.
New scripts expose required inputs via `--help`; full input/source hashes,
current code commit and budgets are captured in each input manifest/report.

Next: fresh matched microturbulence0.9/1.0/1.1km/s, using the already-declared
control amplitudes from `rya1281_mac_checks.json`, now at the solar5772/g4.438
source physics. Rebuild exact local-input coefficient certificates rather than
relabel an xi1 binding; prove producer xi handling and preserve native-only
admission guards. Then fresh O solves/rate checks, followed by wider0.5/1.5km/s
controls if informative. Keep atmosphere thermal structure fixed for this causal
xi test; self-consistent atmosphere response is a separate gate. Do not retune
xi to the Kitt gap. Full EOS/continuum/collision consistency, independent matched1D
O/model-dimensionality benchmark, broader local opacity and continuum anchors
remain open. C/N before Fe abundance or another star; no adoption, merge,
publication, full3D or all-band claim.
