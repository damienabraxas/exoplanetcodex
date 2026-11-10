# RYA-1281 controlled native-Ng-off comparison

Scientific HOLD. The previous [operator correction and qualification](solar_solver_acceleration_review.md)
restored physical rate action and passed the native baseline, but temperature-only
still exhausted its original 600 s budget with repeated excursions. This experiment
changes only the native Ng flag. It does not assume that acceleration explains the
failure or select a physical spectrum from a preferred convergence method.

The predeclared plan is `config/experiments/rya1281_solar_matched_transition_no_ng.json`.
It uses the same isolated transition-local plugin/build receipt, source atom and
collision hypotheses, bare background, 56 depths, xi=1 km/s and fixed original
geometric cap. Every case starts independently from the original completed model-H
fractional populations, renormalized to its own O total; no failed-run seed or
sequential matrix initialization is used. All tolerance, radiation, iteration,
ray and 600 s whole-case limits remain unchanged. Only `controls.ng=false` and
the explicit named numerical control select the difference. The inherited scope
string's earlier generic reference to unchanged Ng is superseded by these explicit
fields and the appended Ng-off declaration; actual constructor options disable Ng.
The plan remains byte-identical to the already dry-checked, pinned input.

All four dry Contexts were qualified before this run. New RT work begins with a
fresh native 5750 K/logg4.5 baseline. Require independently reconstructed final
physical C+transition-R residual <1e-6, conservation <1e-10, converged radiation,
valid emergent/continuum extraction and original 1e-6 endpoint flux/continuum
parity screens. Only after baseline qualification run T-only5772/4.5, gravity-only
5750/4.438 and combined5772/4.438 serially; stop the matrix on a failed case.
The analyzer still requires all four complete cases plus gravity endpoint parity
before the factorial observational comparison. No failed attempt is promoted.

Emergent coverage remains three VIS O777 lines, 4001 AIR samples over7763-7783 A;
full source SE retains6796 wavelengths,41 levels,256 lines and22 continua. No
all-band, physical3D, abundance, source/reference, EOS/chemistry/provider or model
admission follows from numerical convergence. Current observer-frame conversions,
RV/continuum/telluric handling and pinned ORCH replay are unchanged.

Origin refresh remains7f023fd7; library audit remains212 documents/127 accounted/
85 unindexed/12 off-Mac, nonzero exit. All runs are Mac-only, serial, one thread;
no Sirius, CI, default/installed-backend change, production, merge or publication.


## Completed Ng-off control and independent replay

Native baseline qualifies in339outer iterations/396.91995s; fresh independent
physical residual2.48239e-7, conservation6.66e-16. Maximum normalized flux and
relative continuum differences from the original direct endpoint are1.42672e-7
and1.10800e-13; both pass the original1e-6screens. Against the corrected Ng-on
baseline, differences are2.55792e-9and8.88e-15. Initial/LTE populations, source
physical inputs, collision arrays, wavelengths and population-map helper are exact.
The independent comparison recomputes physical incoming/outgoing flow from saved
offdiagonal rates instead of trusting the reported residual.

Ng-on has38population-change increases, exactly on iterations9,14,...194;
Ng-off has none. Ng-off takes more outer steps(339vs197), with nearly the same
single-run wall time(396.92vs398.33s). This is a measured single-atmosphere control,
not a general speed or convergence theorem. The two qualified spectra are almost
identical. A standalone figure and history table preserve those observations.
Additional input/map guards in comparisonv2 leave the original history CSV and
figure byte-identical to v1; no scientific values change.

Temperature-only Ng-off exhausts600.06059s/239outer iterations. Its last population
change is0.000132623, old-radiation rate residual3.29814e-5, J change0.000132868.
These are incomplete in-iteration diagnostics, not a fresh final physical result.
Gravity and combined solar are not started after this failure. The trajectory is
smoother than Ng-on, but Ng removal is insufficient to solve the slow coupled mode.
All original Ng-on/off attempts remain, including initialization-only failure seeds.

## Prospective fresh-rate-screened extrapolation

The next bounded control uses the same corrected transition-local MALI map with
native Ng off, plus log-fraction Anderson(memory5/delay10/rcond1e-10). Each plain
and proposed population receives a NEW fixed-population radiation/cap solution.
Only a proposal with smaller fresh physical rate residual than the direct-map
step is accepted. Rejected proposals restore the plain state by actual radiation
re-evaluation; no radiation/rate snapshot from the other population is reused.
This deliberately costs extra RT calls inside the SAME600swhole-case budget.
All original joint/physical/radiation/iteration/ray thresholds are retained.

A separate explicit backtracking plan reduces an oversized proposal uniformly in
log coordinates by powers of two, normalizing after each trial, until its ACTUAL
normalized log-population step is<=0.05. There is no per-component clipping, floor,
changed physical rate or larger step bound. The older reject-only policy remains
unchanged by default and its dry contexts are retained. Every proposed population
still faces the fresh physical screen and original final admission checks.

Two analytic nonlinear-rate tests distinguish a stale exact SE solution from the
true feedback root, reject the former and verify that the evaluator is restored
to the selected state. The step test verifies positivity, O conservation, actual
log-step bound and improvement along a slow map. An adversarial comparison test
rejects changed background, map helper/method and compute budget before output.
All6targeted tests pass. Four backtracking dry contexts pass; actual baseline
qualification is a separate requirement and is running under the original limits.
The current control has not qualified any new solar-parameter result.


The bounded fresh-screened baseline PASSES51outer iterations/245.05756s, fresh
independent physical residual1.534120685e-7 and conservation7.77e-16. It accepts
4proposals and rejects/restores37; the final field/population changes pass the
original joint threshold. Its spectrum/continuum differ1.45670e-7/1.10134e-13
from the original direct endpoint, and2.99804e-9/9.10e-15from plain Ng-off;
all original1e-6parity screens pass. This qualifies the controller for this
baseline. Temperature, gravity and solar now run serially from the SAME original
common fractional seeds; no stellar-parameter conclusion is claimed yet.

Subsequent milestone: [all four matched columns now pass](solar_full_cap_milestone_review.md). Earlier failures remain preserved; scientific HOLD continues.
