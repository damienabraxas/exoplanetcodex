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
