# RYA-1281: local LTE opacity strengths at frozen anchors

Predeclared2026-11-05 before evaluation. Start3ded8e9d, stacked base4635c742,
draftPR598. Main advanced; exact refreshed SHA is
9a74d843112b526e671f0dce7132d9d2bc3b16a1. No merge into this historical
control branch: maintain the frozen population/atmosphere/source inputs.

Main's model-guided bound is now15percent versus historical3percent; no continuum
fit is executed here, and frozen historical anchors do not claim the current
ORCH rule. The388row region has one changed main row, TiI gf_127774, superseded
NIST-3.053 by Lawler2013-3.180 (RYA-1315). Withhold that row entirely; do not
choose/adopt a competing gf or label the historical catalog current truth.

Manifest config/experiments/rya1281_lte_screen.json owns the180s total Mac budget
and logtau500 nearest0 diagnostic depth. Both frozen columns, all56depths, all
actual historical coefficient-support anchors with source identities checked.
Construct only static continuum/background coefficients with saved O populations,
no formal RT, SE or flux fit. No population/EOS/ne/charge feedback or new opacity
installed. Keep existing modeled41FeI separate from other candidate species.

Require unique original source identity, single component, canonical/source gf
agreement, positive source natural width with logarithmic consistency, isotope0,
full-column partition domain/stage coverage and bound upper energy. Preserve all
withheld reasons, no zeros substituted for missing line physics. Natural-Voigt
widths only: pressure/Stark/isotope/molecular strength and line-wing qualification
remain owed. This screen is not a bound on omitted pressure wings.

Use the same native complete abundance vector and fixed T/ne/nH/xi; independently
verify Kurucz LTE stages via log-space Saha recursion and nuclei closure. Sample
source-local opacity at observer anchors inverse-shifted using the inspected
iSpec relativistic gravitational convention in AIR, then convert to vacuum once.
Macro/rotation/instrument convolution is not local opacity and is not applied.

Record per-line maximum line/continuum-absorption ratios and the depth/wavelength
of each maximum, plus the declared reference-depth ratios. Retain every-depth
species summaries over each footprint, weighted by absolute normalized historical
centre-prediction weights. Those weights describe footprint geometry, not a
formation depth/contribution function. Coefficients cannot predict flux/EW or
resolve the observed reference gap without a source-qualified transfer test.

## Screen results

Plan/module/script committed before evaluation at73d0b916. Both frozen columns
complete:242rows evaluated (41existing FeI plus201other candidates),146withheld.
Per-column disposition:29source identities not unique;42natural widths unqualified;
61upper energies at/above ionization;4canonical/source gf differences;6combined
FeI exclusions;3active-OI candidates;1changed-main TiI row. No arbitrary source,
width, gf, ion stage or population extrapolation substituted. Stage/partition and
isotope guards also apply; they cause no additional refusal in these inputs.

Each case samples318distinct observer-anchor positions and all56native depths.
Independent log-space Saha populations agree with native multiplicative recursion
within3.326e-15 of the total elemental population; nuclei closure4.441e-16. This
is numerical verification of the stated atomic LTE calculation, not atmospheric
chemistry/EOS or NLTE accuracy. Saved oxygen populations remain array-identical.
No formal RT/SE entry point is invoked and no emergent spectrum is emitted.

Background absorption is the static native/bound-Fe background with saved O
bound-free populations, no O bound-bound or added trace-line wrapper. The
41modeledFeI comparison is a natural-only recalculation, not the current model's
full pressure/Stark profile. Line/background coefficient ratios do not predict
observed flux, equivalent widths, continuum corrections or an abundance change.

At the reference native depth (index41), actual logtau500 is-0.0099128 for the
canonical case and-0.0089528 for lower gravity. Absolute-weight footprint means:

| Frozen case | Synthesis-anchor footprint | All other screened natural lines / background absorption | FeII portion |
|---|---|---:|---:|
|canonical|7771|2.176416e-7|1.636470e-8|
|canonical|7774|1.867715e-7|8.065700e-9|
|canonical|7775|0.002250493|0.002247655|
|lower-g|7771|2.252013e-7|1.760178e-8|
|lower-g|7774|1.873699e-7|8.680473e-9|
|lower-g|7775|0.002375557|0.002372575|

The eligible natural-only additions are very small at the7771/7774 synthesis
anchor locations. This cannot exclude omitted pressure wings,42unqualified
natural-width rows,61unbound-upper rows, other withheld candidates, molecules,
continuum/normalization differences or source-model errors. It is not a general
upper bound on missing opacity.

FeII gf_079010 near7777.12 is now a strength-supported transfer-test candidate,
beyond the previous centre-only inventory. In the canonical7775synthesis footprint,
its maximum local ratio at index41 is0.02242763; the all-depth maximum is0.4779355
at logtau500=-4.9152276, near observer AIR7777.119802. The large outer-column
value is not a47.8% emergent flux effect. Across footprint means, FeII dominates
the screened additions there. Local natural widths alone do not qualify its
pressure profile, LTE population physics or actual effect on the O7775 continuum.

## Source and main integration

Main advanced from f8c723a8 to9a74d843 (PR607, RYA-1311 with RYA-1315 updates).
Inspected changed bootstrap/science/decision/provider paths: no change to the
applicable AGENTS/LEDGERS/DEV_CYCLE/DATA_SAFETY/SCIENCE_STANDARDS/decision files
or Bride kernels in that range. Main changed canonical_gf and the model-guided
local-continuum bound. This branch preserves historical sources rather than
silently updating an established differential; the changed TiI row is withheld.
Historical side refusals remain historical and must not be asserted as today's
15percent-rule behavior. No continuum fit or new refusal under that rule is
claimed here. Current production truth stays on main, not this frozen experiment.

The SEQUENCE hook initially reports the newly merged RYA-1311 landing absent
from this review branch's sequence ledger. Backfill only that verified landing,
advance the SEQUENCE stamp to9a74d843, and record the source/rule integration
caveat in registerv215. Do not merge/update raw science inputs, regenerate other
elements, or claim independent main publication approval from these diagnostics.

## Commands, resources and validation

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /usr/bin/nice -n 19 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/screen_bride_lte_opacity.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-lte-opacity-screen-v1
/Users/ryanschmitt/codex/venv/bin/python scripts/analyze_bride_lte_screen.py --source /Users/ryanschmitt/codex/bride_mac_checks/rya1281-lte-opacity-screen-v1 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-lte-opacity-analysis-v1
```

Both final commands exit0. The sandbox denies nice's setpriority request;
computation continues, actual nice level unmeasured. Do not claim nice19 applied.
The one-thread coefficient computation takes0.856615s of its180s allowance;
process-lifetime peak RSS665387008bytes, host load3.39/3.01/2.44. Timing excludes
imports/input intake/output writing and is not a transfer/grid/3D compute budget.
Existing Lightweaver config warning retained; no benchmark/install triggered.
A read-only source introspection initially used the wrong internal class name;
corrected to KuruczPf before analysis, no source file altered.

Immutable registered resources experiments.bride_lte_screen and
experiments.bride_lte_analysis retain the originals; reruns require fresh directories.
Four derived result artifacts/one visually inspected figure copied byte-for-byte
under data/results/rya1281/lte_opacity_screen_v1 and
 data/output/rya1281/lte_opacity_screen_v1. Largest result2350328bytes, below5MiB.
No external raw spectral flux arrays copied. Source/report/partition/atom/structure/
binding and saved-population identities are recorded. Actual writer invocations
registered in GENERATORS.yaml.

[Full-column opacity ratios](../../../data/output/rya1281/lte_opacity_screen_v1/depth_ratios.png).

Focused suite4passed in1.00s. Expanded science-env suite75passed in7.31s with
292existing printUpdate warnings: previous71test command plus
 tests/test_bride_lte_screen_rya1281.py and tests/test_bride_lte_blend_rya1281.py.
Independent Saha and species closure tests, partition-domain refusal, ambiguous
identity/gf/changed-main/width refusal, frequency-integrated line-strength sum rule,
Kirchhoff relation and zero-opacity control pass. No new bug or silent fallback
discovered; prior characterized/fixed ledger rows retained, no empty placeholder.

Self-review: complete source/certificate/population identity checks; changed-main
row withheld, no competing gf selected. All388rows accounted for; absences are
explicit refusals, not zero coefficients. Native composition/ne/geometry and saved
O fixed, no charge/SE/RT feedback; shifted anchors map back to rest AIR then vacuum.
All depth data retained; absolute prediction weights describe geometry only.
Local LTE/natural-only approximations remain explicit and cannot clear physical,
continuum, EOS, native FTS, independent-O, 3D or all-band admission gates.

Next bounded task: qualify a same-transfer FeII gf_079010 differential with the
existing FeI/background/O held fixed, a zero-extra-line identity control, unchanged
radiation/cap accuracy gates, both continuum references and RV controls retained.
No fitted pressure width/gf/abundance, and no expectation that this solves7771/7774.
Then address source-qualified omitted wings/molecules/continuum and independent
O atom/EOS/reference. Broader xi/T/g/composition interactions, additional O lines,
C/N then Fe then reference star and adaptive3D/all-band goals remain OPEN.
Scientific HOLD; RYA-1281 In Progress / RYA-1278 In Review. No new opacity installed,
native/default/production change, adopted result, merge or publication. Mac only.

Final checks: MIRROR targets/v215 pointer/whitespace and reconciled main landing pass. Generator guard retains exactly2436existing diagnostics, zero added/removed versus prior influence audit. Fresh library audit194documents/127accounted/67unindexed/12off-Mac registered paths. Existing historical RYA-587 ledger and gold/tracker debts remain visible; no green CI claim. No newly discovered bug or empty ledger row.
