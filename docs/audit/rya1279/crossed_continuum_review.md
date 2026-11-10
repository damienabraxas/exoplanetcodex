# Crossed continuum controls

2026-11-07, RYA-1281 / RYA-1278, draft PR598, HOLD.
Worktree `/Users/ryanschmitt/codex/rya1279`, branch
`ryandamienschmitt/rya-1279-1280-bride-oxygen-gates`; stacked base4635c742,
starting evidence3f0fb5ab, main459f5614. No unrelated worktree changes.

## Predeclared design

`config/experiments/rya1281_crossed_continuum.json`: both FeII Stark alternatives,
both columns, all three O lines, synthesis and IAG references. For each set,
cross the before/after denominator with before/after support in three designs:
original selection with clipping; original selection without clipping; actual
native coefficient support without clipping. All four corners retained. Final
coefficient support means the pixels used by the polynomial, not the postfit
survivor count (known7774actual73/reported70 distinction remains intact).

Native-clipping diagonal endpoints and final-coefficient diagonal endpoints
must reproduce archived level, slope, indices and observed area **exactly**.
Off-diagonal values are labelled counterfactual diagnostics, never adopted
continua. Finite path averages, interactions and closure are arithmetic accounting
of correlated changes, not physical causal fractions or independent uncertainties.
No new RT, O populations, widths/gf/RV fits, moved windows, reference preference,
source/observed conditioning, product/default/adoption/production/merge/publication.
Historical3percent rule retained by hash; main15percent rule not exercised.
The matched continuum audit uses the existing plus-gravity reference; both RV
variants remain in the preceding frozen transfer analysis, not newly refitted.

## Initial validation

Mac held Bride interpreter, one thread:
`python -m pytest -q tests/test_bride_crossed_continuum_rya1281.py tests/test_bride_continuum_audit_rya1281.py`
10passed in0.14s, exit0. Tests cover finite paths/interaction/closure/refusal,
analytic denominator effects, exact fixed supports, immutable arrays and a
discriminating clipped-versus-fixed control. Initial synthetic outlier fixture
triggered the native minimum-anchor refusal rather than clipping as intended
(1failed/9passed); sparse outliers corrected the test setup, no native-rule change.
No scientific implementation defect or silent fallback discovered.

Precommit plan/kernel/tests/runner before evaluation. Exact next command:

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_crossed_continuum.py --output <fresh-registered-crossed-directory>
```

Results, artifact registrations, review and final guards remain owed at this checkpoint.

## Evaluated controls

Predeclared code/plan commit5c7ee2fd. Audit exits0:72crossed2x2designs,
288fits,96exact native/coefficient endpoint replays, no refused corner and every
finite closure exactly0. The matching level, slope, actual coefficient indices
and integrated observed areas reproduce source traces, including the7774cap
support distinction. IAG holding/reader and delivered Kitt product hashes checked;
all source transfer, analysis and trace receipts checked. No raw-data copy.

For synthesis-guided7775observed area (mA):

| Column / Stark | Native total | Support path average | Denominator path average | Interaction |
|---|---:|---:|---:|---:|
| Canonical / Cowley |-0.137092523|-0.141296433|+0.004203910|-0.000462553|
| Canonical / supplied |-0.137067985|-0.141297849|+0.004229865|-0.000465385|
| Lower-g / Cowley |+0.003694366|0|+0.003694366|0|
| Lower-g / supplied |+0.003718064|0|+0.003718064|0|

Canonical Cowley native crossed corners:

| Support source | Natural denominator | W+Cowley denominator |
|---|---:|---:|
| Natural selection, native clipping (83coefficient pixels) |48.797953500|48.802388687|
| Width selection, native clipping (90coefficient pixels) |48.656888343|48.660860977|

Holding the actual coefficient support fixed with clipping disabled gives
**identical four corners** here: changing the denominator alone does not change
which of those pixels survive in this experiment. Changing the initial selection
(two replacements out of90; symmetric difference4) changes the clipped fit from
83to90pixels. The83-to90coefficient sets have symmetric difference11. This is
the localized estimator response observed here, not a claim about every spectrum
or the current main15percent estimator.

Using original selected pixels **without clipping** instead gives canonical
Cowley total+0.013424557mA, support path+0.009262941 and denominator+0.004161616;
interaction-0.000377964. Supplied total+0.013449095. For Cowley, the unclipped
natural endpoint is48.647436420 rather than the native48.797953500mA; the
width endpoint48.660860977is unchanged. Thus the difference between these two
diagnostic totals is-0.150517080mA. It is not an independent physical clipping
error or a recommendation to disable clipping. The observed data and original
fits remain unchanged.

Lower-g has identical initial and coefficient supports between the two physical
legs: support term0, a small reference-value response remains. Its no-clipping
totals are+0.003683091Cowley/+0.003706717supplied; native/fixed-coefficient totals
remain in the table. All IAG terms are exactly0 in all designs.7771/7774total
responses across all designs remain<=0.000016522mA.

These experiments separate the finite effects of reference denominator and
support choices. They do not establish physical uncertainties, an optimum mask,
the true continuum, a preferred atlas/gravity or solar abundance. No favorable
corner replaces the original measurement. The much larger8.8-9.0mA reference
spread, one-sided7771/7774anchors, structured FeII shoulders, missing neighboring
opacity and O/EOS/reference gates remain open. No new RT/SE/coupled-O result.

## Resources, validation and self-review

Registered Mac resources `experiments.bride_crossed_continuum` and
`experiments.bride_crossed_continuum_plot` retain immutable originals. Audit log
`/private/tmp/rya1281-crossed-continuum.log`. Plot command, exit0:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/plot_bride_crossed_continuum.py --source <registered-crossed-directory> --output <fresh-registered-crossed-plot>
```

3results and1visually inspected figure byte-match their originals at
`data/results/rya1281/crossed_continuum_v1` and
`data/output/rya1281/crossed_continuum_v1`; largest result621193bytes. Actual
writer invocations registered. No old output overwritten; all owned jobs finished.

Full held Bride regression87passed,295existing dependency warnings,9.02s,
exit0; log `/private/tmp/rya1281-crossed-tests.log`. Exact suite is the preceding
84test suite in `feii_width_transfer_review.md` plus all3new crossed tests.
All science/analysis on Mac, one thread; no Sirius, benchmark, installs/downloads,
CI, nice-priority or new RT. Fresh library audit exit1 retains195documents,
127accounted/68unindexed and12off-Mac declarations. None newly used as evidence.
Refinement debt remains10not-fired/5TBD per preceding session report, including
N RYA697/unfiled engine-B atom; no debt cleared by this continuum diagnostic.

Main advanced during checks to de96ae51(PR609/RYA1315): Fe/Ti product feeds,
derived artifacts, bibliography and publication bookkeeping. Changed code is
element_publish/Fe-schema emitter and associated product tests; no changed Bride
kernel, continuum rule or authority in that range. Reconciled SEQUENCE landing,
without merging production changes into historical controls. Statev219 records
the finite finding, not a new scientific admission.

Self-review checked support identity versus reported survivors, unclipped versus
native semantics, before/after identity, all four retained corners, exact replay,
finite closure, correlated-path interpretation, zero-change IAG controls, source
hashes, AIR/vacuum and plus-gravity scope. No new bug or silent fallback found;
prior bug-ledger dispositions retained. The first synthetic fixture failure was
a test setup issue; no physical defect attributed to the existing estimator.

Final guards: generator check exit1 with exactly2436existing diagnostic strings,
zero added or removed versus the preceding phase; no new invocation/registration
failure. MIRROR targets, version pointer, main-landing reconciliation and whitespace
pass. Historical global-ledger1failed/1passed from the preceding checkpoint remains
documented (RYA587blank science_impact; no ledger changes on newer main). Postcommit
freshness and pre-push checks run before handoff; no guard bypass.

## Next bounded work

Replay the same frozen spectra/observations through the actual current-main ORCH
continuum implementation as a separately sourced control. Establish whether its
different model-guided bound changes these cases; do not silently relabel the
historical3percent results as current ORCH behavior. No main/native rewrite.

Qualify the missing neighboring atomic/molecular absorption and local continuum
information rather than selecting the numerically convenient support. Prepare
the existing natural/source-width controls for coupled O feedback, with fresh
independent physical R+C residual and conservation checks; do not equate stable
emergent flux with statistical equilibrium. Additional O diagnostics, C/N then
Fe/reference star and wider parameter/interaction/adaptive3D/all-band gates remain
open. No product/default/production/adoption/merge/publication; HOLD.
