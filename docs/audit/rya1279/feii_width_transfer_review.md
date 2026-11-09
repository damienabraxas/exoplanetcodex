# FeII source-width transfer: predeclared control

2026-11-07; RYA-1281 / RYA-1278, draft PR598, scientific HOLD.
Worktree `/Users/ryanschmitt/codex/rya1279`, branch
`ryandamienschmitt/rya-1279-1280-bride-oxygen-gates`, stacked base
`4635c7421a7dee8671ffefe05d771c68da1cd24e`.
Starting evidence HEAD `e7cf9319eaadb29e5aa61810748d055a0f97c547`.
Refreshed main advanced from9a74d843 to459f5614; changes concern CNO/Na
products, Asplund atmosphere and uncertainty/publication bookkeeping. No changed
Bride transfer kernel, canonical target row or science authority. No merge into
this frozen control. Main's15percent continuum rule remains unexercised.

## Scope and acceptance

Plan `config/experiments/rya1281_feii_width_transfer.json` declares two columns
and two independent Stark alternatives (export-Cowley and supplied). Each pair
first replays the archived natural-FeII transfer and must meet the unchanged
1e-6 absolute normalized-flux guard before its physical-width leg. Eight serial
legs total,180s per leg,1e-8 radiation tolerance,1000iterations maximum. No
new coupled O/charge equilibrium or fresh physical-rate certificate claimed.

Widths use live depth-dependent H0/He0/ne/T and the verified W branch, natural
gamma plus pressure plus the explicit ionic Stark alternative. H2 remains zero.
No mean-width convolution, independent-spectrum multiplication or fitted width.
Source native hashes/commit, canonical line/strength eligibility/current main,
saved O, continuum, inherited41FeI opacity, atmosphere and numerical geometry
are checked. Cache refuses changed collider or atmosphere fields. AIR comparison
and internal vacuum convention remain unchanged. Both continuum references and
both fixed RV conventions must be retained in postanalysis.

Mac only, one thread, installed LW0.16.1; no Sirius or CI. Uptime before work
showed load2.71/2.21/2.14; sandbox denied process enumeration. No process killed.
Fresh output directories preserve old runs. No reference-target tuning, opacity
adoption, abundance change, publication or merge.

## Initial verification and debt

`source config/hosts/mac.env` then
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q tests/test_bride_feii_width_control_rya1281.py tests/test_bride_feii_control_rya1281.py`:
5passed,3existing dependency warnings,2.86s, exit0. Real-backend tests independently
form both total-width equations, verify depth-dependent LTE FeII absorption and
thermal emission, unchanged scattering/outside support, deepcopy, natural/empty
controls and refusal of changed H populations, identities and unsupported modes.

Refinement-debt report exit0:10not-fired and5TBD remain; N I has RYA-697 and
unfiled engine-B model-atom debt. None cleared by this oxygen control. Existing
bug-ledger dispositions retained; no new defect discovered in implementation.

## Execution / handoff

Predeclare and commit the provider, tests, plan and runner before real transfer.
Runner takes `--stark-mode cowley|provided --case canonical_refined|lower_gravity_refined
--output <fresh-directory>`. Analysis, numerical results and final source/ledger
checks were still owed at that pre-run checkpoint (commit3a91383f).

## Completed transfer and interpretation

All eight legs pass. In each pair `zero` means **natural FeII baseline**, and
`feii` means **natural+W pressure+the named ionic Stark alternative**; these
historical directory names do not mean an absent line in this experiment.
All four natural normalized spectra differ from the archived natural run by
exactly0; their CSV receipts match too. Saved O populations and standalone
continuum arrays remain exact. Radiation/cap convergence uses65steps canonical
and67lower-g, final J8.2973e-9/8.6213e-9, boundary1.1184e-9/1.5218e-9.

| Alternative | Canonical max wider flux change | Lower-g max wider flux change | Max direct O area change, both RVs (mA) |
|---|---:|---:|---:|
| W + export-Cowley |0.000440860|0.000425736|0.000135702|
| W + supplied Stark |0.000443279|0.000428169|0.000136473|

This is the change from **natural-only FeII**, not from the old absent-line
control. The two Stark alternatives differ in intrinsic flux by at most
2.419e-6canonical/2.433e-6lower-g. No alternative selected for better agreement.
Maximum intrinsic O-window flux change3.021e-7; direct physical absorption is
tiny in the O windows despite the larger local damping rates. These integrals
are blended fixed-window absorption diagnostics, not isolated oxygen EWs.

Synthesis-guided observed7775area: canonical48.79795350 ->48.66086098mA
(Cowley,-0.13709252) or48.66088552mA(supplied,-0.13706798);
lower-g48.91159408 ->48.91528845(+0.00369437) or48.91531214(+0.00371806).
IAG-reference observed7775area remains57.70338185mA for both columns and all
legs; model residuals still include the direct model change. The remaining
IAG-versus-synthesis spread is about9.0425canonical/8.7881lower-g mA. No atlas
or gravity is preferred from the agreement.7771/7774guided responses are tiny.

All48 native continuum-record/area replays pass (2alternatives x2columns
x2legs x2references x3lines). Canonical7775synthesis selection changes at4
positions and coefficient supports differ at11; actual support83->90,
blue/red29/54->36/54. Lower-g selection and coefficient support86(30/56)
remain exact. All IAG supports and7771/7774supports unchanged. Both denominator
values and selection/clipping can change; no causal selector-only fraction is
claimed without crossed fixed-support accounting. The earlier7774reported70
versus actual73 fifth-pass discrepancy remains characterized, not silently fixed.

The four response figures were visually inspected. They show slight FeII core
shallowing and separate direct-area/continuum-placement effects. They do not
qualify an isolated observed FeII profile; the previously documented structured
atlas shoulders, RU gf uncertainty and missing neighboring opacity remain open.
No new population equilibrium, physical reference parity, full-opacity spectrum,
solar abundance, full3D or all-band result is established.

## Reproduction and resources

Resolve `experiments.bride_feii_width_transfer`, `experiments.bride_feii_width_analysis`
and `experiments.bride_feii_width_anchors` via the Mac host profile. Full immutable
outputs and logs remain there. For each Stark alternative and each declared case:

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/bride_rya1118/.venv/bin/python scripts/run_bride_feii_width_transfer.py --stark-mode cowley --case canonical_refined --output <fresh-transfer>/cowley/canonical_refined
```

Repeat with lower_gravity_refined and with provided. Then for each alternative:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/analyze_bride_feii_width_transfer.py --stark-mode cowley --source <transfer>/cowley --output <fresh-analysis>/cowley
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/trace_bride_feii_anchors.py --transfer <transfer>/cowley --analysis <analysis>/cowley --output <fresh-anchors>/cowley
```

All four runner, two analysis and two tracing commands exit0. Source-selection,
native commit/source, unchanged atom and inherited-kernel checks pass at runtime.
Per-leg17.19-18.71s; pair totals36.47-37.61s, peak process1.284GB. Every solver
uses one thread. A scheduling mistake briefly overlapped the end of lower-g
Cowley with the start of canonical supplied; do not claim fully serial execution.
No Sirius, nice-priority, benchmark, installation/download or CI work.

46results/4figures byte-checked into `data/results/rya1281/feii_width_transfer_v1`
and `data/output/rya1281/feii_width_transfer_v1`; largest result520712bytes.
Actual writers/invocations registered; helper plots/logs stay external. No raw
source spectral copy, old-output overwrite, binary rename or guard exemption.

## Validation, review and exact next step

Focused full regression suite84passed,295existing backend warnings,8.95s,
exit0, one-thread held Bride interpreter. Exact output log:
`/private/tmp/rya1281-feii-width-tests.log`. It includes the previous81test suite
listed in `feii_profile_source_review.md` plus all3new width-provider tests.
Library audit exit1:195documents/127accounted/68unindexed,12off-Mac declarations;
no new science claim depends on an unidentified document. Existing global
run-bug-ledger RYA-587 blank science_impact remains unrelated; no new scientific
defect or silent fallback discovered, prior ledger rows retained.

Self-review checked line/stage identity, ionic versus neutral Stark equations,
m^-3 to cm^-3 density boundary, live neutral fields and all-depth caches,
same-transfer thermal emission, scattering, source receipts, both references/RVs,
fixed O/continuum and AIR/vacuum boundary. No new blocking finding. Scientific
admission remains HOLD. Registerv218 and SEQUENCE record these bounded findings
and the fresh main459f5614/PR608 landing, without merging main into the control.

Final artifact guards: all2436generator diagnostic strings exactly match the
previous phase, zero added/removed (exit1 retained). MIRROR targets, version pointer,
main-landing reconciliation and whitespace checks pass. Initial reconciliation
correctly flagged missing RYA-1315/PR608 before the new landing entry/stamp.
Global ledger tests1failed/1passed: historical RYA-587 science_impact blank; the
first invalid row matches current main459f5614 in every field. An initial typo
`tests/test_run_bug_ledger.py` gave exit4/no tests before the correct
`tests/test_run_bug_ledger_rya926.py` invocation, not an extra passing test.
No unrelated ledger repair or bypass.

Next: cross the frozen natural/source-width reference denominators with fixed
selection and fixed coefficient supports to quantify the7775continuum response,
without moving shoulders or retuning. Qualify missing neighboring opacity, then
test coupled O feedback with independently recomputed physical R+C residuals.
O source/EOS/reference closure and additional oxygen indicators precede C/N,
then Fe abundance and a reference star. Broader parameter interactions and actual
adaptive3D/all-band grids remain open. No default, product, abundance, adoption,
production, merge or publication changed; child In Progress, parent In Review.
