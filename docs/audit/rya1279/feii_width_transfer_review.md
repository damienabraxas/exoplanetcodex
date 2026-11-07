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
checks are still owed at this pre-run checkpoint.
