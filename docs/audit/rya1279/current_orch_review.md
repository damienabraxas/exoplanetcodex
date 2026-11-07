# Latest RYA-1315 ORCH comparison

2026-11-07, RYA-1281/1278. User explicitly requested latest RYA-1315 ORCH.
Fetched origin; main de96ae51fd8c52f79e5457fc0205e4baaa4596fd contains final
1315branch8efedb8e via PR609. No ORCH pipeline/config differences between that
branch and main. Pulled complete main into detached Mac worktree
`/Users/ryanschmitt/codex/rya1281-orch-current-de96ae51`; historical Bride branch
stays isolated at `/Users/ryanschmitt/codex/rya1279`, stacked base4635c742.
Read RYA1315 description and comments; later Ti-II headline decision and final
merged evidence supersede the original Ti-I-primary brief and queued-run prose.

Predeclare `config/experiments/rya1281_current_orch.json`: load actual pinned
`pipeline/local_continuum.py`, verify checkout bytes against Git blob and latest
1315branch, then compare48native fits on frozen spectra, both columns/Stark modes,
all3O windows and synthesis/IAG references. Historical record/continuum/area replay
must remain exact. Report any current change, never assume it is zero. Synthetic
0.92 sag must be refused by historical3percent but admitted by current15percent;
gap stays refused and upper-envelope bound stays3percent. Source is imported
without editing or monkeypatching its constants. Same held reference arrays and
plus-gravity convention; prior transfer retains both RVs. This is actual current
continuum implementation, not a full ORCH engine/abundance/publication rerun.

Current checkout's native continuum tests:8passed in0.20s, exit0 (project Mac
interpreter,1thread). Exact command: `python -m pytest -q tests/test_rya1230_local_continuum.py`.
Precommit plan/runner before observational evaluation. No new RT, source/product
rewrite, adoption/production/merge/publication. Scientific HOLD.
