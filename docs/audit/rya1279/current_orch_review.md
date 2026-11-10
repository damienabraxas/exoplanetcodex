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

Current checkout's native continuum tests:8passed in0.08s, exit0 (project Mac
interpreter,1thread). Exact command: `python -m pytest -q tests/test_rya1230_local_continuum.py`.
Precommit plan/runner before observational evaluation. No new RT, source/product
rewrite, adoption/production/merge/publication. Scientific HOLD.

Initial harness v1 refused before any output directory/observational fit: synthetic canary omitted the required exclude_half_width_A keyword. Corrected both calls to explicit0.5A; native/source data unchanged. Fresh v2 output/log follows. No accepted v1 result.

## Accepted result

Code/plan cce5ebde; API-call fix2cd27ce6 before accepted evaluation.
Accepted v2 command, exit0, Mac project interpreter, one thread:

```sh
source config/hosts/mac.env
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_current_orch.py --orch-worktree /Users/ryanschmitt/codex/rya1281-orch-current-de96ae51 --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-current-orch-v2
```

All48historical native records/areas replay exactly against saved traces. All48
current records (including method/reason), entire applied continuum arrays and
observed areas exactly match historical. No current measurement changed. Current
levels span0.9946038057630487 to1.0045515524639725; all lie within the former
3percent guard, so the15percent model bound is inactive on these frozen cases.
The previous support/clipping response therefore remains present with the actual
current ORCH estimator on these inputs; the1315HARPS blue-sag fix does not close
the O777reference discrepancy here. This says nothing about untested bands or
the broader1315HARPS production changes.

| Synthetic control | Historical model-guided | Current model-guided | Current envelope |
|---|---|---|---|
| Uniform0.92level | Refused | Applied | Refused |
| Uniform0gap | Refused | Refused | Not tested |

Both module thresholds checked directly (.03envelope/.15model); constants were
not edited. Native loaded source bytes equal pinned Git blob and1315branch blob.
Same frozen source spectra/observations and reference arrays, AIR comparison,
internal vacuum spectra, plus-gravity continuum reference, no new telluric or
conditioning operation. Prior transfer still retains both RVs. Current eight
native tests passed0.08s; observational comparison is additional runtime evidence,
not48independent physical validations. No new RT or coupled O/R+C claim.

## Artifacts, QA and handoff

Registered checkout `experiments.bride_current_orch_checkout`; accepted output
`experiments.bride_current_orch`. Two result copies at
`data/results/rya1281/current_orch_v1` byte-match external accepted v2; largest
54542bytes. Actual writer/invocation registered. Original v1 failure log
`/private/tmp/rya1281-current-orch.log` and v2log
`/private/tmp/rya1281-current-orch-v2.log` retained, no output overwrite. No binary
or raw spectral copy. Bug-ledger row `rya1281-current-orch-canary-signature`
records the fixed harness API error; no accepted science affected.

Fresh library audit retains195documents/127accounted/68unindexed/12off-Mac
declarations, exit1. No newly unidentified literature used. No production or
publication, adopted continuum/abundance, source/default/product change, merge,
Sirius, CI, benchmark, installation/download or nice-priority claim. All owned
jobs finished. Complete current ORCH checkout remains available for subsequent
checks; historical experiment code/results were not replaced by a mass merge.

Self-review checked actual source identity/branch ancestry, dirty-source refusal,
both distinct native modules, active-bound canaries, exact records/arrays/areas,
all48identities and explicit plus-gravity scope. Registerv220/SEQUENCE record the
bounded finding and scientific HOLD. Existing ledger/tracker/gold and refinement
debts remain; the new harness defect is fixed, not hidden as a science result.

Final guards: generator check exit1, exactly2436diagnostic strings identical to
preceding phase, zero added/removed; no new invocation/registration failure.
MIRROR/version/SEQUENCE/whitespace pass. Global ledger tests retain historical
1failed/1passed (RYA587blank science_impact); that first invalid row matches current
main in every field, and the newly appended harness row has the required fields.
Postcommit freshness and normal pre-push gates run before handoff, no bypass.

Next: qualify missing neighboring opacity and continuum information, then compare
coupled O feedback using the fixed natural/source-width controls and independently
computed physical R+C/conservation diagnostics. Additional oxygen indicators,
C/N then Fe/reference star, wider parameter interactions and actual adaptive3D/
all-band goals remain open. RYA1281 In Progress, RYA1278 In Review, draftPR598.
