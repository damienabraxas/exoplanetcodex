# synth-v2 reproduces exactly — the harness was the only broken thing — RYA-713

Ryan: *"run `_run_synthesis_v2_mode` and confirm 7.520 still reproduces."*

| | banked (2026-08-08) | reproduction (2026-08-09) |
|---|---|---|
| A(Fe I) | 7.520 | **7.520** |
| A_X_std | 0.066 | **0.066** |
| n_lines | 23 | **23** |
| med_red_χ² | 2.976 | **2.976** |
| nlte_applied | False | **False** |

**Identical in every digit.** No regression. The validated synthesis path works, and the
eleven failing control runs were failures of my harness alone.

## What the run also shows

The fit-quality gate, doing exactly what my adapter was missing until the last iteration:

```
[synth-v2] fit-quality gate: 24 accepted / 61 rejected (χ²ᵣ ≥ 10.0) of 85 fitted
[synth-v2] REJECT(fit) Fe I 4969.917 Å — χ²ᵣ=42.9 ≥ 10.0 (a_synth=7.578 not used)
```

**61 of 85 fitted lines are rejected on merit** — fitted, recorded, and excluded, never
dropped. My control accepted every converged line, which is why its scatter ran 1.181
against this path's 0.066.

## Method note

`run(star_id='solar', engine='synthesis-v2')` could not be used: it dies on an unrelated
RYA-342 guard about Li I 6707.76 (see the Li section of this ticket). This called
`_run_synthesis_v2_mode` directly with **Fe-only linemasks**, never entering `run()`, so
the Li guard was not bypassed — it simply was not on the path.

## Consequence

The synthesis handler should **wrap this driver**, not re-derive it. Line selection
(`last_linemasks`, carried from the EW stage with per-line windows from measured EWs) and
aggregation both belong to the validated path. The handler's legitimate scope is band
policy, coverage checking and quarantine translation — nothing that touches an abundance.
