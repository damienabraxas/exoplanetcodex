# 🛑 BLOCKING — the RYA-553 Fe 1D→3D correction double-applies against gold v3

**RYA-669 (RYA-527 Phase 2). Found 2026-08-07. Reproduced on both machines.**
**This is a §4 STOP condition: "Fe I moves off 7.466 … or double-corrected 7.416 landed".**

Nothing was frozen, promoted or adopted. Gold v3 is byte-unchanged
(`sha256 47ad869e…83421`). The canonical phase_c verdict on the branch is the
**committed** one (Fe 7.466) — the regenerated file is preserved as evidence only, at
`evidence_fe_double_correction/solar_phase_c_verdict_REGENERATED_FE_7416.json`.

---

## What happens

Re-running the verdict generator on current main produces:

```
RYA-553 Fe 1D→3D: A(Fe I) 7.466 (1D-NLTE) -0.050 -> 7.416 (3D-NLTE, Magic 2013)
```

```json
"A_measured": 7.416,
"fe_1d3d_correction": {
  "applied": true, "correction_dex": -0.05,
  "a_1dnlte_pre": 7.466, "a_3dnlte_post": 7.416, "scale": "3D-NLTE"
}
```

The Magic-2013 −0.05 dex correction is applied to a value that **already carries it**.

## Why — one desynchronised cell

`scripts/phase_c_verdict_rya371.py` reads its Fe anchor from the frozen gold reference
(`read_solar_reference('CURRENT')`, RYA-469) and guards against exactly this:

```python
# Idempotent on the gold row's method_scale: applied to a 1D-NLTE anchor, SKIPPED if
# the anchor is already 3D-NLTE (so a re-frozen v3 gold is never double-corrected …)
if '3D' not in scale.upper():
    a_meas = round(a_meas + dex, 3)
```

The guard keys on the gold row's `method_scale` label. And in gold v3 that label
contradicts the value it labels:

| version | `A_X` | `method_scale` | consistent? |
|---|---|---|---|
| v2 | 7.516 | `1D-NLTE (Fe I)` | ✅ 1D value, 1D label → correction applies once → 7.466 |
| **v3** | **7.466** | **`1D-NLTE (Fe I)`** | ❌ **3D value, 1D label → correction applies again → 7.416** |

RYA-665 froze the post-correction number under the pre-correction scale label. `A_X`
moved 7.516 → 7.466; `method_scale` did not move with it. The guard RYA-553 built for
this precise scenario reads the label, sees no `3D`, and re-corrects.

**Idempotency was asserted on the wrong thing.** The correction is idempotent with
respect to a *label a separate process is responsible for maintaining*, not with respect
to the value. Any freeze that updates one without the other silently re-arms it.

## Why nothing caught it

`FE_GATE = [7.410, 7.510]`:

| A(Fe I) | in gate? |
|---|---|
| 7.466 (correct) | INSIDE — passes |
| **7.416 (double-corrected)** | **INSIDE — passes, by 0.006 dex** |
| 7.516 (uncorrected) | OUTSIDE — fails |

The gate was tightened to catch a *dropped* correction, and it does. It cannot catch a
*doubled* one: 7.416 clears the lower bound by 0.006 dex.

This was **measured, not inferred** — the double-corrected verdict was substituted for
the committed one and the gates were executed against it:

```
pytest tests/test_solar_calibration_gate.py   ->  9 passed          (RYA-166, all 9 assertions)
python pipeline/ledger_consistency_guard.py   ->  0 undocumented    (RYA-632)
```

Both green on a verdict whose Fe anchor is 0.05 dex wrong. RYA-632 is structurally
blind here by design — it compares verdicts, tiers and counts, never values.

The defect is latent on main only because nothing has re-run phase_c since the freeze
landed on 2026-08-07. phase_c is the canonical status channel (RYA-654); the next
regeneration propagates 7.416 into the element status tracker, the disposition report,
the gold builder and the differential denominator for every future star — with every
gate green.

## Reproduction

Both machines, independently, bit-identical:

```bash
python scripts/phase_c_verdict_rya371.py
python -c "import json; d=json.load(open('data/audit/cno_synthesis/solar_phase_c_verdict.json')); \
           print([v['A_measured'] for v in d['verdicts'] if v['element']=='Fe'][0])"
# -> 7.416   (committed value: 7.466)
```

Mac (`.venv312`, numpy 2.2.6) and Sirius (`/mnt/codex-data/venv312`) agree.

## What this does NOT mean

- **Not a measurement change.** No spectrum, EW, gf or grid moved. The fresh two-engine
  record is unaffected — it never reads the gold reference.
- **Not a v3 value error.** Gold v3's Fe **7.466 is the correct 3D-scale number**. The
  wrong cell is the `method_scale` label beside it.
- **Not a regression in this branch.** RYA-669 changed nothing on this path; the defect
  arrived with the RYA-665 freeze (PR #184) and has been latent on main since.

## Options — Ryan's call, nothing actioned

Gold v3 is write-once/immutable (RYA-469) and RYA-669 §5 forbids touching it, so
"correct the cell in place" is not available.

1. **Harden the guard (recommended, and independently worth doing).** Stop keying
   idempotency on a label another process must remember to update. The freeze already
   records the pre/post pair; the verdict can assert against the value it expects rather
   than trusting a string. This makes the class of defect unrepresentable instead of
   fixing one instance of it — and it is a code fix, so it needs no re-freeze.
2. **Re-freeze v4 with `method_scale` = `3D-NLTE (Fe I)`.** Corrects the record, but
   costs a ratified freeze cycle and leaves the guard just as fragile for the next one.
3. **Both** — 1 to close the class, 2 to make the frozen record self-describing.

Doing 2 alone would leave the same trap armed for the v4→v5 freeze.

**Also worth Ryan's attention:** `FE_GATE`'s lower bound cannot distinguish "correct"
from "double-corrected", because the correction (0.05) equals the gate half-width
(0.05). That is a gate-design question separate from this defect and is not something
this ticket should decide.
