# Element status tracker — seed-vs-live-verdict drift (RYA-594)

Companion to `data/audit/element_status_tracker.csv`.

- **Seed** = RYA-524 master 27×2 audit, Linear comment **2026-07-14** (§2 table, §3 Cr ion balance,
  §5 batched fix plan) + the Sr correction comment **2026-07-15**.
- **Live verdict** = `docs/audit/solar_phase_c_verdict_rya371.md` at `main` **809542b**
  (file last regenerated **2026-07-18**, by RYA-559 `564824a`).
  Re-verified byte-identical at `main` **35f4359** (RYA-561 PR#164 landed a guard helper only and
  did not regenerate the verdict), so every comparison below still holds at that SHA.

Per RYA-594 §5: divergences are **flagged, not silently reconciled**. The CSV carries the audit
snapshot; nothing below has been folded into it.

---

## A. Seed-internal inconsistency (found while seeding — flag, do not fix)

The RYA-524 comment's **§4 summary tally** does not match its own **§2 per-row table**.

| bucket | §4 summary claims | §2 per-row table gives |
|---|---:|---:|
| WIRED-OK | 9 | **8** |
| GENUINELY-OWED | 13 | **14** |
| DONE-BUT-STALE | 2 | 2 |
| WRONG-SPECIES | 1 | 1 |
| VINTAGE-INFLATED | 1 | 1 |
| NLTE-VOID | 1 | 1 |

Cause: §4's WIRED-OK list is "Fe I, Fe II, C, O, K, Mn, **V-void**, Cu, Li" — it counts **V** in
WIRED-OK *and* again as the NLTE-VOID row, so V is double-counted and GENUINELY-OWED is short by one.
Both sum to 27 only because of the double count.

**Resolution used:** the CSV follows **§2**, the per-row table, which is the per-element authority and
the thing the tracker is a row-for-row transcription of. §4 is a summary. Nothing was edited in Linear.

---

## B. Seed vs live verdict — element by element

**AGREE (13):** C, O, Mg, Si, Co, K, P, Y, Eu, Cr, V, Sc, Cu, Zr — value, blank-ness and disposition
all match. (Vocabulary differs: the audit says `gold PASS` / `gf_floor` / `owed`, the verdict says
`PASS` / `CURATION-OWED`. Not a disagreement.)

**DISAGREE (11 rows):**

| # | row | seed (2026-07-14/15) | live verdict (main 809542b) | read |
|---|---|---|---|---|
| 1 | **Sr I** | 07-15 correction reclassifies Sr → **WIRED-OK**, A(Sr II) ≈ **2.74** (RYA-551, 4077 PRIMARY) | **still blank**, `CURATION-OWED`, "EW present; no independent-gf line survives the graded cull", NLTE grid cited as `Sr_Bergemann2012_INSPECT.csv` | **LOUDEST.** RYA-551 is merged (`4b3280e`, confirmed ancestor of `origin/main`) but the verdict has **not adopted** the Sr II synthesis. The 07-15 reclassification to WIRED-OK is **not** reflected in the live verdict, and the wrong-species Sr I channel is what the verdict still describes. Not reconciled here. |
| 2 | **N I** | 8.202, `NLTE-OWED`, Δ+0.372; action = "apply the wired ≈ **−0.37**" | 8.188, `CURATION-OWED`, Δ+0.358; NLTE wired via RYA-556, per-line δ −0.0115/−0.0145/−0.0154, **mean −0.0138** | The audit's predicted −0.37 NLTE correction is **refuted**. The NLTE debt is cleared; the surviving +0.36 is a Kitt-Peak red-multiplet **gf/data-channel floor** (RYA-161), not an NLTE debt. Seed's `action_needed` is now stale. |
| 3 | **Fe I** | 7.516, Δ+0.056 | **7.466**, Δ+0.006 | RYA-553 applied the tabulated Magic-2013 1D→3D solar Fe correction at the reported layer; `FE_GATE [7.41,7.51]` governs. The anchor moved after the seed. |
| 4 | **S I** | 7.753, Δ+0.633, `DONE-BUT-STALE`, expected ≈ 7.486 | **7.486**, Δ+0.366 | Audit prediction **confirmed** by RYA-557/492. Drift resolved in the live verdict; seed row is simply older. |
| 5 | **Ba II** | blank, `owed`, EW culled | **2.410**, Δ+0.140 (RYA-559 Ba II 5853 HFS synthesis + Korotin2015 δ −0.0285) | Live verdict is **ahead** of the seed. Live also flags the pool EW as blend-inflated by ≈ +0.15 (clean ≈ 2.19–2.23) — deblend owed. |
| 6 | **Mn I** | 5.470, Δ+0.050 | 5.466, Δ+0.046 | Small; RYA-546 ab-initio δ+0.024 re-derivation. |
| 7 | **Ca I** | 6.324, Δ+0.024 | **blank** — "no independent-gf line survives the graded cull" | see note below |
| 8 | **Ti I** | 5.471, Δ+0.501 | **blank** — same cull | see note below |
| 9 | **Ni I** | 6.946, Δ+0.746 | **blank** — same cull | see note below |
| 10 | **Na I** | 6.264, Δ+0.024 | **blank** — same cull | see note below |
| 11 | **Al I** | 7.406, Δ+0.976 | **blank** — same cull | see note below |

**Rows 7–11 (one cause):** the live verdict reports **no value** for Ca/Ti/Ni/Na/Al because the
RYA-398 graded-gf firewall — wired into the default run by RYA-456 — culls every line in those pools
(pool gf is Kurucz/ungraded). The audit's 2026-07-14 numbers are pre-cull. This is a **systematic
channel change, not five independent regressions**, and it moves those five from "owed with a number"
to "owed with no number". Direction of truth not adjudicated here.

**Structural mismatch (not counted above):**

- **Fe II** — the tracker (and the RYA-524 §2 table) carries Fe II as its own row
  (≈7.486, arbiter/diagnostic). The live verdict has **no Fe II row**: Fe is a single row annotated
  "EW: 62 Fe I + 3 Fe II". So the tracker's 27-row element/ion grain is **finer** than the verdict's
  26-row grain. Kept as-is — the ticket's canonical list is 27 and splits Fe.
- **Li I** — value agrees (0.727) but **tier does not**: seed says `upper_limit`, the live verdict has
  no upper-limit tier and files Li under `CURATION-OWED`. RYA-563 merged the guard
  (`engine_selection.is_upper_limit_disposition`); the consumer-side fix rides with the **unmerged**
  RYA-527. Tier drift is expected to close when RYA-527 lands.

---

## C. Stale-doc findings from the seed that are STILL live at 809542b

- `data/curation/nlte_grid_availability.csv` **Ti row still reads `Ti_Bergemann2011_MPIA.csv`** —
  verified at `809542b`. RYA-545 changed `constants.py` (and the verdict now cites
  `Ti_Mallinson2024_PySME.csv`) but not this wiring matrix. The audit flagged exactly this; it has not
  been fixed.
- `docs/science/problem_children_registry.md` (RYA-463) was dated 2026-07-04 at audit time — not
  re-checked in this pass.

---

## D. Live verdict counts at 809542b (for reference)

`PASS 5` (O, C, Fe, Mn, K) · `NLTE-OWED 0` · `CURATION-OWED 21` · `DATA-GAP 0`.
The seed predates the NLTE-OWED → 0 transition (N cleared by RYA-556).
