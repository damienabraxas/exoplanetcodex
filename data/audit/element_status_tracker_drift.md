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
| 7 | **Ca I** | 6.324, Δ+0.024 | **blank** — ~~"no independent-gf line survives the graded cull"~~ **REFUTED, see below** | see note below |
| 8 | **Ti I** | 5.471, Δ+0.501 | **blank** — same | see note below |
| 9 | **Ni I** | 6.946, Δ+0.746 | **blank** — same | see note below |
| 10 | **Na I** | 6.264, Δ+0.024 | **blank** — same | see note below |
| 11 | **Al I** | 7.406, Δ+0.976 | **blank** — same | see note below |

**Rows 7–11 — ADJUDICATED by RYA-596 (2026-07-27). The cause recorded above was WRONG.**

The original read here — "the RYA-398 graded-gf firewall culls every line in those pools" — was
**quoting the verdict's own channel text, which was itself a hardcoded, unverified claim.** It is
refuted on both legs:

1. **The firewall is not over-culling, and not culling to zero.** A fresh grade-restricted cull on
   the live pool leaves survivors for every one of the five:

   | element | pool | kept | gf tiers | cull reasons |
   |---|---|---|---|---|
   | Ca I | 28 | **3** | LOW 16 / MED 10 / HIGH 2 | SAT 13, GRADE 16, HIERR 3 |
   | Ti I | 81 | **10** | MED 46 / LOW 35 | SAT 44, GRADE 35, HIERR 25 |
   | Ni I | 26 | **2** | LOW 20 / MED 6 | GRADE 20, SAT 10, HIERR 6, BLEND 3 |
   | Na I | 4 | **2** | **HIGH 4** | SAT 2, BLEND 2 — **GRADE culls ZERO Na lines** |
   | Al I | 2 | **1** | LOW 1 / MED 1 | GRADE 1 |

   Na is the cleanest refutation: all four Na lines carry a NIST **HIGH** grade, so the RYA-398
   grade cull removes none of them. The firewall could not have blanked Na.

2. **The real cause is the RYA-522 gold tiering — and it is by ratified design, not a bug.**
   `build_solar_reference_v2_rya522.py` freezes a value only for the `gold` / `gf_floor` /
   `upper_limit` tiers; everything else is `owed` and freezes **no value even when one was
   produced** (`a_frozen = a_verdict if conf != "owed" else None`) — Ryan's 2026-07-05
   ratification, "suspect → held, not immortalised". All five sit in `owed`. Gold **v1** still
   carries their curated numbers (Al 7.406, Ca 6.324, Na 6.264, Ni 6.946, Ti 5.471 — exactly the
   seed values), proving the curation produced them.

3. **The bug was the round-trip.** `phase_c_verdict_rya371.py` reads the frozen gold back in
   (RYA-469), sees a blank `A_X`, and fell through to a branch that asserted the firewall cause
   **without checking it**. So a deliberate hold re-entered the verdict disguised as a cull. The
   tell was on the row all along: `n_lines` = 2/10/2/2/1 sitting next to "no line survives".
   **Sr I (row 1) carried the same phantom string** for the same reason.

**Fixed in RYA-596:** the classifier now distinguishes *held at tier* from *zero survivors*, and a
generation-time tripwire (`_assert_blank_cause_is_honest`) refuses to emit the zero-survivor claim
on any row with `n_lines > 0`. Verdict counts are unchanged (5/0/21/0) — this corrects the stated
cause, it does not invent a value.

**Still genuinely zero-survivor (claim verified, text retained):** Mg (0/5), Y (0/3), Zr (0/6),
Eu (0/1). Caveat worth its own ticket: for **Zr and Eu the blank is NOT gf-grade-driven** — all
their lines are MED tier and die on SAT/BLEND/HIERR (the RYA-395 quality cuts), so the surviving
"pool gf is Kurucz/ungraded" wording is right for Mg/Y but imprecise for Zr/Eu.

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

---

## E. Mg I — a REAL measurement disagreement, recorded not reconciled (RYA-592, 2026-07-27)

Reads together with RYA-596 (§B, Mg row): RYA-596 verified that Mg's EW pool is a **genuine**
zero-survivor state (0 of 5 lines kept) — unlike Ca/Ti/Ni/Na/Al, Mg's blank really is a cull, not
a gold-tier hold. That is precisely why Mg has to be measured by synthesis, and what RYA-592 did.

RYA-592 measured the second clean Mg I line (5528.405) by in-window blend-fit synthesis on
Sirius, and re-measured 5711.088 **with the same harness** so the comparison is like-for-like.
Three numbers for the same element now disagree, and per the corollary-1 rule none of them is
silently overwritten:

Per-line Engine-B values (all include the Gerber δ = −0.023), so like is compared with like:

| source | line | A(Mg) | red-χ² | method |
|---|---|---|---|---|
| committed channel | 5711.088 | **7.494** | 4.85 | synth-v2 flux fit |
| committed channel | 6319.237 | **7.734** | **58.5** | synth-v2 flux fit |
| committed channel | **reported = mean of the two** | **7.614** | — | two-engine floor, n=2 |
| RYA-592 harness | 5711.088 | **7.400 / 7.414** (HARPS / IAG) | 1.28 / 1.48 | in-window blend-fit, rest-frame |
| RYA-592 harness | 5528.405 | **7.176 / 7.180** (HARPS / IAG) | 1.37 / 1.59 | in-window blend-fit, rest-frame |

The headline is **line-to-line scatter of ~0.2 dex in Mg, present in BOTH methods**:

1. **5528 vs 5711, RYA-592 harness** — 0.207–0.234 dex, well outside the 0.10 band, and
   *consistent across two independent solar arms* (HARPS and the IAG FTS atlas agree to
   ≤0.015 dex on each line). A line-to-line disagreement, not an arm artefact.
2. **6319.237 vs 5711, committed channel** — **0.240 dex**, the same size and the same
   direction of problem, already present inside the value the verdict carries today.
3. **Same-line, method-to-method (5711)** — only **0.080–0.094 dex**, i.e. the two
   independent profile-fit implementations agree far better with each other than either
   agrees across lines. That localises the problem to the LINES, not to either harness.

Point 3 also means the committed 7.614 is not a robust central value: it is the mean of a
reasonable fit (5711, red-χ² 4.85) and a poor one (6319.237, red-χ² **58.5**) — and 6319.237
is the strong 6318/6319 complex sitting on the Ca I autoionisation feature, i.e. a line the
RYA-592 ticket itself lists among Mg's unusable measurements. It contributes half the
reported value and pulls it **up** by +0.12 dex.

What is already ruled out as the cause (all tested, see the ticket comment):

- **NLTE atom** — both engines reproduce their own committed 5711 δ (Engine-A residual
  −0.0000 vs the committed grid; Engine-B −0.0008 vs the RYA-534 record), and the two atoms
  agree on 5528 to 0.007 dex.
- **gf** — canonical, NIST-verified, and identical in all three stores for both lines.
- **Rest frame** — the solar arms carry a real +0.76 (HARPS) / +0.28 (IAG) km/s residual
  velocity (RYA-309). Fitting it as a nuisance parameter cut red-χ² ~4× (5.5→1.4) but moved
  A by only 0.003–0.005 dex. Not the cause.
- **van der Waals damping / broadening grid** — 5528 is insensitive (7.217 vs 7.220 across
  the VALD and GES-ABO values); 5711 *is* damping-sensitive (7.37 ABO → 7.43 VALD → railed
  at 8.05 with Unsold), which is a lead for (2) but does not explain (1).

**Live suspicion (not yet demonstrated):** 5528 is a strongly saturated line (observed EW
≈340 mÅ, 3.4× the ratified 100 mÅ saturation knee), so it sits on the damping/microturbulence
part of the curve of growth where 1D-LTE is weakest — exactly the regime standard solar Mg
analyses avoid. A ξ/damping/3D systematic on the strong line is the leading explanation, which
would mean 5528 is *measurable but not trustworthy as an abundance indicator*, distinct from
"not measurable" (it comfortably clears the reliability floor at dEW/dA = 130 mÅ/dex).

**Consequence:** Mg does NOT promote. The second line contradicts rather than confirms the
first, and RYA-561 gate 3 fails independently anyway (5528 is saturated → Engine-B by the
ratified clause-3 classifier → dCE stays `None`). Mg stays `CURATION-OWED-with-value`.

**Owed, before the v3 freeze carries Mg's +0.064:**

1. Adjudicate the ~0.2 dex Mg line-to-line scatter (which lines are admissible at all).
2. Re-examine whether **6319.237** belongs in the reported value: red-χ² 58.5 on a known
   bad blend, contributing half of 7.614 and +0.12 dex of it. On 5711 alone the committed
   channel gives 7.494 (+0.056 →) and the RYA-592 harness 7.400 (−0.15).
3. Note the *sign* is at stake, not just the magnitude: drop 6319 and/or adopt the RYA-592
   scale and Mg's offset moves from **+0.064** to somewhere in **−0.06 … −0.37**. Mg's
   "+0.064, comfortably inside the band" framing does not survive either change, so this is
   not a cosmetic difference.
