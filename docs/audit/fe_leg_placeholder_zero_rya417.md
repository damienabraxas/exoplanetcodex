# RYA-417 — MPIA Bergemann Fe-leg placeholder-zero sweep

Extends the RYA-413 Ca 6166 fix to the **separate primary-Fe NLTE leg**
(`_load_mpia_fe_grid` / `apply_fe_nlte_corrections`, grid `Fe_Bergemann_MPIA.csv`) +
the stale `Si_Bergemann_MPIA.csv`. That leg is NOT covered by the RYA-413 registry
guard (`_load_mpia_element_grid`, for Ca/Ti/Cr) — which is exactly why these survived.
Fe is not in `NLTE_CORRECTION_ELEMENTS`; it has its own dedicated path.

## Part A — per-line LIVE MPIA verification (load-bearing)

Candidates (lines identically 0 across ALL nodes in the committed grid): **32** —
**14 Fe I + 17 Fe II + 1 Si I** (the ticket cited "23 Fe + Si 6253"; the live grid carries
**31 Fe**, reported as-found). Each re-queried live against the MPIA tool at 6 representative
nodes spanning the grid box (the same POST the grid was built with):

| class | count | meaning | action |
|-------|-------|---------|--------|
| **TRUE_PLACEHOLDER** | **32 (all)** | MPIA serves the line but returns literal 0.0 at every node | DROP |
| GENUINE_NEAR_ZERO | 0 | small-but-real correction (NLTE-insensitive) | (would KEEP) |
| NOT_SERVED | 0 | only error codes | (n/a) |

**Every one of the 32 is a TRUE placeholder** — MPIA returns `0.000` at every served node
(a handful of Fe II nodes returned MPIA error codes → NaN, but no served node was nonzero).
No genuine near-zero physics exists among them, so dropping deletes nothing real. Per-line
verdict table: `data/audit/fe_leg_placeholders_rya417.csv`
(`scripts/audit_fe_leg_placeholders_rya417.py --offline|live`).

## Part B — symmetric guard

`detect_placeholder_zero_lines` (RYA-413) now also guards the **Fe leg** loader
`_load_mpia_fe_grid` (ion-aware — the grid carries Fe I + Fe II) and the **curation**
target `curate_nonfe_pools.solar_nlte_delta`; both **refuse** a placeholder-carrying grid
loudly (`ValueError [PLACEHOLDER_ZERO]`), matching the registry leg. Build-time drop added
to `build_fe_nlte_grid_rya319.py` (mirrors `build_nlte_grids_mpia`). After this, no
placeholder-zero can silently register or apply on ANY leg. Tests:
`tests/test_fe_leg_placeholder_guard_rya417.py` (ion-aware detect, Fe-leg refuse,
genuine-near-zero kept, committed grid loads clean, curation-leg refuse).

## Part C — calibration impact (re-derived, NOT tuned)

Verified placeholders dropped (`scripts/drop_fe_leg_placeholders_rya417.py`); clean solar
Fe gate re-run:

| quantity | contaminated | clean | shift |
|----------|-------------|-------|-------|
| A(Fe I) NLTE **median** (gate verdict) | 7.516 | 7.516 | **≈0 (<0.001)** |
| A(Fe I) NLTE **mean (A+B)** | 7.5066 | 7.5082 | **+0.0016** |
| A(Fe II) NLTE | 7.7145 | 7.7145 | 0 |
| Fe I scatter (σ) | 0.138 | 0.138 | 0 |
| mean Δ(NLTE) diagnostic | +0.0089 | +0.0105 | +0.0016 |

The per-line abundances are **unchanged** (a placeholder line gets δ=0 whether it is
"NLTE-corrected to 0" or "no correction / 1D-LTE retained"). So the **median** Fe reference
(the gate's primary) does not move, and the **[X/Fe] ripple is ≈0** for the differential
backbone. The only shift is on the **mean** estimator (+0.0016 dex): the 6 low-side fake-zero
Fe I lines were dragging the mean of the NLTE-corrected set down; excluding them lifts it
+0.0016 — the Fe analogue of Ca's +0.0052, but smaller and median-robust. Reported as-derived;
never nudged toward 7.46.

## Part D — Fe I pool intersection × the RYA-407 floor

Of the solar Fe pool, **6 Fe I + 1 Fe II** lines fall within 0.15 Å of a placeholder and
received a fake-zero correction. Their residuals from the Fe I median (7.5095):

| line | A(Fe I) | residual | side |
|------|---------|----------|------|
| 5225.5 | 7.500 | −0.010 | mid |
| 5855.1 | 7.437 | −0.073 | low |
| 6219.3 | 7.449 | −0.061 | low |
| 6609.1 | 7.561 | +0.052 | high-ish |
| 6157.7 / 6430.8 | — | (not in the 62-line gate pool) | — |

**Verdict: the placeholder lines do NOT account for the RYA-407 floor.** The 407 signature is
an asymmetric **high** tail (15 lines > +0.10); the overlapping placeholder lines are
predominantly **low/mid-side**. And the fix is numerically inert on the scatter (σ = 0.138 →
0.138, identical), because a fake-zero and a real-zero produce the same per-line abundance —
there was no hidden correction to restore (MPIA has no nonzero value for these lines). So the
**~0.13 floor stands as bedrock** for 277/282/279/395; it is not partly a fixable placeholder
defect. The RYA-417 fix is about honesty (no fake-zero NLTE flag) + the recurrence guard, and
removes a tiny low-side contaminant from the mean estimator only.
