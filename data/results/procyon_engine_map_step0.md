# RYA-349 — Procyon synth-vs-EW engine map: STEP 0 RECON (hard gate)

Analysis-only. Harness: `scripts/procyon_element_map.py` (models `scripts/procyon_fe_2x2.py`).
Regenerate: `python scripts/procyon_element_map.py`
→ `procyon_step0_linelist.csv`, `procyon_step0_nlte_status.csv`.
STOP-at-verdict; no `STAR_PARAMS`/line-list edits; no merge.

Step 0 is a **hard recon gate**, confirmed LIVE here (never assumed). It surfaced three
findings that supersede the ticket's stale premises and reshape the run.

## Finding 1 — NLTE wiring is WIDER than the ticket assumed (Mg/Si/Na/Al now wired)
The ticket says "we have Fe + Ca/Ti/Cr; confirm Mg, Si, Na, Al." Live, the
`config.constants.NLTE_CORRECTION_ELEMENTS` registry now wires **Ca, Ti, Cr, Na, Mg, Si,
Al, Ba, Sr, Mn, S, K** (Amarsi 2020 PySME per-line grids via RYA-410), plus Fe on its
separate MPIA path. All grid CSVs are on disk. So **Mg/Si/Na/Al are wired**, not missing.

## Finding 2 — …but at Procyon's 6554 K, EVERY non-Fe grid is OUT OF BOUNDS
This is the gate. Grid Teff ceilings: Ca/Ti/Cr/Mn/Ba 6500 K; Na/Mg/Si/Al/S/K 6200 K;
Sr 6000 K. **Procyon Teff = 6554 K sits above all of them.** Only the RYA-319-extended
Fe MPIA grid (Teff→6600+) covers Procyon.

| element | wired | on disk | in-bounds @ 6554 K | status |
|---|---|---|---|---|
| **Fe** | yes | yes | **YES** | NLTE-LIVE |
| Mg, Si, Na, Al | yes | yes | **NO** (ceil 6200) | grid present, out-of-bounds → edge-clamp or LTE |
| Ca, Ti, Cr, Mn, Ba | yes | yes | **NO** (ceil 6500) | grid present, out-of-bounds → edge-clamp or LTE |
| Sr | yes | yes | **NO** (ceil 6000) | grid present, out-of-bounds → edge-clamp or LTE |
| Co, Cu, Sc, V, Li, Eu, Y | no | no | — | not wired → LTE (RYA-242-style grid need) |

**Consequence for the run:** the NLTE arm the ticket anticipated as "a real stress-test"
**cannot fire at Procyon** — the grids don't reach it. For every non-Fe element the {NLTE}
cell is either (a) 1D-LTE (flagged NLTE_unavailable) or (b) a 0th-order Teff-edge-clamp to
6500/6200 K (the Amarsi/Fe precedent, RYA-319/339) — a *finding to log*, not physics. So
**the Procyon "NLTE-critical?" column of the RYA-277 map is decided by coverage, not by
physics: only Fe has live NLTE at 6554 K.** Extending the non-Fe grids past 6550 K is a
RYA-319-style grid-extension need (one ticket per element / a batched warm-edge extension).

## Finding 3 — the RYA-347 gf-source confound is CLOSED
The ticket flags that synth (GES v6) and EW (`linelist_solar` VALD3) carried diverging gf
(|Δ| up to 0.125 dex), so a raw synth−EW gap is partly a gf artifact. Live: **both paths
now import `gf_resolver` and resolve against `canonical_gf.csv` (RYA-353)** — EW path
(`abundances_derive`) and synth path (`cno_synthesis`). Canonical table present. So the
confound is **structurally closed**: a residual synth−EW gap on Procyon is no longer a
gf-source artifact and can be read as engine/physics.

## Step 0.1 — in-list optical lines per species (HARPS 3780–6910 Å)
Full `linelist_procyon.csv` extraction counts (the measured/vetted EW pool is a subset,
staged downstream). Highlights: Mg I 62, Si I 151, Ca I 106, Ti I 536 / Ti II 1166,
Cr I 600, Na I 156, Al I 21; HFS: Mn I 1665, Co I 2598, Cu I 99, Sc II 494, V I 1791;
problem children: Ba II 92, Eu II 228, Y II 132; Li I not present in the extraction window
(the 6707 line is a single injected diagnostic — confirm at EW stage). Fe I 2837 / Fe II
207. Full table: `procyon_step0_linelist.csv`.

## Partial engine map (what Step 0 already fixes for RYA-277, Procyon row)
| axis | Procyon (6554 K) verdict from recon |
|---|---|
| NLTE-critical & available | **Fe only** — every other wired grid is out-of-bounds at 6554 K |
| gf confound | closed (single-source canonical_gf) → synth−EW gap is engine, not gf |
| not-wired (LTE-only here) | Co, Cu, Sc, V, Li, Eu, Y |
The per-species **engine recommendation** (EW vs synth) and **[X/H] vs literature** need
the four-cell matrix run (below).

## Steps 1–4 — GATED (exact resume point)
The four-cell {EW,synth}×{LTE,NLTE} matrix, per-line χ²ᵣ, Arm-B rescue tests, and the
final engine map run through `ad.run('procyon', engine='spectrum')` and `engine='synthesis'`
(both emit per-element A_X + A_X_nlte; NLTE edge-clamp logged per Finding 2). They require
the Procyon **processed** inputs, which are gitignored and not staged in this worktree:
`data/processed/procyon_normalized.csv`, `data/processed/procyon_ew.csv`.

**Resume:** stage them from the 20 raw HARPS frames (external store), then re-run the harness:
1. `pipeline.spectra_normalize` (procyon) → `procyon_normalized.csv`
2. `python -m pipeline.lines_fit --star procyon` → `procyon_ew.csv`
3. `python scripts/procyon_element_map.py` → runs `run_matrix()` (guarded snapshot/restore of
   the production processed files) → `procyon_engine_map_{ew,synth}.csv` + the engine map.

No abundance numbers are reported here — none were computed; Step 0 is the deliverable that
gates them, and it changes the run: at Procyon the NLTE stress-test is Fe-only until the
warm-edge grids land.
