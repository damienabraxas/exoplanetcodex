# RYA-407 — Solar Fe I scatter (σ = 0.138 dex) diagnostic → verdict (c) honest floor

Diagnostic-only (RYA-405 split). Owns **only** the solar Fe I line-to-line scatter; the
Fe II EW-vs-synth arbiter mismatch is RYA-406. On `main` 166cb4b — the correctly-integrated
state per RYA-405, so this is the **live** pool, not a stale snapshot. Fe I is sourced from
the committed GES reference (RYA-330), so this verdict is independent of RYA-406/408.

Reproduce: `python scripts/validate_fe_rya238.py --star solar` (writes `solar_per_line.csv`),
then `python scripts/diag_fe1_scatter_rya407.py`. The ranked table is committed at
`data/audit/fe1_scatter/fe1_per_line_residuals_rya407.csv` (byte-reproducible).

## Step 0 — reproduce + isolate
Gate (verbatim): `Fe I scatter = 0.138 -> FAIL (< 0.1 dex)`, n = 62 Fe I lines,
A(Fe I) NLTE 7.516, median A 7.5095. Diagnostic σ(1D) = 0.1398 (NLTE σ ≈ 1D; Fe I δ_NLTE
is ~uniform so it does not change line-to-line scatter).

## Step 1 — classify the drivers (residual-INDEPENDENT, reusing existing machinery)
Each mechanism is evaluated on its own cited signal, **not** on abundance residual, so the
test cannot cull toward a target σ.

| Mechanism | Signal (single-source) | Result |
|---|---|---|
| Path / COG / excitation systematic | corr(resid, REW/EP/log_gf/snr) | all \|corr\| ≤ 0.13 → **no trend** |
| gf error | per-line `log_gf` vs canonical `gf_synth_ges` (RYA-203/350) | already == GES (Δ ≈ 0) → **ruled out** |
| Saturation / COG | RYA-220 `saturation_score` ≥ 0.8 | drops 0 lines, σ unchanged → **not the driver** |
| Weak-line noise | `ew_snr_score` < 0.5 | σ 0.1398 → 0.1405 → **not the driver** |
| Proximity blend | single-source `vald_proximity_flag` ≥ 0.5 | σ 0.1398 → 0.1276 (drops 10) → **small, real (a)** |
| proximity ≥0.5 **and** snr<0.5 | both | σ → 0.1193 (n=45) — floor, not ~0.05–0.10 |

**Skew:** 15 lines > +0.10 vs 5 lines < −0.10, mean resid **+0.041** — a one-sided **high**
tail, i.e. unrecognized weak blends / 1D-LTE line-formation inadequacy, not symmetric noise.
Top outliers (e.g. 5481.243 +0.494 prox 1.0; 5104.030 +0.369; 5432.948 +0.263 prox 0.885)
are predominantly high-side and proximity-flagged or strong, consistent with the tail.

## Step 2 — verdict + routing
- **(b) path artifact — RULED OUT:** no resid-vs-REW/EP/gf/snr trend; NLTE σ ≈ 1D σ.
- **(a) real pool issue — PARTIAL:** the only cited mechanism that moves σ is the
  proximity-blend cull (0.138 → 0.128; both-cuts → 0.119). Real but small; does not explain
  the gap to the historical ~0.05–0.10.
- **(c) honest floor — PRIMARY:** after all residual-independent cited vetting σ floors at
  **~0.12–0.13**, well above the ~0.10 gate. The current solar Fe I floor on the GES-EW
  1D-LTE pool is genuinely ≈ 0.12–0.13; **no cited mechanism reaches 0.10**. Per the
  ticket's CRITICAL rule, this is reported as the floor — **not** massaged toward a smaller
  number.

**Routing:**
- (c) honest floor → **RYA-277** (per-spectral-type gate threshold — 0.10 is mis-set for this
  pool) + **RYA-282** (σ/uncertainty budget that adopts the measured floor).
- small (a) → **RYA-395** (principled proximity-blend cull of the named lines) +
  **RYA-279** (compute the gate σ over the ceiling'd/vetted pool, not the raw output).

No pool was culled here (analysis-only): the proximity tail is a small, cited effect routed
to RYA-395, not applied, because it does not by itself clear the gate and the dominant result
is the honest floor.
