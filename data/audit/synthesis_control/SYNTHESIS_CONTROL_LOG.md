# Synthesis handler — optical control log — RYA-713

Ryan: *"Lets look at it from both angles"* → *"run it against harps and see if we hit 7.520"*
→ *"make it a thin adapter over `_fit_synth_flux`"*.

## Every run, same question: does the handler reproduce a known answer?

| # | change | spectrum | ref | n | A(Fe I) | offset | scatter |
|---|---|---|---|---|---|---|---|
| 1 | as built, `macro = 0` | KP | 7.466 | 8 | 7.046 | −0.420 | — |
| 2 | vmac 3.8 / vsini 1.8 from resolver | KP | 7.466 | 7 | 7.226 | −0.240 | 0.232 |
| 3 | response-weighted χ² | KP | 7.466 | 7 | 7.226 | −0.240 | 0.280 |
| 4 | target-only EW, `xi` sourced | KP | 7.466 | 9 | 7.266 | −0.200 | 0.278 |
| 5 | envelope in every band | KP | 7.466 | 8 | 7.346 | −0.120 | 0.285 |
| 6 | coverage quarantine, band clipped | KP | 7.466 | 10 | 7.266 | −0.200 | 0.481 |
| 7 | synth-v2 window + constant σ | KP | 7.466 | ? | 7.546 | **+0.026 vs 7.520** | 0.571 |
| 8 | → HARPS, ref 7.520 | HARPS | 7.520 | 9 | 7.956 | +0.436 | 0.446 |
| 9 | envelope re-gated to pseudo bands | HARPS | 7.520 | ? | **7.516** | **−0.004** | 0.428 |
| 10 | **thin adapter over `_fit_synth_flux`** | HARPS | 7.520 | 11 | 7.409 | −0.111 | **1.181** |
| 11 | + χ² acceptance gate | HARPS | 7.520 | 4 | 7.428 | −0.093 | 0.492 |

**Banked comparand: synth-v2 gives A(Fe I) = 7.520, A_X_std = 0.066, on 23 HARPS lines.**

## Six defects found, all in my reimplementation

1. **`macroturbulence = 0`** — inherited from EW-mode defaults, where it is correct because
   EW is broadening-invariant. A flux fit is not. Bypassed `_resolve_broadening` (RYA-288),
   whose docstring says exactly this. Worth −0.42 dex.
2. **Fixed ±1.00 Å window** — synth-v2 uses wing-wide ±0.237 Å at 45 mÅ. 4× too wide dilutes
   the line among blends and continuum.
3. **Response-weighted χ²** — breaks the σ-invariance synth-v2 relies on by design. Did not
   help either (`corr(EW, ΔA)` +0.891 → +0.878).
4. **Envelope in every band** — harmless at 0.400 Å knot spacing, self-cancelling at
   0.095 Å because the *synthetic's* envelope deepens with trial abundance and divides out
   the sensitivity being fitted. Run 8: +0.436.
5. **Circular seed** — `a_start` defaulted to the reference. Three of eight lines returned
   7.466 exactly (`ΔA = −3.6e-15`) and counted as perfect.
6. **No χ² acceptance gate** — synth-v2: *"the gate is FIT QUALITY… acceptance is on
   merit."* Without it, badly-fitted lines entered the median: scatter 1.181 vs 0.066.

## The honest conclusion

**The oscillation is the finding.** Runs 7, 9 and 11 land at +0.026, −0.004 and −0.093 with
scatters of 0.571, 0.428 and 0.492 — a median that moves 0.12 dex between configurations
while the spread stays 7× the banked 0.066. **A median landing near the answer with that
scatter is not a controlled method; it is a small sample averaging out.**

Delegating the fit did not fix it, which localises the remaining problem to **the control
harness, not the fitter**: line selection, the observed-spectrum handling, or the context.
synth-v2 fits `last_linemasks` — a curated set carried from the EW stage with per-line
windows from measured EWs — while this control hand-picks 4–16 pool lines by EW range.

## Recommendation

Stop iterating this harness. Two cleaner options:

1. **Run synth-v2's own driver** (`_run_synthesis_v2_mode`) on the current tree and confirm
   it still returns 7.520. That separates "did something regress" from "is my harness
   wrong", which 11 runs have not separated.
2. **Have the handler wrap that driver** rather than re-deriving the per-line loop. The
   handler's legitimate job is band policy, coverage checking and quarantine translation —
   the line selection and aggregation belong to the validated path too.

## What survives regardless

* **Angle 1 ≈ 0.54–0.72** across every configuration: synthesis attributes ~55–70 % of the
  pool's EW to Fe itself. The profile fitter assigns blended flux to the line; synthesis
  separates it. **That is the argument for synthesis in crowded regions**, and it is a
  result, not a bug.
* **The coverage quarantine.** 4065.381 Å sits below the GES list's 4200 Å blue limit and
  synthesises to EW 0.00 against a pool value of 73.59. Real, and it is the same 4200 Å
  wall that blocks near-UV synthesis and Engine B entirely.
