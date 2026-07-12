# RYA-545 — Ti I gf-scale re-grade + blend-aware ionization-balance gate: STOP (measurement-limited), grid CORROBORATED

**Date:** 2026-07-11 · Branch `ryandamienschmitt/rya-544-...` (RYA-545 work). **NO MERGE. NOT wired.**
Scripts: `scripts/rya545_route1_blended.py` (PySME) + `scripts/rya545_route2_ts.py` (production TS).
Logs: `rya545_route1.log`, `rya545_route2_ts.log`. Runs on Sirius.

## The task

RYA-544 STOPPED because OUR Ti I pool had a +0.6-dex gf-scale problem (MAD 0.5) swamping the +0.05
NLTE correction. RYA-545: fix the Ti I gf-scale **by pre-declared gf-provenance grade** (never by
which line breaks balance), then re-run the ionization-balance gate. Ti II = untouched validator.

## Pre-declared gf-provenance retire criterion (fixed BEFORE inspecting any abundance)

KEEP a Ti I line iff its `canonical_gf.loggf_reference` is a PRIMARY LABORATORY transition
probability of NIST quality: **LGWSC / 2013ApJS** (Lawler et al. 2013 — the definitive Ti I lab gf),
`1982/1983/1986MNRAS` (Blackwell/Oxford lab), `NWL/MFW/BLNP/SK`, or `nist_grade` set. RETIRE
ungraded **K10** (GES-synth-fitted gf) + bare **VALD3** compilation. (`nist_grade` is blank for ~all
Ti — the provenance proxy is `loggf_reference`.) In the gate pool this retired 11 (10 K10 + 1 VALD3),
kept the lab lines. Firewall: by provenance, not by balance.

## First gate (single-line COG, RYA-544 tool) — gf-scale FIXED, but at the noise floor

Full pool: LTE balance **+0.590** → lab-gf pool: LTE balance **−0.003** (Ti I 5.006 = Ti II 5.008).
The gf-scale systematic is gone. NLTE balance +0.050 — at the gate edge, BUT Ti I MAD 0.508 / SEM
0.26 (single-line COG, Unsold vdW, a 6745/6746 blend artifact) ⇒ noise floor, not a verdict.

## Refined gate — TWO independent blend-aware instruments (pre-declared: widen-by-quality, blend-fix, correct-inversion)

Both on MARCS, reference-blind, lab-gf Ti I pool (all EW), Ti II untouched:

- **Route 1 — PySME full-window BLENDED synthesis** (`linelist_full`, proper ABO vdW, differential EW
  = with-target minus without-target, isolates the line on top of the blend):
  Ti I **5.013**, Ti II **5.004**, **LTE balance +0.009**, NLTE +0.061 ± 0.207 (SEM 0.15).
- **Route 2 — production Turbospectrum synth-EW bisection** (RYA-285, iSpec `generate_spectrum` +
  GES linelist; iSpec TS synthesizer compiled here, `make turbospectrum`). The `_bisect` ew_floor
  method mis-converges saturated lines to A=7–8 (returns conv=True on runaway — unlike Route 1's
  rail), so a physical-plausibility bracket A∈[4.07,5.87] (= Route 1's bracket, applied identically)
  drops them: Ti I **5.283**, Ti II **5.301**, **LTE balance −0.018**, NLTE +0.035 ± 0.157 (SEM 0.16).

## Result — the grid is CORROBORATED; the balance is measurement-limited

Both independent instruments agree:
- **LTE ionization balance ≈ 0** (+0.009 / −0.018) — the lab-gf re-grade + blend-aware synthesis
  fixed the systematic (from +0.59). Ti I and Ti II agree in LTE.
- **NLTE balance ≈ +0.05** (+0.061 / +0.035) — consistent with the Mallinson-2024 +0.0506 correction
  tipping Ti I just above Ti II.
- Absolute zero-points differ by ~0.28 (different linelist/gf/synthesizer) but the balance is
  differential and robust.

**Per the pre-declared gate — require SEM ≪ 0.05 THEN |balance| < 0.05 — the verdict is STOP:** SEM
(~0.12–0.16) is NOT ≪ 0.05, so the balance (point estimate ~+0.05, consistent with the grid) cannot
be **certified** to the 0.05 gate. Even with lab gf + blend-aware synthesis and two instruments, the
per-line scatter (~0.2 dex intrinsic: EW-measurement error + synthesis fidelity) over thin reliable
pools (6–7 Ti I, 3 Ti II; saturated lines rail/mis-converge) floors SEM at ~0.12. **Ti I's ionization
balance on OUR data is measurement/precision-limited — a real finding, earned not defaulted.**

## Disposition — STOP, no wire (grid corroborated but not certified)

- **Do NOT wire.** `constants.py` Ti stays `Ti_Bergemann2011_MPIA.csv`; register Ti stays CHECK; the
  RYA-534 Ti strict-xfail is NOT flipped. (Per the pre-declared STOP branch.)
- **What is now known (strengthened):** (1) the correct Ti I NLTE correction is ~+0.05 ab-initio
  (RYA-544, Mallinson-2024, derived not tuned); (2) with proper lab gf + blend-aware synthesis the
  LTE Ti I/Ti II balance is ≈ 0 and the NLTE balance ≈ +0.05 — the Mallinson grid is **corroborated**
  by two independent instruments; (3) the remaining barrier is measurement precision (SEM ≫ the 0.05
  bar), not the grid or the gf-scale.
- The grid remains banked (Zenodo 10753497, M.2, md5-pinned) + derived (+0.0506). Wiring it needs a
  precision the ionization balance can't deliver on our current solar Ti EWs — a better solar Ti EW
  set / more lab lines / a careful line-by-line analysis would be required to certify, OR the wire
  decision rests on the corroboration rather than a certified <0.05 balance (a science call).

## Engineering notes (reusable)
- Compiled iSpec's Turbospectrum synthesizer on Sirius: `cd /srv/codex/engines/ispec_src && make
  turbospectrum` → `synthesizer/turbospectrum/bin/{babsma_lu,bsyn_lu}`; then
  `ispec.is_turbospectrum_support_enabled()` → True. Env: venv312 + `ISPEC_DIR=PYTHONPATH=
  /srv/codex/engines/ispec_src`; `mkdir -p /tmp/ispec_codex_synth`.
- The production MOOG-EW path (`_ew_to_abundance`, EW_BASELINE_CODE='moog') is NOT usable standalone
  for an arbitrary Ti pool: MOOG `abfind` is single-line (not blend-aware; A(Ti II)=8.5 on a blend),
  and the iSpec 42000_VALD region file overlaps only 6/81 of our Ti EWs (0/37 lab). Blend-aware
  requires SYNTHESIS (Route 1/2).
- The pipeline synth-EW (`_bisect_synth_abundance`) needs a NARROW per-line window and a physical
  A bracket — its ew_floor method runs A away on saturated/wide-window lines (returns conv=True).
