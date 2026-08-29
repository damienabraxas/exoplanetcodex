---
name: codex-product-audit
description: Standard per-product audit for any completed Exoplanet Codex abundance product. Use this skill whenever Ryan says "audit this product", "audit each product", "did we do the error right", "why is this value high/off", "sign off this product", or before any product is marked final, combined into a consolidated value, or put on the website. Covers the uncertainty budget (Type A sigma/sqrt(N) vs raw sigma, Type B stellar-parameter systematics, sourced delta_xi, sigma_syst-is-a-total-not-a-floor), label/provenance integrity (do route/gf/atmosphere/departure labels match what was actually computed), the gf graded/deepgraded differential (applied everywhere it should be), the gate-as-flag check, and the standard scientific-interrogation battery (why is the value where it is; continuum; blends; broadening; model domain in UV/IR; line selection; portability to other stars). Always use before signing off a product and whenever a value looks anomalous.
---

# Codex Product Audit Skill

## Purpose

Every completed abundance product gets audited before it is signed off, combined into a consolidated value, or published. Two jobs, always both:

1. **Is the uncertainty budget computed correctly and honestly?**
2. **Do we understand WHY the value sits where it does?**

Two permanent rules frame the whole audit:

- **Report honest error bars — never tune to a reference (RYA-161).** A large bar or an offset is a puzzle to *diagnose* (RCA), not a number to accept blindly nor to tune away.
- **The gate is a FLAG, not a blocker.** A product with an honest bar that exceeds the field target is REPORTED and LABELED, never hidden. Wabi-sabi: the honest bar is the truth of the measurement; a fake-precise tuned number is the ugliness.

Never change a value during an audit. If a defect is found, file a fix ticket.

---

## Part A — Uncertainty budget (did we do the error right?)

**Statistical (Type A):**
- [ ] Reported statistical error is **sigma/sqrt(N)**, not raw line-to-line sigma (raw sigma is QA/diagnostic only) — RYA-282.
- [ ] **Small-N flag:** if sigma/sqrt(N) is inflated by few lines (e.g. NIR KP Synth.1D-NLTE.Bergemann = 7.503 +/- 0.497, n=7), record it and route to Part E — is it irreducible (few clean lines) or a fixable defect?
- [ ] Per-line rejections logged (no silent line drops) — a suppressed outlier hides real scatter.

**Systematic (Type B):**
- [ ] Stellar-parameter systematics (Teff, logg, xi) propagated into sigma_reported.
- [ ] **delta_xi is SOURCED**, not the uncited 0.05 (RYA-1089). The honest delta_xi reflects the real method+selection spread, not a formal error on a xi value we don't actually adopt.
- [ ] Region-specific systematic included (UV / VIS / IR each differ) — see Part E.
- [ ] **sigma_syst is a real total, not a FLOOR** with unmeasured terms (e.g. the Amarsi MLP-network term left uncomputed, RYA-1095). A floor quoted as a total understates the bar.

**No tuning:**
- [ ] No value, threshold, or bar was chosen because it lands near a reference or passes a gate (RYA-161). Every constant carries a citation (single source of truth).

---

## Part B — Label / provenance integrity (do the labels tell the truth?)

The recurring defect class this session: a label asserting a property with no evidentiary basis. **Measure from the artifact, never trust the tag.**
- [ ] `route` (ew / synth) matches the handler that actually ran — not a pinned `_ROUTE_BY_LABEL` or a stranded hardcoded handler (RYA-1104: the "EW" Amarsi product was synth).
- [ ] `gf` tier matches the lines actually used — not a LEGACY constant (RYA-1104: `gf=kurucz` published on an all-lab product).
- [ ] `atmosphere` / `deck` labels match the model actually run.
- [ ] "departures applied" columns reflect reality (RYA-1104: `nlte_delta_dex=0.0` / "no departure applied" on a product that DID apply departures).
- [ ] The error budget does not charge a term the product didn't earn (RYA-1104: +0.0129 dex profile-fitter residual on a synthesis-measured product).

---

## Part C — gf-tier differential (applied everywhere it should be?)

- [ ] The graded / deepgraded gf treatment is applied to THIS product in **every region and holding** where it applies — not just Fe I VIS.
- [ ] Graded lab-gf lines and deepgraded / Kurucz-floored lines are handled per the ratified differential, and the tier is stamped **per line** (lambda+EP dual key, RYA-1037).
- [ ] Where a region is heavily gf-floored (e.g. near-UV ~84% Kurucz, RYA-822), that is REPORTED as a tier fact, not silently folded into the value.
- [ ] Cross-check: does the tier composition explain part of a value offset? (The Elo/gf line-population story — RYA-1104.)

---

## Part D — Gate as flag, not blocker

- [ ] The product is REPORTED with its honest bar regardless of whether it clears the field target (~0.05 — itself to be sourced, RYA-282).
- [ ] The gate renders as a **LABEL** (bar vs target), never suppresses the product.
- [ ] If the bar exceeds the target, that is stated plainly with the reason from Part E — not hidden, not tuned down.

---

## Part E — Scientific interrogation (why is the value where it is?)

Ask ALL of these, every product. A value that "looks fine" still gets asked; a value that looks off gets RCA'd until the cause is **named**.

**Measurement**
- Is the measurement correct? Why is the value where it is — especially if high/low vs the reference?
- Is there something we missed — a line set, a correction, a data product?

**Continuum**
- Was normalization done appropriately for this region — a broad global continuum, or region-local? What is the continuum doing in the **UV** specifically (crowded blue: is there even a true-continuum window)?
- Could continuum placement be biasing the value (the UV-high symptom on both Kitt Peak holdings)?

**Blends / line quality**
- Any blended, telluric-bitten, or bad lines inflating scatter or biasing the value? (Blends resolve by synthesis, not EW — RYA-287/338/581.)

**Broadening**
- If line broadening (xi, macroturbulence, rotation) is driving the measurement, can the model handle it, or is it an unmodeled systematic?

**Model domain (UV / IR)**
- Is the model IN-domain for this region? (Amarsi MLP is Fe I VIS only, encodes wavelength, ~1.82 eV training floor, never saw an IR photon; `atom.fe607a` tops out at 7.505 eV, below IR upper levels.)
- What uncertainty does out-of-domain extrapolation introduce in the UV and the IR — is it quantified and carried, or silently ignored?

**Line selection / portability**
- Is our line set appropriate and complete for this region? Does it miss the reference's lines (the Asplund low-Elo gap, RYA-1106)?
- Does this treatment / line set **travel to other stars**? Bring the reference (e.g. Asplund) line list along for the ride — how does it compare on Alpha Cen, 55 Cnc? What lines did the reference analyses choose for *those* stars?

**RCA verdict per finding:** IRREDUCIBLE (accept + report, maybe add lines later) vs FIXABLE DEFECT (fix, without tuning) vs OPEN (log to OPEN_QUESTIONS.md).

---

## Verdict

| Dimension | Status | Notes |
|---|---|---|
| Type A error is sigma/sqrt(N) | check | |
| Type B systematics + sourced delta_xi | | |
| sigma_syst is a total, not a floor | | |
| Labels match the measurement | | |
| gf graded/deepgraded differential applied | | |
| Gate rendered as a flag (not a block) | | |
| Value interrogated (Part E), cause named | | |
| No tuning / everything sourced | | |

**Verdict: SIGN OFF / SIGN OFF WITH CAVEATS (documented) / HOLD (defect or open question to resolve).**

---

## Output

Post the audit as a Linear comment on the product's audit ticket (or the per-element audit umbrella). Log any OPEN findings to `OPEN_QUESTIONS.md`. Never change a value during an audit — a defect gets its own fix ticket.
