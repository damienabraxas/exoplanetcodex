# RYA-669 — RYA-527 Phase 2 run summary

**REVIEW ARTIFACT. Nothing frozen, promoted or adopted.** Gold v3 is byte-untouched; `data/reference/solar/` was not written to.

## The answer

🛑 **STOPPED — 1 §4 STOP condition(s) tripped.** The artifacts below were still written so the evidence is readable, but no further step was taken.

- Fe 1D→3D correction WILL DOUBLE-APPLY on the next phase_c regeneration: gold CURRENT holds A_X=7.466 (the POST-correction 3D value) but labels it method_scale='1D-NLTE (Fe I)', which carries no '3D'. The idempotency guard in phase_c_verdict_rya371.py keys on that label, so it re-applies 0.0 and the anchor lands at 7.466 — INSIDE FE_GATE [7.410, 7.510], so no gate catches it. See BLOCKING_FINDING_fe_double_correction.md.

## What is actually fresh here

| leg | fresh? | what it means |
|---|---|---|
| two-engine record (`solar_two_engine_records.json`) | **YES — re-computed** | both engines re-driven over real solar data on current main (Sirius): Engine A EW→A(X) per line + production NLTE δ; Engine B a new Turbospectrum synthesis-v2 flux fit per line + Gerber TS-native NLTE δ. Nothing cherry-picked from the July 18 branch. |
| phase_c verdict | fresh run, **derived** input | since RYA-469 phase_c CLASSIFIES the frozen gold (`read_solar_reference('CURRENT')`) plus the dedicated-channel measurements — it does not re-derive A(X) from spectra. Re-running it answers "does the freeze re-classify consistently", NOT "does the measurement reproduce". |
| disposition report | **YES** | same classifier, run over the FRESH record. It was expected to retire the gate-3 staleness RYA-663 flagged; it does not, and the section below shows why that flag cannot clear by re-running. |

Gold compared against: **v3**. Verdict counts: `{'PASS': 7, 'CURATION-OWED': 19}`.

## Elements whose value moved > 0.01 dex vs frozen v3

| element | v3 | v4 | Δ | why |
|---|---|---|---|---|
| Si | 7.888 | 7.64 | -0.248 | two-engine synthesis floor (I, B_synth) |
| Cr | 6.022 | 5.661 | -0.361 | two-engine synthesis floor (I, B_synth) |

## Gate 3 — the two signals, separately (RYA-675)

✅ **Age: fresh.** The two-engine artifact is at least as new as every live input it is read against, so nothing here is provisional on age. This is the state a fresh re-emit produces, and the state the pre-RYA-675 detector could not report because any value difference read as age.

**Cross-channel disagreement: 5 element(s).** Remedy **ADJUDICATE** — these do NOT regenerate away; a two-engine floor's diagnostic leg is designed to be able to disagree with the ratified leg, and someone must pick which channel is the reported truth.

| element | two-engine | live phase_c | delta | why they differ |
|---|---|---|---|---|
| Cr | 5.661 | 6.022 | -0.361 | different legs — Cr I synthesis floor vs the gf_floor EW value |
| Fe | 7.584 | 7.466 | +0.118 | the RATIFIED Fe policy — the two-engine number is a diagnostic that sits above the anchor BY CONSTRUCTION (Ryan, 2026-07-16) |
| O | 8.730 | 8.735 | -0.005 | rounding |
| S | 7.369 | 7.486 | -0.117 | different legs — EW vs the RYA-492 Costa-Silva dedicated synthesis |
| Si | 7.640 | 7.888 | -0.248 | different legs — synthesis floor vs the gf_floor EW value |

Promotes under the ratified three gates: **Ca**

## The four species-adoption decisions — NOT adopted, Ryan decides

**D.1 — Sr: adopt Sr II 2.759 into gold v4?**

- Fresh number: `2.759`
- Recommendation: ADOPT Sr II 2.759 — reproduces within ±0.02 of the RYA-643 corrected 2.759.
- ⚠ near-UV fit red-chi2 78.27 — high, and flagged as the open near-UV item in RYA-643 (red-chi2 78-180 across the blue channels). The line is marked reliable on dEW/dA grounds, not on chi2. Adoption inherits that systematic; it is a known ~0.05-0.1 dex near-UV floor, not noise.

**D.2 — Co: which Co value goes into gold v4?**

- Fresh number: `4.96`
- Recommendation: Phase 2 cannot break this tie — it produced no Co number (see above). On the two existing candidates the split is 0.005 dex, inside every gate and below the reported precision, so it is a provenance choice, not a measurement one: 4.960 comes from the run with the RYA-643 rest-frame/gsig defect fixed, which is the better-founded of the two. Ryan picks.
- ⚠ Co is verdict PASS at tier `owed` ⇒ v3 freezes NO value for it (RYA-665). Whichever number is picked, it stays HELD until the tier moves — adopting a value here does not by itself freeze one.

**D.3 — Ba: fire RYA-581 deblend BEFORE the v4 freeze, or freeze 2.410 HELD-with-caveat?**

- Fresh number: `2.237`
- Recommendation: No recommendation on timing — this is a scheduling call, not a measurement one. The measurement fact: 2.410 is blend-inflated by ~+0.15 against a clean cross-check at 2.187/2.231, and RYA-581 exists to replace it. Freezing 2.410 into v4 would immortalise a number already known to be high; deferring costs one ticket.
- ⚠ Ba is `owed` tier in v3 ⇒ blank A_X, so nothing is frozen today either way. The urgency is about v4, not v3.

**D.4-Ca — Ca: promote Ca to PASS on the fresh gate 3?**

- Fresh number: `6.324`
- Recommendation: PROMOTE Ca — clears all three ratified gates on the FRESH cross-engine delta.

**D.4-Na — Na: promote Na to PASS on the fresh gate 3?**

- Fresh number: `6.264`
- Recommendation: DO NOT promote Na — gate 3: |dCE| = 0.124 > 0.1. This is a FRESH answer: RYA-664 cleared its gate 1, and gate 3 is now decided on a current delta rather than deferred.

## Known defect carried forward (not fixed here)

`data/reference/solar/solar_abundances_v3.csv` holds a **Sr I** row and no Sr II row, while the NLTE registry ratifies Sr as **Sr II** (RYA-551/643). The diff table therefore shows Sr II against a blank v3 cell. This is the RYA-663 defect, unchanged by the v3 freeze; adopting Sr II (D.1) is what would repair it.

