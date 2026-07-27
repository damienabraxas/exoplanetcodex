# RYA-527 — re-emitted solar verdict + PROPOSED gold v3 (two-engine floor)

_REVIEW ARTIFACT — NOT frozen. Freeze = `promote_solar_reference.py --apply` (Ryan). Fe = 7.466 (3D, RYA-553); two-engine 7.580 is the RYA-525 diagnostic only._

Verdict counts: {'PASS': 6, 'NLTE-OWED': 1, 'CURATION-OWED': 19}

## Proposed v3 vs v2 vs Asplund (tiered), with the two-engine per-element record

| El | Asp | v2 gold | v3 proposed | dAsp | verdict | engA | engB | selected | dCE | source |
|----|-----|---------|-------------|------|---------|------|------|----------|-----|--------|
| O I | 8.69 | 8.735 | 8.735 | +0.045 | PASS | - | 8.73 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| C I | 8.46 | 8.491 | 8.491 | +0.031 | PASS | - | 8.491 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Fe I | 7.46 | 7.516 | 7.466 **NEW** | +0.006 | PASS | 7.475 | 7.495 | A_1dnlte,B_synth | 0.084 | EW Fe I ionization-gated, 3D-corrected (RYA-406/407/553); two-engine 7.580 is the RYA-525 cross-engine diagnostic ONLY |
| Ca I | 6.3 | - | 6.372 **NEW** | +0.072 | PASS | 6.344 | 6.327 | A_1dnlte,B_synth | 0.016 | two-engine synthesis floor (I, A_1dnlte,B_synth) — RYA-561 PASS: 534-validated atom [Ca_gerber2023.prov.json], d_ref +0.072 <= 0.1, dCE +0.016 <= 0.1 over 4 lines; 1D-NLTE value vs 3D-NLTE reference; un-applied 3D term folded into the offset |
| Mn I | 5.42 | 5.47 | 5.466 **NEW** | +0.046 | PASS | - | 5.466 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| K I | 5.07 | 5.099 | 5.099 | +0.029 | PASS | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| N I | 7.83 | - | 8.202 **NEW** | +0.372 | NLTE-OWED | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Mg I | 7.55 | - | 7.614 **NEW** | +0.064 | CURATION-OWED | - | 7.614 | B_synth | - | two-engine synthesis floor (I, B_synth) |
| Si I | 7.51 | 7.888 | 7.639 **NEW** | +0.129 | CURATION-OWED | 8.06 | 7.633 | B_synth | -0.756 | two-engine synthesis floor (I, B_synth) |
| S I | 7.12 | - | 7.369 **NEW** | +0.249 | CURATION-OWED | 7.369 | 7.41 | A_1dnlte | 0.041 | two-engine synthesis floor (I, A_1dnlte) |
| Al I | 6.43 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |
| Na I | 6.24 | - | 6.369 **NEW** | +0.129 | CURATION-OWED | 6.37 | 6.248 | A_1dnlte | -0.121 | two-engine synthesis floor (I, A_1dnlte) |
| Ni I | 6.2 | - | 6.253 **NEW** | +0.053 | CURATION-OWED | 7.297 | 6.2 | B_synth | -1.253 | two-engine synthesis floor (I, B_synth) |
| Cr I | 5.62 | 6.022 | 6.022 | +0.402 | CURATION-OWED | 6.036 | 5.638 | B_synth | -1.058 | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| P I | 5.41 | - | 6.61 **NEW** | +1.200 | CURATION-OWED | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Ti I | 4.97 | - | 4.88 **NEW** | -0.090 | CURATION-OWED | 6.401 | 4.862 | B_synth | -1.326 | two-engine synthesis floor (I, B_synth) |
| Co I | 4.94 | - | 6.128 **NEW** | +1.188 | CURATION-OWED | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Cu I | 4.18 | - | 4.345 **NEW** | +0.165 | CURATION-OWED | - | 4.345 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| V I | 3.9 | - | 3.917 **NEW** | +0.017 | CURATION-OWED | - | 3.917 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Sc II | 3.14 | 3.203 | 3.203 | +0.063 | CURATION-OWED | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Sr II | 2.83 | - | 2.769 **NEW** | -0.061 | CURATION-OWED | - | 2.769 | B_synth | - | two-engine synthesis floor (II, B_synth) |
| Zr I | 2.59 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |
| Ba I | 2.27 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |
| Y I | 2.21 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |
| Li I | 1.05 | 0.727 | 0.727 | -0.323 | CURATION-OWED | 0.727 | 1.409 | B_synth | 0.682 | upper-limit disposition (phase_c, RYA-103/458); two-engine synth recorded as DIAGNOSTIC-ONLY, never the reported value |
| Eu I | 0.52 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |

## RYA-561 floor-promotion ledger (ratified 2026-07-27 — STRICT gate 3)

Gates: (1) RYA-534 anchor-validated NLTE atom · (2) |dAsp| <= 0.1 · (3) a REAL cross-engine |dCE| <= 0.1. A MISSING dCE **fails** gate 3 — the atom delta may never substitute for it (that is gate 1 under a second name; validate-don't-tune, RYA-161).

| El | A(X) | dAsp | n_lines | g1 atom | g2 tol | g3 dCE | dCE | outcome |
|----|------|------|---------|---------|--------|--------|-----|---------|
| Ca I | 6.372 | +0.072 | 4 | PASS | PASS | PASS | +0.016 | **PROMOTED -> PASS** |
| Mg I | 7.614 | +0.064 | 2 | PASS | PASS | FAIL | **None** | held (CURATION-OWED) |
| Si I | 7.639 | +0.129 | 15 | PASS | FAIL | FAIL | -0.756 | held (CURATION-OWED) |
| S I | 7.369 | +0.249 | 1 | FAIL | FAIL | PASS | +0.041 | held (CURATION-OWED) |
| Na I | 6.369 | +0.129 | 2 | FAIL | FAIL | FAIL | -0.121 | held (CURATION-OWED) |
| Ni I | 6.253 | +0.053 | 3 | PASS | PASS | FAIL | -1.253 | held (CURATION-OWED) |
| Ti I | 4.88 | -0.090 | 6 | FAIL | PASS | FAIL | -1.326 | held (CURATION-OWED) |
| Sr II | 2.769 | -0.061 | 1 | PASS | PASS | FAIL | **None** | held (CURATION-OWED) |

Reasons:

- **Ca**: PROMOTED: 534-validated atom; |d_ref|=0.072 <= 0.1; |dCE|=0.016 <= 0.1. 1D-NLTE value vs 3D-NLTE reference; un-applied 3D term folded into the offset.
- **Mg**: HELD: gate3 NO cross-engine delta (single-engine record — zero independent confirmation of the value; atom-delta fallback is REJECTED, RYA-561)
- **Si**: HELD: gate2 |d_ref|=0.129 > 0.1; gate3 |dCE|=0.756 > 0.1
- **S**: HELD: gate1 NLTE atom not 534-validated (no RYA-534 Engine-B grid provenance on record for S); gate2 |d_ref|=0.249 > 0.1
- **Na**: HELD: gate1 NLTE atom not 534-validated (Na_gerber2023.prov.json: no verdict recorded); gate2 |d_ref|=0.129 > 0.1; gate3 |dCE|=0.121 > 0.1
- **Ni**: HELD: gate3 |dCE|=1.253 > 0.1
- **Ti**: HELD: gate1 NLTE atom not 534-validated (Ti_gerber2023.prov.json: CHECK (not PASS): |median-anchor|=0.114 > tol 0.06. FINDING — NOT a deck/line failure: departures engaged, 3 directly-MPIA-comparable well-measured lines all ~2x the Bergemann-2011 value => a genuine Ti model-atom difference (Gerber-2023 atom.ti503b/503-level vs Bergemann-2011). Analogous to the Mn ~1/2 finding; a cross-engine model-atom systematic for RYA-525, not a pass. First-pass weak lines (5702/5703, 0.9-5.7 mA) gave +0.228 vs a wrong +0.05 guess.); gate3 |dCE|=1.326 > 0.1
- **Sr**: HELD: gate3 NO cross-engine delta (single-engine record — zero independent confirmation of the value; atom-delta fallback is REJECTED, RYA-561)

## RYA-524 reconciliation (S / Sr / N)

- **S**: STALE (owed, no value / pre-492 gf) -> CURATION-OWED, two-engine 7.369. RYA-492 Costa-Silva gf; gf barely moves it (+0.004 on GES) so the +0.37 vs Asplund is a gf-scale floor (RYA-161), not a line-ID error — stays owed.
- **Sr**: WRONG-SPECIES (Sr I raw-EW +2.13 artifact) -> CURATION-OWED, Sr II 2.769 (two-engine synth, RYA-551). Discarded the Sr I raw-EW artifact; Sr II 4077 synthesis (INSPECT NLTE). ~0.05-0.1 dex near-UV systematic -> refinement owed, not a clean PASS.
- **N**: STALE/unwired (NLTE-OWED, no value) -> 8.202 / NLTE-OWED (Kittpeak); RYA-526 grid registered. N I NLTE grid registered (RYA-369/526). See the N WIRING FLAG: the KP channel does not auto-apply the grid delta in phase_c yet; the +0.36 residual is a data-channel/gf item, not an NLTE debt.

## Honest flags

- N: base phase_c = 8.202 / NLTE-OWED (Kittpeak red multiplets). RYA-526 registered the N I NLTE grid; the phase_c KP channel does not auto-apply its (~-0.014) delta, so N still reads NLTE-OWED here — a properly-wired verdict moves it to ~8.188 CURATION-OWED (data-channel/gf floor, not an NLTE debt). WIRING FLAG.
- Mg: HELD at CURATION-OWED-with-value 7.614 (+0.064 vs Asplund, INSIDE the 0.1 band, and the NLTE atom IS RYA-534 anchor-validated) — held by gate 3 ONLY: the record has NO cross-engine delta because Engine-A is suppressed for Mg BY DESIGN (RYA-520 synth-required + the SAT-culled b-triplet EW pool), so the value rests on 2 Engine-B lines with zero independent confirmation. This is a RATIFIED HOLD (Ryan 2026-07-27), not an omission: the atom delta may NOT stand in for gate 3 (it is gate 1 under a second name — validate-don't-tune, RYA-161). Mg's path to PASS is the real second-line measurement, Mg I 5528 (RYA-592).
- Fe: reported 7.466 (3D, RYA-553). Two-engine 7.580 is the RYA-525 cross-engine diagnostic ONLY (per-line winner-combine biases high) — NOT the verdict (Ryan 2026-07-16).
- S: the committed two-engine record engineA=7.369 is on the PRE-RYA-492 gf (the records were run off 7fb2224). A fresh two-engine run on current main picks up the Costa-Silva-2020 gf (A(S)~7.486); S stays CURATION-OWED (gf-scale floor, RYA-161) either way.
- Ca: PROMOTED CURATION-OWED -> PASS at 6.372 (+0.072) by the ratified RYA-561 three-gate rule — 534-validated atom, d_ref within 0.1, and a REAL cross-engine dCE +0.016 from BOTH engines over 4 lines. SCALE CAVEAT: 1D-NLTE value vs 3D-NLTE reference; un-applied 3D term folded into the offset — 3D is Fe-only so far (Magic-2013, RYA-553); small for weak lines like Ca 6122/6162, not zero (RYA-399/586).
- Cr: reported = Cr I gf-floor 6.022 (+0.402 vs Asplund, the RYA-398 graded-pool CANARY — stays CURATION-OWED at floor, NOT PASS). Cr II 5.676 is DIAGNOSTIC-ONLY (RYA-240 ratified exclusion — COG/saturation artifact; enforced by the engine_selection guard, RYA-558) and the two-engine Cr I synthesis 5.654 sits near-anchor — both are diagnostics, never the reported value. Promotion of Cr II needs clean unsaturated weak lines (future decision), not the blind floor.
- Ti: production NLTE = Engine-A Mallinson-2024 (RYA-545). The Engine-B Gerber Ti (+0.221) ships atom.ti503b and is a strict xfail (RYA-548) — recorded as diagnostic, not applied to the reported value.