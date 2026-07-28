# RYA-527 — re-emitted solar verdict + PROPOSED gold v3 (two-engine floor)

_REVIEW ARTIFACT — NOT frozen. Freeze = `promote_solar_reference.py --apply` (Ryan). Fe = 7.466 (3D, RYA-553); two-engine 7.580 is the RYA-525 diagnostic only._

Verdict counts: {'PASS': 5, 'NLTE-OWED': 1, 'CURATION-OWED': 20}

## Proposed v3 vs v2 vs Asplund (tiered), with the two-engine per-element record

| El | Asp | v2 gold | v3 proposed | dAsp | verdict | engA | engB | selected | dCE | source |
|----|-----|---------|-------------|------|---------|------|------|----------|-----|--------|
| O I | 8.69 | 8.735 | 8.735 | +0.045 | PASS | - | 8.73 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| C I | 8.46 | 8.491 | 8.491 | +0.031 | PASS | - | 8.491 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Fe I | 7.46 | 7.516 | 7.466 **NEW** | +0.006 | PASS | 7.475 | 7.495 | A_1dnlte,B_synth | 0.084 | EW Fe I ionization-gated, 3D-corrected (RYA-406/407/553); two-engine 7.580 is the RYA-525 cross-engine diagnostic ONLY |
| Mn I | 5.42 | 5.47 | 5.466 **NEW** | +0.046 | PASS | - | 5.466 | B_synth | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| K I | 5.07 | 5.099 | 5.099 | +0.029 | PASS | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| N I | 7.83 | - | 8.202 **NEW** | +0.372 | NLTE-OWED | - | - | - | - | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Mg I | 7.55 | - | 7.614 **NEW** | +0.064 | CURATION-OWED | - | 7.614 | B_synth | - | two-engine synthesis floor (I, B_synth) |
| Si I | 7.51 | 7.888 | 7.639 **NEW** | +0.129 | CURATION-OWED | 8.06 | 7.633 | B_synth | -0.756 | two-engine synthesis floor (I, B_synth) |
| S I | 7.12 | - | 7.369 **NEW** | +0.249 | CURATION-OWED | 7.369 | 7.41 | A_1dnlte | 0.041 | two-engine synthesis floor (I, A_1dnlte) |
| Al I | 6.43 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |
| Ca I | 6.3 | - | 6.372 **NEW** | +0.072 | CURATION-OWED | 6.344 | 6.327 | A_1dnlte,B_synth | 0.016 | two-engine synthesis floor (I, A_1dnlte,B_synth) |
| Na I | 6.24 | - | 6.369 **NEW** | +0.129 | CURATION-OWED | 6.37 | 6.248 | A_1dnlte | -0.121 | two-engine synthesis floor (I, A_1dnlte) |
| Ni I | 6.2 | - | 6.253 **NEW** | +0.053 | CURATION-OWED | 7.297 | 6.2 | B_synth | -1.253 | two-engine synthesis floor (I, B_synth) |
| Cr II | 5.62 | 6.022 | 5.676 **NEW** | +0.056 | CURATION-OWED | 8.354 | 5.676 | B_synth | -2.677 | two-engine synthesis floor (II, B_synth) |
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
| Li I | 1.05 | 0.727 | 1.409 **NEW** | +0.359 | CURATION-OWED | 0.727 | 1.409 | B_synth | 0.682 | two-engine synthesis floor (I, B_synth) |
| Eu I | 0.52 | - | - | - | CURATION-OWED | - | - | - | - | phase_c (owed) |

## RYA-524 reconciliation (S / Sr / N)

- **S**: STALE (owed, no value / pre-492 gf) -> CURATION-OWED, two-engine 7.369. RYA-492 Costa-Silva gf; gf barely moves it (+0.004 on GES) so the +0.37 vs Asplund is a gf-scale floor (RYA-161), not a line-ID error — stays owed.
- **Sr**: WRONG-SPECIES (Sr I raw-EW +2.13 artifact) -> CURATION-OWED, Sr II 2.769 (two-engine synth, RYA-551). Discarded the Sr I raw-EW artifact; Sr II 4077 synthesis (INSPECT NLTE). ~0.05-0.1 dex near-UV systematic -> refinement owed, not a clean PASS.
- **N**: STALE/unwired (NLTE-OWED, no value) -> 8.202 / NLTE-OWED (Kittpeak); RYA-526 grid registered. N I NLTE grid registered (RYA-369/526). See the N WIRING FLAG: the KP channel does not auto-apply the grid delta in phase_c yet; the +0.36 residual is a data-channel/gf item, not an NLTE debt.

## Honest flags

- N: base phase_c = 8.202 / NLTE-OWED (Kittpeak red multiplets). RYA-526 registered the N I NLTE grid; the phase_c KP channel does not auto-apply its (~-0.014) delta, so N still reads NLTE-OWED here — a properly-wired verdict moves it to ~8.188 CURATION-OWED (data-channel/gf floor, not an NLTE debt). WIRING FLAG.
- Fe: reported 7.466 (3D, RYA-553). Two-engine 7.580 is the RYA-525 cross-engine diagnostic ONLY (per-line winner-combine biases high) — NOT the verdict (Ryan 2026-07-16).
- S: the committed two-engine record engineA=7.369 is on the PRE-RYA-492 gf (the records were run off 7fb2224). A fresh two-engine run on current main picks up the Costa-Silva-2020 gf (A(S)~7.486); S stays CURATION-OWED (gf-scale floor, RYA-161) either way.
- Ti: production NLTE = Engine-A Mallinson-2024 (RYA-545). The Engine-B Gerber Ti (+0.221) ships atom.ti503b and is a strict xfail (RYA-548) — recorded as diagnostic, not applied to the reported value.