# RYA-669 — RYA-527 Phase 2: proposed gold v4 vs frozen v3 vs Asplund

_**REVIEW ARTIFACT — NOTHING IS FROZEN.** The freeze is a separate, Ryan-ratified ticket running `promote_solar_reference.py --apply`. Gold v3 is untouched by this run._

Verdict counts: {'PASS': 6, 'CURATION-OWED': 20}

`v3 gold` is blank wherever the RYA-522 tier is `owed` — that is the tier holding the value unfrozen, not a missing measurement.

## Proposed v4 vs frozen v3 vs Asplund (tiered), with the two-engine record

| El | Asp | v3 gold | tier | v4 proposed | Δ(v4−v3) | ΔAsp | verdict | engA | engB | selected | dCE | source |
|----|-----|---------|------|-------------|----------|------|---------|------|------|----------|-----|--------|
| O I | 8.69 | 8.735 | gold | 8.735 | +0.000 | +0.045 | PASS | — | 8.73 | B_synth | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| C I | 8.46 | 8.491 | gold | 8.491 | +0.000 | +0.031 | PASS | — | 8.491 | B_synth | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Fe I | 7.46 | 7.466 | gold | 7.466 | +0.000 | +0.006 | PASS | 7.475 | 7.497 | A_1dnlte,B_synth | 0.088 | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Mn I | 5.42 | 5.466 | gold | 5.466 | +0.000 | +0.046 | PASS | — | 5.466 | B_synth | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| K I | 5.07 | 5.099 | gold | 5.099 | +0.000 | +0.029 | PASS | — | — | — | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Co I | 4.94 | — | owed | 4.965 | — | +0.025 | PASS | — | — | — | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| N I | 7.83 | — | owed | 8.188 | — | +0.358 | CURATION-OWED | — | — | — | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Mg I | 7.55 | — | owed | 7.62 | — | +0.070 | CURATION-OWED | — | 7.62 | B_synth | — | two-engine synthesis floor (I, B_synth) |
| Si I | 7.51 | 7.888 | gf_floor | 7.64 **MOVED** | -0.248 | +0.130 | CURATION-OWED | 8.06 | 7.633 | B_synth | -0.756 | two-engine synthesis floor (I, B_synth) |
| S I | 7.12 | — | owed | 7.369 | — | +0.249 | CURATION-OWED | 7.369 | 7.427 | A_1dnlte | 0.058 | two-engine synthesis floor (I, A_1dnlte) |
| Al I | 6.43 | — | owed | — | — | — | CURATION-OWED | — | — | — | — | phase_c (owed) |
| Ca I | 6.3 | — | owed | 6.355 | — | +0.055 | CURATION-OWED | 6.344 | 6.306 | A_1dnlte,B_synth | -0.003 | two-engine synthesis floor (I, A_1dnlte,B_synth) |
| Na I | 6.24 | — | owed | 6.369 | — | +0.129 | CURATION-OWED | 6.37 | 6.246 | A_1dnlte | -0.124 | two-engine synthesis floor (I, A_1dnlte) |
| Ni I | 6.2 | — | owed | 6.252 | — | +0.052 | CURATION-OWED | 7.297 | 6.197 | B_synth | -1.254 | two-engine synthesis floor (I, B_synth) |
| Cr I | 5.62 | 6.022 | gf_floor | 5.661 **MOVED** | -0.361 | +0.041 | CURATION-OWED | 6.036 | 5.645 | B_synth | -1.051 | two-engine synthesis floor (I, B_synth) |
| P I | 5.41 | — | owed | 6.61 | — | +1.200 | CURATION-OWED | — | — | — | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Ti I | 4.97 | — | owed | 4.879 | — | -0.091 | CURATION-OWED | 6.401 | 4.863 | B_synth | -1.326 | two-engine synthesis floor (I, B_synth) |
| Cu I | 4.18 | — | owed | 4.345 | — | +0.165 | CURATION-OWED | — | 4.345 | B_synth | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| V I | 3.9 | — | owed | 3.917 | — | +0.017 | CURATION-OWED | — | 3.917 | B_synth | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Sc II | 3.14 | 3.203 | gold | 3.203 | +0.000 | +0.063 | CURATION-OWED | — | — | — | — | ratified/dedicated channel (phase_c); two-engine = diagnostic |
| Sr II | 2.83 | — | owed | 2.759 | — | -0.071 | CURATION-OWED | — | 2.759 | B_synth | — | two-engine synthesis floor (II, B_synth) |
| Zr I | 2.59 | — | owed | — | — | — | CURATION-OWED | — | — | — | — | phase_c (owed) |
| Ba I | 2.27 | — | owed | 2.41 | — | +0.140 | CURATION-OWED | — | — | — | — | phase_c (owed) |
| Y I | 2.21 | — | owed | — | — | — | CURATION-OWED | — | — | — | — | phase_c (owed) |
| Li I | 1.05 | 0.727 | upper_limit | 0.727 | +0.000 | -0.323 | CURATION-OWED | 0.727 | 1.155 | B_synth | 0.428 | HELD — RYA-563 UPPER_LIMIT disposition — the reference-blind floor may NEVER emit a point value for this element; two-engine value recorded as diagnostic only |
| Eu I | 0.52 | — | owed | — | — | — | CURATION-OWED | — | — | — | — | phase_c (owed) |

## Elements whose value moved > 0.01 dex vs frozen v3

| element | v3 | v4 | Δ | why |
|---|---|---|---|---|
| Si | 7.888 | 7.64 | -0.248 | two-engine synthesis floor (I, B_synth) |
| Cr | 6.022 | 5.661 | -0.361 | two-engine synthesis floor (I, B_synth) |

## The four species-adoption decisions — Ryan decides, nothing adopted

### D.1 — Sr: adopt Sr II 2.759 into gold v4?

- `fresh_value`: 2.759
- `reference_value`: 2.759
- `reproduces_within_0.02`: True
- `near_uv_red_chi2`: 78.27
- `near_uv_reliable`: True
- `current_verdict_channel`: None
- **Recommendation:** ADOPT Sr II 2.759 — reproduces within ±0.02 of the RYA-643 corrected 2.759.
- ⚠ near-UV fit red-chi2 78.27 — high, and flagged as the open near-UV item in RYA-643 (red-chi2 78-180 across the blue channels). The line is marked reliable on dEW/dA grounds, not on chi2. Adoption inherits that systematic; it is a known ~0.05-0.1 dex near-UV floor, not noise.

### D.2 — Co: which Co value goes into gold v4?

- `fresh_value`: None
- `no_fresh_value_because`: Co produced NO two-engine record at all. It has no EW-pool lines and no synth-v2 lines, and `_dedicated_engine_B()` in rya527_two_engine_run.py wires C/O/Mn/Cu/V/Sr/Zr/Mg but NOT the RYA-564 Co red-line synthesis (data/results/co_synthesis_rya564.json). So Phase 2 CANNOT arbitrate this split — it produced no third number to weigh. Wiring Co into the dedicated Engine-B set is the prerequisite, and is not in this ticket's scope.
- `candidates`: {'phase_c (RYA-564)': 4.965, 'RYA-643 corrected re-run': 4.96, 'fresh two-engine': None}
- **Recommendation:** Phase 2 cannot break this tie — it produced no Co number (see above). On the two existing candidates the split is 0.005 dex, inside every gate and below the reported precision, so it is a provenance choice, not a measurement one: 4.960 comes from the run with the RYA-643 rest-frame/gsig defect fixed, which is the better-founded of the two. Ryan picks.
- ⚠ Co is verdict PASS at tier `owed` ⇒ v3 freezes NO value for it (RYA-665). Whichever number is picked, it stays HELD until the tier moves — adopting a value here does not by itself freeze one.

### D.3 — Ba: fire RYA-581 deblend BEFORE the v4 freeze, or freeze 2.410 HELD-with-caveat?

- `fresh_value`: None
- `no_fresh_value_because`: Ba produced NO two-engine record either, for the same reason as Co: no EW-pool or synth-v2 lines, and the RYA-559 Ba II 5853 synthesis (data/results/solar_ba_synthesis_rya559.json) is not wired into `_dedicated_engine_B()`. Ba's gate 3 therefore stays UNEVALUABLE after Phase 2 — which the ticket expected the re-run to fix, and it does not.
- `current_verdict_channel`: 2.41
- `clean_cross_check`: [2.187, 2.231]
- **Recommendation:** No recommendation on timing — this is a scheduling call, not a measurement one. The measurement fact: 2.410 is blend-inflated by ~+0.15 against a clean cross-check at 2.187/2.231, and RYA-581 exists to replace it. Freezing 2.410 into v4 would immortalise a number already known to be high; deferring costs one ticket.
- ⚠ Ba is `owed` tier in v3 ⇒ blank A_X, so nothing is frozen today either way. The urgency is about v4, not v3.

### D.4-Ca — Ca: promote Ca to PASS on the fresh gate 3?

- `value`: 6.324
- `reference`: 6.3
- `gate1_atom_validated`: True
- `gate2_within_tol`: True
- `gate3_state`: OK
- `cross_engine_delta`: -0.003
- `promoted_by_ratified_rule`: True
- `blocker`: none — promotes under the ratified three gates
- `gate3_still_provisional`: True
- **Recommendation:** PROMOTE Ca — clears all three ratified gates on a FRESH cross-engine delta (-0.003) computed this run from both engines over real solar data. The report's blanket PROVISIONAL stamp is spurious here (see the gate-3 section): it comes from cross-CHANNEL disagreements on other elements, not from anything about Ca's delta or the artifact's age.

### D.4-Na — Na: promote Na to PASS on the fresh gate 3?

- `value`: 6.264
- `reference`: 6.24
- `gate1_atom_validated`: True
- `gate2_within_tol`: True
- `gate3_state`: FAILED
- `cross_engine_delta`: -0.124
- `promoted_by_ratified_rule`: False
- `blocker`: gate 3: |dCE| = 0.124 > 0.1
- `gate3_still_provisional`: True
- **Recommendation:** DO NOT promote Na — gate 3: |dCE| = 0.124 > 0.1. This is a FRESH answer: RYA-664 cleared its gate 1, and gate 3 is now decided on a current delta rather than deferred.

