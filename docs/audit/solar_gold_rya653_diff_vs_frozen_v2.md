# Rebuilt solar gold candidate vs frozen v2 — per-cell diff

Candidate: `data/reference/solar/solar_abundances_corrected_candidate_rya653.csv` (rebuilt from `data/audit/cno_synthesis/solar_phase_c_verdict.json`).
Frozen reference: `v2` via `pipeline.data_namespace.read_solar_reference()` (the CURRENT pointer). **This diff is a proposal, not a freeze** — promoting it is `scripts/promote_solar_reference.py --apply` and a ratification decision (RYA-527).

**29 cell(s) differ.**

| element | field | frozen v2 | rebuilt candidate |
|---|---|---|---|
| Fe | A_X | 7.516 | 7.466 |
| Fe | A_X_nlte | 7.516 | 7.466 |
| Mg | note | no independent-gf line survives the graded cull | EW present; no independent-gf line survives the graded cull |
| Ca | note | LOW_CONFIDENCE / thin graded pool (+0.02) | EW: 2 curated line(s); value HELD at gold tier 'owed' (RYA-522) — not a graded-cull blank |
| Ti | note | LOW_CONFIDENCE / thin graded pool (+0.50) | EW: 10 curated line(s); value HELD at gold tier 'owed' (RYA-522) — not a graded-cull blank |
| Ni | note | LOW_CONFIDENCE / thin graded pool (+0.75) | EW: 2 curated line(s); value HELD at gold tier 'owed' (RYA-522) — not a graded-cull blank |
| Na | note | LOW_CONFIDENCE / thin graded pool (+0.02) | EW: 2 curated line(s); value HELD at gold tier 'owed' (RYA-522) — not a graded-cull blank |
| P | note | near-IR multiplet, gf-limited (RYA-460) | near-IR multiplet, gf-limited (RYA-460) — A(X) 6.610 HELD at tier 'owed' (RYA-522), not frozen |
| S | method_scale | EW 1D-LTE/NLTE | synthesis |
| S | note | LOW_CONFIDENCE / thin graded pool (+0.63) | synthesis: Costa-Silva gf (RYA-492) — A(X) 7.486 HELD at tier 'owed' (RYA-522), not frozen |
| N | verdict | NLTE-OWED | CURATION-OWED |
| N | note | N I red multiplets; +0.37 owed NLTE (RYA-369) | N I red multiplets; +0.37 owed NLTE (RYA-369) — A(X) 8.188 HELD at tier 'owed' (RYA-522), not frozen |
| Co | verdict | CURATION-OWED | PASS |
| Co | method_scale | atlas 1D | synthesis |
| Co | n_lines | 1 | 5 |
| Co | note | blue-edge, SNR-limited — not trusted (RYA-460) | synthesis: Co I red HFS (RYA-564) — A(X) 4.965 HELD at tier 'owed' (RYA-522), not frozen |
| Al | note | LOW_CONFIDENCE / thin graded pool (+0.98) | EW: 1 curated line(s); value HELD at gold tier 'owed' (RYA-522) — not a graded-cull blank |
| Ba | method_scale | EW 1D-LTE/NLTE | synthesis |
| Ba | n_lines | 0 | 1 |
| Ba | note | no independent-gf line survives the graded cull | synthesis: Ba II 5853 HFS (RYA-559) — A(X) 2.410 HELD at tier 'owed' (RYA-522), not frozen |
| Y | note | no independent-gf line survives the graded cull | EW present; no independent-gf line survives the graded cull |
| V | note | HFS-resolved synthesis (RYA-411/466/473) | HFS-resolved synthesis (RYA-411/466/473) — A(X) 3.917 HELD at tier 'owed' (RYA-522), not frozen |
| Cu | note | HFS-resolved synthesis (RYA-411/466/473) | HFS-resolved synthesis (RYA-411/466/473) — A(X) 4.345 HELD at tier 'owed' (RYA-522), not frozen |
| Mn | A_X | 5.47 | 5.466 |
| Mn | A_X_nlte | 5.47 | 5.466 |
| Eu | note | no independent-gf line survives the graded cull | EW present; no independent-gf line survives the graded cull |
| Zr | note | no independent-gf line survives the graded cull | EW present; no independent-gf line survives the graded cull |
| Sr | method_scale | EW (SUSPECT) | EW 1D-LTE/NLTE |
| Sr | note | +2.13 — NOT a gf-floor; saturated-line-on-flat-COG signature (the RYA-520 disease) → saturation-trace owed | EW: 1 curated line(s); value HELD at gold tier 'owed' (RYA-522) — not a graded-cull blank |

