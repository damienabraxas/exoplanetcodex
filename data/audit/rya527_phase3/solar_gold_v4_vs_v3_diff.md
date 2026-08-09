# Rebuilt solar gold candidate vs frozen v3 — per-cell diff

Candidate: `data/audit/rya527_phase3/solar_abundances_v4_candidate.csv` (rebuilt from `data/audit/cno_synthesis/solar_phase_c_verdict.json`).
Frozen reference: `v3` via `pipeline.data_namespace.read_solar_reference()` (the CURRENT pointer). **This diff is a proposal, not a freeze** — promoting it is `scripts/promote_solar_reference.py --apply` and a ratification decision (RYA-527).

**57 cell(s) differ.**

| element | field | frozen v3 | rebuilt candidate |
|---|---|---|---|
| Fe | method_scale | 1D-NLTE (Fe I) | 3D-NLTE (Fe I, Magic 2013) |
| Fe | scale_state | <absent> | 3D-NLTE |
| Fe | corrections_applied | <absent> | ["1D_3D_solar_Fe_Magic2013"] |
| Fe | note | our 1D-NLTE runs ~+0.05 above Asplund 3D-true 7.46 (RYA-336) — documented offset, NOT a discrepancy | reported anchor on the true 3D scale — the tabulated Magic-2013 1D→3D solar Fe offset is APPLIED (RYA-553); scale carried as DATA in `scale_state`, never inferred from this prose (RYA-681) |
| C | scale_state | <absent> |  |
| C | corrections_applied | <absent> | [] |
| O | scale_state | <absent> |  |
| O | corrections_applied | <absent> | [] |
| Mg | scale_state | <absent> |  |
| Mg | corrections_applied | <absent> | [] |
| Si | scale_state | <absent> |  |
| Si | corrections_applied | <absent> | [] |
| Ca | scale_state | <absent> |  |
| Ca | corrections_applied | <absent> | [] |
| Ti | scale_state | <absent> |  |
| Ti | corrections_applied | <absent> | [] |
| Ni | scale_state | <absent> |  |
| Ni | corrections_applied | <absent> | [] |
| Na | scale_state | <absent> |  |
| Na | corrections_applied | <absent> | [] |
| P | scale_state | <absent> |  |
| P | corrections_applied | <absent> | [] |
| S | scale_state | <absent> |  |
| S | corrections_applied | <absent> | [] |
| N | scale_state | <absent> |  |
| N | corrections_applied | <absent> | [] |
| Co | scale_state | <absent> |  |
| Co | corrections_applied | <absent> | [] |
| Co | note | synthesis: Co I red HFS (RYA-564) — A(X) 4.965 HELD at tier 'owed' (RYA-522), not frozen | synthesis: Co I red HFS (RYA-564) — A(X) 4.960 HELD at tier 'owed' (RYA-522), not frozen |
| Cr | scale_state | <absent> |  |
| Cr | corrections_applied | <absent> | [] |
| Al | scale_state | <absent> |  |
| Al | corrections_applied | <absent> | [] |
| K | scale_state | <absent> |  |
| K | corrections_applied | <absent> | [] |
| Ba | verdict | CURATION-OWED | PASS |
| Ba | scale_state | <absent> |  |
| Ba | corrections_applied | <absent> | [] |
| Ba | note | synthesis: Ba II 5853 HFS (RYA-559) — A(X) 2.410 HELD at tier 'owed' (RYA-522), not frozen | synthesis: Ba II 5853 in-window deblend (RYA-581) — A(X) 2.237 HELD at tier 'owed' (RYA-522), not frozen |
| Y | scale_state | <absent> |  |
| Y | corrections_applied | <absent> | [] |
| V | scale_state | <absent> |  |
| V | corrections_applied | <absent> | [] |
| Cu | scale_state | <absent> |  |
| Cu | corrections_applied | <absent> | [] |
| Mn | scale_state | <absent> |  |
| Mn | corrections_applied | <absent> | [] |
| Sc | scale_state | <absent> |  |
| Sc | corrections_applied | <absent> | [] |
| Li | scale_state | <absent> |  |
| Li | corrections_applied | <absent> | [] |
| Eu | scale_state | <absent> |  |
| Eu | corrections_applied | <absent> | [] |
| Zr | scale_state | <absent> |  |
| Zr | corrections_applied | <absent> | [] |
| Sr | scale_state | <absent> |  |
| Sr | corrections_applied | <absent> | [] |

