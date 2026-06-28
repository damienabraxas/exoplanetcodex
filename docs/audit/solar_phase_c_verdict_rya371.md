# RYA-371 Phase C — Solar 27-element verdict table (RYA-239 retry)

_Generated 2026-06-27 by scripts/phase_c_verdict_rya371.py. Validate-don't-tune: classification only, no correction fitted to the anchor._

## Verdict counts

- **PASS**: 3
- **NLTE-OWED**: 1
- **CURATION-OWED**: 18
- **DATA-GAP**: 4

## Per-element table

| El | Asplund21 | A(meas) | Delta | sigma | n | NLTE | Verdict | Channel |
|----|----------:|--------:|------:|------:|--:|:----:|:--------|:--------|
| O | 8.69 | 8.735 | +0.045 | 0.01 | 3 | - | **PASS** | synthesis: O I 777 (primary) + [O I] 6300 (cross-check) |
| C | 8.46 | 8.491 | +0.031 | 0.05 | 5 | - | **PASS** | synthesis: CH G-band + C I 5052 + C2 Swan (C I 5380 BAD_FIT-excluded) |
| N | 7.83 | 7.774 | -0.056 | 0.96 | 3 | - | **NLTE-OWED** | synthesis: N I 8216 + CN red; NH 3360 primary |
| Mg | 7.55 |  |  |  |  | wired | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Si | 7.51 | 7.888 | +0.378 | 0.36 | 7 | wired+3D | **CURATION-OWED** | EW: 7 curated line(s), graded-gf (RYA-395/398) |
| Fe | 7.46 | 7.516 | +0.056 | 0.14 | 62 | - | **PASS** | EW: 62 Fe I + 3 Fe II, NLTE-wired (Bergemann MPIA) |
| S | 7.12 | 7.753 | +0.633 | 0.37 | 2 | wired | **CURATION-OWED** | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Al | 6.43 |  |  |  |  | wired | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Ca | 6.30 | 6.324 | +0.024 | 0.12 | 2 | wired | **CURATION-OWED** | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Na | 6.24 |  |  |  |  | wired | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Ni | 6.20 | 6.946 | +0.746 | 0.51 | 2 | - | **CURATION-OWED** | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Cr | 5.62 | 6.022 | +0.402 | 0.60 | 7 | wired+3D | **CURATION-OWED** | EW: 7 curated line(s), graded-gf (RYA-395/398) |
| Mn | 5.42 |  |  |  |  | wired | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| P | 5.41 |  |  |  |  | - | **DATA-GAP** | none on present ground set |
| K | 5.07 |  |  |  |  | - | **DATA-GAP** | no curated solar lines in present set |
| Ti | 4.97 | 5.471 | +0.501 | 0.93 | 10 | wired+3D | **CURATION-OWED** | EW: 10 curated line(s), graded-gf (RYA-395/398) |
| Co | 4.94 |  |  |  |  | - | **DATA-GAP** | no curated solar lines in present set |
| Cu | 4.18 |  |  |  |  | - | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| V | 3.90 |  |  |  |  | - | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Sc | 3.14 |  |  |  |  | - | **DATA-GAP** | no curated solar lines in present set |
| Sr | 2.83 | 4.961 | +2.131 |  | 1 | wired | **CURATION-OWED** | EW: 1 curated line(s), low-confidence (RYA-395/398) |
| Zr | 2.59 |  |  |  |  | - | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Ba | 2.27 |  |  |  |  | wired | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Y | 2.21 |  |  |  |  | - | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |
| Li | 1.05 | 0.727 | -0.323 |  | 1 | - | **CURATION-OWED** | EW: Li I 6707 (single line, UPPER LIMIT) |
| Eu | 0.52 |  |  |  |  | - | **CURATION-OWED** | EW present; no independent-gf line survives the graded cull |

## Remaining-work map

### PASS (3)

- **O** — cross-arm AGREE; O I 777 Amarsi-2019 3D-NLTE, [O I] Caffau-2015 3D anchor — measured 8.74 vs Asplund 8.69 (+0.05). RYA-455.
- **C** — C I Amarsi-2019 3D-NLTE -> 8.46; CH 8.49 (3D-offset-owed); C I 5380 formally excluded ew_integrity=BAD_FIT (RYA-458); surviving cross-arm spread 0.054 on 5 indicators. [C I 5380 EXCLUDED (ew_integrity=BAD_FIT, espresso:CI_5380); cross-arm spread 0.149->0.054 on 5 surviving indicators.]
- **Fe** — A(Fe I) NLTE 7.516 vs Asplund 7.46 (+0.056); ionization-balance gated, scatter 0.139 = honest floor (RYA-407). Documented +0.05 1D/3D scale offset (RYA-336), not the verdict.

### NLTE-OWED (1)

- **N** — N I 8216 LTE 7.99 (NLTE owed -> N I grid, RYA-369, would pull toward 7.83); NH 3360 primary channel is a DATA-GAP (UVES-blue under-SNR). cross-arm FLAGGED-DISAGREEMENT (reported, not averaged) (spread 0.96).

### CURATION-OWED (18)

- **Mg** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Mg_Amarsi2020_PySME.csv).
- **Si** — A(X) 7.888 vs Asplund 7.51 (+0.378), sigma 0.36 — gross offset removed by the blind cull, gf-scale residual survives on the graded pool → escalate to RYA-161/162 (do NOT tune).
- **S** — A(X) 7.753 vs Asplund 7.12 (+0.633) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Al** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Al_Amarsi2020_PySME.csv).
- **Ca** — A(X) 6.324 vs Asplund 6.30 (+0.024) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Na** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Na_Amarsi2020_PySME.csv).
- **Ni** — A(X) 6.946 vs Asplund 6.20 (+0.746) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Cr** — A(X) 6.022 vs Asplund 5.62 (+0.402), sigma 0.60 — gross offset removed by the blind cull, gf-scale residual survives on the graded pool → escalate to RYA-161/162 (do NOT tune).
- **Mn** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Mn_Bergemann_MPIA.csv).
- **Ti** — A(X) 5.471 vs Asplund 4.97 (+0.501), sigma 0.93 — gross offset removed by the blind cull, gf-scale residual survives on the graded pool → escalate to RYA-161/162 (do NOT tune).
- **Cu** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged).
- **V** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged).
- **Sr** — A(X) 4.961 vs Asplund 2.83 (+2.131) on 1 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Zr** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged).
- **Ba** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Ba_Korotin2015.csv).
- **Y** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged).
- **Li** — CN-blended UPPER LIMIT (RYA-103/458, ew_integrity disposition=UPPER_LIMIT); A(Li) 0.73 is a LTE lower bound, not a clean determination. A clean low value here would be a RED FLAG (CN deblend not applied).
- **Eu** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged). Eu II 6645 EW 6.8 mA, ew_integrity disposition=RECOVERED (RYA-102/458 HFS-summing).

### DATA-GAP (4)

- **P** — usable P lines are FUV/near-UV-hard, below the ~300 nm atmospheric floor; needs HST/STIS (RYA-119). Do not force.
- **K** — K not in the canonical solar EW pool; extraction owed. no NLTE grid.
- **Co** — Co not in the canonical solar EW pool; extraction owed. no NLTE grid.
- **Sc** — Sc not in the canonical solar EW pool; extraction owed. no NLTE grid.

## Off-Sun — DEFERRED to RYA-348

A green solar run validates the MACHINERY, not target-transferability. These cannot be validated on the Sun and are explicitly deferred:

- per-star broadening sourcing (vmac/vsini)
- temperature-dependent NLTE (e.g. C I LTE-on-Sun but NLTE-on-Procyon)
- F-star / cool-star line-list adequacy

_Not "55 Cnc-ready" — that is RYA-348 and beyond._
