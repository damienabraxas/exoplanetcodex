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
| C | 8.46 | 8.491 | +0.031 | 0.15 | 6 | - | **PASS** | synthesis: CH G-band + C I 5052/5380 + C2 Swan |
| N | 7.83 | 7.774 | -0.056 | 0.96 | 3 | - | **NLTE-OWED** | synthesis: N I 8216 + CN red; NH 3360 primary |
| Mg | 7.55 | 8.424 | +0.874 | 1.04 | 2 | wired | **CURATION-OWED** | EW: 2 line(s) |
| Si | 7.51 | 8.060 | +0.550 | 0.61 | 15 | wired+3D | **CURATION-OWED** | EW: 15 line(s) |
| Fe | 7.46 | 7.516 | +0.056 | 0.14 | 62 | - | **PASS** | EW: 62 Fe I + 3 Fe II, NLTE-wired (Bergemann MPIA) |
| S | 7.12 | 7.386 | +0.266 |  | 1 | wired | **CURATION-OWED** | EW: 1 line(s) |
| Al | 6.43 |  |  |  |  | wired | **CURATION-OWED** | EW present, not wired into production A(X) |
| Ca | 6.30 | 6.382 | +0.082 | 0.56 | 4 | wired | **CURATION-OWED** | EW: 4 line(s) |
| Na | 6.24 |  |  |  |  | wired | **CURATION-OWED** | EW present, not wired into production A(X) |
| Ni | 6.20 | 7.297 | +1.097 | 0.97 | 3 | - | **CURATION-OWED** | EW: 3 line(s) |
| Cr | 5.62 | 5.906 | +0.286 | 1.09 | 3 | wired+3D | **CURATION-OWED** | EW: 3 line(s) |
| Mn | 5.42 |  |  |  |  | wired | **CURATION-OWED** | EW present, not wired into production A(X) |
| P | 5.41 |  |  |  |  | - | **DATA-GAP** | none on present ground set |
| K | 5.07 |  |  |  |  | - | **DATA-GAP** | no curated solar lines in present set |
| Ti | 4.97 | 6.374 | +1.404 | 0.85 | 6 | wired+3D | **CURATION-OWED** | EW: 6 line(s) |
| Co | 4.94 |  |  |  |  | - | **DATA-GAP** | no curated solar lines in present set |
| Cu | 4.18 |  |  |  |  | - | **CURATION-OWED** | EW present, not wired into production A(X) |
| V | 3.90 |  |  |  |  | - | **CURATION-OWED** | EW present, not wired into production A(X) |
| Sc | 3.14 |  |  |  |  | - | **DATA-GAP** | no curated solar lines in present set |
| Sr | 2.83 |  |  |  |  | wired | **CURATION-OWED** | EW present, not wired into production A(X) |
| Zr | 2.59 |  |  |  |  | - | **CURATION-OWED** | EW present, not wired into production A(X) |
| Ba | 2.27 |  |  |  |  | wired | **CURATION-OWED** | EW present, not wired into production A(X) |
| Y | 2.21 |  |  |  |  | - | **CURATION-OWED** | EW present, not wired into production A(X) |
| Li | 1.05 | 0.727 | -0.323 |  | 1 | - | **CURATION-OWED** | EW: Li I 6707 (single line, upper limit) |
| Eu | 0.52 |  |  |  |  | - | **CURATION-OWED** | EW present, not wired into production A(X) |

## Remaining-work map

### PASS (3)

- **O** — cross-arm AGREE; O I 777 Amarsi-2019 3D-NLTE, [O I] Caffau-2015 3D anchor — measured 8.74 vs Asplund 8.69 (+0.05). RYA-455.
- **C** — C I Amarsi-2019 3D-NLTE -> 8.46; CH 8.49 (3D-offset-owed); ESPRESSO C I 5380 chi2r~103 flagged outlier, NOT averaged. spread 0.149.
- **Fe** — A(Fe I) NLTE 7.516 vs Asplund 7.46 (+0.056); ionization-balance gated, scatter 0.139 = honest floor (RYA-407). Documented +0.05 1D/3D scale offset (RYA-336), not the verdict.

### NLTE-OWED (1)

- **N** — N I 8216 LTE 7.99 (NLTE owed -> N I grid, RYA-369, would pull toward 7.83); NH 3360 primary channel is a DATA-GAP (UVES-blue under-SNR). cross-arm FLAGGED-DISAGREEMENT (reported, not averaged) (spread 0.96).

### CURATION-OWED (18)

- **Mg** — A(X) 8.424 vs Asplund 7.55 (+0.874), sigma 1.04 — gf-/blend-limited pool (RYA-395/398). NLTE grid Mg_Amarsi2020_PySME.csv wired but measured lines fall outside its node coverage (flag NLTE_unavailable).
- **Si** — A(X) 8.060 vs Asplund 7.51 (+0.550), sigma 0.61 — gf-/blend-limited pool (RYA-395/398). NLTE grid Si_Amarsi2020_PySME.csv wired but measured lines fall outside its node coverage (flag NLTE_unavailable).
- **S** — A(X) 7.386 vs Asplund 7.12 (+0.266), sigma nan — gf-/blend-limited pool (RYA-395/398). NLTE grid S_Amarsi2025_PySME.csv wired but measured lines fall outside its node coverage (flag NLTE_unavailable).
- **Al** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). NLTE grid available (Al_Amarsi2020_PySME.csv).
- **Ca** — A(X) 6.382 vs Asplund 6.30 (+0.082), sigma 0.56 — gf-/blend-limited pool (RYA-395/398). NLTE grid Ca_Mashonkina2017.csv wired but measured lines fall outside its node coverage (flag NLTE_unavailable).
- **Na** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). NLTE grid available (Na_Amarsi2020_PySME.csv).
- **Ni** — A(X) 7.297 vs Asplund 6.20 (+1.097), sigma 0.97 — gf-/blend-limited pool (RYA-395/398). no NLTE grid (LTE-flagged).
- **Cr** — A(X) 5.906 vs Asplund 5.62 (+0.286), sigma 1.09 — gf-/blend-limited pool (RYA-395/398). NLTE grid Cr_Bergemann2010_MPIA.csv wired but measured lines fall outside its node coverage (flag NLTE_unavailable).
- **Mn** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). NLTE grid available (Mn_Bergemann_MPIA.csv).
- **Ti** — A(X) 6.374 vs Asplund 4.97 (+1.404), sigma 0.85 — gf-/blend-limited pool (RYA-395/398). NLTE grid Ti_Bergemann2011_MPIA.csv wired but measured lines fall outside its node coverage (flag NLTE_unavailable).
- **Cu** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). no NLTE grid (would be LTE-flagged).
- **V** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). no NLTE grid (would be LTE-flagged).
- **Sr** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). NLTE grid available (Sr_Bergemann2012_INSPECT.csv).
- **Zr** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). no NLTE grid (would be LTE-flagged).
- **Ba** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). NLTE grid available (Ba_Korotin2015.csv).
- **Y** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). no NLTE grid (would be LTE-flagged).
- **Li** — CN-blended upper limit (RYA-103); A(Li) 0.73 is a LTE lower bound, not a clean determination. Curation/3D-NLTE owed for a real value.
- **Eu** — solar EW measured + matched in linelist_solar, but the line drops in the GES synthesis-region match of the EW->A(X) path (RYA-395 curate_nonfe_pools not in default run). no NLTE grid (would be LTE-flagged).

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
