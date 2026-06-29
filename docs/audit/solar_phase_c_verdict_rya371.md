# RYA-371 Phase C — Solar 27-element verdict table (RYA-239 retry)

_Generated 2026-06-29 by scripts/phase_c_verdict_rya371.py. Validate-don't-tune: classification only, no correction fitted to the anchor._

_RYA-460: Kitt Peak Solar Flux Atlas wired in (N + P/K/Co/Sc). Leg validated by overlap: PASS._

## Verdict counts

- **PASS**: 4  (prior 3, diff +1)
- **NLTE-OWED**: 1  (prior 2, diff -1)
- **CURATION-OWED**: 21  (prior 21, diff +0)
- **DATA-GAP**: 0  (prior 0, diff +0)

## RYA-460 overlap cross-check (Kitt Peak leg validation)

| line | Kitt Peak raw | Phase A (arm) | Δ | agree |
|------|--------------:|---------------|---:|:-----:|
| OI_6300 | 8.835 | 8.8 (harps) | +0.035 | YES |
| OI_777 | 8.955 | 8.917 (espresso) | +0.038 | YES |

**N solar cross-indicator map** (NLTE-OWED, not validated): N I red mean **8.202** (spread 0.033) from {'NI_7442_7468': 8.189, 'NI_8216_8223': 8.222, 'NI_8680_8718': 8.196, 'NH_3360': 7.197, 'CN_violet_3883': 9.018}. NH 3360 / CN violet blue-edge FLAGGED (['NH_3360', 'CN_violet_3883']).

## Per-element table

| El | Asplund21 | A(meas) | Delta | sigma | n | NLTE | Verdict | Provenance | Channel |
|----|----------:|--------:|------:|------:|--:|:----:|:--------|:-----------|:--------|
| O | 8.69 | 8.735 | +0.045 | 0.01 | 3 | - | **PASS** | synthesis: harps/espresso (Phase A) | synthesis: O I 777 (primary) + [O I] 6300 (cross-check) |
| C | 8.46 | 8.491 | +0.031 | 0.05 | 5 | - | **PASS** | synthesis: harps/espresso (Phase A) | synthesis: CH G-band + C I 5052 + C2 Swan (C I 5380 BAD_FIT-excluded) |
| N | 7.83 | 8.202 | +0.372 | 0.03 | 3 | - | **NLTE-OWED** | kittpeak-measured | kittpeak: N I red 7442/7468 + 8216/8223 + 8680-8718 (NH/CN blue-edge flagged) |
| Mg | 7.55 |  |  |  |  | wired | **CURATION-OWED** | harps-measured (EW pool) | EW present; no independent-gf line survives the graded cull |
| Si | 7.51 | 7.888 | +0.378 | 0.36 | 7 | wired+3D | **CURATION-OWED** | harps-measured (EW pool) | EW: 7 curated line(s), graded-gf (RYA-395/398) |
| Fe | 7.46 | 7.516 | +0.056 | 0.14 | 62 | - | **PASS** | harps-measured (EW) | EW: 62 Fe I + 3 Fe II, NLTE-wired (Bergemann MPIA) |
| S | 7.12 | 7.753 | +0.633 | 0.37 | 2 | wired | **CURATION-OWED** | harps-measured (EW pool) | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Al | 6.43 | 7.406 | +0.976 |  | 1 | wired | **CURATION-OWED** | harps-measured (EW pool) | EW: 1 curated line(s), low-confidence (RYA-395/398) |
| Ca | 6.30 | 6.324 | +0.024 | 0.12 | 2 | wired | **CURATION-OWED** | harps-measured (EW pool) | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Na | 6.24 | 6.264 | +0.024 | 0.02 | 2 | wired | **CURATION-OWED** | harps-measured (EW pool) | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Ni | 6.20 | 6.946 | +0.746 | 0.51 | 2 | - | **CURATION-OWED** | harps-measured (EW pool) | EW: 2 curated line(s), low-confidence (RYA-395/398) |
| Cr | 5.62 | 6.022 | +0.402 | 0.60 | 7 | wired+3D | **CURATION-OWED** | harps-measured (EW pool) | EW: 7 curated line(s), graded-gf (RYA-395/398) |
| Mn | 5.42 | 5.554 | +0.134 | 0.15 | 3 | wired | **CURATION-OWED** | synthesis: HFS-resolved (RYA-473) | HFS synthesis: Mn I 3 lines (6013/6016/6021, Den Hartog e6S→z6P), gf=Den Hartog+2011 (MED), NLTE vendored MPIA/Bergemann grid (live Amarsi .grd offline) |
| P | 5.41 | 6.610 | +1.200 |  | 2 | - | **CURATION-OWED** | kittpeak-measured | kittpeak: P I 10581/10596 near-IR multiplet |
| K | 5.07 | 5.099 | +0.029 |  | 1 | wired | **PASS** | kittpeak-measured | kittpeak: K I 7699 (clean; 7665 in the telluric O2 A-band) — NLTE-wired |
| Ti | 4.97 | 5.471 | +0.501 | 0.93 | 10 | wired+3D | **CURATION-OWED** | harps-measured (EW pool) | EW: 10 curated line(s), graded-gf (RYA-395/398) |
| Co | 4.94 | 6.128 | +1.188 |  | 1 | - | **CURATION-OWED** | kittpeak-measured | kittpeak: Co I 3845 (blue-edge, SNR-limited) |
| Cu | 4.18 | 4.345 | +0.165 | 0.14 | 5 | - | **CURATION-OWED** | synthesis: HFS-resolved (RYA-466) | HFS synthesis: Cu I 5 lines (5105/5218/5220/5700/5782), gf=Kock&Richter, NLTE vendored RYA-402 b-factor (live .grd offline) |
| V | 3.90 | 3.917 | +0.017 | 0.03 | 6 | - | **CURATION-OWED** | synthesis: HFS-resolved LTE (RYA-466) | HFS synthesis LTE: V I 6 lines — NLTE-VOID (no model atom) |
| Sc | 3.14 | 3.203 | +0.063 |  | 1 | - | **CURATION-OWED** | kittpeak-measured | kittpeak: Sc II 4246 (blue-edge, HFS) |
| Sr | 2.83 | 4.961 | +2.131 |  | 1 | wired | **CURATION-OWED** | harps-measured (EW pool) | EW: 1 curated line(s), low-confidence (RYA-395/398) |
| Zr | 2.59 |  |  |  |  | - | **CURATION-OWED** | harps-measured (EW pool) | EW present; no independent-gf line survives the graded cull |
| Ba | 2.27 |  |  |  |  | wired | **CURATION-OWED** | harps-measured (EW pool) | EW present; no independent-gf line survives the graded cull |
| Y | 2.21 |  |  |  |  | - | **CURATION-OWED** | harps-measured (EW pool) | EW present; no independent-gf line survives the graded cull |
| Li | 1.05 | 0.727 | -0.323 |  | 1 | - | **CURATION-OWED** | harps-measured (EW pool) | EW: Li I 6707 (single line, UPPER LIMIT) |
| Eu | 0.52 |  |  |  |  | - | **CURATION-OWED** | harps-measured (EW pool) | EW present; no independent-gf line survives the graded cull |

## Remaining-work map

### PASS (4)

- **O** — cross-arm AGREE; O I 777 Amarsi-2019 3D-NLTE, [O I] Caffau-2015 3D anchor — measured 8.74 vs Asplund 8.69 (+0.05). RYA-455.
- **C** — C I Amarsi-2019 3D-NLTE -> 8.46; CH 8.49 (3D-offset-owed); C I 5380 formally excluded ew_integrity=BAD_FIT (RYA-458); surviving cross-arm spread 0.054 on 5 indicators. [C I 5380 EXCLUDED (ew_integrity=BAD_FIT, espresso:CI_5380); cross-arm spread 0.149->0.054 on 5 surviving indicators.]
- **Fe** — A(Fe I) NLTE 7.516 vs Asplund 7.46 (+0.056); ionization-balance gated, scatter 0.139 = honest floor (RYA-407). Documented +0.05 1D/3D scale offset (RYA-336), not the verdict.
- **K** — MEASURED Kitt Peak K I 7699 = 5.411 (1D-LTE, +0.34 vs 5.07). K_Amarsi2020_PySME NLTE delta -0.312 APPLIED via the existing interpolation subsystem (RYA-462 wiring; NLTE_Amarsi2020_PySME_1D; validate-don't-tune) -> A(K) 5.099 (+0.029 vs Asplund 5.07). Reconciles within TOL after the cited NLTE correction — the severe negative K I resonance NLTE is real, not tuned.

### NLTE-OWED (1)

- **N** — MEASURED from Kitt Peak N I red — 3 independent multiplets AGREE: 8.189 / 8.222 / 8.196 (mean 8.202, spread 0.033). +0.37 vs Asplund 7.83 is the N I NLTE offset OWED (N I grid RYA-369; NLTE is negative, pulls toward 7.83). NOT validated: Teff-bracket owed (Procyon / aCen B, RYA-369). NH 3360 + CN violet 3883 UNMEASURABLE here — blue-edge no-true-continuum (SNR~28, RYA-451/454) + the Turbospectrum molecular linelist is absent — FLAGGED, not forced. Kitt Peak leg VALIDATED by the [O I]6300/O I 777 overlap cross-check vs HARPS/ESPRESSO (agree within 0.04).

### CURATION-OWED (21)

- **Mg** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Mg_Amarsi2020_PySME.csv).
- **Si** — A(X) 7.888 vs Asplund 7.51 (+0.378), sigma 0.36 — gross offset removed by the blind cull, gf-scale residual survives on the graded pool → escalate to RYA-161/162 (do NOT tune).
- **S** — A(X) 7.753 vs Asplund 7.12 (+0.633) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Al** — A(X) 7.406 vs Asplund 6.43 (+0.976) on 1 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Ca** — A(X) 6.324 vs Asplund 6.30 (+0.024) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Na** — A(X) 6.264 vs Asplund 6.24 (+0.024) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Ni** — A(X) 6.946 vs Asplund 6.20 (+0.746) on 2 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Cr** — A(X) 6.022 vs Asplund 5.62 (+0.402), sigma 0.60 — gross offset removed by the blind cull, gf-scale residual survives on the graded pool → escalate to RYA-161/162 (do NOT tune).
- **Mn** — MEASURED via HFS-resolved synthesis — the EW path SAT-culls the Den Hartog triplet (REW −4.78..−4.82 over the −4.90 knee, HFS-split hfs_n=6; RYA-468 finding: gf graded but saturation is the blocker). Synthesis on the GES HFS line list (6 components/feature = Den Hartog Table 4, cited) measures it. A(Mn)_LTE 5.446 + Mn NLTE +0.108 (vendored MPIA/Bergemann grid (live Amarsi .grd offline)) = 5.554 (+0.134 vs Asplund 5.42; σ 0.154, n=3). The NLTE δ is the MPIA grid's HIGH-EP reference value (+0.107); the low-EP triplet is NOT a grid node and its HFS-resolved Amarsi δ differs (RYA-411, .grd offline) → do NOT certify PASS on the vendored δ. The MEASUREMENT-TOOL blocker is fixed (off no-value); the line-exact NLTE is owed (RYA-411 Amarsi grid).
- **P** — MEASURED from Kitt Peak P I near-IR = 6.61 (+1.20 vs 5.41) — OFF DATA-GAP: the near-IR multiplet is reachable from the ground, no HST/STIS needed (RYA-119 superseded for the Sun). The large +1.2 offset is a gf-scale residual (P I near-IR gf are uncertain) → curation owed RYA-161/162; do NOT tune.
- **Ti** — A(X) 5.471 vs Asplund 4.97 (+0.501), sigma 0.93 — gross offset removed by the blind cull, gf-scale residual survives on the graded pool → escalate to RYA-161/162 (do NOT tune).
- **Co** — Kitt Peak covers Co, but the extracted Co I 3845 sits in the blanketed blue edge (SNR~24, chi2r~3100) → the value 6.128 is NOT trusted (blue-edge per the RYA-451/454 caveat). OFF pure DATA-GAP (a measured reference now exists) but curation owed: extract cleaner red Co I lines (within KP's 1300 nm reach) + HFS. Do NOT force the blue value.
- **Cu** — MEASURED via HFS-resolved synthesis — the EW path could not reach Cu (all 5 lines hyperfine-split, RYA-354 finding); synthesis on the GES HFS line list (4-10 components/feature) measures it. gf adjudicated to GES=Kock&Richter (KR), VALD3 superseded (GES = Kock & Richter (reference_code KR; the NIST-graded authority RYA-354 named) — VALD3 superseded, median |Δ(GES−VALD)|=0.327 dex). A(Cu)_LTE 4.342 + RYA-402 b-factor NLTE +0.003 (vendored RYA-402 b-factor (live .grd offline)) = 4.345 (+0.165 vs Asplund 4.18; σ 0.139, n=5). The +offset survives the small Cu NLTE → 1D-LTE optical Cu I sits high; 3D / fuller-NLTE curation owed (RYA-161/162). The blocker was the MEASUREMENT TOOL, now fixed; the residual is a finding, do NOT tune.
- **V** — MEASURED via HFS-resolved LTE synthesis — V was unmeasurable by the EW path (HFS, RYA-354). A(V)_LTE 3.917 (+0.017 vs Asplund 3.90; σ 0.029, n=6) — tight, AGREES with Asplund, BUT V is the NLTE-VOID (no model atom anywhere, RYA-463) → LTE-only, LOWER CONFIDENCE; a V NLTE model + graded gf (V I Lawler+2014, V II Wood+2014) owed before any PASS. Off no-value, do NOT certify on LTE agreement alone.
- **Sc** — MEASURED from Kitt Peak Sc II 4246 = 3.203 (+0.06 vs 3.14) — OFF DATA-GAP, value close to Asplund BUT single blue-edge HFS line (SNR~180, no true continuum) → LOW_CONFIDENCE; HFS-resolved synthesis + a cleaner Sc II line owed before any PASS.
- **Sr** — A(X) 4.961 vs Asplund 2.83 (+2.131) on 1 graded line(s) — below the stable-mean floor; thin independent-gf pool, differential-survey curation owed (RYA-161/162).
- **Zr** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged).
- **Ba** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). NLTE grid available (Ba_Korotin2015.csv).
- **Y** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged).
- **Li** — CN-blended UPPER LIMIT (RYA-103/458, ew_integrity disposition=UPPER_LIMIT); A(Li) 0.73 is a LTE lower bound, not a clean determination. A clean low value here would be a RED FLAG (CN deblend not applied).
- **Eu** — solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf firewall (now wired into the default run, RYA-456) culls every line — the pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). no NLTE grid (would be LTE-flagged). Eu II 6645 EW 6.8 mA, ew_integrity disposition=RECOVERED (RYA-102/458 HFS-summing).

## Off-Sun — DEFERRED to RYA-348

A green solar run validates the MACHINERY, not target-transferability. These cannot be validated on the Sun and are explicitly deferred:

- per-star broadening sourcing (vmac/vsini)
- temperature-dependent NLTE (e.g. C I LTE-on-Sun but NLTE-on-Procyon)
- F-star / cool-star line-list adequacy

_Not "55 Cnc-ready" — that is RYA-348 and beyond._
