# Solar gold reference v2 — TIERED ratification table (RYA-522)

Verdict-sourced (RYA-521), tiered by row-confidence per Ryan's ratification. Asplund = `SOLAR_ASPLUND2021` (Asplund, Amarsi & Grevesse 2021, A&A 653, A141).

**Tiers:** gold=6 · gf_floor=2 · upper_limit=1 · owed(held, no frozen value)=17.

Scales differ (ours 1D-NLTE/synth vs Asplund 3D-NLTE) — `note` states each row's scale so documented offsets are not misread as disagreement. `owed` rows freeze NO value (the C=10.26 lesson: suspect → held, not immortalised).

| Element | conf | v2 (frozen) | method/scale | v1 (old) | Δ(v2−v1) | Asplund 2021 | Δ(v2−Asp) | note |
|---|---|---|---|---|---|---|---|---|
| Fe | gold | 7.516 | 1D-NLTE (Fe I) | 7.516 | +0.000 | 7.46 | +0.056 | our 1D-NLTE runs ~+0.05 above Asplund 3D-true 7.46 (RYA-336) — documented offset, NOT a discrepancy |
| C | gold | 8.491 | synthesis | 10.260 | -1.769 | 8.46 | +0.031 | CH G-band + C I (RYA-237) |
| O | gold | 8.735 | synthesis | — | — | 8.69 | +0.045 | O I 777 + [O I] 6300 (RYA-237) |
| Mg | owed | owed | EW 1D-LTE/NLTE | — | — | 7.55 | — | no independent-gf line survives the graded cull |
| Si | gf_floor | 7.888 | EW 1D-LTE/NLTE | 7.888 | +0.000 | 7.51 | +0.378 | characterized gf-scale floor (+0.38); 3D not the lever (RYA-398/399) |
| Ca | owed | [6.324 held] | EW 1D-LTE/NLTE | 6.324 | — | 6.30 | +0.024 | LOW_CONFIDENCE / thin graded pool (+0.02) |
| Ti | owed | [5.471 held] | EW 1D-LTE/NLTE | 5.471 | — | 4.97 | +0.501 | LOW_CONFIDENCE / thin graded pool (+0.50) |
| Ni | owed | [6.946 held] | EW 1D-LTE/NLTE | 6.946 | — | 6.20 | +0.746 | LOW_CONFIDENCE / thin graded pool (+0.75) |
| Na | owed | [6.264 held] | EW 1D-LTE/NLTE | 6.264 | — | 6.24 | +0.024 | LOW_CONFIDENCE / thin graded pool (+0.02) |
| P | owed | [6.610 held] | atlas 1D | — | — | 5.41 | +1.200 | near-IR multiplet, gf-limited (RYA-460) |
| S | owed | [7.753 held] | EW 1D-LTE/NLTE | 7.753 | — | 7.12 | +0.633 | LOW_CONFIDENCE / thin graded pool (+0.63) |
| N | owed | [8.202 held] | atlas 1D | — | — | 7.83 | +0.372 | N I red multiplets; +0.37 owed NLTE (RYA-369) |
| Co | owed | [6.128 held] | atlas 1D | — | — | 4.94 | +1.188 | blue-edge, SNR-limited — not trusted (RYA-460) |
| Cr | gf_floor | 6.022 | EW 1D-LTE/NLTE | 6.022 | +0.000 | 5.62 | +0.402 | characterized gf-scale floor (+0.40); 3D not the lever (RYA-398/399) |
| Al | owed | [7.406 held] | EW 1D-LTE/NLTE | 7.406 | — | 6.43 | +0.976 | LOW_CONFIDENCE / thin graded pool (+0.98) |
| K | gold | 5.099 | atlas 1D | — | — | 5.07 | +0.029 | K I 7699 + K NLTE grid (RYA-462) |
| Ba | owed | owed | EW 1D-LTE/NLTE | — | — | 2.27 | — | no independent-gf line survives the graded cull |
| Y | owed | owed | EW 1D-LTE/NLTE | — | — | 2.21 | — | no independent-gf line survives the graded cull |
| V | owed | [3.917 held] | synthesis | — | — | 3.90 | +0.017 | HFS-resolved synthesis (RYA-411/466/473) |
| Cu | owed | [4.345 held] | synthesis | — | — | 4.18 | +0.165 | HFS-resolved synthesis (RYA-411/466/473) |
| Mn | gold | 5.470 | synthesis | — | — | 5.42 | +0.050 | HFS-resolved synthesis (RYA-411/466/473) |
| Sc | gold | 3.203 | atlas 1D | — | — | 3.14 | +0.063 | blue-edge HFS single line (RYA-460) |
| Li | upper_limit | 0.727 | EW (upper limit) | 0.727 | +0.000 | 1.05 | -0.323 | CN-blended, carried as UPPER LIMIT (RYA-103) — a clean low value is a red flag |
| Eu | owed | owed | EW 1D-LTE/NLTE | — | — | 0.52 | — | no independent-gf line survives the graded cull |
| Zr | owed | owed | EW 1D-LTE/NLTE | — | — | 2.59 | — | no independent-gf line survives the graded cull |
| Sr | owed | [4.961 held] | EW (SUSPECT) | 4.961 | — | 2.83 | +2.131 | +2.13 — NOT a gf-floor; saturated-line-on-flat-COG signature (the RYA-520 disease) → saturation-trace owed |

**Headline:** C 10.26 → 8.491 (RYA-520 saturated-C I-5380 fix, −1.77 dex).
**+2.1-tail suspect:** Sr (Δ+2.13) — saturated-line-on-flat-COG signature (RYA-520 class) → held owed, routed to a saturation-trace ticket.

