# Solar gold reference v2 — ratification diff table (RYA-522)

Source: the phase_c **verdict** channel (RYA-521), regenerated — not hand-edited. Asplund column = `SOLAR_ASPLUND2021` (Asplund, Amarsi & Grevesse 2021, A&A 653, A141).

**Scales differ:** our values are 1D-NLTE / synthesis on our stack; Asplund 2021 is 3D-NLTE photospheric. The `note` states each row's scale so documented offsets (e.g. Fe I +0.05, RYA-336) are not misread as disagreement.

| Element | v2 (verdict) | method/scale | v1 (old) | Δ(v2−v1) | Asplund 2021 | Δ(v2−Asp) | verdict | note |
|---|---|---|---|---|---|---|---|---|
| Fe | 7.516 | 1D-NLTE (Fe I) | 7.516 | +0.000 | 7.46 | +0.056 | PASS | our 1D-NLTE runs ~+0.05 above Asplund 3D-true 7.46 (RYA-336) — documented offset, NOT a discrepancy |
| C | 8.491 | synthesis | 10.260 | -1.769 | 8.46 | +0.031 | PASS | CH G-band + C I (RYA-237) |
| O | 8.735 | synthesis | — | — | 8.69 | +0.045 | PASS | O I 777 + [O I] 6300 (RYA-237) |
| Mg | owed | EW 1D-LTE/NLTE | — | — | 7.55 | — | CURATION-OWED | curated EW; low-confidence |
| Si | 7.888 | EW 1D-LTE/NLTE | 7.888 | +0.000 | 7.51 | +0.378 | CURATION-OWED | gf-limited residual floor (+0.38); NOT an Asplund disagreement (RYA-399) |
| Ca | 6.324 | EW 1D-LTE/NLTE | 6.324 | +0.000 | 6.30 | +0.024 | CURATION-OWED | curated EW; low-confidence |
| Ti | 5.471 | EW 1D-LTE/NLTE | 5.471 | +0.000 | 4.97 | +0.501 | CURATION-OWED | gf-limited residual floor (+0.50); NOT an Asplund disagreement (RYA-399) |
| Ni | 6.946 | EW 1D-LTE/NLTE | 6.946 | +0.000 | 6.20 | +0.746 | CURATION-OWED | gf-limited residual floor (+0.75); NOT an Asplund disagreement (RYA-399) |
| Na | 6.264 | EW 1D-LTE/NLTE | 6.264 | +0.000 | 6.24 | +0.024 | CURATION-OWED | curated EW; low-confidence |
| P | 6.610 | atlas 1D | — | — | 5.41 | +1.200 | CURATION-OWED | near-IR multiplet, gf-limited (RYA-460) |
| S | 7.753 | EW 1D-LTE/NLTE | 7.753 | +0.000 | 7.12 | +0.633 | CURATION-OWED | gf-limited residual floor (+0.63); NOT an Asplund disagreement (RYA-399) |
| N | 8.202 | atlas 1D | — | — | 7.83 | +0.372 | NLTE-OWED | N I red multiplets; +0.37 = owed NLTE (RYA-369) |
| Co | 6.128 | atlas 1D | — | — | 4.94 | +1.188 | CURATION-OWED | blue-edge, SNR-limited — value not fully trusted (RYA-460) |
| Cr | 6.022 | EW 1D-LTE/NLTE | 6.022 | +0.000 | 5.62 | +0.402 | CURATION-OWED | gf-limited residual floor (+0.40); NOT an Asplund disagreement (RYA-399) |
| Al | 7.406 | EW 1D-LTE/NLTE | 7.406 | +0.000 | 6.43 | +0.976 | CURATION-OWED | gf-limited residual floor (+0.98); NOT an Asplund disagreement (RYA-399) |
| K | 5.099 | atlas 1D | — | — | 5.07 | +0.029 | PASS | K I 7699 + K NLTE grid (RYA-462) |
| Ba | owed | EW 1D-LTE/NLTE | — | — | 2.27 | — | CURATION-OWED | curated EW; low-confidence |
| Y | owed | EW 1D-LTE/NLTE | — | — | 2.21 | — | CURATION-OWED | curated EW; low-confidence |
| V | 3.917 | synthesis | — | — | 3.90 | +0.017 | CURATION-OWED | HFS-resolved synthesis (RYA-411/466/473) |
| Cu | 4.345 | synthesis | — | — | 4.18 | +0.165 | CURATION-OWED | HFS-resolved synthesis (RYA-411/466/473) |
| Mn | 5.470 | synthesis | — | — | 5.42 | +0.050 | PASS | HFS-resolved synthesis (RYA-411/466/473) |
| Sc | 3.203 | atlas 1D | — | — | 3.14 | +0.063 | CURATION-OWED | blue-edge HFS single line (RYA-460) |
| Li | 0.727 | EW (upper limit) | 0.727 | +0.000 | 1.05 | -0.323 | CURATION-OWED | CN-blended, carried as UPPER LIMIT (RYA-103) — a clean low value is a red flag |
| Eu | owed | EW 1D-LTE/NLTE | — | — | 0.52 | — | CURATION-OWED | curated EW; low-confidence |
| Zr | owed | EW 1D-LTE/NLTE | — | — | 2.59 | — | CURATION-OWED | curated EW; low-confidence |
| Sr | 4.961 | EW 1D-LTE/NLTE | 4.961 | +0.000 | 2.83 | +2.131 | CURATION-OWED | gf-limited residual floor (+2.13); NOT an Asplund disagreement (RYA-399) |

**Headline:** C 10.26 → 8.491 (the RYA-520 saturated-C I-5380 artifact corrected; −1.77 dex).

Freeze is GATED on Ryan's explicit ratification of this table (RYA-522 step 4).
