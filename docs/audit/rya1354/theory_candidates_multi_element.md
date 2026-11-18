# RYA-1354: Theoretical Codex candidates for 15 species, checked against the observed Sun

Status: **CANDIDATE / SCREENING**. Nothing here writes a grade or a sigma into a product. It needs Ryan's ruling under RYA-1352.

This extends [theory_calibration_fe.md](theory_calibration_fe.md) and [theory_calibration_oi.md](theory_calibration_oi.md).

## 1. Kurucz vs laboratory, 15 species

Reproduce with `scripts/rya1354_theory_calibration_kurucz.py`; outputs go to `data/audit/rya1354/theory_calibration_kurucz/`.

- **Theory:** Kurucz's calculations for Si I, Ca I, Ti I/II, V I/II, Cr I/II, Mn I/II, Fe I/II, Co I/II and Ni I (`data/reference/kurucz/`).
- **Lab:** 3,458 primary lab gf, matched by level identity (≤ 13 unmatched per species).
- **Quality class:** spin change × Kurucz strength, as for Fe I.

| Species | Lab lines | Robust σ (all) | Best class (dS0 \| strong): bias / σ68 |
|---|---|---|---|
| Ca I | 77 | 0.13 | +0.04 / 0.08 |
| Ti II | 364 | 0.09 | +0.01 / 0.05 |
| V I | 836 | 0.20 | +0.09 / 0.14 |
| V II | 203 | 0.06 | −0.01 / 0.07 |
| Cr I | 263 | 0.10 | +0.04 / 0.11 |
| Cr II | 183 | 0.10 | −0.01 / 0.04 |
| Mn I | 62 | 0.07 | −0.02 / 0.08 |
| Fe I | 493 | 0.20 | +0.13 / 0.11 |
| Fe II | 131 | 0.06 | +0.04 / 0.04 |
| Co I | 898 | 0.27 | +0.11 / 0.17 |
| Co II | 84 | 0.08 | +0.08 / 0.08 |
| Ni I | 371 | 0.21 | +0.05 / 0.10 |

Pooled over all species:

| Class | σ68 (dex) |
|---|---|
| dS0 \| strong | 0.117 |
| dS0 \| mid | 0.20 |
| jK \| strong | 0.15 |
| dS \| mid | 0.30 |
| dS \| weak | 0.57 |
| dS0 \| weak | 0.74 |

- **Ions are the most accurate:** Kurucz's best class reaches 0.04–0.07 dex for Ti II, V II, Cr II and Fe II.
- **Si I has no usable local lab calibration** (2 lines). Its candidates use the pooled class σ and are tiered T2P, not T2.

**Candidate rule (gap):** VALD's gf equals Kurucz's to within 0.006 dex, meaning no lab value exists. On top of that:
- class σ68 ≤ 0.15 dex;
- solar depth 0.1–0.8 in the VALD estimate;
- not blend-flagged;
- ≥ 300 nm.

This gives 1,243 species-calibrated candidates plus 349 with only the pooled σ.

## 2. Check against the observed solar spectrum

Reproduce with `scripts/rya1354_candidate_spectrum_check.py`; output goes to `data/audit/rya1354/candidate_spectrum_check/checked.csv`.

The O I T1/T2 lines are included. Components within 6 km/s of each other are checked as one feature (multiplets). Each feature is looked at in a normalised, telluric-corrected or telluric-free solar atlas:

| Range | Atlas |
|---|---|
| < 5001 Å | Kitt Peak 1984 |
| 5001–11083 Å | IAG Baker 2020 (telluric-free) |
| 0.98–1.08 µm | CRIRES+ Vesta Y |
| 1.12–1.35 µm | CRIRES+ Vesta J (molecfit-corrected) |
| 1.49–1.76 µm | CRIRES+ Vesta H (molecfit-corrected) |
| 1.95–2.48 µm | CRIRES+ Vesta K (molecfit-corrected) |

A feature is **CLEAN** when:
- its observed depth is ≥ 5%;
- the minimum lies within 2 km/s of the catalogue wavelength;
- the core asymmetry is ≤ 0.25;
- no neighbour deeper than half the line sits within 2–4 HWHM;
- molecfit transmission is ≥ 0.9 in the J/H/K core.

**Positive control:** the 40 Asplund 2021 Fe I lines plus O I 777 pass the offset, asymmetry and telluric tests 43/43; "close neighbour" flags 3/40.

Observed/predicted depth is reported but is **not** a test. The predicted depths are VALD LTE estimates with uncorrected theory gf:
- O I 777 is observed 1.8× deeper than predicted (NLTE).
- H-band Fe I is observed at a median 0.64 of the predicted depth, consistent with Kurucz's gf running too strong.

**Result:** of 1,566 features, 278 are CLEAN, 1,143 FLAGGED and 145 not covered (the 1.35–1.49 and 1.76–1.95 µm water bands, plus the O I 1.8 µm multiplets).
- **New:** 263 of the clean features are new, i.e. in no Codex line set.
- **Already in a set:** the other 15 are mainly Si I IR lines Elgueta 2026 or Bergemann 2013 already use.
- **Near-UV:** almost everything fails there (about 8% of KP-atlas features clean), because the region is crowded.

**Clean, new, ≥ 700 nm (132 features):**

| Species | n | Deepest (Å, observed depth) |
|---|---|---|
| Fe I | 52 | 7780.6 (0.55), 8085.2 (0.53), 8793.3, 8764.0, 8028.3, **16486.7 (0.40)**, 8439.6, **16444.8**, 7244.8, 10145.6 |
| Si I | 25 | 10827.1 (0.68), 10585.1, **15888.4 (0.55)**, 10603.4, 10979.3, 10843.9, **15960.1**, **16680.8**, 13176.9, 7932.3 (pooled σ) |
| Ni I | 15 | 7555.6 (0.52), 7797.6, 7574.0, 7182.0, 8862.5, **16363.1**, **17001.0**, 8965.9, **15555.4**, 10302.6 |
| Cr I | 13 | 9009.8 (0.34), 9021.6, 8947.2, 11015.5, 10905.7, 11398.1, 10486.3 |
| Ti I | 13 | 9675.5 (0.39), 9705.7, 9728.4, 9787.7, 9599.6, 9647.4, 8396.9 |
| Ca I | 8 | **19853.1 (0.37)**, **22651.2**, **22625.0**, **16197.1**, 13134.9, 13033.6 (K and H band) |
| O I | 3 | 777.19/777.42/777.54 nm (Li class B, σ 0.064 dex) |
| Ti II | 3 | 7214.7, 8979.2, 15873.8 |

Every Kurucz candidate carries a bias-corrected log gf and the class σ68 and σ95 (`checked.csv`).

## 3. Caveats and next steps

- **Screening, not measurement.** The geometric tests catch obvious blends and misidentifications. Each line still needs the usual line inventory and synthesis fit before it is graded.
- **Gap test:** "VALD = Kurucz" is the gap test. A recent lab paper VALD does not include would remove a line, so each element's literature check runs before grading.
- **Si I:** it needs its own lab or independent-theory calibration before its candidates can be T2. Deshmukh 2022, Bergemann 2013 and Elgueta 2026 line sets already exist in `data/reference/line_sets/`.
- **Natural first test:** with Ryan's ruling, the obvious first run is to measure a handful of these with the ORCH using the theory gf ± σ. Good starters are the Fe I H-band set (where the H-band lab calibration says bias-corrected Kurucz is lab-quality) and Ca I K band. Checking the abundance they return against each element's lab-graded result tests the whole framework end to end.
