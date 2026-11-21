# RYA-1354: calibrating theoretical Fe I gf against the laboratory

Status: **CANDIDATE / HOLD**. Nothing here writes a grade or a sigma into a product. It needs Ryan's ruling under RYA-1352.

Reproduce with `~/codex/venv/bin/python scripts/rya1354_theory_calibration_fe.py`. The outputs go to `data/audit/rya1354/theory_calibration_fe/`:
- `summary.json`
- `k17_vs_lab.csv`
- `candidates.csv`

This is the companion to the O I calibration ([theory_calibration_oi.md](theory_calibration_oi.md)). Fe I is the opposite case. O I has almost no lab gf, so the theory error had to be inferred from 13 lifetimes. Fe I has **493 modern primary lab gf with per-line errors**, so the error of a calculation is **measured line by line**.

## Inputs

- **Theory:** Kurucz's Fe I calculation (`data/reference/kurucz_fe1/gf2600.pos.gz`). It covers 126,288 lines between known levels, uses semi-empirical least-squares eigenvectors, and carries reference code K17 on every row, with no lab values substituted. VALD falls back on this source when no lab gf exists.
- **Lab:** `data/reference/fe_gf_lab/`. These are the primary papers already curated in RYA-799, every row with a per-line sigma:

  | Paper | Lines | Range |
  |---|---|---|
  | Belmonte 2017 | 120 | UV |
  | Den Hartog 2014 | 203 | |
  | Ruffoni 2014 | 142 | |
  | Ruffoni 2013 | 28 | H band |

- **Matching:** both level energies within 0.25 cm⁻¹. All 493 lab lines match exactly one Kurucz line.

## Results

**Overall:** Kurucz minus lab has a median of +0.075 dex and a robust σ of 0.20 dex. Only 30% of lines fall within 0.1 dex, and the tails are heavy (rms 0.60 dex).

**The errors are predictable.** Kurucz publishes no quality indicator; there is no gauge pair and no cancellation factor. Two properties sort his lines:

| K class (spin change \| Kurucz log gf) | n | Bias (K17 − lab) | σ68 about bias | σ95 | Off by > 0.5 dex |
|---|---|---|---|---|---|
| dS0 \| strong (≥ −1) | 183 | +0.126 | **0.109** | 0.335 | 2% |
| jK \| strong (4f, pair coupling) | 13 | +0.035 | 0.161 | 0.443 | 0% |
| dS0 \| mid (−2.5 … −1) | 147 | +0.023 | 0.212 | 0.595 | 7% |
| dS \| strong | 29 | +0.150 | 0.187 | 0.752 | 10% |
| dS \| mid | 86 | −0.060 | 0.350 | 0.881 | 20% |
| dS \| weak (< −2.5) | 25 | −0.954 | 0.842 | 2.9 | 60% |
| dS0 \| weak | 9 | −1.106 | 1.213 | 2.5 | 67% |

- **The worst failures are spin-changing lines that Kurucz predicts to be weak.** For example, 11316 Å y⁵F → e⁷D: Kurucz −8.26, lab −1.70. The real spin mixing is far stronger than his fit.
- **Upper-level siblings do not help.** Correcting a line by the offset of lab lines sharing its upper level leaves the scatter unchanged (0.205 → 0.22 dex). Kurucz's errors are line-specific.

**IR (> 1 µm, H band, Ruffoni 2013):** for the best class (dS0 | strong, n = 13), Kurucz runs **+0.117 dex** too strong. After that bias is removed, the residual is **fully explained by the lab errors** (extra theory scatter 0.008 dex; robust σ 0.064). In the H band, bias-corrected Kurucz is about as good as the lab for strong, spin-conserving lines. Thirteen lines is a small sample, so the class-wide σ68 of 0.11 dex stays the adopted candidate.

## Proposed Fe I theoretical uncertainty (candidate)

A Kurucz gf used where no lab gf exists:
- is **corrected by its class bias**;
- carries **σ = σ68 of its class**: 0.11 dex for dS0 | strong, 0.16 jK | strong, 0.21 dS0 | mid;
- lists σ95 alongside, because the tails are real.

dS | mid and both weak classes are **not** Theoretical Codex material.

Compare O I with Li's MCDHF: 0.064 dex, with independent theories agreeing to 0.01. Modern ab initio O I is about twice as good as semi-empirical Fe I. Kurucz's failure modes, spin-changing and weak lines, are exactly where an ab initio Fe model would have to earn its place.

## Fe I line candidates (our solar line list)

The starting set is Fe I lines that are measurable (solar depth ≥ 0.05), not blend-flagged, have a unique Kurucz match, have no primary lab gf, and are not in the Asplund set. That gives 12,088 lines, split into tiers:

- **Gap:** VALD's own gf *equals* Kurucz's (within 0.006 dex), i.e. VALD had no lab value.
- **T2 (Theoretical Codex candidate):** a gap line in class dS0 | strong or jK | strong, 0.1 ≤ depth ≤ 0.8 (not saturated), observable from the ground (≥ 300 nm).

| Tier | UV < 300 | 300–400 | 400–700 | 700–1000 | > 1 µm |
|---|---|---|---|---|---|
| T0: VALD has a non-Kurucz gf (not a gap) | 2930 | 3196 | 2443 | 318 | 1074 |
| **T2: Theoretical Codex candidate** | — | 1 | 38 | 35 | **271** |
| T3: gap, lower class or depth | — | 255 | 162 | 47 | 98 |
| T4: gap, poor class or space UV | 592 | 341 | 194 | 38 | 55 |

**The Fe I gap is overwhelmingly infrared.** 271 of the 345 candidates lie beyond 1 µm. The deepest are:

| λ (Å) | Depth | Transition |
|---|---|---|
| 15769.4 | 0.50 | e⁵D → u⁵F |
| 16486.7 | 0.49 | e⁵F → 5G |
| 16444.8 | 0.49 | e⁵F → 5F |
| 14400.5 | 0.48 | e⁷D → ⁷F |
| 14512.2 | 0.48 | e⁷D → n⁷F |
| 16102.4 | 0.47 | e⁵F → 5G |

These are 4s5s → 4s5p and (⁴F)5s → (⁴F)5p transitions, the same family as Ruffoni 2013's lab lines, which is where the IR calibration above was measured. NIR examples: 7780.6, 7389.4, 8085.2, 8793.3 and 8764.0 Å.

## Caveats

- Depth and blend flags come from `linelist_solar.csv` (VALD synthesis estimates), not from the observed atlas. A candidate still has to pass the usual line inventory in the spectra.
- "Gap" means VALD carries Kurucz's value. A lab gf VALD does not know about, such as a recent paper, would remove a line from the gap. The literature check runs per line during grading.
- The 493 lab lines are dominated by 400–1000 nm, log gf ≥ −2. Class sigmas for the UV and for weak lines rest on fewer lines.
