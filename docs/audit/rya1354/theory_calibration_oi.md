# RYA-1354: calibrating theoretical O I gf against the laboratory

Status: **CANDIDATE / HOLD**. Nothing here writes a grade or a sigma into a product. The class thresholds and the uncertainty rule need Ryan's ruling under RYA-1352 (Theoretical Codex governance).

Reproduce with `~/codex/venv/bin/python scripts/rya1354_theory_calibration_oi.py`. The outputs go to `data/audit/rya1354/theory_calibration_oi/`:
- `summary.json`
- `direct_gf.csv`
- `li_vs_nist.csv`
- `candidates.csv`

## The question

Many UV and IR lines have no laboratory gf. Can a calculated gf stand in for one, and how uncertain is it?

To find out, I took the best published O I calculation and checked it three ways:
- against every laboratory measurement that exists;
- against an independent calculation;
- then I ranked the O I lines in our solar line list as Theoretical Codex candidates.

The calculation is Li et al. 2023: GRASP MCDHF/RCI, 81 states, 989 lines. It is the scheme the RYA-1354 pilot copies.

## First finding: O I has almost no laboratory gf

NIST's O I values (Wiese, Fuhr & Deters 1996) are **themselves calculated**, mostly Hibbert et al. 1991. That includes the 777, 844.6, 926 and 615.8 nm lines Asplund uses. The prior curation in `source_curation_v2/evaluated_4d_validation.json` reached the same conclusion: zero experimental gf rows.

Every experimental anchor Li could find:
- **Direct gf:** 12 lines, all vacuum UV (98.9–130.6 nm), with absorption, emission and beam-foil data.
- **Lifetimes:** 13 upper levels. A lifetime fixes the *sum* of the A-values out of a level. For the 3p ³P level that sum is entirely the **844.6 nm** triplet; for 4d ⁵D° it is mostly **615.8 nm** plus 2.65 µm lines.

There is **no** experimental gf or lifetime for the 777 nm triplet's upper level, nor for any 926 nm, 1.13 µm, 1.32 µm or 1.8 µm line.

## Test 1: direct gf (UV)

| | |
|---|---|
| Lines | 12, all Li class A |
| Li minus experiment | mean −0.009 dex, rms 0.046 dex |
| χ² using experimental errors only | 12.1 for 12 lines |
| Extra theory scatter needed | **0** |

Li agrees with the UV measurements within the experiments' own precision (0.01–0.07 dex). The experiments cannot see a theory error smaller than about 0.03 dex.

## Test 2: lifetimes (reaches VIS/IR lines)

| Level | Lines it tests | Li τ | Exp. τ (weighted) | Δlog A (Li − exp) | z |
|---|---|---|---|---|---|
| 3s ³S° | 130.2/130.5/130.6 nm | 1.70 ns | 1.80 ± 0.04 | +0.025 | +2.4 |
| 3s ⁵S° | 135.6/135.9 nm | 209 µs | 181 ± 4 | −0.064 | −6.5 |
| **3p ³P** | **844.6 nm** (100%) | 32.1 ns | 39.0 ± 1.2 | **+0.084** | **+5.7** |
| 4p ⁵P | 2.76 µm, 394.9 nm | 200 ns | 193 ± 9 | −0.015 | −0.8 |
| 4p ³P | 2.89 µm, 436.9 nm | 182 ns | 155 ± 9 | −0.071 | −3.1 |
| **4d ⁵D°** | **615.8 nm** (54%), 2.65 µm | 72.2 ns | 95.8 ± 3.7 | **+0.123** | **+6.5** |

The other seven levels are in `summary.json`. Their experiments disagree with each other (Birge ratio up to 4), or they are single measurements.

Across the levels whose experiments agree, the lifetimes need **0.063 dex of extra scatter** beyond the experimental errors. Theory and experiment differ at the 10–30% level for some VIS/IR upper levels:
- The 844.6 nm level is the clearest case. Its main measurement is two-photon laser-induced fluorescence (Bischel et al.), 39.1 ± 1.4 ns. It says Li, and Hibbert/NIST, which agree with Li, make 844.6 nm about **0.08 dex too strong**.
- The 615.8 nm upper level says about 0.12 dex too strong. The source there is older beam-foil work, where cascades are a known systematic.

These are quotes via Li's table. The primary papers have not been re-read, so this is a lead, not a correction.

## Test 3: independent theory (Li vs NIST/Hibbert), 347 common lines

Per-method scatter is the Li−NIST robust σ divided by √2, assuming the two methods are of equal quality.

| Li class | n | Median Li−NIST | Per-method scatter | Li's own claim |
|---|---|---|---|---|
| A | 61 | −0.008 | 0.009 dex | ≤ 0.013 |
| B | 75 | −0.011 | 0.008 | ≤ 0.041 |
| C | 81 | −0.019 | 0.036 | ≤ 0.097 |
| D | 25 | +0.086 | 0.131 | ≤ 0.176 |
| E | 27 | +0.370 | 0.250 | > 0.18 |

Li's accuracy classes rank lines correctly: independent theories agree to about 0.01 dex for classes A/B and diverge rapidly for D/E. Strong lines (A ≥ 10⁶ s⁻¹) agree at 0.009 dex, while A < 10³ s⁻¹ lines scatter by 0.56 dex. By band, the scatter is 0.012 (IR), 0.013 (VIS) and 0.024 (UV).

**2p³(⁴S)nf lines (1.80/1.82 µm, 1.25 µm).**
- **Components:** the ⁵F and ³F upper levels are 0.094 cm⁻¹ apart and strongly mixed. NIST/VALD (pure LS coupling) and Li (intermediate coupling) put the strength into different components, differing by up to 0.9 dex per feature.
- **Multiplet totals agree** to 0.013 dex (3d ⁵D → 4f) and 0.015 dex (3d ³D → 4f).
- **In the Sun** each multiplet spans ≤ 0.12 nm, comparable to the line width, so these are graded as **one blended feature at the multiplet level**. VALD's per-component split for these lines is coupling-dependent and should not be trusted component by component.

## Proposed theoretical uncertainty (candidate)

σ_theory(class) = √(floor² + per-method scatter²)

- **floor:** 0.063 dex, from the lifetime test. It is the only experiment that reaches VIS/IR upper levels.
- **per-method scatter:** from Test 3.

| Class | Proposed σ | ≈ % in gf |
|---|---|---|
| A | 0.064 dex | 16% |
| B | 0.064 dex | 16% |
| C | 0.073 dex | 18% |
| D | 0.145 dex | 40% |
| E | 0.258 dex | 81% |

For context, Asplund 2021 does not include atomic data in its O budget; Bautista et al. 2022 quote 0.021 dex for 777 from multi-method spread alone. The 0.063 dex floor is driven by three lifetimes (3s ⁵S°, 3p ³P, 4d ⁵D°). If a primary re-read finds systematics in those experiments, the floor drops toward the 0.01–0.03 dex that the theory-to-theory spread and the UV gf allow.

## Candidate lines (our solar line list, 103 O I rows)

| Tier | n | Meaning |
|---|---|---|
| T1 experimental anchor | 18 | direct gf, or an upper level with a lifetime and ≥10% branching |
| T2 Theoretical Codex candidate | 30 | no lab anchor; Li class A/B (or nf multiplet); CF ≥ 0.1; solar depth ≥ 1% |
| T3 | 1 | measurable but low theory class (394.7 nm, class D) |
| T4 | 46 | depth < 1% in the Sun |
| not in Li's set | 8 | [O I] 630.0/636.4 (forbidden, not E1), 121.8 nm, 477/497 nm, 2.11–2.13 µm |

**Calibration anchors:**
- **UV:** 115.2 and 130.2/130.5/130.6 nm, plus 135.6/135.9 nm (class C, intercombination). These need space data.
- **IR:** 844.6 nm (lifetime-anchored).
- **VIS:** 615.8 nm and 436.8 nm, both weak in the Sun.

**Theoretical Codex candidates (ground-observable, IR):**
- 777.2/777.4/777.5 nm (class B)
- 926.1–926.6 nm (A)
- 1128.6–1128.7 nm (A)
- 1129.5–1130.2 nm (B)
- 1316.4/1316.5 nm (A)
- **1589.3 nm** (B; H band)
- the **1.80 and 1.82 µm 3d–4f multiplets**: depth 2–7%, CF 0.5–1.0, graded as multiplets.

## Our own model on the same yardstick

The RYA-1354 FAC models (sealed v41) computed one quantity end to end: the 3s ³S° lifetime, the upper level of 130.2/130.5/130.6 nm. The scores are in `summary.json → own_model_v41`.

| Model | τ (Babushkin) | τ (Coulomb) | Gauge dT | Class on Li's scale | Δlog A vs experiment |
|---|---|---|---|---|---|
| Ours, ground-weighted | 1.50 ns | 1.88 ns | 20% | **C** | +0.080 dex |
| Ours, equal-weighted | 1.51 ns | 3.06 ns | 51% | E | +0.077 dex |
| Li 2023 | 1.70 | 1.69 | 0.6% | A | +0.025 dex |

Experiment: 1.80 ± 0.04 ns.

Our current model is class C: about 0.07 dex on the proposed scale, and 3–4× worse than Li on the cleanest O I line. Before it can fill a gap it has to reach class A/B on Li's lines. That target is concrete: gauge agreement within 10%, matching Li's 989 lines and the 13 lifetimes.

## What this sets up

1. **Our own model:** the RYA-1354 pilot can now be scored against the same three tests. That is the "theoretical lab" check before it is trusted for gaps.
2. **Gaps:** lines outside Li's 81 states have no published high-quality calculation. Our model would have to fill these:
   - 2.11–2.13 µm, K band;
   - higher Rydberg series;
   - 2p³(²D,²P) parents, e.g. 121.8 nm, 477/497 nm.
3. **Fe:** the same framework applies with a much larger lab sample (`data/reference/fe_gf_lab/`). Fe I theory vs lab gives a direct, line-by-line calibration of the floor, which O I cannot provide.
