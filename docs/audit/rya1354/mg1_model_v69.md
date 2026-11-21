# RYA-1354 v69: the Sun as a referee, NIST composition, convergence rule, and a blind isoelectronic test (Al II, Si III)

Status: **CANDIDATE / HOLD** (RYA-1352). Declared on RYA-1354 (2026-11-20 14:55 UTC) before any run. No published data touched.

Outputs:
- `data/audit/rya1354/mg1_model_v69/`
- `al2_model_v69/`
- `si3_model_v69/`
- cores `al3_core_v69/` and `si4_core_v69/`

## Headline: Mg I → Al II → Si III, blind, same code, nothing retuned

| Ion | Blind intercombination line (3s² ¹S₀–3s3p ³P₁) | Model (ion SO, declared) | Lab | Δ (dex) |
|---|---|---|---|---|
| Mg I | 4571 Å | −5.650 | −5.686 (Kwong 1982) | +0.036 |
| **Al II** | **2669 Å** | **−4.988** | **−4.979** (Träbert et al. 1999 storage ring; NIST A+) | **−0.009** |
| **Si III** | **1892 Å** | **−4.569** | **−4.570** (Kwong et al. 1983 RF ion trap, A = 1.67(10)e4 s⁻¹ *) | **+0.002** |

\* The Si III lab value was taken from the secondary citation in A&A 2015 (aa24394-14). The primary paper is not yet read.

**Source check (added after v69, no numbers changed).** The cited sources are in the library and have rows in `data/refs/bibliography.csv`:
- **Si III:** the measuring group's own table (P. L. Smith, IAU Highlights, doi 10.1017/S0251107X00005332) lists 1892 Å A = 16.7(1.0) × 10³ s⁻¹ (Kwong et al. 1983a). Atalay et al. 2019 (A&A 631, A29) states it is the only measurement. The primary paper (Phys. Rev. A 27, 3040) is behind the APS paywall, and its page is confirmed by Crossref as 3040–3043.
- **Al II:** the measuring group's table also gives A = 3.33(0.23) × 10³ s⁻¹ (Johnson, Smith & Parkinson 1986, ApJ 308, 1013, held). Träbert's NIST review gives the storage-ring τ(³P₁) = 305 ± 10 µs, i.e. A = 3.28 × 10³ s⁻¹ and log gf = −4.978, the NIST A+ value we used. Träbert et al. 1999 (J. Phys. B 32, 537) is paywalled.
- **Mg I:** Kwong, Smith & Parkinson 1982 (Phys. Rev. A 25, 2629) is paywalled. The value is read via Bergemann 2017 and KP08 (both held).
- **Paywalled primaries not held:** Kwong 1982, Kwong 1983 and Träbert 1999. Del Zanna 2015 is free, but the A&A site blocks automated download.

The ab initio spin–orbit variant gives −5.824 (Mg I), −5.054 (Al II) and −4.575 (Si III).

The p-shell ion calibration is what matters, and it shrinks along the sequence: p scale 1.22 → 1.08 → 1.01. A one-body ξ becomes more nearly complete as the ion charge rises, which is the expected trend.

Spin-allowed lines vs NIST graded gf (mostly MCHF compilations, so a theory referee, not lab):

| Ion | n | Median | abs68 | Within 1σ / 2σ |
|---|---|---|---|---|
| Mg I (lab, Pehlivan) | 24 | +0.020 | 0.022 | 83 / 96 % |
| Al II (NIST A+/A/B+) | 246 | −0.001 | 0.021 | 90 / 97 % |
| Si III (NIST A/B+) | 107 | −0.004 | 0.022 | 87 / 96 % |

By grade, Al II: A+ (n = 5) 0.018; A (109) 0.020; B+ (132) 0.024 (abs68). The tighter NIST grades agree better, as they should.

Energy gates pass:

| Ion | Max energy error | Core potential fit (rel. err.) |
|---|---|---|
| Mg I | 0.45 % | |
| Al II | 0.38 % | 3 × 10⁻⁵ |
| Si III | 0.91 % | 3 × 10⁻⁵ |

The cores were fitted exactly as Mg II was in v63, with α_d = 0.2649 / 0.1624 a.u. (Mitroy, Safronova & Clark 2010).

Declared rules that carry along the sequence (no theory refit): where there are no lab-independent lines, the floor (0.050) and the CF-bin term (0.155 / 0.34) are taken from Mg I.

Not run for the sister ions: convergence variants (Mg I shows ≤ 0.03 dex at p90).

### Findings

1. **The core d spin–orbit scale is not physical.** In Mg II, Al III and Si IV every nd ²D term is inverted in NIST (J = 5/2 below 3/2). That comes from core polarization, which a one-body ξ cannot produce.
   - The per-l d scale (Mg 1.10, Al 0.39, Si 0.02) is therefore arbitrary.
   - It does not affect the blind tests (set by p).
   - The SO-treatment σ term (|ion − ab initio|) already prices it.
   - **Proposal for v70 (to be declared): an l whose core fine structure is inverted stays unscaled.**
2. **The NIST Si III 1501.150 Å entry (3s3d ³D₃–3s4f ³F₃, ref T3003, B+) looks wrong.**
   - Given A = 2.75e9 s⁻¹, which is stronger than the multiplet's main line 3D₃–3F₄ (2.10e9). LS predicts about a tenth of that.
   - Our model gives log gf −0.25, against NIST +0.81. This looks like a J mis-assignment in the source.
   - Noted for the line-list ticket (RYA-1393 style); no action here.
3. **The worst Al II / Si III misses are high-l Rydberg lines** (ng–nf, 5g–6f; 0.2–0.3 dex). These are the LS-broken configurations, and the NIST values there are LS-coupled theory (T6298c4). So the misses are not clearly ours.

## A1: same-lower-level line ratios in the observed Sun (abundance-free)

**Declared primary (weak limit, Gaussian EW): no power.**
- 25 lines in 8 groups.
- Method null 0.25 dex.
- 0 discriminating lines.

Cause, measured: solar Mg I lines are saturated far below log(W/λ) = −5. In a MARCS LTE curve of growth, 6318 Å (37 mÅ) has slope 0.67 and 7659 Å (80 mÅ) has 0.43. Gaussian EWs also miss the van der Waals wings, so the observed ratios come out compressed (e.g. the 4p ³P–6s ³S triplet at 1:3:5).

**A1b (post-hoc, disclosed): Voigt EW inverted through a per-line LTE curve of growth** (iSpec/SPECTRUM, MARCS 5772/4.44/0, ξ 1.0).
- The 12418 / 12423 / 12433 Å triplet (CRIRES+ J) then reproduces the pure-LS ratios to **≤ 0.03 dex**; without the curve of growth the error was ≈ 0.2.
- Over all groups (11 lines) the null is 0.18. The weakest and blended lines dominate it (7753 Å at 1.9 mÅ).
- **Still no discriminating lines.** Every clean, weak Mg I group in the Sun is LS-pure, and there all gf sources agree within 0.02. So the Sun cannot referee Mg I's relative gf.
- The method is validated and moves to Ca I and Fe I, which have many more lines and real intermediate coupling.

Script: `scripts/rya1354_solar_ratio_referee.py` (`--mode=cog`). No gf is changed.

**By-product:** the codex IR iSpec line lists carry **no van der Waals constants for Mg I** (all zero, Turbospectrum-only). This is input for v70 (broadening from our own wavefunctions).

## A2: composition vs NIST level percentages (report only)

| Ion | n levels | Median model − NIST (pct points) |
|---|---|---|
| Mg I | 16 | −1.0 |
| Al II | 11 | −1.2 |
| Si III | 10 | 0.0 |

The one notable miss is Mg I 3s3d ¹D, 3p² admixture 16.6 vs 22 %.

NIST lists **"3d² ³S" (2 %) in 3s4s ³S, which is Pauli-forbidden** (d² allows only L + S even). That is a basis-notation artifact; the row is flagged and excluded from the median.

## A3: convergence rule (declared; disclosed as following the v68 rcomp35 result)

- A convergence variant enters σ_conv only if its energy error is ≤ primary + 0.10 %.
- New variant rcomp20: energy **2.13 %**, i.e. worse. So RCOMP = 25 a₀ is an optimum, not a convergence knob.
- rcomp35 (0.87 %) and rcomp20 are excluded and reported (p90 |Δ| 0.63 / 0.81 dex).
- More compact functions at fixed box (ncomp18) moves lines by ≤ 0.02 dex. The compact space is complete along the axis that refines the basis.
- **Effect:** σ_conv p90 for n 11–12 goes from 0.34 to 0.030. Mg I class C goes from 1,318 to 1,359. Lab coverage is unchanged (83 / 96 %).

## Also

- RYA-1393 audit rerun on v69: unchanged (10 refuted, 4 confirmed).
- **Tests:** 15 pass (+ configuration parsing, the Pauli rule, the convergence rule, sister-ion labels).

## Open for Ryan

None blocking.

## v70 proposals (to be declared)

1. Inverted-core-FS rule for the d scale.
2. Roadmap item 3: radiative damping, van der Waals broadening (ABO from our ⟨r²⟩ / n*), photoionization.
3. Read Kwong 1983 directly.
