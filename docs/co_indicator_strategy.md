# C/O Indicator Strategy — The Exoplanet Codex

**Method document — carbon and oxygen abundances and the C/O ratio**
Status: draft for review · Owner: Ryan Schmitt · Both the optical and IR abundance arms report into this document.

---

## 0. Why this document exists

C/O is the Codex flagship metric — the bridge from stellar chemistry to rocky-planet/habitability inference. It deserves a *standing method*, not a per-session improvisation. This document fixes three things:

1. **The mission is good C and good O — not CO.** CO is one road to carbon and oxygen, and (see §2) it is the least precise and most model-sensitive one. We pursue a *panel* of indicators per element, not a hero line.
2. **The robust C/O comes from matched-pair construction** (§3), where shared temperature sensitivity cancels systematics in the ratio even when individual abundances carry offsets.
3. **Different indicators own different wavelength regions and temperature regimes** (§4). Agreement across the panel *is* the systematic meter — the "quiltwork of truth." This is the input side of the cross-region presentation architecture (RYA-237 / abundance-presentation): cross-indicator and cross-region spread is flagged and adjudicated, never silently averaged.

The concrete per-line wavelengths, log gf, and broadening are **not** fixed here — they come from the per-target VALD extraction and the cited NLTE/molecular line data (single source of truth). This document fixes the *strategy*: which indicator families, how to combine them, and which dominate for which target.

---

## 1. Reference frame and anchors

Solar anchors, adopted from the weighted-mean-over-panel of Amarsi et al. 2021 (`aa4044521`, Table 3):

| Element | A(X)☉ | Source |
|---|---|---|
| C | 8.46 ± 0.04 | weighted mean of [C I], C I, C₂, CH, CO |
| O | 8.69 ± 0.04 | weighted mean of [O I], O I, OH |
| N | 7.83 ± 0.07 | atomic N I + molecular NH/CN |

These are **validation targets, never tuning targets.** Differential [X/H] uses a consistently NLTE-corrected solar anchor through the same apply path; the C/O ratio itself is largely immune to the anchor (RYA-282).

---

## 2. The indicator panel (per element)

Precision and model-sensitivity from Amarsi et al. 2021 (`aa4044521`, Table 3); error column is the combined statistical + systematic per indicator.

### Carbon — ranked best-behaved → worst

| Indicator | Region | Regime | σ (dex) | Notes |
|---|---|---|---|---|
| C₂ Swan (5135, 5165 Å) | optical | molecular, LTE-OK | **0.029** | Teff-sensitive; strong cross-check |
| [C I] 8727 Å | optical | forbidden, **NLTE-immune** | **0.034** | well-reproduced in 3D-LTE; the clean differential anchor |
| C I (high-excitation) | optical + Y/J/H IR | atomic, **large −NLTE** | **0.038** | NLTE grows toward metal-poor (Fabbian, `Carbon.pdf`) |
| CH (G-band 4300 Å; Δν=1 rovib) | optical / IR | molecular, LTE | 0.043–0.048 | strong in metal-poor |
| CO fundamental (Δν=1, 4.6 µm) | IR (L/M) | molecular | 0.101 | model-sensitive; ground-based access poor |
| **CO first overtone (Δν=2, 2.3 µm)** | IR (K) | molecular | **0.122** | **largest error of any C indicator; ~0.28 dex spread across model atmospheres (3D 8.467 vs HM 8.743); telluric-cursed** |

### Oxygen

| Indicator | Region | Regime | σ (dex) | Notes |
|---|---|---|---|---|
| O I 777 nm triplet | optical | atomic, high-exc, **large −NLTE (~−0.12)** | **0.030** | the FGK workhorse; 3D-NLTE grid from Amarsi 2016 (`stv2608`) |
| [O I] 6300, 6364 | optical | forbidden, **LTE-safe** | 0.051 | weak; 6300 blended with Ni I 6300.34 → Johansson 2003 gf + consistent A(Ni) (RYA-365); the model-independent anchor |
| O I 6158 / 8446 / 9266 | optical/NIR | atomic, high-exc | — | weaker auxiliaries (Amarsi 2021 line set) |
| OH vib-rotation (Δν=0,1,2; ~1.5–2.0 µm; Δν=1 ~3 µm) | IR (H + L) | molecular | 0.06–0.09 | model-sensitive; valuable in cool stars; needs telluric |
| OH electronic (3100–3200 Å) | UV | molecular | — | strong but UV access + crowding (HST/STIS territory) |

**Takeaway:** precise C and O live in the **optical** — C₂/[C I]/C I for carbon, O I 777 + [O I] for oxygen. CO is 3–4× less precise than the optical carbon indicators and the most model-dependent feature in either panel.

---

## 3. C/O construction — the matched-pair principle

C/O is more trustworthy than either abundance alone *because* shared temperature sensitivity cancels in the ratio (Fabbian, `Carbon.pdf`: "the [C/O] ratio should be relatively unaffected, given the similar temperature dependence of the C and O lines we employ"). We exploit this deliberately.

**Published C/O = best C ⊗ best O regardless of region** (per the presentation architecture), with the constructions below as consistency checks:

| Construction | C indicator | O indicator | Property |
|---|---|---|---|
| **High-excitation atomic** (primary, warm stars) | C I | O I 777 | correlated Teff + NLTE response → robust ratio |
| **Forbidden** (model-independent anchor) | [C I] 8727 | [O I] 6300 | both NLTE-immune, LTE-safe; weak but clean |
| **Molecular** (cool stars) | CH / C₂ / CO | OH | accessible where atomic lines fade |

Rules:
- **Matched-pool discipline:** do not mix construction families without flagging; pool-composition artifacts distort spread attribution.
- **Molecular degeneracy:** CO ≈ f(C × O). To extract C from CO, **fix O from an independent indicator (OH / [O I] / O I 777) and propagate its error** (Amarsi 2021 method); symmetric for O from CO. CO never yields C *or* O standalone.
- **Disagreement is data:** cross-construction or cross-region spread beyond errors is flagged and adjudicated, never averaged away — it is the systematic meter.

---

## 4. Per-temperature-regime strategy

Which panel carries the weight shifts with Teff. This sets where each target's effort goes.

| Regime | Targets | Carbon | Oxygen | C/O primary | IR role |
|---|---|---|---|---|---|
| Warm FGK | Sun, **α Cen A** (G2V), Procyon (F5 IV–V) | C I + [C I] + C₂ | O I 777 + [O I] | high-exc atomic pair; forbidden-pair cross-check | optional cross-check + ¹³C bonus |
| Cool K dwarf | **α Cen B** (K1V), **55 Cnc A** (K0 IV–V) | CH + C₂ + (CO if IR clean) | [O I] + OH (IR) | molecular ⊗ molecular, anchored by optical where present | **load-bearing** — atomic lines weaken |
| M dwarf | future targets | CO, OH, CN bands (MARCS.GES) | OH, CO | molecular-only | **essential** |

**Strategic consequence:** the IR arm (RYA-373 / 380 / 390 / 391) pays its biggest dividend on the *cool* targets and the future M-dwarf program — not on the Sun or α Cen A, where optical already delivers a precise, robust C/O. This is the principled justification for building the IR capability now while deferring it off the warm-star critical path.

---

## 5. Calibration and the gold standard

- **3D-NLTE is non-negotiable for the molecular indicators** (CO, OH). The 0.28 dex 3D-vs-1D CO spread (§2) is the proof — CO forms deep and reads the temperature structure. Atomic O I 777 and C I require their NLTE grids (Amarsi 2016 `stv2608` for O; Amarsi 2019 / Fabbian for C).
- **Solar IR gold standard** = telluric-free space atlas (ACE-FTS / Hase 2010, RYA-390) + disk-resolved FTS atlas (Kitt Peak photatl) + a 3D-NLTE solar model, validated by **reproducing Elgueta et al. (`2602.14294`) robust-line results on the shared Vesta CRIRES+ data** (their G-dwarf benchmark is Vesta/Sun, obs 2022-11-22 — the same data we hold).
- **IR telluric recipe:** cr2res + molecfit + per-night GDAS (RYA-380), with loud-fail on standard-atmosphere fallback. Apply a **per-band robustness gate** (depth / non-saturation / purity / goodness-of-fit, evaluated independently per band — Elgueta) to the IR line list.
- **Note on K-band scope:** the field's high-resolution IR abundance work concentrates in **Y/J/H**, not K (Elgueta). K-band CO is a specialist C/O/¹³C diagnostic; broader, better-behaved IR C and O signal lives in the Y/J/H C I and OH lines. "Good IR" ≠ "K-band CO."
- **Benchmark ladder:** Sun → Procyon (α CMi) → α Cen A/B → 55 Cnc A, each read against its GBS reference panel. Procyon and the Sun are both in Elgueta's GBS IR set — direct external validation points.

---

## 6. Standing operational rules

1. Every C and O result reports a **per-indicator panel**, not a single number.
2. C/O is published as **primary (best C ⊗ best O) + matched-pair constructions as consistency checks**, with a cross-region/cross-indicator agreement badge.
3. **CO is cross-check + isotopes (¹²C/¹³C), never the sole C or O.**
4. Molecular C or O is only derived with the partner element fixed from an independent indicator, error propagated.
5. Cross-indicator / cross-region disagreement beyond errors → **flag and adjudicate, never silently average.**
6. Solar anchors are validation targets, never tuning knobs.

---

## References (project library)

- Amarsi et al. 2021, A&A 653, A141 — solar C/N/O reanalysis; the indicator panel and anchors (`aa4044521`)
- Amarsi et al. 2016, MNRAS 455, 3735 — O I 3D-NLTE / LTE abundance-correction grid (`stv2608`)
- Amarsi et al. 2019 (CDS J/A+A/630/A104) — C/O 3D-NLTE and 1D-NLTE line-by-line correction grids
- Fabbian et al. — C I NLTE line formation in late-type stars; C/O matched-temperature robustness (`Carbon.pdf`)
- Elgueta et al. — selecting IR lines for abundance determination; per-band robustness gate; Vesta/Procyon CRIRES+ GBS (`2602.14294`)

*Per-line wavelengths, log gf, broadening, and HFS are sourced per target from VALD + the cited line data — not from this document.*
