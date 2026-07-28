# RYA-400 "The Beast" — per-element physics-regime map (all 27)

The homework that earns the RYA-371 Phase-C serious validation. Before the 27-element
solar abundance set is printed, every element's physics regime is settled — so each
residual in the 371 run is *genuine* (the element got its proper 1D/3D × LTE/NLTE
treatment), not an artifact of an unexamined default.

- **Single source:** [`config/physics_regime_rya400.yaml`](../config/physics_regime_rya400.yaml) — the structured per-element rows (cited).
- **Living check:** `scripts/audit_physics_regime_rya400.py` cross-checks the map against the live code — a `LOCKED` element that is **not actually NLTE-wired loud-fails**, coverage gaps fail, measured-line counts must match `data/measured/sol_ew_results_v1.csv`. **Current status: PASS** (12 LOCKED all wired; 26 symbols mapped; Fe I/II = the 27th species).
- **Principle:** validate-don't-tune. An element still off after its *correct* physics is **real information** — documented and **carried to the next star** (Procyon → α Cen → …), never tuned to close a residual.

## Regime authorities (value-provenance — no regime asserted from memory)

| Tag | Reference | Role |
|-----|-----------|------|
| Asplund2021 | Asplund, Amarsi & Grevesse 2021, **A&A 653, A141** | The 3D non-LTE solar-abundance authority; per-element regime calls (new 3D-NLTE for Na/Mg/K/Ca/Fe + Li/C/N/O/Al/Si/Ba). Also the repo's `SOLAR_ASPLUND2021` source. |
| Amarsi2020 | Amarsi et al. 2020, **A&A 642, A62** (GALAH) | Public NLTE departure-coefficient grids for **H, Li, C, N, O, Na, Mg, Al, Si, K, Ca, Mn, Ba** — the GET-GRID source for Al/K/Li. |
| Amarsi2019 | Amarsi, Nissen & Skúladóttir 2019, A&A 630, A104 | C I / O I (and N) 3D-NLTE grids (repo `nlte_cno`). |
| Repo grids | Lind2011 (Na), Osorio2015 (Mg), Korotin2015 (Ba), Bergemann2010/2011 (Cr/Ti), Bergemann&Gehren2008 (Mn), Mashonkina2017 (Ca), MAFAGS-OS (Si/Mg/Mn/Fe) | The vendored, currency-verified NLTE grids (RYA-165/235/256/319/396). |

## Verdict legend & tally

`LOCKED` (NLTE applied + wired) · `GET-GRID` (NLTE grid exists, vendor it) · `GET-DATA`
(measure/observe first) · `GET-3D` → routes to **RYA-399** (3D arm) · `HARD-carry-forward`
(document the limit, carry to the next star).

**Tally (round-1):** LOCKED **12** · GET-GRID **6** · GET-DATA **3** · HARD-carry-forward **3** · LTE-OK **2**.

> **Round-2 update (RYA-401):** the 6 GET-GRID elements were probed against the live
> sources — **0 were vendorable** (each a distinct documented gap; no fabrication), so they
> are flipped out of GET-GRID. **Current tally:** LOCKED **12** · GET-DATA **4** (Co, Sc,
> P, **Sr**) · HARD-carry-forward **8** (Li, Eu, Zr + **Al, K, S, Cu, V**) · LTE-OK **2**.
> See [`nlte_grid_acquisition_rya401.md`](nlte_grid_acquisition_rya401.md) + the `rya401:`
> blocks in the YAML. The GET-GRID table below is the round-1 audit (superseded).

---

## LOCKED — NLTE applied + verified wired (3D status confirmed/flagged)

| El | Dimensionality | Departure (solar Δ) | Grid (in repo) | Indicators (in data) | Hazards | Verdict |
|----|----------------|---------------------|----------------|----------------------|---------|---------|
| **Fe I/II** | 1D-OK, 3D small `[A21]` | NLTE over-ion. −0.01…−0.06 `[A21; Bergemann]` | `Fe_Bergemann_MPIA` (RYA-319) | many Fe I/II ✓ | Fe II = logg anchor | LOCKED |
| **C** | 3D needed `[A19/A21]` | NLTE up to −0.1 `[A19]` | `nlte_cno` ✓ | C I 5052/5380 + IR; CH | CH/C2 blends | LOCKED |
| **N** | 3D needed `[A19/A21]` | NLTE `[A19/A20]` | `nlte_cno` ✓ | N I 7468/8216; NH 3360; **synthesis, 0 EW** | no clean VIS N I EW | LOCKED |
| **O** | **3D applied** `[A16; Caffau15]` | NLTE: 777 +0.2…+0.3, [OI] small `[A16]` | `nlte_cno` ✓ | [O I] 6300 (+Ni), O I 777 | [OI]6300 Ni blend (RYA-365/367) | LOCKED |
| **Na** | 3D-NLTE `[A21]` | NLTE −0.05…−0.15 `[Lind11]` | `Na_Lind2011_INSPECT` (d=−0.107) | 5682/5688/6154/6160 ✓ | D lines saturated → use subordinate | LOCKED |
| **Mg** | 3D-NLTE `[A21]` | NLTE ~+0.007 `[Osorio15]` | `Mg_Bergemann_MPIA` | 5711/4730 ✓ | RYA-239 +1.72 = **saturation/curation, not NLTE** | LOCKED |
| **Si** | **3D → RYA-399** `[A21]` | NLTE ~−0.004 (clean) `[A21]` | `Si_Bergemann_MPIA` | many Si I (80) ✓ | 3D check owed | LOCKED |
| **Ca** | 3D-NLTE `[A21]` | NLTE ~+0.013 `[Mash17]` | `Ca_Mashonkina2017` | Ca I (29) ✓ | strong lines saturate | LOCKED |
| **Ti** | **3D-dominated → RYA-399** `[A21]` | NLTE ~+0.108 `[Berg11]` | `Ti_Bergemann2011_MPIA` | Ti I (91) + Ti II ✓ | Ti I/II NLTE differ | LOCKED |
| **Cr** | **3D-dominated → RYA-399** `[A21]` | NLTE ~+0.073 `[Berg10]` | `Cr_Bergemann2010_MPIA` | Cr I (80) ✓ | raw LTE reads HIGH; NLTE worsens uncurated → curate first (RYA-395) | LOCKED |
| **Ba II** | 3D-NLTE `[A21]` | NLTE ~−0.03 `[Korotin15]` | `Ba_Korotin2015` (ion II) | 4554/5853/6141/6496 ✓ | resonance saturated; HFS+isotopes | LOCKED |
| **Mn** | 3D-NLTE `[A21]` | NLTE ~+0.10 `[B&G08]` | `Mn_Bergemann_MPIA` | Mn I (9) ✓ | HFS-heavy (VALD HFS-on) | LOCKED |

*Even the LOCKED set carries its 3D status: Si/Ti/Cr are flagged 3D-dominated and routed to **RYA-399** (the 3D-extension arm); O is already 3D.*

## GET-GRID (round-1 audit — SUPERSEDED by RYA-401: all 6 → documented gaps)

| El | Departure (solar Δ) | Grid source | Indicators (in data) | Hazards | Verdict |
|----|---------------------|-------------|----------------------|---------|---------|
| **Al** | NLTE `[Nordlander&Lind17; A20]` | Amarsi2020 GALAH / Nordlander&Lind2017 | 6696/6698 subord. (2) | resonance 3944/3961 strong | GET-GRID |
| **K** ⭐ | NLTE strong (resonance) `[A20/A21]` | Amarsi2020 GALAH (K) | 7665/7699 (**0 measured**) | **triple-threat:** saturated + NLTE + **in O₂ telluric band** → needs RYA-380 telluric | GET-GRID |
| **S** ⭐ | NLTE line-dep: NIR triplet up to −1.1, opt −0.1, 8694 −0.26 `[Korotin09; Takeda05]` | Korotin/Takeda (no INSPECT/MPIA public grid; CDS/author or compute) | opt 6743/6757 vs NIR 9212/9228/9237 (7) | indicator decision: optical weak-but-small-NLTE vs NIR strong-but-−1.1 (Takeda/Korotin disagree ~0.4) | GET-GRID |
| **Cu** | NLTE significant; HFS `[Shi14; Korotin18]` | Korotin2018 / Andrievsky2018 (CDS/author) | 5105/5218/5782 (2) | HFS (VALD HFS-on) | GET-GRID |
| **V** | NLTE moderate; HFS-dominated `[A21]` | gap (author/compute) | V I (7) | HFS likely > NLTE; confirm | GET-GRID (low) |
| **Sr II** | NLTE resonance; LTE ion-disc. ~0.5 `[Bergemann12]` | Bergemann2012 (MPIA / ChETEC-INFRA 3DNLTE) | 4077/4215 resonance (1) | strong + blue-edge | GET-GRID |

## GET-DATA — measurement/observation must precede physics

| El | Note | Indicators | Verdict |
|----|------|-----------|---------|
| **Co** ⭐ | **0 measured lines today.** NLTE+HFS both required: solar NLTE−LTE ≈ +0.14, A(Co)=4.95 `[Bergemann10]`. Grid gettable (Bergemann10) but **measure first**. | Co I (HFS-heavy) — 0 | GET-DATA |
| **P** | FUV needs HST/STIS (**RYA-119**); NIR P I 10511/10581 weak. CHNOPS data-gap. | P I NIR/FUV (2 marginal) | GET-DATA |
| **Sc II** | **0 measured lines** — measure Sc II first; then small-NLTE/LTE-OK `[Zhang08]`; HFS. | Sc II 5031/5526/5657 — 0 | GET-DATA |

## LTE-OK — ionised-majority / well-behaved (small NLTE, confirmed)

| El | Note | Verdict |
|----|------|---------|
| **Ni** | Fe-peak, small NLTE `[A21]`; many Ni I lines (26). LTE acceptable at solar; MPIA Ni I grid optional (low priority). | LTE-OK |
| **Y II** | Ionised majority → small ionisation NLTE `[A21]`; 4883/5087/5200 (3). | LTE-OK |

## HARD-carry-forward — document the limit, carry to the next star

| El | Limiting issue | Verdict |
|----|----------------|---------|
| **Li** | Li I 6707 **CN molecular-blend** contamination (EW 10.2 vs ~2 expected) — **RYA-103**, blend-limited & deferred. NLTE grid exists (Lind09/Amarsi20). | HARD |
| **Eu II** | Eu II 6645 EW **0.3 mÅ — below noise floor** (expected ~8) — **RYA-102**. HFS-dominant; ~LTE. | HARD |
| **Zr II** | Saturation + blue-edge log gf calibration — **RYA-161**. ~LTE. | HARD |

*Carry-forward rationale: a benchmark in a different regime (cooler, enhanced, higher-S/N) may unveil what the Sun won't — the value is across the multi-star arc, not papered over.*

---

## Spotlight lurkers (most likely to bite Phase C)

- **S (sulfur) — CHNOPS-critical, no handling.** NLTE is line-dependent and severe in the NIR: the strong 9212/9228/9237 Å triplet reaches ≈ **−1.1 dex** non-LTE, while the weak high-EP optical lines (6743/6757) are small (≈−0.1). The **indicator decision** matters more than the grid: optical = clean but weak; NIR = strong but huge-NLTE and Takeda 2005 vs Korotin 2009 disagree by ~0.4 dex at 9213. *Decide indicator (RYA-162) + source S NLTE before Phase C trusts a sulfur number.*
- **K (potassium) — triple-threat.** 7665/7699 resonance doublet: **saturated + NLTE + inside the O₂ telluric A/B band**. Needs the RYA-380 telluric recipe applied (red-optical wavelength-gated) **and** the Amarsi 2020 K grid **and** a measurement (0 lines today).
- **Co (cobalt) — measurement precedes physics.** **0 measured lines** today, HFS-heavy. The Bergemann 2010 NLTE+HFS grid is gettable, but nothing can be applied until Co I is measured.

## What this gates / spawns

- **Blocks RYA-371 Phase C** — the serious 27-element validation runs once every element's regime is decided and the gettable grids/data are in hand.
- **Spawns grid-vendoring children (RYA-396 pattern):** Al, K, S, Cu, Sr (and Co after measurement) → INSPECT/MPIA/CDS/Amarsi-GALAH scrapers.
- **Feeds RYA-399 (3D arm):** Si, Ti, Cr (3D-dominated) + any element the audit flags 3D.
- **GET-DATA dependencies:** Co (measure), Sc (measure), P (RYA-119 HST/MAST).
- **HARD-carry-forward → 371 reports honestly, never tuned:** Li (RYA-103), Eu (RYA-102), Zr (RYA-161).

## Citations

- Asplund M., Amarsi A. M., Grevesse N. 2021, A&A 653, A141 — *The chemical make-up of the Sun: A 2020 vision.*
- Amarsi A. M. et al. 2020, A&A 642, A62 — *The GALAH Survey: non-LTE departure coefficients for large spectroscopic surveys.*
- Amarsi A. M., Nissen P. E., Skúladóttir Á. 2019, A&A 630, A104 (C/O); Amarsi et al. 2016, MNRAS 463, 1518 (O I 777).
- Lind K. et al. 2011, A&A 528, A103 (Na); Lind K. et al. 2009, A&A 503, 541 (Li).
- Osorio Y. et al. 2015, A&A 579, A53 (Mg); Mashonkina L. et al. 2017, A&A 606, A147 (Ca).
- Bergemann M. & Cescutti G. 2010, A&A 522, A9 (Cr); Bergemann M. 2011, MNRAS 413, 2184 (Ti); Bergemann M. & Gehren T. 2008, A&A 492, 823 (Mn); Bergemann M. et al. 2010, MNRAS 401, 1334 (Co); Bergemann M. et al. 2012, A&A 546, A90 (Sr).
- Korotin S. A. et al. 2015, A&A 581, A70 (Ba); Korotin S. A. 2009, Astron. Rep. 53, 651 (S); Korotin S. A. et al. 2018, MNRAS 480, 965 (Cu).
- Takeda Y. et al. 2005 (S NIR triplet); Nordlander T. & Lind K. 2017, A&A 607, A75 (Al); Shi J. R. et al. 2014 / Andrievsky S. et al. 2018 (Cu); Zhang H. W. et al. 2008 (Sc II).
- Caffau E. et al. 2015, A&A 579, A88 (solar [O I] 3D, RYA-367).
