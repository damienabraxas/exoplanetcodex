# RYA-531 — per-element NLTE engine map (input for RYA-525's two-engine floor)

**Date:** 2026-07-06 · **Branch:** `ryandamienschmitt/rya-531-...` (stacked on RYA-529 `dcbc741`) · **No merge — Ryan reviews.**

## The architectural finding that reframes this map

RYA-402's "two-family" split (Family A TS-native Gerber / Family B PySME) is a map of **how the
NLTE _delta_ is derived**, not of two production synthesis engines. Investigation of merged
`main` (RYA-531 Step 0) establishes:

- **Production NLTE is universally a pre-computed 1D delta CSV** interpolated over (Teff, logg,
  [Fe/H]) and **added to the Turbospectrum/EW LTE abundance** (`pipeline/nlte_corrections.py`
  `apply_element_nlte_corrections`; `pipeline/nlte_cno.py` for C/O). There is **no bsyn
  in-synthesis NLTE anywhere** — grep for `set_nlte` / departure-feed to bsyn in production = 0.
- Therefore the two engines of RYA-525 are: **Engine A = 1D-NLTE (EW + delta)** and **Engine B =
  Turbospectrum spectral synthesis**. Engine B is **LTE-only for every element today.** Giving
  Engine B an NLTE mode requires the **Gerber TS-native departure machinery** (Gerber, Bergemann
  et al. **2023, A&A 669, A43**; arXiv 2206.00967 — note: the RYA-531 ticket's "A&A 666 A18"
  citation is wrong) — `NLTE .true.` + a bsyn `NLTEINFOFILE` whose per-element rows carry the
  model-atom + departure-coefficient file paths (MODELATOMFILE/DEPARTUREFILE are columns of the
  NLTEINFOFILE, not separate top-level keys, in the released `Turbospectrum_NLTE`). This is
  **NOT built and NOT provisioned** (empty `gerber_ts/` on Sirius; `nlte_bfactor_synth.synth_ew_nlte_vs_lte`
  is the abandoned Family-B TS stub, a `NotImplementedError`).
- The RYA-402 "Family A TS-native" route was described as "the route" but **never executed**:
  Na/Mg/Si were instead re-sourced onto Amarsi-2020 **PySME**-derived deltas (RYA-410), and
  Ca/Ti/Mn/Ba/Sr kept their pre-existing published delta grids. So Family A's NLTE is fully
  served — just not by TS-native synthesis.

**Bottom line for RYA-525:** Engine-B **NLTE** is the genuine unbuilt gap (all 27 elements). It is
a from-scratch RT-engine build, not a re-provision. Engine-A 1D-NLTE is well covered (see below).

## Per-element map (all 27 targets, RYA-109 `TARGET_ELEMENTS`)

Legend — **Engine A source flag:** the production `NLTE_CORRECTION_ELEMENTS` / `nlte_cno` engine.
**Engine B (synth) NLTE:** TS-native NLTE status = **NONE-BUILT** for all (LTE synth only today).

| # | Element | Engine A (1D-NLTE delta) — current production source | RYA-402 family | Engine B (TS synth) NLTE |
|---|---|---|---|---|
| 1 | Fe | MPIA/Amarsi Fe grid, dedicated Fe leg (`Fe_Bergemann_MPIA.csv`, RYA-319) | — (Fe handled pre-402) | NONE-BUILT (LTE synth) |
| 2 | C | `nlte_cno` Amarsi 2019 CDS 1D-NLTE C grid | — (nlte_cno) | NONE-BUILT (LTE synth) |
| 3 | O | `nlte_cno` Amarsi 2019 CDS 1D-NLTE O grid | Family A | **TS-native NLTE VALIDATED** (RYA-534; O 777 δ −0.105 vs Amarsi-2019 1D −0.134) |
| 4 | Mg | `Mg_Amarsi2020_PySME.csv` — **PySME**-derived (RYA-410) | Family A | **TS-native NLTE VALIDATED** (RYA-534; 5711 δ −0.023 vs PySME −0.022) |
| 5 | Si | `Si_Amarsi2020_PySME.csv` — **PySME**-derived (RYA-410) | Family A | **TS-native NLTE VALIDATED** (RYA-534; 5772 δ −0.034 vs PySME −0.013) |
| 6 | Ca | `Ca_Mashonkina2017.csv` — MPIA MAFAGS-OS delta grid | Family A | **TS-native NLTE VALIDATED** (RYA-534; 6122/6162 δ −0.009 vs Mashonkina +0.017) |
| 7 | Ti | `Ti_Bergemann2011_MPIA.csv` — MPIA MAFAGS-OS delta grid | Family A | NONE-BUILT (LTE synth) |
| 8 | Ni | **LTE** (Engine-A: no grid) | Family A | **TS-native NLTE VALIDATED** (RYA-534; 5018/5035 δ +0.018 vs Bergemann +0.02) |
| 9 | Na | `Na_Amarsi2020_PySME.csv` — **PySME**-derived (RYA-410; RYA-529 Sirius gate δ −0.129) | Family A **and** B (dual-validator) | **TS-native NLTE VALIDATED** (RYA-533; 5682/5688 δ −0.068 vs −0.107) |
| 10 | P | **LTE** — no NLTE grid (not an Amarsi/MPIA element) | — | NONE-BUILT (LTE synth) |
| 11 | S | `S_Amarsi2025_PySME.csv` — **PySME**-derived (RYA-402) | Family B | NONE-BUILT (LTE synth) |
| 12 | N | `N_Amarsi2020_PySME.csv` — **PySME**-derived on Sirius (RYA-526) | Family A | NONE-BUILT (LTE synth) |
| 13 | Co | **LTE** (Engine-A: no grid) | Family A | **TS-native NLTE VALIDATED** (RYA-534; 5000/5013 δ +0.099 vs Bergemann +0.10; weak lines) |
| 14 | Cr | `Cr_Bergemann2010_MPIA.csv` — MPIA MAFAGS-OS delta grid | — (MPIA/RYA-235) | NONE-BUILT (LTE synth) |
| 15 | Al | `Al_Amarsi2020_PySME.csv` — **PySME**-derived (RYA-402) | Family B | NONE-BUILT (LTE synth) |
| 16 | K | `K_Amarsi2020_PySME.csv` — **PySME**-derived (RYA-402/462) | Family B | NONE-BUILT (LTE synth) |
| 17 | Ba | `Ba_Korotin2015.csv` — Korotin CDS delta grid | Family A | **TS-native NLTE VALIDATED** (RYA-534; Ba II 4554 δ −0.018 vs Korotin) |
| 18 | Y | **LTE** — no NLTE grid | — | NONE-BUILT (LTE synth) |
| 19 | V | **LTE** — true NLTE void (RYA-526 → build RYA-470) | — | NONE-BUILT (LTE synth) |
| 20 | Cu | **LTE** in production — PySME machinery validated but NOT registered (junk measured lines, RYA-402/395) | Family B | NONE-BUILT (LTE synth) |
| 21 | Mn | `Mn_Bergemann_MPIA.csv` — MPIA MAFAGS-OS delta grid | Family A | **TS-native NLTE VALIDATED** (RYA-534; 6013/6021 δ +0.043 vs Bergemann +0.107 / ~½; RYA-411 xref) |
| 22 | Sc | **LTE** — no NLTE grid | — | NONE-BUILT (LTE synth) |
| 23 | Li | **LTE** — deferred (RYA-103 Li 6707 CN-blend) | Family B (deferred) | NONE-BUILT (LTE synth) |
| 24 | Eu | **LTE** — no NLTE grid | — | NONE-BUILT (LTE synth) |
| 25 | Zr | **LTE** — dominant-ion, no grid needed (RYA-526) | — | NONE-BUILT (LTE synth) |
| 26 | Sr | `Sr_Bergemann2012_INSPECT.csv` — INSPECT delta grid (Sr II 4077/4215) | Family A | **TS-native NLTE VALIDATED** (RYA-534; Sr II 4215 δ −0.013 vs INSPECT ~−0.005) |

*(26 unique symbols = 27 targets; Fe covers Fe I + Fe II.)*

## Summary counts (Engine A / 1D-NLTE)

- **PySME-derived Amarsi delta CSV (8):** Na, Mg, Si, Al, K, S, N — plus Cu (validated, unregistered).
- **Published MPIA/CDS/INSPECT delta grid (6):** Fe, Ca, Ti, Cr, Mn (MPIA); Ba (Korotin); Sr (INSPECT/Bergemann).
- **`nlte_cno` Amarsi-2019 (2):** C, O.
- **LTE (no NLTE, 9):** Ni, P, Co, Y, V, Sc, Li, Eu, Zr (+ Cu in production).

## Engine B (TS-native NLTE synthesis) — status: **DECK BUILT + Na-VALIDATED (RYA-533); 1/11 Family-A provisioned**

**RYA-533 (2026-07-06):** the Engine-B TS-native NLTE deck is **built and validated on Sirius**.
Turbospectrum_NLTE (BSYN v20.1) compiled; the **Na Gerber-2023** departure grid (15.9 GB) + model atom
(`atom.na_qmh`) + aux provisioned (md5-pinned, Sirius-only); bsyn NLTE deck wired (interpol departures →
babsma → bsyn NLTE-vs-LTE COG). **Na gate GREEN: median δ −0.068 vs anchor −0.107 (PASS, model-atom tol
0.05); validate-don't-tune.** Cross-engine Na = INSPECT −0.107 / PySME −0.129 / TS-Gerber −0.068 (the
RYA-525 model-atom-systematic diagnostic). Deck `scripts/ts_gerber_na_gate_rya533.py`; gate test
`tests/test_ts_gerber_nlte_rya533.py`.

Gerber-2023 (A&A 669 A43) provides TS-native NLTE for 11 of the 27: **O, Na, Mg, Si, Ca, Ti, Mn, Co, Ni,
Sr, Ba** (+ Fe, H, Al, Y wider). **Na is provisioned + validated (1/11);** the other 10 are a **mechanical
per-element repeat** — download the element's Gerber grid + atom (Keeper), run the interpolator, swap the
line list, re-run the gate. So for RYA-525's Engine-B-NLTE column: Na = REAL now; O/Mg/Si/Ca/Ti/Mn/Co/Ni/
Sr/Ba = deck-ready, grid-download-pending. `gerber_ts/` on Sirius (empty at RYA-531) now holds the Na set.

**Provisioning feasibility (RYA-531 recon — data IS available, but the build is from-scratch):**
- Gerber 2023 (A&A 669, A43) publicly releases native model atoms + departure-coefficient grids
  (1D-MARCS and avg-3D-STAGGER, 4D over Teff/logg/[Fe/H]/abundance) for **H, O, Na, Mg, Si, Ca,
  Ti, Mn, Fe, Co, Ni, Sr, Ba** — all 11 Family-A elements covered.
- Hosting: **MPG Keeper (Seafile) share** `https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/`
  (subfolders `dep grids`, `model_atoms`). Openly downloadable per-element via TSFitPy's
  `download_nlte_grids.py`. **Caveats vs. our provenance bar: no DOI, mutable share, grids
  unversioned** → must checksum whatever is pulled. Per-element grids are multi-GB (~20–100 GB
  for a handful; ~1 TB all) → Sirius-only (RYA-526).
- Code: `github.com/bertrandplez/Turbospectrum_NLTE` (v20 NLTE) — must be **built from source on
  Sirius** (separate from the current LTE TS build; expect the RYA-517 gcc-flags class of issues).
- **Effort: this is a from-scratch RT-engine build** — compile Turbospectrum_NLTE, download the
  Gerber grids+atoms on Sirius, write the babsma/bsyn NLTE deck + 4D departure interpolation into
  the pipeline (none exists today), then run Na through it and validate vs the anchor. It is NOT a
  download-and-wire re-provision. RYA-402 spent a whole ticket here and deliberately chose PySME
  instead (skipped the TS model-atom route). **Decision to commit this build = pending Ryan
  (RYA-531 Linear thread).** Recommended as its own dedicated build ticket, since Family A is
  already NLTE-covered via Engine A (so this is a NEW Engine-B capability for 525's cross-engine
  comparison, not a coverage-blocking gap).
