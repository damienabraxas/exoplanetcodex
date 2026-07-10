# RYA-535 — Ti same-atom systematic: RCA findings

**Date:** 2026-07-09 · **Branch:** `ryandamienschmitt/rya-535-...` (on the RYA-534 tip 62d7329) · **NO MERGE.**

## The premise (why RYA-534's Ti label was wrong)

RYA-534 (v15) recorded Ti as a "genuine ~2× model-atom difference (Gerber-2023 atom.ti503b vs
Bergemann-2011)". **Refuted by source:** Gerber et al. 2023 (A&A 669 A43), Table 1 lists the Ti NLTE
model atom as **Bergemann (2011)**. `atom.ti503b` is the release's filename for that atom, not a
separate Gerber atom. Engine-A (our banked `Ti_Bergemann2011_MPIA` anchor +0.108) and Engine-B
(TS-Gerber synthesis, +0.221) therefore use the **same atom** → a model-atom difference is
impossible. +0.221 vs +0.108 is a **same-atom systematic**. Part A reverted the label (test xfail,
register v16); Part B is this RCA.

## B1 — Engine-A provenance (inherited, RYA-256/235/244-245)

+0.108 anchor = the MPIA per-line **Bergemann-2011** grid `Ti_Bergemann2011_MPIA.csv`, scraped via
`build_nlte_grids_mpia.py` (RYA-244/245); atmosphere = **MAFAGS-OS**; solar Ti I δ = +0.108 dex;
feh axis relative, unit-convention verified (RYA-256). ⇒ **Decision-tree CASE (ii): SAME atom
(Bergemann 2011), DIFFERENT atmosphere** — Engine-A MAFAGS-OS vs Engine-B MARCS.

## Cheap defect-class checks (both CLEAN — ruled out)

1. **Unit convention** (Fe PR#8 silent-clamp class). The Ti departure grid's abundance axis at the
   solar node spans absolute A(Ti) ∈ [4.0, 5.4] (16 nodes); the deck interpolates at A(Ti)=4.90 —
   **inside range, no clamp**. The deck passes absolute A(Ti) to babsma/bsyn + the interpolator's
   `aburef`. Not the cause.
2. **Grid placeholder/zero/NaN** (RYA-413 class). The solar departure file has **zero** NaN/inf and
   **zero** exact-zeros; **all 503/503 levels thermalise to b→1.0 in the deep layers** (median deep
   b = 1.000) — the unforgeable correctness signature. The grid is fully physical. Not the cause.

## B2 — b-factor dump (Engine-B, MARCS) — `Ti_bfactor_dump.txt`

Solar node echoed from `STAR_PARAMS['solar']` (5772 / 4.438 / 0.0 / 1.0), NOT hardcoded; grid
interpolated at the nearest MARCS node 5750/4.5/0.0, A(Ti)=4.90. (Engine-A is a MAFAGS-OS **scalar**
delta — it has no per-level b-factors to dump; stated, not faked.)

**Dump-tool bug found + fixed (mine, not the deck's):** `interpol_modeles_nlte.f` (line 772) writes
the departures **depth-major** (`do depth: write all n_lev`); my first reader reshaped level-major →
garbage (ground state didn't thermalise, wild b). Fixed to `(ndepth, nlevel).T`; ground state now
thermalises 0.062→1.000, all 503 levels →1 deep. **This is a bug in the RCA dump tool only — bsyn
reads the file via its own `read_departure.f`, so the gate's +0.221 is unaffected.**

Corrected b-factors, line-forming layer (logτ ∈ [−1, 0]) for the 3 gate lines:

| line | ⟨b_lower⟩ | ⟨b_upper⟩ | deep b |
|---|---|---|---|
| 5689.460 | 0.813 | 0.895 | →1.000 |
| 5648.565 | 0.822 | 0.865 | →1.000 |
| 5662.150 | 0.814 | 0.857 | →1.000 |

The departures are **physical and modest** (b_lower ≈ 0.81 — under-populated lower level → line
weaker in NLTE → positive Δ, correct sign). A b_lower ≈ 0.8 is on the LOW side for a +0.22 dex
correction.

## B3 — which branch fired

**Not** unit/clamp (ruled out). **Not** grid placeholder/zero (ruled out; grid thermalises). **Not**
a model-atom difference (same atom, refuted). The Engine-B departures are physical.

Remaining candidates (case ii): **(a)** a MAFAGS-OS-vs-MARCS **atmosphere** effect on Ti I (a trace,
strongly over-ionization-sensitive neutral species), and/or **(b)** a **downstream** Ti-line-specific
effect (line formation / gf / vdW of the GES Ti lines / EW-inversion). Two constraints bear on it:

- The **same deck code** reproduced 10 other elements' anchors within tol — including **O 777
  (−0.105 vs Amarsi-2019 1D −0.134, MARCS-synth vs a non-MARCS anchor, agree to 0.03)**. A generic
  downstream EW-inversion bug would have shown across elements; it did not. This argues **against** a
  generic deck bug and toward a **Ti-line/atmosphere-specific** cause.
- The +0.11-dex Engine-A↔Engine-B gap is far larger than the ~0.03–0.05 cross-atmosphere spread seen
  for the passing elements → Ti is a genuine anomaly, not noise.

**Verdict: RCA has ruled out the cheap defect classes and confirmed physical departures + same atom;
it has NARROWED the +0.11 excess to MAFAGS-OS-vs-MARCS atmosphere and/or a Ti-line-specific effect,
but has NOT fully localized it.** Ti remains **CHECK — NOT registered, NOT a "model-atom difference"**.

## Decisive next step — RESOLVED by RYA-542 (see below)

The deferred discriminator has now been run. See **"RYA-542 resolution"**. Result: **ATMOSPHERE
branch** — an independent TSFitPy MARCS synthesis reproduces the deck's ~+0.20, so +0.221 is NOT a
deck bug; the +0.11 Engine-A↔Engine-B gap is the MAFAGS-OS-vs-MARCS atmosphere systematic. Ti's
Engine-B value is legitimate (atmosphere-flagged), NOT a model-atom difference, NOT a deck defect.

---

# RYA-542 resolution — independent TSFitPy MARCS Bergemann-2011 Ti synthesis

**Date:** 2026-07-09 · **Branch:** `ryandamienschmitt/rya-542-...` (off origin/main 4ff3fb5) · **NO MERGE.**
Script: `scripts/rya542_ti_tsfitpy.py`. Full log: `rya542_tsfitpy_ti_run.log` (this dir). Runs on Sirius only.

## Method (what makes it an independent discriminator)

An **independent driver** on the **same physics inputs** as the deck. It reuses TSFitPy's own
synthesis + EW + abundance inversion (`generate_and_fit_atmosphere` → `TurboSpectrum` class +
`calculate_equivalent_width` + `root_scalar`), sharing with `scripts/ts_gerber_gate.py` **only** the
compiled Fortran (`bsyn_lu`/`babsma_lu`/`interpol_modeles_nlte`) and the physics:

- atmosphere: **MARCS** standard-composition grid (same grid the deck uses),
- model atom: **atom.ti503b (Bergemann 2011)** + the Gerber `NLTEgrid4TS_TI_MARCS_Feb-21-2022.bin`
  (54.8 GB), pulled **once** into the RYA-540 persistent md5-pinned cache (`bin_md5 5677e3a7…`) and **kept**,
- line data: the **verbatim GES Ti I rows** (loggf/χ/level-IDs) via the deck's `ges_lines()`,
- solar node: via the deck's `_solar_node()` = STAR_PARAMS **5772 / 4.438 / 0.0 / 1.0** (no hardcode),
- reference A(Ti) = **4.90** (deck a_sun; [Ti/Fe]=−0.04 vs TSFitPy solar 4.94). Validate-don't-tune.

The suspected "deck line handling" lives entirely in the Python driver — which here is **TSFitPy's,
not ours**. Departures verified to **engage** (bsyn log: `read departure file header`, `NLTE
abundance: 4.90`, 503-level departures, `Ti I NLTE … 1 lines in the interval`) — not the RYA-533
silent-LTE trap. `delta = A_NLTE − A_LTE` (same convention as the deck's +0.221 and Engine-A +0.108).

## Result — three deltas side by side

| line (Ti I) | TSFitPy-MARCS (independent) | deck Engine-B (TS-Gerber, MARCS) | EW_LTE (mA) |
|---|---|---|---|
| 5689.460 | **+0.249** | +0.266 | 14.0 |
| 5648.565 | **+0.193** | +0.207 | 14.2 |
| 5662.150 | **+0.203** | +0.221 | 26.3 |
| **median** | **+0.203** | **+0.221** | — |

Engine-A (MPIA **Bergemann-2011**, **MAFAGS-OS**) = **+0.108**.

`|median − deck +0.221| = 0.017`  ·  `|median − MAFAGS +0.108| = 0.095`.

## Decision — ATMOSPHERE branch

A **fully independent** MARCS Bergemann-2011 synthesis reproduces the deck's large correction
(+0.203 vs +0.221, per-line within ~0.02) and does **NOT** reproduce MAFAGS-OS (+0.108, off by
0.095). Two independent drivers agreeing on **MARCS+Bergemann-2011 → ~+0.20** rules out the
"deck-Ti-line-handling-bug" hypothesis.

⇒ **The deck's +0.221 is the genuine MARCS Ti I NLTE correction, not a defect.** The +0.11
Engine-A↔Engine-B gap is a real **MAFAGS-OS-vs-MARCS atmosphere systematic** (Ti I is a trace,
strongly over-ionization/atmosphere-sensitive neutral). It is **NOT** a model-atom difference (same
Bergemann-2011 atom, confirmed) and **NOT** a deck bug.

**Ti disposition:** moves from "CHECK — unresolved same-atom systematic" to **atmosphere-flagged
(legitimate value), same-atom-confirmed**. The Engine-A(MAFAGS-OS)↔Engine-B(MARCS) reconciliation is
now a **documented atmosphere systematic** to record when Ti is registered (the register/test flip is
the reconciliation follow-on — NOT done on this ticket; NO MERGE). Off the Beta-science critical path.

## Cache confirmation (RYA-540 governing rule)

Ti grid `NLTEgrid4TS_TI_MARCS_Feb-21-2022.bin` (**54.8 GB**, `bin_md5 5677e3a728c2…`, `zip_md5
c4ef399f…`) is in `gerber_ts/_cache_index.json` and **retained on disk after the run — NOT freed**
(no free-after-gate). This was the **first real production pull through the RYA-540 persistent cache**
(prior smoke was synthetic). The last Ti download, ever.

## Engineering notes (fixes made to stand up the independent path)

1. **Compiled the LTE `interpol_modeles`** binary (`gfortran -o interpol_modeles interpol_modeles.f`)
   in `Turbospectrum_NLTE/interpolator/` — only the NLTE interpolator had been built; TSFitPy's LTE
   EW step needs the LTE one. (Engine completion, not a data mutation.)
2. **TSFitPy builds the model-atom path by raw string concat** (`model_atom_path + atom_file`), so
   `model_atom_path` must end in `/` — else bsyn fails opening `gerber_tsatom.ti503b`. Wrapper sets it.
3. Built a dedicated `venv_tsfitpy` (dask) — did NOT pollute the pinned `venv_pysme`.

## Mn parallel (noted, not blocked here)

Mn carries the same same-atom question (Gerber ships **Bergemann 2019** Mn; deck ~½ vs MPIA). The
identical discriminator (independent TSFitPy MARCS Mn synthesis, or a 1D-vs-3D analog) applies. Mn
**passed** its gate so it stays registered and is **not** blocked on this — a TSFitPy Mn cross-check is
the analog to run only if Mn's ~½ is ever contested.

## Note on Mn (related, out of scope here)

Mn's RYA-534 "~½ model-atom" framing carries the **same** question (Gerber ships Bergemann's Mn atom
too). Mn PASSED its gate so it stays registered, but its "model-atom" caveat was softened to
"same-atom provenance question (RYA-411/535)". A parallel Mn RCA is warranted.
