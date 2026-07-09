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

## Decisive next step (deferred; needs a fresh grid pull — Keeper-throttled 55 GB)

Run an **independent MARCS-Bergemann-2011 Ti synthesis via TSFitPy** (the reference pipeline that
ships these exact grids): if it returns ~+0.22 → the MARCS atmosphere genuinely gives it (real
atmosphere effect, same atom; then compare to a 3rd reference e.g. Mallinson 2024 / Sitnova). If it
returns ~+0.11 → the excess is in our Engine-B deck's Ti line handling (localize gf/vdW/blend/EW).
Also map the GES line indices to the `atom.ti503b` level table to confirm the dump reads the exact
transition levels bsyn uses. Until then: Ti CHECK, RCA open, test strict-xfail, not merged.

## Note on Mn (related, out of scope here)

Mn's RYA-534 "~½ model-atom" framing carries the **same** question (Gerber ships Bergemann's Mn atom
too). Mn PASSED its gate so it stays registered, but its "model-atom" caveat was softened to
"same-atom provenance question (RYA-411/535)". A parallel Mn RCA is warranted.
