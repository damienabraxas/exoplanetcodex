# RYA-533 — TS-native Gerber NLTE synthesis deck (Engine-B NLTE): forensic verdict + build status

**Date:** 2026-07-06 · **Branch:** `ryandamienschmitt/rya-533-...` (stacked on RYA-531 `e2c4689` → RYA-529 `dcbc741`) · **No merge — Ryan reviews.**

## Step 0 — forensic verdict: **NOT-BUILT** (evidence on all 5 channels)

| Channel | Finding | Verdict |
|---|---|---|
| 0a test tell | No TS-Gerber reproduction test in `tests/`. Only hit is a docstring in `test_nlte_bfactor_synth_rya402.py:10` saying live TS synth is "gated on the Gerber-2022 TS atom" (= not built). That test exercises only the PySMEGrid reader. | NOT-BUILT |
| 0b git archaeology | No surviving `build/rya-402-bfactor-nlte-synth` branch on origin (402's capability merged via the PySME line — PR#50 rya-409, PR#95 rya-466, commits RYA-410/411/466). No `rya402` worktree on Mac. No 402/gerber/nlte-synth branch on any Sirius checkout. | NOT-BUILT |
| 0c merged-tree | `nlte_bfactor_synth.synth_ew_nlte_vs_lte` = `NotImplementedError`. No driver writes a NLTEINFOFILE/DEPARTUREFILE, does 4D departure interpolation, or invokes babsma+bsyn in NLTE. Production NLTE is 100% delta-CSV (RYA-531). | NOT-BUILT |
| 0d grid/atom store | `/mnt/codex-data/grids/nlte/gerber_ts/` was **empty**; no Turbospectrum_NLTE build on Sirius (engines/ had only ispec/nirps/pysme); zero model_atom/departure/coef artifacts; Mac canonical store empty. | NOT-BUILT |
| 0e registration + Linear | No Family-A element carries a Gerber/TS-native flag in `NLTE_CORRECTION_ELEMENTS` (all PySME/MPIA/Korotin/INSPECT/nlte_cno). **RYA-361** ("Gerber 2023 in-synthesis TS-NLTE — install + assess", RYA-364 epic) is the pre-existing scoping ticket — **Backlog, never started, no branch, scoping-only comments.** | NOT-BUILT (scoped, never executed) |

The TS binary's OWN NLTE hooks have existed since RYA-402 (`read_departure.f`, control keys, unity-departure LTE self-test) — but **hooks ≠ our deck.** Our deck (Gerber grids on Sirius + babsma/bsyn NLTE deck + 4D departure interpolation + a passing Na gate) was never stood up.

**Provenance resolution (2022 vs 2023):** same paper — arXiv 2206.00967 (2022) = Gerber et al. **2023, A&A 669, A43**. The "A&A 666 A18" locator in RYA-531/533 is wrong; RYA-361 already cites A&A 669 A43.

## Step 2 — build progress THIS session (Sirius-only; nothing on the Mac)

**DONE:**
1. **Turbospectrum_NLTE compiled on Sirius.** Cloned `github.com/bertrandplez/Turbospectrum_NLTE` → `/mnt/codex-data/engines/Turbospectrum_NLTE`; `make` in `exec-gf/` (gfortran 15.2.0, `-fconvert=big-endian -O3 -fno-automatic -std=legacy -mcmodel=small`) built `babsma_lu` + `bsyn_lu` (**BSYN v20.1**). `strings bsyn_lu` confirms NLTE mode live: `NLTEINFOFILE`/`MODELATOMFILE`/`DEPARTUREFILE`, "asking for NLTE calculation", "missing NLTEINFOFILE!", "continuing with departure coefficient = 1." (unity-departure = LTE self-test).
2. **Deck format recovered from the repo (not memory).** `COM/script-NLTE-multi-element.com` + `DATA/SPECIES_LTE_NLTE.dat` (the NLTEINFOFILE) + `interpolator/interpol_modeles_nlte.f` (`interpol-nlte.script`). Key facts:
   - NLTEINFOFILE row: `Z  'El'  'nlte'/'lte'  'model_atom_file'  'departure_file'  'ascii'/'binary'`, preceded by a model-atom-path line + a departure-path line.
   - `bsyn_lu` `ABFIND` (EW iteration) is **NOT implemented in v20** → the NLTE−LTE Δ must be obtained by synthesizing NLTE vs LTE and integrating EWs on an LTE COG (same inversion as `pipeline/pysme_nlte`).
   - The departure file per star is produced by `interpol_modeles_nlte.f` (4D interpolation of the grid over Teff/logg/[Fe/H]/A(X) — the aux file supplies those axes + pointers).
3. **Na Gerber-2023 data provisioning (Sirius `/mnt/codex-data/grids/nlte/gerber_ts/`):**
   - `atom.na_qmh` (model atom, 2.51 MB, md5 `8297b392587941eaece3dfc9626f23be`) ✅
   - `auxData_Na_MARCS_Jul-14-2023.dat` (91.9 MB param index, md5 `bea18709a21e9bbfe144436cb435ac88`) ✅
   - `NLTEgrid4TS_Na_MARCS_Jul-14-2023.bin.zip` (1D departure grid, **15.9 GB**, md5 `d1e8b51efd66ad079bdef3377ce164d1`) ✅ downloaded.
   - Provenance JSON committed: `data/nlte_grids/gerber_ts/Na_gerber2023.prov.json`. URLs pinned from TSFitPy `nlte_grids_links.cfg`.

## Resume point (multi-session build — as the ticket anticipates)

Engine + data are being provisioned; the **deck driver + gate are not yet written** (so the register stays Engine-B-OWED; no SETTLED claim — the Na gate has NOT run). To finish:

1. Unzip the Na grid `.bin.zip` (15.9 GB, md5 `d1e8b51efd66ad079bdef3377ce164d1`, already on Sirius) → `.bin` for the interpolator.
2. Get a **solar MARCS `.mod`** atmosphere (from iSpec's MARCS pack or the TS interpolator) at Teff 5772/logg 4.44/[Fe/H] 0.
3. Run `interpol_modeles_nlte.f` to produce the solar Na **departure file** from the grid (4D interp).
4. Build a **Na line list** in TS format for 5682.633 / 5688.205 (the anchor doublet; VALD long / GES format).
5. Write the pipeline **deck driver**: babsma_lu (contopac) → bsyn_lu twice (NLTE-with-departures via NLTEINFOFILE, and LTE), integrate EWs, invert on the LTE COG → Δ = A(NLTE) − A(LTE). No silent LTE fallback (RAISE if departures not engaged).
6. **Na gate:** reproduce the anchor (~−0.107; Amarsi-2020/Lind-2011; allow 0.02–0.05 dex model-atom difference — path-correctness, not a tune). Spot-check one more Family-A element (e.g. Ca via `atom.ca105b`).
7. **Write the test** that codifies the gate (closes the RYA-530 loop — the next forensic sweep must find it).
8. Register **v8 → v9**: flip Family-A Engine-B-NLTE OWED → SETTLED (only on a green gate); update `data/audit/rya531_engineB_map/engineB_per_element.md`.

**STOP discipline:** do NOT register SETTLED until the Na gate passes; a miss → RCA. TS NLTE not engaged (0.0/LTE) → RAISE.
