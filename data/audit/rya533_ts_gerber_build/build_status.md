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

## Step 3 — THE GATE: DONE, GREEN (validate-don't-tune)

The full deck is wired and the Na gate reproduces the anchor. Steps executed (all Sirius):
1. Unzipped the Na grid → `.bin` (8.4 GB unpacked).
2. Provisioned MARCS standard-composition models (141 MB, 15691 models) from the same MPG Keeper share; solar node `p5750_g+4.5_z+0.00`.
3. Compiled + ran `interpol_modeles_nlte` → `sun_Na_coef.dat` departure file (4D interp; **gotcha: it reads the aux ROW COUNT from stdin — pass wc-l-1 = 436255, not "9"**).
4. Na line list = the 5682/5688 lines taken VERBATIM from the bundled GES NLTE line list (carries the level IDs; **the silent-error site — a plain VALD line makes bsyn fall back to departure=1**).
5. Deck driver `scripts/ts_gerber_na_gate_rya533.py`: babsma (contopac, `DATA/` symlinked) → bsyn NLTE (NLTEINFOFILE with exact `# path for ...` markers) + bsyn LTE COG → per-line EW inversion → Δ. RAISES if departures don't engage.

**RESULT (`na_gate_result.txt`):** Na I 5682.633 δ −0.065 (EW 112.5 mA) / 5688.205 δ −0.072 (EW 147.8 mA) → **median −0.068 vs anchor −0.107 ± 0.05 → PASS**. Departures engaged (read_departure 56 depths × 290 levels; no silent-LTE fallback). Cross-engine Na: INSPECT −0.107 / PySME −0.129 (RYA-529) / TS-Gerber −0.068 = the RYA-525 model-atom systematic.

Gate codified in `tests/test_ts_gerber_nlte_rya533.py` (recorded-result test always runs; live test runs on Sirius). Register **v9** Engine-B TS-native NLTE flipped to **SETTLED (Na)**; Engine-B map updated.

## Remaining (mechanical, per-element — beyond this ticket's Na gate)

The other 10 Family-A elements (O/Mg/Si/Ca/Ti/Mn/Co/Ni/Sr/Ba) are a repeat: download the element's Gerber grid + atom from the Keeper (`download_nlte_grids.py 1D <El>`), run the interpolator, swap the line list block, re-run the gate. Deck + engine are done.
