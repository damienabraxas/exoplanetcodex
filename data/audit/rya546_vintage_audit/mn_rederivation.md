# RYA-546 Mn re-derivation — triplet-exact ab-initio δ → solar Mn PASS (closes RYA-411)

**Date:** 2026-07-13 · Branch `ryandamienschmitt/rya-546-...` · **MERGES** (the re-derivation; the audit parts stay no-merge).
Reused RYA-473 HFS machinery + RYA-476 grid wiring; no new synthesis path. Runs on Sirius.

## The reversal

RYA-411 measured the Amarsi-2020 vs Bergemann Mn difference (~2×) and **kept Bergemann +0.107**,
trusting "the dedicated Mn-NLTE literature (Bergemann & Gehren 2008)." **RYA-546's H-collision
vintage audit reverses that reasoning:** B&Gehren-2008 *is* the outdated **scaled-Drawin** inelastic-H
recipe, which inflates the over-ionization correction (the exact Ti defect, RYA-542/544). The
Amarsi-2020 GALAH grid uses **ab-initio (Barklem/Grumer)** H collisions → the smaller δ is *correct*,
not "survey-grade under-correction."

## Gate 0 — Amarsi Mn grid live + non-NaN triplet δ ✓ PASS

The Den Hartog e6S→z6P triplet (6013/16/21) NaN'd on the MPIA high-EP grid (why RYA-476 blocked). On
the **live** `nlte_Mn_scatt_pysme.grd` (Sirius), HFS-resolved identically on both legs:

| Mn I line | live Amarsi HFS δ |
|---|---|
| 6013.510 | +0.0296 |
| 6016.670 | +0.0242 |
| 6021.820 | +0.0171 |
| **median** | **+0.0242** |

**NON-NaN → Gate 0 PASS.** ≪ Engine-A +0.108 (scaled-Drawin); corroborated by Engine-B TS-Gerber
+0.043 (RYA-534). Solar node from `get_star_params('solar')` = 5772/4.438/0.0/1.0 (no hardcode).
In-hull guard fired: grid coverage Teff(2500,8000) logg(−0.5,5.5) [Fe/H](−5,1); solar in-hull.

## Additions (all landed)

- **A — HFS symmetry:** `nlte_delta` derives δ = A_NLTE − A_LTE with the same HFS-split triplet on BOTH
  legs (the harness `_synth_ew` expands the DLSSC hyperfine components for the NLTE EW and the LTE COG
  alike). Apples-to-apples with the +0.043 Engine-B corroboration.
- **B — hard in-hull guard (`pysme_nlte.assert_in_grid_hull`, wired into `nlte_delta`):** hard-fails +
  logs coverage if (Teff,logg,[Fe/H]) is outside the grid box (rectangular ⇒ box = hull). No silent
  extrapolation. Solar in-hull; lands now so 55 Cnc A / α Cen A target runs can't extrapolate.
- **C — solar-anchored per-line report:** triplet δ printed per-line (+0.030/+0.024/+0.017 — varies,
  NOT flat, so no grid/HFS red flag); gate against STAR_PARAMS Asplund-2021 A(Mn)☉.

## Solar gate — A(Mn)☉ ✓ PASS (CURATION-OWED → PASS)

Re-ran `measure_mn_hfs_synthesis_rya473` (HFS flux fit for A_LTE, absolute anchor from the HFS-synth
path NOT the TS synth-EW path; RYA-546 Part B +0.25 offset is irrelevant — Δ=NLTE−LTE within one code
cancels it):

| Mn I line | A(Mn)_LTE | χ²ᵣ | σfit |
|---|---|---|---|
| 6013.51 | 5.447 | 120.6 | 0.006 |
| 6016.67 | 5.177 | 126.9 | 0.011 |
| 6021.82 | 5.442 | 119.4 | 0.012 |
| **median** | **5.442** (mean 5.355, σ 0.154, n=3) | | |

**A(Mn)_LTE 5.442 + live Amarsi Δ +0.024 → A(Mn)_NLTE = 5.466 → +0.046 vs Asplund 5.42.**
`phase_c_verdict_rya371`: **Mn 5.466 +0.046 σ0.15 n=3 → PASS** (HFS synthesis, triplet). Was 5.554
(+0.13) on the +0.107 high-EP mismatch. Verdict **CURATION-OWED → PASS**; the RYA-473 δ-provenance
caveat is resolved; **RYA-411 closed.**

## Files
- `pipeline/pysme_nlte.py` — Addition B `assert_in_grid_hull` + `nlte_delta` guard (MERGES).
- `scripts/rya546_mn_gate0.py` — Gate 0 + triplet-exact δ driver.
- `config/constants.py` — Mn comment reversed (RYA-411→RYA-546); +0.107 marked SUPERSEDED for solar.
- `data/audit/mn_hfs_synthesis/solar_mn_hfs_synthesis_rya473.{json,csv}` — regenerated with the live δ.

## Scope / follow-on
- **Solar-only this pass** (per brief). Fan-out (55 Cnc A / α Cen A) deferred — but Addition B (the
  in-hull guard) is landed so those runs hard-fail rather than extrapolate.
- **RYA-527 fan-out:** repoint the EW-path Mn NLTE grid (`Mn_Bergemann_MPIA.csv` +0.107) to an Amarsi
  ab-initio Mn δ-grid for non-solar stars; fold the corrected A(Mn)☉ 5.466 into the gold v3 re-freeze
  (Mn gold was 5.554 → now 5.466).
- **Env note:** venv_pysme needed h5py/hdf5plugin/lockfile/dill/statsmodels to import iSpec (added);
  run with `venv_pysme + ISPEC_DIR=PYTHONPATH=/srv/codex/engines/ispec_src` for both iSpec+PySME.

## Sources
- Amarsi et al. 2020 (A&A 642 A62), Zenodo 3982506 — GALAH ab-initio Mn departure grid.
- Den Hartog, Lawler, Sobeck, Sneden, Cowan 2011 (ApJS 194, 35) — the e6S→z6P triplet HFS + gf.
- Bergemann & Gehren 2008 (A&A 492, 823) — the scaled-Drawin MPIA Mn atom (superseded for solar).
- Grumer & Barklem 2020 (A&A 637 A68) — ab-initio Mn+H collisions.
