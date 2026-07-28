# RYA-529 — RYA-402 b-factor NLTE synthesis re-provisioned + re-validated on Sirius

**Date:** 2026-07-06 · **Branch:** `ryandamienschmitt/rya-529-...` (off `origin/main` 21d6e6c) · **No merge — Ryan reviews.**

This is the Engine-B-substrate enabler for the two-engine floor (RYA-525), unblocking the
RYA-527 re-freeze. It re-provisions + re-validates the RYA-402 b-factor NLTE-correction
capability on the Sirius-authoritative compute stack (post RYA-506/511/517 migration).

## Step 0 — runtime spec recovered from merged RYA-402 (not from memory)

The merged capability's **validated engine is PySME**, not Turbospectrum `bsyn`:

- `pipeline/nlte_bfactor_synth.py::synth_ew_nlte_vs_lte` (the TS-`bsyn` NLTE deck) is a
  documented `NotImplementedError` — the SME `.grd` grids carry departures + levels but no
  TS model atom (transitions/collisions), and Ryan's 2026-06-21 RYA-402 decision resolved
  Family B (the b-factor grids) onto **PySME as the production Δ-engine** (the native `.grd`
  consumer), keeping Turbospectrum as the LTE spectral-synthesis engine and applying the
  derived Δ = A(NLTE) − A(LTE) to the TS/EW LTE baseline.
- The Na validation command is **`python -m pipeline.pysme_nlte Na`** →
  `pysme_nlte.validate('Na')`, PASS iff `|median Δ − (−0.107)| ≤ 0.03`
  (`_ANCHOR['Na'] = (-0.107, 0.03, 'Lind et al. 2011 INSPECT (Na I 5682/5688)')`).
- Inputs: the Amarsi PySME `.grd` departure grid (`nlte_Na_scatt_pysme.grd`) in
  `data/nlte_grids/amarsi_galah/`, the MARCS `marcs2012.sav` atmosphere (PySME auto-fetches
  from `sme.astro.uu.se/atmos`), plus committed `atmos_*/label_*` axis files + provenance JSON.

## Step 1 — Sirius inventory (gap table)

| Component | On Sirius? | Detail / source |
|---|---|---|
| PySME + venv | ✅ | `/mnt/codex-data/venv_pysme` — py3.12.13, numpy 2.2.6 (`np.trapz` present), scipy 1.18.0, pandas 2.3.3, **PySME 1.0.2** (RYA-526). venv312 reference stack untouched. |
| Amarsi-2020 b-factor grids | ✅ | Base-13 in `/mnt/codex-data/grids/nlte/amarsi_galah/` (Zenodo 3982506; fetched by RYA-477/526): H, Li, C, N, O, Na, Mg, Al, Si, K, Ca, Mn, Ba. Na grid = `nlte_Na_scatt_pysme.grd` (2.76 GB), reads 43680 models × 140 levels. |
| Na grid provenance | ✅ | `nlte_Na_scatt_pysme.tar.gz.prov.json` on Sirius (tar.gz md5 `b4b1408dc8f296fa3d14ad2e2a62452f`, Zenodo 3982506) + repo-committed `Na_amarsi2020_v3.prov.json` (RYA-402). |
| MARCS `marcs2012.sav` | ✅ (auto-fetch) | Not pre-cached; PySME downloaded it from `sme.astro.uu.se/atmos` during the gate run (Sirius has internet — same path RYA-477/526 used for Zenodo). |
| pipeline code | ✅ | Fresh worktree of `origin/main` @ `21d6e6c` at `/mnt/codex-data/codex/rya529`; Na grid symlinked in. |
| Cu (Caliskan) / S (Amarsi-2025) grids | ❌ | Not on Sirius. RYA-402 vendored them on the Mac; NOT re-provisioned here. Acquire on Sirius if the two-engine floor needs Cu/S NLTE-synth. |

**Owed → acquired this ticket:** nothing had to be downloaded — the base-13 grids + venv_pysme
were already on Sirius from RYA-477/526. Provenance already recorded. **All prior downloads/
extractions ran on Sirius, never the Mac** (RYA-526 hard rule; no `.grd` touched the Mac here).

## Step 3 + 4 — smoke + THE GATE (one command covers both)

Run on Sirius from the `rya529` worktree with `venv_pysme`:

```
$ /mnt/codex-data/venv_pysme/bin/python -m pipeline.pysme_nlte Na
{
  "element": "Na",
  "per_line": {
    "5682.633": -0.12138596510804334,
    "5688.205": -0.1373110032239575
  },
  "delta_median": -0.12934848416600042,
  "anchor": -0.107,
  "anchor_tol": 0.03,
  "passed": true,
  "ref": "Lind et al. 2011 INSPECT (Na I 5682/5688)"
}
```

- **Step 3 (NLTE engaged, no silent LTE):** per-line Δ = −0.121 / −0.137 — nonzero + finite.
- **Step 4 (THE GATE):** median Δ **−0.12935** vs anchor **−0.107 ± 0.03 → PASS**
  (|−0.12935 − (−0.107)| = 0.0223 ≤ 0.03). Validate-don't-tune: the anchor is reproduced,
  never fitted. Matches the Mac RYA-402 result (−0.121 / −0.138 / −0.129) to ≤0.001/line —
  the machinery is **stack- and platform-robust** (Mac darwin-arm64 ↔ Sirius linux-x86_64).

## Step 5 — register

`CODEX_STATE_REGISTER.md` v6 → **v7**: added the `TS NLTE synthesis (Engine B substrate)` row
to the NATIVE table (immediately after the two-engine-floor row) + changelog line. The row
truthfully attributes the correction engine to **PySME** (TS-`bsyn` deck = NotImplementedError)
and cites Zenodo **3982506** for the base-13 grids. Section-edit only; no pre-existing row
altered (Solar 5 PASS verdict intact); file grew 208 → 210 lines.
