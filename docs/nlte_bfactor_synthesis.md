# Departure-coefficient (b-factor) NLTE synthesis (RYA-402)

The generalisation of the RYA-396 NLTE machinery from *"apply a pre-computed
delta"* to *"compute the delta from departure coefficients"*. Some elements
(Al, K) are published only as **departure coefficients (b = n_NLTE/n_LTE)**, not
as abundance-correction deltas — so they need NLTE *synthesis* to convert. Unlocks
Al/K (RYA-401 follow-up); complements the RYA-399 3D *delta*-application path
(this is the b-factor-only path; available to 399 wherever 3D data is published
only as departure coefficients).

```
python -m pipeline.nlte_bfactor_synth --status              # engine + data readiness
python -m pipeline.nlte_bfactor_synth --validate-against Na # reproduce INSPECT Na -0.107 (STOP-gate)
python -m pipeline.nlte_bfactor_synth --elements Al,K       # derive (after Na validates)
```

## Engine finding (verified on disk — corrects the ticket premise)

We do **not** need to build a Turbospectrum NLTE engine from scratch:

- **`bsyn_lu` already has the NLTE path compiled in** — `strings` shows
  `read_departure.f` / `read_nlteinfofile.f`, the control keys `NLTEINFOFILE` /
  `MODELATOMFILE` / `DEPARTUREFILE`, and *"no departure file → coefficients = 1.0"*
  (so unity departures reproduce LTE — the engine self-test).
- **iSpec already interpolates departure grids and feeds bsyn** —
  `ispec.atmospheres.interpolate_nlte_departure_coefficients`, with grids read from
  `$ISPEC_DIR/input/dep-grid/{El}_nlte_grid_data.h5` (HDF5, lz4).

So the remaining work is: **(1)** a reader for the Amarsi ASCII b-factor format,
**(2)** a converter to iSpec's dep-grid HDF5 (same MARCS nodes) *or* a direct
`bsyn_lu` driver, **(3)** the b-factor → delta extraction. **(1)** and **(3)** are
built and unit-tested in `pipeline/nlte_bfactor_synth.py`; **(2)** plus the live
synthesis are wired to the verified engine and **fail loud** until the grids are
installed (no silent LTE fallback — the RYA-289 anti-pattern).

## Validate-don't-tune (the critical Step 3)

Before trusting the path for Al/K, reproduce an **already-known** correction: run
the vendored Na grid (INSPECT, δ = −0.107) through the synthesis path and confirm
it lands within tolerance. **If it cannot, STOP** and do RCA — never fit to
Asplund. `--validate-against Na` is that gate; it loud-stops today because the Na
grid isn't fetched yet.

## Data — Ryan fetches (Step 1, pinned)

`data/nlte_grids/amarsi2020_galah/SOURCE_rya402.json` pins the canonical source:
**Amarsi et al. 2020, A&A 642, A62** (DOI 10.1051/0004-6361/202038650; arXiv
2008.09582), departure-coefficient grids on the Zenodo *"Grid/NLTE"* deposit
(Amarsi, Uppsala; v6 = `zenodo.org/records/15062813`). Fetch the **Na** (validation)
plus **Al** and **K** grids — `label_{El}.txt` + `nlte_{El}_*.txt` (ASCII; optional
`*_pysme.grd`) — into that folder. Na first.

## Resume point (after the fetch)

1. `read_amarsi_grid(El)` → convert to `input/dep-grid/{El}_nlte_grid_data.h5`
   (match the iSpec MARCS pack nodes; assert exact parameter coverage).
2. Wire `synth_ew_nlte_vs_lte` to the live Turbospectrum NLTE run.
3. `--validate-against Na` must reproduce −0.107 ± 0.03 → machinery trusted.
4. Derive Al/K, register source-aware in `NLTE_CORRECTION_ELEMENTS`, reproduce the
   solar anchors, flip Al/K in `physics_regime_rya400.yaml` (GET-GRID → handled,
   synthesis provenance) and keep the RYA-400 audit PASS. K: confirm the 7665/7699
   resonance doublet is telluric-recoverable (RYA-380) for the correction to matter.
