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

So the remaining work is: **(1)** a reader for the Amarsi grid, **(2)** feeding the
departures to Turbospectrum (the Gerber-2022 adaptation, or iSpec's dep-grid HDF5),
**(3)** the b-factor → delta extraction. **(1)** the `PySMEGrid` reader is **built
and verified against the real vendored Na grid**; **(3)** is built and unit-tested
in `pipeline/nlte_bfactor_synth.py`; **(2)** plus the live synthesis are wired and
**fail loud** until hooked to the Gerber-2022 TS departure machinery (no silent LTE
fallback — the RYA-289 anti-pattern).

## The grid format — PySME DirectAccessFile (`.grd`), reader done

Intake established the real on-disk format: a PySME **"DirectAccess file Version
1.10"** binary (the legacy IDL/SME save format) — *not* the ASCII layout the Step-1
scaffold assumed. Structure: a 64-byte version string; a `(nblocks u8, dir_length
i2, ndir u8)` header; `ndir` directory entries `(key S256, size[23] i4, pointer
i8)`; arrays memmapped at `pointer` with the IDL (column-major) shape reversed for
NumPy. Keys: axis vectors `teff`/`grav`/`feh`/`abund`; model-atom data
`conf`/`term`/`spec`/`J`/`energy`; and one departure block **per model**, keyed
`t{teff}_g{logg}_m{feh}_a{abund}`, each an `(ndepth, nlevel)` array of
`b = n_NLTE/n_LTE`.

`pipeline.nlte_bfactor_synth.PySMEGrid` reads it (memmapped — never loads the
multi-GB file into RAM). **Format check (the physics):** at the solar node the deep
layers give `b → 1.000` (thermalised to LTE) while the surface departs — verified
live on the Na grid, and asserted in the tests.

## Validate-don't-tune (the critical Step 3)

Before trusting the path for Al/K, reproduce an **already-known** correction: run
the vendored Na grid (INSPECT, δ = −0.107) through the synthesis path and confirm
it lands within tolerance. **If it cannot, STOP** and do RCA — never fit to Asplund.
`--validate-against Na` is that gate (loud-stops until the Gerber-2022 TS wiring is
in).

## Data — intaken (md5-verified, gitignored, provenance committed)

The grids live in `data/nlte_grids/amarsi_galah/` (the rolling Amarsi *"Grid/NLTE"*
Zenodo deposit). Each `.grd` is multi-GB → **gitignored**; only the kB axis/level
files + a provenance JSON (source DOI + Zenodo version + md5 + coverage + refetch)
are committed. Vendored so far:

| El | Zenodo | record | models × levels | [Fe/H] |
|----|--------|--------|-----------------|--------|
| Cu | v6 (Caliskan 2024) | 15062813 | 36855 × 150 | −5…+1 |
| S  | v7 (Amarsi 2025, A&A 703 A35) | 17064337 | 36855 × 181 | −5…+1 |
| Na | v3 (Amarsi 2020, A&A 642 A62) | 3982506 | 43680 × 140 | −5…+1 |
| Al | v3 | 3982506 | 43680 × 118 | −5…+1 |
| K  | v3 | 3982506 | 43680 × 189 | −5…+1 |

All cover [Fe/H] up to +1 (no metal-rich clamp → 55 Cnc). `atmos_*.txt` is shared
across elements (same MARCS grid; byte-identical md5).

## Resume point (grids + reader done)

1. **The TS departure adaptation (the remaining engineering):** interpolate the
   `PySMEGrid` departures to the target (Teff, logg, [Fe/H], abund) and write them as
   a Turbospectrum departure file + `NLTEINFOFILE` + model atom (Gerber 2022, A&A
   666 A18), *or* convert to iSpec's `input/dep-grid/{El}_nlte_grid_data.h5`. Wire
   `synth_ew_nlte_vs_lte` to run `babsma_lu` + `bsyn_lu` NLTE vs LTE.
2. `--validate-against Na` must reproduce −0.107 ± 0.03 → machinery trusted (the
   adaptation is exactly where a silent unit/format error hides; this is the guard).
3. Derive Cu / S / Al / K, register source-aware in `NLTE_CORRECTION_ELEMENTS`,
   reproduce the solar anchors, flip them in `physics_regime_rya400.yaml`
   (GET-GRID → handled, synthesis provenance) keeping the RYA-400 audit PASS. K:
   confirm the 7665/7699 resonance doublet is telluric-recoverable (RYA-380).
