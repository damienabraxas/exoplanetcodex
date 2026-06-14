# RYA-303 — STIS UV de-dup + re-split (α Cen A / B)

Conditioning of the HST/STIS UV set surfaced by the RYA-301 audit, **before**
any UV-inclusive run (RYA-302). Source tree is read-only — nothing was deleted
or moved (DATA_SAFETY.md). Regenerate with:

    python3 scripts/dedup_resplit_stis_uv_rya303.py

## Step 1 — de-dup (SHA-256)

| metric | value |
|---|---|
| STIS UV FITS on disk (both folders) | 222 |
| **Unique frames** | **111** |
| Duplicate rows collapsed | 111 |
| Unique frames present in BOTH "Alpha Centauri A" & "B" folders | 111 (all) |

The two star folders are **byte-identical copies of the same mixed 111-frame
set** — folder name carries no information about the true star. Full hash →
path mapping: `stis_uv_dedup_map.csv`.

## Step 2 — re-split by OBJECT (not folder)

Attribution is on the **base HD number in `TARGNAME`** (the OBJECT/target
header): `HD128620`→A, `HD128621`→B. STIS visit/pointing-revision suffixes
(`-1`, `-2`, `-NEW`, `-COPY`) are **not** different stars and are ignored.

| component | unique UV frames |
|---|---|
| α Cen A (HD 128620) | **48** |
| α Cen B (HD 128621) | **63** |
| unresolved | 0 |

Matches the RYA-301 audit split (A 48 / B 63). Per-frame inventory:
`stis_uv_unique_inventory.csv`. Per-star/per-mode loader file lists:
`alpha_cen_{A,B}_{e140m,e140h,e230h,e230m}_files.txt`.

**Coordinate-named frames:** 0 — every frame carries an explicit HD identifier,
so OBJECT alone resolves attribution; the PM-cross-match fallback was not
triggered.

**Proper-motion note (documented limitation):** α Cen A/B is a tight
common-proper-motion binary whose internal separation is dominated by the
~80-yr orbit. Linear-PM propagation of catalog J2000 positions cannot model the
orbit, so an *absolute* "nearest component" PM verdict is degenerate (collapses
to one component) and is **not** used to override OBJECT. Instead the HD labels
are corroborated geometrically and orbit-independently: each frame's
`RA_TARG/DEC_TARG` vs the nearest-in-time *opposite*-component frame is
separated by a median ≈5″ (two distinct pointing clusters that drift together
under the system's large PM across 1999→2018) — i.e. the labels track two
physically distinct targets.

## Step 3 — product type, modes, coverage

All 111 are extracted 1D echelle science products (`_x1d.fits`) — no raw/2D.

| star | mode | band | frames | λ coverage (Å, **vacuum**) |
|---|---|---|---|---|
| A | E140M | FUV | 25 | 1140–1730 |
| A | E140H | FUV | 3  | 1141–1688 |
| A | E230H | NUV | 20 | 1629–3159 |
| B | E140M | FUV | 17 | 1140–1730 |
| B | E140H | FUV | 19 | 1196–1587 |
| B | E230H | NUV | 23 | 2327–3000 |
| B | E230M | NUV | 4  | 1607–2511 |

⚠️ **Wavelengths are VACUUM.** The pipeline works in **air Å for λ ≥ 2000 Å**
(VALD air/vacuum convention). The STIS loader **must apply vacuum→air
conversion** at the spectrum-matching boundary. The `wavelength_frame` column in
the inventory is set to `VACUUM` as a standing flag.

## Files in this directory

- `stis_uv_dedup_map.csv` — every on-disk frame: sha256, folder, path, size.
- `stis_uv_unique_inventory.csv` — 111 unique frames with attribution, mode,
  coverage, PM diagnostics, flags.
- `alpha_cen_{A,B}_{mode}_files.txt` — per-star/per-mode absolute-path lists for
  the loader.
