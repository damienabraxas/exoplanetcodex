# data/linelists

Atomic line data for the Exoplanet Codex spectral analysis pipeline.

## Files

### `vald_55cnc_raw.txt`
Raw VALD3 long-format output for 55 Cancri A stellar parameters.
- Source: VALD3 Extract Stellar (vald.astro.uu.se)
- Teff: 5196 K · log g: 4.41 · Vmicro: 0.9 km/s
- Wavelength range: 3780–6910 Å (full HARPS range)
- Detection threshold: 0.001 (0.1% central depth minimum)
- Lines extracted: 125,615 (includes HFS components for Ba, Eu, Li, Mn)
- Retrieved: 2026-05-30

### `vald_55cnc_nir_raw.txt`
Raw VALD3 long-format output for 55 Cancri A stellar parameters, NIR extension.
- Source: VALD3 Extract Stellar (vald.astro.uu.se) — request RyanSchmitt.019389
- Teff: 5196 K · log g: 4.41 · Vmicro: 0.9 km/s (same as `vald_55cnc_raw.txt`)
- Wavelength range: 5000–30000 Å
- Lines extracted: 21,312 — not truncated (019387 covering 1150–30000 Å hit the 100K cap)
- Retrieved: 2026-06-07

### `vald_solar_raw.txt`
Raw VALD3 long-format output for solar parameters (Teff 5777 K, log g 4.44, Vmicro 1.0 km/s).
Used for solar calibration EW measurements. Same threshold and format as vald_55cnc_raw.txt.
- Wavelength range: 3780–6910 Å (optical) · Detection threshold: 0.001 · HFS: on

### Solar non-optical extensions — `vald_solar_{fuv,nearuv,redopt,ir_*}_hfson_raw.txt` (RYA-381)
Six VALD3 Extract Stellar deliveries that extend the solar list beyond the optical, at the
SAME solar parameters (Teff 5777 / log g 4.44 / Vmicro 1.0; applied [M/H]=0.00) with **HFS
on**. Gated by `scripts/intake_solar_nonoptical_rya381.py` (all 6 ACCEPT: complete span,
not truncated, solar composition, `hfs:` reference tags present).

| File | Range (Å) | VALD id | Lines | Threshold |
|------|-----------|---------|-------|-----------|
| `vald_solar_fuv_1150_2000_hfson_raw.txt`    | 1150–2000  | 019691 | 28,111 | 0.05 |
| `vald_solar_nearuv_2000_3780_hfson_raw.txt` | 2000–3780  | 019692 | 83,434 | 0.05 |
| `vald_solar_redopt_6910_9500_hfson_raw.txt` | 6910–9500  | 019693 |  1,521 | 0.05 |
| `vald_solar_ir_9500_14000_hfson_raw.txt`    | 9500–14000 | 019694 |    989 | 0.05 |
| `vald_solar_ir_14000_18000_hfson_raw.txt`   | 14000–18000| 019695 |  1,696 | 0.05 |
| `vald_solar_ir_18000_25000_hfson_raw.txt`   | 18000–25000| 019696 |  1,019 | 0.05 |

> **⚠ Detection-threshold heterogeneity (RYA-381):** these wings were extracted at central-depth
> threshold **0.05**, the optical core (`vald_solar_raw.txt`) at **0.001** — 50× shallower. The
> assembled `linelist_solar.csv` is therefore NOT depth-homogeneous; the wings carry only the
> strong lines. Weak-line / blend / n-capture work beyond the optical (e.g. the Elgueta IR
> windows — see `scripts/elgueta_ir_crosscheck_rya381.py`) needs a **0.001 re-extraction** of the
> non-optical chunks. Per-source provenance: `linelist_solar_provenance_rya381.json`.

### `linelist_master.csv`
Processed master line list for 55 Cnc A. 125,617 rows.

Built from `vald_55cnc_raw.txt` by `scripts/build_linelist.py`. Priority is assigned from
`data/config/elements_master.json` (27 target elements). Two NIST-only science lines (O I,
Ni I) are injected unconditionally for the O I/Ni I 6300 Å blend analysis.

### `linelist_full.csv`
Extended line list for 55 Cnc A covering optical + NIR. 140,483 rows.

Built from `vald_55cnc_raw.txt` (3780–6910 Å) merged with `vald_55cnc_nir_raw.txt`
(6910–30000 Å) via `scripts/build_linelist.py --vald2`. Same schema as `linelist_master.csv`.
Covers 3780.038–29994.710 Å; 97,899 target-element lines (27 elements, 74 species total).
UV extension (1150–3780 Å) pending a separate VALD request.

### `linelist_solar.csv`
Solar line list used by the solar EW pipeline (`pipeline/lines_fit.py run(star_id='solar')`).

**RYA-381: extended from optical-only (3780–6910 Å) to the full 1150–24985 Å span** by
appending the six non-optical solar deliveries above via
`scripts/assemble_solar_linelist_rya381.py`. 225,741 rows. The curated optical core
(108,971 rows — RYA-365/368 gf adjudications, NIST injections, vetted `blend_flag`) is
preserved byte-for-byte; only wing rows are added. `vald_proximity_flag` is recomputed over
the full union. See the depth-heterogeneity warning above before using the wings for
weak-line work.

#### Column schema

| Column | Description |
|--------|-------------|
| `element` | Element symbol (Fe covers both Fe I and Fe II; distinguished by `ion`) |
| `ion` | Ionization stage (I = neutral, II = singly ionized) |
| `wavelength_air_A` | Air wavelength (Å) |
| `excitation_potential_eV` | Lower-level excitation potential (eV) |
| `log_gf` | Log oscillator strength × statistical weight |
| `loggf_source` | `VALD3` or `NIST` |
| `nist_grade` | NIST ASD accuracy grade: A+ (<0.3%), A (<1%), B (<3%), C (<10%) |
| `damping_rad` | Radiative damping constant (log) |
| `damping_stark` | Stark damping constant (log) |
| `damping_vdW` | van der Waals damping constant (log) |
| `central_depth` | Predicted central depth (0–1) from VALD3 |
| `vald_proximity_flag` | Continuous proximity-contamination score 0–1 (RYA-209; VALD neighbour density × relative depth, 0.5 Å window) |
| `blend_flag` | Vetted spectroscopic exclusion (RYA-209/358): True only for literature/synthesis-confirmed non-separable blends, not raw proximity |
| `priority` | 1 = science critical, 2 = tracers (iron-peak/s-r-process/bio), 3 = supplementary |
| `notes` | Special handling notes |

### `nist_crosscheck.csv`
NIST ASD cross-validation for Tier 1 + Tier 2 science-critical lines.
Use to verify VALD3 log gf values and confirm NIST grades before abundance analysis.

Tier 1 (anchors + blend-critical): Fe I 5576/6065/6136, Fe II 6149/6247,
O I 6300.304, C I 5380, Li I 6707, P I 6034/6043

Tier 2: Na I, Mg I, Si I, S I, Ca I, Ti I, Ni I 6300.336, Ba II, Eu II

### `nist_reference.csv`
NIST A and A+ grade reference lines for pipeline QA. Use as Type B uncertainty
anchors. Highest-accuracy wavelength/log gf values available.

### `loader.py`
Pipeline module for loading `linelist_master.csv` into the analysis workflow.
See `tests/test_linelist_loader.py` for usage examples.

## Rebuilding the line lists

`linelist_master.csv`, `linelist_solar.csv`, and `linelist_full.csv` are produced by
`scripts/build_linelist.py`. Run from the repo root whenever the VALD3 extract or element
targets change:

```bash
# 55 Cnc A optical only (HARPS range)
python3 scripts/build_linelist.py \
    --star 55cnc \
    --vald data/linelists/vald_55cnc_raw.txt \
    --qa

# Solar calibration
python3 scripts/build_linelist.py \
    --star solar \
    --vald data/linelists/vald_solar_raw.txt \
    --out data/linelists/linelist_solar.csv \
    --qa

# 55 Cnc A full range (optical + NIR, 3780–30000 Å)
python3 scripts/build_linelist.py \
    --vald  data/linelists/vald_55cnc_raw.txt \
    --vald2 data/linelists/vald_55cnc_nir_raw.txt \
    --vald2-min-wave 6910.0 \
    --out   data/linelists/linelist_full.csv \
    --qa

# When UV extract arrives (1150–3780 Å), extend to full range:
# python3 scripts/build_linelist.py \
#     --vald  data/linelists/vald_55cnc_uv_raw.txt \
#     --vald2 data/linelists/vald_55cnc_nir_raw.txt \
#     --vald2-min-wave 3780.0 \
#     --out   data/linelists/linelist_full.csv \
#     --qa
```

Inputs consumed:
- `data/linelists/vald_<star>_raw.txt` — raw VALD3 long-format extract (in Git LFS)
- `data/linelists/nist_crosscheck.csv` — NIST grades for Tier 1+2 lines
- `data/config/elements_master.json` — 27 target elements with priorities

## Special handling notes

| Line | Issue | Action |
|------|-------|--------|
| O I 6300.304 | Forbidden line; Ni I 6300.336 blend | Ni contribution subtracted via COG from clean Ni I lines (Allende Prieto+2001) |
| Ni I 6300.336 | O I blend partner | Measured jointly with O I; EW predicted from COG |
| Ba II 5853 | 23 HFS components in 0.015 Å | Fit single Voigt profile = total HFS EW |
| Eu II 6645 | 30 HFS components in 0.10 Å | Narrow ±0.15 Å window; single profile = total HFS EW |
| Li I 6707 | CN blend; HFS doublet | ±0.25 Å window at doublet midpoint; flagged upper limit |
| C I 5380.337 | High excitation (7.685 eV) | Sensitive to Teff errors ±50 K |
| P I 6034, 6043 | Weak lines, depth < 1% | Require S/N > 300; NIST-injected into linelist |
| N I 7468.31 | Outside HARPS range | Near-IR only — skip for this dataset |
| Cu I 5105/5218 | Only 2 usable lines | Best effort; upper limit acceptable |
| Sr II 4077/4215 | Blue HARPS edge; crowded | Best effort; large continuum uncertainty |

## Grade definitions (NIST ASD)

| Grade | Accuracy | Pipeline use |
|-------|----------|--------------|
| A+ | < 0.3% | Yes |
| A  | < 1%   | Yes |
| B  | < 3%   | Yes (default minimum) |
| C  | < 10%  | Optional (P I lines only) |
| D  | < 25%  | Excluded |

Default minimum grade: `B` (set in `config/constants.py` → `PIPELINE['min_nist_grade']`)

## `reference_sets/` — imported reference line sets (the `line_set` provenance axis)

Line sets **published by someone else** and transcribed here so a replication can be run on
*their* selection rather than ours. The `line_set` vocabulary is closed and declared in
`pipeline.model_registry.LINE_SETS`; RYA-1111 wires these to measurements.

### `gbs_solar_fe_rya1110.csv` — Gaia FGK Benchmark Stars, solar Fe (RYA-1110)

Jofré et al. 2014, A&A **564**, A133 (`jofre2014`). `line_set=gbs`, band VIS, 159 rows —
the Sun's full published selected set, reproducing Table 3's N(Fe I)=150 / N(Fe II)=9.

* **The replication set is `rew_class == "pass"`: 142 lines** (133 Fe I + 9 Fe II),
  λ 4787.83–6820.37 Å, E<sub>low</sub> 0.11–5.10 eV, log(EW/λ) −5.934 … −4.820.
  `excluded` (14) and `ambiguous` (3) are kept in the file so the cut can be audited.
* **Two gf columns, both shipped, neither preferred:** `log_gf_gbs` (Jofré Tables 4/5, on
  121 of the 142) and `log_gf_ours` (`canonical_gf`, on all 142), plus `gf_synth_ges`.
  Choosing between them is an **open decision for Ryan** — see
  `docs/design/rya1110_gbs_reference_lineset.md`.
* **Per-line gf provenance is decoded** (`gf_source_per_line`, populated on all 159 rows)
  from two sources whose ordering matters: Heiter+2021's per-line `r_loggf`
  (`data/reference/heiter2021_ges/`, VizieR `J/A+A/645/A106`) where its log gf **equals**
  the published GBS value, and Jofré's own Table 4/5 footnote where it does not — because
  Heiter is GES **v6** and Jofré used **v3**, so a revised value's source is not the GBS
  value's source. `gf_source_basis` says which route answered.
* 🔴 **`gf_source_firewalled` flags three Fe II lines** (5414.07, 5425.26, 6432.68) whose
  GBS gf is Meléndez & Barbuy 2009 — RYA-161-firewalled, partly solar-fitted. On two of
  them our own adopted value is the identical number.
* 🔴 **The paper's stated log(EW/λ) ≤ −4.8 cut does not select the set it publishes** —
  14 of the Sun's 159 published lines violate it on every method (713 of 4252 across all
  34 benchmark stars). Both facts are carried; neither is overridden.
* Built by `scripts/rya1110_build_gbs_fe_lineset.py` from the pinned holding
  `data/reference/jofre2014_gbs/`; guarded by `tests/test_gbs_lineset_rya1110.py`, which
  rebuilds and compares byte-for-byte.

### `gbs_solar_fe_coverage_rya1110.csv`

The 142 selected lines against every holding in `measure_band_ew._INSTRUMENT_HOLDINGS`.
All five HARPS/Kitt Peak holdings reach all 142; the two IAG holdings split them 132 + 10;
no line falls in a registered telluric band (the bluest is O₂ B at 6867 Å, the reddest GBS
line is 6820.37 Å).
