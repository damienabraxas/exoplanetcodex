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

### `vald_solar_raw.txt`
Raw VALD3 long-format output for solar parameters (Teff 5777 K, log g 4.44, Vmicro 1.0 km/s).
Used for solar calibration EW measurements. Same threshold and format as vald_55cnc_raw.txt.

### `linelist_master.csv`
Processed master line list for 55 Cnc A. 125,617 rows.

Built from `vald_55cnc_raw.txt` by `scripts/build_linelist.py`. Priority is assigned from
`data/config/elements_master.json` (27 target elements). Two NIST-only science lines (O I,
Ni I) are injected unconditionally for the O I/Ni I 6300 Å blend analysis.

### `linelist_solar.csv`
Same schema as `linelist_master.csv`, built from `vald_solar_raw.txt` at solar parameters.
Used by the solar EW pipeline (`pipeline/lines_fit.py run(star_id='solar')`).

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
| `blend_flag` | True if another line is within 0.10 Å |
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

`linelist_master.csv` and `linelist_solar.csv` are produced by `scripts/build_linelist.py`.
Run from the repo root whenever the VALD3 extract or element targets change:

```bash
# 55 Cnc A
python scripts/build_linelist.py \
    --star 55cnc \
    --vald data/linelists/vald_55cnc_raw.txt \
    --qa

# Solar calibration
python scripts/build_linelist.py \
    --star solar \
    --vald data/linelists/vald_solar_raw.txt \
    --out data/linelists/linelist_solar.csv \
    --qa
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
