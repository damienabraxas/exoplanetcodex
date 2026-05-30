# data/linelists

Atomic line data for the Exoplanet Codex spectral analysis pipeline.

## Files

### `vald_55cnc_raw.txt`
Raw VALD3 long-format output for 55 Cancri A stellar parameters.
- Source: VALD3 Extract Stellar (vald.astro.uu.se)
- Teff: 5196 K · log g: 4.41 · Vmicro: 0.9 km/s
- Wavelength range: 3780–6910 Å (full HARPS range)
- Detection threshold: 0.01 (1% central depth minimum)
- Lines extracted: 17,926
- Retrieved: 2026-05-29

### `linelist_master.csv`
Processed master line list used by the pipeline. 17,931 rows.

Built from `vald_55cnc_raw.txt` plus 5 critical lines injected from NIST ASD
that fall below the 1% VALD depth threshold but are required for science:
- O I 6300.304 Å (forbidden line — primary oxygen indicator)
- Ni I 6300.336 Å (O I blend partner — must model jointly)
- C I 5380.337 Å (high excitation carbon line)
- P I 6034.04 Å (phosphorus tracer)
- P I 6043.12 Å (phosphorus tracer)

#### Column schema

| Column | Description |
|--------|-------------|
| `element` | Element symbol |
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
| `blend_flag` | True if flagged as blend |
| `priority` | 1 = science critical, 2 = alpha elements, 3 = iron-peak / bio-elements |
| `notes` | Special handling notes |

### `nist_crosscheck.csv`
NIST ASD cross-validation for the 20 Tier 1 + Tier 2 science-critical lines
defined in RYA-64. Use to verify VALD3 log gf values and confirm NIST grades
before abundance analysis.

Tier 1 (anchors + blend-critical): Fe I 5576/6065/6136, Fe II 6149/6248,
O I 6300.304, C I 5380, Li I 6707, P I 6034/6043

Tier 2: Na I, Mg I, Si I, S I, Ca I, Ti I, Ni I 6300.336, Ba II, Eu II

### `nist_reference.csv`
NIST A and A+ grade reference lines for pipeline QA. Use as Type B uncertainty
anchors. Highest-accuracy wavelength/log gf values available.

### `loader.py`
Pipeline module for loading `linelist_master.csv` into the analysis workflow.
See `tests/test_linelist_loader.py` for usage examples.

## Rebuilding the line list

`linelist_master.csv` is produced by `scripts/build_linelist.py`. Run it from
the repo root whenever the VALD3 extract or element targets change:

```bash
python scripts/build_linelist.py \
    --star 55cnc \
    --vald data/linelists/vald_55cnc_raw.txt \
    --qa
```

Inputs consumed:
- `data/linelists/vald_55cnc_raw.txt` — raw VALD3 long-format extract
- `data/linelists/nist_crosscheck.csv` — NIST grades for Tier 1+2 lines
- `data/config/elements_master.json` — 24 target elements with priorities

The 5 NIST-only science lines (O I, Ni I, C I, P I ×2) are injected automatically
regardless of `--min-depth`. Full process documented in `docs/linelist_pipeline.md`.

## Special handling notes

| Line | Issue | Action |
|------|-------|--------|
| O I 6300.304 | Forbidden line; Ni I 6300.336 blend | Subtract Ni contribution before EW measurement |
| Ni I 6300.336 | O I blend partner | Model jointly with O I |
| Li I 6707.76 | Strong NLTE effects | LTE abundance is lower bound by ~0.1–0.3 dex |
| C I 5380.337 | High excitation (7.685 eV) | Sensitive to Teff errors ±50 K |
| P I 6034, 6043 | Weak lines, depth < 1% | Require S/N > 300; absent from VALD extract |
| N I 7468.31 | Outside HARPS range | Near-IR only — skip for this dataset |

## Grade definitions (NIST ASD)

| Grade | Accuracy | Pipeline use |
|-------|----------|--------------|
| A+ | < 0.3% | Yes |
| A  | < 1%   | Yes |
| B  | < 3%   | Yes (default minimum) |
| C  | < 10%  | Optional (P I lines only) |
| D  | < 25%  | Excluded |

Default minimum grade: `B` (set in `config/constants.py` → `PIPELINE['min_nist_grade']`)
