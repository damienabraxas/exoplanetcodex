# Line List Pipeline: Repeatable Process

This document describes how to build `linelist_master.csv` (and `linelist_solar.csv`) from
scratch for any target star. The process is fully repeatable.

## Overview

```
VALD3 (online)          NIST ASD (online)
      │                       │
      ▼                       ▼
vald_<star>_raw.txt    nist_crosscheck.csv    elements_master.json
      │                       │                       │
      └───────────────────────┴───────────────────────┘
                              │
                    build_linelist.py
                              │
                              ▼
                    linelist_master.csv
                              │
                    data/linelists/loader.py
                              │
                    pipeline/lines_fit.py
```

---

## Step 1: Obtain VALD3 Extract

1. Register (free) at https://vald.astro.uu.se
2. Log in and choose **Extract Stellar** (not Extract All)
3. Enter the target star's atmospheric parameters:

   | Parameter | 55 Cnc A | Solar | Description |
   |-----------|----------|-------|-------------|
   | Teff      | 5196     | 5777  | Effective temperature (K) |
   | log g     | 4.41     | 4.44  | Surface gravity (cgs) |
   | Vmicro    | 0.9      | 1.0   | Microturbulence (km/s) |
   | Wav start | 3780     | 3780  | Start wavelength (Å) |
   | Wav end   | 6910     | 6910  | End wavelength (Å) |
   | Detection | **0.001**| **0.001** | Central depth threshold (0.1%) |
   | Format    | **Long** | **Long** | Must be Long format, not Short |
   | HFS       | **Yes**  | **Yes**  | Enable hyperfine structure |

4. Element list — all 27 targets (use the 26 unique symbols):
   `Fe, C, O, Mg, Si, Ca, Ti, Ni, Na, P, S, N, Co, Cr, Al, K, Ba, Y, V, Cu, Mn, Sc, Li, Eu, Zr, Sr`

5. Submit and wait (typically 5–60 minutes for full optical range)
6. Download via FTP (HTTP truncates large requests) and save to
   `data/linelists/vald_<star>_raw.txt`

> **Format note:** VALD Long format outputs 4 lines per transition. The first line is
> the data line. `build_linelist.py` reads only data lines; reference lines are ignored.

---

## Step 2: Run `build_linelist.py`

From the repo root:

```bash
# For a science target
python scripts/build_linelist.py \
    --star 55cnc \
    --vald data/linelists/vald_55cnc_raw.txt \
    --qa

# For solar calibration
python scripts/build_linelist.py \
    --star solar \
    --vald data/linelists/vald_solar_raw.txt \
    --out data/linelists/linelist_solar.csv \
    --qa
```

The script:
1. Parses all transitions from the VALD3 long-format file
2. Assigns priorities from `data/config/elements_master.json` (27 targets)
3. Flags blend candidates (lines within 0.10 Å of each other)
4. Applies NIST grades from `nist_crosscheck.csv`
5. Injects 2 NIST-only science lines (O I 6300.304, Ni I 6300.336)
6. Writes the sorted master CSV

Expected output for 55 Cnc A (current):
```
[1/6] Parsing VALD3: vald_55cnc_raw.txt
      125,615 transitions parsed
[2/6] Loading elements: elements_master.json
      26 target elements
[3/6] Filtering to target elements
      87,972 target-element lines; 37,643 non-target kept for blend coverage
[4/6] Flagging blends (threshold: 0.10 Å)
      123,090 blend flags set
[5/6] Applying NIST grades: 72 lines graded
[5/6] Injecting NIST-only science lines: 2 lines injected
[6/6] No depth filter applied
Wrote 125,617 lines → data/linelists/linelist_master.csv
```

---

## Step 3: Validate the Output

```python
from data.linelists.loader import load_linelist, summarize_linelist

# All priority-1 lines, no blend filter
df = load_linelist(priority=1, exclude_blends=False, min_nist_grade=None)
summarize_linelist(df)

# Verify O I / Ni I blend pair
o  = df[(df['element']=='O')  & df['wavelength_air_A'].between(6300.3, 6300.31)]
ni = df[(df['element']=='Ni') & df['wavelength_air_A'].between(6300.33, 6300.34)]
assert len(o) >= 1,  "O I 6300.304 missing"
assert len(ni) >= 1, "Ni I 6300.336 missing"
print("O I / Ni I blend pair OK")

# Verify alpha elements at priority 1
for sym in ['Mg','Si','Ca','Ti']:
    lines = df[df['element']==sym]
    assert len(lines) > 0, f"{sym} has no priority-1 lines"
    print(f"{sym}: {len(lines)} priority-1 lines OK")

# Verify new additions
for sym in ['Cu','Sr']:
    n = len(df[df['element']==sym])
    print(f"{sym}: {n} lines (priority may be 2 or 3)")
```

---

## Step 4: Commit

```bash
git add data/linelists/linelist_master.csv
git add data/linelists/vald_<star>_raw.txt  # tracked via Git LFS
git commit -m "RYA-XX: Rebuild linelist for <star>"
git push
```

---

## Adding a New Star

1. Obtain VALD3 extract with new star parameters (Step 1)
2. Save as `data/linelists/vald_<new_star>_raw.txt`
3. Run `build_linelist.py` with `--star <new_star> --out data/linelists/linelist_<new_star>.csv`
4. Same `elements_master.json` and `nist_crosscheck.csv` apply to all solar-type FGK stars

---

## Files Reference

| File | Role |
|------|------|
| `data/linelists/vald_55cnc_raw.txt` | VALD3 extract for 55 Cnc A (125,615 lines, Git LFS) |
| `data/linelists/vald_solar_raw.txt` | VALD3 extract for Sun (108,969 lines, Git LFS) |
| `data/linelists/nist_crosscheck.csv` | Hand-curated NIST grades for Tier 1+2 lines |
| `data/linelists/nist_reference.csv` | NIST A/A+ reference lines for pipeline QA |
| `data/linelists/linelist_master.csv` | Built master list for 55 Cnc A (125,617 rows) |
| `data/linelists/linelist_solar.csv` | Built master list for solar calibration (108,971 rows) |
| `data/config/elements_master.json` | 27 target elements with priorities |
| `scripts/build_linelist.py` | Pipeline builder CLI |
| `data/linelists/loader.py` | Pipeline loader (used by all analysis scripts) |
| `config/constants.py` | Stellar parameters and pipeline settings |

---

## Known Issues and Special Cases

### O I 6300.304 (forbidden line + Ni blend)
Primary oxygen indicator. Two complications:
1. **Forbidden transition** — log gf = −9.717; VALD depth < 0.1% threshold; NIST-injected.
2. **Ni I 6300.336 blend** (Δλ = 0.032 Å) — Ni contribution predicted via linear COG from
   clean Ni I lines and subtracted (Allende Prieto et al. 2001, ApJ 556, L63).

### Ba II 5853 / Eu II 6645 / Li I 6707 (HFS)
All have unresolved hyperfine structure. Treated as single absorption features; the EW
from a single Voigt fit = total HFS EW. See `pipeline/lines_fit.py` → `LINE_WINDOWS`.

### P I 6034/6043 (weak)
Depth < 0.1% at 55 Cnc parameters — NIST-injected into all linelists. Requires S/N > 300.

### C I 5380 (high excitation)
χ = 7.685 eV; 50 K error in Teff → ~0.1 dex uncertainty. Use only after Teff anchored.
