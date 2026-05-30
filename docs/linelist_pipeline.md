# Line List Pipeline: Repeatable Process

This document describes how to build `linelist_master.csv` from scratch for any target star.
The process is fully repeatable: running the same VALD3 input through `build_linelist.py`
always produces the same output.

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
                    pipeline/03_linelist.py
```

---

## Step 1: Obtain VALD3 Extract

1. Register (free) at https://vald.astro.uu.se
2. Log in and choose **Extract Stellar** (not Extract All)
3. Enter the target star's atmospheric parameters:

   | Parameter | 55 Cnc A value | Description |
   |-----------|---------------|-------------|
   | Teff      | 5196          | Effective temperature (K) |
   | log g     | 4.41          | Surface gravity (cgs) |
   | [Fe/H]    | +0.32         | Metallicity |
   | Vmicro    | 0.9           | Microturbulence (km/s) |
   | Wav start | 3780          | Start wavelength (Å) |
   | Wav end   | 6910          | End wavelength (Å) |
   | Detection | 0.01          | Central depth threshold (1%) |
   | Format    | **Long**      | Must be Long format, not Short |

4. Submit and wait (typical: 5–30 minutes for a full optical range extract)
5. Download the result and save to `data/linelists/vald_<star>_raw.txt`

> **Format note:** VALD Long format outputs 4 lines per transition. The first line is
> the data line (`'El Num', wavelength, log_gf, E_low, ...`). The other three lines are
> configuration labels and bibliographic references. `build_linelist.py` reads only the
> data lines; the reference lines are ignored.

---

## Step 2: Update NIST Crosscheck (if adding a new star)

`data/linelists/nist_crosscheck.csv` contains hand-curated NIST ASD values for the
20 Tier 1 + Tier 2 science-critical lines (defined in RYA-64). These grades are applied
by `build_linelist.py` to override whatever VALD3 lists.

For a new star, verify:
- The same 20 lines are present and measurable at the new star's metallicity/Teff
- NIST ASD grades have not been revised since last retrieval (check NIST ASD version)

Source: https://physics.nist.gov/PhysRefData/ASD/lines_form.html

---

## Step 3: Run `build_linelist.py`

From the repo root:

```bash
python scripts/build_linelist.py \
    --star 55cnc \
    --vald data/linelists/vald_55cnc_raw.txt \
    --qa
```

The script will:
1. Parse all transitions from the VALD3 long-format file
2. Filter to the 24 target elements in `data/config/elements_master.json`
3. Flag blend candidates (lines within 0.10 Å of each other)
4. Apply NIST grades from `nist_crosscheck.csv` to matching lines
5. Inject 5 NIST-only science lines (O I, Ni I, C I, P I ×2) that fall below the 1% VALD threshold
6. Write the sorted master CSV to `data/linelists/linelist_master.csv`

Expected output for 55 Cnc A:
```
[1/6] Parsing VALD3: vald_55cnc_raw.txt
      17,926 transitions parsed
[2/6] Loading elements: elements_master.json
      24 target elements
[3/6] Filtering to target elements
      17,926 lines retained
[4/6] Flagging blends (threshold: 0.10 Å)
      N blend flags set
[5/6] Applying NIST grades from: nist_crosscheck.csv
      20 lines graded
[5/6] Injecting NIST-only science lines
      5 lines injected (O I, Ni I, C I, P I x2)
[6/6] No depth filter applied
Wrote 17,931 lines → data/linelists/linelist_master.csv
```

---

## Step 4: Validate the Output

### 4a. Count check

```bash
wc -l data/linelists/linelist_master.csv
# Expected: 17,932 lines (1 header + 17,931 data rows)
```

### 4b. NIST injection check

```bash
grep "^O,I,6300" data/linelists/linelist_master.csv
grep "^C,I,5380" data/linelists/linelist_master.csv
grep "^P,I,603"  data/linelists/linelist_master.csv
```

All five NIST-injected lines must be present with `loggf_source=NIST`.

### 4c. Python QA

```python
from data.linelists.loader import load_linelist, summarize_linelist

# Load all science-critical lines (priority 1, no blend filter)
df = load_linelist(priority=1, exclude_blends=False, min_nist_grade=None)
summarize_linelist(df)

# Verify O I and Ni I blend pair present
o_line = df[(df['element'] == 'O') & (df['wavelength_air_A'].between(6300.3, 6300.31))]
ni_line = df[(df['element'] == 'Ni') & (df['wavelength_air_A'].between(6300.33, 6300.34))]
assert len(o_line) == 1, "O I 6300.304 missing"
assert len(ni_line) == 1, "Ni I 6300.336 missing"
assert o_line.iloc[0]['blend_flag'] == True
assert ni_line.iloc[0]['blend_flag'] == True
print("O I / Ni I blend pair OK")

# Verify Li doublet
li = df[df['element'] == 'Li']
assert len(li) == 2, f"Expected 2 Li lines, got {len(li)}"
print("Li I doublet OK")
```

### 4d. Special lines quick-check

| Line | Check |
|------|-------|
| O I 6300.304 | `blend_flag=True`, `loggf_source=NIST`, `nist_grade=A` |
| Ni I 6300.336 | `blend_flag=True`, `loggf_source=NIST`, `nist_grade=B` |
| C I 5380.337 | `loggf_source=NIST`, `excitation_potential_eV=7.685` |
| Li I 6707.76 | `nist_grade=A+`, NLTE note in `notes` column |
| P I 6034.04 | `nist_grade=C`, `loggf_source=NIST` |
| Fe I 5576 | `nist_grade=A` (from crosscheck) |

---

## Step 5: Commit

```bash
git add data/linelists/linelist_master.csv
git add data/linelists/vald_<star>_raw.txt  # only on first add; file is ~7 MB
git commit -m "RYA-XX: Rebuild linelist_master for <star> from VALD3 <date>"
git push origin main
```

---

## Adding a New Star

1. Obtain VALD3 extract with new star parameters (Step 1)
2. Save as `data/linelists/vald_<new_star>_raw.txt`
3. Run:
   ```bash
   python scripts/build_linelist.py \
       --star <new_star> \
       --vald data/linelists/vald_<new_star>_raw.txt \
       --out data/linelists/linelist_<new_star>.csv \
       --qa
   ```
4. The same `elements_master.json` and `nist_crosscheck.csv` apply to all solar-type stars.
   Update `nist_crosscheck.csv` if the new star's parameters significantly change which
   NIST lines are detectable (e.g., very different Teff → different excitation sensitivities).
5. Update `config/constants.py` to add the new star's parameters to a `STAR_<NAME>` dict.

---

## Files Reference

| File | Role |
|------|------|
| `data/linelists/vald_55cnc_raw.txt` | Raw VALD3 long-format extract (17,926 transitions) |
| `data/linelists/nist_crosscheck.csv` | Hand-curated NIST grades for 20 Tier 1+2 lines |
| `data/linelists/nist_reference.csv` | NIST A/A+ reference lines for pipeline QA |
| `data/linelists/linelist_master.csv` | Built master list (17,931 rows) |
| `data/config/elements_master.json` | 24 target elements with priorities |
| `scripts/build_linelist.py` | Pipeline builder (this step's CLI tool) |
| `data/linelists/loader.py` | Pipeline loader (used by all analysis scripts) |
| `config/constants.py` | Stellar parameters and pipeline settings |

---

## Known Issues and Special Cases

### O I 6300.304 (forbidden line)
This [O I] line is the primary oxygen indicator but has two complications:
1. It is **forbidden** (magnetic dipole transition) — the transition probability is 5.63×10⁻³ s⁻¹,
   three orders of magnitude weaker than typical permitted lines. VALD's 1% depth threshold
   excludes it; it is injected from NIST.
2. It is **blended** with Ni I 6300.336 (Δλ = 0.032 Å). The Ni I contribution must be
   subtracted before measuring the O I equivalent width. See `pipeline/04_ew_measure.py`.

### Li I 6707 doublet (NLTE)
The Li I resonance doublet is the primary age/activity diagnostic but LTE abundances are
systematically too low by ~0.1–0.3 dex (increasing toward lower Teff and higher [Fe/H]).
Apply NLTE corrections from Lind et al. 2009 (A&A 503, 541) post-EW measurement.

### P I 6034/6043 (weak lines)
Both phosphorus lines have predicted central depth < 0.6% at 55 Cnc parameters.
Equivalent widths are measurable only at S/N ≥ 300 per pixel. NIST grades are C (<10%
accuracy on log gf), introducing ~0.04 dex systematic uncertainty in the abundance.

### C I 5380.337 (high excitation)
The lower level at 7.685 eV makes this line extremely sensitive to Teff. A 50 K error
in Teff propagates to ~0.1 dex uncertainty in [C/H]. Use only after Teff is anchored
by Fe I excitation equilibrium.
