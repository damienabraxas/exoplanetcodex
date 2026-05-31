# Exoplanet Codex — Data Safety Rules

## THIS DOCUMENT IS INVIOLABLE
These rules apply to ALL code written by Mr. Code (Claude Code) and ALL pipeline runs on any machine (MacBook, Sirius AI, any future node). No exceptions without explicit written instruction from Ryan in the chat.

---

## The Golden Rule

> **Source data is sacred. Pipeline code reads it. Pipeline code never writes to it, moves it, or deletes it.**

---

## Directory Structure & Permissions

```
~/Documents/Exoplanet Codex/
├── data/
│   ├── spectra/          ← READ ONLY. NEVER WRITE HERE. NEVER DELETE HERE.
│   │   ├── Solar Calibration/     (Dumusque HARPS, ESO 1102.D-0954(A))
│   │   ├── 55_Cancri_A/           (88 HARPS S1D files)
│   │   ├── Alpha_Cen_A/           (75 HARPS SNR>1000 + HST downloading)
│   │   └── [future stars]/
│   ├── results/          ← Pipeline outputs go here ONLY
│   │   └── [star]/[YYYY-MM-DD]/
│   ├── output/           ← Plots, diagnostics, intermediate files
│   ├── test/             ← Synthetic/sample data for testing only
│   └── config/
│       ├── elements_master.json
│       └── [other config files]
├── pipeline/
│   ├── spectra_normalize.py
│   ├── lines_fit.py
│   ├── abundances_derive.py
│   ├── uncertainty_calc.py
│   └── build_linelist.py
└── linelist_master.csv
    linelist_solar.csv
```

---

## Hard Rules for Mr. Code

### NEVER do these without explicit Ryan approval:
- `os.remove()`, `shutil.rmtree()`, `rm`, or any deletion command on any file
- Write, overwrite, or move any file inside `data/spectra/`
- Rename or reorganize any FITS file
- Run `git clean`, `git reset --hard`, or any destructive git operation
- Delete or overwrite a `data/results/` folder from a previous run
- Empty, truncate, or overwrite `linelist_master.csv` or `linelist_solar.csv`

### ALWAYS do these:
- Write pipeline outputs to `data/results/[star]/[YYYY-MM-DD]/`
- Write plots and diagnostics to `data/output/`
- Open FITS files read-only: `fits.open(path, mode='readonly')`
- Check that output directories exist before writing; create them if not
- Print a summary of what files were read and what files were written at the end of every run
- Never assume a previous results folder can be overwritten — use a new datestamped subfolder

---

## Backup Architecture

### Current State (MacBook only)
- Time Machine: should be running — verify this is on
- Manual backup: periodically copy `data/spectra/` to an external drive

### Target State (when Sirius AI is online)
```
MacBook (R&D)
    ↕ SSH / rsync
Sirius AI SSD (primary working copy — fast pipeline access)
    ↕ scheduled rsync
External Hard Drive (backup — source data + results)
```

- Source FITS data lives on Sirius SSD for fast pipeline access
- External drive holds a full backup of source data + all results
- MacBook retains a copy of pipeline code and config (via git)
- Sirius AI hardware upgrade (RAM + SSD) required before data migration — tracked in RYA-79/RYA-113

### Backup Checklist (before any major pipeline run)
- [ ] Time Machine last backup: < 24 hours ago
- [ ] `data/spectra/` file count matches expected (Solar: 10, 55 Cnc A: 88, Alpha Cen A: 75+)
- [ ] Previous results folder preserved (not overwritten)

---

## Pipeline Code Safety Patterns

Mr. Code must follow these patterns in all pipeline scripts:

```python
# CORRECT: Open FITS read-only
with fits.open(filepath, mode='readonly') as hdul:
    data = hdul[0].data.copy()  # copy data out before closing

# CORRECT: Write results to datestamped output dir
from datetime import date
output_dir = Path(f"data/results/{star_name}/{date.today().isoformat()}")
output_dir.mkdir(parents=True, exist_ok=True)

# CORRECT: End-of-run summary
print(f"[DONE] Read {n_files} FITS files from data/spectra/{star_name}/")
print(f"[DONE] Results written to {output_dir}")
print(f"[DONE] No source files were modified.")

# NEVER: Write to spectra directory
# with open("data/spectra/some_file.fits", "w") as f:  ← FORBIDDEN

# NEVER: Delete anything
# os.remove(some_file)  ← FORBIDDEN without explicit Ryan approval
# shutil.rmtree(some_dir)  ← FORBIDDEN without explicit Ryan approval
```

---

## If Something Goes Wrong

If Mr. Code accidentally modifies or deletes a file:
1. Stop immediately — do not run anything else
2. Tell Ryan exactly what happened and what files were affected
3. Do not attempt to "fix" it by running more code
4. Check Time Machine for recovery
5. Open a Linear bug ticket documenting what happened

---

## Data Inventory (as of May 2026)

| Dataset | Files | Location | Status |
|---------|-------|----------|--------|
| Solar HARPS (Dumusque) | 10 S1D | data/spectra/Solar Calibration/ | ✅ Complete |
| 55 Cancri A HARPS | 88 S1D | data/spectra/55_Cancri_A/ | ✅ Complete |
| Alpha Cen A HARPS | 75 S1D SNR>1000 | data/spectra/Alpha_Cen_A/ | ✅ Complete |
| Alpha Cen A HST/MAST | ~111 obs | data/spectra/Alpha_Cen_A/HST/ | 🔄 Downloading |
| 55 Cnc A HST | G140M, G750L, G430L, E230M | data/spectra/55_Cancri_A/HST/ | ✅ Complete |
