# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         audit_procyon_hst.py
# Module:       scripts (one-off audit tools)
# Description:  Audits the ~400 HST MAST files downloaded for Procyon
#               (HD 61421). Filters to science-ready STIS 1D extracted
#               spectra, inventories by grating and programme, checks FITS
#               headers, and outputs a filtered file list ready for pipeline
#               ingestion. Non-destructive — reads only, no files moved.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-11
# Last modified: 2026-06-11
# Linear issue: RYA-222 — AUDIT: Procyon HST data filter
#
# -----------------------------------------------------------------------------
# SCIENTIFIC CONTEXT
# -----------------------------------------------------------------------------
# Pipeline stage: Pre-pipeline data audit (before normalization)
# Target star:    Procyon (HD 61421), F5 IV-V, Teff~6530K, [Fe/H]=-0.04
# Instrument:     HST/STIS (multiple gratings)
# Science goal:   UV abundance calibration — C I, O I, Si II, S I, Fe II, Mg II
#                 Cross-instrument validation before 55 Cnc A pipeline run
#
# -----------------------------------------------------------------------------
# KEY REFERENCES
# -----------------------------------------------------------------------------
# Woodgate et al. 1998 — PASP 110, 1183 — STIS instrument description
# Allende Prieto et al. 2002 — ApJ 567, 544 — Procyon reference parameters
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
# External: astropy, numpy, pandas
# Data:     MAST download directory (path set by DATA_DIR below)
# =============================================================================

import os
import sys
from pathlib import Path
from astropy.io import fits
import pandas as pd
import numpy as np
# Standalone-script bootstrap (RYA-313): put the REPO ROOT on sys.path BEFORE
# importing config/pipeline, so this runs from any cwd. Derived from __file__.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

# =============================================================================
# CONFIGURATION — update DATA_DIR to match your local download path
# =============================================================================
DATA_DIR = codex_path('data.spectra_local') / 'Procyon' / 'Procyon HST'
OUTPUT_DIR = Path("data/audit/procyon_hst")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Science-ready gratings for abundance analysis (in priority order)
SCIENCE_GRATINGS = [
    'E140H', 'E140M',   # UV echelle — primary for C, O, Si, S
    'E230H', 'E230M',   # Near-UV echelle — Fe II, Mg II
    'G140M',            # UV first-order — lower R, wider context
    'G230MB',           # Near-UV first-order
]

# Products to keep (1D extracted spectra only)
# _cspec.fits = HASP (HST Advanced Spectral Products) co-added spectra —
# highest S/N combined product; primary science target for pipeline ingestion.
KEEP_SUFFIXES = ['_x1d.fits', '_sx1.fits', '_cspec.fits']

# Products to skip (2D, calibration, imaging)
SKIP_SUFFIXES = [
    '_raw.fits', '_flt.fits', '_crj.fits', '_wav.fits',
    '_tag.fits', '_corrtag.fits', '_acq.fits', '_drz.fits',
    '_asn.fits', '_spt.fits', '_jit.fits', '_jif.fits',
]

# Expected target name variants in FITS headers
# HRS (Goddard High Resolution Spectrograph, pre-1997) files present in archive.
# Not included in KEEP_SUFFIXES — different data format, requires separate loader.
# 24 files detected. Future scope: HRS_AUDIT issue if UV coverage needed pre-1997.

TARGET_NAMES = ['PROCYON', 'HD61421', 'HD 61421', 'HD061421', 'ALPHA CMI', 'ALPHACMI',
                'ALFCMI', 'ALF CMI']  # * alf CMi = HASP Bayer designation for Procyon

# =============================================================================
# AUDIT FUNCTIONS
# =============================================================================

def is_science_product(filepath):
    """
    Check if a file is a 1D extracted spectrum (not calibration/imaging).

    Returns True only for _x1d.fits and _sx1.fits products — the pipeline-
    reduced 1D spectra ready for science analysis.
    """
    name = filepath.name.lower()
    return any(name.endswith(s.lower()) for s in KEEP_SUFFIXES)


def read_header_safe(filepath):
    """
    Safely read primary FITS header, returning None on any error.

    STIS files can have multiple extensions; we check [0] for instrument
    metadata and [1] for science data keywords.
    """
    try:
        with fits.open(filepath, memmap=False) as hdul:
            h0 = hdul[0].header
            # Some STIS products store keywords in extension 1
            h1 = hdul[1].header if len(hdul) > 1 else {}
            return h0, h1
    except Exception as e:
        return None, None


def get_header_value(h0, h1, key, default='UNKNOWN'):
    """Check header h0 first, then h1, then return default."""
    val = h0.get(key, h1.get(key, default) if h1 else default)
    return str(val).strip()


def opt_elem_from_filename(filename):
    """
    Fallback grating extraction for HASP _cspec.fits products.

    HASP single-grating filenames embed the grating in lowercase:
      hst_12278_stis_hd61421_e140h_obkk_cspec.fits → E140H
    Multi-grating combined products (e140h-e230h) return 'UNKNOWN'
    to avoid miscategorisation — they span multiple settings.
    Returns 'UNKNOWN' if zero or multiple gratings found.
    """
    name_lower = filename.lower()
    matches = [g for g in SCIENCE_GRATINGS if g.lower() in name_lower]
    return matches[0] if len(matches) == 1 else 'UNKNOWN'


def audit_file(filepath):
    """
    Audit a single FITS file. Returns a dict with all relevant metadata.

    Checks instrument, grating, target, date, exptime, wavelength range,
    and flags any issues for review.
    """
    record = {
        'filepath': str(filepath),
        'filename': filepath.name,
        'is_science_product': is_science_product(filepath),
        'instrume': 'UNKNOWN',
        'opt_elem': 'UNKNOWN',   # grating/prism
        'cenwave': 'UNKNOWN',    # central wavelength
        'object': 'UNKNOWN',
        'date_obs': 'UNKNOWN',
        'exptime': 0.0,
        'proposid': 'UNKNOWN',   # HST programme ID
        'pr_inv_l': 'UNKNOWN',   # PI last name
        'detector': 'UNKNOWN',
        'wl_min_A': None,
        'wl_max_A': None,
        'target_confirmed': False,
        'is_science_grating': False,
        'flags': [],
    }

    h0, h1 = read_header_safe(filepath)
    if h0 is None:
        record['flags'].append('FITS_READ_ERROR')
        return record

    record['instrume']  = get_header_value(h0, h1, 'INSTRUME')
    record['opt_elem']  = get_header_value(h0, h1, 'OPT_ELEM')
    # HASP _cspec.fits products store grating in filename, not OPT_ELEM header
    if record['opt_elem'] == 'UNKNOWN' and filepath.name.endswith('_cspec.fits'):
        record['opt_elem'] = opt_elem_from_filename(filepath.name)
    record['cenwave']   = get_header_value(h0, h1, 'CENWAVE')
    record['object']    = get_header_value(h0, h1, 'TARGNAME',
                          h0.get('OBJECT', 'UNKNOWN'))
    record['date_obs']  = get_header_value(h0, h1, 'DATE-OBS')
    record['proposid']  = get_header_value(h0, h1, 'PROPOSID')
    record['pr_inv_l']  = get_header_value(h0, h1, 'PR_INV_L')
    record['detector']  = get_header_value(h0, h1, 'DETECTOR')

    try:
        record['exptime'] = float(h0.get('EXPTIME', h1.get('EXPTIME', 0)))
    except (ValueError, TypeError):
        record['exptime'] = 0.0

    # Target confirmation — Procyon has many name variants in HST headers
    obj_upper = record['object'].upper().replace(' ', '').replace('-', '')
    record['target_confirmed'] = any(
        t.replace(' ', '').replace('-', '') in obj_upper
        for t in TARGET_NAMES
    )

    # Science grating check
    record['is_science_grating'] = record['opt_elem'] in SCIENCE_GRATINGS

    # Wavelength range from data extension if available
    try:
        with fits.open(filepath, memmap=False) as hdul:
            if len(hdul) > 1 and 'WAVELENGTH' in hdul[1].columns.names:
                wl = hdul[1].data['WAVELENGTH'].flatten()
                wl = wl[wl > 0]
                if len(wl) > 0:
                    record['wl_min_A'] = float(np.nanmin(wl))
                    record['wl_max_A'] = float(np.nanmax(wl))
                    # STIS wavelengths are in Å — sanity check
                    if record['wl_max_A'] < 100:
                        # Might be in nm — flag it
                        record['flags'].append('WAVELENGTH_UNITS_CHECK_nm?')
    except Exception:
        pass  # Wavelength read is best-effort; not a failure condition

    # Flag issues
    if not record['target_confirmed']:
        record['flags'].append(f"TARGET_MISMATCH: {record['object']}")
    # HASP _cspec.fits co-added products don't populate EXPTIME — suppress flag
    if record['exptime'] <= 0 and not filepath.name.endswith('_cspec.fits'):
        record['flags'].append('ZERO_EXPTIME')
    if record['instrume'] not in ('STIS', 'COS', 'FOS', 'GHRS'):
        record['flags'].append(f"NON_SPECTROSCOPIC: {record['instrume']}")

    return record


# =============================================================================
# MAIN AUDIT
# =============================================================================

def run_audit():
    """
    Run the full HST data audit for Procyon.

    Inventories all FITS files, filters to science products, groups by
    grating, checks headers, and writes summary reports.
    """
    print(f"\n{'='*70}")
    print(f"PROCYON HST DATA AUDIT")
    print(f"{'='*70}")
    print(f"Data directory: {DATA_DIR}")

    # Collect all FITS files
    all_fits = sorted(DATA_DIR.rglob("*.fits"))
    print(f"\nTotal FITS files found: {len(all_fits)}")

    if len(all_fits) == 0:
        print(f"ERROR: No FITS files found in {DATA_DIR}")
        print("Check DATA_DIR path at top of script.")
        sys.exit(1)

    # Audit all files
    print(f"Auditing headers...")
    records = [audit_file(f) for f in all_fits]
    df = pd.DataFrame(records)

    # ==========================================================================
    # SUMMARY REPORT
    # ==========================================================================
    print(f"\n{'─'*70}")
    print("PRODUCT TYPE BREAKDOWN")
    print(f"{'─'*70}")

    # Count by suffix
    for suffix in KEEP_SUFFIXES + SKIP_SUFFIXES:
        count = df['filename'].str.lower().str.endswith(suffix.lower()).sum()
        tag = '✅ KEEP' if suffix in KEEP_SUFFIXES else '   skip'
        if count > 0:
            print(f"  {tag}  {suffix:<20} {count:>4} files")

    # Science products only
    science_df = df[df['is_science_product']].copy()
    print(f"\nScience products (_x1d + _sx1): {len(science_df)} files")

    print(f"\n{'─'*70}")
    print("INSTRUMENT BREAKDOWN (science products only)")
    print(f"{'─'*70}")
    for instrume, grp in science_df.groupby('instrume'):
        print(f"  {instrume}: {len(grp)} files")
        for grating, g2 in grp.groupby('opt_elem'):
            flag = '✅' if grating in SCIENCE_GRATINGS else '  '
            print(f"    {flag}  {grating:<12} {len(g2):>3} files")

    print(f"\n{'─'*70}")
    print("STIS SCIENCE GRATINGS — usable for abundance analysis")
    print(f"{'─'*70}")
    stis_science = science_df[
        (science_df['instrume'] == 'STIS') &
        (science_df['is_science_grating'])
    ].copy()

    if len(stis_science) == 0:
        print("  ⚠️  NO STIS science grating files found!")
        print("  Check that _x1d.fits files are present in the download.")
    else:
        for grating, grp in stis_science.groupby('opt_elem'):
            wl_str = ""
            if grp['wl_min_A'].notna().any():
                wmin = grp['wl_min_A'].min()
                wmax = grp['wl_max_A'].max()
                wl_str = f"  λ {wmin:.0f}–{wmax:.0f} Å"
            print(f"\n  {grating} — {len(grp)} files{wl_str}")
            # Show date range
            dates = sorted(grp['date_obs'].dropna().unique())
            if dates:
                print(f"    Date range: {dates[0]} → {dates[-1]}")
            # Show programmes
            props = grp[['proposid', 'pr_inv_l']].drop_duplicates()
            for _, row in props.iterrows():
                print(f"    Programme {row['proposid']} (PI: {row['pr_inv_l']})")
            # Exptime stats
            exptimes = grp['exptime']
            print(f"    Exptime: min={exptimes.min():.0f}s  "
                  f"max={exptimes.max():.0f}s  "
                  f"total={exptimes.sum():.0f}s")

    print(f"\n{'─'*70}")
    print("TARGET CONFIRMATION")
    print(f"{'─'*70}")
    n_confirmed = science_df['target_confirmed'].sum()
    n_mismatch = (~science_df['target_confirmed']).sum()
    print(f"  Target confirmed: {n_confirmed} files")
    if n_mismatch > 0:
        print(f"  ⚠️  Target mismatch: {n_mismatch} files — check manually")
        mismatch = science_df[~science_df['target_confirmed']][['filename','object']]
        print(mismatch.to_string(index=False))

    print(f"\n{'─'*70}")
    print("FLAGS / ISSUES")
    print(f"{'─'*70}")
    flagged = df[df['flags'].apply(len) > 0]
    if len(flagged) == 0:
        print("  No issues found ✅")
    else:
        for _, row in flagged.iterrows():
            print(f"  {row['filename']}: {', '.join(row['flags'])}")

    # ==========================================================================
    # OUTPUT FILES
    # ==========================================================================

    # Full inventory
    inventory_path = OUTPUT_DIR / "procyon_hst_full_inventory.csv"
    df.to_csv(inventory_path, index=False)
    print(f"\n{'─'*70}")
    print(f"Full inventory saved: {inventory_path}")

    # Science-ready filtered list (STIS science gratings, target confirmed)
    ready_df = stis_science[stis_science['target_confirmed']].copy()
    ready_path = OUTPUT_DIR / "procyon_hst_science_ready.csv"
    ready_df.to_csv(ready_path, index=False)
    print(f"Science-ready file list saved: {ready_path}  ({len(ready_df)} files)")

    # Per-grating file lists
    for grating, grp in ready_df.groupby('opt_elem'):
        grating_path = OUTPUT_DIR / f"procyon_hst_{grating.lower()}_files.txt"
        grp['filepath'].to_csv(grating_path, index=False, header=False)
        print(f"  {grating} file list: {grating_path}  ({len(grp)} files)")

    print(f"\n{'='*70}")
    print("AUDIT COMPLETE")
    print(f"Science-ready STIS files: {len(ready_df)} / {len(all_fits)} total")
    print(f"{'='*70}\n")

    return df, ready_df


if __name__ == '__main__':
    df, ready_df = run_audit()
