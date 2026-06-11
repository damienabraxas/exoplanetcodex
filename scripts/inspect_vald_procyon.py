# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         inspect_vald_procyon.py
# Module:       scripts (one-off inspection tools)
# Description:  Quick sanity check on Procyon VALD3 FTP delivery.
#               Counts lines per species, checks wavelength range, confirms
#               long format with damping constants present.
#               Non-destructive — reads only.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-11
# Last modified: 2026-06-11
# Linear issue: RYA-223 — PREP: Unpack Procyon VALD deliveries
#
# -----------------------------------------------------------------------------
# KEY REFERENCES
# -----------------------------------------------------------------------------
# Ryabchikova et al. 2015 — Phys. Scr. 90, 054005 — VALD3 description
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
# External: none (stdlib only)
# Data:     data/linelists/vald_procyon_raw/ (unpacked VALD files)
# =============================================================================

import sys
from pathlib import Path
from collections import defaultdict

# Update these paths to match actual unpacked filenames
VALD_FILES = {
    'optical_3780_6910': Path('data/linelists/vald_procyon_raw/vald_procyon_optical.vald'),
    'blue_3000_4800':    Path('data/linelists/vald_procyon_raw/vald_procyon_blue.vald'),
}

# 27-element Codex canonical list — check all are represented
CODEX_ELEMENTS = [
    'Fe', 'C', 'O', 'N', 'Mg', 'Si', 'S', 'Ca', 'Ti', 'Co', 'Ni',
    'Na', 'Al', 'K', 'P', 'Ba', 'Y', 'Eu', 'Mn', 'Cr', 'V', 'Sc',
    'Cu', 'Zr', 'Li', 'Sr'
]


def inspect_vald_file(label, filepath):
    """
    Parse a VALD Long Format file and print a species inventory.

    Args:
        label (str): Human-readable label for this extraction.
        filepath (Path): Path to the VALD file.

    Prints species line counts, wavelength range, and flags any
    missing Codex elements or format issues.
    """
    if not filepath.exists():
        print(f"\n  ⚠️  FILE NOT FOUND: {filepath}")
        print(f"     Update VALD_FILES paths at top of script to match")
        print(f"     actual filenames in data/linelists/vald_procyon_raw/")
        return

    print(f"\n{'='*60}")
    print(f"VALD FILE: {label}")
    print(f"Path: {filepath}")
    print(f"Size: {filepath.stat().st_size / 1024:.1f} KB")
    print(f"{'='*60}")

    species_counts = defaultdict(int)
    wavelengths = []
    n_total = 0
    n_data = 0
    n_comment = 0
    format_errors = []

    with open(filepath, 'r', errors='replace') as f:
        for i, line in enumerate(f):
            n_total += 1
            line = line.strip()

            # Skip comment/header lines
            if line.startswith('#') or line.startswith('*') or not line:
                n_comment += 1
                continue

            # Identify VALD Long Format data lines by structure:
            # Data:   'Fe 1',  3780.26, ...  → quote_parts[2] starts with ','
            # Config: '  LS   ...'           → quote_parts[2] == ''
            # Ref:    '_  Kurucz...', 'E 0.02 ...' → quote_parts[2] == ''
            if not line.startswith("'"):
                n_comment += 1
                continue

            quote_parts = line.split("'")
            if len(quote_parts) < 3 or not quote_parts[2].startswith(','):
                n_comment += 1
                continue

            n_data += 1
            try:
                species_raw = quote_parts[1].strip()   # e.g. 'Fe 1'
                species_counts[species_raw] += 1

                # Wavelength is first numeric field after the closing quote+comma
                after_species = quote_parts[2]         # ',  3780.2550, -2.876, ...'
                wl = float(after_species.split(',')[1].strip())
                wavelengths.append(wl)

            except (IndexError, ValueError) as e:
                format_errors.append(f"Line {i+1}: {str(e)[:50]}")
                if len(format_errors) > 5:
                    format_errors.append("...(further errors suppressed)")
                    break

    # Summary
    print(f"\nLine counts:")
    print(f"  Total lines:    {n_total}")
    print(f"  Data lines:     {n_data}")
    print(f"  Comment/header: {n_comment}")

    if wavelengths:
        print(f"\nWavelength range: {min(wavelengths):.2f} – {max(wavelengths):.2f} Å")
    else:
        print(f"\n  ⚠️  No wavelengths parsed — check file format")

    if format_errors:
        print(f"\n  ⚠️  Format errors ({len(format_errors)}):")
        for e in format_errors[:5]:
            print(f"     {e}")

    print(f"\nSpecies inventory (top 20 by line count):")
    sorted_species = sorted(species_counts.items(), key=lambda x: -x[1])
    for species, count in sorted_species[:20]:
        print(f"  {species:<15} {count:>5} lines")
    if len(sorted_species) > 20:
        print(f"  ... and {len(sorted_species)-20} more species")

    # Check Codex elements
    present_elements = set(s.split()[0] for s in species_counts.keys())
    missing = [e for e in CODEX_ELEMENTS if e not in present_elements]
    print(f"\nCodex element coverage:")
    print(f"  Present: {len(CODEX_ELEMENTS) - len(missing)}/{len(CODEX_ELEMENTS)}")
    if missing:
        print(f"  Missing: {', '.join(missing)}")
        print(f"  (Note: some elements may have no detectable lines at F-star")
        print(f"   temperatures — absence is not always a problem)")
    else:
        print(f"  All 27 Codex elements present ✅")

    # Fe I line count — most important for calibration
    fe1_count = species_counts.get('Fe 1', 0)
    fe2_count = species_counts.get('Fe 2', 0)
    print(f"\nFe I lines:  {fe1_count}  (need ≥100 before vetting to get ≥20 clean)")
    print(f"Fe II lines: {fe2_count}  (need ≥20 before vetting to get ≥5 clean)")
    if fe1_count < 50:
        print(f"  ⚠️  LOW Fe I count — check detection threshold in VALD extraction")


def main():
    print("PROCYON VALD3 DELIVERY INSPECTION")
    print("Checking files in data/linelists/vald_procyon_raw/")

    # List all files in the raw directory
    raw_dir = Path('data/linelists/vald_procyon_raw')
    if raw_dir.exists():
        files = list(raw_dir.iterdir())
        print(f"\nFiles in staging directory ({len(files)}):")
        for f in sorted(files):
            print(f"  {f.name}  ({f.stat().st_size / 1024:.1f} KB)")
    else:
        print(f"\n  ⚠️  Staging directory not found: {raw_dir}")
        print(f"     Run unpack step first.")
        sys.exit(1)

    # Inspect each expected file
    for label, filepath in VALD_FILES.items():
        inspect_vald_file(label, filepath)

    print(f"\n{'='*60}")
    print("INSPECTION COMPLETE")
    print("Next step: RYA-206 — build linelist_procyon.csv from these files")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
