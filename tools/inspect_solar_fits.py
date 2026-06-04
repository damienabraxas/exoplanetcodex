#!/usr/bin/env python3
"""
inspect_solar_fits.py — Survey HARPS FITS headers to identify Ceres vs Vesta
Usage: python inspect_solar_fits.py
Run from the folder containing your FITS files.
"""

import os
import sys
import csv
import glob
import argparse
from astropy.io import fits

parser = argparse.ArgumentParser(description='Survey HARPS FITS headers')
parser.add_argument('--path', default='.', help='Directory to search for FITS files')
args = parser.parse_args()

def get_snr(header):
    """Try several common SNR keyword variants."""
    for key in ['SNR', 'ESO DRS SPE EXT SN', 'ESO QC SNR', 'SN', 'BSS SNR']:
        try:
            val = header.get(key)
            if val is not None:
                return str(round(float(val), 1))
        except:
            pass
    # Try order 50 SNR (mid-range, common in HARPS pipelines)
    for key in [f'HIERARCH ESO DRS SPE EXT SN{i}' for i in [50, 40, 30]]:
        try:
            val = header.get(key)
            if val is not None:
                return str(round(float(val), 1))
        except:
            pass
    return 'N/A'

def classify(object_name, prog_id=''):
    if not object_name or object_name == 'N/A':
        obj = ''
    else:
        obj = object_name.lower()
    if 'ceres' in obj: return 'Ceres'
    if 'vesta' in obj: return 'Vesta'
    if 'sol' in obj or 'sun' in obj: return 'Sol/Sun'
    if '1102.D-0954' in prog_id: return 'Direct Solar (Dumusque)'
    return 'Unknown'

fits_files = sorted(
    glob.glob(os.path.join(args.path, '**/*.fits'), recursive=True) +
    glob.glob(os.path.join(args.path, '**/*.fits.gz'), recursive=True)
)

if not fits_files:
    print("No FITS files found in current directory.")
    sys.exit(1)

print(f"Found {len(fits_files)} FITS files. Inspecting headers...\n")

results = []
counts = {'Ceres': 0, 'Vesta': 0, 'Sol/Sun': 0, 'Unknown': 0}

for fname in fits_files:
    try:
        with fits.open(fname, memmap=True) as hdul:
            h = hdul[0].header
            obj      = str(h.get('OBJECT', 'N/A')).strip()
            date     = str(h.get('DATE-OBS', 'N/A'))
            prog_id  = str(h.get('HIERARCH ESO OBS PROG ID', h.get('PROG-ID', 'N/A')))
            pi       = str(h.get('HIERARCH ESO OBS PI-COI NAME', 'N/A'))
            targ     = str(h.get('HIERARCH ESO OBS TARG NAME', 'N/A'))
            exptime  = str(h.get('EXPTIME', 'N/A'))
            snr      = get_snr(h)
            cls      = classify(obj, prog_id=prog_id)
            counts[cls] += 1
            results.append({
                'filename': fname,
                'OBJECT': obj,
                'DATE-OBS': date,
                'PROG-ID': prog_id,
                'PI': pi,
                'TARGET': targ,
                'EXPTIME': exptime,
                'SNR': snr,
                'classification': cls
            })
    except Exception as e:
        results.append({
            'filename': fname, 'OBJECT': 'ERROR', 'DATE-OBS': '',
            'PROG-ID': '', 'PI': '', 'TARGET': '',
            'EXPTIME': '', 'SNR': '', 'classification': f'ERROR: {e}'
        })

# Print table
hdr = f"{'File':<45} {'OBJECT':<12} {'DATE-OBS':<22} {'PROG-ID':<20} {'SNR':<8} {'Class'}"
print(hdr)
print('-' * len(hdr))
for r in results:
    print(f"{r['filename']:<45} {r['OBJECT']:<12} {r['DATE-OBS']:<22} {r['PROG-ID']:<20} {r['SNR']:<8} {r['classification']}")

# Summary
print(f"\n--- Summary ({len(results)} files) ---")
for k, v in counts.items():
    if v > 0:
        print(f"  {k}: {v}")

# Unique program IDs
prog_ids = set(r['PROG-ID'] for r in results if r['PROG-ID'] not in ('N/A', ''))
print(f"\nUnique Program IDs: {', '.join(sorted(prog_ids))}")

# Save CSV
out = 'inspect_solar_fits_results.csv'
with open(out, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)
print(f"\nResults saved to {out}")
