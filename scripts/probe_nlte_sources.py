"""
scripts/probe_nlte_sources.py
==============================
Reconnaissance script — probe available NLTE correction sources for
Ca, Ti, and Cr. Reports what grids are available and their parameter
coverage. Does not modify any pipeline files.

Run this first. Post full output to Linear before proceeding to Step 2.
"""

import requests

MPIA_BASE = "https://nlte.mpia.de"

# 1. Check MPIA main page for available elements
print("=== MPIA NLTE Grid Server ===")
r = requests.get(MPIA_BASE, timeout=30)
print(f"Status: {r.status_code}")
print(r.text[:3000])  # Print first 3000 chars — look for element list

# 2. Try common MPIA API patterns for Ca, Ti, Cr
# The MPIA server uses element-specific endpoints; try plausible patterns:
for element in ['Ca', 'Ti', 'Cr']:
    for path in [f'/api/{element}', f'/grids/{element}', f'/{element}',
                 f'/nlte/{element}', f'/corrections/{element}']:
        try:
            resp = requests.get(f"{MPIA_BASE}{path}", timeout=10)
            print(f"{element} {path}: HTTP {resp.status_code}  len={len(resp.text)}")
            if resp.status_code == 200:
                print(resp.text[:500])
                break
        except Exception as e:
            print(f"{element} {path}: ERROR {e}")

# 3. Check correct Mashonkina Ca VizieR ID
print("\n=== VizieR check: J/A+A/601/A96 ===")
vizier_url = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/601/A96/"
r2 = requests.get(vizier_url, timeout=30)
print(f"Status: {r2.status_code}")
if r2.status_code == 200:
    print(r2.text[:1000])
else:
    print("Not found at this path — try alternate:")
    # Try alternate CDS FTP path
    alt = "https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/601/A96"
    r3 = requests.get(alt, timeout=30)
    print(f"Alternate status: {r3.status_code}")
    print(r3.text[:500])
