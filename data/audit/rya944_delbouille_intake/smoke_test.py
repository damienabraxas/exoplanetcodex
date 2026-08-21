"""RYA-944 smoke test, run against the atlas AS DELIVERED (no vac->air conversion).

The spec's medium assertion is reproduced verbatim AND replaced with one that can
actually fail -- see findings.md section 2 for why the original cannot.
"""
import csv, gzip, hashlib, sys
from pathlib import Path
import numpy as np

src = Path(sys.argv[1] if len(sys.argv) > 1 else
           "/mnt/codex-data/solar_reference/delbouille_bass2000/delbouille_bass2000_visible.csv.gz")
sha = hashlib.sha256(src.read_bytes()).hexdigest()
w, f = [], []
with gzip.open(src, "rt") as fh:
    for i, r in enumerate(csv.reader(fh)):
        if i == 0: continue
        w.append(float(r[0])); f.append(float(r[1]))
wav, flux = np.array(w), np.array(f)

assert 3000 <= wav.min() and wav.max() >= 9999, "coverage gap"

# --- the spec's assertion, verbatim -----------------------------------------
idx = np.argmin(np.abs(wav - 6430.846))
assert abs(wav[idx] - 6430.846) < 0.002, "vac->air FAILED"
print(f"spec assertion      : PASS (grid point {wav[idx]:.4f})")
print("                      ^ passes in ANY medium; it tests grid spacing, not medium")

# --- the assertion that can actually fail -----------------------------------
def centroid(pos, half=0.09):
    m = (wav > pos - half) & (wav < pos + half)
    ww, ff = wav[m], flux[m]
    d = np.clip(ff.max() - ff, 0, None)
    return float((ww * d).sum() / d.sum())

def n_air(l):
    s2 = (1e4 / l) ** 2
    return 1 + (0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))

air = 6430.846
c = centroid(air)
vac = air * n_air(air)
assert abs(c - air) < abs(c - vac), f"medium is NOT air: centroid {c:.4f}"
assert abs(c - air) < 0.020, f"air centroid off by {(c-air)*1e3:.1f} mA"
print(f"medium (centroid)   : AIR  ({c:.4f}, {(c-air)*1e3:+.1f} mA vs air, "
      f"{(c-vac)*1e3:+.1f} mA vs vacuum)")

med = np.median(flux[(wav > 6425) & (wav < 6428)])
assert 0.97 < med < 1.03, f"normalization off: {med}"
print(f"normalization       : PASS (median {med:.4f} in 6425-6428 A)")
print(f"continuum ceiling   : max flux over {wav.size} points = {flux.max():.4f} (never 1.0)")
print(f"\nSHA-256: {sha} | coverage: {wav.min()} - {wav.max()} | geometry: disk-center intensity mu=1.0")
