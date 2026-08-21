"""RYA-944: assemble the BASS2000 visible arm into one atlas and audit it.

Every assertion here is a RELATIONSHIP, not a memorised constant, so the audit
travels if BASS2000 re-reduces the atlas.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, sys
from pathlib import Path
import numpy as np

SCRATCH = Path("/private/tmp/claude-501/-Users-ryanschmitt/07bd9c77-6066-4985-9818-751ff0c31e9c/scratchpad")
CHUNKS = SCRATCH / "delb"
OUT = SCRATCH / "delbouille_bass2000_visible.csv.gz"

# ---------------------------------------------------------------- assemble
waves, fluxes, ident = [], [], 0
files = sorted(CHUNKS.glob("chunk_*.csv"))
if not files:
    sys.exit("no chunks")
for f in files:
    for i, r in enumerate(csv.reader(f.read_text(encoding="utf-8", errors="replace").splitlines())):
        if i == 0 or len(r) < 2:
            continue
        try:
            w = float(r[0]); v = float(r[1])
        except ValueError:
            continue
        waves.append(w); fluxes.append(v)
        if len(r) > 2 and r[2].strip():
            ident += 1
w = np.asarray(waves, float); f = np.asarray(fluxes, float) / 1e4

order = np.argsort(w, kind="stable")
w, f = w[order], f[order]
keep = np.ones(w.size, bool); keep[1:] = np.diff(w) > 0     # drop chunk-seam duplicates
dupes = int((~keep).sum())
w, f = w[keep], f[keep]
assert np.all(np.diff(w) > 0), "wavelengths not strictly increasing after dedupe"

# ---------------------------------------------------------------- write
with gzip.open(OUT, "wt", newline="") as fh:
    out = csv.writer(fh)
    out.writerow(["wavelength_A_air", "intensity_norm"])
    for a, b in zip(w, f):
        out.writerow([f"{a:.4f}", f"{b:.5f}"])
sha = hashlib.sha256(OUT.read_bytes()).hexdigest()

# ---------------------------------------------------------------- audit
def n_air(l):
    s2 = (1e4 / l) ** 2
    return 1 + (0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))

def centroid(pos, halfwin=0.09):
    """Flux-weighted centroid of the absorption core -- sharper than argmin."""
    m = (w > pos - halfwin) & (w < pos + halfwin)
    if m.sum() < 5:
        return None
    ww, ff = w[m], f[m]
    d = np.clip(ff.max() - ff, 0, None)
    return float((ww * d).sum() / d.sum()) if d.sum() > 0 else None

# air rest wavelengths of clean, strong, unblended solar lines
LINES = [("Fe I 4045", 4045.812), ("Fe I 4383", 4383.545), ("Mg I b2 5172", 5172.684),
         ("Mg I b1 5183", 5183.604), ("Na D2 5889", 5889.951), ("Na D1 5895", 5895.924),
         ("Fe I 6430", 6430.846), ("H-alpha 6562", 6562.797), ("Fe I 8688", 8688.624)]
medium = []
for name, air in LINES:
    c = centroid(air)
    if c is None:
        continue
    vac = air * n_air(air)
    medium.append({"line": name, "air_rest_A": air, "vac_rest_A": round(vac, 3),
                   "measured_A": round(c, 4), "resid_vs_air_mA": round((c - air) * 1e3, 1),
                   "resid_vs_vac_mA": round((c - vac) * 1e3, 1),
                   "verdict": "AIR" if abs(c - air) < abs(c - vac) else "VACUUM"})

# normalization: continuum should sit at 1.0 in clean windows
CLEAN = [(4700, 4703), (5240, 5243), (6425, 6428), (7480, 7483), (8710, 8713)]
norm = []
for a, b in CLEAN:
    m = (w > a) & (w < b)
    if m.any():
        norm.append({"window_A": [a, b], "median": round(float(np.median(f[m])), 4),
                     "p95": round(float(np.percentile(f[m], 95)), 4)})

steps = np.diff(w)
audit = {
    "source_url": "https://bass2000.obspm.fr/php/getSolarSpectrumDB.php?WL=<start>&DW=<width>&resol=0.002&fmt=txt",
    "atlas": "BASS2000 visible arm - Jungfraujoch (Delbouille, Neven & Roland 1972)",
    "geometry": "disk-center intensity (mu=1.0)",
    "n_points": int(w.size), "seam_duplicates_dropped": dupes,
    "coverage_A": [round(float(w.min()), 4), round(float(w.max()), 4)],
    "sampling_A": {"min": round(float(steps.min()), 5), "median": round(float(np.median(steps)), 5),
                   "max": round(float(steps.max()), 5)},
    "largest_gap_A": round(float(steps.max()), 5),
    "intensity_range": [round(float(f.min()), 5), round(float(f.max()), 5)],
    "rows_with_line_identification": ident,
    "sha256": sha, "bytes": OUT.stat().st_size,
    "wavelength_medium_test": medium,
    "normalization_test": norm,
}
(SCRATCH / "rya944_audit.json").write_text(json.dumps(audit, indent=2))
print(json.dumps(audit, indent=2))
