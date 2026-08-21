"""RYA-944: is the Delbouille atlas telluric-contaminated, and how badly?

An absence is a hypothesis, never a conclusion (RYA-833). So Delbouille is measured
BETWEEN two atlases whose telluric state is already established in
pipeline/telluric_policy.py:

    Kitt Peak flux   tellurics RETAINED  -> flux driven to zero in the bands
    IAG FTS          tellurics CORRECTED -> sits at continuum in the same bands

If Delbouille lands with KP it retains tellurics; if it lands with IAG it does not.
The clean-continuum window is the control that proves the metric discriminates at all.
"""
from __future__ import annotations
import csv, glob, gzip, os
import numpy as np

WINDOWS = [("O2 B-band",      6867.0, 6884.0),
           ("O2 A-band",      7594.0, 7685.0),
           ("H2O 8100-8400",  8100.0, 8400.0),
           ("H2O 9280-9600",  9280.0, 9600.0),
           ("clean continuum",6425.0, 6525.0)]

def n_air(l):
    s2 = (1e4/l)**2
    return 1 + (0.0000834254 + 0.02406147/(130-s2) + 0.00015998/(38.9-s2))

def load_delbouille():
    w, f = [], []
    with gzip.open("/mnt/codex-data/solar_reference/delbouille_bass2000/"
                   "delbouille_bass2000_visible.csv.gz", "rt") as fh:
        for i, r in enumerate(csv.reader(fh)):
            if i == 0: continue
            w.append(float(r[0])); f.append(float(r[1]))
    return np.array(w), np.array(f)

def load_kp():
    base = "/mnt/codex-data/spectra/Solar Calibration/Kitt Peak Flux Atlas"
    w, f = [], []
    for p in sorted(glob.glob(os.path.join(base, "lm[0-9]*"))):
        if os.path.basename(p).startswith("._"): continue
        try:
            d = np.loadtxt(p)
        except Exception:
            continue
        if d.ndim != 2 or d.shape[1] < 2: continue
        w.append(d[:, 0] * 10.0); f.append(d[:, 1])       # col0 nm -> A (air)
    w = np.concatenate(w); f = np.concatenate(f)
    o = np.argsort(w); return w[o], f[o]

def load_iag():
    d = np.loadtxt("/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz")
    vac = 1e8 / d[:, 0]                                    # col0 VACUUM wavenumber cm-1
    air = vac / n_air(vac)
    o = np.argsort(air); return air[o], d[o, 1]

def stats(w, f, lo, hi):
    m = (w >= lo) & (w <= hi)
    if m.sum() < 50: return None
    y = f[m]
    return {"n": int(m.sum()), "min": float(np.nanmin(y)), "mean": float(np.nanmean(y)),
            "pct50": float(np.nanmean(y < 0.5) * 100)}

atlases = []
for name, loader in (("Delbouille", load_delbouille), ("KittPeak", load_kp), ("IAG", load_iag)):
    try:
        atlases.append((name, *loader()))
        print(f"loaded {name}: {atlases[-1][1].size} pts "
              f"{atlases[-1][1].min():.1f}-{atlases[-1][1].max():.1f} A", flush=True)
    except Exception as e:
        print(f"FAILED {name}: {e}", flush=True)

print(f'\n{"window":<18}{"atlas":<12}{"n":>9}{"min":>9}{"mean":>9}{"%<0.5":>9}')
for lbl, lo, hi in WINDOWS:
    for name, w, f in atlases:
        s = stats(w, f, lo, hi)
        if s is None:
            print(f'{lbl:<18}{name:<12}{"-- not covered --":>36}')
        else:
            print(f'{lbl:<18}{name:<12}{s["n"]:>9}{s["min"]:>9.3f}{s["mean"]:>9.3f}{s["pct50"]:>9.2f}')
    print()

# ---- appended: the ACTUALLY telluric-corrected IAG (Baker+2020 "telfree") ----
from astropy.io import fits

def load_baker():
    d = fits.open("/mnt/codex-data/solar_reference/iag_baker2020/"
                  "iag_telfree_solaratlas.fits")[1].data
    vac = 1e8 / d["v"].astype(float)
    air = vac / n_air(vac)
    o = np.argsort(air)
    return air[o], d["s"].astype(float)[o]

wb, fb = load_baker()
print("\nloaded IAG-Baker2020(telfree): %d pts %.1f-%.1f A" % (wb.size, wb.min(), wb.max()))
hdr = "%-18s%-24s%9s%9s%9s%9s" % ("window", "atlas", "n", "min", "mean", "%<0.5")
print("\n" + hdr)
rename = {"IAG": "IAG-Reiners2016 RAW", "Delbouille": "Delbouille (NEW)"}
for lbl, lo, hi in WINDOWS:
    for name, w, f in atlases + [("IAG-Baker2020 telfree", wb, fb)]:
        s = stats(w, f, lo, hi)
        tag = rename.get(name, name)
        if s is None:
            print("%-18s%-24s%36s" % (lbl, tag, "-- not covered --"))
        else:
            print("%-18s%-24s%9d%9.3f%9.3f%9.2f" % (lbl, tag, s["n"], s["min"], s["mean"], s["pct50"]))
    print()
