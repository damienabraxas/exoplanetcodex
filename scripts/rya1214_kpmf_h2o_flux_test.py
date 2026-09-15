"""RYA-1214 — is the KP molecfit 1984 composite contaminated at 9250-9700 A?

Direct flux test, no synthesis: the same solar lines on three holdings. K05 (Kurucz 2005
residual, telluric-corrected at source) and IAG (FTS, corrected) are the references.
Suspect windows: the O I 926 triplet and C I 9603 / 9658, where KP-mf fits went +1 to +3.8
dex. Controls: windows where KP-mf C I fits were normal (10683, 10691, 11330 A) and a
red-optical window (C I 8335) — a holding-wide offset would show there too.
Per window: median |KP-mf - K05| and |KP-mf - IAG|, the minimum flux on each holding (a
telluric core drives flux down), and the K05-vs-IAG difference as the reference-pair null.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from measure_band_ew import load_window_ex  # noqa: E402

HOLD = {"kpmf": ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
        "k05": ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
        "iag": ("iag_fts_solar_atlas", "solar_iag")}
WINDOWS = [("SUSPECT O I 926.1", 9260.87, 1.0), ("SUSPECT O I 926.3", 9262.68, 1.0),
           ("SUSPECT O I 926.6", 9265.92, 1.0), ("SUSPECT C I 9603", 9603.03, 1.4),
           ("SUSPECT C I 9658", 9658.44, 1.4),
           ("CONTROL C I 8335", 8335.15, 1.1), ("CONTROL C I 10683", 10683.08, 1.4),
           ("CONTROL C I 10691", 10691.25, 1.4), ("CONTROL C I 11330", 11330.28, 1.4)]


def grid(name, c, pad):
    inst, hold = HOLD[name]
    try:
        w = load_window_ex(inst, c, pad, holding=hold)
    except Exception as e:  # noqa: BLE001
        return None, f"{type(e).__name__}: {str(e)[:70]}"
    x, y = np.asarray(w.wave, float), np.asarray(w.flux, float)
    m = np.isfinite(x) & np.isfinite(y)
    return (x[m], y[m]), None


rows = []
print(f"{'window':22s} {'min kpmf':>8s} {'min k05':>8s} {'min iag':>8s} "
      f"{'|mf-k05|':>9s} {'|mf-iag|':>9s} {'|k05-iag|':>9s}")
for label, c, pad in WINDOWS:
    g = {k: grid(k, c, pad) for k in HOLD}
    ok = {k: v[0] for k, v in g.items() if v[0] is not None and len(v[0][0]) > 5}
    common = np.linspace(c - pad, c + pad, 400)
    f = {k: np.interp(common, *ok[k]) for k in ok}
    mins = {k: (f"{f[k].min():.3f}" if k in f else "n/a") for k in HOLD}

    def med(a, b):
        return f"{np.median(np.abs(f[a] - f[b])):.4f}" if a in f and b in f else "n/a"
    print(f"{label:22s} {mins['kpmf']:>8s} {mins['k05']:>8s} {mins['iag']:>8s} "
          f"{med('kpmf', 'k05'):>9s} {med('kpmf', 'iag'):>9s} {med('k05', 'iag'):>9s}"
          + "".join(f"  [{k}: {v[1]}]" for k, v in g.items() if v[1]))
    rows.append(dict(window=label, centre_A=c, half_width_A=pad,
                     **{f"min_{k}": (round(float(f[k].min()), 4) if k in f else None) for k in HOLD},
                     **{f"median_abs_{a}_minus_{b}": (round(float(np.median(np.abs(f[a] - f[b]))), 4)
                                                     if a in f and b in f else None)
                        for a, b in (("kpmf", "k05"), ("kpmf", "iag"), ("k05", "iag"))},
                     unavailable="; ".join(f"{k}: {v[1]}" for k, v in g.items() if v[1])))

out = ROOT / "data/audit/rya1214_cno_products/kpmf_h2o_flux_test.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"wrote {out.relative_to(ROOT)}")
