"""RYA-929 three-way telluric diagnostic: 1984 KPNO vs Kurucz 2005 vs IAG.

This emits window diagnostics only. It never replaces a science product or reports an
abundance. Kurucz irradiance is locally continuum-normalized with a broad max filter;
the result is explicitly a smoke comparison until the published continuum prescription
and wavelength medium are fully resolved.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
from pipeline.wavelength_util import vac_to_air


WINDOWS = (
    ("O2_B", 6867.0, 6884.0),
    ("H2O_7160_7340", 7160.0, 7340.0),
    ("O2_A", 7594.0, 7685.0),
    ("H2O_8100_8400", 8100.0, 8400.0),
    ("H2O_9280_9600", 9280.0, 9600.0),
    ("clean_6600_6650", 6600.0, 6650.0),
    ("clean_6690_6702", 6690.0, 6702.0),
    ("clean_7690_7705", 7690.0, 7705.0),
)
LINES = (
    ("Al_I_6631", 6631.218), ("Al_I_6696", 6696.185),
    ("Fe_I_6872", 6872.1617), ("Fe_I_6875", 6875.44503),
    ("Si_I_6876", 6876.35917), ("Fe_I_6881", 6881.44175),
    ("K_I_7665", 7665.0), ("K_I_7699", 7698.964), ("N_I_8216", 8216.0),
)


def read_ascii_windows(path: Path, windows, scale=1.0, skip=0):
    out = {name: ([], []) for name, _, _ in windows}
    with path.open(errors="replace") as fh:
        for _ in range(skip): next(fh, "")
        for line in fh:
            q = line.split()
            if len(q) < 2: continue
            try: w, f = float(q[0]) * scale, float(q[1])
            except ValueError: continue
            for name, lo, hi in windows:
                if lo <= w <= hi: out[name][0].append(w); out[name][1].append(f)
    return {k: (np.asarray(v[0]), np.asarray(v[1])) for k, v in out.items()}


def read_kpno(directory: Path, windows):
    out = {name: ([], []) for name, _, _ in windows}
    for p in sorted(directory.glob("lm[0-9]*")):
        chunk = read_ascii_windows(p, windows, scale=10.0)
        for name in out:
            out[name][0].extend(chunk[name][0]); out[name][1].extend(chunk[name][1])
    return {k: (np.asarray(v[0]), np.asarray(v[1])) for k, v in out.items()}


def read_iag(path: Path, windows):
    out = {name: ([], []) for name, _, _ in windows}
    with gzip.open(path, "rb") as fh:
        for raw in fh:
            q = raw.split()
            if len(q) < 2: continue
            try: wa, f = 1e8 / float(q[0]), float(q[1]); w = float(vac_to_air(np.array([wa]))[0])
            except (ValueError, TypeError): continue
            for name, lo, hi in windows:
                if lo <= w <= hi: out[name][0].append(w); out[name][1].append(f)
    return {k: (np.asarray(v[0]), np.asarray(v[1])) for k, v in out.items()}


def read_iag_fits(path: Path, windows):
    from astropy.io import fits
    out = {name: ([], []) for name, _, _ in windows}
    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        wn = np.asarray(d["v"], float)
        fl = np.asarray(d["s"], float)
    w = np.asarray(vac_to_air(1e8 / wn))
    for name, lo, hi in windows:
        m = (w >= lo) & (w <= hi) & np.isfinite(fl)
        out[name] = (w[m], fl[m])
    return out


def interp_grid(w, f, lo, hi):
    m = (w >= lo) & (w <= hi) & np.isfinite(f)
    return w[m], f[m]


def normalized_kurucz(w, f):
    # Streaming-safe preliminary continuum: robust high percentile per diagnostic
    # window. Formal continuum prescription remains part of method intake.
    return f / np.nanpercentile(f, 99.5)


def metrics(w, f, lo, hi, normalized=False):
    x, y = interp_grid(w, f, lo, hi)
    if len(x) < 10:
        return {"n": len(x), "fraction_below_0p5": np.nan, "min_flux": np.nan}
    if not normalized:
        y = normalized_kurucz(x, y)
    return {"n": int(len(y)),
            "fraction_below_0p5": float(np.mean(y < 0.5)),
            "min_flux": float(np.nanmin(y))}


def line_metric(w, f, center, normalized=False):
    m = (w >= center - 0.5) & (w <= center + 0.5) & np.isfinite(f)
    if m.sum() < 5:
        return {"n": int(m.sum()), "min_wave_A": np.nan, "depth": np.nan}
    x, y = w[m], f[m]
    if not normalized:
        y = y / np.nanpercentile(y, 99.5)
    j = int(np.nanargmin(y))
    return {"n": int(len(y)), "min_wave_A": float(x[j]),
            "depth": float(1.0 - y[j])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kpno", type=Path, required=True)
    ap.add_argument("--kurucz", type=Path, required=True)
    ap.add_argument("--iag", type=Path, required=True)
    ap.add_argument("--iag-fits", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--lines-out", type=Path, required=True)
    args = ap.parse_args()
    data = {"kpno_1984": read_kpno(args.kpno, WINDOWS),
            "kurucz_2005": read_ascii_windows(args.kurucz, WINDOWS, scale=10.0, skip=9),
            "iag": (read_iag_fits(args.iag, WINDOWS) if args.iag_fits
                    else read_iag(args.iag, WINDOWS))}
    rows = []
    for name, lo, hi in WINDOWS:
        row = {"region": name, "lo_A": lo, "hi_A": hi,
               "status": "DIAGNOSTIC_ONLY",
               "note": "Kurucz local 5-A max-filter continuum; not an abundance product"}
        for label, by_window in data.items():
            w, f = by_window[name]
            m = metrics(w, f, lo, hi, normalized=(label != "kurucz_2005"))
            row[f"{label}_n"] = m["n"]
            row[f"{label}_fraction_below_0p5"] = m["fraction_below_0p5"]
            row[f"{label}_min_flux"] = m["min_flux"]
        rows.append(row)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    for row in rows:
        print(row)
    line_rows = []
    for name, center in LINES:
        row = {"line": name, "center_A": center, "status": "DIAGNOSTIC_ONLY"}
        for label, by_window in data.items():
            w_all = np.concatenate([v[0] for v in by_window.values()])
            f_all = np.concatenate([v[1] for v in by_window.values()])
            m = line_metric(w_all, f_all, center, normalized=(label != "kurucz_2005"))
            row[f"{label}_n"] = m["n"]
            row[f"{label}_min_wave_A"] = m["min_wave_A"]
            row[f"{label}_depth"] = m["depth"]
        line_rows.append(row)
    args.lines_out.parent.mkdir(parents=True, exist_ok=True)
    with args.lines_out.open("w", newline="") as fh:
        fields = list(line_rows[0]); writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader(); writer.writerows(line_rows)
    for row in line_rows: print(row)


if __name__ == "__main__":
    main()
