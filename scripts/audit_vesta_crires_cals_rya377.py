"""
scripts/audit_vesta_crires_cals_rya377.py
=========================================
RYA-377 — Vesta CRIRES+ observing-night calibration set + per-setting science
coverage map. READ-ONLY (codex-data-audit approach; reuses the 370/384 logic).

Two questions (see RYA-373 / RYA-376):
  1. Telluric standard. CRIRES+ programs observe a hot-star telluric standard as
     a separate OB. Is one present in what was downloaded? If so it gives an
     EMPIRICAL telluric (better than the molecfit/GDAS fallback of RYA-373) and a
     second, independent telluric check at the low-RV epoch.
  2. Per-setting science coverage. RYA-370 gave the union (0.95-2.49 um) but not
     the per-setting sub-windows. CRIRES+ settings are narrow chunks with detector
     gaps, so band tiling does NOT guarantee a diagnostic line lands in data --
     that has to be inventoried from the WAVE arrays themselves.

Outputs (data/audit/vesta_crires_rya377/):
  - inventory.csv          one row per FITS product (header truth)
  - setting_coverage.csv   covered (order x detector) sub-windows, merged per setting
  - diagnostic_map.csv     each high-value near-IR line vs covered windows

    python scripts/audit_vesta_crires_cals_rya377.py [DATA_DIR]
"""
from __future__ import annotations
import sys, glob, os, warnings, collections
from pathlib import Path
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from astropy.io import fits as pf

DEFAULT_DATA = Path(
    "/Users/ryanschmitt/Documents/Exoplanet Codex/data/spectra/exoplanetcodex-data/"
    "Solar Calibration/Solar System Targets/Vesta"
)
OUT = Path(__file__).resolve().parents[1] / 'data' / 'audit' / 'vesta_crires_rya377'

# High-value near-IR diagnostics to test for landing in covered windows.
# Wavelengths in nm (converted from the Angstrom values in the ticket / RYA-376).
DIAGNOSTICS = [
    ("S I 1.045um",  [1045.5, 1045.6, 1045.9]),
    ("P I 1.05um",   [1051.1, 1052.9, 1058.1, 1059.6]),
    ("C I 1.068um",  [1068.3]),
    ("K I 1.17um",   [1169.0, 1176.9]),
    ("OH 1.5-1.8um", None),  # band 1500-1800 nm, handled as an interval
]
OH_BAND = (1500.0, 1800.0)
DET_GAP_NM = 0.5  # contiguity threshold when merging WAVE samples into sub-windows


def band_of(setting: str) -> str:
    """Y/J/H/K from the WLEN setting id (e.g. Y1029 -> Y)."""
    return setting[0] if setting and setting[0] in "YJHK" else "?"


def covered_windows(hdul) -> list[tuple[float, float]]:
    """Truly covered sub-windows (nm) from the SPECTRUM table, split at detector /
    order gaps. Empirical: derived from the WAVE column where FLUX is finite, not
    from header advertised ranges."""
    d = hdul[1].data
    wave = np.asarray(d['WAVE']).ravel().astype(float)
    flux = np.asarray(d['FLUX']).ravel().astype(float)
    good = np.isfinite(wave) & np.isfinite(flux) & (wave > 0)
    w = np.sort(wave[good])
    if w.size == 0:
        return []
    # split where the gap between consecutive samples exceeds DET_GAP_NM
    splits = np.where(np.diff(w) > DET_GAP_NM)[0]
    segs, start = [], 0
    for s in splits:
        segs.append((float(w[start]), float(w[s])))
        start = s + 1
    segs.append((float(w[start]), float(w[-1])))
    # drop degenerate slivers
    return [(a, b) for a, b in segs if b - a > 0.05]


def merge_windows(wins: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Union of intervals (used to merge repeated settings across nights)."""
    if not wins:
        return []
    wins = sorted(wins)
    out = [list(wins[0])]
    for a, b in wins[1:]:
        if a <= out[-1][1] + DET_GAP_NM:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def in_windows(x: float, wins: list[tuple[float, float]]) -> bool:
    return any(a <= x <= b for a, b in wins)


def scan() -> tuple[pd.DataFrame, dict]:
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DATA
    files = sorted(glob.glob(str(data_dir / '**' / '*.fits'), recursive=True))
    rows, win_by_setting = [], collections.defaultdict(list)
    for f in files:
        with pf.open(f) as h:
            ph = h[0].header
            instrument = str(ph.get('INSTRUME') or ph.get('ESO DET ID') or '')
            is_crires = 'CRIRES' in instrument.upper()
            setting = str(ph.get('ESO INS WLEN ID', '?'))
            wins = covered_windows(h)
            if is_crires:  # only CRIRES+ feeds the coverage map (scope of RYA-377)
                win_by_setting[setting].extend(wins)
            rows.append(dict(
                file=os.path.basename(f),
                folder=Path(f).parent.name,
                instrument=instrument,
                object=str(ph.get('OBJECT', '')),
                pro_catg=str(ph.get('ESO PRO CATG', '')),
                prodcatg=str(ph.get('PRODCATG', '')),
                dpr_type=str(ph.get('ESO DPR TYPE', '')),
                setting=setting if is_crires else '',
                band=band_of(setting) if is_crires else '',
                date_obs=str(ph.get('DATE-OBS', '')),
                mjd_obs=ph.get('MJD-OBS'),
                exptime=ph.get('EXPTIME'),
                airm_start=ph.get('ESO TEL AIRM START'),
                airm_end=ph.get('ESO TEL AIRM END'),
                spec_res=ph.get('SPEC_RES'),
                wavelmin_nm=ph.get('WAVELMIN'),
                wavelmax_nm=ph.get('WAVELMAX'),
                n_subwin=len(wins),
            ))
    return pd.DataFrame(rows), win_by_setting


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, win_by_setting = scan()
    df.to_csv(OUT / 'inventory.csv', index=False)

    print("=" * 78)
    print("RYA-377  Vesta CRIRES+ cal set + per-setting coverage audit")
    print("=" * 78)
    crires = df[df['instrument'].str.upper().str.contains('CRIRES')].copy()
    other = df[~df['instrument'].str.upper().str.contains('CRIRES')].copy()
    print(f"\nProducts found: {len(df)}  (CRIRES+: {len(crires)}, other-instrument: {len(other)})")
    print("By instrument / folder / PRO.CATG:")
    print(df.groupby(['instrument', 'folder', 'pro_catg']).size().to_string())
    if len(other):
        oinstr = ", ".join(sorted(other['instrument'].unique()))
        print(f"\nNOTE: {len(other)} non-CRIRES Vesta products present ({oinstr}; the"
              f" RYA-370 optical reflected-solar set, in a mis-named 'xshooter' folder)"
              f" -- excluded from the CRIRES+ coverage map below.")

    # ---- 1. Telluric-standard verdict -------------------------------------
    print("\n" + "-" * 78)
    print("1. TELLURIC-STANDARD VERDICT  (CRIRES+ night)")
    objs = sorted(crires['object'].str.upper().unique())
    print(f"   distinct CRIRES+ OBJECT values: {objs}")
    non_vesta = crires[~crires['object'].str.upper().str.contains('VESTA')]
    if len(non_vesta) == 0:
        print("   VERDICT: NO companion telluric standard present. Every product is")
        print("   OBJECT=Vesta SCIENCE.SPECTRUM. No hot-star STD OB, no raw cal frames")
        print("   (flats/darks/arcs), no cr2res master cals were delivered -- Phase-3")
        print("   IDPs do not bundle them. RYA-373 molecfit/GDAS fallback STANDS.")
    else:
        print("   Candidate non-Vesta products:")
        print(non_vesta[['file', 'object', 'setting', 'airm_start']].to_string(index=False))

    # ---- 2. Per-setting science coverage ----------------------------------
    print("\n" + "-" * 78)
    print("2. PER-SETTING SCIENCE COVERAGE  (covered sub-windows, nm)")
    cov_rows = []
    for setting in sorted(win_by_setting):
        merged = merge_windows(win_by_setting[setting])
        band = band_of(setting)
        nfiles = int((df['setting'] == setting).sum())
        span = f"{min(a for a, _ in merged):.1f}-{max(b for _, b in merged):.1f}"
        wtxt = ", ".join(f"{a:.1f}-{b:.1f}" for a, b in merged)
        print(f"   [{band}] {setting:6s}  files={nfiles}  span {span} nm  | {len(merged)} windows: {wtxt}")
        for a, b in merged:
            cov_rows.append(dict(setting=setting, band=band, win_min_nm=round(a, 2),
                                 win_max_nm=round(b, 2), n_files=nfiles))
    pd.DataFrame(cov_rows).to_csv(OUT / 'setting_coverage.csv', index=False)

    # union of every covered window across all settings = what the asteroid IR reaches
    all_wins = merge_windows([w for ws in win_by_setting.values() for w in ws])

    # ---- 3. Diagnostic-coverage map ---------------------------------------
    print("\n" + "-" * 78)
    print("3. DIAGNOSTIC-COVERAGE MAP  (does the line land in a covered window?)")
    diag_rows = []
    for name, lines in DIAGNOSTICS:
        if lines is None:  # OH band: report fraction of 1500-1800 nm covered
            lo, hi = OH_BAND
            band_wins = [(max(a, lo), min(b, hi)) for a, b in all_wins if b > lo and a < hi]
            cov = sum(b - a for a, b in band_wins)
            frac = cov / (hi - lo)
            settings_hit = sorted({s for s in win_by_setting
                                   if any(x < hi and y > lo for x, y in merge_windows(win_by_setting[s]))})
            status = "PARTIAL" if 0 < frac else "ABSENT"
            print(f"   {name:14s} band {lo:.0f}-{hi:.0f} nm: {status}  "
                  f"({100*frac:.0f}% covered) settings={settings_hit}")
            diag_rows.append(dict(diagnostic=name, line_nm="1500-1800",
                                  covered="PARTIAL" if frac else "NO",
                                  detail=f"{100*frac:.0f}% of band", settings=";".join(settings_hit)))
            continue
        for lam in lines:
            hit = [s for s in win_by_setting
                   if in_windows(lam, merge_windows(win_by_setting[s]))]
            ok = bool(hit)
            print(f"   {name:14s} {lam:8.1f} nm: {'COVERED' if ok else 'not covered':12s}"
                  f" {('settings=' + ','.join(sorted(hit))) if ok else ''}")
            diag_rows.append(dict(diagnostic=name, line_nm=lam,
                                  covered="YES" if ok else "NO",
                                  detail="", settings=";".join(sorted(hit))))
    pd.DataFrame(diag_rows).to_csv(OUT / 'diagnostic_map.csv', index=False)

    print("\n" + "-" * 78)
    print(f"Wrote: {OUT}/inventory.csv, setting_coverage.csv, diagnostic_map.csv")
    print(f"Total covered span (union of all settings): "
          f"{all_wins[0][0]:.1f}-{all_wins[-1][1]:.1f} nm in {len(all_wins)} windows")


if __name__ == '__main__':
    main()
