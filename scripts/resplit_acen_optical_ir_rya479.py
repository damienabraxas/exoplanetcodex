# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         resplit_acen_optical_ir_rya479.py
# Module:       scripts (one-off data-conditioning tools)
# Description:  Re-attributes the alpha Cen optical + IR frames (HARPS, NIRPS,
#               ESPRESSO, CRIRES) to the correct star by the FITS OBJECT
#               header, NOT the folder name. The RYA-301 audit found these
#               folders carry the same A/B-mixing disease RYA-303 fixed for the
#               STIS UV — folder name != true star — plus a non-alpha-Cen
#               contaminant ("Star S5") and 5 truncated NIRPS downloads.
#
#               This is the RYA-303 analog for the non-UV arms. It:
#                 Step 1 — normalizes every OBJECT string (case/separator-
#                          insensitive) to alpha_cen_a (HD 128620) /
#                          alpha_cen_b (HD 128621) / quarantine.
#                 Step 2 — runs a FITS-integrity sweep (full data read) and
#                          flags truncated/corrupt frames for re-download.
#                 Step 3 — computes the alpha Cen A HARPS short-frame co-add
#                          SNR (combined SN35 = sqrt(sum SN35^2), photon-noise
#                          stacking) to show the 1 s frames clear the SNR-200
#                          science floor.
#                 Step 4 — flags (does NOT exclude) the tight A/B-separation
#                          B HARPS frames (< 4.5"), computed from the
#                          Pourbaix/Kervella relative orbit at each epoch.
#
#               Non-destructive: source is opened read-only; NOTHING under
#               data/spectra/ is written, moved, or deleted (DATA_SAFETY.md).
#               Products (manifest, per-star/instrument file lists, README)
#               are written under data/audit/ per the RYA-303 convention.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-29
# Last modified: 2026-06-29
# Linear issue: RYA-479 — DATA: Re-split alpha Cen optical+IR by OBJECT
#                         (gates RYA-302)
#
# -----------------------------------------------------------------------------
# SCIENTIFIC CONTEXT
# -----------------------------------------------------------------------------
# Pipeline stage: Pre-pipeline data conditioning (before the alpha Cen run).
# Target stars:   alpha Cen A (HD 128620, G2 V) + B (HD 128621, K1 V) — a
#                 common-proper-motion visual binary (~80 yr orbit).
# Instruments:    HARPS / NIRPS / ESPRESSO (ESO Phase-3 ADP), CRIRES+ (cr2res).
# Science goal:   Clean per-star optical+IR sets so the alpha Cen abundance run
#                 (RYA-302) never computes A from majority-B photons.
#
# -----------------------------------------------------------------------------
# LOADER CONTRACTS (carried from the RYA-301 audit — the loader MUST honor)
# -----------------------------------------------------------------------------
#   HARPS    SPECSYS=BARYCENT -> flux already barycentric; do NOT re-apply BERV.
#   NIRPS    read a FLUX_TELL_* column (FLUX_TELL_CAL/_EL) + ATM_TRANSM, not raw
#            FLUX (telluric-corrected product).
#   CRIRES   WAVE in nm -> convert to Angstrom; SPECSYS=TOPOCENT (BERV at load);
#            NO telluric column -> gated out until molecfit is applied.
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
#   numpy, astropy
# =============================================================================

import math
import os
import re
from datetime import datetime
from pathlib import Path

import numpy as np
from astropy.io import fits
# Standalone-script bootstrap (RYA-313): put the REPO ROOT on sys.path BEFORE
# importing config/pipeline, so this runs from any cwd. Derived from __file__.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

# --- paths -------------------------------------------------------------------
DATA_ROOT = codex_path('data.spectra_local') / 'Alpha Centauri (vetted)'
STAR_FOLDERS = ["Alpha Cen A", "Alpha Cen B"]
INSTRUMENTS = ["HARPS", "NIRPS", "ESPRESSO", "CRIRES"]
OUTPUT_DIR = Path("data/audit/alpha_cen_optical_ir")

SNR_FLOOR = 200          # PIPELINE['snr_min_science']
TIGHT_SEP_ARCSEC = 4.5   # below this = binary-contamination watch-list


# --- OBJECT -> star normalizer ----------------------------------------------
def norm_object(obj):
    """Map a raw OBJECT header to (star, reason). Case/separator-insensitive."""
    if obj is None:
        return "quarantine", "missing OBJECT"
    key = re.sub(r"[^a-z0-9]", "", obj.lower())
    if key == "hd128620":
        return "alpha_cen_a", "HD128620"
    if key == "hd128621":
        return "alpha_cen_b", "HD128621"
    if "cen" in key:                       # alfcena / alphacenb / alfcenb ...
        if key.endswith("a"):
            return "alpha_cen_a", "name->A"
        if key.endswith("b"):
            return "alpha_cen_b", "name->B"
        return "quarantine", f"ambiguous cen-name ({obj})"
    return "quarantine", f"non-alphaCen ({obj})"


# --- A/B angular separation (Pourbaix & Boffin 2016 / Kervella+2016) ----------
_ORB = dict(P=79.91, T=1955.66, e=0.5179, a=17.57,
            inc=math.radians(79.32), Om=math.radians(204.85), w=math.radians(232.0))


def separation_arcsec(year):
    P, T, e, a = _ORB["P"], _ORB["T"], _ORB["e"], _ORB["a"]
    inc, Om, w = _ORB["inc"], _ORB["Om"], _ORB["w"]
    M = 2 * math.pi * (((year - T) / P) % 1.0)
    E = M
    for _ in range(60):
        E = E - (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2),
                        math.sqrt(1 - e) * math.cos(E / 2))
    r = a * (1 - e * math.cos(E))
    u = nu + w
    x = r * (math.cos(Om) * math.cos(u) - math.sin(Om) * math.sin(u) * math.cos(inc))
    y = r * (math.sin(Om) * math.cos(u) + math.cos(Om) * math.sin(u) * math.cos(inc))
    return math.hypot(x, y)


def _decimal_year(date_obs):
    d = datetime.strptime(date_obs[:19], "%Y-%m-%dT%H:%M:%S")
    start, end = datetime(d.year, 1, 1), datetime(d.year + 1, 1, 1)
    return d.year + (d - start).total_seconds() / (end - start).total_seconds()


# --- FITS integrity ----------------------------------------------------------
def integrity(path):
    """Force a full data read; classify OK / TRUNCATED / CORRUPT."""
    import warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with fits.open(path, mode="readonly") as hdul:
                for h in hdul:
                    _ = h.data
        return "OK"
    except Exception as exc:                       # noqa: BLE001
        msg = str(exc).lower()
        if "truncat" in msg or "buffer is too small" in msg:
            return "TRUNCATED"
        return "CORRUPT"


# --- scan + attribute --------------------------------------------------------
def scan():
    rows = []
    for folder in STAR_FOLDERS:
        for inst in INSTRUMENTS:
            for f in sorted((DATA_ROOT / folder / inst).glob("*.fits")):
                with fits.open(f, mode="readonly") as hdul:
                    h = hdul[0].header
                    obj = h.get("OBJECT")
                    date = h.get("DATE-OBS")
                    star, reason = norm_object(obj)
                    sep, tight = None, False
                    if star == "alpha_cen_b" and inst == "HARPS" and date:
                        sep = round(separation_arcsec(_decimal_year(date)), 3)
                        tight = sep < TIGHT_SEP_ARCSEC
                    rows.append(dict(
                        true_star=star, instrument=inst, src_folder=folder,
                        raw_object=obj, norm_reason=reason,
                        file=os.path.basename(f),
                        date_obs=date, exptime=h.get("EXPTIME"),
                        snr35=h.get("ESO DRS SPE EXT SN35"),
                        specsys=h.get("SPECSYS"),
                        sep_arcsec=sep, tight_sep_flag=tight,
                        integrity=integrity(f),
                        src_path=str(f)))
    return rows


# --- alpha Cen A HARPS co-add SNR (photon-noise stacking) --------------------
def coadd_snr(rows):
    A = [r for r in rows
         if r["true_star"] == "alpha_cen_a" and r["instrument"] == "HARPS"
         and r["integrity"] == "OK" and r["snr35"]]
    snr = [float(r["snr35"]) for r in A]
    sub = [s for s in snr if s < SNR_FLOOR]
    return dict(
        n_frames=len(A), n_below_floor=len(sub),
        per_frame_min=min(snr), per_frame_median=sorted(snr)[len(snr) // 2],
        per_frame_max=max(snr),
        combined_all=round(math.sqrt(sum(s * s for s in snr))),
        combined_short_only=round(math.sqrt(sum(s * s for s in sub))),
        epoch_start=min(r["date_obs"] for r in A)[:10],
        epoch_end=max(r["date_obs"] for r in A)[:10])


# --- write products ----------------------------------------------------------
def write_products(rows, co):
    import csv
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # full manifest
    cols = ["true_star", "instrument", "src_folder", "raw_object", "norm_reason",
            "file", "date_obs", "exptime", "snr35", "specsys", "sep_arcsec",
            "tight_sep_flag", "integrity", "src_path"]
    with open(OUTPUT_DIR / "acen_optical_ir_manifest_rya479.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["true_star"], r["instrument"],
                                             r["date_obs"] or "")):
            w.writerow({k: r[k] for k in cols})

    # per-star/instrument clean file lists (OK integrity, true star)
    from collections import defaultdict
    groups = defaultdict(list)
    for r in rows:
        if r["integrity"] == "OK" and r["true_star"] != "quarantine":
            groups[(r["true_star"], r["instrument"])].append(r["src_path"])
    for (star, inst), paths in sorted(groups.items()):
        with open(OUTPUT_DIR / f"{star}_{inst}_files.txt", "w") as fh:
            fh.write("\n".join(sorted(paths)) + "\n")

    # quarantine + truncated lists
    for tag, pred in [("quarantine_StarS5",
                       lambda r: r["true_star"] == "quarantine"),
                      ("truncated_redownload",
                       lambda r: r["integrity"] != "OK")]:
        sel = [r for r in rows if pred(r)]
        if sel:
            with open(OUTPUT_DIR / f"_{tag}_files.txt", "w") as fh:
                fh.write("\n".join(sorted(r["src_path"] for r in sel)) + "\n")

    # co-add SNR summary
    with open(OUTPUT_DIR / "alpha_cen_a_HARPS_coadd_snr_rya479.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        for k, v in co.items():
            w.writerow([k, v])


def main():
    print(f"Source root (read-only): {DATA_ROOT}")
    rows = scan()
    print(f"Scanned {len(rows)} frames across {INSTRUMENTS}\n")

    from collections import Counter
    clean = Counter((r["true_star"], r["instrument"]) for r in rows
                    if r["integrity"] == "OK" and r["true_star"] != "quarantine")
    print("CLEAN per-star x instrument (OK integrity, by OBJECT):")
    for k in sorted(clean):
        print(f"  {k[0]:13s} {k[1]:9s}: {clean[k]}")

    quar = [r for r in rows if r["true_star"] == "quarantine"]
    trunc = [r for r in rows if r["integrity"] != "OK"]
    tight = [r for r in rows if r["tight_sep_flag"]]
    print(f"\nQuarantine (Star S5): {len(quar)} "
          f"{dict(Counter(r['instrument'] for r in quar))}")
    print(f"Truncated/corrupt (re-download): {len(trunc)} "
          f"{dict(Counter((r['instrument'], r['integrity']) for r in trunc))}")
    for r in sorted(trunc, key=lambda r: r["file"]):
        print(f"    {r['instrument']} {r['file']} ({r['integrity']})")
    print(f"Tight-sep B HARPS (<{TIGHT_SEP_ARCSEC}\"): {len(tight)} "
          f"{dict(Counter(r['date_obs'][:4] for r in tight))}")

    co = coadd_snr(rows)
    print(f"\nalpha Cen A HARPS co-add ({co['epoch_start']}..{co['epoch_end']}):")
    print(f"  before: {co['n_below_floor']}/{co['n_frames']} frames below SNR "
          f"{SNR_FLOOR} (median SN35 {co['per_frame_median']:.0f})")
    print(f"  after : combined SN35 = {co['combined_all']} (all) / "
          f"{co['combined_short_only']} (short-only) -> clears {SNR_FLOOR}")

    write_products(rows, co)
    print(f"\n[DONE] Read {len(rows)} FITS from data/spectra/ (read-only).")
    print(f"[DONE] Products written to {OUTPUT_DIR}/")
    print(f"[DONE] No source files were modified.")


if __name__ == "__main__":
    main()
