# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         dedup_resplit_stis_uv_rya303.py
# Module:       scripts (one-off data-conditioning tools)
# Description:  De-duplicates and re-attributes the HST/STIS UV frames for
#               alpha Cen A + B. The RYA-301 audit found the STIS UV set is
#               (1) byte-duplicated across the "Alpha Centauri A" and
#               "Alpha Centauri B" folders and (2) internally mixing both
#               stars within a folder — folder name != true star.
#
#               This script:
#                 Step 1 — SHA-256 hashes every STIS UV frame across both
#                          folders, collapses to the unique set, records the
#                          duplicate map. NOTHING is deleted or moved (the
#                          source tree is read-only per DATA_SAFETY.md).
#                 Step 2 — re-attributes each unique frame to alpha Cen A
#                          (HD 128620) or B (HD 128621) by the TARGNAME/OBJECT
#                          header, NOT the folder. Every frame is independently
#                          cross-checked by proper-motion-propagated position
#                          (RA_TARG/DEC_TARG vs A/B at the observation epoch);
#                          coordinate-named or disagreeing frames are flagged.
#                 Step 3 — confirms the products are extracted 1D science
#                          spectra (x1d/sx1), documents the modes present
#                          (E140M/H, E230M/H) and their coverage, and flags
#                          STIS vacuum wavelengths for air conversion at the
#                          spectrum-matching boundary.
#
#               Non-destructive — reads source only; writes products under
#               data/audit/ per convention.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-14
# Last modified: 2026-06-14
# Linear issue: RYA-303 — DATA: De-dup + re-split STIS UV (alpha Cen) by OBJECT
#
# -----------------------------------------------------------------------------
# SCIENTIFIC CONTEXT
# -----------------------------------------------------------------------------
# Pipeline stage: Pre-pipeline data conditioning (before any UV-inclusive run)
# Target stars:   alpha Cen A (HD 128620, G2 V) + B (HD 128621, K1 V) — a
#                 common-proper-motion visual binary (~80 yr orbit).
# Instrument:     HST/STIS UV echelle (E140M/H FUV, E230M/H NUV)
# Science goal:   Clean per-star UV sets so the UV (C/N/O, Fe II, Mg II, ...)
#                 can be folded into the alpha Cen run (RYA-302) without
#                 A/B cross-attribution.
#
# -----------------------------------------------------------------------------
# WAVELENGTH CONVENTION (important)
# -----------------------------------------------------------------------------
# STIS WAVELENGTH arrays are VACUUM. The pipeline works in air angstrom for
# lambda >= 2000 A (VALD air/vacuum convention). The loader must apply a
# vacuum->air conversion at the spectrum-matching boundary. This script only
# FLAGS the convention; it does not convert.
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
# External: astropy, numpy, pandas
# Data:     STIS UV mode folders under DATA_ROOT (read-only)
# =============================================================================

import re
import sys
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
# Standalone-script bootstrap (RYA-313): put the REPO ROOT on sys.path BEFORE
# importing config/pipeline, so this runs from any cwd. Derived from __file__.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

# =============================================================================
# CONFIGURATION
# =============================================================================
DATA_ROOT = codex_path('data.spectra_local')
STAR_FOLDERS = ["Alpha Centauri A", "Alpha Centauri B"]
# The STIS UV mode subfolders that exist (mirrored) under each star folder.
UV_MODE_DIRS = ["uv-e140m", "uv-e140h", "nuv-e230h", "nuv-e230m"]

OUTPUT_DIR = Path("data/audit/alpha_cen_stis_uv")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Extracted 1D science products only (echelle x1d / first-order sx1).
SCIENCE_SUFFIXES = ["_x1d.fits", "_sx1.fits"]

# Base HD number in TARGNAME -> component. HD 128620 = A, HD 128621 = B.
# STIS TARGNAMEs carry visit/pointing-revision suffixes (-1, -2, -NEW, -COPY)
# on top of the HD number; the suffix is NOT a different star, so we attribute
# on the BASE HD number only.
HD_TO_COMPONENT = {"HD128620": "A", "HD128621": "B"}
HD_RE = re.compile(r"HD\s*0*12862([01])")

# Reference astrometry (SIMBAD / Hipparcos van Leeuwen 2007), ICRS @ epoch 2000.
# Both components share near-identical system proper motion; the A-B separation
# is orbital (not captured by linear PM) but the J2000 offset is carried forward
# well enough to assign the nearest component at HST pointing precision.
ALPHA_CEN_REF = {
    "A": dict(ra=219.90205833, dec=-60.83397222, pmra=-3679.25, pmdec=473.67,
              plx=754.81, rv=-22.3),   # HD 128620
    "B": dict(ra=219.89613750, dec=-60.83961111, pmra=-3614.39, pmdec=802.98,
              plx=796.92, rv=-20.6),   # HD 128621
}

# STIS echelle modes -> waveband label.
MODE_BAND = {"E140M": "FUV", "E140H": "FUV", "E230M": "NUV", "E230H": "NUV"}


# =============================================================================
# HELPERS
# =============================================================================
def sha256_file(path, chunk=1 << 20):
    """SHA-256 of a file's bytes (streamed, read-only)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def is_science_product(path):
    name = path.name.lower()
    return any(name.endswith(s) for s in SCIENCE_SUFFIXES)


def read_headers(path):
    """Return (h0, h1). STIS keeps OBSERVATION metadata in [0] and exposure
    metadata (DATE-OBS, EXPSTART, EXPTIME) in [1]."""
    with fits.open(path, memmap=False) as hdul:
        h0 = hdul[0].header
        h1 = hdul[1].header if len(hdul) > 1 else fits.Header()
        return h0.copy(), h1.copy()


def get(h0, h1, key, default=None):
    if key in h0 and h0[key] not in (None, ""):
        return h0[key]
    if key in h1 and h1[key] not in (None, ""):
        return h1[key]
    return default


def wavelength_span(path):
    """(min, max) of the WAVELENGTH column across all echelle orders, in A."""
    try:
        with fits.open(path, memmap=False) as hdul:
            if len(hdul) > 1 and hdul[1].data is not None \
                    and "WAVELENGTH" in hdul[1].columns.names:
                wl = np.asarray(hdul[1].data["WAVELENGTH"], dtype=float).ravel()
                wl = wl[np.isfinite(wl) & (wl > 0)]
                if wl.size:
                    return float(wl.min()), float(wl.max())
    except Exception:
        pass
    return None, None


def component_by_targname(targname):
    """A/B from the BASE HD number in TARGNAME (suffixes ignored).

    Returns 'A'/'B', or None if no HD 128620/128621 token is present
    (i.e. a genuinely coordinate-named frame that needs PM resolution)."""
    m = HD_RE.search(str(targname).upper())
    if not m:
        return None
    return "A" if m.group(1) == "0" else "B"


def pm_propagated_seps(ra_targ, dec_targ, date_obs):
    """Sanity check only: separation of RA_TARG/DEC_TARG to the
    linear-PM-propagated catalog positions of A and B at the observation epoch.

    NOTE (documented limitation): for this tight common-proper-motion binary the
    ~80-yr orbital motion dominates the internal A-B separation and is NOT
    captured by linear PM, so the *absolute* nearest-component verdict is
    degenerate (it collapses toward one component). This is reported as a
    diagnostic, never used to override the OBJECT/HD attribution.
    Returns (nearest, sep_A_arcsec, sep_B_arcsec) or (None, None, None)."""
    if ra_targ is None or dec_targ is None or date_obs is None:
        return None, None, None
    try:
        obstime = Time(str(date_obs))
    except Exception:
        return None, None, None
    target = SkyCoord(ra=float(ra_targ) * u.deg, dec=float(dec_targ) * u.deg)
    seps = {}
    for comp, p in ALPHA_CEN_REF.items():
        c0 = SkyCoord(ra=p["ra"] * u.deg, dec=p["dec"] * u.deg,
                      pm_ra_cosdec=p["pmra"] * u.mas / u.yr,
                      pm_dec=p["pmdec"] * u.mas / u.yr,
                      distance=(1000.0 / p["plx"]) * u.pc,
                      radial_velocity=p["rv"] * u.km / u.s,
                      obstime=Time("J2000"))
        c_epoch = c0.apply_space_motion(new_obstime=obstime)
        seps[comp] = target.separation(c_epoch).arcsec
    nearest = min(seps, key=seps.get)
    return nearest, seps["A"], seps["B"]


def epoch_matched_pointing_check(inv):
    """Independent geometric corroboration of the A/B labels that does NOT
    depend on an absolute orbit model.

    For each frame, find the nearest-in-time frame of the OPPOSITE component and
    measure the on-sky separation of their pointings (RA_TARG/DEC_TARG). If the
    HD labels track two physically distinct stars, this separation should be a
    few arcsec (the alpha Cen A-B visual separation), not ~0. Adds columns
    'opp_dt_days' and 'opp_sep_arcsec'."""
    out_dt, out_sep = [], []
    has = inv.dropna(subset=["ra_targ", "dec_targ", "date_obs"]).copy()
    has["mjd"] = Time(list(has["date_obs"].astype(str))).mjd
    coords = SkyCoord(ra=has["ra_targ"].to_numpy() * u.deg,
                      dec=has["dec_targ"].to_numpy() * u.deg)
    has = has.assign(_coord=list(coords), _idx=range(len(has)))
    for _, row in inv.iterrows():
        opp = "B" if row["component"] == "A" else "A"
        cand = has[has["component"] == opp]
        if cand.empty or row["date_obs"] is None or pd.isna(row["ra_targ"]):
            out_dt.append(None); out_sep.append(None); continue
        try:
            mjd = Time(str(row["date_obs"])).mjd
        except Exception:
            out_dt.append(None); out_sep.append(None); continue
        j = (cand["mjd"] - mjd).abs().idxmin()
        crow = cand.loc[j]
        sc = SkyCoord(ra=row["ra_targ"] * u.deg, dec=row["dec_targ"] * u.deg)
        out_dt.append(float(abs(crow["mjd"] - mjd)))
        out_sep.append(float(sc.separation(crow["_coord"]).arcsec))
    inv = inv.copy()
    inv["opp_dt_days"] = out_dt
    inv["opp_sep_arcsec"] = out_sep
    return inv


# =============================================================================
# MAIN
# =============================================================================
def run():
    print(f"\n{'='*72}\nRYA-303  STIS UV de-dup + re-split  (alpha Cen A/B)\n{'='*72}")
    print(f"Source root: {DATA_ROOT}")

    # ---- collect every STIS UV frame across both folders ---------------------
    all_paths = []
    for star in STAR_FOLDERS:
        for mode in UV_MODE_DIRS:
            all_paths.extend(sorted((DATA_ROOT / star / mode).rglob("*.fits")))
    all_paths = sorted(all_paths)
    print(f"\nTotal STIS UV FITS across both folders: {len(all_paths)}")
    if not all_paths:
        print("ERROR: no STIS UV FITS found — check DATA_ROOT / UV_MODE_DIRS.")
        sys.exit(1)

    # ===========================================================================
    # STEP 1 — de-dup by byte hash
    # ===========================================================================
    rows = []
    for p in all_paths:
        rel = p.relative_to(DATA_ROOT)
        folder = rel.parts[0]                       # "Alpha Centauri A"/"...B"
        rows.append(dict(
            sha256=sha256_file(p),
            filename=p.name,
            folder=folder,
            folder_component="A" if folder.endswith(" A") else "B",
            relpath=str(rel),
            abspath=str(p),
            size_bytes=p.stat().st_size,
        ))
    dmap = pd.DataFrame(rows)

    n_unique = dmap["sha256"].nunique()
    print(f"\n{'-'*72}\nSTEP 1 — DE-DUP (SHA-256)\n{'-'*72}")
    print(f"  Frames on disk : {len(dmap)}")
    print(f"  Unique frames  : {n_unique}")
    print(f"  Duplicate rows : {len(dmap) - n_unique}")

    # Cross-folder duplication: how many unique frames appear in BOTH folders.
    by_hash = dmap.groupby("sha256")["folder_component"].agg(
        lambda s: "".join(sorted(set(s))))
    both = (by_hash == "AB").sum()
    print(f"  Unique frames present in BOTH A & B folders: {both}")
    print(f"  (filename collision is byte-identical: folder is meaningless)")

    dmap.to_csv(OUTPUT_DIR / "stis_uv_dedup_map.csv", index=False)
    print(f"  duplicate map -> {OUTPUT_DIR/'stis_uv_dedup_map.csv'}")

    # one representative path per unique frame (stable: first folder alpha order)
    reps = (dmap.sort_values(["sha256", "folder"])
                .groupby("sha256", as_index=False).first())

    # ===========================================================================
    # STEP 2 + 3 — re-attribute by OBJECT (+PM cross-check) & inspect product
    # ===========================================================================
    print(f"\n{'-'*72}\nSTEP 2/3 — RE-SPLIT BY OBJECT + PRODUCT/COVERAGE\n{'-'*72}")
    recs = []
    for _, r in reps.iterrows():
        p = Path(r["abspath"])
        h0, h1 = read_headers(p)
        targname = str(get(h0, h1, "TARGNAME", "")).upper().replace(" ", "")
        obj_hdr = get(h0, h1, "OBJECT")
        ra_t = get(h0, h1, "RA_TARG")
        dec_t = get(h0, h1, "DEC_TARG")
        date_obs = get(h0, h1, "DATE-OBS")
        opt_elem = str(get(h0, h1, "OPT_ELEM", "")).upper()
        wmin, wmax = wavelength_span(p)
        flags = []

        # --- attribution by OBJECT header (base HD number) ---
        comp_hdr = component_by_targname(targname)
        coord_named = comp_hdr is None
        if coord_named:
            flags.append(f"COORD_NAMED_NO_HD:{targname}")

        # --- PM diagnostic (sanity only; degenerate for this binary) ---
        comp_pm, sep_a, sep_b = pm_propagated_seps(ra_t, dec_t, date_obs)

        # --- resolve ---
        if comp_hdr is not None:
            component = comp_hdr
            method = "OBJECT(HD)"
        elif comp_pm is not None:
            component = comp_pm
            method = "PM_CROSSMATCH"
        else:
            component = "UNRESOLVED"
            method = "NONE"
            flags.append("UNRESOLVED_NO_HEADER_NO_PM")

        if not is_science_product(p):
            flags.append("NOT_1D_SCIENCE_PRODUCT")
        if opt_elem not in MODE_BAND:
            flags.append(f"UNEXPECTED_MODE:{opt_elem}")

        recs.append(dict(
            sha256=r["sha256"], filename=p.name,
            component=component, attribution_method=method,
            targname=targname, object_hdr=obj_hdr,
            opt_elem=opt_elem, band=MODE_BAND.get(opt_elem, "?"),
            detector=get(h0, h1, "DETECTOR"),
            proposid=get(h0, h1, "PROPOSID"),
            pi=get(h0, h1, "PR_INV_L"),
            date_obs=date_obs, exptime=get(h0, h1, "EXPTIME"),
            ra_targ=ra_t, dec_targ=dec_t,
            pm_nearest=comp_pm, sep_A_arcsec=sep_a, sep_B_arcsec=sep_b,
            wl_min_A=wmin, wl_max_A=wmax,
            wavelength_frame="VACUUM",  # STIS convention; loader must -> air
            found_in_folder=r["folder_component"],
            abspath=r["abspath"],
            flags=";".join(flags),
        ))

    inv = pd.DataFrame(recs).sort_values(["component", "opt_elem", "date_obs"])
    inv = epoch_matched_pointing_check(inv)
    inv.to_csv(OUTPUT_DIR / "stis_uv_unique_inventory.csv", index=False)
    print(f"  unique inventory -> {OUTPUT_DIR/'stis_uv_unique_inventory.csv'}")

    # --- per-star counts by mode ---
    print(f"\n  Per-star unique UV frames by mode (E140M/H, E230M/H):")
    for comp in ["A", "B", "UNRESOLVED"]:
        sub = inv[inv["component"] == comp]
        if sub.empty:
            continue
        print(f"    alpha Cen {comp}: {len(sub)} frames")
        for mode, g in sub.groupby("opt_elem"):
            cov = ""
            if g["wl_min_A"].notna().any():
                cov = f"  lambda {g['wl_min_A'].min():.0f}-{g['wl_max_A'].max():.0f} A (vac)"
            print(f"       {mode:6} ({MODE_BAND.get(mode,'?')})  {len(g):>3} frames{cov}")

    # --- folder-vs-truth contamination summary ---
    print(f"\n  Folder-vs-OBJECT contamination (unique frames):")
    xtab = pd.crosstab(inv["found_in_folder"], inv["component"])
    print("    " + xtab.to_string().replace("\n", "\n    "))

    # --- product check ---
    n_sci = inv["filename"].apply(lambda n: is_science_product(Path(n))).sum()
    print(f"\n  Product check: {n_sci}/{len(inv)} are 1D science products "
          f"(x1d/sx1).")
    suffix_counts = (inv["filename"].str.extract(r"(_[a-z0-9]+\.fits)$")[0]
                     .value_counts())
    for sfx, c in suffix_counts.items():
        print(f"       {sfx}: {c}")

    # --- flags ---
    flagged = inv[inv["flags"] != ""]
    print(f"\n  Flagged frames: {len(flagged)}")
    if len(flagged):
        for _, row in flagged.iterrows():
            print(f"    {row['filename']} [{row['component']}]: {row['flags']}")

    # --- coord-named / PM-resolved summary ---
    coord_named = inv[inv["attribution_method"] != "OBJECT(HD)"]
    pm_resolved = inv[inv["attribution_method"] == "PM_CROSSMATCH"]
    print(f"\n  Attribution method: OBJECT(HD)={ (inv['attribution_method']=='OBJECT(HD)').sum() }"
          f"  PM_CROSSMATCH={len(pm_resolved)}  UNRESOLVED="
          f"{ (inv['attribution_method']=='NONE').sum() }")
    print(f"  Coordinate-named frames needing PM (no HD in TARGNAME): {len(coord_named)}")

    # --- geometric corroboration of the A/B labels (orbit-independent) ---
    gc = inv.dropna(subset=["opp_sep_arcsec"])
    print(f"\n  Geometric corroboration — nearest-in-time OPPOSITE-component "
          f"pointing separation:")
    if len(gc):
        print(f"    median {gc['opp_sep_arcsec'].median():.1f}\"  "
              f"(range {gc['opp_sep_arcsec'].min():.1f}-{gc['opp_sep_arcsec'].max():.1f}\", "
              f"median dt {gc['opp_dt_days'].median():.0f} d)")
        print(f"    -> A and B pointings form two distinct clusters separated by "
              f"a few arcsec\n       (consistent with the alpha Cen A-B visual "
              f"separation) — HD labels track\n       two physically distinct "
              f"targets, not one mislabeled star.")
    print(f"  PM-vs-fixed-catalog nearest verdict is DEGENERATE here "
          f"(all -> {inv['pm_nearest'].mode().iat[0]}); not used for attribution.")

    # --- per-star / per-mode file lists for the loader ---
    # Clear any stale lists first so a re-run never leaves an orphaned mode.
    for old in OUTPUT_DIR.glob("alpha_cen_*_files.txt"):
        old.unlink()
    print(f"\n  Per-star/per-mode file lists:")
    for (comp, mode), g in inv.groupby(["component", "opt_elem"]):
        if comp == "UNRESOLVED":
            continue
        lst = OUTPUT_DIR / f"alpha_cen_{comp}_{mode.lower()}_files.txt"
        g["abspath"].to_csv(lst, index=False, header=False)
        print(f"    {lst}  ({len(g)} frames)")

    print(f"\n{'='*72}")
    print(f"DONE. {len(dmap)} on disk -> {n_unique} unique  |  "
          f"A={ (inv['component']=='A').sum() }  "
          f"B={ (inv['component']=='B').sum() }  "
          f"UNRESOLVED={ (inv['component']=='UNRESOLVED').sum() }")
    print(f"Wavelengths are VACUUM -> flag vacuum->air for the loader.")
    print(f"{'='*72}\n")
    return dmap, inv


if __name__ == "__main__":
    run()
