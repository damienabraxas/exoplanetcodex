"""
pipeline/telluric/gdas_fetch.py
===============================
RYA-380 — per-night GDAS atmospheric-profile retrieval for molecfit telluric
correction. The single reusable retriever behind the standing wavelength-gated
telluric recipe: every red-optical (λ ≳ 6800 Å) + IR dataset gets the real
observation-night GDAS profile (the night's T/P/humidity over the site) before
molecfit runs.

LOUD-FAIL is the whole point. molecfit's `GDAS_PROFILE=auto` silently falls back to
a generic standard atmosphere when it can't locate the right hourly profile — that is
the RYA-373 CRITICAL bug (the conditioned Vesta CO came out telluric-dominated). This
module returns a REAL per-night profile or raises `GDASUnavailable`. There is
deliberately NO standard-atmosphere fallback path here; a silent fallback can never
originate from this code.

Retrieval order (the standing recipe):
  1. cache          — a previously fetched/staged per-night profile (FITS)
  2. ESO GDAS       — the per-site 3-hourly tarball shipped with ESO telluriccorr,
                      `…/share/molecfit/data/profiles/gdas/gdas_profiles_<loc>.tar.gz`.
                      Mechanic (RYA-373 commit 8b0b551): molecfit auto-mode requests
                      odd hours (T01/T02) absent from the 3-hourly tarball → fallback;
                      we instead extract the NEAREST real 3-hourly profile for the obs
                      MJD and convert the ASCII profile to the FITS molecfit's CFITSIO
                      `GDAS_PROF` loader requires.
  3. NOAA ARL READY — manual-pull fallback: GDAS1 by site lat-lon + datetime from the
                      NOAA Air Resources Lab READY archive
                      (https://www.ready.noaa.gov/READYamet.php ; archive
                      ftp://arlftp.arlhq.noaa.gov/pub/archives/gdas1/). There is no
                      clean unauthenticated REST endpoint to call from code, so the
                      operator stages the extracted ASCII (same 4-column P/HGT/T/RELHUM
                      format) into the cache and we ingest it. We do NOT fabricate an
                      HTTP endpoint.

Sources / formats (cited):
  • ESO GDAS tarball: telluriccorr 4.3.3 (`/opt/homebrew/Cellar/telluriccorr/4.3.3_4/
    share/molecfit/data/profiles/gdas/gdas_profiles_C-70.4-24.6.tar.gz`, Paranal),
    3-hourly ASCII entries `C<lon><lat>D<YYYY-MM-DD>T<HH>.gdas` with a header line then
    columns `P[hPa]  HGT[m]  T[K]  RELHUM[%]`. Discovered by glob (version-agnostic)
    with `MOLECFIT_GDAS_DIR` override.
  • molecfit GDAS_PROF FITS columns: `press[hPa], height[km], temp[K], relhum[%]`.

Site coordinates come from `config.constants.SITES` (single source, RYA-388) — never
hardcoded here. Generic over site / lat-lon / datetime so the red-optical arms and the
55 Cnc / α Cen CRIRES+ pulls reuse it unchanged.
"""
from __future__ import annotations

import glob
import os
import tarfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from config.constants import PATHS, get_site


class GDASUnavailable(RuntimeError):
    """No real per-night GDAS profile could be retrieved (cache + ESO + NOAA all
    empty). Raising this is correct behaviour — the caller MUST NOT fall back to a
    standard atmosphere (the RYA-373 failure mode). Either stage a manual NOAA ARL
    pull into the cache, or treat the dataset as GDAS-blocked."""


# ESO telluriccorr GDAS tarballs (per-site, 3-hourly). The candidate directories come
# from pipeline.telluric.esorex_runtime.gdas_dirs(), which puts the REGISTERED
# `eso_pipelines` prefix first. RYA-963: this module used to glob only Homebrew and
# /usr/share/esopipes, so on Sirius — where the ESO source kit lives under
# /srv/codex/eso/molecfit and ships the same Paranal tarball — the per-night profile was
# unreachable and GDASUnavailable would have fired on a profile sitting on the disk. That
# is the worst kind of loud failure: correct behaviour for a wrong reason. Override the
# whole directory with MOLECFIT_GDAS_DIR.
from pipeline.telluric.esorex_runtime import gdas_dirs as _gdas_dirs

_GDAS_3H = 3  # GDAS profiles are 3-hourly (00/03/06/…/21 UTC)


# ── public API ────────────────────────────────────────────────────────────────

def gdas_cache_path(cache_dir, gdas_loc: str, slot: datetime) -> Path:
    """Canonical cache filename for a resolved 3-hourly GDAS slot."""
    return Path(cache_dir) / f"gdas_{gdas_loc}_{slot:%Y-%m-%dT%H}.gdas.fits"


def nearest_3hourly(night, mjd: "float | None" = None) -> datetime:
    """Resolve the nearest real 3-hourly GDAS slot (UTC) for an observation. Prefer the
    exact obs MJD; else the `night` date (a `date`/`datetime`/`YYYY-MM-DD` str) at 00 UT.
    GDAS is 3-hourly, so the per-night profile is the slot nearest the actual exposure."""
    if mjd is not None:
        from astropy.time import Time
        t = Time(float(mjd), format='mjd').to_datetime()
    elif isinstance(night, datetime):
        t = night
    else:
        from datetime import date
        if isinstance(night, date):
            t = datetime(night.year, night.month, night.day)
        else:
            t = datetime.strptime(str(night)[:10], "%Y-%m-%d")
    hh = int(round((t.hour + t.minute / 60.0) / _GDAS_3H)) * _GDAS_3H
    base = t.replace(hour=0, minute=0, second=0, microsecond=0)
    return base + timedelta(hours=hh)


def fetch_gdas(site: str, night=None, *, mjd: "float | None" = None,
               cache_dir=None, work_dir=None,
               lat: "float | None" = None, lon: "float | None" = None,
               gdas_loc: "str | None" = None) -> Path:
    """Return a Path to a REAL per-night GDAS profile (FITS, molecfit-ready) for `site`
    on the observation night. Order: cache → ESO GDAS tarball → NOAA ARL manual pull →
    raise GDASUnavailable. NO standard-atmosphere fallback.

    `site` resolves coordinates + the ESO tarball id from config.constants.SITES; pass
    `lat`/`lon`/`gdas_loc` explicitly to override (generic use). One of `night`/`mjd`
    must be given (mjd preferred — exact exposure time)."""
    if night is None and mjd is None:
        raise ValueError("fetch_gdas needs one of night=/mjd= (the observation night).")
    # Generic path: fully-specified lat/lon/gdas_loc bypass the SITES registry (e.g. an
    # ad-hoc site). Otherwise resolve coordinates + tarball id from constants.SITES.
    if not (gdas_loc and lat is not None and lon is not None):
        rec = get_site(site)
        gdas_loc = gdas_loc or rec['gdas_loc']
        lat = rec['lat'] if lat is None else lat
        lon = rec['lon'] if lon is None else lon
    cache_dir = Path(cache_dir) if cache_dir else Path(PATHS['gdas_cache'])
    work_dir = Path(work_dir) if work_dir else cache_dir
    cache_dir.mkdir(parents=True, exist_ok=True)

    slot = nearest_3hourly(night, mjd)

    cached = gdas_cache_path(cache_dir, gdas_loc, slot)
    if cached.exists() and cached.stat().st_size > 0:
        return cached

    prof = (_try_eso_gdas(gdas_loc, slot, work_dir, cached)
            or _try_noaa_ready(gdas_loc, lat, lon, slot, cache_dir, cached))
    if prof is None:
        raise GDASUnavailable(
            f"No GDAS profile for site={site} ({gdas_loc}) near {slot:%Y-%m-%dT%H}:00 UT "
            f"(ESO tarball + NOAA ARL manual-pull both empty). A standard-atmosphere "
            f"fallback is forbidden (RYA-373 failure mode). Stage a NOAA ARL READY GDAS1 "
            f"pull for this night into {cache_dir} as "
            f"noaa_{gdas_loc}_{slot:%Y-%m-%dT%H}.txt (cols: P[hPa] HGT[m] T[K] RELHUM[%]).")
    return prof


# ── ESO GDAS tarball backend ───────────────────────────────────────────────────

def _eso_gdas_dir() -> "Path | None":
    """Locate the installed ESO telluriccorr GDAS profile directory."""
    override = os.environ.get('MOLECFIT_GDAS_DIR')
    if override and Path(override).is_dir():
        return Path(override)
    for pat in _gdas_dirs():
        hits = sorted(glob.glob(pat))
        if hits:
            return Path(hits[-1])          # newest install version
    return None


def _try_eso_gdas(gdas_loc: str, slot: datetime, work_dir: Path,
                  out_fits: Path) -> "Path | None":
    """Extract the nearest real 3-hourly profile for `slot` from the per-site ESO
    tarball and convert ASCII→FITS at `out_fits`. Returns None if the tarball/profile is
    absent (→ caller tries NOAA, then loud-fails)."""
    gdir = _eso_gdas_dir()
    if gdir is None:
        return None
    tarball = gdir / f"gdas_profiles_{gdas_loc}.tar.gz"
    if not tarball.exists():
        return None
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(tarball) as tf:
            for dh in (0, -3, 3, -6, 6):     # nearest real slot, then neighbours
                cand = slot + timedelta(hours=dh)
                name = f"{gdas_loc}D{cand:%Y-%m-%d}T{cand.hour:02d}.gdas"
                try:
                    member = tf.getmember(name)
                except KeyError:
                    continue
                tf.extract(member, path=str(work_dir))
                return _gdas_ascii_to_fits(work_dir / name, out_fits)
    except (OSError, tarfile.TarError):
        return None
    return None


# ── NOAA ARL READY manual-pull backend ─────────────────────────────────────────

def _try_noaa_ready(gdas_loc: str, lat: float, lon: float, slot: datetime,
                    cache_dir: Path, out_fits: Path) -> "Path | None":
    """Manual-pull fallback: ingest an operator-staged NOAA ARL READY GDAS1 profile.

    The READY archive (https://www.ready.noaa.gov/READYamet.php, GDAS1 by lat-lon +
    datetime; archive ftp://arlftp.arlhq.noaa.gov/pub/archives/gdas1/) is an interactive
    web request — there is no clean code endpoint to call, so we do NOT fabricate one.
    The operator extracts the GDAS1 column profile for this night to the 4-column ASCII
    `P[hPa] HGT[m] T[K] RELHUM[%]` and drops it in the cache as one of:
      noaa_<loc>_<YYYY-MM-DD>T<HH>.txt   |   noaa_<loc>_<YYYY-MM-DD>.txt
    We parse + convert it to the molecfit FITS. Returns None if no manual pull staged."""
    for cand in (cache_dir / f"noaa_{gdas_loc}_{slot:%Y-%m-%dT%H}.txt",
                 cache_dir / f"noaa_{gdas_loc}_{slot:%Y-%m-%d}.txt"):
        if cand.exists() and cand.stat().st_size > 0:
            return _gdas_ascii_to_fits(cand, out_fits)
    return None


# ── ASCII → molecfit FITS conversion ───────────────────────────────────────────

def _gdas_ascii_to_fits(ascii_path, out_fits) -> Path:
    """Convert a 4-column GDAS ASCII profile (`# P[hPa] HGT[m] T[K] RELHUM[%]`) to the
    FITS bintable molecfit's GDAS_PROF (CFITSIO) loader expects: columns
    `press[hPa], height[km], temp[K], relhum[%]` (metres → km)."""
    from astropy.io import fits
    ascii_path, out_fits = Path(ascii_path), Path(out_fits)
    rows = []
    for ln in ascii_path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        parts = ln.split()
        if len(parts) < 4:
            continue
        p, hgt_m, t, rh = (float(x) for x in parts[:4])
        rows.append((p, hgt_m / 1000.0, t, rh))      # m → km
    if not rows:
        raise ValueError(f"GDAS ASCII profile {ascii_path} parsed to zero rows "
                         f"(unexpected format).")
    arr = np.array(rows)
    out_fits.parent.mkdir(parents=True, exist_ok=True)
    fits.BinTableHDU.from_columns([
        fits.Column(name='press', format='1D', array=arr[:, 0]),
        fits.Column(name='height', format='1D', array=arr[:, 1]),
        fits.Column(name='temp', format='1D', array=arr[:, 2]),
        fits.Column(name='relhum', format='1D', array=arr[:, 3])],
    ).writeto(out_fits, overwrite=True)
    return out_fits
