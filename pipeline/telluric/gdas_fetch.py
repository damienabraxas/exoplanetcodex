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
  3. NOAA ARL archive — AUTOMATED (RYA-983). GDAS1 by site lat-lon + datetime, pulled
                      straight from the NOAA Air Resources Lab archive at
                      https://www.ready.noaa.gov/data/archives/gdas1/ .
                      🔴 THIS DOCSTRING USED TO SAY THERE WAS NO ENDPOINT TO CALL. That
                      is true of the INTERACTIVE READY form (READYamet.php) and false of
                      the archive, which is a plain file server: 1286 weekly files over
                      ordinary HTTPS. Believing the note cost RYA-931 and RYA-973 a
                      manual pull per night, and left every La Silla instrument — HARPS,
                      FEROS — telluric-blocked on a wall that was not there.
                      💡 AND WE DO NOT DOWNLOAD THE WEEKLY FILE. ARL is a fixed-record
                      format, so one 3-hour slot is a computable byte range: 164 records
                      x 65210 B = 10.7 MB out of 599 MB, a 56x saving (tau Ceti's 29
                      nights: ~310 MB instead of ~17 GB). The server sends
                      `Accept-Ranges: bytes`, and the slice's own ARL header is checked
                      against the slot we asked for before it is used.
  4. NOAA ARL manual — the operator-staged ASCII, kept as the last resort for a night
                      the archive cannot serve (same 4-column P/HGT/T/RELHUM format).

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
from datetime import date, datetime, timedelta
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

# ── NOAA ARL GDAS1 archive (RYA-983) ──────────────────────────────────────────
_ARL_BASE = "https://www.ready.noaa.gov/data/archives/gdas1"
#: 🔴 A HARD FLOOR, NOT A BACKLOG ITEM. The archive begins here; a night before it can
#: NEVER get a real per-night profile from this source, and under the RYA-380 no-fallback
#: rule those frames are permanently un-correctable by this route. Three of tau Ceti's 32
#: HARPS nights (2004-09-19, 2004-10-02, 2004-10-03) sit below it.
ARL_ARCHIVE_START = date(2004, 12, 1)
#: ARL is fixed-record: 50-byte header + a 360x181 byte grid, 164 records per 3-hour slot.
#: Verified against the live archive — 598888640 B / 65210 = 9184 records = exactly 56
#: slots = 7 days x 8. These mirror pipeline.telluric.arl_gdas and are asserted against it.
_ARL_RECORD_LEN = 50 + 360 * 181
_ARL_RECORDS_PER_SLOT = 164
_ARL_SLOT_BYTES = _ARL_RECORD_LEN * _ARL_RECORDS_PER_SLOT


def arl_week_file(when: datetime) -> str:
    """The weekly archive file holding `when`. ARL weeks are calendar-anchored: w1 =
    days 1-7, w2 = 8-14, w3 = 15-21, w4 = 22-28, w5 = 29-end (so w5 is short)."""
    return f"gdas1.{when.strftime('%b').lower()}{when.strftime('%y')}.w{(when.day - 1) // 7 + 1}"


def arl_slot_offset(when: datetime) -> int:
    """Byte offset of `when`'s 3-hour slot inside its weekly file. Slots run sequentially
    from the first day of the week at 00 UT."""
    week_start_day = ((when.day - 1) // 7) * 7 + 1
    slot = (when.day - week_start_day) * 8 + when.hour // _GDAS_3H
    return slot * _ARL_SLOT_BYTES


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
            or _try_noaa_arl_archive(gdas_loc, lat, lon, slot, cache_dir, cached)
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


# ── NOAA ARL archive backend (RYA-983, automated) ──────────────────────────────

class ArlBeforeArchive(GDASUnavailable):
    """The night predates the ARL archive (2004-12-01). This is a PERMANENT data gap,
    not a retryable miss, and it is a distinct exception so a caller can tell "nobody has
    this and nobody ever will" from "the fetch failed, try again"."""


def _arl_fetch_slot(slot: datetime, work_dir: Path, timeout: float = 180.0) -> Path:
    """Download ONLY the 3-hour slot's byte range from the weekly ARL file, and verify
    the slice is what we asked for before returning it.

    The verification is the point. The record layout is a property of the ARL era, not a
    law, so the computed offset is a hypothesis: the first record of the slice carries
    its own timestamp, and if that disagrees with the slot we wanted we must NOT decode
    it — a silently-wrong atmosphere is exactly the failure the whole no-fallback rule
    exists to prevent (RYA-373)."""
    import urllib.error
    import urllib.request

    name = arl_week_file(slot)
    off = arl_slot_offset(slot)
    url = f"{_ARL_BASE}/{name}"
    req = urllib.request.Request(
        url, headers={'Range': f'bytes={off}-{off + _ARL_SLOT_BYTES - 1}',
                      'User-Agent': 'exoplanetcodex/RYA-983 (telluric GDAS retrieval)'})
    work_dir.mkdir(parents=True, exist_ok=True)
    out = work_dir / f"{name}.slot{off}"
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            if r.status not in (200, 206):
                raise GDASUnavailable(f"ARL {url} returned HTTP {r.status}")
            if r.status == 200:
                # No partial content: the server ignored Range and is sending 599 MB.
                # Refuse rather than quietly pulling the whole file 32 times.
                raise GDASUnavailable(
                    f"ARL {url} ignored the Range request (HTTP 200, not 206). Refusing "
                    f"to pull {_ARL_SLOT_BYTES / 1e6:.0f} MB as {598:.0f} MB.")
            data = r.read()
    except urllib.error.HTTPError as exc:
        raise GDASUnavailable(f"ARL {url}: HTTP {exc.code} {exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise GDASUnavailable(f"ARL {url}: {exc.reason}") from exc

    if len(data) != _ARL_SLOT_BYTES:
        raise GDASUnavailable(
            f"ARL {name}: asked for {_ARL_SLOT_BYTES} B at {off}, got {len(data)} B")
    if len(data) % _ARL_RECORD_LEN:
        raise GDASUnavailable(f"ARL {name}: slice is not record-aligned")

    head = data[:50].decode('ascii', 'replace')
    try:
        got = datetime(2000 + int(head[0:2]), int(head[2:4]), int(head[4:6]),
                       int(head[6:8]))
    except ValueError as exc:
        raise GDASUnavailable(f"ARL {name}: slice header is not an ARL header "
                              f"({head[:14]!r})") from exc
    if got != slot:
        raise GDASUnavailable(
            f"ARL {name}: the slice at offset {off} carries {got:%Y-%m-%dT%H} but "
            f"{slot:%Y-%m-%dT%H} was requested. The record layout differs from the "
            f"assumed {_ARL_RECORDS_PER_SLOT} records/slot for this era — refusing to "
            f"decode an atmosphere from the wrong slot.")
    out.write_bytes(data)
    return out


def _try_noaa_arl_archive(gdas_loc: str, lat: float, lon: float, slot: datetime,
                          cache_dir: Path, out_fits: Path) -> "Path | None":
    """Automated NOAA ARL pull → molecfit FITS. Returns None only when the caller should
    fall through to the manual-pull branch; a PERMANENT gap raises instead."""
    if slot.date() < ARL_ARCHIVE_START:
        raise ArlBeforeArchive(
            f"{slot:%Y-%m-%d} predates the NOAA ARL GDAS1 archive "
            f"({ARL_ARCHIVE_START}). No real per-night profile exists for this night "
            f"from this source and none ever will; under the RYA-380 no-fallback rule "
            f"these frames are permanently un-correctable by this route.")
    if os.environ.get('CODEX_NO_NETWORK'):
        return None
    from pipeline.telluric.arl_gdas import extract_profile, write_molecfit_ascii
    try:
        raw = _arl_fetch_slot(slot, cache_dir / '_arl')
    except GDASUnavailable:
        return None                     # let the manual-pull branch try
    try:
        profile = extract_profile(raw, slot, lat, lon)
        ascii_path = write_molecfit_ascii(
            profile, cache_dir / f"noaa_{gdas_loc}_{slot:%Y-%m-%dT%H}.txt")
    finally:
        # The slice is 10.7 MB of a 599 MB file we deliberately did not keep; the
        # ASCII and the FITS are the artifacts worth caching.
        raw.unlink(missing_ok=True)
    return _gdas_ascii_to_fits(ascii_path, out_fits)


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


# ── CLI (RYA-963 STEP 0: the gate is reported before any molecfit runs) ───────
def main(argv=None) -> int:
    """Report whether a REAL per-night GDAS profile resolves for a site/date, and for
    which exposure times. The whole point of running this first is that the one external
    dependency of a telluric leg is data availability, not compute — findable up front,
    and a clean STOP if it is absent."""
    import argparse
    ap = argparse.ArgumentParser(description=main.__doc__)
    ap.add_argument('--site', default='paranal')
    ap.add_argument('--date', required=True, help='observation night, YYYY-MM-DD')
    ap.add_argument('--times', default=None,
                    help='comma-separated UTC HH:MM of the exposures')
    ap.add_argument('--verify', action='store_true',
                    help='also open the profile and report its physical ranges')
    a = ap.parse_args(argv)

    times = [t.strip() for t in a.times.split(',')] if a.times else ['00:00']
    slots = {}
    for t in times:
        hh, mm = (int(x) for x in t.split(':'))
        when = datetime.strptime(a.date, '%Y-%m-%d').replace(hour=hh, minute=mm)
        slot = nearest_3hourly(when)
        slots.setdefault(f"{slot:%Y-%m-%dT%H}", []).append(t)

    print(f"GDAS gate: site={a.site} night={a.date} exposures={','.join(times)}")
    rc = 0
    for slot_key, ts in sorted(slots.items()):
        hh, mm = (int(x) for x in ts[0].split(':'))
        when = datetime.strptime(a.date, '%Y-%m-%d').replace(hour=hh, minute=mm)
        try:
            path = fetch_gdas(a.site, night=when)
        except GDASUnavailable as exc:
            print(f"  slot {slot_key}  STOP — {exc}")
            rc = 1
            continue
        print(f"  slot {slot_key}  {Path(path).name}   (exposures {', '.join(ts)})")
        if a.verify:
            from astropy.io import fits
            d = fits.getdata(str(path))
            print(f"      {len(d)} levels  press {d['press'].min():g}-{d['press'].max():g} hPa"
                  f"  height {d['height'].min():.3f}-{d['height'].max():.3f} km"
                  f"  temp {d['temp'].min():g}-{d['temp'].max():g} K"
                  f"  relhum {d['relhum'].min():g}-{d['relhum'].max():g} %")
    print(f"  {len(times)} exposure(s) -> {len(slots)} distinct 3-hourly profile(s); "
          f"standard-atmosphere fallback: forbidden (RYA-373)")
    return rc


if __name__ == '__main__':
    raise SystemExit(main())
