"""RYA-983 — the automated NOAA ARL GDAS1 pull.

These run OFFLINE. The arithmetic and the floor are the parts that must never drift, and
both are pure; the one thing that needs the network (a real fetch) is exercised by hand
and recorded in the ticket, not run in CI on a runner with no route out.
"""
from __future__ import annotations

import os
import sys
from datetime import date, datetime
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.telluric import arl_gdas                              # noqa: E402
from pipeline.telluric.gdas_fetch import (                          # noqa: E402
    ARL_ARCHIVE_START, ArlBeforeArchive, GDASUnavailable, _ARL_RECORD_LEN,
    _ARL_RECORDS_PER_SLOT, _ARL_SLOT_BYTES, _try_noaa_arl_archive, arl_slot_offset,
    arl_week_file)


# ── the record layout, which the byte-range trick depends on ────────────────
def test_record_layout_matches_the_decoder_it_feeds():
    """gdas_fetch computes byte offsets; arl_gdas decodes what lands there. If the two
    ever disagree about the record size the offsets point into the middle of a record and
    the header check fires — but they should not be able to disagree in the first place."""
    assert _ARL_RECORD_LEN == 50 + arl_gdas.NX * arl_gdas.NY
    assert _ARL_RECORDS_PER_SLOT == arl_gdas.RECORDS_PER_SLOT
    assert _ARL_SLOT_BYTES == _ARL_RECORD_LEN * _ARL_RECORDS_PER_SLOT


def test_a_weekly_file_is_exactly_56_slots():
    """598888640 B / 65210 = 9184 records = 56 slots = 7 days x 8. Verified against the
    live archive; pinned because the whole offset scheme rests on it."""
    assert 598888640 % _ARL_RECORD_LEN == 0
    assert 598888640 // _ARL_SLOT_BYTES == 56
    assert 56 == 7 * 8


# ── week/offset arithmetic, checked against live-verified values ────────────
@pytest.mark.parametrize("when,name", [
    (datetime(2023, 12, 1, 0), 'gdas1.dec23.w1'),
    (datetime(2023, 12, 7, 21), 'gdas1.dec23.w1'),
    (datetime(2023, 12, 8, 0), 'gdas1.dec23.w2'),
    (datetime(2023, 12, 29, 0), 'gdas1.dec23.w5'),
    (datetime(2004, 12, 1, 0), 'gdas1.dec04.w1'),
])
def test_week_file_naming(when, name):
    assert arl_week_file(when) == name


@pytest.mark.parametrize("when,offset", [
    (datetime(2023, 12, 1, 0), 0),
    (datetime(2023, 12, 3, 6), 192_499_920),      # probed live: header '2312 3 6'
    (datetime(2023, 12, 7, 21), 588_194_200),     # probed live: header '2312 721 3'
])
def test_slot_offsets_match_the_live_archive(when, offset):
    assert arl_slot_offset(when) == offset


def test_every_slot_of_a_week_lands_inside_the_file():
    for day in range(1, 8):
        for hour in range(0, 24, 3):
            off = arl_slot_offset(datetime(2023, 12, day, hour))
            assert 0 <= off and off + _ARL_SLOT_BYTES <= 598888640


# ── the permanent floor ─────────────────────────────────────────────────────
def test_before_the_archive_is_a_PERMANENT_gap_not_a_retryable_miss(tmp_path):
    """A distinct exception, so a caller can tell 'nobody has this and nobody ever will'
    from 'the fetch failed, try again'. Three of tau Ceti's 32 HARPS nights are here."""
    assert ARL_ARCHIVE_START == date(2004, 12, 1)
    assert issubclass(ArlBeforeArchive, GDASUnavailable)
    for night in (datetime(2004, 9, 19, 0), datetime(2004, 11, 30, 21)):
        with pytest.raises(ArlBeforeArchive, match='permanently un-correctable'):
            _try_noaa_arl_archive('C-70.7-29.3', -29.26, -70.73, night,
                                  tmp_path, tmp_path / 'x.fits', site_registered=True)


def test_the_first_night_of_the_archive_is_not_refused(tmp_path, monkeypatch):
    """The boundary must be inclusive; verified live that 2004-12-01 fetches."""
    monkeypatch.setenv('CODEX_NO_NETWORK', '1')
    assert _try_noaa_arl_archive('C-70.7-29.3', -29.26, -70.73,
                                 datetime(2004, 12, 1, 0), tmp_path,
                                 tmp_path / 'x.fits', site_registered=True) is None


def test_no_network_short_circuits_without_pretending_to_succeed(tmp_path, monkeypatch):
    monkeypatch.setenv('CODEX_NO_NETWORK', '1')
    assert _try_noaa_arl_archive('C-70.7-29.3', -29.26, -70.73,
                                 datetime(2023, 12, 3, 6), tmp_path,
                                 tmp_path / 'x.fits', site_registered=True) is None


# ── the docstring that cost two tickets a manual pull each ──────────────────
def test_the_docstring_no_longer_says_there_is_no_endpoint():
    """It used to state there was 'no clean unauthenticated REST endpoint to call from
    code'. True of the interactive READY form, false of the archive — and believing it
    left every La Silla instrument telluric-blocked on a wall that was not there."""
    from pipeline.telluric import gdas_fetch
    doc = gdas_fetch.__doc__
    assert 'no clean unauthenticated REST endpoint' not in doc
    assert 'ready.noaa.gov/data/archives/gdas1' in doc


# ── the safety property RYA-380's test caught ───────────────────────────────
def test_the_arl_backend_refuses_an_UNREGISTERED_site(tmp_path):
    """🔴 The ESO tarball backend is keyed by SITE, so an unknown site gets nothing and
    RYA-380's 'unknown site must fail loud' guarantee holds. The ARL archive is keyed by
    COORDINATE and is GLOBAL — it returns a real, correct profile for any lat/lon on
    Earth, including lat 0 / lon 0 for a site that does not exist, which is the middle of
    the Atlantic. The number is not fabricated; it is the atmosphere above the WRONG
    PLACE, which is worse, because nothing about it looks wrong.

    RYA-380's own test caught this the first time the backend ran. Registered sites only.
    """
    assert _try_noaa_arl_archive('C-999.9-99.9', 0.0, 0.0, datetime(2023, 12, 3, 6),
                                 tmp_path, tmp_path / 'x.fits',
                                 site_registered=False) is None


def test_an_unregistered_site_still_raises_end_to_end(tmp_path):
    """The whole point: fetch_gdas must still refuse a site nobody has, even though the
    ARL archive could physically answer for its coordinates."""
    from pipeline.telluric.gdas_fetch import fetch_gdas
    with pytest.raises(GDASUnavailable):
        fetch_gdas('nowhere', night='2023-12-03', cache_dir=tmp_path,
                   lat=0.0, lon=0.0, gdas_loc='C-999.9-99.9')


def test_registered_sites_are_the_ones_that_may_use_it():
    """la_silla must be reachable — unblocking it is the whole ticket — and the check is
    the registry, not a name allow-list here."""
    from config.constants import SITES, get_site
    assert 'la_silla' in SITES and 'paranal' in SITES
    assert get_site('la_silla')['gdas_loc'] == 'C-70.7-29.3'
    with pytest.raises(Exception):
        get_site('nowhere')
