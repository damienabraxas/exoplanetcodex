"""RYA-1074 — verification must resolve the run's ACTUAL output, never rebuild it.

THE DEFECT, LIVE ON 2026-08-27. The coverage guard narrowed a requested 4200-6910 A HARPS
band to the arm's real extent and wrote every product under `FeI_4200_6908_harps_...`. A
verification step globbed the REQUESTED stem `FeI_4200_6910_harps_*`, matched files a
PREVIOUS run had left there, and reported "EXACT reproduction". It could not have failed:
it compared the store against files that run never wrote.

The trim is correct science and is NOT touched here (a window half outside the data is not
a measurement). What is fixed is that the rename was SILENT, so any check keyed on the
requested band was vacuous.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.run_output_resolver import (  # noqa: E402
    StemResolutionError, requested_stem, resolve_products, resolve_stem,
)

REQ = dict(element="Fe", ion="I", lo=4200.0, hi=6910.0, instrument="harps")
REQ_STEM = "FeI_4200_6910_harps"
EFF_STEM = "FeI_4200_6908_harps"


def _runinfo(d: Path, stem: str, req_stem: str, files: list[str], when="2026-08-27T12:00:00Z"):
    (d / f"{stem}_runinfo.json").write_text(json.dumps({
        "ticket": "RYA-1074", "written_utc": when,
        "element": "Fe", "ion": "I", "instrument": "harps", "holding": None,
        "requested_band_A": [4200.0, 6910.0], "effective_band_A": [4200.0, 6908.98],
        "band_trimmed": True, "stem": stem, "requested_stem": req_stem,
        "files_written": files, "base_sha": "deadbeef",
    }) + "\n")
    for f in files:
        (d / f).write_text("element,ion,A\nFe,I,7.512\n")


def test_the_requested_stem_is_not_what_gets_written():
    """The premise: int() on a trimmed edge changes the filename."""
    assert requested_stem(**REQ) == REQ_STEM
    assert requested_stem(element="Fe", ion="I", lo=4200.0, hi=6908.98,
                          instrument="harps") == EFF_STEM


def test_a_trimmed_run_resolves_to_its_effective_stem(tmp_path):
    _runinfo(tmp_path, EFF_STEM, REQ_STEM, [f"{EFF_STEM}_products.csv"])
    assert resolve_stem(tmp_path, **REQ) == EFF_STEM


def test_POSITIVE_CONTROL_a_stale_requested_stem_file_does_not_satisfy_the_check(tmp_path):
    """🔴 THE TICKET'S REQUIRED CONTROL, AND THE EXACT SHAPE OF THE ORIGINAL BUG.

    Stage a STALE file at the requested stem — as a previous run really had left one —
    then record a run that wrote the TRIMMED stem. The resolver must return the trimmed
    stem and must never offer the stale file, because matching it is how a check passes
    while comparing nothing the run produced.
    """
    stale = tmp_path / f"{REQ_STEM}_products.csv"
    stale.write_text("element,ion,A\nFe,I,9.999\n")      # a wrong value, deliberately
    _runinfo(tmp_path, EFF_STEM, REQ_STEM, [f"{EFF_STEM}_products.csv"])

    assert resolve_stem(tmp_path, **REQ) == EFF_STEM, "resolved to the stale stem"
    got = resolve_products(tmp_path, **REQ)
    for p in got.values():
        assert p.name.startswith(EFF_STEM), f"comparison reached a stale file: {p.name}"
        assert "9.999" not in p.read_text(), "the stale value entered the comparison"
    assert stale.exists(), "the control is only meaningful while the stale file is present"


def test_no_runinfo_means_REFUSE_never_fall_back_to_the_requested_stem(tmp_path):
    """🔴 NO FALLBACK, BY DESIGN. Falling back is the bug: a stale file at the requested
    path makes the fallback SUCCEED, silently, against the wrong data."""
    (tmp_path / f"{REQ_STEM}_products.csv").write_text("element,ion,A\nFe,I,9.999\n")
    with pytest.raises(StemResolutionError) as e:
        resolve_stem(tmp_path, **REQ)
    assert "runinfo" in str(e.value).lower()


def test_a_run_that_asked_for_something_else_does_not_answer_this_request(tmp_path):
    _runinfo(tmp_path, "FeI_3000_3780_kpno", "FeI_3000_3780_kpno",
             ["FeI_3000_3780_kpno_products.csv"])
    with pytest.raises(StemResolutionError) as e:
        resolve_stem(tmp_path, **REQ)
    assert "no run recorded writing" in str(e.value)


def test_an_unreadable_manifest_raises_rather_than_being_skipped(tmp_path):
    """Skipping a corrupt manifest would drop the caller back to globbing — the defect."""
    _runinfo(tmp_path, EFF_STEM, REQ_STEM, [f"{EFF_STEM}_products.csv"])
    (tmp_path / "broken_runinfo.json").write_text("{not json")
    with pytest.raises(StemResolutionError) as e:
        resolve_stem(tmp_path, **REQ)
    assert "unreadable" in str(e.value)


def test_only_files_the_run_RECORDED_enter_the_comparison(tmp_path):
    """A file that appeared afterwards by any other route is not this run's output."""
    _runinfo(tmp_path, EFF_STEM, REQ_STEM, [f"{EFF_STEM}_products.csv"])
    (tmp_path / f"{EFF_STEM}_ENGINE-Z_products.csv").write_text("element,ion,A\nFe,I,0.0\n")
    got = resolve_products(tmp_path, **REQ)
    assert "ENGINE-Z" not in got, "an unrecorded file entered the comparison"


def test_the_newest_run_wins_when_a_request_was_served_twice(tmp_path):
    _runinfo(tmp_path, "FeI_4200_6907_harps", REQ_STEM,
             ["FeI_4200_6907_harps_products.csv"], when="2026-08-26T09:00:00Z")
    _runinfo(tmp_path, EFF_STEM, REQ_STEM, [f"{EFF_STEM}_products.csv"],
             when="2026-08-27T12:00:00Z")
    assert resolve_stem(tmp_path, **REQ) == EFF_STEM


def test_the_deriver_still_trims_and_still_announces_it():
    """⚠️ DO NOT 'FIX' THIS BY REMOVING THE TRIM. The narrowing is correct science; only
    its silence was the defect. Pin both the trim and the new announcement."""
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    assert "a.lo, a.hi = float(_lo_new), float(_hi_new)" in src, "the trim was removed"
    assert "a.requested_lo, a.requested_hi" in src, "the request is not captured"
    assert "BAND TRIMMED" in src, "the trim no longer announces the effective stem"
    assert "write_runinfo(out, stem, a)" in src, "the run no longer records what it wrote"


def test_files_written_excludes_files_that_predate_the_run(tmp_path, monkeypatch):
    """🔴 THE HOLE THE FIRST FIX LEFT OPEN, caught live on 2026-08-27.

    `write_runinfo` globbed `{stem}_*` at write time. That is safe when the coverage guard
    RENAMED the stem — nothing else can be sitting there — and wrong whenever it did not:
    a 14:37 run recorded four 07:34 files from a previous run as its own output. Stale-file
    absorption, one level down from the defect this module exists to prevent, and inside
    the very manifest verification is told to trust.

    `files_written` must now carry only files whose mtime is at or after the run's start;
    anything older is reported under `pre_existing_at_stem` rather than dropped silently.
    """
    import os
    import time as _t
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_dbp", ROOT / "scripts" / "derive_band_products.py")

    stem = "FeI_4200_6910_kpno"
    old = tmp_path / f"{stem}_ENGINE-B-NLTE_products.csv"
    old.write_text("element,ion,A\nFe,I,7.497\n")
    os.utime(old, (_t.time() - 7200, _t.time() - 7200))       # two hours stale

    t0 = _t.time()
    fresh = tmp_path / f"{stem}_products.csv"
    fresh.write_text("element,ion,A\nFe,I,7.447\n")

    class A:                                                   # minimal stand-in
        element, ion, instrument, holding = "Fe", "I", "kpno", None
        lo, hi = 4200.0, 6910.0
        requested_lo, requested_hi = 4200.0, 6910.0
        band_trimmed = False
        _run_t0 = t0

    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    assert "_run_t0" in src and "pre_existing_at_stem" in src, "the mtime filter is gone"

    import json as _j
    ns: dict = {}
    exec(compile(src[src.index("def write_runinfo"):src.index("def _base_sha")],
                 "<write_runinfo>", "exec"),
         {"Path": Path, "json": _j, "time": _t, "print": lambda *a, **k: None,
          "_base_sha": lambda: "test"}, ns)
    ns["_base_sha"] = lambda: "test"
    ns["write_runinfo"](tmp_path, stem, A())

    info = _j.loads((tmp_path / f"{stem}_runinfo.json").read_text())
    assert f"{stem}_products.csv" in info["files_written"]
    assert f"{stem}_ENGINE-B-NLTE_products.csv" not in info["files_written"], \
        "a file predating the run was recorded as this run's output"
    assert f"{stem}_ENGINE-B-NLTE_products.csv" in info["pre_existing_at_stem"], \
        "the pre-existing file was dropped silently instead of being named"
