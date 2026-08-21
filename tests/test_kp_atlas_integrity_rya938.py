"""RYA-938 — a corrupt atlas segment must not present as missing coverage.

`lm0840` in both staged copies of the 1984 Kitt Peak atlas was a saved HTTP 500
page. `measure_band_ew.kp_segments()` swallowed the parse failure in a bare
`except Exception: continue`, so the file vanished from the inventory and
8400-8440 A answered "no Kitt Peak segment covers 8420.000 A" -- a coverage
message for a data fault, which is the RYA-833 shape.

These tests pin what a segment IS and what the loader must do about one that
isn't, so they keep holding when the example file is repaired (it now is).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pipeline.kp_atlas_integrity import (
    KpAtlasCorrupt, inspect_atlas, inspect_segment, require_parseable)

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "data" / "audit" / "rya938_kp_crosscheck"

HTML_500 = (
    '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">\n'
    "<html><head>\n<title>500 Internal Server Error</title>\n</head><body>\n"
    "<h1>Internal Server Error</h1>\n</body></html>\n")


def _segment(path: Path, start_nm: float, n: int = 400, *, flux: float = 0.98) -> Path:
    wave = start_nm + np.arange(n) * 0.001
    data = np.column_stack([wave, np.full(n, flux), np.full(n, 100.0)])
    np.savetxt(path, data, fmt="%12.5f %12.7f %10.2f")
    return path


def test_a_saved_error_page_is_not_data_and_says_why(tmp_path):
    """The reason must name the real fault, not merely 'unparseable'."""
    bad = tmp_path / "lm0840"
    bad.write_text(HTML_500)
    report = inspect_segment(bad)
    assert report.ok is False
    assert "html" in report.reason.lower()
    assert "download" in report.reason.lower()


def test_a_good_segment_passes_and_reports_its_real_span(tmp_path):
    """Positive control: the test must be able to say yes, or it proves nothing."""
    report = inspect_segment(_segment(tmp_path / "lm0840", 840.0))
    assert report.ok is True, report.reason
    assert report.monotonic is True
    assert report.n_columns == 3
    assert report.lo_A == pytest.approx(8400.0, abs=0.1)
    assert report.stem_matches_data is True


def test_non_monotonic_wavelength_is_refused(tmp_path):
    path = tmp_path / "lm0840"
    wave = np.r_[840.0 + np.arange(200) * 0.001, 840.0 + np.arange(200) * 0.001]
    np.savetxt(path, np.column_stack([wave, np.full(400, 0.98), np.full(400, 100.0)]))
    report = inspect_segment(path)
    assert report.ok is False
    assert "increasing" in report.reason


def test_filename_that_disagrees_with_its_data_is_refused(tmp_path):
    """The lm#### stem is a claim about content; a mismatch is a mislabelled file."""
    report = inspect_segment(_segment(tmp_path / "lm0840", 900.0))
    assert report.ok is False
    assert "filename" in report.reason


def test_absurd_flux_is_refused(tmp_path):
    """Residual flux lives near unity; a raw-irradiance column here is a wrong file."""
    report = inspect_segment(_segment(tmp_path / "lm0840", 840.0, flux=1.0e4))
    assert report.ok is False
    assert "flux" in report.reason


def test_require_parseable_names_every_bad_file(tmp_path):
    good = inspect_segment(_segment(tmp_path / "lm0836", 836.0))
    (tmp_path / "lm0840").write_text(HTML_500)
    bad = inspect_segment(tmp_path / "lm0840")
    require_parseable([good])                       # silent when all good
    with pytest.raises(KpAtlasCorrupt) as excinfo:
        require_parseable([good, bad])
    message = str(excinfo.value)
    assert "lm0840" in message
    assert "NOT a coverage gap" in message, "the message must not read as missing data"


def test_inspect_atlas_separates_corruption_from_a_coverage_gap(tmp_path):
    """The two failure modes are different questions and must report differently."""
    _segment(tmp_path / "lm0836", 836.0, n=4000)     # 8360-8400 A
    _segment(tmp_path / "lm0844", 844.0, n=4000)     # 8440-8480 A
    report = inspect_atlas(tmp_path)
    assert report["n_corrupt"] == 0
    assert report["n_gaps"] == 1, "a real hole between segments is a GAP"
    (tmp_path / "lm0840").write_text(HTML_500)
    report = inspect_atlas(tmp_path)
    assert report["n_corrupt"] == 1, "an unreadable file is CORRUPTION, not a gap"


def test_the_staged_atlas_is_intact_and_seamless():
    """Recorded state of the repaired atlas: every segment parses, no holes."""
    report = json.loads((EVIDENCE / "atlas_integrity.json").read_text())
    assert report["n_corrupt"] == 0
    assert report["n_gaps"] == 0
    assert report["n_files"] == report["n_ok"]
    summary = report["segments_summary"]
    assert summary["all_monotonic"] and summary["all_stem_consistent"]
