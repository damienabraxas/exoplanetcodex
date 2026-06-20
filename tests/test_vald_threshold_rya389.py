"""
tests/test_vald_threshold_rya389.py
===================================
RYA-389 item 3 — the VALD extraction-threshold consistency check: every delivery's
detection depth must match the synthesis-era canonical 0.001 (RYA-285/387); the
EW-era 0.05 is under-deep and drops blends + trace species (RYA-381).

  * `vald_parse.effective_extraction_threshold` / `verify_extraction_threshold` —
    the reusable intake-gate check.
  * `check_stewardship.check_vald_threshold` — the CI invariant (tracked per star).
"""
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / 'data' / 'linelists'))
import vald_parse as vp                                        # noqa: E402
import scripts.check_stewardship as sc                         # noqa: E402


def _vald_file(tmp_path, name, central_depths):
    """Minimal VALD long-format delivery: header + one data line per central depth."""
    lines = [' 3780.00000, 6910.00000, 10, 100, 1.0 Wavelength region, lines selected, '
             'lines processed, Vmicro\n']
    for i, cd in enumerate(central_depths):
        wl = 5000.0 + i
        lines.append(f"'Fe 1', {wl:.4f}, -1.000, 1.0000, 1.0, 2.0000, 1.0, "
                     f"1.000, 1.000, 1.000, 8.000,-6.000,-7.000, {cd:.5f},\n")
    p = tmp_path / name
    p.write_text(''.join(lines))
    return p


class TestEffectiveThreshold:
    def test_min_central_depth_is_the_threshold(self, tmp_path):
        p = _vald_file(tmp_path, 'vald_x_raw.txt', [0.5, 0.1, 0.001, 0.02])
        assert vp.effective_extraction_threshold(p) == pytest.approx(0.001)

    def test_ew_era_delivery(self, tmp_path):
        p = _vald_file(tmp_path, 'vald_x_raw.txt', [0.9, 0.2, 0.05])
        assert vp.effective_extraction_threshold(p) == pytest.approx(0.05)


class TestVerifyThreshold:
    def test_accepts_synthesis_grade(self, tmp_path):
        p = _vald_file(tmp_path, 'vald_solar_raw.txt', [0.5, 0.001])
        verdict, msg, eff = vp.verify_extraction_threshold(p)
        assert verdict == 'ACCEPT' and eff == pytest.approx(0.001)

    def test_flags_under_deep_0p05(self, tmp_path):
        p = _vald_file(tmp_path, 'vald_procyon_raw.txt', [0.5, 0.05])
        verdict, msg, eff = vp.verify_extraction_threshold(p)
        assert verdict == 'FLAG' and 'under-deep' in msg and eff == pytest.approx(0.05)

    def test_rejects_unparseable(self, tmp_path):
        p = tmp_path / 'vald_empty_raw.txt'
        p.write_text('no data lines here\n')
        verdict, msg, eff = vp.verify_extraction_threshold(p)
        assert verdict == 'REJECT' and eff is None

    def test_canonical_constant(self):
        assert vp.THRESHOLD_CANONICAL == 0.001


class TestStewardshipInvariant:
    def test_current_deliveries_all_tracked(self):
        # The committed under-deep deliveries (αCen/Procyon/55 Cnc at 0.05) are tracked
        # against their re-extraction tickets → CI green, never UNTRACKED.
        v = sc.check_vald_threshold()
        assert all(x.ticket is not None for x in v)
        assert {x.ticket for x in v} <= {'RYA-382', 'RYA-384', 'RYA-385', 'RYA-387'}

    def test_optical_cores_pass(self):
        # solar + 55 Cnc optical are 0.001 → not flagged.
        flagged = {x.locus for x in sc.check_vald_threshold()}
        assert 'vald_solar_raw.txt' not in flagged
        assert 'vald_55cnc_raw.txt' not in flagged

    def test_fires_untracked_on_unknown_star(self, tmp_path, monkeypatch):
        d = tmp_path / 'data' / 'linelists'
        d.mkdir(parents=True)
        _vald_file(d, 'vald_newstar_raw.txt', [0.5, 0.05])     # under-deep, no ticket
        monkeypatch.setattr(sc, '_REPO', tmp_path)
        v = sc.check_vald_threshold()
        assert len(v) == 1 and v[0].ticket is None             # UNTRACKED → loud FAIL
        assert v[0].invariant == 'vald_threshold'

    def test_quarantine_files_skipped(self, tmp_path, monkeypatch):
        d = tmp_path / 'data' / 'linelists'
        d.mkdir(parents=True)
        _vald_file(d, 'vald_x_hfsoff_quarantine.txt', [0.5, 0.05])   # under-deep but quarantined
        monkeypatch.setattr(sc, '_REPO', tmp_path)
        assert sc.check_vald_threshold() == []                 # out-of-band, not flagged


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-v']))
