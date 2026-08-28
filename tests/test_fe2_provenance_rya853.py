"""
tests/test_fe2_provenance_rya853.py — RYA-853, the Fe II half of scopes 1 and 4
===============================================================================
`canonical_gf` labelled Fe II 6149.246 / 6247.557 `NIST ASD v5.11 grade B`. NIST — queried
in AIR, EP-matched on both sides, uniquely — says -2.854 acc **E** and -2.444 acc **D**.
Wrong on both axes, and `B` sits in NIST_GRADE_HIGH ("trusted, <=0.041 dex") while E and D
sit in NIST_GRADE_CULL.

The claim was removed and the VALUE HELD, because no source can adjudicate it: Den Hartog
2019 stops at 4584 A, Melendez & Barbuy 2009 has both but flagged S (solar-fitted,
firewalled by RYA-161), and NIST itself is the scale RYA-853 scope 3 measured as the low
one. These tests pin that the claim is gone, the number did not move, and the gap is
declared rather than papered over.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
SOLAR = ROOT / "data" / "linelists" / "linelist_solar.csv"
NIST_REF = ROOT / "data" / "linelists" / "nist_reference.csv"
NIST_XC = ROOT / "data" / "linelists" / "nist_crosscheck.csv"
GAP = ROOT / "data" / "results" / "rya853" / "rya853_fe2_declared_gap.json"

#: (wavelength, EP, the value that must NOT move)
TARGETS = [(6149.246, 3.8892, -2.724), (6247.557, 3.8916, -2.329)]


@pytest.fixture(scope="module")
def canon():
    d = pd.read_csv(CANON, comment="#", low_memory=False)
    return d[d.species.astype(str) == "Fe II"]


def _row(df, wl, ep):
    m = df[((df.wavelength_air_A - wl).abs() <= 0.05)
           & ((df.excitation_potential_eV - ep).abs() <= 0.05)]
    assert len(m) == 1, f"{len(m)} rows at {wl} — the match itself is ambiguous"
    return m.iloc[0]


def test_the_fabricated_nist_grade_is_gone(canon):
    """🔴 The defect that opened the ticket. `B` is in NIST_GRADE_HIGH; NIST says E and D."""
    for wl, ep, _ in TARGETS:
        r = _row(canon, wl, ep)
        assert str(r.nist_grade) in ("nan", "", "None"), (
            f"{wl} still carries nist_grade={r.nist_grade!r}")
        assert "grade B" not in str(r.loggf_reference)


def test_the_value_did_not_move(canon):
    """🔴 CRITICAL. This was a provenance fix, not a re-measurement. Adopting NIST's value
    would have moved the ionization arbiter ~0.13 dex onto the scale scope 3 measured as
    the LOW one, on grade E/D data."""
    for wl, ep, held in TARGETS:
        assert float(_row(canon, wl, ep).log_gf) == held


def test_the_reference_states_what_nist_actually_says(canon):
    """A corrected label that merely goes quiet would lose the finding. The row has to
    carry the number it is NOT."""
    for wl, ep, _ in TARGETS:
        ref = str(_row(canon, wl, ep).loggf_reference)
        assert "NOT NIST" in ref and "RYA-853" in ref
        assert ("-2.854" in ref) or ("-2.444" in ref)


def test_no_cited_sigma_was_invented(canon):
    """The honest state for a line with no graded source is NO cited sigma — the band
    budget then charges the ungraded blanket. Attaching NIST's grade-E dex to a value that
    is not NIST's would be a second fabrication."""
    for wl, ep, _ in TARGETS:
        assert not np.isfinite(float(_row(canon, wl, ep).gf_sigma_dex or np.nan))


def test_the_stamped_grade_is_gone_from_the_upstream_files_too():
    """🔴 The defect lived in THREE files: the extracts are the source, crosscheck_nist()
    stamped the grade onto linelist_solar's VALD3 row, and migrate_gf_single_source built
    canonical_gf from them. Fixing only canonical_gf would be undone by the next build."""
    ll = pd.read_csv(SOLAR, low_memory=False)
    for wl, ep, _ in TARGETS:
        r = _row(ll[ll.element.astype(str) == "Fe"], wl, ep)
        assert str(r.nist_grade) in ("nan", "", "None"), f"linelist_solar {wl}: {r.nist_grade!r}"
        assert str(r.loggf_source) == "VALD3"
    for path in (NIST_REF, NIST_XC):
        d = pd.read_csv(path, comment="#")
        d = d[(d.element.astype(str) == "Fe") & (d.ion.astype(str) == "II")]
        for _, r in d.iterrows():
            if any(abs(float(r.wavelength_air_A) - wl) <= 0.05 for wl, _, _ in TARGETS):
                assert str(r.nist_grade).strip() in ("nan", "", "None"), path.name
                assert "RYA-853" in str(r.notes)


def test_the_declared_gap_names_what_was_checked():
    """RYA-833: an absence is a hypothesis to check. The gap record must say which sources
    were looked at and why each fails, not merely that none was found."""
    g = json.loads(GAP.read_text())
    checked = g["checked_not_assumed"]
    assert "DenHartog2019" in checked and "4584" in checked["DenHartog2019"]
    assert "FIREWALL" in checked["MelendezBarbuy2009"].upper()
    assert "NIST_ASD" in checked
    assert "HELD" in g["value_disposition"]


# ── scope 4: the two wavelength-only matchers ────────────────────────────────

def _load_build_linelist():
    spec = importlib.util.spec_from_file_location(
        "build_linelist_rya853", ROOT / "scripts" / "build_linelist.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _nist_csv(tmp_path, rows):
    p = tmp_path / "nist.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "ion", "wavelength_air_A", "excitation_potential_eV",
                    "log_gf", "nist_grade", "tier", "notes"])
        w.writerows(rows)
    return p


def test_control_a_wrong_level_inside_the_window_is_NOT_graded(tmp_path):
    """POSITIVE CONTROL, scope 4. Two levels of one species can sit inside 0.010 A. The old
    matcher took the first in FILE ORDER; a grade then lands on a line NIST never graded.
    RYA-853 measured this directly: on wavelength alone, 6149.246 picks up EP 13.436 eV
    instead of 3.889."""
    m = _load_build_linelist()
    nist = _nist_csv(tmp_path, [["Fe", "II", 6149.246, 13.436, -4.983, "B", 1, "wrong level"]])
    lines = [{"element": "Fe", "ion": "II", "wavelength_air_A": 6149.246,
              "excitation_potential_eV": 3.8892, "nist_grade": ""}]
    out = m.crosscheck_nist(lines, nist)
    assert out[0]["nist_grade"] == "", "a wrong-EP row was still graded"


def test_control_the_right_level_IS_still_graded(tmp_path):
    """The guard must not simply refuse everything — a matcher that never matches reads as
    'no graded source' and is a silent false absence."""
    m = _load_build_linelist()
    nist = _nist_csv(tmp_path, [["Fe", "II", 6149.246, 3.889, -2.724, "C", 1, "right level"]])
    lines = [{"element": "Fe", "ion": "II", "wavelength_air_A": 6149.246,
              "excitation_potential_eV": 3.8892, "nist_grade": ""}]
    assert m.crosscheck_nist(lines, nist)[0]["nist_grade"] == "C"


def test_control_two_candidates_on_both_axes_are_REFUSED(tmp_path):
    """Ambiguity is refused, never resolved by proximity. An argmin over rows the data
    cannot separate manufactures a provenance nobody measured."""
    m = _load_build_linelist()
    nist = _nist_csv(tmp_path, [
        ["Fe", "II", 6149.244, 3.889, -2.724, "B", 1, "a"],
        ["Fe", "II", 6149.248, 3.890, -2.700, "D", 1, "b"]])
    lines = [{"element": "Fe", "ion": "II", "wavelength_air_A": 6149.246,
              "excitation_potential_eV": 3.8892, "nist_grade": ""}]
    assert m.crosscheck_nist(lines, nist)[0]["nist_grade"] == ""


def test_control_a_line_with_no_EP_is_not_graded(tmp_path):
    """No EP on our side means the level cannot be verified, so no grade may be claimed."""
    m = _load_build_linelist()
    nist = _nist_csv(tmp_path, [["Fe", "II", 6149.246, 3.889, -2.724, "B", 1, "x"]])
    lines = [{"element": "Fe", "ion": "II", "wavelength_air_A": 6149.246, "nist_grade": ""}]
    assert m.crosscheck_nist(lines, nist)[0]["nist_grade"] == ""


def test_the_rya347_matcher_carries_an_EP_guard_and_the_EP_to_use_it():
    """The other instance of the same defect. Source-level, because running rya347 needs
    the Fe II fit products. Both halves are asserted: the guard, AND that step0 actually
    emits the ep_eV it consumes — a guard fed NaN would refuse everything and read as
    'no graded source'."""
    src = (ROOT / "scripts" / "rya347_fe2_atomic_data_audit.py").read_text()
    assert "ep_eV=float(m.excitation_potential_eV)" in src, "step0 no longer emits ep_eV"
    assert "EP_TOL_EV" in src and "excitation_potential_eV" in src
    assert ".abs() < 0.1]" not in src, "the wavelength-only +/-0.1 A match came back"
