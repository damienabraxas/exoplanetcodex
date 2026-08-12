"""RYA-794 — CRIRES+ Vesta solar Y band.

These pin the two things that were actually hard here and that a future edit could break
silently: the fixed-width byte offsets into Elgueta's atomicy.dat, and the two conflicting
wavelength-unit conventions inside one VizieR ReadMe.

The Sirius-only tests skip cleanly on the Mac, where the 33 MB sp/ spectra are not staged.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "scripts" / "normalize_vesta_ir.py"


def _mod():
    s = importlib.util.spec_from_file_location("normalize_vesta_ir", SPEC)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


m = pytest.importorskip("numpy") and _mod()
has_atomicy = m.ATOMICY.exists()
sirius = pytest.mark.skipif(not has_atomicy, reason="Elgueta VizieR pull is Sirius-only")


@sirius
def test_atomicy_is_fixed_width_805():
    """Every offset in this reader assumes Lrecl 805 blank-padded records."""
    lens = {len(l) for l in m.ATOMICY.read_text(errors="replace").splitlines() if l.strip()}
    assert lens == {805}, lens


@sirius
def test_gdrob_certifies_zero_fe_i_for_a_g_dwarf():
    """THE headline finding, and the one most likely to be mistaken for a bug.

    Elgueta's own G-dwarf robust flag selects no Fe I line at all in the Y band. Pinned
    with its corroborating counts so that a future reader sees a deliberate result rather
    than an empty list to 'fix'.
    """
    assert len(m.elgueta_y_lines("FeI", tier="all")) == 141
    assert len(m.elgueta_y_lines("FeI", tier="assessed")) == 25
    assert m.elgueta_y_lines("FeI", tier="robust") == []
    # ...while other species DO get certified, so the flag itself is being read correctly
    robust_species = {r["element"] for sp in ("SiI", "CI", "CrI", "TiI")
                      for r in m.elgueta_y_lines(sp, tier="robust")}
    assert robust_species == {"SiI", "CI", "CrI", "TiI"}
    assert len(m.elgueta_y_lines("SiI", tier="robust")) == 18


@sirius
def test_gdsat_is_never_n_so_blank_means_uncertified():
    """The reason Fe I fails: GDSat is 'Y' or blank, never 'N'.

    If this ever shows an 'N', the flag semantics changed and the candidate tier below
    (which treats blank as *uncertified*, not *failed*) needs revisiting.
    """
    vals = {l[m.GD["sat"]] for l in m.ATOMICY.read_text(errors="replace").splitlines()
            if len(l) == 805}
    assert vals == {"Y", " "}


@sirius
def test_blank_gdsat_among_depth_passers_is_iron_only():
    """Blank GDSat means FAILED, not UNEVALUATED — and this is the evidence.

    GDSat never takes 'N', so blank is ambiguous in isolation. But of the lines that pass
    the depth cut, the ones with blank GDSat are iron and only iron. A merely sparse column
    would not single out one element, so the candidate tier is Elgueta-rejected science we
    choose to carry, not a hole in their table.
    """
    lines = [l for l in m.ATOMICY.read_text(errors="replace").splitlines() if len(l) == 805]
    depth_pass = [l for l in lines if l[m.GD["depth"]] == "Y"]
    blank_sat = [l for l in depth_pass if l[m.GD["sat"]] != "Y"]
    assert len(depth_pass) == 71
    assert {l[13:17].strip() for l in blank_sat} == {"FeI", "FeII"}
    assert len(blank_sat) == 13


@sirius
def test_candidate_tier_is_disjoint_from_robust():
    cand = m.elgueta_y_lines("FeI", tier="candidate")
    assert len(cand) == 7
    assert all(r["gd_robust"] != "Y" for r in cand)
    assert all(r["gd_depth"] == r["gd_purity"] == r["gd_gof"] == "Y" for r in cand)
    assert all(r["tier"] == "candidate" for r in cand)


@sirius
def test_the_two_wavelength_conventions_do_not_get_swapped():
    """One ReadMe, one label ('0.1nm'), two different actual units.

    atomicy.dat really is Angstrom; sp/ is nanometres and needs x10. Reading either by the
    other's convention yields a confident wrong answer, so both are asserted.
    """
    fe = m.elgueta_y_lines("FeI", tier="all")
    assert 9800.0 < fe[0]["wave_A"] < 10800.0, "atomicy is already Angstrom"
    assert m.NM_TO_A == 10.0
    if m.SP_Y.exists():
        w, _ = m.load_y_band()
        assert 9796.0 < w.min() < 9797.0 and 10796.0 < w.max() < 10797.0


@pytest.mark.skipif(not m.SP_Y.exists(), reason="sp/ spectra are Sirius-only")
def test_selected_lines_are_present_at_their_stated_depth():
    """End-to-end: the 5 candidates land in OUR spectrum at Elgueta's published depths.

    This is the test that would catch a unit slip, a wavelength-solution error or a
    continuum error, because the reference depths come from the paper, not from us.
    """
    import numpy as np
    w, f = m.load_y_band()
    sel = (w >= 10280.0) & (w <= 10680.0)
    norm, _cont, order = m.normalize_reflectance(w[sel], f[sel])
    ws = w[sel][order]
    rows = [r for r in m.elgueta_y_lines("FeI", tier="candidate")
            if 10280.0 <= r["wave_A"] <= 10680.0]
    assert len(rows) == 5
    m.verify_line_depths(ws, norm, rows)
    for r in rows:
        assert abs(r["meas_dlam_A"]) < 0.05, r
        assert 0.8 < r["depth_ratio"] < 1.25, r


@sirius
def test_no_abundance_product_is_written():
    """The ticket's CRITICAL constraint, enforced rather than remembered.

    Single exposure, zero certified lines: this script must never emit an A(Fe).
    """
    src = SPEC.read_text()
    assert "NO ABUNDANCE IS QUOTED" in src, "the refusal must stay visible in the output"
    if m.OUT_DIR.exists():
        stray = [p.name for p in m.OUT_DIR.iterdir()
                 if "abund" in p.name.lower() or p.name.endswith("_products.csv")]
        assert not stray, f"RYA-794 must not emit an abundance product: {stray}"
