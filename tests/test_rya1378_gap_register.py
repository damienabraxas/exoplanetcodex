"""RYA-1378: the Atomic Data Gap Register lists every line the grading keeps out of the
reference pool, for every element, from the ORCH's own pool definitions."""
import pandas as pd
import pytest

from pipeline import gap_register as g

HEAD = ["line_id", "physical_id", "species", "wavelength_air_A", "excitation_potential_eV",
        "log_gf", "loggf_reference", "gf_tier", "nist_grade", "gf_sigma_dex", "in_synth",
        "lab_source_tag"]
ROWS = [
    # VIS: one LAB line (reference) -> a graded route exists; two list lines held
    ("a", "pa", "Xx I", 5000.00, 2.0, -1.0, "LAB1", "LAB", "", 0.02, True, "lab"),
    ("b", "pb", "Xx I", 5100.00, 2.0, -1.5, "K16", "KURUCZ", "", None, True, ""),
    ("c", "pc", "Xx I", 5200.00, 2.0, -2.0, "VALD3", "VALD3", "", None, False, ""),
    # graded but not lab: codex pool only
    ("d", "pd", "Xx I", 5300.00, 2.0, -2.0, "NIST", "NIST-C+", "B", None, True, ""),
    # H band: no graded line at all -> priority 1 even out of the synthesis list
    ("e", "pe", "Xx I", 16000.0, 5.0, -0.5, "K16", "KURUCZ", "", None, False, ""),
    # a LAB line the registry excludes
    ("f", "pf", "Xx II", 5400.00, 3.0, -1.0, "LAB2", "LAB", "", 0.02, True, "lab"),
    # a higher ion stage: not an ORCH ion
    ("g", "pg", "Xx III", 5500.00, 9.0, -1.0, "K16", "KURUCZ", "", None, False, ""),
]


@pytest.fixture
def reg(tmp_path, monkeypatch):
    canon = tmp_path / "canonical_gf.csv"
    pd.DataFrame(ROWS, columns=HEAD).to_csv(canon, index=False)
    known = tmp_path / "known.csv"
    pd.DataFrame([
        ("Xx", "Xx I", 5000.00, "GF_DISCREPANCY", "lab vs list", "T-1"),   # has a LAB row
        ("Xx", "Xx I", 5100.00, "GF_DISCREPANCY", "theory vs list", "T-2"),
        ("Xx", "Xx II", 4000.00, "GF_DISCREPANCY", "lab line absent", "T-3"),
    ], columns=["element", "species", "wavelength_air_A", "flags", "evidence", "source"]
    ).to_csv(known, index=False)
    monkeypatch.setattr(g, "CANON", canon)
    monkeypatch.setattr(g, "KNOWN", known)
    moves = tmp_path / "moves.csv"
    pd.DataFrame([("b", "Xx I", 5100.0, 2.0, -1.5, "KURUCZ", -0.9, "lab2099", 0.6)],
                 columns=["line_id", "species", "wavelength_air_A", "ep_eV", "canonical_loggf",
                          "canonical_tier", "candidate_loggf", "candidate_source", "d_dex"]
                 ).to_csv(moves, index=False)
    monkeypatch.setattr(g, "HELD_MOVES", moves)
    monkeypatch.setattr(g, "PROBLEMS", tmp_path / "{star}_problem_lines.csv")
    monkeypatch.setattr(g, "DIAG", tmp_path / "{star}_{el}_diag.csv")
    monkeypatch.setattr(g, "_excluded_by_registry",
                        lambda sym, rows: rows.wavelength_air_A.astype(float).eq(5400.0))
    monkeypatch.setattr(g, "_solar_value", lambda sym: 5.0)
    return g.build("Xx")


def row(t, w):
    r = t[t.wavelength_air_A.astype(float) == w]
    assert len(r) == 1, w
    return r.iloc[0]


def test_reference_lines_are_not_listed_and_higher_ions_are_dropped(reg):
    ws = set(reg.wavelength_air_A.astype(float))
    assert 5000.0 not in ws and 5500.0 not in ws


def test_held_reason_follows_the_pool_definitions(reg):
    assert row(reg, 5100.0).held_reason == "UNGRADED_GF"
    assert row(reg, 5300.0).held_reason == "GRADED_NOT_LAB"
    ex = row(reg, 5400.0)
    assert ex.held_reason == "REGISTRY_EXCLUDED" and "UNRATED_GF" not in ex["flags"]


def test_absence_is_judged_against_every_row_including_lab(reg):
    missing = reg[reg["flags"].str.contains("MISSING_GF")]
    assert list(missing.wavelength_air_A.astype(float)) == [4000.0]
    assert row(reg, 4000.0).held_reason == "NOT_IN_CANONICAL"


def test_flags_coexist(reg):
    r = row(reg, 5100.0)
    assert {"UNRATED_GF", "GF_DISCREPANCY"} <= set(r["flags"].split(";")) and r.priority == 2
    assert "PROVENANCE_GAP" in row(reg, 5200.0)["flags"]


def test_priority_one_is_a_band_without_a_graded_line_whatever_in_synth(reg):
    r = row(reg, 16000.0)
    assert r.band == "H" and not r.in_synth and r.no_graded_route_in_band and r.priority == 1
    assert row(reg, 5200.0).priority == 4      # not in synthesis, band has a graded route


def test_the_register_is_deterministic(reg, monkeypatch):
    assert g.build("Xx").equals(reg)


def test_band_of_takes_the_narrowest_band():
    assert g.band_of(16000.0) == "H" and g.band_of(10000.0) == "NIR"


def test_a_held_large_lab_gf_move_is_a_gf_discrepancy(reg):
    r = row(reg, 5100.0)
    assert r["flags"].count("GF_DISCREPANCY") == 1          # curated + held move: one flag
    assert "lab2099" in r.evidence and "theory vs list" in r.evidence and r.priority == 2
