"""
tests/test_al_level_reach_rya773.py
===================================
RYA-773 — the Al departure extract now serves the lines the pipeline actually measures.

WHAT THIS PINS
--------------
The extract carried 6696/6698 while Al's clean, multi-arm-corroborated lines are
7835/7836 + 8772/8773 (RYA-708/716): Al's best data and Al's NLTE coverage did not
overlap, so three of four in-aggregate lines were dispositioned ENGINE-A UNCOVERED and
the Engine-A product rested on a single line.

The fix was cheap because the levels were already there. These tests pin the reasoning
that established that — the level identification, which is the step where a silent
error would live — plus the extract's new coverage. The identification logic is pure
(energies and selection rules), so it is testable without PySME or the 2.3 GB `.grd`;
the delta VALUES are a Sirius derivation and are recorded on the ticket, not asserted
here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(ROOT / "scripts"))
import rya773_al_level_reach as LR                     # noqa: E402

EXTRACT = ROOT / "data" / "nlte_grids" / "Al_Amarsi2020_PySME.csv"


# ── the atom's level table ───────────────────────────────────────────────────

def test_level_table_parses_the_variable_width_configuration():
    """`3s2.nd  y 2D` carries a series letter, so the configuration field is not a
    fixed number of tokens. A plain split() shears that row into the wrong columns."""
    lv = LR.read_levels()
    odd = lv[lv.conf.str.contains("nd")]
    assert len(odd) == 2                                    # 3s2.nd y 2D, J=1.5 and 2.5
    assert set(odd.term) == {"y 2D"}
    assert odd.E_eV.min() == pytest.approx(4.826632, abs=1e-6)
    # and the ordinary rows are unharmed
    assert lv.loc[lv.idx == 7, "conf"].iloc[0] == "3s2.3d"
    assert lv.loc[lv.idx == 7, "term"].iloc[0] == "2D"


def test_the_doublet_levels_are_in_the_atom():
    """The whole cheap-vs-expensive verdict. 3d 2D, 5f 2F* and 6f 2F* must all exist."""
    lv = LR.read_levels()
    for conf, term, n in [("3s2.3d", "2D", 2), ("3s2.5f", "2F*", 2), ("3s2.6f", "2F*", 2)]:
        got = lv[(lv.conf == conf) & (lv.term == term)]
        assert len(got) == n, f"{conf} {term} missing from the model atom"


# ── the E1 selection rule, which is what disambiguates the upper level ───────

def test_orbital_l_reads_the_last_subshell():
    assert LR.orbital_l("3s2.3d") == 2
    assert LR.orbital_l("3s2.6f") == 3
    assert LR.orbital_l("3s2.5g") == 4
    assert LR.orbital_l("3s2.3p") == 1


def test_2d_to_2f_allowed_but_2d_to_2g_forbidden():
    """This is the tie-breaker: 6g 2G sits only 3.4 meV from 6f 2F*, so an energy match
    alone cannot separate them — but l=2 -> l=4 is not an electric-dipole transition."""
    lv = LR.read_levels()
    lo = lv[(lv.conf == "3s2.3d") & (lv.J == 2.5)].iloc[0]
    f_up = lv[(lv.conf == "3s2.6f") & (lv.J == 3.5)].iloc[0]
    g_up = lv[(lv.conf == "3s2.6g") & (lv.J == 3.5)].iloc[0]
    assert LR.e1_allowed(lo, f_up) is True
    assert LR.e1_allowed(lo, g_up) is False
    assert abs(float(f_up.E_eV) - float(g_up.E_eV)) < 0.004      # the near-degeneracy


def test_parity_and_delta_j_are_enforced():
    lv = LR.read_levels()
    lo = lv[(lv.conf == "3s2.3d") & (lv.J == 1.5)].iloc[0]       # even, J=1.5
    same_parity = lv[(lv.conf == "3s2.4s")].iloc[0]              # even -> forbidden
    too_far_J = lv[(lv.conf == "3s2.6f") & (lv.J == 3.5)].iloc[0]  # dJ = 2 -> forbidden
    assert LR.e1_allowed(lo, same_parity) is False
    assert LR.e1_allowed(lo, too_far_J) is False


def test_match_level_picks_the_allowed_term_over_the_nearer_forbidden_one():
    """With a tolerance deliberately wider than the 6f/6g gap, the selection rule — not
    the tolerance — must be what returns 2F*."""
    lv = LR.read_levels()
    lo = lv[(lv.conf == "3s2.3d") & (lv.J == 1.5)].iloc[0]
    up, resid, margin = LR.match_level(lv, 5.6034263, lower=lo, tol_eV=0.02)
    assert up is not None and up.conf == "3s2.6f" and up.term == "2F*"
    assert resid < 1e-5
    assert margin == pytest.approx(0.00335, abs=5e-4)            # to the 6g competitor


def test_match_level_refuses_when_nothing_is_within_tolerance():
    lv = LR.read_levels()
    up, resid, _ = LR.match_level(lv, 99.0, tol_eV=0.006)
    assert up is None and resid > 90


# ── what the extract now serves ──────────────────────────────────────────────

def _served():
    return np.unique(pd.read_csv(EXTRACT).wave_A.to_numpy())


def test_extract_serves_the_rya773_pairs_and_names_new_uncovered_lines():
    """RYA-925 expanded the LTE pool after RYA-773.

    Coverage is a property of each model product, not a veto on the LTE measurement
    pool.  Pin the newly exposed reduced coverage so later code cannot silently claim
    Amarsi served those four lines or silently drop them from the LTE product.
    """
    meas = LR.measured_lines()
    served = _served()
    uncovered = [float(w) for w in meas[meas.in_aggregate].wavelength_air_A
                 if np.min(np.abs(served - float(w))) > 0.15]
    assert uncovered == [5557.062, 6783.638, 7361.568, 7362.296]
    for wave in (6696.0137, 6698.6715, 7835.309, 7836.134, 8772.865, 8773.8975):
        assert np.min(np.abs(served - wave)) <= 0.15


def test_the_banked_pair_is_untouched():
    """The 22 RYA-402 rows are anchor-validated. Extending must not move them — a
    re-derivation that silently shifted a registered correction is the failure this
    guards (the new rows were appended; the old ones were never recomputed)."""
    d = pd.read_csv(EXTRACT)
    old = d[d.wave_A.isin([6696.023, 6698.673])]
    assert len(old) == 22
    solar = old[(old.teff_K == 5772) & (old.logg == 4.44) & (old.feh == 0.0)]
    assert dict(zip(solar.wave_A, solar.delta_nlte)) == {6696.023: -0.0275,
                                                         6698.673: -0.0171}


def test_new_deltas_are_small_and_negative_like_the_published_atom():
    """Nordlander & Lind 2017: subordinate Al I corrections are near-zero to slightly
    negative. Reproduced, never tuned — so this is a SIGN-and-ORDER check, deliberately
    loose, not a re-assertion of the derived numbers."""
    d = pd.read_csv(EXTRACT)
    new = d[~d.wave_A.isin([6696.023, 6698.673])]
    assert len(new) == 44
    assert (new.delta_nlte < 0).all(), "NLTE strengthens these lines -> delta negative"
    assert new.delta_nlte.min() > -0.10, "an order larger than the published band"


def test_every_registered_al_diagnostic_is_in_the_extract():
    """NLTE_LINES and the emitted extract must not drift apart — a line registered as a
    diagnostic but absent from the CSV is a correction the pipeline cannot apply."""
    from pipeline.pysme_nlte import NLTE_LINES
    served = _served()
    for row in NLTE_LINES["Al"]:
        assert np.min(np.abs(served - row[0])) <= 0.01, f"{row[0]} not in the extract"


# ── the derivation machinery the extension needed ────────────────────────────

def test_line_subsetting_selects_by_wavelength_and_by_row():
    from pipeline.pysme_nlte import _select_lines, NLTE_LINES
    all_al = NLTE_LINES["Al"]
    assert _select_lines("Al", None) == list(all_al)
    got = _select_lines("Al", [8772.865, 8773.897])
    assert [r[0] for r in got] == [8772.865, 8773.897]
    assert _select_lines("Al", [all_al[0]]) == [all_al[0]]


def test_line_subsetting_raises_on_an_unregistered_wavelength():
    """A typo'd wavelength must not quietly derive a smaller set — that looks exactly
    like a successful derivation of the lines you asked for."""
    from pipeline.pysme_nlte import _select_lines
    with pytest.raises(KeyError, match="no registered NLTE diagnostic"):
        _select_lines("Al", [8772.865, 1234.5])


def test_components_may_carry_their_own_upper_level():
    """7836.134 and 8773.9 are each two J-components from one lower level, and the grid
    says their upper levels do NOT share departures (5f 2F* J=2.5 vs 3.5 differ by 1.8%
    in b). So each component must be emitted with its own upper label rather than
    inheriting the dominant one through the HFS channel."""
    from pipeline.pysme_nlte import _linelist_rows, NLTE_LINES
    feature = [r for r in NLTE_LINES["Al"] if r[0] == 7836.134][0]
    rows = _linelist_rows("Al", [feature])
    assert len(rows) == 2
    assert sorted(r["j_up"] for r in rows) == [2.5, 3.5]
    assert {r["term_upper"] for r in rows} == {"3s2.6f 2F*"}
    assert {r["term_lower"] for r in rows} == {"3s2.3d 2D"}
    assert sorted(r["gflog"] for r in rows) == [-1.795, -0.494]


def test_plain_hfs_components_still_inherit_the_feature_label():
    """The 2-tuple component form is unchanged — true hyperfine structure splits a level
    by ~ueV, so those components DO share a departure coefficient."""
    from pipeline.pysme_nlte import _linelist_rows
    line = (5782.130, -1.488, 1.642, 1.5, 3.786, 0.5, '3d9.4s2 2D', '3d10.4p 2P*', -7.79,
            [(5782.10, -1.8), (5782.16, -1.9)])
    rows = _linelist_rows("Cu", [line])
    assert len(rows) == 2
    assert all(r["term_upper"] == "3d10.4p 2P*" and r["j_up"] == 0.5 for r in rows)


def test_grid_labels_decode_without_nul_padding():
    """The `.grd` stores level labels as fixed-width NUL-padded fields and str.strip()
    does not remove NUL. Measured on Sirius, PySME normalises the padding itself, so
    this is hygiene rather than a fix — pinned so the decoder does not regress to
    emitting labels no human comparison will match."""
    from pipeline.pysme_nlte import decode_grid_label
    assert decode_grid_label(b"3s2.5f\x00\x00") == "3s2.5f"
    assert decode_grid_label(b"2F*\x00") == "2F*"
    assert decode_grid_label(b"y 2D  ") == "y 2D"
    assert decode_grid_label(np.bytes_(b"2D\x00\x00")) == "2D"
