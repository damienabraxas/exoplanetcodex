"""
tests/test_fe_ir_gf_adjudication_rya780.py
==========================================
RYA-780 — the Fe I IR gf-source adjudication.

WHAT IS WORTH PINNING HERE
--------------------------
The scientific result (14 of 14 disputed lines read high against a validated primary
scale) depends on a measured band product that does not live on main, so it is recorded
on the ticket and in the emitted artifact rather than asserted here. What these tests
pin is the machinery that could silently produce a WRONG result and still look right:

  * transition IDENTITY — matching a primary measurement to a pool line on wavelength
    alone pairs different transitions. That is not hypothetical: it manufactured a
    2.85 dex "discrepancy" at 8876 Å before the EP filter went in.
  * the SOURCE classification — if FMW were ever classed as primary, the reference scale
    would absorb the offset under test and the whole comparison would quietly collapse.
  * the VALD parser — it is what supplies both the level identity and VALD's independent
    source attribution, and it reads a fiddly 4-record format.
  * the emitted artifact's shape and firewall — the RYA-161 rule that a quarantine reason
    cites primary-source provenance and never delta-vs-anchor.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import rya780_fe_ir_gf_adjudication as ADJ                     # noqa: E402

DISPOSITION = ROOT / "data" / "results" / "fe_ir_gf" / "FeI_IR_gf_disposition.csv"


# ── source classification ────────────────────────────────────────────────────

def test_the_disputed_compilations_are_not_classed_as_primary():
    """If FMW counted as primary it would enter the reference scale, absorb the offset
    under test, and the comparison would silently return ~0."""
    for bad in ("FMW", "GESB82c", "GESB82c+BW"):
        assert not bad.startswith(ADJ.PRIMARY_PREFIXES), f"{bad} must not be primary"
    import re
    for bad in ("FMW", "GESB82c", "GESB82c+BW"):
        assert re.search(ADJ.DISPUTED_PATTERN, bad)


def test_the_laboratory_sources_are_classed_as_primary():
    for good in ("BWL", "BK", "BKK", "BK+BWL", "BWL+GESHRL", "2014MNRAS.", "GESHRL14"):
        assert good.startswith(ADJ.PRIMARY_PREFIXES), f"{good} must be primary"


def test_k07_is_neither_primary_nor_disputed():
    """Kurucz semi-empirical is explicitly out of scope (RYA-709), so it must neither
    define the reference scale nor be adjudicated as if it were FMW."""
    import re
    assert not "K07".startswith(ADJ.PRIMARY_PREFIXES)
    assert not re.search(ADJ.DISPUTED_PATTERN, "K07")


# ── transition identity ──────────────────────────────────────────────────────

def _prim(rows):
    return pd.DataFrame(rows, columns=["wave_A", "ep_eV", "loggf", "e_loggf",
                                       "prev", "prev_ref"])


def test_a_primary_match_requires_the_lower_level_too():
    """The real 8876 case: two transitions 0.018 A apart with lower levels 0.44 eV apart.
    Wavelength alone accepts it and invents a 2.85 dex discrepancy."""
    prim = _prim([[8876.0241, 5.020, -1.05, 0.03, np.nan, ""]])
    assert ADJ.match_primary(prim, 8876.0059, 5.020) is not None      # right level
    assert ADJ.match_primary(prim, 8876.0059, 4.584) is None          # our K07 line


def test_a_primary_match_requires_the_wavelength_too():
    prim = _prim([[7016.3920, 4.154, -1.07, 0.03, np.nan, ""]])
    assert ADJ.match_primary(prim, 7016.3920, 4.154) is not None
    assert ADJ.match_primary(prim, 7443.0217, 4.154) is None


def test_the_nearest_primary_wins_when_several_qualify():
    prim = _prim([[7016.30, 4.154, -1.20, 0.03, np.nan, ""],
                  [7016.39, 4.154, -1.07, 0.03, np.nan, ""]])
    got = ADJ.match_primary(prim, 7016.392, 4.154)
    assert got is not None and got.loggf == pytest.approx(-1.07)


# ── the vendored primary sources ─────────────────────────────────────────────

def test_den_hartog_parses_with_energies_in_eV():
    """Den Hartog tabulates the lower level in cm^-1; the pool uses eV. A unit slip here
    would break every identity check silently — the match would just always fail."""
    d = ADJ.read_primary("denhartog2014_FeI_6900_9250.tsv")
    assert len(d) >= 30
    assert d.wave_A.between(6900, 9250).all()
    assert d.ep_eV.between(2.0, 6.0).all(), "lower levels should be a few eV, not cm^-1"
    assert d.loggf.between(-4, 1).all()


def test_ruffoni_parses_with_energies_already_in_eV():
    d = ADJ.read_primary("ruffoni2014_FeI_table5.tsv")
    assert len(d) >= 30
    assert d.ep_eV.between(0.0, 6.0).all()


def test_no_primary_covers_any_disputed_transition():
    """The finding that decides the whole ticket: RECOVERED is unavailable because the
    measurements do not exist, not because of a judgement call. If a future pull changes
    that, this test SHOULD fail — that is a recovery, and it is the good outcome."""
    prims = {n: ADJ.read_primary(f) for n, f in
             (("DenHartog2014", "denhartog2014_FeI_6900_9250.tsv"),
              ("Ruffoni2014", "ruffoni2014_FeI_table5.tsv"))}
    d = pd.read_csv(DISPOSITION)
    for _, r in d.iterrows():
        for nm, p in prims.items():
            assert ADJ.match_primary(p, r.wave_A, r.ep_eV) is None, (
                f"{nm} now covers {r.wave_A} — re-source it (RYA-354 RECOVERED)")


# ── the emitted disposition ──────────────────────────────────────────────────

def test_disposition_covers_the_disputed_set_and_nothing_else():
    d = pd.read_csv(DISPOSITION)
    assert len(d) == 14
    assert d.pool_source.str.contains("FMW|GESB82c").all()
    assert d.wave_A.is_unique


def test_vald_is_not_an_independent_opinion_on_these_lines():
    """The reason a catalogue cannot referee this: VALD names FMW for every disputed line
    and carries the same number. Pinned because it is the premise of the whole approach."""
    d = pd.read_csv(DISPOSITION)
    assert (d.vald_source == "FMW").all()
    assert (d.pool_loggf - d.vald_loggf).abs().max() < 0.005


def test_the_two_gesb82c_lines_are_fmw_in_vald():
    """RYA-760 reported GESB82c as a separate (underpowered, n=2) population. It is not —
    both lines are FMW in VALD, so the disputed set is one population of 14."""
    d = pd.read_csv(DISPOSITION)
    g = d[d.pool_source.str.contains("GESB82c")]
    assert len(g) == 2 and (g.vald_source == "FMW").all()


def test_every_reason_cites_provenance_and_never_the_anchor():
    """The RYA-161 firewall. A quarantine justified by distance from 7.466 would be
    tuning; these must cite the primary-source comparison."""
    d = pd.read_csv(DISPOSITION)
    for _, r in d.iterrows():
        assert "7.466" not in r.reason and "anchor" not in r.reason.lower()
        assert ("primary" in r.reason.lower()) or ("Den Hartog" in r.reason)


def test_quarantines_are_the_lines_that_read_high_not_an_arbitrary_set():
    d = pd.read_csv(DISPOSITION)
    q = d[d.disposition == "QUARANTINED-SCALE-EVIDENCE"]
    keep = d[d.disposition == "NO-PRIMARY-NO-EVIDENCE-AGAINST"]
    assert len(q) and len(keep)
    assert q.n_sigma.min() >= 2.0
    assert keep.n_sigma.max() < 2.0


def test_no_line_was_recovered_and_the_artifact_says_why():
    d = pd.read_csv(DISPOSITION)
    assert not (d.disposition == "RECOVERED").any()
    assert d.primary_source.isna().all() or (d.primary_source.astype(str) == "").all() \
        or (d.primary_source.fillna("") == "").all()
    assert d.reason.str.contains("No primary laboratory measurement").all()


def test_the_whole_disputed_population_reads_high():
    """The result that survives the per-line threshold being arbitrary: every one of the
    14 residuals is positive, which is a 6e-05 sign test on its own."""
    d = pd.read_csv(DISPOSITION)
    assert (d.residual_dex > 0).all()
    assert d.residual_dex.median() > 0.2
