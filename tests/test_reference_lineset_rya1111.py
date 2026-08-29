"""RYA-1111 - the reference line-set ingest, its vocabulary, and its refusals.

The measuring run itself needs Sirius (atlases + a synthesiser; MOOGSILENT is a Linux ELF),
so what is guarded here is everything up to the fit: the ingest, the normalisation, the
lambda+EP join at a SOURCE-DERIVED tolerance, the loud refusals, and the promise that none
of this touches an existing product.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import line_match, model_registry  # noqa: E402
from pipeline import reference_lineset as rls  # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"


# ── the vocabulary ────────────────────────────────────────────────────────────

def test_the_vocabulary_has_exactly_one_definition():
    """`reference_lineset` imports LINE_SETS; it does not restate it. Two copies of a
    vocabulary are two things free to disagree -- the defect the registry exists to end."""
    assert rls.LINE_SETS is model_registry.LINE_SETS


def test_the_spec_spelling_won_over_my_earlier_guess():
    """RYA-1101 opened the column before the axis had an owner and guessed
    `asplund-graded`. RYA-1111 owns it and says `asplund`."""
    assert "asplund" in rls.LINE_SETS
    assert "asplund-graded" not in rls.LINE_SETS
    for v in ("gbs", "our-graded", "our-deep-graded"):
        assert v in rls.LINE_SETS


def test_every_declared_set_uses_a_vocabulary_name():
    for name, spec in rls.SETS.items():
        assert name == spec.name
        assert spec.name in rls.LINE_SETS


# ── ingest ────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name,n,species", [("asplund", 53, {"Fe I": 40, "Fe II": 13}),
                                            ("gbs", 159, {"Fe I": 150, "Fe II": 9})])
def test_both_sets_load_into_one_schema(name, n, species):
    d = rls.load(name)
    assert len(d) == n
    assert dict(d.species.value_counts()) == species
    for c in ("line_set", "native_line_set", "species", "wavelength_air_A", "elo_eV",
              "loggf", "gf_source", "source", "gf_missing"):
        assert c in d.columns
    assert set(d.line_set.unique()) == {name}


def test_the_native_tag_is_recorded_not_rewritten():
    """The RYA-1109 artifact's own column says `asplund_agss21`. That is carried as
    `native_line_set`; the file is never relabelled."""
    d = rls.load("asplund")
    assert set(d.native_line_set.unique()) == {"asplund_agss21"}
    assert set(d.line_set.unique()) == {"asplund"}


def test_an_unknown_set_is_refused_not_guessed():
    with pytest.raises(rls.ReferenceLineSetError, match="unknown line set"):
        rls.load("asplund2005")


def test_a_mislabelled_file_is_refused(tmp_path, monkeypatch):
    """If the adapter and the artifact disagree about what the file IS, refuse -- do not
    relabel it into agreement."""
    import pandas as pd
    d = pd.read_csv(rls.SETS["asplund"].path)
    d["line_set"] = "something_else"
    p = tmp_path / "mislabelled.csv"
    d.to_csv(p, index=False)
    bad = rls.ReferenceSet(**{**rls.SETS["asplund"].__dict__, "path": p})
    monkeypatch.setitem(rls.SETS, "asplund", bad)
    with pytest.raises(rls.ReferenceLineSetError, match="disagree about"):
        rls.load("asplund")


# ── the tolerance: the RYA-1109 trap ──────────────────────────────────────────

def test_each_set_declares_a_derived_tolerance_with_its_basis():
    """🔴 THE RYA-1109 TRAP. The module default (0.005 A) suits a table printed to 0.01 A.
    AGSS21 prints nm to 2 dp -- 0.1 A -- so on the default its Fe I overlap read 2/40
    against a true 19/40. A tolerance is a property of the SOURCE."""
    asp, gbs = rls.SETS["asplund"], rls.SETS["gbs"]
    assert asp.match_tol_A == 0.05 != line_match.MATCH_TOL_A
    assert gbs.match_tol_A == 0.015 != line_match.MATCH_TOL_A
    assert asp.match_tol_A != gbs.match_tol_A, "two sources, two precisions"
    for s in rls.SETS.values():
        assert s.tol_basis.strip(), "a tolerance without a stated basis is a magic number"


def test_coverage_ships_a_plateau_so_the_count_can_be_checked():
    """A real overlap plateaus; a count of coincidences keeps climbing. Shipping the sweep
    is what lets a reader see which one they are looking at."""
    import pandas as pd
    pool = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    pool = pool[pool.species == "Fe I"].dropna(
        subset=["wavelength_air_A", "excitation_potential_eV"])
    cov = rls.coverage("asplund", pool, pool_label="canonical_gf Fe I")
    assert cov["match_tol_A"] == 0.05
    pl = cov["plateau"]
    assert pl["0.005"] < pl["0.05"], "the default really does undercount this source"
    assert pl["0.05"] == pl["0.1"] == pl["0.25"], "the count must plateau"


# ── loud refusals ─────────────────────────────────────────────────────────────

def test_gbs_lines_without_a_published_gf_are_flagged_not_dropped():
    """🔴 21 GBS lines read 'NOT PUBLISHED IN TABLES 4/5'. They stay in the set, flagged,
    and are refused at measurement. A `heiter2021_log_gf` IS staged for all 21 and taking
    it silently would turn 'GBS's own scale' into 'GBS where published, GES elsewhere' and
    report one number for the mixture. That adoption is RYA-1110's OPEN decision."""
    gap = rls.gf_gap("gbs")
    assert gap["n_ref"] == 159
    assert gap["n_without_published_gf"] == 21
    assert gap["n_measurable"] == 138
    assert gap["by_species"] == {"Fe I": 20, "Fe II": 1}
    d = rls.load("gbs")
    assert len(d) == 159, "flagged, NOT dropped"
    assert len(rls.measurable(d)) == 138
    assert "OPEN DECISION" in gap["disposition"] or "OPEN" in gap["disposition"]


def test_asplund_has_no_gf_gap():
    assert rls.gf_gap("asplund")["n_without_published_gf"] == 0


def test_a_partial_gf_override_is_refused(monkeypatch):
    """Measuring some lines on their scale and the rest on ours, then reporting one number
    for the mixture, is the confound a replication exists to remove (RYA-429)."""
    import numpy as np
    ll = np.array([(5000.0, 1.0, "Fe 1", -1.0)],
                  dtype=[("wave_A", "f8"), ("lower_state_eV", "f8"),
                         ("element", "U8"), ("loggf", "f8")])
    targets = rls.load("asplund").head(3)          # none of these is at 5000 A
    with pytest.raises(rls.ReferenceLineSetError, match="PARTIAL override"):
        rls.apply_gf_override(ll, targets, rls.SETS["asplund"])


# ── the feed axis ─────────────────────────────────────────────────────────────

def test_the_axis_is_derived_for_our_own_products_over_the_whole_live_feed():
    """Our products already say their pool in `tier`; storing it again would create a
    second source of truth free to disagree with the first."""
    feed = json.loads(FEED.read_text())
    got = [rls.line_set_for_product(p) for p in feed["products"]]
    assert set(got) == {"our-graded", "our-deep-graded"}
    assert len(got) == len(feed["products"])


def test_an_unrecognised_tier_is_refused_rather_than_defaulted():
    """`treatment_axes` warns that a derivation correct today goes silently wrong the day
    the correlation breaks. The guard is to fail loudly, not to pick the nearest value."""
    with pytest.raises(rls.ReferenceLineSetError, match="refusing to guess"):
        rls.line_set_for_product({"tier": "CONSISTENT"})


def test_a_stored_line_set_wins_and_must_be_in_the_vocabulary():
    assert rls.line_set_for_product({"line_set": "asplund", "tier": "GRADED"}) == "asplund"
    with pytest.raises(rls.ReferenceLineSetError, match="not in the vocabulary"):
        rls.line_set_for_product({"line_set": "asplund-graded"})


def test_tag_product_stamps_the_set_and_its_provenance():
    out = rls.tag_product({"A": 7.46}, "asplund")
    assert out["line_set"] == "asplund"
    assert "Asplund" in out["line_set_source"] and out["line_set_ticket"] == "RYA-1109"
    assert out["A"] == 7.46, "tagging must not touch the value"


# ── the promise: existing products are untouched ──────────────────────────────

def test_the_entrypoint_dry_run_changes_no_existing_product():
    """RYA-161 / spec item 5: this ADDS the asplund/gbs line sets, it does not change our
    products. Guarded by running the entrypoint and asking git what moved -- not by
    reading the code and believing it."""
    before = FEED.read_text()
    r = subprocess.run([sys.executable, "scripts/measure_reference_lineset.py",
                        "--line-set", "asplund", "--dry-run"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert FEED.read_text() == before, "the entrypoint modified the live feed"

    dirty = subprocess.run(
        ["git", "status", "--porcelain",
         "data/products", "data/results/band_products", "data/linelists"],
        cwd=ROOT, capture_output=True, text=True).stdout.strip()
    assert dirty == "", f"the dry run touched protected paths:\n{dirty}"


def test_nothing_here_writes_canonical_gf():
    src = (ROOT / "pipeline" / "reference_lineset.py").read_text()
    entry = (ROOT / "scripts" / "measure_reference_lineset.py").read_text()
    for text in (src, entry):
        assert "to_csv" not in text or "canonical_gf" not in text.split("to_csv")[0][-200:]
    assert "canonical_gf.csv is NEVER written" in src or "never written" in src.lower()
