"""
tests/test_element_status_tracker_generator_rya654.py
====================================================
RYA-654 — the element status tracker is a GENERATED product (RYA-436 Move B, results
side), not a hand-maintained file. What these tests defend:

  * the committed CSV equals a fresh regeneration (the hand-edit detector);
  * the STATUS columns come from phase_c and nowhere else — in particular NOT from the
    frozen gold reference, which owns values, not status (RYA-654 §1), and which
    currently lags the live channel on Co and N;
  * the generator refuses to run off an UNCOMMITTED phase_c (the Mac-canary guard);
  * it loud-fails rather than emitting a blank cell when a source is missing.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import generate_element_status_tracker_rya654 as G  # noqa: E402

TRACKER = ROOT / G.TRACKER_REL


def _read_tracker():
    return pd.read_csv(TRACKER, comment="#")


def test_committed_tracker_equals_regeneration():
    """The contract that makes hand-editing a build break rather than a shortcut."""
    assert G.main(["--check"]) == 0


def test_check_mode_catches_a_hand_edit(tmp_path, monkeypatch):
    tampered = TRACKER.read_text(encoding="utf-8").replace(
        "CURATION-OWED", "PASS", 1)
    target = tmp_path / "element_status_tracker.csv"
    target.write_text(tampered, encoding="utf-8")
    monkeypatch.setattr(G, "TRACKER_PATH", target)
    assert G.main(["--check"]) == 1


def test_provenance_header_pins_the_committed_phase_c_run():
    head = [ln for ln in TRACKER.read_text(encoding="utf-8").splitlines()
            if ln.startswith("#")]
    text = "\n".join(head)
    assert "GENERATED from phase_c" in text and "do not hand-edit" in text
    prov = G.phase_c_provenance()
    assert prov["blob"] in text, "the header must pin the exact phase_c blob"
    assert prov["commit"] in text, "…and the commit that ratified it"


def test_generator_refuses_an_uncommitted_phase_c(monkeypatch):
    """The Sirius-only lever we can actually enforce: a local re-run sitting dirty in the
    working tree must not reach the ledger."""
    monkeypatch.setattr(G, "_git",
                        lambda *a: " M data/audit/cno_synthesis/solar_phase_c_verdict.json"
                        if a[0] == "status" else "deadbeef")
    with pytest.raises(G.SourceError, match="uncommitted"):
        G.phase_c_provenance()


def test_status_columns_match_phase_c_exactly():
    phase_c = G.load_phase_c()
    by_el = {r["element"]: r for r in phase_c["verdicts"]}
    df = _read_tracker()
    rows = df[df["tier"] != "diagnostic"]          # the Fe II arbiter row has no verdict
    assert len(rows) == len(by_el) == phase_c["summary"]["n_elements"]
    for _, row in rows.iterrows():
        rec = by_el[row["element"]]
        assert row["verdict"] == rec["verdict"]
        assert int(row["n_lines"]) == int(rec["n_lines"])
        if rec.get("A_measured") is None:
            assert pd.isna(row["verdict_value"])
        else:
            assert row["verdict_value"] == pytest.approx(rec["A_measured"], abs=5e-4)


def test_gold_is_not_a_source_for_status():
    """Co and N are the live proof: gold v2 still says CURATION-OWED / NLTE-OWED while
    phase_c says PASS / CURATION-OWED. The tracker must follow phase_c. If this ever
    flips, status has been re-sourced from a frozen snapshot — the exact defect that made
    the hand-maintained tracker stale."""
    from pipeline.data_namespace import read_solar_reference
    gold, _version = read_solar_reference()
    gold_verdict = {str(r["element"]).strip(): str(r["verdict"]).strip()
                    for _, r in gold.iterrows()}
    df = _read_tracker()
    tracker_verdict = {r["element"]: r["verdict"]
                       for _, r in df[df["tier"] != "diagnostic"].iterrows()}
    assert gold_verdict["Co"] == "CURATION-OWED" and tracker_verdict["Co"] == "PASS"
    assert gold_verdict["N"] == "NLTE-OWED" and tracker_verdict["N"] == "CURATION-OWED"


def test_tier_comes_from_the_rya522_ratified_tiers():
    from build_solar_reference_v2_rya522 import confidence_of
    df = _read_tracker()
    for _, row in df[df["tier"] != "diagnostic"].iterrows():
        assert row["tier"] == confidence_of(row["element"])


def test_fe_ii_arbiter_row_comes_from_the_constant_not_a_hand_typed_number():
    from config.constants import FE_IONIZATION_SYNTH_ARBITER
    df = _read_tracker()
    fe2 = df[df["tier"] == "diagnostic"]
    assert len(fe2) == 1
    row = fe2.iloc[0]
    assert row["element"] == "Fe" and row["ion"] == "II"
    assert row["verdict_value"] == pytest.approx(
        FE_IONIZATION_SYNTH_ARBITER["solar"]["fe2_synth"], abs=5e-4)
    assert pd.isna(row["verdict"]), "the arbiter row asserts no element verdict"


def test_missing_editorial_entry_loud_fails_instead_of_blanking(monkeypatch):
    editorial = G.load_editorial()
    editorial.pop("Ti")
    monkeypatch.setattr(G, "load_editorial", lambda: editorial)
    with pytest.raises(G.SourceError, match="no editorial entry for 'Ti'"):
        G.build_rows()


def test_element_missing_from_the_regime_map_loud_fails(monkeypatch):
    regime = G.load_regime()
    regime.pop("Ba")
    monkeypatch.setattr(G, "load_regime", lambda: regime)
    with pytest.raises(G.SourceError, match="not in config/physics_regime_rya400.yaml"):
        G.build_rows()


def test_generated_and_editorial_columns_do_not_overlap():
    """Each column has exactly one source. An editorial key that shadows a generated
    column would let a hand-edit win over phase_c silently."""
    editorial = G.load_editorial()
    for key, row in editorial.items():
        clash = set(row) & G.GENERATED_COLUMNS
        assert not clash, f"editorial entry {key!r} tries to set generated column(s) {clash}"


def test_runnable_as_a_script():
    r = subprocess.run([sys.executable, str(ROOT / "scripts" /
                        "generate_element_status_tracker_rya654.py"), "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "== regeneration" in r.stdout
