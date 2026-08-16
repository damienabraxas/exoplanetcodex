"""
tests/test_redward_gf_extension_rya834.py — RYA-834
===================================================
canonical_gf now reaches 12934.67 A, so RYA-762's near-IR products stop being
gf-gated. What these pin is not the extension itself but the three judgements
inside it, each of which could have gone the other way and looked fine.

1. THE REFEREE IS PRIMARY LAB, NOT NIST — and that differs from RYA-822 deliberately.
   Blueward, NIST ASD was the only graded source reaching the band. Redward it is not
   independent: matched on wavelength AND EP against this band's own linelist it agrees
   to a MEDIAN 0.0003 dex over 34 lines. That is the same number arriving twice
   (RYA-760: FMW is a NIST compilation, VALD copies it), so NIST is recorded and never
   adopted.

2. THE 11316-12935 A ABSENCE IS CHECKED, NOT ASSUMED (RYA-833). It is also not where
   the ticket said: the ticket put the gap at 10864 A, but our vendored Ruffoni pull
   reaches 11316.06 A and Den Hartog reaches 11013.24 A.

3. AN UNADJUDICATED LINE STAYS ON THE KURUCZ FLOOR AND SAYS SO. 70 Fe I lines above the
   last lab measurement keep `single_source`, even though NIST has values for 16 of them
   — because adopting a compilation echo would look like adjudication while adding
   nothing.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
AUDIT = ROOT / "data" / "audit" / "rya834_redward_gf"
ADJ = AUDIT / "redward_fe1_gf_adjudication.csv"
SUMMARY = AUDIT / "redward_gf_summary.json"
LAB = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"

BAND_LO_A = 9199.90
BAND_HI_A = 12935.0
LAST_LAB_A = 11316.064


def _canon() -> list[dict]:
    with open(CANON, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _adj() -> list[dict]:
    with open(ADJ, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _summary() -> dict:
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def _true(v) -> bool:
    return str(v).strip().lower() == "true"


# ── the extension itself ─────────────────────────────────────────────────────
def test_canonical_gf_reaches_the_top_of_the_762_band():
    w = [float(r["wavelength_air_A"]) for r in _canon() if r.get("wavelength_air_A")]
    assert max(w) > 12934.0, "the band tops out at 12934.67 A"
    assert min(w) < 3001.0, "RYA-822's blue edge must not regress"


def test_the_band_is_populated_not_merely_spanned():
    """A span is not coverage — a single row at 12934 would pass a max() check."""
    w = [float(r["wavelength_air_A"]) for r in _canon() if r.get("wavelength_air_A")]
    assert sum(1 for x in w if BAND_LO_A < x <= BAND_HI_A) == 571


def test_the_resolver_now_serves_every_line_in_the_band():
    """THE UNBLOCK. Before this, `apply_to_synth_array` raised above 9199.90 A, so
    none of RYA-762's three product legs could be derived at all."""
    import pandas as pd
    from pipeline.gf_resolver import apply_to_linelist_df

    a = pd.read_csv(ADJ)
    d = pd.DataFrame({
        "element": a.element_raw.str.split().str[0],
        "ion": a.ion.astype(int),
        "wavelength_air_A": a.wavelength_air_A.astype(float),
        "excitation_potential_eV": a.excitation_potential_eV.astype(float),
        "log_gf": a.gf_linelist_vald.astype(float),
    })
    out = apply_to_linelist_df(d)          # raises GfResolutionError if any line misses
    assert len(out) == len(d) == 571


# ── judgement 1: NIST is recorded, never adopted ─────────────────────────────
def test_no_band_row_adopts_nist_as_its_gf():
    band = [r for r in _canon()
            if r.get("wavelength_air_A")
            and BAND_LO_A < float(r["wavelength_air_A"]) <= BAND_HI_A]
    assert band
    for r in band:
        assert "NIST" not in (r.get("loggf_reference") or ""), (
            f"{r['wavelength_air_A']} adopted a NIST value; in this band NIST agrees "
            f"with the linelist to ~0.0003 dex and cannot referee it")
        assert r.get("adjudication_status") in ("lab_rya834", "single_source")


def test_the_summary_records_why_nist_was_not_adopted():
    s = _summary()
    assert s["nist_adopted"] is False
    assert "0.0003" in s["nist_note"], "the measurement, not just the assertion"
    assert "RYA-760" in s["nist_note"]


def test_adjudicated_rows_cite_a_primary_laboratory_source():
    band = [r for r in _canon()
            if r.get("wavelength_air_A")
            and BAND_LO_A < float(r["wavelength_air_A"]) <= BAND_HI_A
            and r.get("adjudication_status") == "lab_rya834"]
    assert len(band) == 28
    ok = {"Ruffoni2014", "DenHartog2014", "Belmonte2017"}
    for r in band:
        ref = r["loggf_reference"]
        assert ref.startswith("PRIMARY LAB ")
        assert ref.split("PRIMARY LAB ")[1].strip() in ok


# ── judgement 2: the absence is checked, and it is not where the ticket said ──
def test_the_lab_sources_reach_further_than_the_ticket_stated():
    """The ticket puts Ruffoni's top at 10864 A. Our vendored pull reaches 11316.06."""
    with open(LAB, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    top = {}
    for r in rows:
        s = r["source"]
        top[s] = max(top.get(s, 0.0), float(r["wavelength_air_A"]))
    assert top["Ruffoni2014"] == pytest.approx(11316.064, abs=0.01)
    assert top["DenHartog2014"] == pytest.approx(11013.237, abs=0.01)
    assert top["Ruffoni2014"] > 10864.0, (
        "the ticket's stated 3526-10864 A span is short; the gap it asked about "
        "opens at 11316.06 A, not 10864")


def test_the_absence_above_the_last_lab_line_is_recorded_as_measured():
    s = _summary()
    assert s["lab_span_A"][1] == pytest.approx(LAST_LAB_A, abs=0.01)
    assert s["lines_above_last_lab_line"] == 70
    assert s["nist_rows_above_last_lab_line"] == 16
    assert "11316.06" in s["ticket_span_correction"]


def test_every_line_above_the_last_lab_measurement_stays_single_source():
    """The Kurucz floor, stated. 16 of these have a NIST value; none may adopt it."""
    band = [r for r in _canon()
            if r.get("wavelength_air_A")
            and float(r["wavelength_air_A"]) > LAST_LAB_A
            and float(r["wavelength_air_A"]) <= BAND_HI_A]
    assert band
    assert all(r["adjudication_status"] == "single_source" for r in band)


# ── judgement 3: the matcher is controlled ───────────────────────────────────
def test_the_match_rate_beats_a_randomised_null():
    """Without this the rate is uninterpretable — dense bands match by accident."""
    s = _summary()
    assert s["adjudicated_on_lab"] == 28
    assert s["ambiguous"] == 0, "an ambiguous line must be left unadjudicated, not guessed"
    assert s["randomised_null_rate"] == 0.0
    assert s["match_rate"] > s["randomised_null_rate"]


def test_matching_used_both_wavelength_and_excitation_potential():
    """Wavelength alone manufactured a 2.85 dex false discrepancy in RYA-780."""
    s = _summary()
    assert s["tolerances"]["wavelength_A"] > 0
    assert s["tolerances"]["excitation_eV"] > 0


def test_the_adopted_values_carry_a_real_per_line_sigma():
    adj = [r for r in _adj() if _true(r["adjudicated"])]
    assert len(adj) == 28
    sig = [float(r["lab_sigma_dex"]) for r in adj]
    assert all(0.0 < x <= 0.5 for x in sig)
    sig.sort()
    assert sig[len(sig) // 2] < 0.20, (
        "the point of a primary source is a sigma tighter than the 0.20 blanket")


def test_the_synthesis_linelist_is_the_SECOND_wall_and_it_has_not_moved():
    """Step 4's blocker, pinned — extending canonical_gf was necessary, not sufficient.

    `canonical_gf` says what gf to USE for a line. The synthesis linelist decides which
    lines EXIST for the engines at all. Both stopped at 9199.9 A and only one has moved,
    so all three product legs return `NOT-IN-SYNTH-LINELIST` and n=0.

    Asserted on the COMMITTED GES linelist rather than on a live `_load_synth_resources()`
    call, so it runs without iSpec. When this test starts failing, the synthesis input has
    been extended redward and RYA-834 step 4 became runnable — that is the signal to go
    derive the products, not to delete the test.
    """
    from config.constants import codex_path
    ges = codex_path('engines.ges_nlte_linelist')
    if not ges.exists():
        pytest.skip("GES linelist not mounted (engine volume offline)")
    top = 0.0
    for raw in ges.read_text(errors="replace").splitlines():
        if raw.startswith("'") or not raw.startswith(" "):
            continue
        try:
            top = max(top, float(raw.split()[0]))
        except (ValueError, IndexError):
            continue
    assert top < 9200.0, (
        f"the synthesis linelist now reaches {top:.2f} A — RYA-834 step 4 is unblocked; "
        f"derive the 762 band products (two of them: red-optical 9203-10000 and NIR "
        f"10000-12935, per RYA-712)")


def test_canonical_gf_reaches_past_the_synthesis_linelist_on_purpose():
    """The two spans are ALLOWED to differ, and the direction matters.

    canonical_gf leading is harmless — a gf nobody asks for is inert. The reverse would
    be the failure: a synthesis line with no canonical gf raises GfResolutionError.
    """
    w = [float(r["wavelength_air_A"]) for r in _canon() if r.get("wavelength_air_A")]
    assert max(w) > 12934.0 > 9200.0


def test_adjudication_moved_gf_for_exactly_the_adopted_lines():
    """28 changed, and only those — an extension must not perturb untouched lines."""
    changed = [r for r in _adj()
               if _true(r["adjudicated"])
               and abs(float(r["lab_log_gf"]) - float(r["gf_linelist_vald"])) > 1e-9]
    assert len(changed) >= 25, "near-identical lab/VALD values are possible but not 28 of 28"
    unadj = [r for r in _adj() if not _true(r["adjudicated"])]
    assert all(r["lab_log_gf"] in ("", "nan") or r["lab_log_gf"] is None for r in unadj)
