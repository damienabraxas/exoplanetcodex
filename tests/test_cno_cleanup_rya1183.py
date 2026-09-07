"""RYA-1183 — an empty bin, a missing blend partner, and two artifacts disagreeing.

RYA-1142 A6b and B3, bundled:

  A6b(1)  the atomic census is neutrals only and NOTHING recorded whether that is the
          source's shape or a filter we applied.
  A6b(2)  [O I] 6300.300 is carried, but the Ni I 6300.34 blend — the best-known
          contaminant of the most-used solar oxygen diagnostic — was retained nowhere.
  B3      summary.json said canonical_matched=0 while intake_verdict.json reported 364
          accepted, two artifacts in one directory 364 rows apart.

No gf value changes: the census gains 3 columns and one clearly-labelled component row.
"""
from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

AUDIT = ROOT / "data/audit/rya1136_cno_intake"
CENSUS = AUDIT / "atomic_source_census.csv"
VERDICT = AUDIT / "intake_verdict.json"
SUMMARY = AUDIT / "summary.json"
SCOPE = AUDIT / "atomic_census_scope_rya1183.json"

ACCEPTED = {"PRIMARY_TUPLE_MATCH", "PRIMARY_UNRESOLVED_SUM_MATCH", "PHYSICAL_TUPLE_MATCH"}


def _rows(p: Path) -> list[dict]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


# ── A6b(1): the scope, measured rather than asserted ─────────────────────────
def test_the_neutrals_only_scope_separates_source_absence_from_our_filter():
    """🔴 The two causes are NOT the same fact, and the raw table proves it: it carries
    CI, OI and FeII. So ionised CNO is absent from the SOURCE, while Fe II is present and
    dropped by US."""
    d = json.loads(SCOPE.read_text())
    codes = d["species_codes_present_in_source"]
    assert set(codes) == {"CI", "OI", "FeII"}, codes
    assert not any(c in codes for c in ("CII", "NII", "OII"))

    assert d["ionised_cno"]["cause"] == "ABSENT_FROM_SOURCE"
    assert "NOT a filter" in d["ionised_cno"]["statement"]
    assert d["fe_ii"]["cause"] == "EXCLUDED_BY_THIS_BUILDER"
    assert d["fe_ii"]["present_in_source"] == 142 and d["fe_ii"]["present_in_census"] == 0
    # and N I is neither: it does not come from this source at all
    assert d["n_i"]["cause"] == "DIFFERENT_SOURCE"
    assert "no nitrogen" in d["n_i"]["statement"]


def test_no_ionised_cno_row_was_fabricated():
    species = {r["species"] for r in _rows(CENSUS)}
    assert not species & {"C II", "N II", "O II"}, "an ionised row was invented"


# ── A6b(2): the blend component ──────────────────────────────────────────────
def test_the_ni_blend_at_6300_is_retained_as_a_component():
    """The [O I] 6300 diagnostic is not usable without it — `pipeline.cno_synthesis`
    already registers the feature as a forbidden_blend requiring joint synthesis."""
    rows = _rows(CENSUS)
    ni = [r for r in rows if r["species"] == "Ni I"]
    assert len(ni) == 1
    ni = ni[0]
    assert ni["wavelength_air_A"] == "6300.342"
    assert ni["blend_role"] == "CONTAMINANT_COMPONENT"
    assert ni["blend_partner_A"] == "6300.300"
    assert ni["use_status"] == "BLEND_COMPONENT_NOT_A_SOURCE_LINE", \
        "it must not read as an AGSS21 adopted line"
    assert "Johansson" in ni["gf_source"] and "2003" in ni["gf_source"]
    # the pairing is stated from BOTH sides
    oi = [r for r in rows if r["species"] == "O I" and r["wavelength_air_A"].startswith("6300.3")]
    assert oi and oi[0]["blend_role"] == "PRIMARY_DIAGNOSTIC_IN_A_BLEND"
    assert oi[0]["blend_partner_A"] == "6300.342"


def test_the_ni_values_come_from_our_own_store_not_from_memory():
    """Every number on that row must be the canonical_gf row it cites."""
    ni = next(r for r in _rows(CENSUS) if r["species"] == "Ni I")
    with (ROOT / "data/linelists/canonical_gf.csv").open(newline="") as fh:
        canon = next(r for r in csv.DictReader(fh) if r["line_id"] == ni["canonical_line_id"])
    assert canon["species"] == "Ni I"
    assert float(canon["wavelength_air_A"]) == float(ni["wavelength_air_A"])
    assert float(canon["log_gf"]) == float(ni["codex_loggf"])
    assert float(canon["excitation_potential_eV"]) == float(ni["codex_EP_eV"])
    assert canon["loggf_reference"] == ni["gf_source"]


def test_an_empty_blend_column_is_not_a_clean_bill():
    """Only [O I] 6300 is populated because it is the only feature the pipeline registers
    as a blend. An empty cell means 'no blend recorded', never 'checked and clean'."""
    import build_cno_intake_rya1136 as B
    assert B.blend_fields("C", 5052.0)["blend_role"] == ""
    assert B.blend_fields("O", 6363.77)["blend_role"] == ""
    assert B.blend_fields("O", 6300.30)["blend_role"] == "PRIMARY_DIAGNOSTIC_IN_A_BLEND"


# ── B3: the verdict ──────────────────────────────────────────────────────────
def test_summary_and_verdict_agree_on_the_match_count():
    """🔴 They were 364 rows apart, in one directory, and the stale one read clean."""
    s, v = json.loads(SUMMARY.read_text()), json.loads(VERDICT.read_text())
    status = v["molecular"]["join_status"]
    accepted = sum(n for k, n in status.items() if k in ACCEPTED)
    assert s["canonical_matched"] == accepted == 364
    assert s["canonical_matched"] + s["crossmatch_review"] == sum(status.values()) == 408
    assert s["verdict"] == v["verdict"]
    assert s["frozen_ready"] == v["frozen_ready_for_measurement"] is False
    assert s["reconciled_with"].startswith("intake_verdict.json")


def test_the_verdict_label_is_explained_rather_than_left_ambiguous():
    s = json.loads(SUMMARY.read_text())
    assert s["verdict"] == "INTAKE_COMPLETE_REVIEW_REQUIRED"
    note = s["verdict_label_note"]
    assert "BLOCKED_MOLECULAR_DATA" in note, "the ticket's expected string must be addressed"
    assert "census is complete" in note and "frozen_ready_for_measurement" in note


def test_the_safety_line_is_true_now_rather_than_aspirational():
    """It claimed 'no wavelength-only join admitted' while the N I join was wavelength-only.
    RYA-1143 fixed the join; RYA-1179's scope-aware guard confirms it does not flag it."""
    v = json.loads(VERDICT.read_text())
    assert "no wavelength-only join admitted" in v["safety"]
    assert "RYA-1143" in v["safety"] and "EP-aware" in v["safety"]
    disp = v["blocker_disposition"]["n_i_wavelength_only_admission"]
    assert disp["state"] == "RESOLVED" and "RYA-1179" in disp["addressed_by"]


def test_a_worked_on_blocker_is_not_recorded_as_a_gone_blocker():
    """🔴 RYA-1179/1180/1181 have all landed, and two of those blockers still block:
    making a defect VISIBLE is not removing the condition."""
    v = json.loads(VERDICT.read_text())
    disp = v["blocker_disposition"]
    assert disp["exomol_co_redistribution"]["state"] == "STILL_BLOCKING"
    assert "no Li 2015 primary" in disp["exomol_co_redistribution"]["detail"]
    assert disp["uv_scope"]["state"] == "STILL_BLOCKING"
    assert v["frozen_ready_for_measurement"] is False
    assert len(v["blocking_findings"]) >= 7


# ── no value moved ───────────────────────────────────────────────────────────
def test_the_census_gained_columns_and_a_component_but_moved_nothing():
    before = subprocess.run(["git", "show", f"HEAD:{CENSUS.relative_to(ROOT)}"],
                            cwd=ROOT, capture_output=True, text=True)
    if before.returncode != 0 or not before.stdout:
        pytest.skip("pre-change revision unavailable")
    old = list(csv.DictReader(io.StringIO(before.stdout)))
    new = _rows(CENSUS)
    if set(old[0]) == set(new[0]) and len(old) == len(new):
        pytest.skip("HEAD already contains the change")
    assert len(new) == len(old) + 1, "exactly one row added (the Ni I component)"
    assert set(new[0]) - set(old[0]) == {"blend_role", "blend_partner_A", "blend_note"}
    for k in old[0]:
        assert all(a[k] == b[k] for a, b in zip(old, new[:len(old)])), f"{k} moved"


def test_canonical_gf_was_not_touched():
    r = subprocess.run(["git", "status", "--porcelain", "data/linelists/canonical_gf.csv"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.stdout.strip() == ""
