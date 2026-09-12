"""RYA-1214 — the CNO gf adjudication, and the two refusals it is built on.

Three things have to stay true, and each of them was a live defect before this ticket:

1. **No CNO line is ever tiered LAB.** C/N/O have no primary laboratory gf in this repo
   and no table in `gf_grades.LAB_TABLES`. Opacity Project and MCHF are calculations
   (RYA-1172); tiering them LAB would repeat RYA-1005.
2. **The pull's non-physical rows stay out.** 39 C II rows carry f-values of 8e5-5e7.
   Three of them sit on a real canonical line.
3. **A component is not a blend.** NIST serves fine-structure components separately;
   `canonical_gf` holds the blend with `hfs_n_components`. Adopting on a unique
   wavelength+EP match alone took the component for the blend and biased 13 rows LOW.

The tests below are written against the ARTIFACTS and the STORE, not against the
script's internals, so a rewrite of the script that reintroduces any of the three fails
here rather than passing on its own terms.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

from pipeline.cno_gf_pedigree import pedigree_for  # noqa: E402
from pipeline.gf_grades import NIST_ACC_PCT, nist_sigma_dex  # noqa: E402

CANON = ROOT / "data/linelists/canonical_gf.csv"
AUDIT = ROOT / "data/audit/rya1214_cno_products"
ADJ = AUDIT / "cno_gf_adjudication.csv"
BY_SPECIES = AUDIT / "cno_gf_adjudication_by_species.csv"
NONPHYS = AUDIT / "nist_pull_nonphysical_rows.csv"
OI777 = AUDIT / "oi_tachiev_vs_nist.csv"
PROV = AUDIT / "cno_gf_adjudication.prov.json"

ADJUDICATED_SPECIES = ("C I", "C II", "N I", "N II", "O I", "O II")
STATUS = "nist_rya1214"


def _rows(p: Path) -> list[dict]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def canon() -> list[dict]:
    return _rows(CANON)


@pytest.fixture(scope="module")
def cno(canon) -> list[dict]:
    return [r for r in canon if r["species"].split()[0] in ("C", "N", "O")]


# ── 1. the rung ──────────────────────────────────────────────────────────────
def test_no_cno_row_is_ever_tiered_lab(cno):
    """`LAB` means a PRIMARY LABORATORY MEASUREMENT. CNO has none, in any species."""
    lab = [r for r in cno if r["gf_tier"] == "LAB"]
    assert not lab, (
        f"{len(lab)} CNO row(s) tiered LAB — e.g. "
        f"{[(r['species'], r['wavelength_air_A']) for r in lab[:5]]}. The Opacity "
        f"Project and MCHF are CALCULATIONS (RYA-1172); a critically-evaluated "
        f"theoretical gf is the NIST-C+ rung, never the laboratory one (RYA-1005).")


def test_adjudicated_rows_carry_the_compilation_tier_and_a_priced_sigma(cno):
    adj = [r for r in cno if r["adjudication_status"] == STATUS]
    assert adj, "no row carries the RYA-1214 adjudication status — the ingest is missing"
    for r in adj:
        assert r["gf_tier"] == "NIST-C+", (
            f"{r['line_id']} is {STATUS} but tiered {r['gf_tier']!r}")
        assert r["nist_grade"] in NIST_ACC_PCT, (
            f"{r['line_id']} carries grade {r['nist_grade']!r}, which is not on NIST's "
            f"ladder — an unmapped class becomes a NaN sigma wearing a tier")
        # A grade with no sigma is a hidden bar (RYA-799/713).
        assert r["gf_sigma_dex"], f"{r['line_id']} is graded but carries no gf_sigma_dex"
        assert math.isclose(float(r["gf_sigma_dex"]), nist_sigma_dex(r["nist_grade"]),
                            rel_tol=1e-4), (
            f"{r['line_id']}: sigma {r['gf_sigma_dex']} does not price grade "
            f"{r['nist_grade']} (log10(1 + pct/100))")
        assert "NIST ASD" in r["loggf_reference"], (
            f"{r['line_id']} is NIST-adjudicated but cites {r['loggf_reference']!r}")


def test_only_species_with_a_stated_pedigree_were_adjudicated(cno):
    """C III-C V / N III-N V / O III-O VI have no held authority. Unknown stays unknown."""
    for r in cno:
        if r["adjudication_status"] != STATUS:
            continue
        assert r["species"] in ADJUDICATED_SPECIES, (
            f"{r['line_id']} ({r['species']}) was adjudicated, but "
            f"cno_gf_pedigree states no authority for it — a pedigree must never be "
            f"extended from a neighbouring species (RYA-1072)")
        assert pedigree_for(r["species"]).is_laboratory is False


# ── 2. the physicality gate ──────────────────────────────────────────────────
def test_the_nonphysical_pull_rows_are_quarantined_and_none_reached_the_store():
    """The Ladenburg identity gi*fik = 1.4992e-16 * lambda^2 * gk * Aki must return a
    statistical weight. Where it returns 1e9 the row is not a datum."""
    bad = _rows(NONPHYS)
    assert bad, "the quarantine file is empty — the gate has stopped firing"
    for r in bad:
        gk = float(r["gk_implied"])
        assert not (0.5 <= gk <= 50.0), (
            f"{r['species']} {r['wavelength_A']} has a plausible implied gk {gk:.3f} "
            f"and should not be quarantined")
    # Every quarantined row is REFUSED, never repaired into a value.
    for r in bad:
        assert float(r["log_gf"]) > 1.5, (
            "a quarantined row with a physical log gf means the gate is catching "
            "something other than the delivery defect it was measured against")


#: The three canonical lines the defective rows would have landed on, with the value each
#: still holds. Named explicitly rather than re-derived, because the point of the test is
#: that THESE rows were not touched -- a re-derivation would pass by reproducing the bug.
_WOULD_HAVE_BEEN_CORRUPTED = {
    "gf_002524": (-2.005, "C II", 5039.639),
    "gf_002540": (-3.063, "C II", 5207.982),
    "gf_002609": (-2.393, "C II", 5717.984),
}


def test_the_three_lines_the_defect_would_have_hit_are_untouched(canon):
    """+8.6 to +9.7 dex, with full NIST provenance attached. The gate is what stopped it.

    ⚠️ This is deliberately NOT a blanket `log gf < 1.5` bound on the store. That bound
    is FALSE: `canonical_gf` legitimately holds log gf up to 1.783 on multi-component
    blends (Fe III 8235.397 over 15 components, O V 4930.271 over 6, O IV 7713.272 over
    3), because a blend's gf is the sum of its components and the highly-ionised species
    carry large statistical weights. A bound written from the neutral species would have
    failed on real data. The pull's own physical rows top out at log gf 1.215 and the
    largest value adopted here is 0.958.
    """
    by_id = {r["line_id"]: r for r in canon}
    for line_id, (value, species, wave) in _WOULD_HAVE_BEEN_CORRUPTED.items():
        r = by_id[line_id]
        assert r["species"] == species and abs(float(r["wavelength_air_A"]) - wave) < 1e-3
        assert math.isclose(float(r["log_gf"]), value, abs_tol=1e-6), (
            f"{line_id} ({species} {wave}) now reads {r['log_gf']} — it must still hold "
            f"{value}. The NIST row for it implies gk ~ 1e9 and is not a datum.")
        assert r["adjudication_status"] != STATUS, (
            f"{line_id} was adjudicated from a quarantined row")


def test_no_rya1214_adopted_row_exceeds_the_pulls_own_physical_ceiling(cno):
    """Every value written here came from a row that passed the Ladenburg gate, and those
    top out at log gf 1.215 across all six species. A higher one means a defective row
    got through."""
    over = [r for r in cno
            if r["adjudication_status"] == STATUS and float(r["log_gf"]) > 1.25]
    assert not over, (
        f"{len(over)} RYA-1214 row(s) above the pull's physical ceiling: "
        f"{[(r['line_id'], r['log_gf']) for r in over[:5]]}")


# ── 3. a component is not a blend ────────────────────────────────────────────
def test_admission_requires_the_multiplicities_to_agree():
    """A unique wavelength+EP match is NOT sufficient, and this is the rule that says so."""
    adj = _rows(ADJ)
    adopted = [r for r in adj if r["verdict"].startswith("ADOPT")]
    assert adopted, "nothing was adopted"
    for r in adopted:
        assert int(r["n_candidates"]) == int(r["hfs_n_components"]), (
            f"{r['line_id']} adopted with {r['n_candidates']} NIST candidate(s) against "
            f"hfs_n_components={r['hfs_n_components']} — that is a component standing in "
            f"for a blend (or the reverse)")
    refused = [r for r in adj if r["verdict"].startswith("MULTIPLICITY_MISMATCH")]
    assert refused, (
        "no row was refused on multiplicity. The mismatch class is real — 13 rows have a "
        "single NIST candidate against hfs >= 2 — so an empty refusal set means the "
        "check has been removed or defeated.")


def test_the_summed_blends_reproduce_the_store_they_replaced():
    """The admission rule is stated on multiplicity ALONE. This is its validation: where
    it fires on a blend it must land on the value the store already held, because both
    sides are summing the same components. Nothing here SELECTS by agreement."""
    adj = _rows(ADJ)
    summed = [r for r in adj if r["verdict"].startswith("ADOPT_SUM")]
    assert len(summed) >= 20, f"only {len(summed)} summed blends — too few to validate on"
    deltas = sorted(abs(float(r["delta_nist_minus_store"])) for r in summed)
    median = deltas[len(deltas) // 2]
    assert median < 0.01, (
        f"median |NIST_sum - store| over {len(summed)} summed blends is {median:.4f} dex. "
        f"Above ~0.01 the two catalogues are not summing the same components and the "
        f"multiplicity rule is not identifying lines.")


# ── the refusals that are NOT this script's to make ──────────────────────────
def test_the_oi_777_offset_is_stated_and_not_resolved():
    """Tachiev and NIST differ systematically on the dominant oxygen indicator. AGSS21
    names Tachiev for N i only and states no O i source, so the choice is Ryan's."""
    rows = _rows(OI777)
    triplet = [r for r in rows if r["wavelength_A"].startswith("777")]
    assert len(triplet) == 3, f"expected the 777 triplet, got {len(triplet)} row(s)"
    for r in triplet:
        d = float(r["tachiev_minus_nist_dex"])
        assert d < 0, "the published offset is negative on all three components"
        assert 0.014 < abs(d) < 0.018, (
            f"777 offset {d:+.4f} dex is outside the RYA-1160 measurement "
            f"(-0.0156/-0.0160/-0.0162) — restate the measurement, do not widen the test")
        assert "OPEN" in r["decision"], "the decision must stay referred, not resolved"
    prov = json.loads(PROV.read_text())
    assert prov["oi_777_decision"].startswith("OPEN")


def test_n_i_agreement_is_recorded_as_one_source_seen_twice():
    """NIST's N I values ARE Tachiev (TP T7370). Agreement there is not confirmation."""
    prov = json.loads(PROV.read_text())
    assert "T7370" in prov["n_i_caveat"]
    assert "not independent" in prov["n_i_caveat"] or "twice" in prov["n_i_caveat"]


def test_every_species_beats_its_own_randomised_null():
    """A match rate is uninterpretable without the rate the SAME matcher gets on
    randomised wavelengths."""
    for r in _rows(BY_SPECIES):
        null = float(r["randomised_null_rate"])
        rate = float(r["match_rate"])
        assert rate > 10 * max(null, 1e-6), (
            f"{r['species']}: match rate {rate:.4f} against a randomised null of "
            f"{null:.4f} — not a decisive separation")
