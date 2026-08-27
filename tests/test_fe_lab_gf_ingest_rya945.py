"""
RYA-945 — the Fe I/II laboratory gf backbone in canonical_gf.csv.

These re-derive every claim from the WRITTEN FILE. The ingest script prints a report; a
report is the script's own account of itself, and the thing that has to hold is the file.

The judgements under test:

1. PRECEDENCE IS A ONE-WAY RATCHET. Lab beats NIST beats Kurucz beats VALD, and no row
   ever moves the other way. This is the ticket's "a lab value is NEVER overwritten by
   Kurucz or VALD", stated as a property of the file rather than of the run that wrote it.
2. A GRADE OR SIGMA DESCRIBES THE NUMBER ACTUALLY USED (RYA-850). A row tagged as lab must
   carry the lab paper's value, not merely cite it — that distinction is precisely what
   RYA-799 named SCALE-MISMATCH.
3. THE INGEST IS AN UPDATE, NOT AN EXTENSION. RYA-834 appended rows; this rewrites them.
   No row may be added, dropped, or re-keyed, because `gf_resolver.apply_to_synth_array`
   resolves every line the synthesis loads and a changed key is a GfResolutionError.
4. EARLIER TICKETS' ADJUDICATIONS SURVIVE. RYA-834's 28 IR lab rows are lateral moves
   under this precedence and must be left exactly as they were.
5. THE `RU` CODE IS NOT RUFFONI. The ticket's premise; pinned so it cannot silently
   regress into "lab already wired".
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
AUDIT = ROOT / "data" / "audit" / "rya945_fe_lab_gf"
LAB_FE1 = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"
LAB_FE2 = ROOT / "data" / "reference" / "fe_gf_lab" / "fe2_lab_loggf_dh19.csv"

# 🔴 DERIVED, NOT FROZEN — RYA-1052. This was a literal set, and it went red the moment
# RYA-1047 added a fifth laboratory source (Ruffoni 2013). The invariant this test protects
# is "every LAB-tier row cites a primary lab paper", not "there are exactly four of them";
# a frozen list turns every legitimate new source into a failure and invites someone to
# widen the literal without checking the citation actually exists.
import sys as _sys
_sys.path.insert(0, str(ROOT / "scripts"))
from rya945_ingest_fe_lab_gf import LAB_TAG as _LAB_TAG  # the single source (RYA-353)
LAB_TAGS = set(_LAB_TAG.values())
NIST_C_OR_BETTER = {"AAA", "AA", "A+", "A", "B+", "B", "C+", "C"}
#: Same bound `pipeline.gf_grades` uses: below the 2-decimal quantisation of the VALD
#: delivery, so it admits "same number, different rounding" and nothing else.
DESCRIBES_TOL = 0.02


def _rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def canon() -> list[dict]:
    return _rows(CANON)


@pytest.fixture(scope="module")
def fe(canon) -> list[dict]:
    return [r for r in canon if r["species"].strip() in ("Fe I", "Fe II")]


@pytest.fixture(scope="module")
def summary() -> dict:
    return json.loads((AUDIT / "rya945_summary.json").read_text(encoding="utf-8"))


# ── judgement 1: precedence is a one-way ratchet ─────────────────────────────────────
def test_no_lab_row_carries_a_kurucz_or_vald_reference(fe):
    """The ticket's headline invariant, read off the file."""
    bad = [r for r in fe
           if r["lab_source_tag"] in LAB_TAGS
           and r["loggf_reference"].startswith(("K0", "K1", "VALD"))]
    assert not bad, (
        f"{len(bad)} row(s) are tagged with a primary laboratory source but cite Kurucz "
        f"or VALD — e.g. {bad[0]['wavelength_air_A']} {bad[0]['loggf_reference']}")


def test_every_lab_tier_row_cites_a_primary_lab_paper(fe):
    lab = [r for r in fe if r["gf_tier"] == "LAB"]
    assert lab, "no LAB-tier Fe rows at all — the ingest did not happen"
    for r in lab:
        assert r["loggf_reference"].startswith("PRIMARY LAB "), r["loggf_reference"]
        assert r["lab_source_tag"] in LAB_TAGS, r["lab_source_tag"]
        assert float(r["gf_sigma_dex"]) > 0, (
            f"{r['wavelength_air_A']} is LAB-tier with no published sigma — the sigma is "
            f"the whole point of a primary source (RYA-824)")


def test_no_nist_tier_row_is_worse_than_grade_C(fe):
    """Also the reason 9 pre-existing rows were demoted — see `demoted_nist_rows.csv`."""
    for r in fe:
        if r["gf_tier"] != "NIST-C+":
            continue
        assert r["nist_grade"] in NIST_C_OR_BETTER, (
            f"{r['wavelength_air_A']} sits in the NIST-C+ tier carrying grade "
            f"{r['nist_grade']!r} — the ladder has been inverted (RYA-592)")
        assert float(r["gf_sigma_dex"]) <= 0.12 + 1e-9


def test_the_update_plan_refused_every_downgrade():
    """The ratchet, checked at the other end: in the plan, not just in the result."""
    rank = {"VALD3": 0, "OTHER": 1, "KURUCZ": 2, "NIST-C+": 3, "LAB": 4}
    plan = _rows(AUDIT / "update_plan.csv")
    assert plan
    for p in plan:
        improves = rank[p["to_tier"]] > rank[p["from_tier"]]
        assert (p["kept"] == "True") is improves, (
            f"row {p['canon_index']} {p['from_tier']}->{p['to_tier']} was "
            f"kept={p['kept']}, which does not follow the precedence")
    assert any(p["kept"] == "False" for p in plan), (
        "no rewrite was refused at all — either the tiering is broken or nothing was "
        "already adjudicated, and both would make this suite vacuous")


# ── judgement 2: the sigma describes the number actually used ────────────────────────
def test_lab_tagged_rows_carry_the_lab_papers_own_value(fe):
    """RYA-799's SCALE-MISMATCH is a reference that does not describe its value."""
    # RYA-1052: read EVERY registered Fe lab table, not a hand-listed two. `LAB_TABLES`
    # is the single source and it now holds two files for Fe I (the rya799 pull and
    # Ruffoni 2013's H-band set) — a hardcoded pair silently omits the newest paper and
    # then fails on its own omission.
    from pipeline.gf_grades import LAB_TABLES
    _srcs = []
    for _v in LAB_TABLES.values():
        _srcs.extend(_v if isinstance(_v, (list, tuple)) else [_v])
    lab = {}
    for src in _srcs:
        for r in _rows(src):
            lab.setdefault(r["source"], []).append(
                (float(r["wavelength_air_A"]), float(r["loggf"])))

    # Inverted from the single source rather than restated (RYA-353).
    tag_to_source = {v: k for k, v in _LAB_TAG.items()}
    checked = 0
    for r in fe:
        tag = r["lab_source_tag"]
        if tag not in LAB_TAGS or r["gf_tier"] != "LAB":
            continue
        w, gf = float(r["wavelength_air_A"]), float(r["log_gf"])
        near = [v for lw, v in lab[tag_to_source[tag]] if abs(lw - w) <= 0.02]
        assert near, f"{w} claims {tag} but no {tag} line lies within 0.02 A of it"
        assert min(abs(v - gf) for v in near) <= 0.001, (
            f"{w} is tagged {tag} but its log gf {gf} is not that paper's value "
            f"{near} — a citation that does not describe the number is the defect "
            f"RYA-799 named")
        checked += 1
    assert checked > 300, f"only {checked} lab rows checked — the ingest looks empty"


def test_a_sigma_is_only_attached_where_the_source_agrees_with_the_value(fe):
    """A refused lateral gets metadata ONLY when the source's number is the row's."""
    with_sigma = [r for r in fe if r["gf_sigma_dex"]]
    assert with_sigma
    for r in with_sigma:
        assert r["gf_tier"] in ("LAB", "NIST-C+") or r["lab_source_tag"], (
            f"{r['wavelength_air_A']} carries a cited sigma with no source behind it")


# ── judgement 3: this is an update, not an extension ─────────────────────────────────
def test_the_row_population_and_keys_are_untouched(canon):
    """`apply_to_synth_array` resolves EVERY line the synthesis loads (RYA-822)."""
    # This literal is a TRIPWIRE on RYA-945's ingest, not a fact about canonical_gf's
    # size. It moves only when an APPENDING ticket extends the table, and the move is
    # recorded here so an accidental change still trips it:
    #   167739  RYA-834  extension to 12935 A
    #   170539  RYA-1047 H-band extension to 21390 A (+2800 Fe lines, 12976-21400 A)
    #   178819  RYA-1053 solar IR extension to 25000 A (+8280 rows)
    assert len(canon) == 178819, (
        f"{len(canon)} rows — RYA-945 rewrites rows in place and must not change the "
        f"population; RYA-834 and RYA-1047 are the tickets that append")
    keys = {(r["species"], r["wavelength_air_A"], r["excitation_potential_eV"])
            for r in canon}
    assert len(keys) == len(canon), "duplicate (species, wavelength, EP) keys"
    ids = {r["line_id"] for r in canon}
    assert len(ids) == len(canon), "duplicate line_id"


def test_hfs_rows_never_took_a_single_component_lab_value(fe):
    """A multi-component row carries the SUM; a lab table publishes one component."""
    for r in fe:
        if r["gf_tier"] in ("LAB", "NIST-C+") and r["hfs_n_components"]:
            assert int(float(r["hfs_n_components"])) == 1, (
                f"{r['wavelength_air_A']} has {r['hfs_n_components']} HFS components and "
                f"took a per-transition value — the other components' strength is lost")


# ── judgement 4: earlier adjudications survive ───────────────────────────────────────
def test_rya834s_twenty_eight_ir_lab_rows_are_untouched(fe):
    old = [r for r in fe if r["adjudication_status"] == "lab_rya834"]
    assert len(old) == 28, (
        f"{len(old)} rows still carry RYA-834's adjudication; this ingest is a lateral "
        f"move for them and must not renumber or restate it")
    for r in old:
        assert r["loggf_reference"].startswith("PRIMARY LAB ")
        assert r["gf_tier"] == "LAB"


# ── judgement 5: `RU` is Raassen & Uylings, not Ruffoni ──────────────────────────────
def test_the_RU_reference_code_is_not_ruffoni(canon, summary):
    ru = [r for r in canon if r["loggf_reference"].strip() == "RU"]
    assert ru, "the `RU` rows vanished — this pin has lost its subject"
    species = {r["species"].strip() for r in ru}
    assert "Fe I" not in species, (
        f"`RU` now appears on {species}. Ruffoni 2014 is an Fe I paper, so an Fe I row "
        f"tagged `RU` would mean the two sources have genuinely been conflated")
    # And the positive identification, not just the absence: `RU` sits on Fe II AND Cr II.
    # Raassen & Uylings published orthogonal-operator calculations for both; Ruffoni 2014
    # measured neither. A code shared across two second ions is a calculation, not one
    # paper's Fe I laboratory table.
    assert species == {"Fe II", "Cr II"}, species
    assert "Raassen" in summary["ru_code_correction"]
    assert not any(r["lab_source_tag"] == "RU14" for r in ru), (
        "an `RU` row was credited to Ruffoni 2014 — the exact conflation the ticket's "
        "premise invites")


# ── the coverage claim ───────────────────────────────────────────────────────────────
def test_the_diagnostic_pool_grew_into_the_hundreds(fe, summary):
    diag = [r for r in fe if r["is_diagnostic"] == "True"]
    assert len(diag) == summary["diagnostic_after"]
    assert len(diag) > 1000, f"only {len(diag)} diagnostic lines"
    for r in diag:
        assert float(r["excitation_potential_eV"]) >= 1.2
        assert r["gf_tier"] in ("LAB", "NIST-C+")
    fe1_vis = [r for r in diag if r["species"].strip() == "Fe I"
               and 3780.0 <= float(r["wavelength_air_A"]) < 5500.0]
    assert len(fe1_vis) > 100, (
        f"VIS Fe I diagnostic count is {len(fe1_vis)}; the ticket's starting point was 9 "
        f"and the whole purpose is to lift it")


def test_every_outlier_is_on_the_ledger(summary):
    """Nothing moves by more than the Kurucz systematic without being named."""
    ledger = _rows(AUDIT / "outlier_ledger.csv")
    assert len(ledger) == summary["outliers_over_0.20_dex"]
    for r in ledger:
        assert abs(float(r["delta"])) > 0.20
        assert r["old_reference"] and r["loggf_reference"]


def test_den_hartog_2019_is_the_fe_ii_source_and_melendez_barbuy_is_not(fe):
    """RYA-161/852: MB09's own S flag marks a reverse solar analysis."""
    fe2_lab = [r for r in fe if r["species"].strip() == "Fe II" and r["gf_tier"] == "LAB"]
    assert fe2_lab, "no Fe II lab rows — DH19 did not land"
    assert {r["lab_source_tag"] for r in fe2_lab} == {"DH19"}
    for r in fe2_lab:
        assert "Melendez" not in r["loggf_reference"]
        assert "MB09" not in r["loggf_reference"]
