"""RYA-1112 — the VIS Fe uncertainty audit. READ-ONLY, and these tests pin its findings.

An audit that is not pinned rots into a claim. Each test here asserts the MECHANISM behind
a finding, not just the number it produced, so the test goes green the day the defect is
actually fixed rather than the day someone edits the artifact.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import rya1112_vis_fe_uncertainty_audit as A          # noqa: E402


@pytest.fixture(scope="module")
def doc():
    return A.audit()


def test_the_audit_reproduces_and_every_dig_in_product_has_a_named_cause(doc):
    assert A.check(doc) == []
    over = [r for r in doc["products"] if r["over_dig_in"]]
    assert len(over) == 9
    assert all(r["rca_verdict"] for r in over)
    # the skill demands a NAMED verdict, not the word "OPEN" with nothing behind it
    assert all(len(r["rca_note"]) > 60 for r in over)


# ── Part A ───────────────────────────────────────────────────────────────────
def test_type_A_is_a_standard_error_on_every_product(doc):
    assert doc["part_A_type_A"]["all_sigma_stat_are_standard_error"]


def test_F1_the_band_budget_HAS_NO_MICROTURBULENCE_TERM():
    """🔴 THE MECHANISM, not the number. `error_budget.build()` is what every published
    `sigma_syst` comes from; if it never mentions microturbulence, no product can carry it.
    This test reads the SOURCE, because a claim about a budget that only checks the budget's
    output cannot tell "absent" from "small"."""
    src = (ROOT / "pipeline" / "error_budget.py").read_text().lower()
    for token in ("vmic", "microturb", "delta_xi"):
        assert token not in src, (
            f"error_budget.py now mentions {token!r} — if the stellar-parameter terms have "
            f"been wired in, RYA-1112's finding F1 is FIXED and this test should be "
            f"replaced by one asserting the term is present and sourced")


def test_F1_the_stellar_parameter_budget_is_read_by_nothing_in_the_product_path():
    """The other half of F1: the sourced term exists and does not travel."""
    r = subprocess.run(["git", "grep", "-l", "solar_uncertainty_rya158", "--",
                        "pipeline", "scripts"], cwd=ROOT, capture_output=True, text=True)
    readers = set(r.stdout.split())
    # ⚠️ THE INSTRUMENT IS NOT THE MEASUREMENT. This audit reads the budget too, and once
    # it was tracked it showed up as a third "reader" and turned this test red. An AUDIT
    # reading a budget is not the budget reaching a product, so the auditor excludes
    # itself — explicitly and by name, never by a pattern that could hide a real reader.
    readers -= {"scripts/rya1112_vis_fe_uncertainty_audit.py"}
    # what is left is only its own two stamping scripts; nothing builds or publishes
    assert readers == {"scripts/rya1088_record_sigma_params.py",
                       "scripts/rya1089_stamp_honest_delta_xi.py"}, readers


def test_F1_the_omitted_term_is_sourced_and_large_enough_to_matter(doc):
    f = doc["findings"]["F1_stellar_parameter_terms_absent_from_published_sigma_syst"]
    assert f["sourced_delta_xi_kms"] == 0.2912          # RYA-1089, not the uncited 0.05
    assert f["reference_magnitude_FeI"] == 0.0699
    # folding it in would move 7 more products over the dig-in threshold
    assert f["products_over_dig_in_if_FeI_vmic_folded_in"] > f["products_over_dig_in_now"]


# ── Part B ───────────────────────────────────────────────────────────────────
def test_F2_the_gf_axis_is_pinned_by_the_TOKEN_not_by_the_lines():
    """The mechanism behind `gf=kurucz` on a lab-gf pool."""
    from pipeline import treatment_axes
    assert treatment_axes.LEGACY["1D-LTE"]["gf"] == "kurucz"
    assert treatment_axes.LEGACY["ENGINE-A"]["gf"] == "kurucz"
    # ...while the selector picks lines on the LAB tier
    sel = (ROOT / "scripts" / "derive_band_products.py").read_text()
    assert 'gf_tier.astype(str).str.contains("LAB"' in sel


def test_F2_a_graded_product_really_is_all_lab_gf():
    """Measured, not argued: the pool the label calls kurucz is 100% primary-laboratory."""
    import numpy as np
    import pandas as pd
    from pipeline import line_match
    f = (ROOT / "data" / "results" / "band_products"
         / "FeI_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_GRADED_1D-LTE_lines.csv")
    if not f.exists():
        pytest.skip(f"{f.name} not staged (band_products is gitignored on some checkouts)")
    d = pd.read_csv(f)
    d = d[d.in_aggregate == True] if "in_aggregate" in d else d       # noqa: E712
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    lab = cg[(cg.species == "Fe I") & cg.gf_tier.astype(str).str.contains("LAB", na=False)]
    r = line_match.match(d.wavelength_air_A.values, lab.wavelength_air_A.values,
                         want_ep=d.ep_eV.values,
                         src_ep=lab.excitation_potential_eV.values,
                         tol_A=0.005, require_ep=True)
    assert r.n_resolved == len(d), f"{r.n_resolved}/{len(d)} lab"
    # and the artifact still calls it kurucz
    # ...and the artifact still calls it kurucz. NOTE the products sibling drops the
    # treatment suffix that the lines file carries — the two names are NOT symmetric.
    prod = ROOT / "data" / "results" / "band_products" / (
        "FeI_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_GRADED_products.csv")
    assert set(pd.read_csv(prod)["gf"]) == {"kurucz"}
    # the repo contradicts itself in its own committed budget for the same product
    budget = (ROOT / "data" / "results" / "band_products"
              / "FeI_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_GRADED_budgets.txt")
    assert "every one of the 67 Fe I lines is GF-LAB" in budget.read_text()


# ── Part D, and the xi learning ──────────────────────────────────────────────
def test_F4_no_product_record_carries_a_gate_or_flag_field(doc):
    feed = json.loads((ROOT / "data" / "products" / "solar" / "Fe.json").read_text())
    fields = {k for p in feed["products"] for k in p}
    for forbidden in ("gate", "flag", "target", "exceeds_target", "honest_bar"):
        assert forbidden not in fields
    # the gate constant exists only inside AUDIT scripts — never in the product path
    r = subprocess.run(["git", "grep", "-l", "SOLAR_GATE_DEX", "--", "pipeline", "scripts"],
                       cwd=ROOT, capture_output=True, text=True)
    holders = set(r.stdout.split()) - {"scripts/rya1112_vis_fe_uncertainty_audit.py"}
    # only AUDIT scripts hold the constant — never the product path
    assert sorted(holders) == ["scripts/rya1089_stamp_honest_delta_xi.py",
                               "scripts/rya1093_xi_robustness_audit.py"], holders


def test_the_xi_applicability_follows_the_line_set_not_the_element(doc):
    """RYA-1093's learning, carried per product: 3D has no xi at all, DEEPGRADED is the
    saturated pool where xi bites, GRADED is the weaker one."""
    by = {}
    for r in doc["products"]:
        by.setdefault(r["xi_applicability"].split(" —")[0], []).append(r)
    assert "NOT APPLICABLE" in " ".join(by)
    assert doc["n_3d_products_where_xi_is_not_applicable"] == 8
    for r in doc["products"]:
        if "3D" in r["treatment"]:
            assert r["xi_applicability"].startswith("NOT APPLICABLE")


def test_the_audit_changed_nothing():
    """RYA-161: an audit reports. It must not touch a value, a product or a line.

    🔴 THIS MEASURES `audit()`, NOT THE WORKING TREE, AND THE DIFFERENCE IS THE WHOLE TEST.
    The first version ran `git status --porcelain` over data/products, data/linelists,
    data/audit and data/results/band_products. That passed on the Mac and went RED on
    Sirius — where the grids are mounted, a grid-guarded test that SKIPS locally actually
    RUNS and regenerates `data/audit/crires_co_conditioned/rya390_co_validation.json`. My
    test then read another test's side effect and reported it as "the audit modified data".

    A claim about what THIS function does has to be measured around THIS function. Hashing
    the inputs across the call does that, and is independent of whatever else the suite is
    doing to the tree.
    """
    import hashlib

    watched = sorted(
        p for d in ("data/products/solar", "data/linelists", "data/audit/uncertainty")
        for p in (ROOT / d).rglob("*") if p.is_file())
    watched += sorted((ROOT / "data" / "results" / "band_products").glob("Fe*"))
    assert len(watched) > 50, f"only {len(watched)} input files watched — too few to mean much"

    def digest():
        h = hashlib.sha256()
        for f in watched:
            h.update(f.name.encode())
            h.update(f.read_bytes())
        return h.hexdigest()

    before = digest()
    A.audit()
    assert digest() == before, "audit() wrote to one of the inputs it may only read"
