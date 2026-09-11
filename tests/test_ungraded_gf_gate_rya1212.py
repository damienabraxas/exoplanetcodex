"""RYA-1212 — the 0.17 blanket is a publication GATE, not a published value.

`UNGRADED_GF_SYSTEMATIC_DEX` is a placeholder for a gf quality nobody measured. Published,
it stops reading as a placeholder: 0.17 dex renders in the same column and on the same
error-bar forest as a cited-lab 0.041 somebody did measure, and the reader cannot tell
"wide" from "unknown". Ryan, 2026-09-10: "there is no reason at all why a published
product should ever have 0.17 stamped."

These tests pin BOTH halves, because either alone is a different and worse thing: the
budget must still ASSEMBLE with the placeholder (or the defect becomes undiagnosable — it
is the budget RYA-1211 read to discover its own `systematic:K07 x15` was a lookup miss),
and the PUBLISH path must refuse it.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from pipeline import error_budget as eb

ROOT = Path(__file__).resolve().parents[1]
BLANKET_SRC = ROOT / "data/results/rya1106/kpno_kurucz2005"
CLEAN_SRC = (ROOT / "data/results/band_products/"
             "FeI_9199_12976_kpno_solar_atlas_solar_kpno_molecfit_corrected_SYNTH_GRADED")


def _budget(**kw):
    base = dict(scatter_dex=0.012, gf_graded=False, harness_residual_dex=0.0,
                handler="SynthesisHandler")
    base.update(kw)
    return eb.build("Fe", 5500.0, 21, **base)


# ── the budget must still be assemblable ──────────────────────────────────────

def test_the_blanket_term_still_builds_so_the_defect_stays_diagnosable():
    """🔴 NOT A `raise` IN `gf_term`. RYA-1211 found its pool's `systematic:K07 x15` was a
    lookup miss BY READING THIS BUDGET. A build-time refusal would have deleted the
    evidence instead of the practice."""
    b = _budget()
    assert b.total()[1] == pytest.approx(eb.UNGRADED_GF_SYSTEMATIC_DEX)
    assert eb.carries_ungraded_gf(b.describe())


def test_the_placeholder_says_it_is_not_publishable():
    """A term that blocks publication must say so where it is read, not only in a gate."""
    term = eb.gf_term(graded=False)
    assert "NOT PUBLISHABLE" in term.source
    assert "RYA-1212" in term.source


def test_the_label_has_exactly_one_home():
    """🔴 A GATE KEYED ON A COPY OF A LABEL STOPS GATING THE DAY THE LABEL IS REWORDED.
    `publish_product` must import the constant, never spell the string (RYA-845)."""
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    assert "error_budget.UNGRADED_GF_TERM_LABEL" in src
    assert "error_budget.carries_ungraded_gf" in src
    literal = re.findall(r'["\']gf scale \(UNGRADED\)["\']', src)
    assert not literal, f"the label is written out in publish_product: {literal}"


# ── the gate ──────────────────────────────────────────────────────────────────

def _publish(tmp_path, src_dir, stem, tier, holding):
    d = tmp_path / "art"
    d.mkdir(exist_ok=True)
    products = next(p for p in Path(src_dir).parent.glob(Path(src_dir).name + "*")
                    if p.name.endswith("_products.csv")) \
        if stem is None else Path(str(src_dir) + "_products.csv")
    budget = Path(str(products)[: -len("_products.csv")] + "_budgets.txt")
    shutil.copy(products, d / "T_products.csv")
    shutil.copy(budget, d / "T_budgets.txt")
    r = subprocess.run(
        [sys.executable, "scripts/publish_product.py", "--from", str(d / "T_products.csv"),
         "--holding", holding, "--tier", tier, "--dry-run"],
        cwd=ROOT, capture_output=True, text=True)
    return r


def test_publishing_a_blanket_product_is_refused(tmp_path):
    """End-to-end through the real entry point, not a restatement of the branch."""
    r = _publish(tmp_path, BLANKET_SRC / "asplund_lines", "x", "ALL",
                 "solar_kpno_kurucz2005_corrected")
    assert r.returncode == 8, f"expected refusal, got {r.returncode}\n{r.stdout}{r.stderr}"
    assert "RYA-1212" in r.stderr and "not publishable" in r.stderr


def test_the_refusal_names_what_blocks_the_pool_and_how_to_resolve_it(tmp_path):
    """A refusal that does not say what to do is an outage. It must carry the decider's
    own verdict and the resolution paths — including drop/re-source, which is what the
    near-UV pool actually needed (one line of 57)."""
    r = _publish(tmp_path, BLANKET_SRC / "asplund_lines", "x", "ALL",
                 "solar_kpno_kurucz2005_corrected")
    for expected in ("MIXED POOL", "cited lab", "NIST grade", "per-line", "drop/re-source"):
        assert expected in r.stderr, f"refusal does not mention {expected!r}"


def test_a_clean_product_still_publishes(tmp_path):
    """🔴 THE CONTROL. A gate that refuses everything is not a gate, and this one sits in
    front of every publish in the repo."""
    r = _publish(tmp_path, CLEAN_SRC, "x", "GRADED", "solar_kpno_molecfit_corrected")
    assert r.returncode == 0, f"clean product refused:\n{r.stdout}{r.stderr}"


# ── RYA-968 wired into the production budget ──────────────────────────────────

def test_the_empirical_per_line_route_is_reachable_from_the_budget():
    """Part B. The pool `gf_graded` calls ungraded takes the per-line route instead, and
    the blanket is GONE from the budget — not merely outvoted in the quadrature."""
    b = _budget(empirical_gf_sigma_dex=0.0475,
                empirical_gf_provenance="cited x17, fallback x4 (RYA-968)")
    assert b.total()[1] == pytest.approx(0.0475, abs=5e-4)
    assert not eb.carries_ungraded_gf(b.describe())
    assert "empirical, per-line" in b.describe()


def test_the_empirical_route_is_not_gated_on_the_switch_it_replaces():
    """It must work with `gf_graded=False` — that is the pool it exists for. Requiring the
    switch to be True first would put the blanket in front of its own replacement."""
    b = _budget(gf_graded=False, empirical_gf_sigma_dex=0.05,
                empirical_gf_provenance="self-reported (n=4)")
    assert not eb.carries_ungraded_gf(b.describe())


def test_the_two_gf_routes_are_exclusive():
    with pytest.raises(ValueError, match="both describe the gf term"):
        _budget(gf_graded=True, cited_gf_sigma_dex=0.05, cited_gf_source="DenHartog2014",
                empirical_gf_sigma_dex=0.05, empirical_gf_provenance="cited")


def test_an_empirical_sigma_must_say_where_each_line_got_it():
    """RYA-968's rule: never 'assumed'. An unsourced term is how a fallback becomes a
    measurement."""
    with pytest.raises(ValueError, match="provenance"):
        _budget(empirical_gf_sigma_dex=0.05)


# ── Part C ────────────────────────────────────────────────────────────────────

def test_the_refused_set_is_exactly_the_feeds_blanket_carriers():
    """The gate's answer must match the feed's own. Both sides measured, neither asserted."""
    doc = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    carriers = {(r["holding"], r["treatment"]) for r in doc["products"]
                if abs(float(r.get("sigma_syst") or 0) - eb.UNGRADED_GF_SYSTEMATIC_DEX) < 5e-4}
    report = json.loads(
        (ROOT / "data/results/rya1212/rya1212_refused_products.json").read_text())
    refused = {(p["holding"], p["treatment"]) for p in report["products"]}
    assert refused == carriers
    assert len(carriers) == 4, f"expected the 4 Reference Grade products, got {carriers}"
    # ⚠️ keyed on the VALUE, not (holding, treatment): that pair is not unique across
    # bands and line sets, and keying on it swept a Codex Grade product into the set.
    assert {r["grade"] for r in doc["products"]
            if abs(float(r.get("sigma_syst") or 0)
                   - eb.UNGRADED_GF_SYSTEMATIC_DEX) < 5e-4} == {"Reference Grade"}


def test_every_other_published_product_is_already_clear_of_the_blanket():
    """The gate is not a cliff the whole feed falls off: 88 of 92 already pass it."""
    doc = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    over = [r["holding"] for r in doc["products"]
            if float(r.get("sigma_syst") or 0) >= eb.UNGRADED_GF_SYSTEMATIC_DEX
            and r.get("line_set") != "asplund"]
    assert not over, f"unexpected blanket carriers outside the Reference set: {over}"


def test_no_published_value_was_changed_by_this_ticket():
    """A gate changes what MAY be published, never what IS. The four refused products keep
    their values and their (wrong-by-construction) bars until they are resolved."""
    doc = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    ref = {r["holding"]: (r["A"], r["sigma_syst"], r["n_lines"]) for r in doc["products"]
           if r.get("line_set") == "asplund"}
    assert ref == {
        "solar_harps_molecfit_corrected": (7.4771, 0.17, 21),
        "solar_iag": (7.5039, 0.17, 21),
        "solar_kpno_kurucz2005_corrected": (7.4713, 0.17, 21),
        "solar_kpno_molecfit_corrected": (7.5037, 0.17, 21)}
