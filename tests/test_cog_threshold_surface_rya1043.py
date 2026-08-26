"""RYA-1043 — the per-line COG measurement, and the surface a threshold is set FROM.

This ticket delivers a MEASUREMENT. It must not deliver a number: RYA-1041 Step 2 is
Ryan's, and `Thresholds` has no defaults precisely because a value nobody consciously
chose is indistinguishable downstream from one somebody measured (RYA-161).
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_empirical import Thresholds  # noqa: E402

AUDIT = ROOT / "data" / "audit" / "rya1041"
SURFACE = AUDIT / "rya1043_threshold_surface.json"
SCRIPT = ROOT / "scripts" / "rya1043_threshold_surface.py"
COG_SCRIPT = ROOT / "scripts" / "rya1043_perline_cog.py"


@pytest.fixture(scope="module")
def surf() -> dict:
    return json.loads(SURFACE.read_text())


def test_the_surface_sets_nothing(surf):
    """🔴 The whole point. A surface that named a winner would be setting the admission
    gate for an entire grid column under cover of reporting one."""
    assert surf["sets_no_threshold"] is True
    assert Thresholds().min_d_rew_dA is None
    for banned in ("recommended", "chosen", "adopted", "default"):
        assert banned not in json.dumps(surf).lower(), (
            f"the surface names a {banned} threshold — it must only report consequences")


def test_the_script_contains_no_candidate_value_as_a_default():
    """The surface enumerates candidate cuts; none may leak into code as a default."""
    tree = ast.parse(SCRIPT.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "add_argument":
            for kw in node.keywords:
                if kw.arg == "default" and isinstance(kw.value, ast.Constant):
                    assert not isinstance(kw.value.value, float), (
                        "a float default appeared on the surface CLI — thresholds are "
                        "declared by a person, never defaulted")


def test_the_cog_writes_its_provenance():
    """🔴 A DERIVATIVE IS MEANINGLESS WITHOUT ITS STEP. `d_rew_dA` is two-sided over
    `delta`; runs at 0.10 and 0.40 write the same columns and different numbers. The
    original script wrote the CSV and nothing else — the RYA-1006 shape, and worse here,
    because `min_d_rew_dA` is set FROM this distribution."""
    src = COG_SCRIPT.read_text()
    assert "provenance.json" in src
    assert '"delta_dex"' in src, "the provenance must record the step the derivative used"
    assert '"sets_no_threshold": True' in src


def test_the_shipped_measurement_is_flagged_as_uncertified(surf):
    """The completed 4200-6910 A run predates the provenance fix, so its delta is not
    recorded. Inferring it from the script default would be exactly the assumption this
    codebase keeps getting burned by — the surface must carry the caveat, not bury it."""
    assert surf["source_provenance"] is None, (
        "provenance now exists — drop this test and certify the delta instead")


def test_the_purity_coupling_is_reported(surf):
    """⚠️ A stricter saturation gate yields a MORE blend-dominated pool. Whoever sets the
    threshold must see that raising it buys linearity by spending purity."""
    c = surf["coupling"]
    assert c["corr_d_rew_dA_blend_fraction"] > 0.4, (
        "the positive blend coupling has vanished — re-derive before trusting the surface")
    assert c["corr_d_rew_dA_rew"] < -0.5
    assert surf["purity_cost"], "the cost of each cut must be reported, not just the count"
    strict = [r for r in surf["purity_cost"] if r["threshold"] == 0.60][0]
    assert strict["kept_blend_gt_0.5"] / strict["kept"] > 0.5, (
        "at a strict cut the majority of ADMITTED lines are mostly blend — if this no "
        "longer holds the coupling changed and the threshold advice is stale")


def test_the_unphysical_lines_are_counted(surf):
    """dREW/dA < 0 is EW falling as abundance rises — impossible for a real line, so it
    marks the blend artifacts rather than a saturation regime."""
    assert surf["n_negative_derivative"] >= 1
    assert surf["quantiles"]["p0"] < 0


def test_the_measurement_itself_is_intact():
    """617 lines with a usable derivative, spanning linear to saturated."""
    df = pd.read_csv(AUDIT / "rya1041_perline_cog_4200_6910.csv")
    ok = df[df.status == "ok"]
    d = pd.to_numeric(ok.d_rew_dA, errors="coerce").dropna()
    assert len(d) == 617
    assert d.max() > 1.0 and d.min() < 0.0
