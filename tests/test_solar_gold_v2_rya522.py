"""
tests/test_solar_gold_v2_rya522.py
RYA-522 — the solar gold reference is v2 (verdict-sourced), C is fixed, owed rows
freeze no value, and v1 is retained-immutable-but-superseded.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import data_namespace as ns  # noqa: E402


def test_current_is_v4_and_prior_versions_are_retained_immutable():
    # RYA-811 froze v4 from v3, changing ONE cell — the Fe I row's `method_scale` label,
    # 1D-NLTE -> 3D-NLTE. The VALUE 7.466 is unchanged; v3 had frozen a post-Magic-2013 3D
    # number under a pre-correction 1D label, which made `solar_scale_provenance` refuse to
    # load gold at all (the RYA-669 defect, stored). v3 is NOT deleted or edited — it is
    # retained-immutable-but-superseded, exactly as v1 was at the v2 freeze and v2 at v3.
    # The RYA-522 invariants below are version-agnostic and keep reading CURRENT; only the
    # pointer identity moved.
    _, ver = ns.read_solar_reference("CURRENT")
    assert ver == "v4"
    for prior in ("v2", "v3"):
        _, got = ns.read_solar_reference(prior)
        assert got == prior
    ns.assert_frozen_references()          # v1..v4 all match the committed manifest


def test_carbon_fixed_and_gold():
    df, _ = ns.read_solar_reference("CURRENT")
    c = df[df["element"].astype(str) == "C"].iloc[0]
    assert abs(float(c["A_X"]) - 8.491) < 1e-3     # RYA-520 fix, not 10.26
    if "confidence" in df.columns:
        assert c["confidence"] == "gold"


def test_owed_rows_freeze_no_value():
    import pandas as pd
    df, _ = ns.read_solar_reference("CURRENT")
    if "confidence" not in df.columns:
        return
    owed = df[df["confidence"] == "owed"]
    # a held-owed element (e.g. Sr, the +2.1 saturation suspect) carries NO frozen value
    assert owed["A_X"].isna().all(), "owed tier must not freeze an authoritative value"
    assert "Sr" in set(owed["element"].astype(str))


def test_v1_retained_and_superseded():
    d = ROOT / "data" / "reference" / "solar"
    assert (d / "solar_abundances_v1.csv").exists()            # historical, immutable
    assert (d / "solar_abundances_v1.SUPERSEDED.md").exists()  # marker
    prov = ns.read_provenance("v2")
    assert prov.get("supersedes") == "v1"
