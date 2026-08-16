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


def test_current_is_v5_and_prior_versions_are_retained_immutable():
    # RYA-819/831 froze v5 from v4, changing ONE cell — the Fe I row's `note`. The VALUE
    # 7.466 is unchanged AND SO IS THE SCALE LABEL; what was withdrawn is the stated ROUTE
    # onto that scale ("1D-NLTE 7.516 minus Magic-2013's 0.05"), because the measured
    # 3D-atmosphere term is -0.013..+0.018 by line sample and is linear in excitation
    # potential, so no scalar can represent it. The -0.05 landed the right answer by
    # absorbing the RYA-161 gf zero point.
    #
    # Each freeze ADDS its predecessor to the immutability check rather than replacing it,
    # so the retained set only ever grows: v1 was checked at the v2 freeze, v2 at v3, v3 at
    # v4 (RYA-811), and v4 here.
    _, ver = ns.read_solar_reference("CURRENT")
    assert ver == "v5"
    for prior in ("v2", "v3", "v4"):
        _, got = ns.read_solar_reference(prior)
        assert got == prior
    ns.assert_frozen_references()          # v1..v5 all match the committed manifest


def test_v5_moved_the_provenance_and_NOT_the_value_or_the_scale():
    """The whole point of the RYA-819 re-freeze: a note changed, nothing else.

    A regression here would mean a provenance correction had quietly moved the anchor —
    which is the failure mode gold's write-once rule (RYA-469) exists to make impossible,
    and which RYA-669 showed can pass every other gate.
    """
    import pandas as pd
    v4, _ = ns.read_solar_reference("v4")
    v5, _ = ns.read_solar_reference("v5")
    key = ["element", "ion"]
    m = v4.merge(v5, on=key, suffixes=("_v4", "_v5"))
    assert len(m) == len(v4) == len(v5), "the freeze added or dropped a row"
    for col in ("A_X", "A_X_nlte", "verdict", "confidence", "method_scale", "n_lines"):
        a, b = m[f"{col}_v4"], m[f"{col}_v5"]
        # NaN-aware: 17 of the 26 rows are `owed` and carry NO value, and NaN != NaN would
        # report every one of them as a change. Comparing as strings does not save it —
        # the two frames can render a missing cell differently ('nan' vs '<NA>') depending
        # on dtype, which is exactly what this test hit first time round.
        same = (a.isna() & b.isna()) | (a.astype(str) == b.astype(str))
        bad = m.loc[~same, ["element", "ion", f"{col}_v4", f"{col}_v5"]]
        assert same.all(), (f"v5 changed {col!r}, which a provenance re-freeze must not:\n"
                            f"{bad.to_string(index=False)}")
    fe = v5[(v5.element == "Fe") & (v5.ion.astype(str).str.upper() == "I")].iloc[0]
    assert abs(float(fe["A_X"]) - 7.466) < 1e-9
    assert "3D-NLTE" in str(fe["method_scale"])     # the SCALE was never what was wrong
    assert "ROUTE CORRECTED" in str(fe["note"])     # the note carries the withdrawal


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
