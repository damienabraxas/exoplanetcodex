"""
tests/test_fe_ir_band_inventory_rya762.py — RYA-762
===================================================
The 9199–13000 Å Fe band: what the inventory says, and the two claims about it
that are easy to get wrong in opposite directions.

CLAIM 1 — "Engine B cannot follow past 9199.9 Å."  FALSE.
    The ticket's product table said so. 9199.9 Å is where the GES NLTE linelist
    ENDS, not where the physics does: `atom.fe607a` is cut at 20000 Å and VALD
    carries (J, energy) for both endpoints natively. Measured: 187 of 239.

CLAIM 2 — "Fe I 12807/12808 sit beside an unmodelled Pa-beta."  ALSO FALSE.
    This one was ours, and it is retracted here. Excluding H I from the ATOMIC
    line list never meant hydrogen was unsynthesised: iSpec appends Turbospectrum's
    own `DATA/Hlinedata` on every call, and that file carries Pa-beta 12818.077
    (nblo 3, nbup 5). It is the same mistake RYA-759's v75 already retracted for
    the Balmer series, made a second time one band over.

The band's PRODUCTS are separately gated on RYA-379 — `canonical_gf.csv` stops at
9199.90 Å and `apply_to_synth_array` raises rather than defaulting — so the gate is
pinned here as a measured fact rather than left as prose.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BAND_ENGINES = ROOT / "data" / "audit" / "rya762" / "band_engines.csv"
CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"

BAND_LO_A = 9199.9


def _rows() -> list[dict]:
    with open(BAND_ENGINES, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ── the inventory ────────────────────────────────────────────────────────────
def test_the_band_is_entirely_above_the_ges_wall():
    w = [float(r["wave_air_A"]) for r in _rows()]
    assert len(w) == 239
    assert min(w) > BAND_LO_A, "a line at or below 9199.9 belongs to the existing product"
    assert max(w) < 13000.0


def test_engine_b_availability_is_what_was_measured():
    """187/239 — and the two non-available states stay DISTINCT.

    REFUSED means the atom genuinely lacks that level. NO-VALD-LEVELS means the
    accounting pool and the VALD extract disagree about the line. Collapsing them
    into one 'unavailable' bucket would hide a bookkeeping gap inside a physics
    result.
    """
    rows = _rows()
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["engine_b"]] = counts.get(r["engine_b"], 0) + 1
    assert counts == {"AVAILABLE": 187, "REFUSED": 41, "NO-VALD-LEVELS": 11}
    assert sum(counts.values()) == len(rows) == 239


def test_engine_b_availability_splits_by_band_as_measured():
    rows = _rows()
    per = {}
    for r in rows:
        key = r["band"]
        per.setdefault(key, [0, 0])
        per[key][1] += 1
        if r["engine_b"] == "AVAILABLE":
            per[key][0] += 1
    assert per["NIR"] == [107, 149]
    assert per["red-optical"] == [80, 90]


def test_the_denominator_is_named_not_a_bare_percentage():
    """187/239 is against a VALD extract cut at depth >= 0.05.

    RYA-763's 0.276 for 6910-9199 A was against GES's weak-line tail. The two are
    not comparable, and a bare '78%' invites exactly that comparison.
    """
    rows = _rows()
    assert all(r.get("predicted_depth") for r in rows), (
        "every row must carry the depth that defines the denominator")
    assert min(float(r["predicted_depth"]) for r in rows) >= 0.05


# ── the products are gated, and the gate is measurable ───────────────────────
def test_the_gf_gate_on_this_band_is_CLEARED():
    """This test used to assert the OPPOSITE, and its failure is what retired it.

    It was written to pin an absence: zero canonical gf rows above 9199.9 A, so no
    product leg could be derived until RYA-379 extended the table. Its failure message
    said that a failure would mean the gate had opened, not that something had broken.

    RYA-834 opened it — canonical_gf now runs to 12934.67 A and covers all 571 of this
    band's physical lines — so the assertion is inverted here rather than deleted. The
    band's history stays legible: the gf gate existed, and this is the commit where it
    stopped existing.
    """
    with open(CANONICAL_GF, encoding="utf-8") as fh:
        gf = [float(r["wavelength_air_A"]) for r in csv.DictReader(fh)
              if r.get("wavelength_air_A")]
    assert gf, "canonical_gf.csv did not parse"
    band = [float(r["wave_air_A"]) for r in _rows()]
    assert max(gf) >= max(band), (
        f"canonical_gf reaches {max(gf):.2f} A but the band runs to {max(band):.2f} A — "
        f"RYA-834's redward extension has regressed and the products are gf-gated again")
    assert sum(1 for x in gf if BAND_LO_A < x <= 12935.0) == 571


def test_the_REMAINING_gate_is_the_synthesis_linelist_not_the_gf_table():
    """Which gate is binding is itself a fact worth pinning.

    RYA-834 ran the products after extending canonical_gf and all three legs still
    returned n=0 — `NOT-IN-SYNTH-LINELIST`. canonical_gf says what gf to USE for a line;
    the synthesis linelist decides which lines EXIST for the engines. Only the first has
    moved, so this band's blocker is now RYA-837, not RYA-379.
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
        f"the synthesis linelist now reaches {top:.2f} A — RYA-837 has landed and this "
        f"band's products are derivable; go finish RYA-762")


# ── the H I exclusion, and the retraction ────────────────────────────────────
def test_the_h_exclusion_is_stated_as_an_engine_property():
    from pipeline.nearuv_linelist import ENGINE_REJECTED_SPECIES
    reason = ENGINE_REJECTED_SPECIES["H"]
    assert "ENGINE PROPERTY" in reason.upper()
    assert "PASCHEN" in reason.upper(), "the IR evidence must travel with the rule"
    assert "12800-12830" in reason, "the head-to-head window that established it"


def test_the_unmodelled_pabeta_caveat_is_retracted_not_merely_dropped():
    """A withdrawn claim must say it was withdrawn, or it gets re-derived."""
    from pipeline.nearuv_linelist import ENGINE_REJECTED_SPECIES
    reason = ENGINE_REJECTED_SPECIES["H"]
    assert "RETRACTED" in reason.upper()
    assert "12818" in reason, "name the line that disproves it"
    assert "Hlinedata" in reason
    assert "unmodelled Pa-beta" not in reason.replace("NOT sit beside an "
                                                      "unmodelled Pa-beta", "")


# ── against the real engine data, when mounted ───────────────────────────────
def _hlinedata() -> Path | None:
    try:
        from config.constants import codex_path
        p = codex_path('engines.ispec') / "synthesizer" / "turbospectrum" / "DATA" / "Hlinedata"
        return p if p.exists() else None
    except Exception:
        return None


@pytest.mark.skipif(_hlinedata() is None, reason="iSpec/TS DATA not mounted")
def test_hlinedata_really_does_carry_pabeta():
    """The measurement the retraction rests on."""
    lines = [l.split() for l in _hlinedata().read_text(errors="replace").splitlines()[2:]
             if l.strip()]
    waves = {}
    for f in lines:
        try:
            waves[round(float(f[0]), 1)] = (int(f[1]), int(f[2]))
        except (ValueError, IndexError):
            continue
    assert 12818.1 in waves, "Pa-beta must be present, or the retraction is wrong"
    assert waves[12818.1] == (3, 5)              # Paschen: nblo 3 -> nbup 5
    assert 10938.1 in waves and 10049.4 in waves  # Pa-gamma, Pa-delta
    series = {nblo for nblo, _ in waves.values()}
    assert {1, 2, 3} <= series, "Lyman, Balmer and Paschen are all carried"


@pytest.mark.skipif(_hlinedata() is None, reason="iSpec/TS DATA not mounted")
def test_the_fe_lines_the_caveat_named_have_pabeta_modelled_beside_them():
    band = [float(r["wave_air_A"]) for r in _rows()]
    near = [w for w in band if 12800.0 <= w <= 12815.0]
    assert near, "the 12807/12808 pair should be in the inventory"
    text = _hlinedata().read_text(errors="replace")
    assert "12818.077" in text
