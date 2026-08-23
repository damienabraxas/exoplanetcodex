"""The artifact stem carries the OBSERVED-SPECTRUM CONDITIONING axis — RYA-1006.

🔴 THIS TEST FILE EXISTS BECAUSE THE FAILURE ALREADY HAPPENED, not because it might.
On 2026-08-23 two RYA-1000 `--local-renorm` runs wrote
`FeI_4200_6910_*_SYNTH_DEEPGRADED_*` — the exact filenames
`pipeline.anchor_pools.ANCHORS['rya984_graded_163']` names — and the anchor's Kitt Peak
value moved 7.417 -> 7.337 under its own name, with a BYTE-IDENTICAL provenance file.

The invariant, not the example (RYA-870): every flag that changes what was measured must
change the stem, and the guard must refuse a part it cannot verify.
"""
from __future__ import annotations

import sys
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from derive_band_products import (conditioning_id, conditioning_tag,   # noqa: E402
                                  _conditioning_note, _selector_tag, synthesis_stem)
from pipeline.anchor_pools import ANCHORS, AnchorPool, _assert_conditioning  # noqa: E402


def _args(**kw) -> Namespace:
    base = dict(local_renorm=False, degrade_to_R=None,
                lines_deep_graded=False, lines_from_ew=None, lines_tier="all",
                element="Fe", ion="I", lo=4200, hi=6910,
                instrument="kpno_solar_atlas", holding=None)
    base.update(kw)
    return Namespace(**base)


# ------------------------------------------------------------------ the stem itself

def test_an_unconditioned_product_keeps_the_name_it_has_always_had():
    """Every committed artifact predates this axis. Renaming them all would break every
    anchor, every register row and every path in the ledgers to fix a defect they do not
    have — the default MUST stay empty (RYA-984's own rule)."""
    assert conditioning_tag(_args()) == ""
    assert _conditioning_note(_args()) == ""


def test_local_renorm_cannot_share_a_filename_with_the_unconditioned_product():
    """🔴 THE EXACT COLLISION THAT DESTROYED THE ANCHOR."""
    anchor = _args(lines_deep_graded=True)
    renorm = _args(lines_deep_graded=True, local_renorm=True)
    assert _selector_tag(anchor) == _selector_tag(renorm) == "_DEEPGRADED", (
        "the premise: these two runs agree on every axis the stem used to carry")
    assert conditioning_tag(anchor) != conditioning_tag(renorm)


def test_degrade_to_R_is_the_same_axis_and_is_also_in_the_stem():
    """RYA-995's resolution-degradation flag changes the observed spectrum too. It was
    open at the same time and by the same mechanism; closing only the one that happened
    to fire would leave the class open."""
    assert conditioning_tag(_args(degrade_to_R=60000)) == "_R60000"
    assert conditioning_tag(_args(degrade_to_R=60000.4)) == "_R60000"


def test_two_conditionings_compose_deterministically():
    both = conditioning_tag(_args(local_renorm=True, degrade_to_R=90000))
    assert both == "_LOCALRENORM_R90000"
    assert both == conditioning_tag(_args(degrade_to_R=90000, local_renorm=True))


def test_every_distinct_conditioning_gets_a_distinct_tag():
    combos = [_args(), _args(local_renorm=True), _args(degrade_to_R=60000),
              _args(degrade_to_R=90000), _args(local_renorm=True, degrade_to_R=60000)]
    tags = [conditioning_tag(c) for c in combos]
    assert len(set(tags)) == len(tags), tags


def test_a_conditioned_product_says_so_in_its_provenance():
    """The overwrite was undetectable because the provenance was byte-identical. A reader
    holding only the file must be able to see what was done to the spectrum."""
    note = _conditioning_note(_args(local_renorm=True))
    assert "DIAGNOSTIC" in note and "95th-percentile" in note and "RYA-1000" in note
    assert "RYA-995" in _conditioning_note(_args(degrade_to_R=60000))


# ------------------------------------------------- the machine-readable row value

def test_native_is_spelled_out_not_left_blank():
    """A consumer must distinguish 'not conditioned' from 'predates the column' (RYA-833).
    A blank cannot carry that distinction; `anchor_pools` depends on it."""
    assert conditioning_id(_args()) == "native"
    assert conditioning_id(_args(local_renorm=True)) == "localrenorm"
    assert conditioning_id(_args(local_renorm=True, degrade_to_R=60000)) == "localrenorm+r60000"


def test_the_row_value_tracks_the_stem():
    for c in (_args(), _args(local_renorm=True), _args(degrade_to_R=75000),
              _args(local_renorm=True, degrade_to_R=75000)):
        tag, rid = conditioning_tag(c), conditioning_id(c)
        assert (tag == "") == (rid == "native")
        if tag:
            assert rid == tag.lstrip("_").replace("_", "+").lower()


# ------------------------------------------------------------- the anchor's refusal

_POOL = AnchorPool(name="t", parts=("p.csv",), species="Fe I", note="", conditioning="native")


def _rows(value):
    d = pd.DataFrame({"wavelength_air_A": [5000.0, 5001.0], "abundance": [7.4, 7.5]})
    if value is not _MISSING:
        d["observed_conditioning"] = value
    return d


_MISSING = object()


def test_the_anchor_accepts_the_conditioning_it_declares():
    _assert_conditioning(_rows("native"), part="p.csv", pool=_POOL)


def test_the_anchor_refuses_a_differently_conditioned_product():
    """🔴 THE GUARD FOR THE ACTUAL INCIDENT: same filename, different measurement."""
    with pytest.raises(SystemExit, match="DIFFERENT measurement"):
        _assert_conditioning(_rows("localrenorm"), part="p.csv", pool=_POOL)


def test_a_product_that_cannot_say_is_REFUSED_not_assumed_native():
    """⚠️ The corrupted files were silent too. Reading silence as 'native' would re-admit
    exactly the products this guard was written for (RYA-833)."""
    with pytest.raises(SystemExit, match="UNVERIFIABLE"):
        _assert_conditioning(_rows(_MISSING), part="p.csv", pool=_POOL)


def test_a_blank_cell_is_UNKNOWN_not_native():
    with pytest.raises(SystemExit, match="UNKNOWN"):
        _assert_conditioning(_rows([None, "native"]), part="p.csv", pool=_POOL)


def test_a_mixed_pool_is_refused_rather_than_averaged():
    with pytest.raises(SystemExit, match="DIFFERENT measurement"):
        _assert_conditioning(_rows(["native", "localrenorm"]), part="p.csv", pool=_POOL)


def test_the_shipped_anchor_declares_its_conditioning():
    """Every anchor row must state it; a default nobody looked at is not a declaration."""
    for name, pool in ANCHORS.items():
        assert pool.conditioning, name
        assert pool.conditioning == pool.conditioning.lower()


# ------------------------------------- the SHIPPED stem, on the exact runs that collided

#: The two RYA-1000 commands, verbatim from `ps` at 2026-08-23 07:23, and the RYA-984
#: anchor run they overwrote. Driving `synthesis_stem` rather than a re-typed f-string:
#: a reconstructed stem agrees with itself while the route writes something else.
_ANCHOR_KP = _args(lines_deep_graded=True)
_ANCHOR_HARPS = _args(lines_deep_graded=True, instrument="harps",
                      holding="solar_harps_molecfit_corrected")
_RENORM_KP = _args(lines_deep_graded=True, local_renorm=True)
_RENORM_HARPS = _args(lines_deep_graded=True, local_renorm=True, instrument="harps",
                      holding="solar_harps_molecfit_corrected")


def test_the_anchor_run_still_writes_the_anchors_filename():
    """The committed anchor name must survive the fix, or the fix breaks what it protects."""
    assert (synthesis_stem(_ANCHOR_KP)
            == "FeI_4200_6910_kpno_solar_atlas_SYNTH_DEEPGRADED")
    assert (synthesis_stem(_ANCHOR_HARPS)
            == "FeI_4200_6910_harps_solar_harps_molecfit_corrected_SYNTH_DEEPGRADED")


def test_the_two_runs_that_collided_no_longer_collide():
    """🔴 THE REGRESSION. Before RYA-1006 both sides of each pair returned the same stem
    and the second run silently overwrote the anchor."""
    for anchor, renorm in ((_ANCHOR_KP, _RENORM_KP), (_ANCHOR_HARPS, _RENORM_HARPS)):
        assert synthesis_stem(anchor) != synthesis_stem(renorm)
        assert synthesis_stem(renorm) == synthesis_stem(anchor) + "_LOCALRENORM"


def test_the_stem_separates_every_axis_pairwise():
    """One axis moving must move the stem, whichever axis it is."""
    base = _args()
    variants = {
        "holding": _args(holding="solar_harps_molecfit_corrected"),
        "instrument": _args(instrument="harps"),
        "band": _args(hi=7000),
        "ion": _args(ion="II"),
        "selection": _args(lines_deep_graded=True),
        "conditioning": _args(local_renorm=True),
        "resolution": _args(degrade_to_R=60000),
    }
    stems = {"base": synthesis_stem(base)}
    stems.update({k: synthesis_stem(v) for k, v in variants.items()})
    assert len(set(stems.values())) == len(stems), stems


def test_the_FROMEW_products_the_live_runs_target_also_separate():
    """The 08:07 runs were `--lines-from-ew ... --local-renorm`, aimed at the RYA-986
    pools RYA-992's arm scale was derived on — the same collision one selector over."""
    plain = _args(lines_from_ew="x.csv")
    renorm = _args(lines_from_ew="x.csv", local_renorm=True)
    assert synthesis_stem(plain).endswith("_SYNTH_FROMEW")
    assert synthesis_stem(renorm).endswith("_SYNTH_FROMEW_LOCALRENORM")
