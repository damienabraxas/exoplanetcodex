"""
tests/test_term_aware_level_keying_rya823.py — RYA-823
======================================================
A super-level has no J, and must never be matched as though it had one.

RYA-776 already guarded the total case: if NOTHING in a level table resolves, the
key is wrong rather than the physics absent, so it reports REACH-UNKNOWN instead of
UNCOVERED. That guard is correct and stays. What it could not catch is the PARTIAL
case — `J = (g-1)/2` on a super-level yields a number, and a number occasionally
matches something. Cr I VIS resolved 4 of 5353 lines that way (0.1%), which is
enough to make `n_either == 0` false, defeat the guard, and let 0.1% be published as
a measured reach.

So the fix is two-sided and both sides are pinned here:
  * J is NaN for a super-level, so the miss is CLEAN and the guard can see it, and
  * the term label is added as a second key, so super-levels resolve for real.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.model_atom import (  # noqa: E402
    read_gerber_atom, resolvable_j, resolve_level)
from scripts.rya763_level_mapping import read_gerber_atom as reader_763  # noqa: E402

MIXED_ATOM = """Ti
    4.95   48.000
    4     9     3     0
*  EC          G       LABEL                   ION
       0.000      9.0  '    Level    1 =       a5D4 '   1
     170.150      7.0  '    Level    2 =       a5D3 '   1
    8090.208     25.0  '    Level    3 =        a5D '   1
   54580.617      6.0  '    Level    4 =        a6S '   2
"""


TERM_ATOM = """Cr
    5.74   52.000
    3     9     2     0
*  EC          G       LABEL                   ION
       0.000      7.0  '    Level    1 =        a7S '   1
    8090.208     25.0  '    Level    2 =        a5D '   1
   24079.555     33.0  '    Level    3 =        a3H '   1
"""


@pytest.fixture
def mixed(tmp_path) -> Path:
    p = tmp_path / "atom.ti4"
    p.write_text(MIXED_ATOM)
    return p


@pytest.fixture
def term_resolved(tmp_path) -> Path:
    p = tmp_path / "atom.cr3"
    p.write_text(TERM_ATOM)
    return p


# ── J is withheld ONLY where the fabricated one is all there is ──────────────
def test_a_term_resolved_atom_offers_no_J_at_all(term_resolved):
    """Cr's case: every level is a term, so (g-1)/2 cannot be right for any of them."""
    lv = read_gerber_atom(term_resolved)
    assert resolvable_j(lv).notna().sum() == 0


def test_a_MIXED_atom_KEEPS_its_super_level_J_deliberately(mixed):
    """⚠️ Restraint, not oversight.

    Withholding J here costs measured reach (Mn II VIS 171 -> 93, Ti I VIS
    1534 -> 1506) and those lost matches cannot be shown to be WRONG: a low-g
    super-level's (g-1)/2 can coincide with a real J, the energy still had to agree
    within 1 meV, and the coefficient returned is the term's either way. Settling it
    needs RYA-818's raw-label route. Until then this change stays MONOTONE — the
    label union only ADDS.
    """
    lv = read_gerber_atom(mixed)
    j = resolvable_j(lv)
    assert float(j[lv["term"] == "a5D4"].iloc[0]) == 4.0     # fine structure kept
    assert float(j[lv["term"] == "a5D"].iloc[0]) == 12.0     # super-level ALSO kept


def test_the_fabricated_J_would_have_been_12(mixed):
    """The exact number the old key produced — pinned so the regression is legible."""
    lv = read_gerber_atom(mixed)
    g = float(lv.loc[lv["term"] == "a5D", "g"].iloc[0])
    assert g == 25.0 and (g - 1) / 2 == 12.0


def test_the_763_reader_follows_the_same_rule(term_resolved, mixed):
    """RYA-776 imports THIS reader for the whole Engine-B deck."""
    assert reader_763(term_resolved)["J"].notna().sum() == 0
    lv = reader_763(mixed)
    assert float(lv.loc[lv["term"] == "a5D4", "J"].iloc[0]) == 4.0


def test_the_763_reader_stopped_calling_everything_Fe(mixed):
    lv = reader_763(mixed)
    assert set(lv["species"]) == {"Ti 1", "Ti 2"}
    assert not any(str(s).startswith("Fe") for s in lv["species"])


# ── the union: each key addresses what it can ────────────────────────────────
def test_j_key_resolves_a_fine_structure_level(mixed):
    lv = read_gerber_atom(mixed)
    e = float(lv.loc[lv["term"] == "a5D4", "energy_eV"].iloc[0])
    assert resolve_level(lv, energy_eV=e, j=4.0)[0] == "UNIQUE"


def test_j_key_alone_cannot_reach_a_super_level_in_a_term_resolved_atom(term_resolved):
    """The failure this ticket fixes, in the atom where it actually bites."""
    lv = read_gerber_atom(term_resolved)
    e = float(lv.loc[lv["term"] == "a5D", "energy_eV"].iloc[0])
    for j in (0.0, 1.0, 2.0, 3.0, 4.0, 12.0):
        assert resolve_level(lv, energy_eV=e, j=j)[0] != "UNIQUE"
    # ...and the label route reaches it
    assert resolve_level(lv, term="a5D")[0] == "UNIQUE"


def test_label_key_reaches_the_super_level(mixed):
    lv = read_gerber_atom(mixed)
    verdict, idx, _ = resolve_level(lv, term="a5D")
    assert verdict == "UNIQUE" and idx == 3


def test_label_key_does_not_hijack_a_fine_structure_level(mixed):
    """`a5D4` is J-resolved, so the label route must not claim it."""
    lv = read_gerber_atom(mixed)
    assert resolve_level(lv, term="a5D4")[0] == "ABSENT"


def test_union_uses_whichever_key_the_level_supports(mixed):
    lv = read_gerber_atom(mixed)
    e_fs = float(lv.loc[lv["term"] == "a5D4", "energy_eV"].iloc[0])
    assert resolve_level(lv, energy_eV=e_fs, j=4.0, term="a5D")[0] == "AMBIGUOUS", (
        "two keys naming DIFFERENT levels is a contradiction, not an answer")
    assert resolve_level(lv, energy_eV=e_fs, j=4.0, term="")[0] == "UNIQUE"
    assert resolve_level(lv, energy_eV=float("nan"), j=float("nan"),
                         term="a5D")[0] == "UNIQUE"


def test_no_key_is_distinct_from_absent(mixed):
    lv = read_gerber_atom(mixed)
    assert resolve_level(lv)[0] == "NO-KEY"
    assert resolve_level(lv, term="none")[0] == "NO-KEY"
    assert resolve_level(lv, energy_eV=999.0, j=1.0)[0] == "ABSENT"


def test_a_J_query_against_a_term_resolved_table_is_NO_KEY_not_ABSENT(term_resolved):
    """The distinction the whole guard rests on.

    A term-resolved atom has NaN for every J, so a (J, energy) query cannot match.
    Reporting that as ABSENT reads as a measured absence — 'the model does not carry
    this line' — when the truth is that the key does not apply to this table. Cr is
    the live case: reported as ABSENT it yielded `SERVED, reach 3 of 5353` for a
    species RYA-818 measures at 88.3%.
    """
    lv = read_gerber_atom(term_resolved)
    assert resolve_level(lv, energy_eV=1.0, j=2.0)[0] == "NO-KEY"
    assert resolve_level(lv, term="a5D")[0] == "UNIQUE"   # the label still works


def test_a_mixed_table_does_offer_J_so_the_same_query_is_a_real_absence(mixed):
    """The other side of the distinction — NO-KEY must not swallow real absences."""
    lv = read_gerber_atom(mixed)
    assert resolve_level(lv, energy_eV=999.0, j=2.0)[0] == "ABSENT"


# ── Engine A must be untouched ───────────────────────────────────────────────
def test_engine_a_frames_without_kind_still_resolve_by_their_own_J():
    """`label_*.txt` decks carry no `g`/`kind`; the union must not demand them."""
    lab = pd.DataFrame(dict(index=[1, 2], species=["Ti 1", "Ti 1"],
                            term=["a5D4", "a5D3"], J=[4.0, 3.0],
                            energy_eV=[0.0, 0.0211], ion=[1, 1]))
    assert resolve_level(lab, energy_eV=0.0, j=4.0)[0] == "UNIQUE"
    assert resolve_level(lab, energy_eV=0.0211, j=3.0)[0] == "UNIQUE"
    assert resolve_level(lab, energy_eV=0.0, j=1.0)[0] == "ABSENT"
    # and a term arrives harmlessly: without `kind` there is no super-level route
    assert resolve_level(lab, energy_eV=0.0, j=4.0, term="a5D4")[0] == "UNIQUE"


# ── against the real atoms, when mounted ─────────────────────────────────────
def _atom(name: str):
    try:
        from config.constants import codex_path
        p = codex_path('grids.gerber_ts') / name
        return p if p.exists() else None
    except Exception:
        return None


needs_grids = pytest.mark.skipif(_atom("atom.cr374") is None,
                                 reason="grid volume not mounted")


@needs_grids
def test_cr_offers_no_J_at_all_so_the_776_guard_can_fire():
    """Cr is term-resolved throughout: with J now NaN, NOTHING resolves by J, which
    is what lets RYA-776's `n_either == 0` guard report REACH-UNKNOWN instead of a
    0.1% reach dressed as a measurement."""
    lv = reader_763(_atom("atom.cr374"))
    assert lv["J"].notna().sum() == 0


@needs_grids
def test_the_mixed_atoms_keep_their_fine_structure_J():
    """Ti/Mn keep every J-resolved level — the fix must not cost them reach."""
    for name, at_least in (("atom.ti503b", 200), ("atom.mn281kbc", 60)):
        lv = reader_763(_atom(name))
        assert lv["J"].notna().sum() >= at_least, name
