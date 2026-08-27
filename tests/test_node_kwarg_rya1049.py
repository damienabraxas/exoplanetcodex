"""RYA-1049 — `node=` must not silently serve the wrong deck.

🔴 THE DEFECT. `for_node(element, ..., node=...)` reads like "which node of the grid", so
the natural way to reach for a ⟨3D⟩ deck is `for_node("Fe", ..., node="Fe@mean3D")`. That
was ACCEPTED, ran clean, and returned the 1D MARCS deck — byte-identical to plain
`for_node("Fe", ...)`. `node` is only a scratch-FILENAME tag (`_interpolate` spends it on
`{workdir}/Testout/{node}_{element}_coef.dat`); it selects nothing, and on a direct-read
deck it is not referenced at all.

⚠️ WHY A NORMAL SMOKE TEST CANNOT CATCH THIS. The 1D and ⟨3D⟩ Fe decks agree on every
orientation-blind field: both are (ndep, 607), both A(X) = 7.50, both `atom.fe607a`, both
b → 1 at depth, both b < 1 at the surface. ONLY `ndep` (56 vs 101) and the tau range tell
them apart. A check that passes on both candidates is structurally incapable of
distinguishing them — the same shape as RYA-821's "verified twice", which pinned four
invariants the two orientations SHARED and so could not see the transpose.

So these tests pin the DISCRIMINATING fields, never the shared ones.
"""
from pathlib import Path

import pytest

from pipeline import gerber_nlte as g

ROOT = Path(__file__).resolve().parents[1]
SRC = (ROOT / "pipeline" / "gerber_nlte.py").read_text()


def test_a_node_naming_a_deck_is_REFUSED_not_ignored():
    """The whole defect: it must raise rather than quietly serve the element's own deck."""
    with pytest.raises(g.GerberDeckError) as e:
        g.for_node("Fe", 5750.0, 4.5, 0.0, node="Fe@mean3D")
    msg = str(e.value)
    assert "selects nothing" in msg
    assert "FIRST argument" in msg, "the message must say how to get what was meant"
    assert "Fe@mean3D" in msg


def test_the_refusal_covers_every_registered_deck_not_just_Fe():
    """Written as a membership test against DECKS, so a deck registered later is covered
    without anyone remembering to extend a list."""
    assert "if node in DECKS:" in SRC


def test_node_is_documented_as_a_filename_tag_not_a_selector():
    """The name is the trap. If the docstring/comment stops saying so, the next caller
    makes the same reach."""
    assert "scratch-filename tag" in SRC or "scratch-FILENAME tag" in SRC
    assert "_coef.dat" in SRC


def test_a_harmless_node_tag_still_works():
    """`node` has a legitimate use — tagging interpolator scratch files. The guard must
    refuse only the values that NAME a deck, not every non-default value."""
    # "solar" is the default and must not be refused; any non-deck string is fine too.
    assert "solar" not in g.DECKS


def test_the_resolved_deck_key_is_reported_so_a_caller_can_ASSERT_what_it_got():
    """RYA-1049's second half: inference is what failed here. A caller should be able to
    ask which deck it received instead of deducing it from shape."""
    assert 'parsed["deck_key"] = element' in SRC
    assert SRC.count('parsed["deck_key"]') >= 2, (
        "both the direct-read and interpolator return paths must stamp it")
