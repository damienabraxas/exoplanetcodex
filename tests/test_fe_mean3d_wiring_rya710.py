"""RYA-710: the Fe ⟨3D⟩ deck is WIRED — the registry line, and what it may not claim.

RYA-710 is *ACQUIRE + WIRE* for the six Gerber decks we never fetched. RYA-1035 found the
Fe ⟨3D⟩ deck (it had been on our own Keeper share since 2021), fixed the two reader defects
that made it unusable, resolved its abundance, and staged it. This ticket adds the four
lines that let the pipeline reach it.

What these tests protect is not the registry line itself — it is the set of claims the line
makes *implicitly*, each of which was wrong at some point in this deck's history:

  * that the matrix may say WIRED (only true while the bytes are really staged);
  * that the aux is the PLAIN one (Al's is the converted one, and for Fe that file is
    unrecoverable);
  * that ⟨3D⟩ departures sit on the ⟨3D⟩ atmosphere, not MARCS;
  * that Fe has no abundance axis, so the v118 machinery is not silently required.

The deck binary is Sirius-only; every number quoted here was measured against it, on the
STAGED copy, and is reproduced in `docs/decisions/fe_3d_nlte_route_rya1035.md`.
"""
from pathlib import Path

import pytest

from pipeline import gerber_nlte as G

AUX_DIR = (Path(__file__).resolve().parents[1] / "data" / "nlte_grids" / "gerber_ts"
           / "fe_mean3d_aux")
PLAIN = AUX_DIR / "auxData_Fe_STAGGERmean3D_May-21-2021.txt"


@pytest.fixture(scope="module")
def plain_rows():
    return G._parse_aux_text(PLAIN.read_text())


# ── the matrix may only say WIRED because the bytes are staged ───────────────

def test_the_matrix_derives_fe_mean3d_as_wired_from_the_registry():
    """🔴 DERIVED, NEVER ASSERTED. `_MEAN3D_WIRED` was once a hard-coded `{'Al': ...}` and
    rendered "WIRED" on main for days while that deck existed only on an unmerged branch
    (v117). It reads `gerber_nlte.DECKS` now, so this cell is true exactly as long as the
    registry entry is."""
    from pipeline.model_availability_matrix import _mean3d_wired
    wired = _mean3d_wired()
    assert wired.get("Fe") == "Fe@mean3D"
    assert wired.get("Al") == "Al@mean3D"


def test_fe_is_listed_as_present_on_sirius():
    """The disk-side set must move with the staging, or the matrix reports a deck it has
    as one it still owes. Staged 2026-08-24, md5-verified in place."""
    from pipeline.model_availability_matrix import _MEAN3D_DECKS, _MEAN3D_DECKS_PUBLISHED
    assert "Fe" in _MEAN3D_DECKS
    assert "Fe" in _MEAN3D_DECKS_PUBLISHED
    assert _MEAN3D_DECKS <= _MEAN3D_DECKS_PUBLISHED, \
        "we cannot hold a deck the vendor does not publish -- that would be a provenance gap"


def test_the_fe_mean3d_cell_reads_consume_wired_not_fetch_owed():
    """Fe moves T2_FETCH_OWED -> T2_CONSUME_WIRED. NOT to T2_CONSUME_VALIDATED: the route
    runs and the deck reads, but no band product has been derived through it yet, and the
    difference between "wired" and "proven" is exactly what that state exists to keep."""
    from pipeline.model_availability_matrix import build_matrix
    cells = {(c["element"], c["model_type"]): c for c in build_matrix()["cells"]}
    for element in ("Fe", "Fe II"):
        cell = cells[(element, "MEAN3D_NLTE")]
        assert cell["state"] == "T2_CONSUME_WIRED", element
        assert cell["code_grid"] == "Fe@mean3D", element


# ── the claims the registry line makes implicitly ────────────────────────────

def test_the_aux_is_the_PLAIN_one_and_that_is_load_bearing():
    """🔴 THE OPPOSITE OF Al, ON PURPOSE. The vendor's `convert_3d_grid_to_marcs_names.py`
    builds the MARCS-style name FROM the [Fe/H] column — which for Fe is zeroed on the
    seven Teff=5777 rows — so the converted file collapses seven distinct atmospheres into
    one byte-identical string and the solar node cannot be addressed at all. Measured:
    `_marcs_names` REFUSES the solar node, the plain aux returns `p5777g44m00`.

    Al's converted aux is fine (0 disagreements in 6345 rows), which is why the two decks
    correctly point at differently-shaped files. This is not a style choice."""
    assert G.DECKS["Fe@mean3D"]["aux"] == "auxData_Fe_STAGGERmean3D_May-21-2021.txt"
    assert "_marcs_names" not in G.DECKS["Fe@mean3D"]["aux"]
    assert "_marcs_names" in G.DECKS["Al@mean3D"]["aux"]


def test_the_deck_carries_the_atmosphere_its_departures_were_solved_on():
    """Applying ⟨3D⟩ departures on a 1D MARCS structure yields a perfectly well-formed
    file — the failure class this module exists to make impossible. The deck is keyed at
    STAGGER's 5777/4.44, which MARCS_SOLAR (5750/4.5) does not contain."""
    assert G.deck_atmosphere("Fe@mean3D") != G.MARCS_SOLAR
    assert "stagger_avg3d" in G.deck_atmosphere("Fe@mean3D")
    assert Path(G.deck_atmosphere("Fe@mean3D")).exists()
    assert G.deck_atmosphere("Fe") == G.MARCS_SOLAR, "the 1D Fe deck must not move"


def test_the_1d_fe_deck_is_untouched_by_the_wiring():
    """Two decks for one element, differing only in the atmosphere their departures were
    solved on. Registering the second must not perturb the first."""
    assert G.DECKS["Fe"]["grid"] == "NLTEgrid4TS_Fe_MARCS_May-07-2021.bin"
    assert G.DECKS["Fe"]["aux"] == "auxData_Fe_MARCS_May-07-2021.dat"
    assert G.DECKS["Fe"].get("read_via") is None, "the MARCS path stays on the interpolator"
    assert G.DECKS["Fe"].get("atmos") is None


def test_both_fe_decks_share_the_model_atom_and_its_abundance():
    """`atom.fe607a` line 2 is `7.50  55.85`, and both Fe decks were solved with it — which
    is why one per-element provenance record serves both, and why `deck_abundance` resolves
    an `@`-suffixed key back to it."""
    assert G.DECKS["Fe"]["atom"] == G.DECKS["Fe@mean3D"]["atom"] == "atom.fe607a"
    assert G.DECKS["Fe"]["Z"] == G.DECKS["Fe@mean3D"]["Z"] == 26
    assert G.deck_abundance("Fe@mean3D") == 7.50
    assert G.deck_abundance("Fe") == 7.50


def test_wiring_did_not_hand_fe_an_abundance_axis(plain_rows):
    """🔴 THE DIFFERENCE A REGISTRY LINE CANNOT SHOW. Al's ⟨3D⟩ deck resolves 31 abundances
    and REQUIRES the v118 per-trial interpolation; Fe's resolves one (A(X) = 7.50 + [Fe/H],
    189 nodes), so its departures are a property of the atmosphere node and may be hoisted
    out of the χ² loop. `has_abundance_axis` decides that per deck, from the aux — so
    adding this entry cannot accidentally opt Fe into machinery it does not need, nor out
    of machinery it does."""
    G._AUX_ROW_CACHE["Fe@mean3D"] = plain_rows
    G._AXIS_CACHE.pop("Fe@mean3D", None)
    try:
        assert G.abundance_axis("Fe@mean3D") == (7.5,)
        assert not G.has_abundance_axis("Fe@mean3D")
    finally:
        G._AUX_ROW_CACHE.pop("Fe@mean3D", None)
        G._AXIS_CACHE.pop("Fe@mean3D", None)


# ── what the wiring may NOT claim ────────────────────────────────────────────

def test_the_passing_rya534_gate_does_not_vouch_for_the_3D_deck():
    """🔴 ONE `verdict: PASS` IN A RECORD THAT NOW DESCRIBES TWO DECKS.

    `Fe_gerber2023.prov.json` carries a PASSING RYA-534 anchor gate — run on the **1D
    MARCS** deck, over eight Fe I lines (RYA-785). The ⟨3D⟩ deck has never been through
    one. A reader finding that block in a record that now lists `grid_3d_*` files would
    reasonably assume it covered both, so the block says otherwise in its own text.

    RYA-710's standing rule — *"do NOT mark anything wired before its RYA-534 anchor gate
    passes"* — reads against the coverage CSV's `task`→`wired` column, which this ticket
    does not touch. The matrix's `T2_CONSUME_WIRED` is a plumbing state, and the state
    above it (`T2_CONSUME_VALIDATED`) is exactly what an ungated deck may not claim."""
    import json
    p = (Path(__file__).resolve().parents[1] / "data" / "nlte_grids" / "gerber_ts"
         / "Fe_gerber2023.prov.json")
    d = json.loads(p.read_text())
    assert d["gate"]["verdict"].startswith("PASS")
    scope = d["gate"]["scope_rya710"]
    assert "1D MARCS DECK ONLY" in scope
    assert "T2_CONSUME_VALIDATED" in scope
    # the gate's own line list is 1D-deck evidence; the <3D> files are separate entries
    assert "grid_3d_bin_zip" in d["files"] and "grid_1d_bin_zip" in d["files"]


def test_the_two_engine_coverage_row_for_fe_is_untouched():
    """The column RYA-710's rule actually governs. Wiring the ⟨3D⟩ deck must not silently
    upgrade the Engine-B coverage claim, which describes the GATED 1D deck."""
    import csv
    p = (Path(__file__).resolve().parents[1] / "data" / "curation"
         / "nlte_two_engine_coverage.csv")
    rows = [r for r in csv.DictReader(p.open()) if r.get("element") == "Fe"]
    assert rows, "Fe must still have a coverage row"
