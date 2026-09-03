"""RYA-1055 — Fe II NLTE is a STATED CAPABILITY LIMIT, and it must be stated where read.

🔴 THE FINDING. `atom.fe607a` — the model atom BOTH registered Gerber decks load (`Fe` and
`Fe@mean3D`) — declares 607 levels and 12,635 bound-bound transitions, and **every single
transition is Fe I → Fe I**. The 58 Fe II levels carry none: they are a pure ionisation
reservoir (targets of Fe I photoionisation, ionising on to Fe III). bsyn applies departures
PER LINE and falls back to departure = 1 for any line whose levels are unidentified, so
every Fe II line synthesises in LTE whatever the label says, and **no line list can change
that**.

⚠️ WHY THIS IS A TEST AND NOT A COMMENT. RYA-1055's first version drew the opposite
conclusion — *"the deck is not the limit, the LINE LIST is"* — and scoped a project to
label a VALD Fe II list against this atom. There was nothing to label against. The
distinction between a line-list gap and a deck limit decides what work gets done, so the
diagnosis is now LOOKED UP from a measured table rather than asserted in prose, and these
tests pin both the table and the places it has to reach.

WHAT IS DELIBERATELY *NOT* ASSERTED HERE
-----------------------------------------
* That Fe II NLTE is unavailable full stop. It is unavailable **on this deck**. The
  ENGINE-A leg reads the MPIA/Bergemann per-line delta grid, which carries 6,400 Fe II
  rows, and its Fe II corrections are real. Over-widening the claim would relabel an
  honest product, so there is a test below that this stays narrow.
* Any abundance. Nothing in this ticket re-derives a value (RYA-161).
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import gerber_nlte as G  # noqa: E402

REACH = ROOT / "data/results/rya1055/atom_ion_reach.json"
FEED = ROOT / "data/products/solar/Fe.json"


# ── the measurement, and the table that quotes it ────────────────────────────────────

def test_the_committed_measurement_says_zero_fe_ii_bound_bound():
    """The number the whole disposition rests on, read from the artifact, not remembered."""
    d = json.loads(REACH.read_text())
    assert d["atom"] == "atom.fe607a"
    # md5 is what `Fe_gerber2023.prov.json` and RYA-1035 pin, so this ties the measurement
    # to the STAGED bytes rather than to some other copy of the file.
    assert d["md5"] == "d08dc8232ed68eec65f9bb6631e82ea8"
    # A SHORT READ would report zero for every stage and look exactly like the finding.
    assert d["n_bb_parsed"] == d["n_bb_declared"] == 12635
    by = {s["ion_stage"]: s for s in d["stages"]}
    assert by[1]["n_bb_both"] == 12635 and by[1]["n_levels"] == 548
    assert by[2]["n_levels"] == 58
    assert by[2]["n_bb_both"] == 0, "the finding"
    assert by[2]["n_bb_any"] == 0, "not even one level of an Fe II transition"
    assert d["n_cross_stage_bb"] == 0
    # The independent corroboration: the highest level index appearing anywhere in the
    # transition block is BELOW the first Fe II level, so the zero is not a parse artefact.
    assert d["highest_level_index_in_any_bb"] < by[2]["level_index_lo"] == 549


def test_the_module_table_and_the_committed_artifact_cannot_drift():
    """🔴 A constant and an artifact describing the same fact drift silently, and the STALE
    side is the one that passes (RYA-1084/1092 class). Hold them against each other."""
    d = json.loads(REACH.read_text())
    measured = {s["ion_stage"]: s["nlte_capable"] for s in d["stages"]}
    assert G.ATOM_ION_REACH[d["atom"]] == measured
    assert G.ATOM_ION_REACH_ARTIFACT == "data/results/rya1055/atom_ion_reach.json"
    assert (ROOT / G.ATOM_ION_REACH_ARTIFACT).exists()


def test_the_generator_is_registered_and_reproduces_the_artifact_where_the_deck_is_staged():
    """RYA-686: the claim must be re-runnable, not merely committed.

    ⚠️ The atom is a VENDOR deck, not a repo file. Where it is not staged this SKIPS —
    and a skipped test is not a passing test, so it says so out loud rather than
    quietly counting as coverage.
    """
    gen = ROOT / "scripts/rya1055_atom_ion_reach.py"
    assert gen.exists()
    manifest = (ROOT / "data/results/GENERATORS.yaml").read_text()
    assert "rya1055/atom_ion_reach.json" in manifest
    assert "scripts/rya1055_atom_ion_reach.py" in manifest

    from config.constants import codex_path
    atom = codex_path("grids.gerber_ts") / "atom.fe607a"
    if not atom.exists():
        pytest.skip(f"atom.fe607a not staged at {atom} — measurement NOT re-run here")
    r = subprocess.run([sys.executable, str(gen), "--check"], cwd=ROOT,
                       capture_output=True, text=True)
    assert r.returncode == 0, f"re-measurement DIFFERS from the artifact:\n{r.stderr}"


# ── the accessor ─────────────────────────────────────────────────────────────────────

def test_the_verdict_is_stage_local_not_element_global():
    """Fe I keeps its NLTE leg. Sweeping the whole element up would delete a real product."""
    assert G.nlte_ion_capability("Fe", "I")[0] is True
    assert G.nlte_ion_capability("Fe", "II")[0] is False


def test_both_fe_decks_answer_the_same_because_they_load_the_same_atom():
    """Reach hangs off the ATOM, so the 1D and <3D> decks CANNOT disagree — which is what
    makes one entry cover model-registry members 4 and 6 at once."""
    assert G.DECKS["Fe"]["atom"] == G.DECKS["Fe@mean3D"]["atom"] == "atom.fe607a"
    assert G.atom_ion_reach(G.DECKS["Fe"]["atom"]) is \
           G.atom_ion_reach(G.DECKS["Fe@mean3D"]["atom"])


def test_an_unmeasured_atom_reads_UNMEASURED_and_never_False():
    """🔴 None is not False. Al's atom has never been measured; answering False would state
    a physical fact we have not established, and the guard branches on the difference."""
    assert G.atom_ion_reach("atom.al_qmh") is None
    assert G.nlte_ion_capability("Al", "I")[0] is None
    assert G.nlte_ion_capability("Mg", "II")[0] is None      # no deck registered at all
    assert G.nlte_ion_capability("Fe", "nonsense")[0] is None


def test_the_accessor_takes_the_ion_in_every_spelling_the_repo_uses():
    """⚠️ `parse_ion` returns the stage as an INT, not a roman string. Reading it as roman
    made the first draft answer None for every input while looking correct — so the
    spellings are pinned rather than assumed."""
    for spelling in ("II", "2", 2, 2.0):
        assert G.nlte_ion_capability("Fe", spelling)[0] is False, spelling


def test_the_limit_names_its_consequence_its_bound_and_its_exception():
    """A flag must carry its consequence. Naming the atom is not a statement about the
    number, and 'unavailable' without the ENGINE-A carve-out is a second wrong claim."""
    t = G.FE_II_NLTE_LIMIT
    assert "departure = 1" in t                 # the mechanism
    assert "Reported LTE" in t                  # what the number IS
    assert "no line list can enable it" in t    # why labelling cannot fix it
    assert "Lind" in t and "Amarsi" in t        # the literature bound, cited
    assert "ENGINE-A" in t                      # the carve-out
    assert "<3D>-LTE IS available" in t         # the axis Fe II CAN have
    assert "RYA-1055" in t


# ── the guard: the diagnosis must be looked up, not hardcoded ────────────────────────

def _linelist(ion_label):
    import numpy as np
    dt = [("turbospectrum_species", "f8"), ("element", "U8"), ("nlte", "U2"),
          ("nlte_label_low", "U8"), ("nlte_label_up", "U8"), ("wave_A", "f8")]
    return np.array([(26.0, ion_label, "F", "none", "none", 5000.0)], dtype=dt)


def test_an_unsupported_stage_is_diagnosed_as_a_DECK_limit_not_a_linelist_gap():
    """🔴 THE SENTENCE THAT WAS WRONG. This message used to close, unconditionally, with
    'This is a LINE-LIST coverage gap, not a deck failure — the deck is fine'. Acting on
    that is what scoped a VALD Fe II labelling project against an atom with nothing to
    label against."""
    with pytest.raises(G.GerberDeckError) as e:
        G.assert_linelist_supports_nlte(_linelist("Fe 2"), 26, "Fe", 4200.0, 6910.0,
                                        ion="II")
    msg = str(e.value)
    assert "NO LINE LIST CAN FIX THIS" in msg
    assert "ZERO Fe II bound-bound transitions" in msg
    assert "LINE-LIST coverage gap, not a deck failure" not in msg


def test_a_supported_stage_still_reads_as_a_linelist_gap():
    """The control. Fe I genuinely IS a labelling gap, and the old sentence was right for
    it — a fix that reported 'deck limit' for everything would be no better."""
    with pytest.raises(G.GerberDeckError) as e:
        G.assert_linelist_supports_nlte(_linelist("Fe 1"), 26, "Fe", 4200.0, 6910.0,
                                        ion="I")
    msg = str(e.value)
    assert "LINE-LIST coverage gap, not a deck failure" in msg
    assert "NO LINE LIST CAN FIX THIS" not in msg


def test_an_unmeasured_species_refuses_to_name_which_side_it_is():
    """Where reach is unknown the message must say UNKNOWN. Guessing here is how the
    original wrong diagnosis got written in the first place."""
    with pytest.raises(G.GerberDeckError) as e:
        G.assert_linelist_supports_nlte(_linelist("Fe 1"), 26, "Fe", 4200.0, 6910.0,
                                        ion=None)
    msg = str(e.value)
    assert "is NOT established for this species" in msg
    assert "LINE-LIST coverage gap, not a deck failure" not in msg
    assert "NO LINE LIST CAN FIX THIS" not in msg


# ── item 1: the stamp reaches every live Fe II product ───────────────────────────────

def test_every_live_fe_ii_product_carries_the_capability_limit():
    """Ryan, 2026-09-03: *stamp the limit ... in every Fe II product's
    science_provenance*. It is a property of the DECK WE SHIP, so it belongs on every
    Fe II product regardless of which treatment produced that number."""
    feed = json.loads(FEED.read_text())
    fe2 = [p for p in feed["products"] if p.get("ion") == "II"]
    assert fe2, "no live Fe II products — this test has stopped measuring anything"
    for p in fe2:
        cap = (p.get("science_provenance") or {}).get("nlte_capability")
        assert cap, f"unstamped: {p['band']}/{p['holding']}/{p['treatment']}"
        assert cap["fe_ii_nlte_available_on_gerber_deck"] is False
        assert cap["limit"] == G.FE_II_NLTE_LIMIT, "the stamp must quote the one source"
        assert cap["measurement"].startswith("data/results/rya1055/atom_ion_reach.json")


def test_the_stamp_stays_off_fe_I_and_names_the_right_nlte_source_per_product():
    """🔴 THE OVER-WIDENING TEST. Fe II ENGINE-A takes its departures from the
    MPIA/Bergemann per-line grid, NOT from atom.fe607a, so its NLTE label is honest.
    A blanket 'Fe II NLTE unavailable' on that product would replace one wrong statement
    with another."""
    feed = json.loads(FEED.read_text())
    for p in feed["products"]:
        cap = (p.get("science_provenance") or {}).get("nlte_capability")
        if p.get("ion") != "II":
            assert cap is None, f"Fe I product stamped: {p['treatment']}"
            continue
        assert cap["this_product_takes_nlte_from_the_gerber_deck"] is False, (
            "no live Fe II product should be taking departures from the Gerber deck — "
            "the RYA-1050 pool guard refuses to emit one")
        src = cap["nlte_source_for_this_product"]
        if p["treatment"] == "ENGINE-A":
            assert "Bergemann" in src and "NOT atom.fe607a" in src
        elif p["treatment"] in ("1D-LTE", "ENGINE-B"):
            assert src.startswith("n/a")


def test_no_live_fe_ii_product_is_on_an_nlte_scale_it_cannot_have():
    """Item 3 — the label audit. The per-line layer is the sharper check: `scale` is
    1D-NLTE on 159 Fe I rows and on ZERO Fe II rows, so the Gerber route never touched
    Fe II in anything we ship."""
    import csv
    rows = list(csv.DictReader(
        (l for l in (ROOT / "data/products/solar/Fe_perline.csv").read_text().splitlines()
         if not l.startswith("#"))))
    fe2 = [r for r in rows if r["ion"] == "II"]
    assert fe2, "no Fe II per-line rows — this test has stopped measuring anything"
    assert not [r for r in fe2 if "NLTE" in r["scale"].upper()], (
        "an Fe II per-line row on an NLTE scale: the Gerber deck cannot produce one, so "
        "either a phantom departure was applied or a label is wrong")
    assert [r for r in rows if r["ion"] == "I" and "NLTE" in r["scale"].upper()], (
        "the control: Fe I DOES carry 1D-NLTE rows, so the assertion above is measuring "
        "the ion split and not an empty column")


# ── item 2: the two published cells are annotated, not deleted ───────────────────────

def test_the_two_published_engine_b_nlte_fe_ii_cells_are_still_present():
    """ANNOTATE, DO NOT DELETE (Ryan, 2026-08-30 / 2026-09-03). They are the only Fe II
    Engine-B numbers we hold, and deleting them destroys the record that the question was
    asked."""
    import csv
    rows = list(csv.DictReader(
        open(ROOT / "data/results/rya783/fe_product_matrix.csv")))
    cells = {(r["band"], round(float(r["A"]), 3)): r for r in rows
             if r["ion"] == "II" and r["treatment"] == "ENGINE-B-NLTE"}
    assert ("VIS", 7.470) in cells and ("red-optical", 7.461) in cells
    assert int(cells[("VIS", 7.470)]["n_lines"]) == 8
    assert int(cells[("red-optical", 7.461)]["n_lines"]) == 2


def test_the_matrix_generator_flags_exactly_those_two_cells():
    """The annotation is DERIVED from the same accessor the guard and the stamp read, so
    the report cannot state a reach the pipeline disagrees with. Run against the committed
    matrix it must flag the two Fe II cells and leave both Fe I ones alone."""
    import pandas as pd
    from scripts.rya783_fe_matrix_report import _annotate_capability
    d = _annotate_capability(
        pd.read_csv(ROOT / "data/results/rya783/fe_product_matrix.csv"))
    flagged = d[d.capability_note.astype(str) != ""]
    assert len(flagged) == 2, f"expected 2 flagged cells, got {len(flagged)}"
    assert set(flagged.ion) == {"II"}
    assert set(flagged.treatment) == {"ENGINE-B-NLTE"}
    assert sorted(round(float(x), 3) for x in flagged.A) == [7.461, 7.470]
    for note in flagged.capability_note:
        assert note.startswith("CANNOT BE SHOWN TO BE NLTE")
        assert "Retained, not deleted" in note
    # the control: Fe I ENGINE-B-NLTE is REAL NLTE and must not be swept up
    fe1 = d[(d.ion == "I") & (d.treatment == "ENGINE-B-NLTE")]
    assert len(fe1) == 2 and (fe1.capability_note == "").all()


def test_the_annotation_record_exists_and_is_registered():
    md = ROOT / "data/results/rya783/CAPABILITY_ANNOTATIONS.md"
    assert md.exists()
    t = md.read_text()
    assert "7.470" in t and "7.461" in t
    assert "cannot be shown to be nlte" in t.lower()
    # It must NOT overclaim: the atmosphere explains the observed per-line difference.
    assert "not of applied departures" in t
    assert "rya783/CAPABILITY_ANNOTATIONS.md" in \
        (ROOT / "data/results/GENERATORS.yaml").read_text()


# ── the limit reaches the model registry, the other place a reader looks ─────────────

def test_the_gerber_nlte_registry_members_state_the_limit():
    """*"beside the deck registration"* is two places: `DECKS` in code, and the model
    registry a product's `model_grid` is resolved from."""
    import csv
    rows = {r["model_id"]: r for r in
            csv.DictReader(open(ROOT / "data/catalog/model_registry.csv"))}
    for mid in ("4", "6"):                       # the two Gerber NLTE members
        assert rows[mid]["model_family"] == "gerber" and "NLTE" in rows[mid]["scale"]
        assert "RYA-1055" in rows[mid]["notes"]
        assert "LTE-EQUIVALENT" in rows[mid]["notes"]
    # model 5 is the axis Fe II CAN have, and it must not be confused for the other one
    assert "ION-AGNOSTIC" in rows["5"]["notes"]
    assert "NEVER label it 3D-NLTE" in rows["5"]["notes"]
    # the standing phantom-departure check rides on the mandatory pair
    assert "model 5 == model 6 on Fe II BY CONSTRUCTION" in rows["6"]["notes"]
