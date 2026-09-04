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


# ── item 3: the live-label confirmation, and the scale it is quoted on ───────────────

AUDIT = ROOT / "data/results/rya1055/fe2_label_audit.json"


def test_the_live_label_audit_is_clean_and_reproduces():
    """Item 3, run rather than believed. The audit exits non-zero if it finds a problem,
    so this is the guard as well as the record."""
    # `--check` re-derives from the LIVE feed and compares; it writes nothing, so the
    # verifier cannot dirty the artifact it verifies on a timestamp alone.
    r = subprocess.run([sys.executable, "scripts/rya1055_fe2_label_audit.py", "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    a = json.loads(AUDIT.read_text())
    assert a["problems"] == []
    assert a["n_live_fe_ii_products"] == a["n_stamped"] == 10
    assert a["n_taking_nlte_from_the_gerber_deck"] == 0
    assert a["per_line_rows_on_an_nlte_scale"] == {"I": 159, "II": 0}


def test_the_balance_is_matched_on_the_full_identity_not_a_looser_key():
    """🔴 Without `selector`, HARPS VIS Fe I offers TWO partners 0.196 dex apart
    (DEEPGRADED 7.535, DEEPGRADED-LOCALRENORM 7.339) and the reader picks whichever suits.
    Pin that the key carries every field, and that no pair came back ambiguous."""
    src = (ROOT / "scripts/rya1055_fe2_label_audit.py").read_text()
    for field in ("band", "instrument", "holding", "tier", "selector", "route",
                  "treatment"):
        assert f'"{field}"' in src.split("KEY = (")[1].split(")")[0], field
    a = json.loads(AUDIT.read_text())
    seen = [(b["holding"], b["band"], b["treatment"]) for b in
            a["ionisation_balance_scale_matched"]]
    assert len(seen) == len(set(seen)), "a cell matched more than one Fe I partner"


def test_the_audit_refuses_to_quote_the_balance_against_the_3d_anchor():
    """⚠️ The Fe II kurucz2005 VIS cell is 7.466 — the anchor's own digits — on a
    DIFFERENT scale. The record must carry the caveat, not the coincidence."""
    a = json.loads(AUDIT.read_text())
    assert "3D-NLTE" in a["anchor_caveat"] and "1D-NLTE" in a["anchor_caveat"]
    assert "RYA-669" in a["anchor_caveat"]
    # and the balance itself must be stated scale-matched, with its weakness admitted
    assert "one scale by construction" in a["balance_note"]
    assert "does not CONTRADICT" in a["balance_note"]


# ── the literature bound is a MEASURED claim, so it must be re-derivable ─────────────

def test_the_mpia_fe_ii_numbers_in_the_limit_are_reproduced_from_the_grid():
    """🔴 THE ERROR THIS TEST EXISTS TO PREVENT, BECAUSE IT WAS MADE. A first draft of the
    limit said the MPIA Fe II corrections were "±0.002 dex across the grid at solar" — read
    off the first five rows of a truncated `describe()` and generalised. The tail is
    +0.016, eight times that, and the wrong number was on its way into every published
    Fe II product's provenance.

    So the numbers the limit quotes are re-derived here from the committed grid. A prose
    bound with no test beside it is a second source of truth for a measured fact.
    """
    import pandas as pd
    df = pd.read_csv(ROOT / "data/nlte_grids/Fe_Bergemann_MPIA.csv")
    node = (df.teff_K == 5800) & (df.logg.isin([4.3, 4.5])) & (df.feh == 0.0)
    fe2 = df[(df.ion == "II") & node]
    fe1 = df[(df.ion == "I") & node]
    assert len(fe2) == 160 and len(fe1) == 504

    assert round(float(fe2.delta_nlte.median()), 3) == 0.000
    assert round(float(fe2.delta_nlte.max()), 3) == 0.016
    assert round(float(fe2.delta_nlte.min()), 3) == -0.002
    assert int((fe2.delta_nlte.abs() <= 0.005).sum()) == 146
    # ⚠️ dropna() FIRST. 8 of the 160 node rows are NaN and 3 lines carry no value at
    # all; a `sort_values().tail()` puts those last and makes the third-from-last look
    # like the maximum — which is how the wrong tail got into a first draft of this test.
    tail = fe2.groupby("wave_A").delta_nlte.mean().dropna().nlargest(3).index.tolist()
    assert sorted(round(w, 3) for w in tail) == [4923.932, 4924.921, 5169.033]
    # the Fe I control — the comparison that makes "small" mean something
    assert round(float(fe1.delta_nlte.median()), 3) == 0.011
    assert round(float(fe1.delta_nlte.max()), 3) == 0.040

    t = G.FE_II_NLTE_LIMIT
    for s in ("+0.000 dex at the solar node", "146 of 160", "+-0.005", "+0.016",
              "4923.932", "4924.921", "5169.033", "median +0.011", "+0.040"):
        assert s in t, f"the limit no longer quotes {s!r}"
    # and it must not claim a BOUND it does not have
    assert "median statement and not a bound" in t
    assert "+-0.002 dex" not in t


def test_what_the_live_fe_ii_products_ACTUALLY_applied_is_read_from_their_own_artifacts():
    """🔴 QUOTE THE PRODUCT'S OWN ARTIFACT, NOT A NEIGHBOURING ONE.

    A draft of the stamp cited "+0.001 dex on 6147.7341 / 6238.3859 / 6247.5570" from
    `Fe_perline.csv`. Those numbers are real — and they belong to a DIFFERENT POOL, the
    RYA-489 replication product's 11-line 5256-6456 A Fe II set. The LIVE VIS band
    product fits 4233.162 / 4303.170 / 4583.829 and applies -0.001/-0.002: same element,
    same ion, same treatment, OPPOSITE SIGN. Both pools are pinned here so the two can
    never again be confused for each other.
    """
    import pandas as pd
    bp = ROOT / "data/results/band_products"
    served = {}
    for f in sorted(bp.glob("FeII_*_ENGINE-A_lines.csv")):
        d = pd.read_csv(f)
        s = d[d.in_aggregate == True]                                   # noqa: E712
        served[f.name] = (len(d), len(s),
                          sorted({round(float(x), 4) for x in s.nlte_delta_dex}))
    assert len(served) == 5, sorted(served)
    for name, (n, k, deltas) in served.items():
        assert deltas and all(-0.0021 <= x <= -0.0009 for x in deltas), (name, deltas)
        assert k < n, f"{name}: MPIA served every line — the n-drop confound is gone?"
    # the VIS pools: 3 of 9 served, at -0.001/-0.002
    vis = [v for k, v in served.items() if "_4200_6910_" in k]
    assert len(vis) == 3 and all(v[0] == 9 and v[1] == 3 for v in vis), vis
    # the near-UV pools: 7 of 12 served, at -0.001 — RYA-1113's n=7-vs-12
    nuv = [v for k, v in served.items() if "_3000_3780_" in k]
    assert len(nuv) == 2 and all(v[0] == 12 and v[1] == 7 for v in nuv), nuv

    # and the OTHER pool, so the two stay distinguishable
    import csv
    rows = list(csv.DictReader(
        (l for l in (ROOT / "data/products/solar/Fe_perline.csv").read_text().splitlines()
         if not l.startswith("#"))))
    fe2 = [r for r in rows if r["ion"] == "II" and r["arm"] == "VIS"]
    lte = {r["wavelength_air_A"]: r["A_X_line"] for r in fe2 if r["engine"] == "1D-LTE"}
    eng = {r["wavelength_air_A"]: r["A_X_line"] for r in fe2
           if r["engine"] == "ENGINE-A" and r["status"] == "in_aggregate"}
    # compare as FLOATS: the CSV writes 6247.557, not 6247.5570.
    assert sorted(round(float(w), 4) for w in eng) == [6147.7341, 6238.3859, 6247.557]
    for w, a in eng.items():
        d = float(a) - float(lte[w])
        assert abs(d - 0.001) < 5e-5, f"{w}: ENGINE-A minus 1D-LTE = {d:+.6f}"
    # THE POINT: the two pools do not overlap at all.
    live_vis = {4233.162, 4303.170, 4583.829}
    assert live_vis & {round(float(w), 3) for w in eng} == set()


def test_the_rya1113_contradiction_is_resolved_by_the_per_line_artifacts():
    """RYA-1113 asked, on this ticket: if Fe II NLTE is structurally unavailable, how does
    a live n=7 Fe II NLTE leg exist at all? It offered two answers — the leg is not real
    NLTE, or the atom has changed. Both are wrong, and the per-line artifact RYA-908 has
    since emitted says which: the source is the MPIA grid, not `atom.fe607a`.

    ⚠️ RYA-1113's OTHER half stands: n=7 against an LTE sibling's n=12 IS a pool change,
    so the published near-UV Fe II NLTE delta is confounded. That is not fixed here.
    """
    import pandas as pd
    f = (ROOT / "data/results/band_products/FeII_3000_3780_kpno_solar_atlas_"
         "solar_kpno_kurucz2005_corrected_SYNTH_DEEPGRADED_ENGINE-A_lines.csv")
    d = pd.read_csv(f)
    served = d[d.in_aggregate == True]                                  # noqa: E712
    assert len(d) == 12 and len(served) == 7, "RYA-1113's n=7 of 12"
    srcs = set(served.nlte_source.dropna())
    assert srcs and all("Bergemann MPIA" in s for s in srcs), srcs
    assert not any("erber" in s or "fe607a" in s for s in set(d.nlte_source.dropna()))
    assert sorted({round(float(x), 4) for x in served.nlte_delta_dex}) == [-0.001]
    # the 5 that dropped out did so for SERVICE coverage, not physics
    dropped = d[d.in_aggregate != True]                                 # noqa: E712
    assert len(dropped) == 5
    assert all("NOT-SERVED" in str(r) for r in dropped.excluded_reason), \
        sorted(set(dropped.excluded_reason))
    a = json.loads(AUDIT.read_text())
    assert "RESOLVED" in a["rya1113_contradiction"]
    assert "confounded with a 5-line pool change" in a["rya1113_contradiction"]


def test_no_fe_ii_band_product_names_the_gerber_deck_as_its_nlte_source():
    """The strongest form of the label audit: `nlte_source` is written at the point the
    correction is applied (RYA-880), so this is the RUN speaking, not a label being read.
    Every LTE leg must say so, and no leg may name the deck."""
    a = json.loads(AUDIT.read_text())
    per = a["band_product_per_line_nlte_sources"]
    #: 15 -> 16: RYA-1135 added the Fe II <3D>-LTE leg
    #: (`..._synth-mean3D-LTE-gerber-stagger_lines.csv`). It is a NEW artifact, not a
    #: changed one — and it is exactly the kind this test exists to police, so it is swept
    #: in rather than excluded: the per-row assertions below run on it too, and it must
    #: report `none — LTE, no departure applied` like every other non-ENGINE-A leg. The
    #: <3D>-mean atmosphere is ion-agnostic; only the departures are ion-specific, and Fe II
    #: takes none.
    assert len(per) == 16, len(per)
    assert any("synth-mean3D-LTE-gerber-stagger" in r["artifact"] for r in per), (
        "the RYA-1135 Fe II <3D>-LTE leg must be audited like every other Fe II product")
    for r in per:
        for s in r["nlte_source"]:
            assert "erber" not in s and "fe607a" not in s, (r["artifact"], s)
        if "_ENGINE-A_" not in r["artifact"]:
            assert r["n_with_a_nonzero_departure"] == 0, r["artifact"]
            assert r["nlte_source"] == ["none — LTE, no departure applied"], r["artifact"]


def test_the_two_disjoint_vis_fe_ii_pools_are_recorded_not_dropped():
    """🔴 TWO PUBLISHED "solar VIS Fe II DEEPGRADED" POOLS, ZERO OVERLAP, same window.

    `Fe_perline.csv` carries 11 lines at 5256.9-6456.4 A (RYA-870, sourced from
    rya847+rya877 — and it is RYA-877's pool that this ticket's own headline "0 of 11
    labelled" was measured on). The live band products carry 9 at 4233.2-4583.8 A. Both
    sit inside 4200-6910 A and they share NOT ONE line.

    It does NOT weaken the finding — a deck with zero Fe II bound-bound transitions is
    zero for ANY pool, which is the strength of a deck-level result over a line-list one.
    It is recorded because it is why the ticket's headline number describes a pool the
    live products no longer use, and because it is FOR RYA TO DISPOSITION, not for this
    ticket to change.
    """
    a = json.loads(AUDIT.read_text())["two_disjoint_vis_fe_ii_pools"]
    perline = a["Fe_perline.csv VIS Fe II (RYA-870, sourced rya847+rya877)"]
    live = a["live band product FeII_4200_6910 kpno molecfit DEEPGRADED"]
    assert len(perline) == 11 and len(live) == 9
    assert a["overlap"] == [], "the pools now overlap — re-read this finding"
    assert max(live) < min(perline), "the two windows no longer separate cleanly"
    assert "FOR RYA TO DISPOSITION" in a["note"]


def test_the_mpia_grid_fill_rates_are_recorded_with_their_fe_i_control():
    """⚠️ Context for how much the live Fe II NLTE label carries, invisible unless counted:
    the MPIA grid's Fe II half is ~120x more sparsely populated than its Fe I half (7.1%
    NaN vs 0.06%), and three Fe II lines are empty at the very node our products
    interpolate at. Not a defect in anything this ticket changes — recorded, with the
    Fe I control beside it so "sparse" means something."""
    g = json.loads(AUDIT.read_text())["mpia_grid_fill"]
    assert g["Fe I"]["rows"] == 20160 and g["Fe I"]["rows_nan"] == 12
    assert g["Fe II"]["rows"] == 6400 and g["Fe II"]["rows_nan"] == 456
    assert g["Fe II"]["nan_fraction"] > 50 * g["Fe I"]["nan_fraction"]
    assert g["Fe I"]["lines_empty_at_the_solar_node"] == []
    assert g["Fe II"]["lines_empty_at_the_solar_node"] == [4319.68, 4577.9, 5722.56]
    assert "FOR RYA TO DISPOSITION" in g["note"]
