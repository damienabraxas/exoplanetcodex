"""RYA-1222: the matrix expands from ledgers, and every cell gets a terminal status.

What is pinned here is not that the orchestrator runs -- it mostly does not, because
most cells are honestly blocked -- but that it never answers by omission. Each test
below corresponds to a way a cell could have disappeared, and two of them are ways a
cell DID disappear while this ticket was being built.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module", autouse=True)
def _kp(tmp_path_factory):
    """The RYA-1064 friction #3 guard, as `test_rya767_run_descriptor` sets it.

    `measure_band_ew` resolves the Kitt Peak atlas AT IMPORT. `setdefault` so a
    machine that really has the atlas staged uses it and CI does not need it.
    """
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))


@pytest.fixture(scope="module")
def rm():
    from pipeline import run_matrix
    return run_matrix


# ── the axes come from the ledgers, never from a literal ─────────────────────

def test_engine_axis_is_bound_to_the_model_registry(rm):
    """Every live treatment is claimed by a deck or named un-dispatchable WITH a reason."""
    rm.assert_registry_covered()
    live = rm.live_treatments()
    claimed = {t for toks in rm.DECK_EMITS.values() for t in toks}
    assert live == claimed | set(rm.NOT_DISPATCHABLE)
    for token, why in rm.NOT_DISPATCHABLE.items():
        assert len(why) > 40, f"{token} is excluded without a reason worth reading"


def test_a_new_live_treatment_fails_loud_rather_than_never_running(rm, monkeypatch):
    """The whole point of the binding: an unclaimed token must not be a silent gap."""
    real = rm.live_treatments()
    monkeypatch.setattr(rm, "live_treatments",
                        lambda: real | {"brand-new-deck-token"})
    with pytest.raises(rm.MatrixError, match="brand-new-deck-token"):
        rm.assert_registry_covered()


def test_a_deck_claiming_a_retired_treatment_fails_loud(rm, monkeypatch):
    """The other direction: a deck dispatched for a product nothing accepts."""
    survivors = set(rm.DECK_EMITS["ts-lte"]) | set(rm.NOT_DISPATCHABLE)
    monkeypatch.setattr(rm, "live_treatments", lambda: survivors)
    with pytest.raises(rm.MatrixError, match="no longer lists as live"):
        rm.assert_registry_covered()


def test_an_unknown_element_is_refused_against_elements_master(rm):
    with pytest.raises(rm.MatrixError, match="canonical targets"):
        rm.validate_element("Unobtainium")


def test_the_ion_axis_is_the_graded_pool_not_the_periodic_table(rm):
    """Si III exists in nature and in canonical_gf; it has no GRADED line, so no cell.
    (Si II gained graded lines from the RYA-1233 NIST ASD intake.)"""
    assert rm.graded_ions("Si") == ["I", "II"]
    assert rm.graded_ions("Fe") == ["I", "II"]


def test_every_holding_on_the_axis_is_a_registered_holding_of_that_star(rm):
    registered = {h["holding_id"] for h in rm.holdings_for("solar")}
    assert {d.holding for d in rm.expand("solar", "Si")} <= registered


def test_an_unregistered_system_is_refused_not_guessed(rm):
    with pytest.raises(rm.MatrixError, match="refuse, do not guess"):
        rm.holdings_for("betelgeuse")


# ── the two cells that vanished while this was being built ───────────────────

def test_harps_vis_is_on_the_matrix(rm):
    """🔴 REGRESSION. It was not, and it is where the published Fe anchor lives.

    `bands_for('harps')` returns VIS as 3780-6910 A from the instrument catalogue;
    `solar_harps` declares its span as 3782.6-6910.0. `HoldingSpec.covers()` demands
    TOTAL coverage, so 2.4 A at the blue edge made the single most productive
    holding in the repo silently absent from its own element's matrix.
    """
    cells = [d for d in rm.expand("solar", "Fe", engines=["ts-lte"])
             if d.holding == "solar_harps_molecfit_corrected" and d.band == "VIS"]
    assert cells, "solar_harps VIS produced no cell"
    for d in cells:
        assert d.lo_A == pytest.approx(3782.6), "the window must be clipped, not dropped"


def test_the_window_is_clipped_to_what_the_holding_serves(rm):
    """IAG serves 5001.1 A upward, so its VIS cell starts there -- as RYA-1218's own
    artifact stem (`SiI_5002_6910_iag_...`) independently says it did."""
    iag = [d for d in rm.expand("solar", "Si", ions=["I"], engines=["ts-lte"], methods=["profile-fit"])
           if d.holding == "solar_iag" and d.band == "VIS"]
    assert len(iag) == 1
    assert iag[0].lo_A == pytest.approx(5001.1)
    assert iag[0].hi_A == pytest.approx(6910.0)


def test_a_holding_that_inventories_its_own_coverage_is_not_clipped(rm):
    """`span_A is None` means the reader knows its segments; there is nothing to clip to."""
    assert rm._clip_to_holding(None, 3780.0, 6910.0) == (3780.0, 6910.0)


# ── rule 2: every cell gets a terminal status, none is dropped ───────────────

@pytest.fixture(scope="module")
def dry(rm, tmp_path_factory):
    return rm.run("solar", "Si", dry_run=True, interpreter=sys.executable,
                  ispec_dir="/nonexistent/ispec", echo=False,
                  report_dir=tmp_path_factory.mktemp("report"))


def test_the_counts_account_for_every_cell(dry, rm):
    assert sum(dry["counts"].values()) == dry["cells_total"] == len(dry["cells"])
    assert set(dry["counts"]) == set(rm.STATUSES)


def test_no_cell_is_silent(dry, rm):
    """A status with no reason is a cell that fell through a branch -- a build defect."""
    for c in dry["cells"]:
        assert c["status"] in rm.STATUSES
        if c["status"] != rm.DONE:
            assert c["reason"].strip(), f"{c} carries a status and no reason"


def test_resume_is_exactly_the_cells_whose_work_does_not_exist(dry, rm):
    expected = [c for c in dry["cells"] if c["status"] not in rm.TERMINAL_OK]
    assert len(dry["resume"]) == len(expected)


# RYA-1233: `test_the_band_edge_disagreement_blocks_rather_than_picks_a_winner` lived here
# and asked to be deleted, with a reason, the day the two tables agreed. They agree:
# pipeline/band_policy.py now reads its optical/NIR edges from config/synth_bands.yaml
# (`test_band_policy_edges_are_synth_bands_edges`). `_band_of` keeps the refusal, so a
# future second edge table still blocks rather than tiebreaks.


def test_the_report_is_written_and_reloadable(rm, tmp_path):
    doc = rm.run("solar", "Si", dry_run=True, engines=["ts-lte"], bands=["VIS"],
                 echo=False, report_dir=tmp_path)
    latest = tmp_path / "solar_Si_latest.json"
    assert latest.exists()
    assert json.loads(latest.read_text())["cells_total"] == doc["cells_total"]
    assert len(list(tmp_path.glob("solar_Si_2*.json"))) == 1, "no timestamped copy"


def test_render_is_ascii_only(dry, rm):
    """Cloudflare's WAF eats non-ASCII in the Linear comment these tables are pasted into."""
    rm.render(dry).encode("ascii")


# ── rule 4: the idempotency key ──────────────────────────────────────────────

def test_a_cell_with_no_recorded_hash_is_not_current(rm):
    """'We cannot tell what it was built from' is not 'yes' (no silent fallbacks)."""
    ok, why = rm.is_current([{"A": 7.5}], "abc", None)
    assert not ok and "no inputs_hash was ever recorded" in why


def test_a_stale_recorded_hash_is_not_current(rm):
    ok, why = rm.is_current([{"A": 7.5}], "new", {"inputs_hash": "old"})
    assert not ok and "recorded inputs_hash" in why


def test_no_published_product_is_not_current(rm):
    """A ledger entry for a product someone quarantined must not skip the cell."""
    ok, why = rm.is_current([], "abc", {"inputs_hash": "abc"})
    assert not ok and "no published product" in why


def test_a_matching_recorded_hash_with_products_present_is_current(rm):
    ok, _ = rm.is_current([{"A": 7.5}], "abc", {"inputs_hash": "abc"})
    assert ok


def test_the_ledger_key_separates_two_decks_over_one_window(rm):
    """The product key does not carry the deck; the cell key must."""
    from pipeline.run_descriptor import RunDescriptor
    mk = lambda deck: rm.ledger_key("solar", RunDescriptor(  # noqa: E731
        element="Fe", ion="I", instrument="harps", holding="solar_harps",
        lo_A=3782.6, hi_A=6910.0, engine_deck=deck), "VIS")
    assert mk("ts-lte") != mk("gerber-nlte")


def test_an_unreadable_ledger_costs_a_rerun_never_a_skip(rm, tmp_path, monkeypatch, capsys):
    bad = tmp_path / "inputs_hashes.json"
    bad.write_text("{not json")
    monkeypatch.setattr(rm, "LEDGER", bad)
    assert rm.load_ledger() == {}
    assert "unreadable" in capsys.readouterr().err


def test_the_recorded_hash_round_trips(rm, tmp_path, monkeypatch):
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")
    rm.record_inputs_hash("solar|Fe|I|VIS|harps|solar_harps|ts-lte", "H",
                          product_keys=["k1"], code_commit="abc")
    led = rm.load_ledger()
    assert led["solar|Fe|I|VIS|harps|solar_harps|ts-lte"]["inputs_hash"] == "H"
    assert rm.is_current([{"A": 1}], "H",
                         led["solar|Fe|I|VIS|harps|solar_harps|ts-lte"])[0]


def test_the_orchestrator_never_writes_into_the_science_feed(rm):
    """🔴 THE SPEC SAID TO, AND THE RYA-587 PUBLICATION CONTRACT REFUSES IT.

    `publish_product.write_feed` re-validates every product that is not
    BYTE-IDENTICAL to the row already on disk. The committed products are legacy
    rows with no `uncertainty` block and pass only via that retention clause, so
    adding one bookkeeping field to one of them gets the whole feed rejected. This
    pins the measurement, so that if the contract ever changes, whoever changes it
    finds out that this decision depended on it.
    """
    import copy
    from pipeline.uncertainty_contract import assert_publication_feed, UncertaintyError
    feed = json.loads((ROOT / "data" / "products" / "solar" / "Fe.json").read_text())
    previous = copy.deepcopy(feed)
    assert_publication_feed(copy.deepcopy(feed), previous=previous)   # unchanged: passes
    feed["products"][0]["inputs_hash"] = "HASH"
    with pytest.raises(UncertaintyError):
        assert_publication_feed(feed, previous=previous)
    assert not hasattr(rm, "stamp_inputs_hash"), \
        "the feed-writing path came back; it cannot work (see LEDGER)"


def _fixture_run(rm):
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element="Fe", ion="I", instrument="harps",
                      holding="solar_harps", lo_A=3782.6, hi_A=6910.0)
    return d, resolve(d, interpreter=sys.executable, ispec_dir="/x/ispec")


def test_the_inputs_hash_is_deterministic(rm):
    d, r = _fixture_run(rm)
    a = rm.inputs_hash(d, r, manifest_path=None)
    b = rm.inputs_hash(d, r, manifest_path=None)
    assert a == b and len(a) == 64


def test_the_inputs_hash_moves_when_a_stage_script_moves(rm, monkeypatch):
    """The `code_commit of the stage modules` term, done on the bytes that will run.

    Repo HEAD would move for edits to files this cell never loads, and would not
    move at all for an uncommitted edit that changes the answer.
    """
    d, r = _fixture_run(rm)
    before = rm.inputs_hash(d, r, manifest_path=None)
    real = rm._file_fingerprint
    monkeypatch.setattr(rm, "_file_fingerprint",
                        lambda p: "EDITED" if p.name.endswith("derive_band_products.py")
                        else real(p))
    assert rm.inputs_hash(d, r, manifest_path=None) != before


def test_the_inputs_hash_moves_with_the_descriptor(rm):
    """Two runs differing only in ion are two runs, and must not share a key."""
    from pipeline.run_descriptor import RunDescriptor, resolve
    d, r = _fixture_run(rm)
    other = RunDescriptor(element="Fe", ion="II", instrument="harps",
                          holding="solar_harps", lo_A=3782.6, hi_A=6910.0)
    ro = resolve(other, interpreter=sys.executable, ispec_dir="/x/ispec")
    assert rm.inputs_hash(d, r, manifest_path=None) != \
        rm.inputs_hash(other, ro, manifest_path=None)


def test_the_inputs_hash_moves_with_the_holdings_own_artifact(rm):
    """A re-staged spectrum invalidates the cell even when nothing else moved."""
    d, r = _fixture_run(rm)
    a = rm.inputs_hash(d, r, manifest_path="data/catalog/instrument_catalog.csv")
    b = rm.inputs_hash(d, r, manifest_path="data/catalog/model_registry.csv")
    assert a != b


def test_a_missing_input_is_named_not_silently_skipped(rm):
    assert rm._file_fingerprint(ROOT / "no" / "such" / "file").startswith("ABSENT:")


def test_the_cell_key_carries_the_holding_not_just_the_instrument(rm):
    from pipeline.run_descriptor import RunDescriptor
    mk = lambda h: rm.cell_prefix_key(RunDescriptor(  # noqa: E731
        element="Fe", ion="I", instrument="harps", holding=h,
        lo_A=3782.6, hi_A=6910.0))
    assert mk("solar_harps") != mk("solar_harps_molecfit_corrected")


def test_cell_products_match_only_the_decks_own_treatments(rm):
    from pipeline.run_descriptor import RunDescriptor
    d = RunDescriptor(element="Fe", ion="I", instrument="harps", holding="solar_harps",
                      lo_A=3782.6, hi_A=6910.0, engine_deck="gerber-nlte", pool="codex")
    feed = {"products": [
        {"element": "Fe", "ion": "I", "band": "VIS", "instrument": "harps",
         "holding": "solar_harps", "treatment": "ENGINE-B-NLTE", "A": 7.5,
         "route": "PROFILEFIT", "selector": "GRADED"},
        # RYA-1233: same everything, another POOL -- another cell's product
        {"element": "Fe", "ion": "I", "band": "VIS", "instrument": "harps",
         "holding": "solar_harps", "treatment": "ENGINE-B-NLTE", "A": 7.7,
         "route": "PROFILEFIT", "selector": "REFERENCE"},
        {"element": "Fe", "ion": "I", "band": "VIS", "instrument": "harps",
         "holding": "solar_harps", "treatment": "1D-LTE", "A": 7.4, "route": "PROFILEFIT"},
        # RYA-1233: same treatment on the OTHER route is another cell's product
        {"element": "Fe", "ion": "I", "band": "VIS", "instrument": "harps",
         "holding": "solar_harps", "treatment": "ENGINE-B-NLTE", "A": 7.6, "route": "SYNTH"},
    ]}
    hit = rm.cell_products(feed, d, "VIS")
    assert [(p["treatment"], p["route"]) for p in hit] == [("ENGINE-B-NLTE", "PROFILEFIT")]


# ── rules 1 and 3: refusals, and loud-fail-CONTINUE ──────────────────────────

def test_the_mean3d_decks_are_a_mandatory_pair(rm):
    """RYA-1040/542: differencing the survivor against 1D-LTE would report the
    1D -> mean-3D ATMOSPHERE shift as non-LTE physics."""
    with pytest.raises(rm.MatrixError, match="MANDATORY PAIR"):
        rm.expand("solar", "Si", engines=["gerber-mean3d"])
    rm.expand("solar", "Si", engines=["gerber-mean3d", "gerber-mean3d-lte"])


def test_an_unknown_deck_is_refused(rm):
    with pytest.raises(rm.MatrixError, match="unknown engine deck"):
        rm.expand("solar", "Si", engines=["ts-nlte-turbo-9000"])


def test_a_dry_run_executes_nothing(rm, monkeypatch, tmp_path):
    monkeypatch.setattr(rm, "_run_step",
                        lambda *a, **k: pytest.fail("a dry run dispatched a stage"))
    doc = rm.run("solar", "Si", dry_run=True, interpreter=sys.executable,
                 ispec_dir="/x/ispec", echo=False, report_dir=tmp_path)
    assert doc["counts"][rm.DONE] == 0


def test_one_failing_cell_does_not_stop_the_matrix(rm, monkeypatch, tmp_path):
    """Rule 3. An aborted matrix reports nothing about the cells it never reached,
    which is the same silent gap as a dropped cell wearing a crash."""
    monkeypatch.setattr(rm, "_run_step", lambda *a, **k: (False, "boom"))

    reaches_the_executor(rm, monkeypatch)
    doc = rm.run("solar", "Si", engines=["ts-lte"], bands=["VIS"],
                 interpreter=sys.executable, ispec_dir="/x/ispec",
                 echo=False, report_dir=tmp_path)
    assert doc["counts"][rm.FAILED] >= 2, "the matrix stopped at the first failure"
    assert sum(doc["counts"].values()) == doc["cells_total"]
    assert all(c["reason"] == "boom" for c in doc["cells"] if c["status"] == rm.FAILED)


def _fake_stage(rm, tmp_path, monkeypatch, calls):
    """Stages that 'succeed' by writing the artifact derive_band_products would write."""
    monkeypatch.setattr(rm, "ARTIFACT_ROOT", tmp_path)
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")

    def run_step(step, star, **k):
        calls.append(step["name"])
        if step["name"] == "derive_products":
            a = step["args"]
            out = tmp_path / a[a.index("--out") + 1]
            out.mkdir(parents=True, exist_ok=True)
            route = "SYNTH" if "--force-synthesis" in a else "PROFILEFIT"
            tag = ("_DEEPGRADED" if "--lines-deep-graded" in a else
                   "_" + {"reference": "REFERENCE", "graded": "GRADED"}[
                       a[a.index("--lines-tier") + 1]])
            stem = (f"{a[a.index('--element') + 1]}{a[a.index('--ion') + 1]}_4200_6908_"
                    f"{a[a.index('--instrument') + 1]}_{a[a.index('--holding') + 1]}_{route}"
                    f"{tag}")
            (out / f"{stem}_products.csv").write_text("treatment,A\n1D-LTE,7.5\n")
        return True, "ok"
    monkeypatch.setattr(rm, "_run_step", run_step)


def test_a_successful_run_that_published_nothing_is_not_DONE_but_is_recorded(
        rm, monkeypatch, tmp_path):
    """Publication is human-gated, so an unpublished build is not DONE -- but RYA-1233
    (Ryan 2026-10-02): it IS recorded, with what it was built from and what it wrote."""
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": []})
    reaches_the_executor(rm, monkeypatch)
    calls: list = []
    _fake_stage(rm, tmp_path, monkeypatch, calls)
    doc = rm.run("solar", "Si", ions=["I"], engines=["ts-lte"], bands=["VIS"],
                 methods=["synthesis"], pools=["codex"],
                 instruments=["solar_harps_molecfit_corrected"],
                 interpreter=sys.executable, ispec_dir="/x/ispec",
                 echo=False, report_dir=tmp_path)
    assert doc["counts"][rm.UNPUBLISHED] == 1 and doc["counts"][rm.DONE] == 0
    assert "publish_product.py" in doc["cells"][0]["reason"]
    led = json.loads((tmp_path / "inputs_hashes.json").read_text())
    (entry,) = led.values()
    assert entry["products"] == [] and entry["artifacts"], entry
    assert any(r["kind"] == "spectrum" for r in entry["inputs"]), "the spectrum is an input"


def test_an_unchanged_unpublished_build_is_not_redone(rm, monkeypatch, tmp_path):
    """THE point of RYA-1233's last change: nothing changed, so nothing re-runs."""
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": []})
    reaches_the_executor(rm, monkeypatch)
    calls: list = []
    _fake_stage(rm, tmp_path, monkeypatch, calls)
    kw = dict(ions=["I"], engines=["ts-lte"], bands=["VIS"], methods=["synthesis"], pools=["codex"],
              instruments=["solar_harps_molecfit_corrected"],
              interpreter=sys.executable, ispec_dir="/x/ispec", echo=False,
              report_dir=tmp_path)
    assert rm.run("solar", "Si", **kw)["counts"][rm.UNPUBLISHED] == 1
    calls.clear()
    second = rm.run("solar", "Si", **kw)
    assert calls == [], f"an unchanged unpublished build was redone: {calls}"
    assert second["counts"][rm.SKIP] == 1
    assert "NOT YET PUBLISHED" in second["cells"][0]["reason"]
    # ...and a deleted artifact is rebuilt rather than trusted
    for f in (tmp_path / "data").rglob("*_products.csv"):
        f.unlink()
    assert rm.run("solar", "Si", **kw)["counts"][rm.UNPUBLISHED] == 1


def test_the_spectrum_a_holding_reads_is_fingerprinted(rm):
    for hid in ("solar_kpno_molecfit_corrected", "solar_harps_molecfit_corrected",
                "solar_iag", "solar_kpno_kurucz2005_corrected"):
        files, why = rm.holding_source_files(hid)
        assert files, f"{hid}: {why}"


def test_the_loop_killer_a_second_run_does_zero_work(rm, monkeypatch, tmp_path):
    """THE acceptance property: run, then run again, and the second run executes nothing.

    Driven through `run()` twice rather than asserted on `is_current` alone, because
    the thing that can break is the WIRING -- a hash computed one way on write and
    another on read would pass every unit test and re-run the whole matrix forever.
    The stages are stubbed and the feed is stubbed; what is under test is the
    orchestrator's own bookkeeping, which is the only part of this that is ours.
    """
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")
    published = [{"element": "Si", "ion": "I", "band": "VIS", "instrument": "harps",
                  "holding": "solar_harps_molecfit_corrected", "treatment": "1D-LTE", "A": 7.51,
                  "route": "PROFILEFIT", "selector": "GRADED",
                  "n_lines": 12, "tier": "GRADED"}]
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": published})

    reaches_the_executor(rm, monkeypatch)

    calls = []
    monkeypatch.setattr(rm, "_run_step",
                        lambda step, *a, **k: (calls.append(step["name"]), (True, "ok"))[1])

    kw = dict(ions=["I"], engines=["ts-lte"], bands=["VIS"],
              instruments=["solar_harps_molecfit_corrected"],
              methods=["profile-fit"],
              interpreter=sys.executable, ispec_dir="/x/ispec", echo=False,
              report_dir=tmp_path)
    first = rm.run("solar", "Si", **kw)
    assert first["counts"][rm.DONE] == 1, first["cells"]
    assert calls, "the first run executed nothing, so this proves nothing"
    assert first["cells"][0]["A"] == 7.51 and first["cells"][0]["n_lines"] == 12

    calls.clear()
    second = rm.run("solar", "Si", **kw)
    assert second["counts"][rm.SKIP] == 1
    assert second["counts"][rm.DONE] == 0
    assert calls == [], f"the second run dispatched {calls} -- idempotency is broken"
    assert second["resume"] == [], "a current cell must not be listed for resume"


def test_a_deck_independent_step_is_built_once_per_run(rm, monkeypatch, tmp_path):
    """🔴 `measure_ew` ran 15 times to produce 3 files on the first full Si run.

    `descriptor.key` carries no deck, so all five decks over one (holding, band)
    name the same EW artifact. Twelve of those fifteen runs rewrote bytes that were
    already correct -- the repeat-work this ticket exists to end, inside a single
    invocation of the thing meant to end it.
    """
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": []})

    reaches_the_executor(rm, monkeypatch)
    ran = []
    monkeypatch.setattr(rm, "_run_step",
                        lambda step, *a, **k: (ran.append(step["name"]), (True, "ok"))[1])

    doc = rm.run("solar", "Fe", ions=["I"], bands=["VIS"], instruments=["solar_harps_molecfit_corrected"],
                 methods=["profile-fit"],
                 interpreter=sys.executable, ispec_dir="/x/ispec",
                 echo=False, report_dir=tmp_path)
    n_decks = len(rm.DECK_EMITS)
    assert doc["cells_total"] == n_decks, "expected one cell per deck"
    assert ran.count("measure_ew") == 1, \
        f"the deck-independent EW step ran {ran.count('measure_ew')} times, not once"
    # 🔴 And the deck-DEPENDENT step must still run for every deck. The resolver
    # declares the same `produces` path for all five, so a reuse rule keyed on that
    # string skipped four of them silently. This is the half that caught it.
    assert ran.count("derive_products") == n_decks, \
        f"derive_products ran {ran.count('derive_products')} times for {n_decks} decks"


def test_a_reused_step_is_only_reused_after_it_SUCCEEDED(rm, monkeypatch, tmp_path):
    """A failed step must not poison the run by looking built to every later cell."""
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": []})

    reaches_the_executor(rm, monkeypatch)
    ran = []
    monkeypatch.setattr(rm, "_run_step",
                        lambda step, *a, **k: (ran.append(step["name"]), (False, "boom"))[1])
    doc = rm.run("solar", "Fe", ions=["I"], bands=["VIS"], instruments=["solar_harps_molecfit_corrected"],
                 methods=["profile-fit"],
                 interpreter=sys.executable, ispec_dir="/x/ispec",
                 echo=False, report_dir=tmp_path)
    assert ran.count("measure_ew") == doc["cells_total"], \
        "a FAILED step was cached as though it had produced its artifact"


def test_a_moved_input_re_runs_the_cell(rm, monkeypatch, tmp_path):
    """The other half: SKIP must be earned, not sticky. A stage script edit re-runs it."""
    monkeypatch.setattr(rm, "LEDGER", tmp_path / "inputs_hashes.json")
    monkeypatch.setattr(rm, "load_feed", lambda *a, **k: {"products": [
        {"element": "Si", "ion": "I", "band": "VIS", "instrument": "harps",
         "holding": "solar_harps_molecfit_corrected", "treatment": "1D-LTE", "A": 7.51,
                  "route": "PROFILEFIT", "selector": "GRADED",
         "n_lines": 12, "tier": "GRADED"}]})

    reaches_the_executor(rm, monkeypatch)
    monkeypatch.setattr(rm, "_run_step", lambda *a, **k: (True, "ok"))

    kw = dict(ions=["I"], engines=["ts-lte"], bands=["VIS"],
              instruments=["solar_harps_molecfit_corrected"],
              methods=["profile-fit"],
              interpreter=sys.executable, ispec_dir="/x/ispec", echo=False,
              report_dir=tmp_path)
    assert rm.run("solar", "Si", **kw)["counts"][rm.DONE] == 1
    assert rm.run("solar", "Si", **kw)["counts"][rm.SKIP] == 1

    real = rm._file_fingerprint
    monkeypatch.setattr(rm, "_file_fingerprint",
                        lambda p: "EDITED" if p.name == "derive_band_products.py" else real(p))
    again = rm.run("solar", "Si", **kw)
    assert again["counts"][rm.SKIP] == 0
    assert again["counts"][rm.DONE] == 1
    assert "stage script moved" in again["cells"][0]["reason"] or \
        again["cells"][0]["status"] == rm.DONE


class _AllGo(dict):
    """Every (holding, band) lookup answers GO -- so cells reach the executor."""

    def __init__(self, row):
        super().__init__()
        self._row = row

    def get(self, _key, _default=None):
        return self._row


def reaches_the_executor(rm, monkeypatch):
    """Clear the two gates that stand between a cell and dispatch, for tests about
    what happens AFTER dispatch.

    🔴 BOTH stubs are load-bearing, and the second one was learned the hard way.

    `_readiness_index` is stubbed because RYA-1069's verdict is its own module's
    business and these tests are about the orchestrator's bookkeeping.

    `verify_numpy_ceiling` is stubbed because it asks the RUNNING INTERPRETER for
    its numpy version, and that made these tests depend on the machine. On this Mac
    numpy is 1.26.4 and they passed; on Sirius CI the interpreter carries numpy
    above the RYA-682 2.3 ceiling, so `run()` correctly BLOCKED every cell and six
    tests failed for a reason that has nothing to do with what they assert. The
    guard was right and the tests were wrong. The ceiling has its own tests
    (`test_the_numpy_ceiling_is_verified_on_the_pinned_interpreter` and
    `test_an_interpreter_that_cannot_run_is_a_refusal_not_a_crash`), which is where
    that behaviour belongs -- not as an invisible precondition of every other test.
    """
    class _Row:
        measurement_ready, blocking_gate, reader_wired = "GO", "", True
    monkeypatch.setattr(rm, "_readiness_index", lambda *a, **k: (_AllGo(_Row()), ""))
    monkeypatch.setattr(rm, "verify_numpy_ceiling",
                        lambda interpreter: (True, f"stubbed for {interpreter}"))
    # RYA-1233: the governing process (steps 1-9) has its own tests below; these tests
    # are about what happens after it has let a cell through.
    monkeypatch.setattr(rm, "process_steps", lambda *a, **k: [])
    monkeypatch.setattr(rm, "cell_process_hold", lambda *a, **k: "")


# ── the RYA-682 ceiling is CHECKED, not asserted ─────────────────────────────

def test_the_numpy_ceiling_is_verified_on_the_pinned_interpreter(rm):
    ok, why = rm.verify_numpy_ceiling(sys.executable)
    assert isinstance(ok, bool) and why


def test_an_interpreter_that_cannot_run_is_a_refusal_not_a_crash(rm):
    ok, why = rm.verify_numpy_ceiling("/nonexistent/python")
    assert not ok and "could not be run" in why


# ── back-compat: the whole-star path is untouched ────────────────────────────

def test_matrix_flags_without_element_are_refused(capsys):
    """`--dry-run` alone would have run the REAL whole-star chain while the operator
    believed nothing was executing."""
    import argparse
    import run_pipeline
    p = argparse.ArgumentParser()
    args = argparse.Namespace(element=None, ion=None, band=None,
                              instrument=None, engine=None, dry_run=True)
    with pytest.raises(SystemExit):
        run_pipeline._refuse_matrix_flags_without_element(p, args)


def test_the_whole_star_stages_are_still_the_same_six_in_the_same_order():
    """The lazy import must not have changed WHICH stages run or their order."""
    import inspect
    import run_pipeline
    src = inspect.getsource(run_pipeline._whole_star_stages)
    assert [n for n in ("spectra_normalize", "lines_fit", "params_stellar",
                        "abundances_derive", "uncertainty_stack", "ratios_interpret")
            if n in src] == ["spectra_normalize", "lines_fit", "params_stellar",
                             "abundances_derive", "uncertainty_stack", "ratios_interpret"]


def test_the_driver_imports_without_loading_any_science_stage():
    """🔴 It did not. `abundances_derive` evaluates `str | None` at class-body time,
    so on Python 3.9 every mode of the driver died on a stage it might never call."""
    import subprocess
    out = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); import run_pipeline; "
         "print(sorted(m for m in sys.modules if m.startswith('pipeline.')))" % str(ROOT)],
        capture_output=True, text=True, cwd=str(ROOT))
    assert out.returncode == 0, out.stderr[-2000:]
    assert "pipeline.abundances_derive" not in out.stdout


# ── RYA-1233: raw is never used; a corrected spectrum is frozen, not re-gated ──

def test_no_raw_ground_holding_is_ever_on_the_matrix(rm):
    """Ryan 2026-10-01/02: raw spectra are not used any more -- the corrected holding IS the data."""
    raw = set(rm.raw_holdings("solar"))
    assert {"solar_kpno", "solar_harps"} <= raw
    cells = {d.holding for d in rm.expand("solar", "Si")}
    assert not (cells & raw), f"raw holdings on the matrix: {sorted(cells & raw)}"


def test_every_corrected_holding_is_kept(rm):
    for hid in ("solar_kpno_molecfit_corrected", "solar_kpno_kurucz2005_corrected",
                "solar_harps_molecfit_corrected"):
        assert hid not in rm.raw_holdings("solar")


def test_a_corrected_holding_is_not_refused_on_a_stale_label(rm, monkeypatch, tmp_path):
    """Both fully-corrected KP holdings were NOT_READY on `evidence_state=audited` /
    a findings.md manifest while raw KP was GO. The element run asks only the reader."""
    class _Row:
        measurement_ready, blocking_gate, reader_wired = (
            "NO-GO", "evidence:audited;telluric:applied-unverified", True)
    monkeypatch.setattr(rm, "_readiness_index", lambda *a, **k: (_AllGo(_Row()), ""))
    monkeypatch.setattr(rm, "verify_numpy_ceiling", lambda i: (True, "stub"))
    doc = rm.run("solar", "Si", instruments=["solar_kpno_molecfit_corrected"],
                 engines=["ts-lte"], bands=["VIS"], dry_run=True,
                 interpreter=sys.executable, ispec_dir="/x/ispec", echo=False,
                 report_dir=tmp_path)
    assert doc["counts"][rm.NOT_READY] == 0, doc["cells"]
    assert "solar_kpno" in doc["excluded_raw_holdings"]


def test_band_policy_edges_are_synth_bands_edges():
    """RYA-1233: two edge tables disagreed (3780/3800, 9199/10000) and BLOCKED every IAG
    and Kurucz-2005 NIR cell. One table now."""
    from config.synth_bands import SYNTH_BANDS
    from pipeline import band_policy
    pol = {p.name: p for p in band_policy.POLICIES}
    for name in ("near-UV", "VIS", "red-optical"):
        assert pol[name].hi_A == SYNTH_BANDS[name].hi_A, name
    for name in ("VIS", "red-optical", "NIR"):
        assert pol[name].lo_A == SYNTH_BANDS[name].lo_A, name


def test_no_cell_is_blocked_by_a_band_edge_disagreement(rm):
    for d in rm.expand("solar", "Si"):
        assert not rm._band_of(d, "solar").disagreement, d


def test_every_deck_writes_its_own_products_file():
    """RYA-1233: five deck runs on one holding left ONE products file (each rewrote it)."""
    from pipeline.run_descriptor import RunDescriptor, deck_out_dir, BAND_PRODUCTS_DIR
    from pipeline.run_matrix import DECK_EMITS
    dirs = [deck_out_dir(d) for d in DECK_EMITS]
    assert len(set(dirs)) == len(dirs)
    assert deck_out_dir(RunDescriptor.__dataclass_fields__["engine_deck"].default) \
        == BAND_PRODUCTS_DIR, "the production deck's directory must not move"
    src = (ROOT / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")
    assert 'OUT = ROOT / ' + ' / '.join(f'"{p}"' for p in BAND_PRODUCTS_DIR.split("/")) in src


def test_the_resolver_passes_each_deck_its_own_out_dir():
    from pipeline.run_descriptor import RunDescriptor, resolve, deck_out_dir
    for deck in ("ts-lte", "gerber-1d-lte", "gerber-nlte"):
        d = RunDescriptor(element="Si", ion="I", instrument="harps",
                          holding="solar_harps_molecfit_corrected",
                          lo_A=3782.6, hi_A=6910.0, engine_deck=deck)
        step = [s for s in resolve(d, interpreter=sys.executable, ispec_dir="/x").steps
                if s["name"] == "derive_products"][0]
        assert step["args"][step["args"].index("--out") + 1] == deck_out_dir(deck)
        assert step["produces"].startswith(deck_out_dir(deck) + "/")


def test_no_two_decks_claim_the_same_treatment(rm):
    seen: dict = {}
    for deck, toks in rm.DECK_EMITS.items():
        for t in toks:
            assert t not in seen, f"{t} emitted by both {seen.get(t)} and {deck}"
            seen[t] = deck


def test_the_profilefit_route_labels_the_gerber_1d_deck_with_its_own_token():
    """It wrote "ENGINE-B" on ATLAS9 -- a byte-identical re-run of ts-lte under ts-lte's label."""
    src = (ROOT / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")
    main_route = src[src.index("if not a.skip_engine_b:\n        # RYA-1040"):]
    head = main_route[:main_route.index("RYA-880.")]
    assert 'elif a.engine_b_deck == "gerber-1d-lte":' in head
    assert "taxes.GERBER1D_LTE_MARCS.token" in head


def test_a_holding_with_no_reader_owns_no_cell_and_is_named(rm, tmp_path, monkeypatch):
    assert "elgueta2026_vizier" in rm.unwired_holdings("solar")
    assert "elgueta2026_vizier" not in {d.holding for d in rm.expand("solar", "Si")}


def test_an_unimportable_harness_excludes_nothing(rm, monkeypatch):
    p = rm.preflight()
    monkeypatch.setattr(p, "harness", lambda: None)
    assert rm.unwired_holdings("solar") == []


def test_a_band_that_permits_both_routes_gets_a_cell_for_each(rm):
    """RYA-1233: 150 of 160 live solar Fe products are SYNTH. The matrix ran only the
    policy's first method (profile-fit) in VIS / red-optical, so it could rebuild almost
    none of them -- and its synthesis leg only fit lines the EW gates had passed, which
    drops exactly the strong lines synthesis exists for."""
    cells = [d for d in rm.expand("solar", "Si", ions=["I"], engines=["ts-lte"], bands=["VIS"],
                                  instruments=["solar_harps_molecfit_corrected"])]
    assert sorted((d.method, d.pool) for d in cells) == [
        ("profile-fit", "codex"), ("synthesis", "codex"), ("synthesis", "deep"),
        ("synthesis", "reference"), ("synthesis", "set:SI_AGSS21")]
    near_uv = [d for d in rm.expand("solar", "Fe", engines=["ts-lte"], bands=["near-UV"])]
    assert near_uv and {d.method for d in near_uv} == {"synthesis"}


def test_a_synthesis_cell_in_a_profile_fit_band_forces_synthesis():
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element="Si", ion="I", instrument="harps",
                      holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                      method="synthesis")
    r = resolve(d, interpreter=sys.executable, ispec_dir="/x")
    assert [s["name"] for s in r.steps] == ["derive_products"], "no EW step on SYNTH"
    assert "--force-synthesis" in r.steps[0]["args"]
    assert r.steps[0]["produces"].endswith("_SYNTH_products.csv")


def test_a_gerber_deck_the_element_does_not_have_is_blocked_before_dispatch():
    """Every solar Si gerber cell FAILED inside the stage: gerber_nlte has decks for Al and
    Fe only. A registry fact is answered by the resolver, not discovered by a run."""
    from pipeline.run_descriptor import RunDescriptor, resolve
    from pipeline import gerber_nlte
    for deck in ("gerber-nlte", "gerber-1d-lte", "gerber-mean3d", "gerber-mean3d-lte"):
        d = RunDescriptor(element="Si", ion="I", instrument="harps",
                          holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                          engine_deck=deck, method="synthesis")
        r = resolve(d, interpreter=sys.executable, ispec_dir="/x")
        if (("Si@mean3D" if "mean3d" in deck else "Si") not in gerber_nlte.DECKS):
            assert r.blocked_reason and "no Gerber deck" in r.blocked_reason, deck
    fe = RunDescriptor(element="Fe", ion="I", instrument="harps",
                       holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                       engine_deck="gerber-nlte", method="synthesis")
    assert "no Gerber deck" not in (resolve(fe, interpreter=sys.executable,
                                            ispec_dir="/x").blocked_reason or "")


def test_a_synthesis_cell_is_clipped_to_its_list_not_failed():
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element="Si", ion="I", instrument="harps",
                      holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                      method="synthesis")
    args = resolve(d, interpreter=sys.executable, ispec_dir="/x").steps[0]["args"]
    assert "--clip-to-synthesis-list" in args
    src = (ROOT / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")
    assert '"--clip-to-synthesis-list"' in src


# ── RYA-1233: Ryan's governing process -- steps 1-9 before measurement ─────────

def test_an_element_without_its_literature_is_held_not_measured(rm, tmp_path, monkeypatch):
    """Step 5. Si was measured on 2026-10-02 with no litscan; under the process it waits."""
    monkeypatch.setattr(rm, "LITSCAN_DIR", tmp_path)          # no <El>.yaml anywhere
    monkeypatch.setattr(rm, "_run_step",
                        lambda *a, **k: pytest.fail("a held cell was dispatched"))
    doc = rm.run("solar", "Si", engines=["ts-lte"], bands=["VIS"], interpreter=sys.executable,
                 ispec_dir="/x", echo=False, report_dir=tmp_path)
    assert doc["counts"][rm.HELD] == doc["cells_total"] > 0
    assert all("PROCESS step 5" in c["reason"] for c in doc["cells"])
    step5 = [r for r in doc["process_steps"] if r["step"] == 5][0]
    assert step5["ok"] is False


def test_iron_first(rm, monkeypatch):
    """Step 9: a non-Fe element waits for a published Fe product."""
    monkeypatch.setattr(rm, "load_feed", lambda star, el: {"products": []})
    rows = {r["step"]: r for r in rm.process_steps("solar", "Si")}
    assert rows[9]["ok"] is False and "[Fe/H]" in rows[9]["evidence"]
    assert 9 not in {r["step"] for r in rm.process_steps("solar", "Fe")}


def test_a_band_with_no_lab_graded_line_is_held_at_step_7(rm):
    from pipeline.run_descriptor import RunDescriptor
    d = RunDescriptor(element="Fe", ion="II", instrument="crires_plus",
                      holding="solar_crires_plus_h_rya1094", lo_A=15007.11, hi_A=17493.69,
                      method="synthesis", pool="reference")
    assert "PROCESS step 7" in rm.cell_process_hold(d)


def test_no_cell_measures_the_ungraded_pool(rm):
    """Step 7: every dispatched cell names a graded pool; `--lines-tier all` never runs."""
    from pipeline.run_descriptor import resolve
    cells = rm.expand("solar", "Fe", engines=["ts-lte"], bands=["VIS"])
    assert cells and all(d.pool in ("reference", "codex", "deep") for d in cells)
    for d in cells[:4]:
        args = resolve(d, interpreter=sys.executable, ispec_dir="/x").steps[-1]["args"]
        assert "--lines-deep-graded" in args or (
            "--lines-tier" in args and args[args.index("--lines-tier") + 1] != "all")


def test_published_line_sets_become_synthesis_pools(rm):
    """RYA-1233 step 7: Asplund's set in VIS, Elgueta's / Bergemann's in the IR."""
    cells = rm.expand("solar", "Si", ions=["I"], engines=["ts-lte"])
    pools = {(d.band, d.pool) for d in cells if d.pool and d.pool.startswith("set:")}
    assert ("VIS", "set:SI_AGSS21") in pools
    assert ("NIR", "set:SI_ELGUETA2026") in pools and ("NIR", "set:SI_BERGEMANN2013") in pools
    assert all(d.method == "synthesis" for d in cells if d.pool.startswith("set:"))


def test_a_set_pool_dispatches_lines_from_set():
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element="Si", ion="I", instrument="harps",
                      holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                      method="synthesis", pool="set:SI_AGSS21")
    args = resolve(d, interpreter=sys.executable, ispec_dir="/x").steps[-1]["args"]
    i = args.index("--lines-from-set")
    assert args[i + 1] == "SI_AGSS21=data/reference/line_sets/si_agss21_SiI.csv"


def test_a_set_pool_reads_its_own_ions_file():
    """SI_AGSS21 is two registry rows (Si I, Si II). Keyed on the name alone, every Si II
    cell was handed the Si I file and refused "7 of 7 lines not in the synthesis list"."""
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element="Si", ion="II", instrument="harps",
                      holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                      method="synthesis", pool="set:SI_AGSS21")
    args = resolve(d, interpreter=sys.executable, ispec_dir="/x").steps[-1]["args"]
    assert args[args.index("--lines-from-set") + 1] == \
        "SI_AGSS21=data/reference/line_sets/si_agss21_SiII.csv"


def test_an_empty_pool_refusal_is_held_not_failed(rm, monkeypatch, tmp_path):
    reaches_the_executor(rm, monkeypatch)
    monkeypatch.setattr(rm, "_run_step", lambda *a, **k: (
        False, "derive_products: exit 1: no graded line in this band sits above the EW depth gate"))
    doc = rm.run("solar", "Si", ions=["I"], engines=["ts-lte"], bands=["VIS"],
                 methods=["synthesis"], pools=["deep"], interpreter=sys.executable,
                 ispec_dir="/x", echo=False, report_dir=tmp_path)
    assert doc["counts"][rm.FAILED] == 0 and doc["counts"][rm.HELD] >= 1


def test_a_set_file_is_part_of_its_cells_fingerprint(rm):
    """A corrected line set must re-run its cells (RYA-1233)."""
    from pipeline.run_descriptor import RunDescriptor, resolve
    d = RunDescriptor(element="Si", ion="I", instrument="harps",
                      holding="solar_harps_molecfit_corrected", lo_A=3782.6, hi_A=6910.0,
                      method="synthesis", pool="set:SI_AGSS21")
    r = resolve(d, interpreter=sys.executable, ispec_dir="/x")
    kinds = {row["kind"]: row["name"] for row in rm.input_fingerprints(d, r, manifest_path=None)}
    assert kinds.get("line_set") == "data/reference/line_sets/si_agss21_SiI.csv"


def test_every_agss21_si_line_matches_the_synthesis_wavelength():
    """Amarsi's Table 1 prints 0.01 A; three lines sat > 0.005 A from the list and were lost."""
    import csv
    rows = list(csv.DictReader(open(ROOT / "data/reference/line_sets/si_agss21_SiI.csv")))
    cg = {round(float(r["wavelength_air_A"]), 3) for r in
          csv.DictReader(open(ROOT / "data/linelists/canonical_gf.csv")) if r["species"] == "Si I"}
    assert len(rows) == 9
    assert all(round(float(r["wavelength_air_A"]), 3) in cg for r in rows)
