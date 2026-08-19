"""
tests/test_preflight_check_rya905.py
====================================
RYA-905 — the pre-flight readiness reconciliation.

WHAT THESE TESTS ARE FOR. The module under test emits NEGATIVES ("this holding is not
reachable", "this line is not grid-served"), and a negative is unfalsifiable without a
control at BOTH ends. So the central test here is not "does it warn about HARPS" — a
check that warns about everything would pass that. It is: **the same input with the
defect removed must stop warning**, and an instrument we genuinely do not hold must
produce INFO rather than WARN, because a survey lacking an arm is the normal state and a
check that cannot tell the two apart is the thing RYA-905 exists to replace.

Nothing here asserts a magic constant (RYA-845): the assertions are on RELATIONSHIPS —
severity flips, one gap counted once, exit status advisory — not on a WARN count that
would pin today's registry contents into a test.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts import preflight_check as pf   # noqa: E402


# ── The dispatch reader and its controls ─────────────────────────────────────

def test_dispatch_reader_controls_pass_on_the_real_harness():
    """The reader must SEE the arm we know is wired, and not see one that cannot exist."""
    d = pf.read_dispatch()
    assert d.controls_ok, d.control_note
    assert "kpno_solar_atlas" in d.instruments
    assert pf.DISPATCH_SENTINEL not in d.instruments


def test_dispatch_reader_reports_control_failure_rather_than_a_clean_report(tmp_path):
    """A reader that cannot parse the dispatch must say so, not report zero arms.

    This is the failure mode that would be most dangerous in production: silently
    returning "no instruments dispatched" turns every holding into a WARN, and a page of
    false WARNs is how a check gets ignored. The controls exist to make that state
    LOUD (ERROR) instead of merely wrong.
    """
    fake = tmp_path / "measure_band_ew.py"
    fake.write_text("def load_window(instrument, centre, pad):\n    return None\n")
    d = pf.read_dispatch(fake)
    assert not d.controls_ok
    assert "CONTROL FAILED" in d.control_note


def _harness_with(instruments) -> str:
    branches = "\n".join(
        f'    if instrument == "{i}":\n        return 1, 2, 3' for i in instruments)
    return (
        "PRE_NORMALISED = {" + ", ".join(f'"{i}": True' for i in instruments) + "}\n"
        '_LOADER_HOLDING = {"kpno_solar_atlas": "solar_kpno"}\n'
        "def load_window(instrument, centre, pad):\n"
        + branches
        + "\n    raise LookupError('no window loader')\n")


def test_dispatch_reader_discriminates_wired_from_unwired(tmp_path):
    """BOTH ENDS. The same reader, two sources, opposite answers about `harps`."""
    without = tmp_path / "without.py"
    without.write_text(_harness_with(["kpno_solar_atlas", "iag_fts_solar_atlas"]))
    with_harps = tmp_path / "with.py"
    with_harps.write_text(
        _harness_with(["kpno_solar_atlas", "iag_fts_solar_atlas", "harps"]))

    assert "harps" not in pf.read_dispatch(without).instruments
    assert "harps" in pf.read_dispatch(with_harps).instruments
    assert pf.read_dispatch(without).controls_ok
    assert pf.read_dispatch(with_harps).controls_ok


# ── Check 1: the RYA-897 retro-catch, and its control ────────────────────────

@pytest.fixture(scope="module")
def solar_fe():
    return pf.run("sun", "Fe", "I", [pf.ROOT])


def _findings(results, number):
    return [f for r in results if r.number == number for f in r.findings]


def test_check1_no_longer_warns_on_harps_because_the_loader_now_EXISTS(solar_fe):
    """🔴 THIS TEST WAS INVERTED BY RYA-911, and that is the correct outcome.

    It was RYA-905's acceptance proof: `solar_harps` is verified and `load_window` has no
    branch for it — the RYA-897 defect, LIVE ON MAIN at the time. RYA-911 wires the HARPS
    loader (porting RYA-897's build onto the RYA-904 holding table), so the premise is
    now false and the original assertion would be pinning a defect in place.

    An acceptance proof for a defect is TIME-LIMITED BY CONSTRUCTION. Keeping it green by
    keeping the defect is the failure mode; the honest move is to invert it and say why.
    Its counterpart, `test_check1_control_the_warning_goes_away_when_the_loader_exists`,
    always tested this end synthetically — this now tests it for real.
    """
    st, results = solar_fe
    warns = [f for f in _findings(results, 1)
             if f.severity == pf.WARN and f.subject == "solar_harps"]
    assert not warns, (
        "check 1 still WARNs that the harness cannot load solar_harps, but RYA-911 "
        f"wired it: {[f.message for f in warns]}")
    oks = [f for f in _findings(results, 1)
           if f.severity == pf.OK and f.subject == "solar_harps"]
    assert oks, "solar_harps is wired but check 1 reports it neither as OK nor as a WARN"
    assert "harps" in oks[0].message


def test_check1_still_discriminates_on_a_harness_that_LACKS_the_loader(monkeypatch,
                                                                      tmp_path):
    """THE OTHER END, kept. Inverting the test above would be worthless if check 1 had
    simply stopped WARNing about anything. Point it at a harness with no harps branch and
    the WARN must come back."""
    patched = tmp_path / "no_harps_harness.py"
    patched.write_text(
        'PRE_NORMALISED = {"kpno_solar_atlas": True, "harps": False}\n'
        '_LOADER_HOLDING = {"kpno_solar_atlas": "solar_kpno"}\n'
        "def load_window(instrument, centre, pad):\n"
        '    if instrument == "kpno_solar_atlas":\n        return 1, 2, 3\n'
        "    raise LookupError('no window loader')\n")
    monkeypatch.setattr(pf, "HARNESS_PY", patched)
    _, results = pf.run("sun", "Fe", "I", [pf.ROOT])
    warns = [f for f in _findings(results, 1)
             if f.severity == pf.WARN and f.subject == "solar_harps"]
    assert warns, "check 1 stopped detecting an unreachable holding altogether"
    assert warns[0].suggested_ticket, "a WARN must carry a suggested-ticket stub"


def test_check1_control_the_warning_goes_away_when_the_loader_exists(solar_fe, monkeypatch,
                                                                     tmp_path):
    """THE OTHER END. Give the harness a `harps` branch and the WARN must disappear.

    Without this, "check 1 warns about HARPS" is compatible with "check 1 warns about
    everything", and the check would be worthless the day the loader lands.
    """
    st, _ = solar_fe
    patched = tmp_path / "patched_harness.py"
    patched.write_text(_harness_with(
        list(st.dispatch.instruments) + ["harps"]))
    monkeypatch.setattr(pf, "HARNESS_PY", patched)

    _, results = pf.run("sun", "Fe", "I", [pf.ROOT])
    warns = [f for f in _findings(results, 1)
             if f.severity == pf.WARN and f.subject == "solar_harps"]
    assert not warns, "the WARN survived a harness that CAN load harps — it is not "\
                      "responding to reachability at all"


def test_check1_an_instrument_we_do_not_hold_is_INFO_not_WARN(solar_fe):
    """CRITICAL per the ticket: a genuinely-absent instrument is not a showstopper.

    No UVES for the Sun is the normal state of a survey. If that produced a WARN, every
    element would emit twenty of them and the report would be noise.
    """
    st, results = solar_fe
    infos = [f for f in _findings(results, 1) if f.severity == pf.INFO]
    not_held = [f for f in infos if f.subject == "not held"]
    assert not_held, "instruments with no holding were not reported at all"
    assert "uves" in not_held[0].message
    for f in _findings(results, 1):
        if f.severity == pf.WARN:
            assert f.subject in set(st.holdings.holding_id.astype(str)), (
                f"{f.subject} WARNed but is not a holding we have — an expected absence "
                f"was reported as a silent gap")


# ── Severity discipline ──────────────────────────────────────────────────────

@pytest.mark.parametrize("element", ["Fe", "Al"])
def test_every_warn_carries_a_discriminator_and_a_ticket_stub(element):
    """A WARN whose author cannot say why it is not an expected-absence has not made the
    distinction this module exists to make."""
    _, results = pf.run("sun", element, "I", [pf.ROOT])
    for r in results:
        for f in r.findings:
            if f.severity == pf.WARN:
                assert f.discriminator.strip(), f"{f.subject}: WARN with no discriminator"
                assert f.suggested_ticket, f"{f.subject}: WARN with no ticket stub"


@pytest.mark.parametrize("element", ["Fe", "Al"])
def test_all_six_checks_run_and_none_is_silent(element):
    _, results = pf.run("sun", element, "I", [pf.ROOT])
    assert [r.number for r in results] == [1, 2, 3, 4, 5, 6]
    for r in results:
        assert r.findings, f"check {r.number} ({r.name}) emitted nothing at all"


def test_it_is_advisory_and_exits_zero_even_with_warnings(capsys):
    """Advisory by construction. A survey lacking an arm must never fail a run."""
    rc = pf.main(["--star", "sun", "--element", "Fe"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "WARN=" in out
    assert "advisory" in out


def test_json_output_round_trips(tmp_path):
    out = tmp_path / "report.json"
    assert pf.main(["--star", "sun", "--element", "Al", "--json", str(out)]) == 0
    payload = json.loads(out.read_text())
    assert payload["ticket"] == "RYA-905"
    assert payload["star"] == "solar"
    assert {c["number"] for c in payload["checks"]} == {1, 2, 3, 4, 5, 6}


def test_an_unknown_star_is_refused_rather_than_reported_empty():
    """A star id nobody registered must not come back as 'no holdings'.

    That is a MANUFACTURED absence (RYA-833): the report would look like a finding about
    the survey when it is a typo in the argument.
    """
    with pytest.raises(SystemExit):
        pf.run("betelgeuse", "Fe", "I", [pf.ROOT])


# ── Check 3: the ion key, which is where a false absence came from ───────────

def test_grid_ion_key_is_normalised_not_string_compared():
    """The Amarsi-2020 Al grid spells the neutral stage `1`; the register spells it `I`.

    A string comparison between the two returns an empty selection and the check reports
    "the grid serves no line of this ion" — a false negative that looks exactly like a
    real coverage gap. The first cut of this module did precisely that.
    """
    _, results = pf.run("sun", "Al", "I", [pf.ROOT])
    c3 = _findings(results, 3)
    assert c3, "check 3 emitted nothing for Al"
    assert not any("serve no Al I wavelength" in f.message for f in c3), (
        "check 3 still reads the Al grid's numeric ion column as covering nothing")
    served = [f for f in c3 if f.severity == pf.OK and "served by the grid" in f.message]
    assert served, "check 3 found no grid-served Al lines at all"


def test_check3_reports_the_join_delta_it_used():
    """A wavelength join across two catalogues must show its own slack, or the match is
    an assertion rather than a measurement (RYA-871: lines keyed on wavelength alone)."""
    _, results = pf.run("sun", "Al", "I", [pf.ROOT])
    served = [f for f in _findings(results, 3) if "join delta" in f.message]
    assert served, "check 3 matched lines to a grid without reporting the delta"


# ── One gap is counted once ──────────────────────────────────────────────────

def test_a_single_gap_is_not_counted_by_two_checks(solar_fe):
    """One defect, one WARN — or the suggested-ticket list becomes three tickets per fix.

    ⚠️ RE-SUBJECTED BY RYA-911, and then CORRECTED. This was written on `solar_harps`,
    which was unreachable (check 1) AND absent from the product (check 5). RYA-911 wired
    HARPS, so that holding stopped being an example of anything.

    🔴 My first re-subjecting asserted the rule over EVERY holding that WARNs — "no
    subject appears under two check numbers" — and that is a DIFFERENT, FALSE rule. CI
    caught it: `solar_crires_plus_y_rya794` is WARNed by check 5 (it is reachable,
    telluric-clear, in wavelength reach, and still appears in no rendered product) and by
    check 6 (its `telluric_applied=applied` disagrees with what the catalogue says the
    arm's basis is). Those are TWO DEFECTS about one holding, not one defect counted
    twice, and fixing either leaves the other standing. Generalising from "subject" to
    "defect" was the error — a subject is not a defect.

    The real invariant, and the one the original was reaching for: **check 5 must never
    fire on a holding check 1 has already declared unreachable.** "It appears in no
    product" is not an independent finding about something nothing can read — that is one
    gap wearing two tickets. Checks 5 and 6 report genuinely separate properties and may
    co-fire; checks 1 and 5 may not.
    """
    _, results = solar_fe
    def _warned_by(number: int) -> set[str]:
        return {f.subject for r in results if r.number == number
                for f in r.findings if f.severity == pf.WARN}
    unreachable, absent = _warned_by(1), _warned_by(5)
    both = unreachable & absent
    assert not both, (
        f"one gap counted twice: {sorted(both)} are WARNed as unreachable (check 1) AND "
        f"as absent from the product (check 5). An unreachable holding cannot also owe a "
        f"product — that is one defect billed to two checks.")


def test_check5_names_the_check_that_owns_each_accounted_absence(solar_fe):
    _, results = solar_fe
    infos = [f for f in _findings(results, 5) if f.severity == pf.INFO]
    accounted = [f for f in infos if "accounted for" in f.message]
    assert accounted, "check 5 accounted for nothing — it cannot be discriminating"
    for f in accounted:
        assert "check " in f.message, f"{f.subject}: absence accounted without a pointer"


# ── Read-only ────────────────────────────────────────────────────────────────

def test_preflight_writes_nothing_it_was_not_asked_to_write(tmp_path, monkeypatch):
    """Read-only reconciliation. It must not touch intake artifacts or the harness."""
    watched = [pf.HOLDINGS_CSV, pf.INSTRUMENTS_CSV, pf.HARNESS_PY]
    before = {p: p.read_bytes() for p in watched}
    pf.run("sun", "Fe", "I", [pf.ROOT])
    for p in watched:
        assert p.read_bytes() == before[p], f"pre-flight modified {p}"
