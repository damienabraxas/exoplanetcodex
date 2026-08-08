"""
tests/test_two_engine_reliable_contract_rya680_691.py
=====================================================
RYA-680 (Co I / Ba II wiring) + RYA-691 (the `reliable` contract), one file because
they are one function: `scripts/rya527_two_engine_run.py::_dedicated_engine_B`.

What these pin, and why each one exists
---------------------------------------
**RYA-680 — two ratified measurements the orchestrator never read.** RYA-673's Engine-B
wiring audit classified Co I and Ba II `NO_HARNESS_INVOCATION`: atom, grid and a merged,
ratified synthesis result all present, and no call site. The consequence was not a wrong
number but NO number — neither species entered the two-engine species set, so RYA-525's
loud-fail could not see them either (it iterates the union of the three coverage sources,
and a species in none of them never enters the loop). Gate 3 read UNEVALUABLE for both,
and no measurement quality could move it.

**The Ba trap, pinned deliberately.** Two Ba artifacts exist. `solar_ba_deblend_rya581.json`
holds A(Ba) 2.237 (in-window blend fit, merged PR #190). `solar_ba_synthesis_rya559.json`
holds 2.410 — superseded, inflated by a `blend_flag=True` pool EW that an EW inversion
cannot deblend. RYA-673's own audit map points at the RYA-559 file, so following the map
naively re-wires the stale value and silently undoes RYA-581. `test_ba_is_the_581_deblend`
and `test_driver_does_not_reference_the_superseded_559_artifact` exist so that cannot
happen again by accident.

**RYA-691 — the flag was consulted for two reads of eight.** Every artifact was assumed
reliable, including the one that carried an explicit flag (Sr II) six lines above a Zr
block whose comment reads "RELIABILITY-GATED throughout". Under either rchi2 ceiling
RYA-679 considered, Sr II 4077 would have been demoted — and this consumer would have
used it anyway, with the demotion recorded in the artifact and invisible downstream.

**The `or` fallback.** `v = m.get('A_nlte') or m.get('A_lte_median')` recorded an LTE
median under an NLTE provenance tag. It is live for V I (`nlte_void: true`), and `0.0`
being falsy made it wrong for a second reason.

The values are pinned too (`test_no_emitted_value_moved`). Both tickets forbid moving a
number: this was a contract fix and a wiring fix, and a changed value is a science
decision for Ryan.
"""
import ast
import importlib.util
import io
import json
import sys
import tokenize
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DRIVER = ROOT / 'scripts' / 'rya527_two_engine_run.py'


def _source_without_comments(path):
    """The driver's source with `#` comments removed.

    Needed because the defects these tests forbid are also DOCUMENTED in the driver —
    the Ba block explains in a comment exactly which artifact it refuses and why. A
    naive substring search cannot tell an explanation from a call site, and a test that
    forbids naming the problem would be a test against writing it down.
    """
    out = []
    for tok in tokenize.generate_tokens(io.StringIO(path.read_text()).readline):
        if tok.type != tokenize.COMMENT:
            out.append(tok.string)
    return '\n'.join(out)

#: Every dedicated Engine-B value the map emitted BEFORE RYA-680/691, verbatim, plus the
#: two species RYA-680 adds. Read as: the contract fix moved nothing, and the wiring fix
#: added exactly Co and Ba.
EXPECTED_MAP = {
    ('C', 'I'): 8.491,     # unchanged — nlte_cno cross-arm primary (RYA-491/237)
    ('O', 'I'): 8.730,     # unchanged
    ('Mn', 'I'): 5.466,    # unchanged — gold v3 PASS
    ('Cu', 'I'): 4.345,    # unchanged
    ('V', 'I'): 3.917,     # unchanged — and now labelled LTE, which it always was
    ('Sr', 'II'): 2.759,   # unchanged — the one flagged artifact that was ungated
    ('Co', 'I'): 4.960,    # RYA-680 NEW — RYA-564 median of 5 reliable red HFS lines
    ('Ba', 'II'): 2.237,   # RYA-680 NEW — RYA-581 deblend, NOT RYA-559's 2.410
}

#: Species deliberately absent: their reliability/concordance gate is shut today.
EXPECTED_ABSENT = {('Zr', 'II'), ('Mg', 'I')}


@pytest.fixture(scope='module')
def drv():
    spec = importlib.util.spec_from_file_location('rya527_two_engine_run', DRIVER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def dedicated(drv):
    return drv._dedicated_engine_B()


# ── RYA-680: Co and Ba are wired ─────────────────────────────────────────────

def test_co_and_ba_are_wired(dedicated):
    """The whole ticket in one assertion: both species now produce an Engine-B record."""
    for key in (('Co', 'I'), ('Ba', 'II')):
        assert key in dedicated, (
            f"{key} produces NO dedicated Engine-B record — this is the RYA-673 "
            f"NO_HARNESS_INVOCATION state RYA-680 removed. Gate 3 cannot see a value "
            f"the orchestrator never reads.")


def test_co_is_the_rya564_reliable_median(dedicated):
    value, source, basis = dedicated[('Co', 'I')]
    summary = json.loads((ROOT / 'data' / 'results' /
                          'co_synthesis_rya564.json').read_text())['_summary']
    # Single-sourced with phase_c's `_co_reclassify`, which reads the same field. Two
    # consumers computing their own Co from the same per-line table is the RYA-669 shape.
    assert value == pytest.approx(float(summary['A_Co'])), (
        "Co's Engine-B value must BE the harness's `_summary.A_Co`, not a re-derivation")
    assert 'RYA-564' in source and 'RYA-679 reliability-gated' in basis
    assert summary['n_reliable'] == 5


def test_ba_is_the_581_deblend(dedicated):
    """The single most important assertion in this file.

    2.237 is RYA-581's in-window deblend. 2.410 is RYA-559's EW->COG value, superseded
    and blend-inflated. Wiring 2.410 would silently undo a merged result.
    """
    value, source, basis = dedicated[('Ba', 'II')]
    deblend = json.loads((ROOT / 'data' / 'results' /
                          'solar_ba_deblend_rya581.json').read_text())
    assert deblend['ticket'] == 'RYA-581'
    assert value == pytest.approx(float(deblend['A_nlte']))
    assert value == pytest.approx(2.237)
    assert value != pytest.approx(2.410), "the SUPERSEDED RYA-559 value is wired"
    assert 'RYA-581' in source
    assert 'reliability-gated' in basis, "the RYA-581 artifact carries `reliable`; use it"


def test_driver_does_not_reference_the_superseded_559_artifact():
    """The driver must not know the RYA-559 path exists — nothing to fall back to.

    `solar_ba_synthesis_rya559.json` stays in the tree (its clean-EW cross-checks are the
    evidence for the supersession), so a future edit could reach for it as a "fallback"
    and silently restore 2.410 the first time the deblend artifact went missing.
    """
    assert 'solar_ba_synthesis_rya559' not in _source_without_comments(DRIVER), (
        "the two-engine driver references the SUPERSEDED RYA-559 Ba artifact")


def test_ba_route_refuses_a_non_581_artifact(drv, tmp_path, monkeypatch):
    """Positive control for the trap: hand the route the RYA-559 file and it must RAISE."""
    stale = ROOT / 'data' / 'results' / 'solar_ba_synthesis_rya559.json'
    monkeypatch.setattr(drv, 'BA_DEBLEND_JSON', stale)
    with pytest.raises(drv.DedicatedEngineBError) as exc:
        drv._dedicated_engine_B()
    assert 'RYA-581' in str(exc.value) and '2.410' in str(exc.value)


def test_co_route_raises_rather_than_dropping_the_element(drv, tmp_path, monkeypatch):
    """RYA-564's ratified rule: no reliable red line -> NO VALUE, and the blue-edge 3845
    artifact is never a fall-back. RYA-680 adds: nor is silence."""
    empty = tmp_path / 'co_empty.json'
    empty.write_text(json.dumps({'_summary': {'A_Co': None, 'n_reliable': 0,
                                              'reason': 'no line cleared the floor'}}))
    monkeypatch.setattr(drv, 'CO_JSON', empty)
    with pytest.raises(drv.DedicatedEngineBError) as exc:
        drv._dedicated_engine_B()
    assert 'no line cleared the floor' in str(exc.value)


def test_co_route_raises_when_the_summary_contradicts_its_own_per_line_flags(
        drv, tmp_path, monkeypatch):
    """Two views of one fact that can disagree is the RYA-669 defect shape."""
    co = json.loads((ROOT / 'data' / 'results' / 'co_synthesis_rya564.json').read_text())
    co['4813.476']['harps']['reliable'] = False       # summary still lists it
    forged = tmp_path / 'co_contradictory.json'
    forged.write_text(json.dumps(co))
    monkeypatch.setattr(drv, 'CO_JSON', forged)
    with pytest.raises(drv.DedicatedEngineBError) as exc:
        drv._dedicated_engine_B()
    assert 'contradicts itself' in str(exc.value).lower()


# ── RYA-691: the `reliable` contract ─────────────────────────────────────────

@pytest.mark.parametrize('attr,mutate', [
    ('BA_DEBLEND_JSON', lambda d: d.__setitem__('reliable', False)),
    ('SR2_JSON', lambda d: d['4077.709']['harps'].__setitem__('reliable', False)),
])
def test_unreliable_value_cannot_reach_the_map(drv, tmp_path, monkeypatch, attr, mutate):
    """THE regression test for RYA-691.

    Before the fix these two reads took the value with no reference to the flag, so an
    artifact that had recorded its own demotion flowed straight into the map. Now the
    read raises. Demonstrated failing on the pre-fix implementation (both returned the
    demoted value silently) and passing here.
    """
    data = json.loads(getattr(drv, attr).read_text())
    mutate(data)
    forged = tmp_path / f'{attr.lower()}_unreliable.json'
    forged.write_text(json.dumps(data))
    monkeypatch.setattr(drv, attr, forged)
    with pytest.raises(drv.DedicatedEngineBError) as exc:
        drv._dedicated_engine_B()
    msg = str(exc.value)
    assert 'reliable=False' in msg and 'RYA-679' in msg


def test_every_read_records_a_reliability_basis(dedicated):
    """No read may be silent about what cleared it. Two legitimate answers — gated by the
    RYA-679 flag, or explicitly UNGATED with the reason — and no third."""
    for key, (value, source, basis) in dedicated.items():
        assert basis, f"{key} records no reliability basis"
        assert basis.startswith('UNGATED — ') or 'RYA-679 reliability-gated' in basis, (
            f"{key} records an uninterpretable reliability basis: {basis!r}")
        assert source, f"{key} records no Engine-B source"


def test_cno_is_documented_as_different_not_forced_under_the_flag(dedicated):
    """RYA-691 §3A explicitly permits C/O to differ, and they do.

    `role == 'primary'` selects WHICH indicator speaks for the element; it says nothing
    about whether the fit was trustworthy. The RYA-491/237 cross-arm artifact carries no
    `reliable` key on any indicator, and cannot: it is a multi-indicator cross-arm
    reconciliation, not an in-window profile fit, so it has neither `dEW_dA` nor `railed`
    and the RYA-679 rule is not computable over it. Fabricating a check would assert a
    gate that does not exist — so the artifact's own quality statement is recorded
    instead.
    """
    cross_arm = json.loads((ROOT / 'data' / 'audit' / 'cno_synthesis' /
                            'solar_phase_a_cross_arm.json').read_text())['cross_arm']
    for el in ('C', 'O'):
        for ind in cross_arm[el]['indicators']:
            assert 'reliable' not in ind, (
                f"the CNO artifact now carries a reliability flag for {el} — the "
                f"documented exemption no longer holds and the read must gate on it")
        basis = dedicated[(el, 'I')][2]
        assert basis.startswith('UNGATED — ')
        assert 'role=' in basis and 'verdict=' in basis


def test_o_records_that_it_takes_the_first_of_two_primaries(dedicated):
    """O has TWO `primary` indicators (OI_6300 8.73, OI_777 8.74) and the read has always
    taken whichever came first in the artifact. That is a real selection decision made by
    file order; it is not changed here (changing it would move O's value), but it is no
    longer invisible."""
    assert 'FIRST of 2 primary indicators' in dedicated[('O', 'I')][2]


# ── RYA-691 §2/§3C: the `or` fallback ────────────────────────────────────────

def test_or_fallback_is_gone_from_the_source():
    """Structural, not textual: an `or` whose left operand fetches an abundance key.

    Matched on the AST so the driver can still DESCRIBE the retired construct in prose —
    which it must, because the reason it was wrong is not obvious from the fix.
    """
    ABUND_KEYS = {'A_nlte', 'A_lte_median', 'A_LTE', 'A_NLTE', 'A_Co', 'A'}

    def _fetches_abundance(node):
        if not isinstance(node, ast.Call):
            return False
        f = node.func
        if not (isinstance(f, ast.Attribute) and f.attr == 'get'):
            return False
        return bool(node.args) and isinstance(node.args[0], ast.Constant) \
            and node.args[0].value in ABUND_KEYS

    offenders = [
        n.lineno for n in ast.walk(ast.parse(DRIVER.read_text()))
        if isinstance(n, ast.BoolOp) and isinstance(n.op, ast.Or)
        and any(_fetches_abundance(v) for v in n.values)]
    assert not offenders, (
        f"the silent NLTE->LTE `or` substitution is back at line(s) {offenders}: it "
        f"records an LTE median under an NLTE provenance tag with no trace the engine "
        f"changed, and `or` treats a genuine 0.0 as missing")


def test_nlte_engine_is_named_not_assumed(drv):
    nlte, tag = drv._nlte_or_lte({'A_nlte': 5.466, 'A_lte_median': 5.442}, 'x', 'T')
    assert (nlte, tag) == (5.466, 'T (NLTE)')
    lte, tag = drv._nlte_or_lte({'A_lte_median': 3.917, 'nlte_void': True}, 'x', 'T')
    assert lte == 3.917 and 'LTE FALL-BACK' in tag and 'nlte_void' in tag


def test_zero_is_not_treated_as_missing(drv):
    """`0.0` is falsy, so `or` would have substituted the LTE median for it.

    NOT reachable with today's inputs — on the A(X) = 12 + log(N_X/N_H) scale a literal
    0.0 means N_X/N_H = 1e-12, nine orders below the rarest species measured here, and no
    artifact in the repo carries it. Fixed regardless: "unreachable today" is a property
    of the data, not of the code, and this is what `is None` buys.
    """
    value, tag = drv._nlte_or_lte({'A_nlte': 0.0, 'A_lte_median': 5.0}, 'x', 'T')
    assert value == 0.0 and tag == 'T (NLTE)'


def test_v_is_labelled_lte_because_it_always_was(dedicated):
    """The live instance of the silent substitution. V I has no NLTE grid
    (`nlte_void: true`, RYA-466), so its emitted value has always been the LTE median
    wearing the same tag as Cu's genuine NLTE value. The number does not move; the label
    stops lying."""
    value, source, _ = dedicated[('V', 'I')]
    assert value == pytest.approx(3.917)
    assert 'LTE FALL-BACK' in source and 'nlte_void' in source
    assert '(NLTE)' in dedicated[('Cu', 'I')][1]


# ── the STOP condition: no value may move ────────────────────────────────────

def test_no_emitted_value_moved(dedicated):
    got = {k: round(v[0], 3) for k, v in dedicated.items()}
    assert got == {k: round(v, 3) for k, v in EXPECTED_MAP.items()}


def test_gated_shut_species_stay_out(dedicated):
    """Zr II (no line over the dEW/dA floor) and Mg I 5528 (reliable but DISCORDANT with
    5711 by 0.21-0.23 dex) are held out on purpose. Wiring is not adoption."""
    for key in EXPECTED_ABSENT:
        assert key not in dedicated


def test_record_carries_the_engine_b_provenance():
    """`src` used to be unpacked into a local and dropped, so the emitted record said
    nothing about which harness produced an Engine-B value or what cleared it."""
    src = DRIVER.read_text()
    assert 'engineB_source=b_source' in src and 'engineB_reliability=b_reliability' in src
