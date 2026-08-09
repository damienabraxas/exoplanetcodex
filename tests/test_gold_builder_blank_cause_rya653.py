"""
tests/test_gold_builder_blank_cause_rya653.py
=============================================
RYA-653 — the RYA-596 blank-cause tripwire, now SHARED, applied to the gold
reference builder.

THE BUG THIS LOCKS OUT
----------------------
`scripts/build_solar_reference_v2_rya522.py::_scale_and_note` carried the same
unchecked fallthrough RYA-596 killed in the verdict: any row without a delta got
stamped "no independent-gf line survives the graded cull", a cull the builder
never checked. It is lying in frozen gold v2's **Ba** row — against a Ba that
RYA-559 has since measured (A(Ba) 2.410, Ba II 5853 HFS + Korotin2015 NLTE) —
and, rebuilt against today's verdict, it would have fabricated the SAME claim on
five more rows (Ca/Ti/Ni/Na/Al, carrying 2/10/2/2/1 graded survivors): RYA-596's
exact five, re-entering through the other file.

THE INVARIANT (shared, one module, both call sites):
    no stage emits or inherits a causal/provenance claim it did not itself
    establish.
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from pipeline import provenance_honesty as H  # noqa: E402
import build_solar_reference_v2_rya522 as B   # noqa: E402
import phase_c_verdict_rya371 as P            # noqa: E402

VERDICT = ROOT / 'data/audit/cno_synthesis/solar_phase_c_verdict.json'
CANDIDATE = ROOT / 'data/reference/solar/solar_abundances_corrected_candidate_rya653.csv'


# ── A. the module is genuinely shared — not a second copy ────────────────────

def test_both_call_sites_use_the_one_tripwire():
    """A forked copy would be the very duplication the invariant forbids, so
    identity — not similarity — is what this asserts."""
    assert P._assert_blank_cause_is_honest is H.assert_blank_cause_is_honest
    assert B.assert_blank_cause_is_honest is H.assert_blank_cause_is_honest
    assert P.ZERO_SURVIVOR_CHANNEL is H.ZERO_SURVIVOR_CHANNEL
    src = (ROOT / 'scripts/build_solar_reference_v2_rya522.py').read_text()
    assert H.ZERO_SURVIVOR_CLAIM not in src, (
        "the gold builder re-hardcodes the zero-survivor claim — it must come "
        "from pipeline.provenance_honesty")


def test_claim_is_recognized_in_both_spellings():
    """The verdict writes the claim behind an 'EW present; ' prefix, the gold
    builder wrote it bare. One claim, two spellings — both must be caught."""
    assert H.claims_zero_survivors(H.ZERO_SURVIVOR_CLAIM)
    assert H.claims_zero_survivors(H.ZERO_SURVIVOR_CHANNEL)
    assert H.claims_zero_survivors(f"{H.ZERO_SURVIVOR_CLAIM} — A(X) 2.410 HELD")
    assert not H.claims_zero_survivors("EW: 2 curated line(s); value HELD")
    assert not H.claims_zero_survivors(None)


# ── B. the tripwire itself ───────────────────────────────────────────────────

def test_tripwire_raises_on_survivors():
    with pytest.raises(AssertionError, match='RYA-596'):
        H.assert_blank_cause_is_honest('Ca', H.ZERO_SURVIVOR_CHANNEL, n_lines=2)


def test_tripwire_raises_on_a_measured_value():
    """The Ba shape: a measured A(X) makes the zero-survivor claim unrepresentable
    even if the row's line count is blank."""
    with pytest.raises(AssertionError, match='2.410'):
        H.assert_blank_cause_is_honest('Ba', H.ZERO_SURVIVOR_CLAIM, n_lines=0,
                                       a_measured=2.410)
    with pytest.raises(AssertionError, match='RYA-596'):
        H.assert_blank_cause_is_honest('Ba', H.ZERO_SURVIVOR_CLAIM, n_lines=1,
                                       a_measured=2.410)


def test_tripwire_is_quiet_on_the_legitimate_shapes():
    """Must NOT over-correct: Mg/Y/Zr/Eu really do end the graded cull empty."""
    H.assert_blank_cause_is_honest('Mg', H.ZERO_SURVIVOR_CHANNEL, n_lines=0)
    H.assert_blank_cause_is_honest('Y', H.ZERO_SURVIVOR_CHANNEL, n_lines=None)
    H.assert_blank_cause_is_honest('Zr', H.ZERO_SURVIVOR_CHANNEL, n_lines=float('nan'),
                                   a_measured=float('nan'))
    H.assert_blank_cause_is_honest('Ca', "EW: 2 curated line(s); value HELD", n_lines=2,
                                   a_measured=6.324)


# ── C. the builder, end to end ───────────────────────────────────────────────

def _verdict_doc():
    return json.loads(VERDICT.read_text())


def _write_verdict(tmp_path, doc):
    p = tmp_path / 'verdict.json'
    p.write_text(json.dumps(doc))
    return p


def _build(tmp_path, verdict_path):
    return B.main(['--verdict', str(verdict_path),
                   '--out-csv', str(tmp_path / 'cand.csv'),
                   '--out-md', str(tmp_path / 'cand.md')])


def test_builder_raises_on_a_ba_like_phantom_row(tmp_path):
    """THE SMOKE TEST. Feed the builder the Ba shape that gold v2 actually froze —
    a zero-survivor cause on a row with survivors AND a measured value — and it
    must refuse to write the reference rather than immortalise the claim."""
    doc = _verdict_doc()
    for r in doc['verdicts']:
        if r['element'] == 'Ba':
            r['channel'] = H.ZERO_SURVIVOR_CHANNEL   # the phantom, as frozen in v2
            assert r['A_measured'] is not None and r['n_lines'] > 0
    with pytest.raises(AssertionError, match='RYA-596'):
        _build(tmp_path, _write_verdict(tmp_path, doc))
    assert not (tmp_path / 'cand.csv').exists(), (
        "the builder wrote a reference carrying a phantom cause — it must raise first")


def test_builder_raises_on_a_held_row_blamed_on_the_cull(tmp_path):
    """The five RYA-596 rows: blank + survivors is a RYA-522 tier HOLD, never a cull."""
    doc = _verdict_doc()
    for r in doc['verdicts']:
        if r['element'] == 'Ti':
            r['channel'] = H.ZERO_SURVIVOR_CHANNEL
            assert r['A_measured'] is None and r['n_lines'] == 10
    with pytest.raises(AssertionError, match='n_lines=10'):
        _build(tmp_path, _write_verdict(tmp_path, doc))


def test_rebuild_from_the_live_verdict_is_honest(tmp_path):
    """Rebuilt against the live verdict, NO row may claim a cull it cannot support.
    Pre-fix this build did not merely lie about Ba — it fabricated the claim for
    Ca/Ti/Ni/Na/Al too, and crashed outright on Sr."""
    _build(tmp_path, VERDICT)
    df = pd.read_csv(tmp_path / 'cand.csv')
    for _, r in df.iterrows():
        if H.claims_zero_survivors(str(r['note'])):
            assert int(r['n_lines']) == 0, f"{r['element']}: cull claim with survivors"
            assert pd.isna(r['A_X']) and pd.isna(r['A_X_nlte']), \
                f"{r['element']}: cull claim on a row carrying a value"


#: The ticket whose Ba measurement the gold row must cite. RYA-653 wrote RYA-559 (the
#: EW->COG value 2.410). RYA-581 then re-measured Ba II 5853 with an in-window blend fit
#: and got 2.237, because the pool EW carries blend_flag=True and an EW inversion cannot
#: deblend — it charged ~10 mA of neighbouring absorption to Ba. RYA-680 pointed both the
#: two-engine floor and phase_c at the deblend, so the gold row follows.
#:
#: Pinned as a constant rather than relaxed to "cites some ticket": the point of this
#: test is that the row names the measurement it actually came from, and that property
#: survives the supersession. If Ba is re-measured again this constant moves again.
BA_MEASUREMENT_TICKET = 'RYA-581'


def test_ba_row_is_sourced_from_its_measurement_ticket(tmp_path):
    """B. The corrected Ba row states its measurement and cites it — and the
    value is READ from the verdict artifact, never typed in."""
    a_ba = {r['element']: r for r in _verdict_doc()['verdicts']}['Ba']['A_measured']
    _build(tmp_path, VERDICT)
    row = pd.read_csv(tmp_path / 'cand.csv').set_index('element').loc['Ba']
    note = str(row['note'])
    assert not H.claims_zero_survivors(note), "Ba still blames the graded cull"
    assert BA_MEASUREMENT_TICKET in note, f"Ba row does not cite its measurement: {note!r}"
    assert f"{a_ba:.3f}" in note, f"Ba row does not carry the measured value: {note!r}"
    # The ratified `owed` tier still freezes NO value (RYA-522) — held, not
    # immortalised. Promoting it is a re-ratification, not this ticket's call.
    assert row['confidence'] == 'owed'
    assert pd.isna(row['A_X']) and pd.isna(row['A_X_nlte'])
    assert int(row['n_lines']) == 1


def test_synthesis_rows_never_borrow_another_tickets_citation(tmp_path):
    """Ba's note used to read 'HFS-resolved synthesis (RYA-411/466/473)' — the Mn/Cu/V
    tickets, a provenance claim this stage never established for Ba."""
    _build(tmp_path, VERDICT)
    rows = pd.read_csv(tmp_path / 'cand.csv').set_index('element')
    assert 'RYA-411' not in str(rows.loc['Ba', 'note'])
    assert 'RYA-466' not in str(rows.loc['Ba', 'note'])
    # the map still covers the rows it was written for
    assert B.SYNTHESIS_NOTES['C'] in str(rows.loc['C', 'note'])
    assert B.SYNTHESIS_NOTES['Mn'] in str(rows.loc['Mn', 'note'])


# ── D. the committed artifact ────────────────────────────────────────────────

def test_committed_corrected_candidate_carries_no_phantom():
    """C. The sibling sweep, as a standing check on the committed artifact."""
    df = pd.read_csv(CANDIDATE)
    claimants = {r['element']: r for _, r in df.iterrows()
                 if H.claims_zero_survivors(str(r['note']))}
    # Mg/Y/Zr/Eu are the verified-genuine zero-survivor rows (RYA-596); Ba is NOT
    # among them any more.
    assert set(claimants) == {'Mg', 'Y', 'Zr', 'Eu'}, sorted(claimants)
    for el, r in claimants.items():
        assert int(r['n_lines']) == 0 and pd.isna(r['A_X']), el


def test_committed_candidate_is_reproducible(tmp_path):
    """No hand-editing: the committed artifact is exactly what the builder emits."""
    _build(tmp_path, VERDICT)
    assert (tmp_path / 'cand.csv').read_text() == CANDIDATE.read_text()


def test_builder_runs_as_a_script(tmp_path):
    """The argv plumbing the tests use must not have broken the CLI."""
    r = subprocess.run([sys.executable, str(ROOT / 'scripts/build_solar_reference_v2_rya522.py'),
                        '--verdict', str(VERDICT),
                        '--out-csv', str(tmp_path / 'c.csv'),
                        '--out-md', str(tmp_path / 'c.md')],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / 'c.csv').exists()
