"""
tests/test_two_engine_coverage_rya526.py
========================================
RYA-526 refresh — the two-engine coverage ledger
(`data/curation/nlte_two_engine_coverage.csv`) is the pre-declared exception
list RYA-525's loud-fail guard consumes. It must be COMPLETE (every one of the 26
TARGET_ELEMENTS, both engine columns populated, no blanks), use only the allowed
state vocabulary, be internally consistent (disposition ⇔ the two engine states),
and its `wired` claims must match disk/registry ground truth (never memory).

This is ADDITIVE: it does not touch `nlte_grid_availability.csv`, so the RYA-543
registry↔disk anti-drift test (`test_grid_coverage_rya526.py`) stays green.
"""
import csv
from pathlib import Path

import pytest

from config.constants import TARGET_ELEMENTS, NLTE_CORRECTION_ELEMENTS

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / 'data' / 'curation' / 'nlte_two_engine_coverage.csv'
GERBER = ROOT / 'data' / 'nlte_grids' / 'gerber_ts'          # Engine-B provenance stubs (Mac)
GRIDS = ROOT / 'data' / 'nlte_grids'

ENGINE_A_STATES = {'wired', 'task', 'LTE-only-by-design'}
#: RYA-695 adds `none-published`. The vocabulary previously had one word, `task`, for
#: two situations the Codex kept confusing:
#:
#:   * the Engine-B model atom EXISTS upstream and was never pulled  -> `task`
#:   * NO Engine-B model atom is published at all                    -> `none-published`
#:
#: Collapsing them cost real time. The registry said "Gerber-2023 includes Al —
#: gettable via RYA-540" and the RYA-673 audit said Al is NO_MODEL_ATOM; both were
#: `task`-shaped and nobody could tell which was true without re-deriving it. RYA-695
#: settled it against the upstream Gerber/TSFitPy catalog
#: (`utilities/nlte_grids_links.cfg`, the file the RYA-534 provenance JSONs cite as
#: `url_source`): SEVENTEEN elements are published — Al, Ba, Ca, Co, Cr, Eu, Fe, H, Mg,
#: Mn, Na, Ni, O, Si, Sr, Ti, Y — and ELEVEN are staged. C, N, K, P, S, Sc, Li, V, Zr
#: and Cu are simply not in it, so their Engine-B TS-NLTE "task" can never be completed
#: and must not be filed as though it could.
#:
#: See `data/curation/engine_b_deck_availability.csv` (generated, Sirius-only).
ENGINE_B_STATES = {'wired', 'task', 'LTE-only-by-design', 'none-published'}
DISPOSITIONS = {'wired-both', 'wired-one', 'acquire-task', 'build-task', 'LTE-only-by-design'}

# Engine-A production coverage = the NLTE correction registry + the Fe leg + C/O (nlte_cno).
ENGINE_A_WIRED_TRUTH = set(NLTE_CORRECTION_ELEMENTS) | {'Fe', 'C', 'O'}


def _rows():
    with open(CSV, newline='') as fh:
        return list(csv.DictReader(fh))


def test_every_target_element_present_exactly_once():
    rows = _rows()
    els = [r['element'] for r in rows]
    assert set(els) == set(TARGET_ELEMENTS), (
        f"ledger vs TARGET_ELEMENTS mismatch: "
        f"missing {set(TARGET_ELEMENTS) - set(els)}, extra {set(els) - set(TARGET_ELEMENTS)}")
    assert len(els) == len(set(els)) == len(TARGET_ELEMENTS)   # no dupes, all 26


def test_no_blanks_and_valid_vocab():
    for r in _rows():
        el = r['element']
        assert r['engineA_1dnlte'] in ENGINE_A_STATES, f"{el}: bad engineA {r['engineA_1dnlte']!r}"
        assert r['engineB_synth_grid'] in ENGINE_B_STATES, f"{el}: bad engineB {r['engineB_synth_grid']!r}"
        assert r['disposition'] in DISPOSITIONS, f"{el}: bad disposition {r['disposition']!r}"
        assert r['reason'].strip(), f"{el}: empty reason"
        assert r['governing_ticket'].strip(), f"{el}: empty governing_ticket"


def test_disposition_is_consistent_with_engine_states():
    for r in _rows():
        a, b, disp, el = r['engineA_1dnlte'], r['engineB_synth_grid'], r['disposition'], r['element']
        wired = {'wired'}
        if a in wired and b in wired:
            assert disp == 'wired-both', f"{el}: both wired but disposition={disp}"
        elif a == 'LTE-only-by-design' and b == 'LTE-only-by-design':
            assert disp == 'LTE-only-by-design', f"{el}: both LTE-by-design but disposition={disp}"
        elif a in wired or b in wired:
            assert disp == 'wired-one', f"{el}: one engine wired but disposition={disp}"
        else:
            # neither engine wired: grid gettable -> acquire-task; genuine void -> build-task
            assert disp in {'acquire-task', 'build-task'}, (
                f"{el}: neither wired/by-design but disposition={disp}")


def test_none_published_is_backed_by_the_upstream_catalog_ledger():
    """RYA-695 — `none-published` is a factual claim about the upstream catalog, so it
    is checked against the generated deck ledger rather than trusted as prose.

    The ledger enumerates every element the Gerber/TSFitPy catalog publishes. An
    element claiming `none-published` must be ABSENT from it; an element claiming
    `task` must be PRESENT (otherwise the task can never be completed, which is the
    conflation this state was added to end).
    """
    ledger = ROOT / 'data' / 'curation' / 'engine_b_deck_availability.csv'
    if not ledger.exists():
        pytest.skip('engine_b_deck_availability.csv not generated (Sirius-only)')
    with open(ledger, newline='') as fh:
        upstream = {r['element'] for r in csv.DictReader(fh)}
    for r in _rows():
        el, b = r['element'], r['engineB_synth_grid']
        if b == 'none-published':
            assert el not in upstream, (
                f"{el}: registry says engineB=none-published but the upstream Gerber "
                f"catalog DOES publish it — that makes it an acquisition `task`")
        if b == 'task':
            assert el in upstream, (
                f"{el}: registry says engineB=task (an acquisition owed) but the "
                f"upstream Gerber catalog does not publish it — the task can never be "
                f"completed; it is `none-published`")


def test_engineA_wired_matches_registry_truth():
    # a 'wired' Engine-A claim must be backed by the NLTE registry (+ Fe/C/O); never memory
    for r in _rows():
        if r['engineA_1dnlte'] == 'wired':
            assert r['element'] in ENGINE_A_WIRED_TRUTH, (
                f"{r['element']}: engineA=wired but absent from NLTE_CORRECTION_ELEMENTS + Fe/C/O")


def test_engineB_wired_matches_gerber_provenance_on_disk():
    # a 'wired' Engine-B claim must have its RYA-534 TS-Gerber provenance stub on disk
    for r in _rows():
        if r['engineB_synth_grid'] == 'wired':
            prov = GERBER / f"{r['element']}_gerber2023.prov.json"
            assert prov.exists(), f"{r['element']}: engineB=wired but missing {prov.name}"


def test_v_is_the_only_build_task():
    builds = [r['element'] for r in _rows() if r['disposition'] == 'build-task']
    assert builds == ['V'], f"expected V the only NLTE_VOID build-task, got {builds}"


def test_state_counts_match_the_declared_ledger():
    from collections import Counter
    c = Counter(r['disposition'] for r in _rows())
    assert c['wired-both'] == 9
    assert c['wired-one'] == 9
    assert c['acquire-task'] == 3
    assert c['build-task'] == 1
    assert c['LTE-only-by-design'] == 4
