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
ENGINE_B_STATES = {'wired', 'task', 'LTE-only-by-design'}
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
