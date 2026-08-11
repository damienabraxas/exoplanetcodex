"""
tests/test_species_row_order_rya768.py
======================================
RYA-768 — the two-engine species artifact must emit rows in a DETERMINISTIC order.
Part of the RYA-166 solar regression-anchor family.

WHY THIS TEST EXISTS
--------------------
RYA-765's proof (a) ran the two-engine floor twice with the tracer disabled in both
runs and got two artifacts that differed. The record SETS were equal, so no number had
moved — only the row ORDER. Root cause: the emit sorted on ``(-A_asplund, element)``,
which is not a total order, because ``A_asplund`` is per-ELEMENT. Fe I and Fe II tie on
the whole key; so do Cr I/Cr II and Ti I/Ti II. Python's sort is stable, so a tie is
resolved by INPUT order, and the input was a ``set`` of ``(element, ion)`` tuples whose
iteration order varies per process under hash randomisation.

An artifact that never diffs clean is how a real change hides. Every diff-based check
downstream — the results-ledger consistency guard, any solar regression anchor that
compares artifacts, the Al finalisation proofs — is only as sharp as this ordering.

WHAT IS AND IS NOT ASSERTED
---------------------------
This guards the ORDERING CONTRACT, not the science: no assertion here reads or compares
an abundance. It also does NOT assert that the COMMITTED artifacts are in canonical
order — they are not, and deliberately so: they were emitted by the old key (the phase-3
artifact carries 'Si II' before 'Si I', a frozen instance of the very tie this fixes)
and rewriting historical evidence is not this ticket's job. They become canonical the
next time they are emitted.
"""
import itertools
import json
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import SOLAR_ASPLUND2021          # noqa: E402
from pipeline.engine_selection import species_row_sort_key  # noqa: E402

#: A committed emit, used ONLY as a realistic supply of (element, ion) pairs — never as
#: an expected order. Reading the species set from a real artifact keeps the fixture
#: honest (no hand-written species list to drift from reality).
_COMMITTED = ROOT / 'data' / 'audit' / 'rya527_phase3' / 'solar_two_engine_records.json'

#: The defective key, reproduced verbatim so the bite test below exercises the real
#: thing rather than a paraphrase of it. Do not "fix" this — it is the negative control.
def _legacy_key(element, ion):
    return (-SOLAR_ASPLUND2021.get(element, -9), element)


def _species_pairs():
    if not _COMMITTED.exists():
        pytest.skip(f"no committed two-engine artifact at {_COMMITTED}")
    records = json.loads(_COMMITTED.read_text())['records']
    return [(str(r['element']), str(r['ion'])) for r in records]


def test_key_is_a_total_order_over_the_real_species_set():
    """No two rows may share a key — that is what 'total order' buys."""
    pairs = _species_pairs()
    keys = [species_row_sort_key(e, i) for e, i in pairs]
    dupes = [k for k in keys if keys.count(k) > 1]
    assert not dupes, (
        f"{len(dupes)} row(s) share a sort key, so their relative order is decided by "
        f"input order: {sorted(set(map(str, dupes)))}")
    assert len(set(keys)) == len(pairs)


def test_every_key_pair_is_strictly_comparable():
    """Totality the strict way: for any two distinct rows, exactly one precedes."""
    keys = [species_row_sort_key(e, i) for e, i in _species_pairs()]
    for a, b in itertools.combinations(keys, 2):
        assert (a < b) != (b < a), f"{a} and {b} are not strictly ordered"


def test_sort_is_invariant_to_input_order():
    """THE regression. Shuffle the input; the emitted order must not move."""
    pairs = _species_pairs()
    rng = random.Random(20260810)
    canonical = sorted(pairs, key=lambda p: species_row_sort_key(*p))
    for _ in range(200):
        shuffled = pairs[:]
        rng.shuffle(shuffled)
        assert sorted(shuffled, key=lambda p: species_row_sort_key(*p)) == canonical


def test_guard_bites_the_legacy_key():
    """Negative control: the OLD key fails the test above.

    Without this, a future refactor could weaken the key back toward per-element and
    the suite would stay green because nothing proved the check discriminates.
    """
    pairs = _species_pairs()
    rng = random.Random(20260810)
    orders = {tuple(sorted(rng.sample(pairs, len(pairs)),
                           key=lambda p: _legacy_key(*p)))
              for _ in range(200)}
    assert len(orders) > 1, (
        "the legacy key produced ONE stable order over 200 shuffles — the fixture no "
        "longer contains a two-ion element, so this guard has stopped discriminating")

    legacy_keys = [_legacy_key(e, i) for e, i in pairs]
    assert len(set(legacy_keys)) < len(pairs), "legacy key unexpectedly had no ties"


def test_mixed_ion_representations_rank_identically():
    """RYA-345: 'I' / '1' / 1 are ONE stage. Sorting them as raw text would put '10'
    before '2' and 'II' before 'I', reintroducing the nondeterminism being removed."""
    assert (species_row_sort_key('Fe', 'I')
            == species_row_sort_key('Fe', '1')
            == species_row_sort_key('Fe', 1))
    assert species_row_sort_key('Fe', 'I') < species_row_sort_key('Fe', 'II')
    assert species_row_sort_key('Fe', 'II') < species_row_sort_key('Fe', 'X')  # 2 < 10


def test_primary_order_is_still_abundance_descending():
    """ORDER-ONLY: the fix may resolve ties, never re-rank distinct elements."""
    pairs = _species_pairs()
    ordered = sorted(pairs, key=lambda p: species_row_sort_key(*p))
    abunds = [SOLAR_ASPLUND2021.get(e, -9) for e, _ in ordered]
    assert abunds == sorted(abunds, reverse=True)


def test_unparseable_ion_sorts_last_without_raising():
    """A sort on an emit path must not turn a data defect into a crash after the
    science is done; it must still be deterministic."""
    k = species_row_sort_key('Fe', 'not-an-ion')
    assert k > species_row_sort_key('Fe', 'II')
    assert k == species_row_sort_key('Fe', 'not-an-ion')


def test_emitter_uses_the_shared_key_not_a_local_copy():
    """Single source of truth: the emit site must import the helper, so the ordering
    cannot drift from this test's subject."""
    try:
        # The driver pulls in abundances_derive -> iSpec, which is not vendored; skip
        # rather than fail where the engines are absent (CI runs on Sirius, which has them).
        sys.path.insert(0, str(ROOT / 'scripts'))
        import rya527_two_engine_run as driver
    except Exception as exc:                        # pragma: no cover - env-dependent
        pytest.skip(f"two-engine driver not importable in this environment: {exc}")
    assert driver.species_row_sort_key is species_row_sort_key
