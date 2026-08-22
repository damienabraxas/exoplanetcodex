"""RYA-935 — the tracker's telluric section must render each metric under its own head.

The page is a VIEW over live_status.json and holds no data of its own, so these tests
assert the contract between the two: every metric the collector can emit has somewhere
to land, and nothing renders under a heading that would misdescribe it.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT / 'data' / 'results' / 'rya935' / 'live_tracker.html'
JSON = ROOT / 'data' / 'results' / 'rya935' / 'live_status.json'


def _rows():
    return json.loads(JSON.read_text()).get('telluric', [])


def _metric(t):
    return t.get('metric', 'pct_below_0.5')


def test_every_metric_the_collector_emits_has_a_group_to_render_it():
    """🔴 The regression this prevents: the section had ONE renderer, keyed on
    before_pct_below_0.5. A D1 row does not carry that field, so six alpha Cen frames
    whose gates all PASS at 80-91% improvement rendered as 'not corrected' — the
    dashboard stating the exact opposite of the truth, in the one place a reader goes to
    find out what is corrected."""
    html = HTML.read_text()
    seen = {_metric(t) for t in _rows()}
    targets = {'pct_below_0.5': 'tell', 'd1_residual': 'tellstar',
               'state_only__NOT_corrected': 'tellstate'}
    for m in seen:
        assert m in targets, f"metric {m!r} is emitted but the page has no group for it"
        assert f"#{targets[m]}" in html or f'id="{targets[m]}"' in html


def test_the_three_metrics_are_not_stacked_in_one_column():
    """They answer different questions and are not comparable as numbers (RYA-873:
    report a value under its DERIVED name). Separate headings encode that."""
    html = HTML.read_text()
    for anchor in ('id="tell"', 'id="tellstar"', 'id="tellstate"'):
        assert anchor in html
    assert 'not comparable as numbers' in html


def test_a_state_determination_never_renders_as_a_correction():
    """A STATE row measures what the product IS; it carries no before/after because
    nothing was corrected. Showing it in a before->after column would read as a
    correction that never happened."""
    html = HTML.read_text()
    assert 'NOT corrected' in html
    for t in _rows():
        if _metric(t).startswith('state_only'):
            assert 'before_d1_residual' not in t
            assert t.get('before_pct_below_0.5') is None


def test_a_contested_star_id_cannot_render_as_a_settled_one():
    """alpha Cen's A/B verdict rests on the acen_orbit branch assignment RYA-963 left
    CONTESTED. not-contested and nobody-recorded-it are different states and the page
    must not collapse them."""
    html = HTML.read_text()
    assert 'CONTESTED' in html
    assert 'unrecorded' in html


def test_the_stellar_rows_carry_what_the_group_renders():
    for t in _rows():
        if _metric(t) == 'd1_residual':
            for k in ('before_d1_residual', 'after_d1_residual', 'gate_passed', 'window'):
                assert k in t, f"stellar row missing {k}"


def test_the_page_holds_no_data_of_its_own():
    """GENERATORS.yaml registers the page against the JSON's generator precisely because
    the page must stay a view. A number hard-coded here would drift silently."""
    html = HTML.read_text()
    # no long decimal literals outside the CSS/JS scaffolding
    body = html.split('<script>')[-1]
    numbers = re.findall(r'\b\d+\.\d{3,}\b', body)
    assert not numbers, f"the page appears to hard-code data: {numbers[:5]}"
