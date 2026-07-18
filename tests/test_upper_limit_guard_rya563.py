"""
tests/test_upper_limit_guard_rya563.py
======================================
RYA-563 — Lithium (Li I) UPPER_LIMIT disposition guard.

Li I 6707.84 carries the registry disposition `required_treatment=upper_limit`
(data/registry/problem_children.csv; RYA-103/458: "CN-blended; a clean low value
is a RED FLAG. Carried as UPPER_LIMIT, never a point value."). The reference-blind
two-engine floor previously fell through to the synthesis branch and emitted the
Engine-B point value 1.409 (+0.359), overriding the disposition.

These pin the LAW: the registry-sourced `is_upper_limit_disposition` helper
(single source of truth = problem_children.csv, no hardcoded element list), and
that the re-emitted Li verdict reports the phase_c UPPER-LIMIT value (0.727) with
1.409 demoted to a DIAGNOSTIC-ONLY record — never the reported value.
"""
import json
import subprocess
import sys
from pathlib import Path

from pipeline import engine_selection as es

ROOT = Path(__file__).resolve().parent.parent


def test_li_is_upper_limit_disposition():
    # Li I carries required_treatment=upper_limit in the registry.
    assert es.is_upper_limit_disposition('Li') is True


def test_non_upper_limit_element_is_false():
    # Fe has no upper_limit disposition — the helper must not over-match.
    assert es.is_upper_limit_disposition('Fe') is False


def test_reemit_li_reports_upper_limit_not_synth():
    """The re-emitted Li verdict must report the phase_c UPPER-LIMIT value (0.727),
    NEVER the two-engine synthesis point value 1.409 (the RYA-103 red flag)."""
    subprocess.run([sys.executable, 'scripts/rya527_reemit_verdict.py'],
                   cwd=ROOT, check=True)
    payload = json.loads((ROOT / 'data' / 'audit' / 'rya527_reemit'
                          / 'proposed_gold_v3_diff.json').read_text())
    li = next(r for r in payload['diff_table'] if r['element'] == 'Li')

    # CRITICAL FAILURE: the synth point value must not survive as the reported value.
    assert li['v3_proposed'] != 1.409
    assert li['v3_proposed'] == 0.727

    # 1.409 is recorded, but only as a DIAGNOSTIC-ONLY species, never the value field.
    ter = li['two_engine_record'] or {}
    assert ter.get('reported') is None
    diag = ter.get('diagnostic_only_species') or []
    assert any(d.get('value') == 1.409 and d.get('DIAGNOSTIC_ONLY') for d in diag)
