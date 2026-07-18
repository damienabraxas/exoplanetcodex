"""
tests/test_upper_limit_guard_rya563.py
======================================
RYA-563 — Lithium (Li I) UPPER_LIMIT disposition guard (durable helper).

Li I 6707.84 carries the registry disposition `required_treatment=upper_limit`
(data/registry/problem_children.csv; RYA-103/458: "CN-blended; a clean low value
is a RED FLAG. Carried as UPPER_LIMIT, never a point value."). The reference-blind
two-engine floor may NOT emit a synthesis point value for such an element.

These pin the LAW of the registry-sourced `is_upper_limit_disposition` helper —
single source of truth = problem_children.csv, no hardcoded element list. The
consumer-side enforcement in the RYA-527 re-emit ladder is exercised separately on
the RYA-527 stack (where that script lives).
"""
from pipeline import engine_selection as es


def test_li_is_upper_limit_disposition():
    # Li I carries required_treatment=upper_limit in the registry.
    assert es.is_upper_limit_disposition('Li') is True


def test_non_upper_limit_element_is_false():
    # Fe has no upper_limit disposition — the helper must not over-match.
    assert es.is_upper_limit_disposition('Fe') is False
