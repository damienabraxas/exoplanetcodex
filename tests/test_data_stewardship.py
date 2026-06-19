"""
tests/test_data_stewardship.py
==============================
RYA-355 — CI wrapper for the data-stewardship invariant check
(scripts/check_stewardship.py).

The contract enforced here:
  • UNTRACKED stewardship violations FAIL the build loudly (the whole point — a
    NEW duplicated-and-divergent canonical value, with no remediation ticket, goes
    red immediately).
  • The two documented backlog defects are gated as known-issues so CI stays green
    while the defect stays visible:
        - gf duplication       → RYA-353 (single-source gf migration)
        - STAR_PARAMS mirrors  → RYA-298 (legacy adapter-dict removal)
    Each is an xfail: it fails today (documenting the open defect) and flips to a
    pass when its remediation lands — at which point the xfail is removed and the
    invariant hard-passes.
  • A canonical source that cannot be parsed raises loudly; it is never skipped.
"""
import pytest

import scripts.check_stewardship as sc


# Run the full check once; share the result across tests (the gf match reads the
# 141k-line synth list, so do it a single time).
@pytest.fixture(scope='module')
def violations():
    return sc.run_all()


# ── The CI gate ───────────────────────────────────────────────────────────────
def test_no_untracked_violations(violations):
    """The hard gate: every stewardship violation must be registered against a
    remediation ticket. An untracked one means a NEW duplication slipped in."""
    untracked = [v for v in violations if not v.tracked]
    assert untracked == [], (
        "Untracked data-stewardship violation(s) — a canonical value was "
        "duplicated-and-divergent with no remediation ticket:\n"
        + "\n".join(f"  {v.invariant}: {v.quantity} @ {v.locus} = {v.value}"
                    for v in untracked))


def test_provenance_complete(violations):
    """Every canonical value carries a non-empty, non-placeholder source. This
    invariant is UNTRACKED — a missing citation fails now, not someday."""
    prov = [v for v in violations if v.invariant == 'provenance']
    assert prov == [], (
        "Canonical value(s) with empty/placeholder provenance:\n"
        + "\n".join(f"  {v.locus}: source={v.source}" for v in prov))


# ── Remediated invariants (RYA-353 + RYA-298 landed; now HARD gates) ───────────
def test_gf_tables_agree(violations):
    """RYA-353 landed: both paths resolve gf from canonical_gf via gf_resolver, so the
    matched Δgf collapses to 0. Was xfail (3,299 divergent) → now a hard pass; a new
    divergence or an orphan line (outside the single source) fails the build."""
    gf = [v for v in violations if v.invariant == 'gf']
    assert gf == [], f"{len(gf)} line(s) carry divergent/orphan gf vs the single source"


def test_star_params_single_source(violations):
    """RYA-298 landed: legacy adapter dicts no longer carry fundamental Teff/logg/feh
    literals (they live only in stars.yaml). Was xfail → now a hard pass."""
    sp = [v for v in violations if v.invariant == 'star_params']
    assert sp == [], f"{len(sp)} mirror copy/divergence(s) of STAR_PARAMS values"


# ── Behavioural guarantees of the checker itself ──────────────────────────────
def test_gf_anchor_selfcheck_runs(violations):
    """The gf match reproduced the RYA-347 Fe II anchors (a wrong match rule would
    have raised inside run_all), so the scope numbers are trustworthy."""
    summary = sc._GF_SUMMARY.get('synth-vs-solar')
    assert summary is not None
    # Large clean 1:1 cross-match (≈15k with the shared physical-line clustering);
    # the anchor self-check inside _match_gf would have raised on a wrong match rule.
    assert summary['clean'] > 13000
    # Post-RYA-353 both paths resolve to the single source → 0 material divergence,
    # and the gf-violation count tracks it.
    assert summary['div_material'] == 0
    assert summary['div_material'] == sum(
        1 for v in violations if v.invariant == 'gf')


def test_untracked_violation_fails_loudly(monkeypatch, tmp_path):
    """If any invariant returns an untracked violation, main() must exit non-zero —
    the guard cannot warn-and-continue."""
    fake = sc.Violation(
        invariant='provenance', quantity='log gf', locus='synthetic',
        value='-1.0', source="''", detail='injected for test', ticket=None)
    monkeypatch.setattr(sc, 'check_gf_pairs', lambda out_dir=None: [])
    monkeypatch.setattr(sc, 'check_star_params', lambda: [])
    monkeypatch.setattr(sc, 'check_provenance', lambda: [fake])
    assert sc.main(['--out', str(tmp_path)]) == 1


def test_tracked_violations_keep_ci_green(monkeypatch):
    """A purely-tracked violation set is a PASS (exit 0) — known backlog stays green."""
    tracked = sc.Violation(
        invariant='gf', quantity='log gf', locus='synthetic', value='x',
        source='y', detail='injected', ticket='RYA-353')
    monkeypatch.setattr(sc, 'check_gf_pairs', lambda out_dir=None: [tracked])
    monkeypatch.setattr(sc, 'check_star_params', lambda: [])
    monkeypatch.setattr(sc, 'check_provenance', lambda: [])
    assert sc.main([]) == 0


def test_parse_error_is_loud(monkeypatch, tmp_path):
    """An unparseable canonical source raises StewardshipParseError → exit 2, never
    a silent skip."""
    bad = sc.ProvenanceCheck('bad', 'csv', path=tmp_path / 'missing.csv',
                             value_col='log_gf', source_col='loggf_source')
    monkeypatch.setattr(sc, 'PROVENANCE_CHECKS', [bad])
    monkeypatch.setattr(sc, 'check_gf_pairs', lambda out_dir=None: [])
    monkeypatch.setattr(sc, 'check_star_params', lambda: [])
    with pytest.raises(sc.StewardshipParseError):
        sc.check_provenance()
    # and through main() it becomes a loud non-zero exit (2)
    assert sc.main([]) == 2
