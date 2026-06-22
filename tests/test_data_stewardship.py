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


# ── RYA-358: blend_flag definition pin ────────────────────────────────────────
def test_blend_flag_definition_holds(violations):
    """blend_flag in linelist_solar matches the re-run vetted builder (RYA-209), the
    vetted blends carry provenance, and nothing flips silently → 0 violations."""
    bf = [v for v in violations if v.invariant == 'blend_flag']
    assert bf == [], f"{len(bf)} blend_flag violation(s): {[v.locus for v in bf]}"
    assert sc._BLEND_SUMMARY.get('mismatch') == 0


def test_vetted_blends_carry_provenance():
    """Each VETTED_BLENDS entry has a (element, ion, wl, source) shape with a real
    citation — the provenance the guard enforces."""
    from pipeline.build_linelist import VETTED_BLENDS
    for entry in VETTED_BLENDS:
        assert len(entry) >= 4, f"VETTED_BLENDS entry missing source: {entry}"
        assert not sc._is_placeholder(entry[3]), f"placeholder source: {entry}"


def test_blend_flag_pin_fails_on_silent_redefinition(monkeypatch, tmp_path):
    """The guard's whole point: if the file's blend_flag diverges from the vetted
    builder (a silent redefinition or a builder swap back to proximity), it must
    produce an UNTRACKED violation that fails the build."""
    import pandas as pd
    bad = tmp_path / 'linelist_mutated.csv'
    pd.DataFrame({
        'element': ['Fe', 'Ca'],
        'ion': ['I', 'I'],
        'wavelength_air_A': [4970.496, 5000.000],
        'blend_flag': [False, True],   # vetted builder says True/False → both mismatch
    }).to_csv(bad, index=False)
    monkeypatch.setattr(sc, '_LL_SOLAR', bad)
    monkeypatch.setattr(sc, '_SOLAR_EW', tmp_path / 'no_ew.csv')  # propagation skipped
    v = sc.check_blend_flag()
    defn = [x for x in v if x.quantity == 'blend_flag definition']
    assert len(defn) == 2, f"expected 2 definition mismatches, got {len(defn)}"
    assert all(not x.tracked for x in defn), "a silent redefinition must be UNTRACKED"
    # and it fails the build loudly through main()
    monkeypatch.setattr(sc, 'check_gf_pairs', lambda out_dir=None: [])
    monkeypatch.setattr(sc, 'check_star_params', lambda: [])
    monkeypatch.setattr(sc, 'check_provenance', lambda: [])
    assert sc.main(['--out', str(tmp_path)]) == 1


# ── RYA-408: solar-EW canonical input source ──────────────────────────────────
def test_solar_ew_canonical_invariant_clean(violations):
    """On a clean checkout (no staging file) the canonical is present + well-formed
    and the gate loader reads it → 0 violations from this invariant."""
    sec = [v for v in violations if v.invariant == 'solar_ew_canonical']
    assert sec == [], f"{len(sec)} solar_ew_canonical violation(s): " \
                      f"{[(v.quantity, v.locus) for v in sec]}"


def test_gate_loader_reads_canonical_not_staging():
    """IDENTITY: _load_solar_ews must read the committed canonical and must NOT read
    the gitignored staging file as its EW pool (the RYA-406 incident)."""
    import inspect
    src = inspect.getsource(sc._ad._load_solar_ews)
    assert 'solar_ew_canonical' in src
    assert "read_csv(str(PATHS['solar_ew']))" not in src


def test_repoint_to_runtime_is_untracked_break(monkeypatch):
    """If a regression re-points the loader at the gitignored staging file, the
    identity guard fires an UNTRACKED violation."""
    def _fake_loader(ew_override=None):
        import pandas as pd
        return pd.read_csv(str(PATHS['solar_ew']))   # the forbidden read pattern
    monkeypatch.setattr(sc._ad, '_load_solar_ews', _fake_loader)
    v = sc.check_solar_ew_canonical()
    ident = [x for x in v if x.quantity == 'gate EW input source']
    assert len(ident) == 1 and not ident[0].tracked


def test_staging_drift_fires_untracked(monkeypatch, tmp_path):
    """DRIFT: a present staging file diverging from the canonical on a measured EW is
    an UNTRACKED loud failure (stale / different-run staging masquerading as source)."""
    import pandas as pd
    canon = pd.read_csv(sc._SOLAR_EW_CANONICAL, low_memory=False)
    stage = canon.copy()
    stage.loc[stage.index[0], 'ew_mA'] = float(stage.iloc[0]['ew_mA']) + 9.0  # >0.5 mÅ
    staging = tmp_path / 'solar_ew.csv'
    stage.to_csv(staging, index=False)
    monkeypatch.setitem(sc._const.PATHS, 'solar_ew', staging)
    v = sc.check_solar_ew_canonical()
    drift = [x for x in v if x.quantity == 'staging↔canonical EW drift']
    assert len(drift) >= 1 and all(not x.tracked for x in drift)


def test_matching_staging_no_drift(monkeypatch, tmp_path):
    """A staging file that agrees with the canonical (a resolved subset of the same
    run) produces NO drift violation — curated blend_flags are not compared here."""
    import pandas as pd
    canon = pd.read_csv(sc._SOLAR_EW_CANONICAL, low_memory=False)
    stage = canon.copy()
    stage['blend_flag'] = False   # raw staging legitimately lacks the curated flags
    staging = tmp_path / 'solar_ew.csv'
    stage.to_csv(staging, index=False)
    monkeypatch.setitem(sc._const.PATHS, 'solar_ew', staging)
    drift = [x for x in sc.check_solar_ew_canonical()
             if x.quantity == 'staging↔canonical EW drift']
    assert drift == []
