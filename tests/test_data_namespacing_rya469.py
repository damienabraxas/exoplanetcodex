"""
tests/test_data_namespacing_rya469.py
=====================================
RYA-469 — per-star output namespacing + the frozen/versioned gold-standard solar
reference + immutability guard. Covers the ticket's Section-3 smoke test:

  * two stars -> two distinct namespaced files, no collision (Deliverable B);
  * solar_abundances_v1 exists, provenance-stamped, frozen; editing it trips the
    immutability guard (Deliverables C + E);
  * the differential denominator reads a pinned version + a target stamp records it
    (Deliverable D);
  * promote_solar_reference creates v2 without touching v1, Fe anchor 7.516 (C).

Promotion/immutability tests run against a TEMP reference dir (monkeypatched) so they
never mutate the committed gold.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import data_namespace as ns  # noqa: E402


# ── Deliverable B: per-star namespacing, no collision ────────────────────────
def test_two_stars_get_distinct_namespaced_paths():
    sol = ns.output_path('solar', 'abundances.csv', create=False)
    pro = ns.output_path('procyon', 'abundances.csv', create=False)
    assert sol != pro
    assert sol.parent != pro.parent                 # star is in the PATH
    assert sol.name == 'solar_abundances.csv'
    assert pro.name == 'procyon_abundances.csv'
    # every product type stays under the star's own dir
    for name in ('abundances.csv', 'per_line.csv', 'ew_integrity.csv', 'verdict.json'):
        assert ns.output_path('procyon', name, create=False).parent.name == 'procyon'
    assert ns.diagnostics_dir('procyon', create=False).parent.name == 'procyon'


def test_star_slug_normalises_sun_aliases():
    assert ns.star_slug('Sun') == 'solar'
    assert ns.star_slug('sol') == 'solar'
    assert ns.star_slug('Tau Boo') == 'tau_boo'
    with pytest.raises(ValueError):
        ns.star_slug('   ')


def test_output_path_does_not_double_prefix():
    assert ns.output_path('solar', 'solar_abundances.csv', create=False).name == 'solar_abundances.csv'


# ── Deliverable C: v1 exists, frozen, provenance-stamped ─────────────────────
def test_v1_exists_with_provenance_and_final_verdict():
    # RYA-522 froze gold v2 from the verdict channel and moved the CURRENT pointer to it
    # (v1 retained immutable + SUPERSEDED — the C=10.26 RYA-520 saturated-C-I-5380
    # artifact); RYA-665 then froze v3 from the RYA-653 corrected candidate and moved the
    # pointer again. So CURRENT is no longer v1. What this test owns is unchanged: v1
    # still EXISTS, is still frozen, and still carries its own provenance and verdict.
    # Read it by explicit version rather than through the moving pointer — this assertion
    # only pins that the pointer has MOVED OFF v1, not where it currently sits.
    assert ns.current_version() != 'v1'
    df, v = ns.read_solar_reference('v1')
    assert v == 'v1'
    prov = ns.read_provenance('v1')
    assert prov['version'] == 'v1'
    assert 'PASS=4' in prov['verdict'] and 'NLTE-OWED=1' in prov['verdict']
    assert 'CURATION-OWED=21' in prov['verdict'] and 'DATA-GAP=0' in prov['verdict']
    fe1 = df[(df['element'] == 'Fe') & (df['ion'] == 'I')].iloc[0]
    assert abs(float(fe1['A_X_nlte_absolute']) - 7.516) < 1e-6


def test_v1_passes_the_immutability_guard_as_committed():
    # the committed manifest must match the committed v1 file
    results = dict((fn, ok) for fn, ok, _ in ns.verify_frozen_references())
    assert results.get('solar_abundances_v1.csv') is True
    ns.assert_frozen_references()                     # does not raise


# ── Deliverable E: editing a frozen version trips the guard ──────────────────
def _seed_v1(tmp_dir, monkeypatch):
    """Point the namespace at a temp reference dir holding a copy of the real v1."""
    df, _ = ns.read_solar_reference('v1')             # real committed data FIRST
    monkeypatch.setattr(ns, 'SOLAR_REFERENCE_DIR', tmp_dir)
    monkeypatch.setattr(ns, 'CURRENT_POINTER', tmp_dir / 'CURRENT')
    monkeypatch.setattr(ns, 'HASH_MANIFEST', tmp_dir / 'hash_manifest.json')
    tmp_dir.mkdir(parents=True, exist_ok=True)
    # write into the temp dir via the production writer (header + hash recorded)
    ns.write_reference_version('v1', df, {'changelog': 'temp copy for test'})
    ns.set_current('v1')
    return tmp_dir


def test_editing_a_frozen_version_fails_the_guard(tmp_path, monkeypatch):
    _seed_v1(tmp_path, monkeypatch)
    ns.assert_frozen_references()                      # clean to start
    # tamper: append a byte to the frozen file
    target = ns.reference_path('v1')
    target.write_text(target.read_text() + '\n# tampered\n')
    with pytest.raises(ns.ImmutableReferenceError):
        ns.assert_frozen_references()


def test_unrecorded_version_on_disk_fails_the_guard(tmp_path, monkeypatch):
    _seed_v1(tmp_path, monkeypatch)
    # drop a v2 file on disk that was never recorded in the manifest
    ns.reference_path('v2').write_text('# rogue\nelement,ion,A_X\nFe,I,7.5\n')
    bad = [fn for fn, ok, _ in ns.verify_frozen_references() if not ok]
    assert 'solar_abundances_v2.csv' in bad


# ── Deliverable C: promotion bumps the version, never overwrites ──────────────
def test_promotion_creates_v2_without_touching_v1(tmp_path, monkeypatch):
    _seed_v1(tmp_path, monkeypatch)
    v1_hash_before = ns.load_manifest()['solar_abundances_v1.csv']
    df, _ = ns.read_solar_reference('v1')
    # promote a v2 (same data is fine for the structural test)
    ns.write_reference_version('v2', df, {'changelog': 'second baseline', 'supersedes': 'v1'})
    ns.set_current('v2')
    assert ns.current_version() == 'v2'
    assert ns.list_versions() == ['v1', 'v2']
    # v1 untouched: hash stable, still passes the guard
    assert ns.load_manifest()['solar_abundances_v1.csv'] == v1_hash_before
    ns.assert_frozen_references()


def test_write_refuses_to_overwrite_existing_version(tmp_path, monkeypatch):
    _seed_v1(tmp_path, monkeypatch)
    df, _ = ns.read_solar_reference('v1')
    with pytest.raises(ns.ReferenceVersionExists):
        ns.write_reference_version('v1', df, {'changelog': 'illegal overwrite'})


# ── Deliverable D: differential denominator pins + stamps a version ──────────
def test_differential_denominator_returns_pinned_version():
    df, v = ns.differential_denominator('CURRENT')
    assert v == ns.current_version()
    assert not df.empty


def test_stamp_records_solar_ref_version_on_a_target():
    target = pd.DataFrame({'element': ['C', 'O'], 'ion': ['I', 'I'], 'A_X': [8.4, 8.7]})
    stamped = ns.stamp_solar_ref_version(target, 'CURRENT')
    assert (stamped['solar_ref_version'] == ns.current_version()).all()
    assert 'solar_ref_version' not in target.columns   # original untouched


def test_phase_c_records_solar_ref_version_in_the_verdict():
    # the committed namespaced solar verdict carries the pinned denominator version
    import json
    p = ns.output_path('solar', 'verdict.json', create=False)
    if not p.exists():
        pytest.skip("solar verdict not generated in this checkout")
    summ = json.loads(p.read_text())['summary']
    assert summ['solar_ref_version'] == ns.current_version()
