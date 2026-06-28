"""
tests/test_problem_children_rya463.py
=====================================
RYA-463 — guard the master problem-children registry: the curated seed, the
RYA-458 auto-aggregate ingest + idempotent upsert, and the prediction layer
reproducing the hand-assembled Sun / Procyon / 55 Cnc watch-lists. Catalog only —
a test asserts it never touches a measured abundance.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import problem_children as P                      # noqa: E402


def _curated_only():
    return P.build_registry(star_files={})


# ── curated layer ─────────────────────────────────────────────────────────────
def test_curated_seed_has_the_charter_and_architectural_rows():
    df = _curated_only()
    sci = set(zip(df['species'], df['lambda_or_scope']))
    assert ('O I', '[O I] 6300') in sci          # the O-indicator risk (RYA-455)
    assert ('C I', '5380.34') in sci             # RYA-458 charter
    assert ('Li I', '6707.84') in sci            # upper-limit charter
    assert ('Eu II', '6645.13') in sci           # HFS-summing charter
    assert ('N I', '7442-8718 (red multiplets)') in sci
    assert any(s.startswith('CH/CN/C2') for s in df['species'])   # molecular carve-out


def test_curated_classes_and_treatments_are_in_vocab():
    df = _curated_only()
    for _, r in df.iterrows():
        # curated problem_class may be a '+'-combo; each part must be a known class
        for part in str(r['problem_class']).split('+'):
            assert part in P.CURATED_CLASSES, part
        assert r['required_treatment'] in P.TREATMENTS
        assert r['status'] in P.STATUSES
        assert r['severity'] in P.SEVERITIES


# ── auto-aggregate + idempotency ──────────────────────────────────────────────
def _fake_ew_integrity(tmp_path):
    df = pd.DataFrame([
        # element, ion, wave, ew, ew_integrity, ew_disposition
        ('C', 'I', 5380.337, 149.5, 'COG_FLAG,BAD_FIT', 'BAD_FIT'),
        ('Fe', 'I', 4924.301, 120.0, 'COG_FLAG,ABUND_OUTLIER', None),
        ('Fe', 'I', 5000.000, 130.0, 'COG_FLAG', None),       # routine saturation
        ('Ti', 'I', 4533.000, 110.0, 'COG_FLAG', None),       # routine saturation
        ('Na', 'I', 5889.950, 40.0, '', None),                # clean -> not flagged
    ], columns=['element', 'ion', 'wavelength_air_A', 'ew_mA', 'ew_integrity', 'ew_disposition'])
    p = tmp_path / 'solar_ew_integrity.csv'
    df.to_csv(p, index=False)
    return p


def test_aggregate_keeps_genuine_perline_collapses_cog(tmp_path):
    rows = P.aggregate_ew_integrity('solar', _fake_ew_integrity(tmp_path))
    by = {(r['species'], r['lambda_or_scope']): r for r in rows}
    # genuine per-line outliers kept
    assert ('C I', '5380.337') in by and 'BAD_FIT' in by[('C I', '5380.337')]['problem_class']
    assert ('Fe I', '4924.301') in by                      # ABUND_OUTLIER kept per-line
    # routine COG saturation collapsed per (species, ion) — Fe gets ONE summary row
    fe_cog = [r for r in rows if r['species'] == 'Fe I' and 'strong-line pool' in r['lambda_or_scope']]
    assert len(fe_cog) == 1
    assert all(r['population_source'] == 'auto_ew_integrity' for r in rows)


def test_upsert_is_idempotent(tmp_path):
    path = _fake_ew_integrity(tmp_path)
    df1 = P.upsert([], P.aggregate_ew_integrity('solar', path))
    df2 = P.upsert(df1, P.aggregate_ew_integrity('solar', path))   # re-ingest same star
    assert len(df1) == len(df2)                            # no duplicate rows
    # observed_in not duplicated
    assert all(s.count('solar') <= 1 for s in df2['observed_in'])


def test_upsert_appends_new_star_to_observed_in(tmp_path):
    path = _fake_ew_integrity(tmp_path)
    df = P.upsert([], P.aggregate_ew_integrity('solar', path))
    df = P.upsert(df, P.aggregate_ew_integrity('procyon', path))   # same lines, new star
    ci = df[(df['species'] == 'C I') & (df['lambda_or_scope'] == '5380.337')].iloc[0]
    assert 'solar' in ci['observed_in'] and 'procyon' in ci['observed_in']


# ── prediction layer (the payoff) ─────────────────────────────────────────────
def test_predict_sun_returns_the_solar_set():
    df = _curated_only()
    preds = P.predict(5777, 0.0, df, star_name='Sun')
    cls = {(p['species'], p['problem_class']) for p in preds}
    assert ('O I', 'CONTINUUM_LIMITED') in cls
    assert ('N I', 'NLTE_OWED') in cls
    assert any(s == 'Cr I' and c == 'BAD_GF' for s, c in cls)


def test_predict_procyon_is_the_fstar_watch_list():
    df = _curated_only()
    preds = P.predict(6554, 0.01, df, star_name='Procyon')
    amp = {p['species'] for p in preds if p['status'] == 'amplified'}
    # F-star: N I NLTE-owed + bad-gf metals amplify with Teff↑
    assert 'N I' in amp
    assert {'Cr I', 'Ti I', 'Si I'} <= amp
    # [O I] 6300 is present but only a WATCH at this [Fe/H] (not amplified)
    o = next(p for p in preds if p['species'] == 'O I')
    assert o['status'] == 'watch'
    # the COG knee heads-up shifts up for the F-star
    assert P.estimated_cog_ceiling_mA(6554) > P.SOLAR_SAT_CEILING_MA


def test_predict_55cnc_amplifies_metal_rich_set():
    df = _curated_only()
    preds = P.predict(5196, 0.32, df, star_name='55 Cnc')
    amp = {p['species'] for p in preds if p['status'] == 'amplified'}
    # very metal-rich: [O I] 6300 amplifies; cool: molecular bands amplify
    o = next(p for p in preds if p['species'] == 'O I')
    assert o['status'] == 'amplified'
    assert any(s.startswith('CH/CN/C2') for s in amp)


def test_headsup_calls_out_oxygen_indicator():
    df = _curated_only()
    hp = P.predict_headsup(6554, 0.01, df, star_name='Procyon')
    assert '[O I] 6300' in hp and '777' in hp


# ── cardinal rule: catalog only, no abundance mutation ────────────────────────
def test_registry_is_catalog_only_no_abundance_columns():
    df = _curated_only()
    forbidden = {'A_X', 'a_lte', 'A_X_nlte', 'abundance', 'ew_mA'}
    assert not (set(df.columns) & forbidden)        # never carries a measured value
    assert list(df.columns) == P.SCHEMA_COLUMNS
