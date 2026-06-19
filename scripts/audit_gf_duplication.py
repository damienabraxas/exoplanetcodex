"""
scripts/audit_gf_duplication.py
===============================
RYA-350 Step 0 — SYSTEMIC SCOPE AUDIT (read-only).

The synth path (iSpec GES v6 `atomic_lines.tsv`, consumed by the flux fit) and the
EW path (`linelist_solar.csv`, VALD3) carry an independent `loggf` for the same
physical line. Because gf is degenerate with the floated abundance, a cross-file
Δgf is a PURE differential abundance bias between the two engines (RYA-347 proved
χ²ᵣ is invariant under a gf change). This script quantifies the blast radius BEFORE
any fix: is the divergence Fe II-only, or systemic across species?

READ-ONLY. Nothing on disk is modified. No gf is adopted, approximated, or
substituted here — this only measures the existing divergence.

Method
------
The two lists handle hyperfine structure (HFS) DIFFERENTLY: some lines are HFS-split
into many components on one side and carried as a single total-gf line on the other
(e.g. Na I / Sc II / La II — synth holds one line at the total gf, `linelist_solar`
holds N VALD3 components each at a fraction). Comparing a component gf to a total gf
is meaningless and would manufacture a huge spurious Δgf. So we compare the TOTAL
oscillator strength per physical line:

  1. AGGREGATE each side to physical lines: cluster same-`species_key` rows that
     share EP (|ΔEP| ≤ EPTOL) and sit within HFS_GAP Å of each other; the cluster's
     gf is log10(Σ 10^gf) (total) and its λ is the gf-weighted centroid. `n_comp`
     records how many components were summed (HFS multiplicity).
  2. MATCH physical lines across the two lists on:
        species : RYA-345 `species_key` (collapses 'Fe 2' / 'Fe','II' / int-ion /
                  molecule 'T'/'F' to one key) — same vocabulary on both sides.
        λ       : air Å, |Δλ| ≤ WTOL.
        EP      : |ΔEP| ≤ EPTOL — separates distinct transitions sharing a λ
                  (e.g. the two Fe II components at 6247.557 / 6247.576).
     gf (total) is then the only quantity that may differ. A pair where either side
     summed >1 component (`hfs` True) is reported separately — even with total-gf
     aggregation, one-sided splitting leaves more matching uncertainty there, so the
     clean 1:1 (`n_comp`==1 both sides) set is the trustworthy headline.

Anchor self-check: the script asserts it reproduces the five RYA-347 Fe II Δgf
values (a wrong match rule would not), so the scope numbers are trustworthy.

Outputs (under data/audit/gf_duplication/):
    gf_divergence_matched.csv   every matched pair + Δgf + flags
    gf_divergence_by_species.csv per-(species) summary
Console: |Δgf| band counts, per-species table, catalog sign pattern, verdict.

Usage:
    python3 scripts/audit_gf_duplication.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from pipeline.species import species_key, z_symbol, MOLECULE  # noqa: E402
import pipeline.abundances_derive as ad  # noqa: E402

# ── match tolerances ─────────────────────────────────────────────────────────
WTOL = 0.02     # Å — physical-line match window (centroids after HFS aggregation)
EPTOL = 0.02    # eV — same lower level; separates distinct transitions at one λ
HFS_GAP = 0.10  # Å — max gap between HFS components grouped into one physical line
BANDS = (0.01, 0.05, 0.10)

_SYNTH = ad._SYNTH_LINELIST_FILE
_LL_SOLAR = _REPO / 'data' / 'linelists' / 'linelist_solar.csv'
_OUT = _REPO / 'data' / 'audit' / 'gf_duplication'

# RYA-347 anchors: ll_solar VALD3 gf, synth GES v6 gf, expected Δ (synth − ll).
_ANCHORS = {
    5234.623: (-2.23, -2.180, +0.050),
    5991.371: (-3.54, -3.647, -0.107),
    6084.102: (-3.78, -3.881, -0.101),
    6247.557: (-2.31, -2.435, -0.125),
    6456.380: (-2.10, -2.185, -0.085),
}


def _species_label(key) -> str:
    if key[0] == MOLECULE:
        return key[1]
    roman = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI'}
    return f"{z_symbol(key[0])} {roman.get(key[1], key[1])}"


def _load():
    syn = pd.read_csv(_SYNTH, sep='\t', low_memory=False)
    ll = pd.read_csv(_LL_SOLAR, low_memory=False)

    syn['key'] = syn.apply(
        lambda r: species_key(r['element'], r['ion'], r['molecule']), axis=1)
    ll['key'] = ll.apply(
        lambda r: species_key(r['element'], r['ion']), axis=1)

    syn = syn.rename(columns={'wave_A': 'wl', 'loggf': 'gf_synth',
                              'lower_state_eV': 'ep_synth'})
    ll = ll.rename(columns={'wavelength_air_A': 'wl', 'log_gf': 'gf_ll',
                            'excitation_potential_eV': 'ep_ll'})
    return (syn[['key', 'wl', 'gf_synth', 'ep_synth', 'reference_code']],
            ll[['key', 'wl', 'gf_ll', 'ep_ll', 'nist_grade', 'loggf_source']])


def _aggregate(df: pd.DataFrame, gf_col: str, ep_col: str, extra: dict) -> pd.DataFrame:
    """Collapse HFS components to physical lines: same key, |ΔEP|≤EPTOL, λ-gap
    ≤HFS_GAP → one line carrying total gf = log10(Σ10^gf) at the gf-weighted
    centroid. `extra` maps output-col → source-col carried from the strongest
    component (for provenance, e.g. nist_grade / source)."""
    # Cluster via the ONE shared physical-line grouper (RYA-353) so the audit, the
    # canonical build, and the synth reroute all agree. EP-first + wl-gap; order-robust.
    from pipeline.gf_resolver import cluster_physical_lines
    d = df.reset_index(drop=True)
    keys = d['key'].tolist()
    wls = d['wl'].to_numpy(float)
    eps = d[ep_col].to_numpy(float)
    gfs = d[gf_col].to_numpy(float)
    out = []
    for cl in cluster_physical_lines(keys, wls, eps):
        w = 10.0 ** gfs[cl]
        rec = {'key': keys[cl[0]],
               'wl': float((wls[cl] * w).sum() / w.sum()),
               'gf': float(np.log10(w.sum())),
               'ep': float(eps[cl].mean()), 'n_comp': len(cl)}
        strongest = cl[int(np.argmax(gfs[cl]))]
        rec.update({k: d.iloc[strongest][v] for k, v in extra.items()})
        out.append(rec)
    return pd.DataFrame(out)


def _match(syn: pd.DataFrame, ll: pd.DataFrame) -> pd.DataFrame:
    """Aggregate each side to physical lines (total gf), then match physical line
    to physical line on (key, λ, EP). Read-only; pure measurement."""
    sa = _aggregate(syn, 'gf_synth', 'ep_synth', {'synth_ref': 'reference_code'})
    la = _aggregate(ll, 'gf_ll', 'ep_ll',
                    {'nist_grade': 'nist_grade', 'll_source': 'loggf_source'})
    rows = []
    common = sorted(set(sa['key']) & set(la['key']), key=str)
    for key in common:
        s = sa[sa['key'] == key].sort_values('wl').reset_index(drop=True)
        l = la[la['key'] == key].sort_values('wl').reset_index(drop=True)
        m = pd.merge_asof(l, s, on='wl', direction='nearest', tolerance=WTOL,
                          suffixes=('_l', '_s'))
        m = m[m['gf_s'].notna() if 'gf_s' in m else m['gf'].notna()]
        m = m.rename(columns={'gf_l': 'gf_ll', 'gf_s': 'gf_synth',
                              'ep_l': 'ep', 'n_comp_l': 'n_comp_ll',
                              'n_comp_s': 'n_comp_synth'})
        m = m[(m['ep'] - m['ep_s']).abs() <= EPTOL]
        if m.empty:
            continue
        for _, r in m.iterrows():
            rows.append({
                'key': key, 'species': _species_label(key), 'wl': r['wl'],
                'gf_ll': r['gf_ll'], 'gf_synth': r['gf_synth'],
                'dgf': r['gf_synth'] - r['gf_ll'], 'ep': r['ep'],
                'n_comp_ll': int(r['n_comp_ll']), 'n_comp_synth': int(r['n_comp_synth']),
                'nist_grade': r['nist_grade'], 'll_source': r['ll_source'],
                'synth_ref': r['synth_ref'],
                'hfs': bool(r['n_comp_ll'] > 1 or r['n_comp_synth'] > 1),
            })
    return pd.DataFrame(rows)


def _anchor_selfcheck(matched: pd.DataFrame) -> None:
    print("\n── anchor self-check (RYA-347 Fe II) ──────────────────────────────")
    ok = True
    fe2 = matched[matched['species'] == 'Fe II']
    for wl, (gf_ll, gf_syn, dexp) in _ANCHORS.items():
        hit = fe2[fe2['wl'].sub(wl).abs() < 0.01]
        if hit.empty:
            print(f"  {wl:9.3f}  MISSING from match  ✗"); ok = False; continue
        r = hit.iloc[0]
        good = (abs(r['gf_ll'] - gf_ll) < 1e-6 and abs(r['gf_synth'] - gf_syn) < 1e-6
                and abs(r['dgf'] - dexp) < 1e-6)
        ok &= good
        print(f"  {wl:9.3f}  ll={r['gf_ll']:+.3f} synth={r['gf_synth']:+.3f} "
              f"Δ={r['dgf']:+.3f} (exp {dexp:+.3f})  {'✓' if good else '✗'}")
    if not ok:
        raise SystemExit("Anchor self-check FAILED — match rule is wrong; "
                         "scope numbers untrustworthy. STOP.")
    print("  all 5 anchors reproduced ✓")


def _report(matched: pd.DataFrame) -> None:
    clean = matched[~matched['hfs']]
    print("\n── match summary ──────────────────────────────────────────────────")
    print(f"  matched physical lines        : {len(matched)}")
    print(f"    clean 1:1 (n_comp==1 both)  : {len(clean)}")
    print(f"    HFS-aggregated (flagged)    : {matched['hfs'].sum()}")

    def bands(df, label):
        n = len(df)
        print(f"\n  |Δgf| bands — {label} (n={n}):")
        for b in BANDS:
            c = (df['dgf'].abs() > b).sum()
            print(f"    > {b:.2f} dex : {c:6d}  ({100*c/n:5.1f}%)" if n else "")
        print(f"    median |Δgf| : {df['dgf'].abs().median():.4f}   "
              f"median Δgf : {df['dgf'].median():+.4f}   "
              f"(synth−ll; +%={100*(df['dgf']>0).mean():.0f})")
    bands(matched, "ALL matched")
    bands(clean, "CLEAN 1:1 only")

    # per-species table (clean only, sorted by # divergent >0.05)
    print("\n── per-species divergence (CLEAN 1:1 matches) ─────────────────────")
    print(f"  {'species':<8} {'n':>6} {'>0.01':>6} {'>0.05':>6} {'>0.10':>6} "
          f"{'medΔ':>8} {'sign':>5}")
    g = []
    for sp, df in clean.groupby('species'):
        g.append((sp, len(df),
                  (df['dgf'].abs() > 0.01).sum(),
                  (df['dgf'].abs() > 0.05).sum(),
                  (df['dgf'].abs() > 0.10).sum(),
                  df['dgf'].median()))
    for sp, n, c1, c5, c10, med in sorted(g, key=lambda x: -x[3]):
        sign = '+' if med > 0.001 else ('-' if med < -0.001 else '0')
        flag = '  <--' if c5 else ''
        print(f"  {sp:<8} {n:>6} {c1:>6} {c5:>6} {c10:>6} {med:>+8.3f} {sign:>5}{flag}")

    # catalog sign pattern
    print("\n── catalog sign pattern (CLEAN, |Δgf|>0.05) ───────────────────────")
    big = clean[clean['dgf'].abs() > 0.05]
    print(f"  divergent lines: {len(big)}")
    print(f"    synth stronger (Δ>0): {(big['dgf']>0).sum()}")
    print(f"    synth weaker   (Δ<0): {(big['dgf']<0).sum()}")
    if len(big):
        print(f"    mean Δgf among divergent: {big['dgf'].mean():+.4f}")

    # verdict
    fe2_big = (clean[(clean['species'] == 'Fe II') & (clean['dgf'].abs() > 0.05)])
    other_big = big[big['species'] != 'Fe II']
    print("\n── VERDICT ────────────────────────────────────────────────────────")
    print(f"  Fe II divergent (>0.05) : {len(fe2_big)}")
    print(f"  non-Fe II divergent     : {len(other_big)} "
          f"across {other_big['species'].nunique()} species")
    if len(other_big) > len(fe2_big):
        print("  => SYSTEMIC: divergence is broader than Fe II. PAUSE for scope "
              "decision before any pipeline-wide policy.")
    elif len(other_big):
        print("  => MOSTLY Fe II but non-trivial spillover; report and decide.")
    else:
        print("  => Fe II-ONLY: local edit suffices.")


def main() -> None:
    print(f"synth : {_SYNTH}")
    print(f"ll    : {_LL_SOLAR}")
    syn, ll = _load()
    print(f"loaded synth={len(syn)} ll_solar={len(ll)}  "
          f"(common species keys={len(set(syn['key']) & set(ll['key']))})")
    matched = _match(syn, ll)
    _anchor_selfcheck(matched)
    _report(matched)

    _OUT.mkdir(parents=True, exist_ok=True)
    matched.sort_values(['species', 'wl']).to_csv(
        _OUT / 'gf_divergence_matched.csv', index=False)
    by_sp = (matched.groupby('species')
             .agg(n=('dgf', 'size'),
                  n_gt01=('dgf', lambda x: (x.abs() > 0.01).sum()),
                  n_gt05=('dgf', lambda x: (x.abs() > 0.05).sum()),
                  n_gt10=('dgf', lambda x: (x.abs() > 0.10).sum()),
                  median_dgf=('dgf', 'median'),
                  mean_dgf=('dgf', 'mean'))
             .sort_values('n_gt05', ascending=False))
    by_sp.to_csv(_OUT / 'gf_divergence_by_species.csv')
    print(f"\nwrote {_OUT/'gf_divergence_matched.csv'}")
    print(f"wrote {_OUT/'gf_divergence_by_species.csv'}")


if __name__ == '__main__':
    main()
