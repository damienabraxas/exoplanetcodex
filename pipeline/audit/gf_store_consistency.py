"""
pipeline/audit/gf_store_consistency.py
======================================
RYA-368 — all-stores gf consistency audit.

RYA-353 single-sourced gf: the canonical table (`data/linelists/canonical_gf.csv`)
is the one authoritative value per physical line; each gf store resolves against it.
This audit reports EVERY gf store's overlap with canonical and its RAW divergence,
and asserts the two lines that triggered RYA-367 ([O I] 6300.304, Ni I 6300.34)
resolve to canonical (−9.717 / −2.11).

The stores (RYA-350/353):
  * #1 synth      atomic_lines.tsv (GES v6)          resolve-at-load: apply_to_synth_array
  * #3 regions    iSpec 42000_VALD ...abundances...  resolve-at-load: apply_to_regions
  * #2 diagnostic linelist_solar.csv (VALD3)         per-line resolve in consumers (RYA-368)

KEY FINDING (RYA-368): store-#2 is an INDEPENDENT VALD3 extraction and does NOT map
cleanly to the GES-based canonical table line-by-line — ~6.7k physical-line clusters
carry a different gf, 247 of them by >1 dex (up to 6.7). A bulk reroute/regeneration
to canonical would CORRUPT the store. Its gf is therefore treated as a raw VALD3
reference: it is NOT load-bearing for abundances (those resolve through store #3),
and the lone gf CONSUMER (the COG A(Fe) diagnostic) resolves per-line via the
resolver. The hard guarantee enforced here is: every store line has a canonical home
(no orphans), and the resolver returns canonical for the science-critical lines.

Smoke test:  python -m pipeline.audit.gf_store_consistency
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd

from pipeline import gf_resolver as gr
from pipeline.species import species_key

_REPO = Path(__file__).resolve().parents[2]
_CANON = _REPO / 'data' / 'linelists' / 'canonical_gf.csv'

# materiality threshold — same as the RYA-355 stewardship gf invariant
GF_DIVERGENCE_DEX = 0.05

# Science-critical target lines (the RYA-367 trigger): must resolve to canonical.
TARGETS = [
    ('[O I] 6300.304', (8, 1), 6300.304, 0.0, -9.717),
    ('Ni I 6300.34',  (28, 1), 6300.337, 4.266, -2.11),
]


@dataclass
class GfStore:
    label: str
    path: Path
    kind: str                # 'csv' | 'atomic_tsv' | 'regions'
    contract: str            # how the store reaches canonical at load
    is_load_bearing_gf: bool # True if its gf feeds an abundance/synthesis result
    ticket: Optional[str] = None   # remediation/tracking ticket for raw divergence


def _load_store_lines(store: GfStore) -> list[tuple]:
    """Return [(key, wl, ep, gf), ...] physical-component rows for a store."""
    rows = []
    if store.kind == 'csv':
        df = pd.read_csv(store.path, comment='#', low_memory=False)
        for el, ion, wl, ep, gf in zip(df['element'], df['ion'],
                                       df['wavelength_air_A'],
                                       df['excitation_potential_eV'], df['log_gf']):
            try:
                rows.append((species_key(str(el), str(ion)), float(wl), float(ep), float(gf)))
            except (ValueError, TypeError):
                continue
    elif store.kind == 'atomic_tsv':
        from pipeline.abundances_derive import ispec
        arr = ispec.read_atomic_linelist(str(store.path))
        for r in arr:
            try:
                rows.append((species_key(r['element'], r['ion'], r['molecule']),
                             float(r['wave_A']), float(r['lower_state_eV']),
                             float(r['loggf'])))
            except (ValueError, TypeError):
                continue
    elif store.kind == 'regions':
        from pipeline.abundances_derive import ispec
        reg = ispec.read_line_regions(str(store.path))
        names = reg.dtype.names or ()
        if 'loggf' not in names:
            return rows
        for r in reg:
            try:
                rows.append((species_key(str(r['note'])), float(r['wave_A']),
                             float(r['lower_state_eV']), float(r['loggf'])))
            except (ValueError, TypeError):
                continue
    else:
        raise ValueError(f"unknown store kind {store.kind!r}")
    return rows


def store_report(store: GfStore) -> dict:
    """Aggregate a store to physical lines, resolve each against canonical, and
    report overlap / raw-divergence / orphans (lines with no canonical home)."""
    rows = _load_store_lines(store)
    if not rows:
        return dict(store=store, n_lines=0, n_overlap=0, n_orphan=0,
                    n_divergent=0, max_abs_dgf=0.0, orphans=[], divergent=[])
    keys = [r[0] for r in rows]
    wls = np.array([r[1] for r in rows])
    eps = np.array([r[2] for r in rows])
    gfs = np.array([r[3] for r in rows])
    n_overlap = n_orphan = n_div = 0
    max_abs = 0.0
    orphans, divergent = [], []
    for cl in gr.cluster_physical_lines(keys, wls, eps):
        comp = gfs[cl]
        raw_total = float(np.log10(np.sum(10.0 ** comp)))
        w = 10.0 ** comp
        centroid = float(np.sum(wls[cl] * w) / w.sum())
        ep_mean = float(np.mean(eps[cl]))
        try:
            canon = gr.resolve(keys[cl[0]], centroid, ep_mean)
        except gr.GfResolutionError:
            n_orphan += 1
            if len(orphans) < 20:
                orphans.append((keys[cl[0]], centroid, ep_mean, raw_total))
            continue
        n_overlap += 1
        dgf = raw_total - canon
        if abs(dgf) > max_abs:
            max_abs = abs(dgf)
        if abs(dgf) > GF_DIVERGENCE_DEX:
            n_div += 1
            if len(divergent) < 20:
                divergent.append((keys[cl[0]], centroid, raw_total, canon, dgf))
    return dict(store=store, n_lines=len(rows),
                n_overlap=n_overlap, n_orphan=n_orphan, n_divergent=n_div,
                max_abs_dgf=max_abs, orphans=orphans, divergent=divergent)


def check_targets() -> list[dict]:
    """Resolve the RYA-367 trigger lines against canonical; each must hit its value."""
    out = []
    for label, key, wl, ep, expect in TARGETS:
        got = gr.resolve(key, wl, ep)
        out.append({'label': label, 'expect': expect, 'got': got,
                    'ok': abs(got - expect) < 0.01})
    return out


def run(verbose: bool = True) -> dict:
    """Print the all-stores report + the target resolution; assert no orphans and
    that the trigger lines resolve to canonical. Returns findings."""
    reports = [store_report(s) for s in STORES]
    targets = check_targets()

    if verbose:
        print("=" * 84)
        print("  RYA-368 — all-stores gf consistency audit (canonical = canonical_gf.csv)")
        print("=" * 84)
        print(f"\n  {'store':<22}{'kind':<11}{'lines':>8}{'overlap':>9}{'orphan':>8}"
              f"{'raw_div':>9}{'max|Δ|':>9}  contract")
        for r in reports:
            s = r['store']
            print(f"  {s.label:<22}{s.kind:<11}{r['n_lines']:>8}{r['n_overlap']:>9}"
                  f"{r['n_orphan']:>8}{r['n_divergent']:>9}{r['max_abs_dgf']:>9.2f}"
                  f"  {s.contract}")
        print("\n  TARGET LINES (must resolve to canonical):")
        for t in targets:
            print(f"    {t['label']:<16} resolver={t['got']:+.3f}  "
                  f"expect={t['expect']:+.3f}  {'OK' if t['ok'] else 'FAIL'}")

        print("\n  store-#2 (linelist_solar.csv) raw-divergence note:")
        s2 = next(r for r in reports if r['store'].label == 'linelist_solar.csv')
        print(f"    {s2['n_divergent']} physical lines diverge >{GF_DIVERGENCE_DEX} dex "
              f"(max {s2['max_abs_dgf']:.2f}) — INDEPENDENT VALD3, not line-by-line "
              f"reconcilable to GES canonical.")
        print(f"    gf is NOT load-bearing (abundances resolve via store #3); the lone "
              f"gf consumer (COG A(Fe) diagnostic) resolves per-line (RYA-368).")

    # ── hard assertions (the durable guarantees) ───────────────────────────────
    problems = []
    for r in reports:
        if r['n_orphan'] > 0:
            problems.append(f"{r['store'].label}: {r['n_orphan']} orphan line(s) with "
                            f"no canonical home (e.g. {r['orphans'][:3]})")
    for t in targets:
        if not t['ok']:
            problems.append(f"target {t['label']} resolves to {t['got']:+.3f}, "
                            f"expected {t['expect']:+.3f}")
    if verbose:
        print("\n" + "-" * 84)
        if problems:
            print("  RESULT: FAIL")
            for p in problems:
                print(f"    ! {p}")
        else:
            print("  RESULT: PASS — every store line resolves to a canonical home; "
                  "trigger lines resolve to canonical.")
        print("=" * 84)
    if problems:
        raise AssertionError("; ".join(problems))
    return {'reports': reports, 'targets': targets, 'pass': True}


# ── store registry (adding a store is a one-line append) ──────────────────────
def _build_stores() -> list[GfStore]:
    from pipeline.abundances_derive import (
        _SYNTH_LINELIST_FILE, _LINE_REGIONS_ALL)
    from config.constants import PATHS
    return [
        GfStore('atomic_lines.tsv', Path(_SYNTH_LINELIST_FILE), 'atomic_tsv',
                'resolve-at-load (apply_to_synth_array, RYA-353)',
                is_load_bearing_gf=True, ticket='RYA-353'),
        GfStore('iSpec line regions', Path(_LINE_REGIONS_ALL), 'regions',
                'resolve-at-load (apply_to_regions, RYA-353)',
                is_load_bearing_gf=True, ticket='RYA-353'),
        GfStore('linelist_solar.csv', Path(PATHS['linelist_solar']), 'csv',
                'raw VALD3 ref; per-line resolve in consumers (RYA-368)',
                is_load_bearing_gf=False, ticket='RYA-368'),
    ]


STORES: list[GfStore] = _build_stores()


if __name__ == '__main__':
    run()
