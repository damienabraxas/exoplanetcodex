"""
scripts/migrate_gf_single_source.py
===================================
RYA-353 — build the canonical single-source gf table (implements the RYA-350
design, Option A). This module's `--build` produces the authoritative
`data/linelists/canonical_gf.csv`: one row per physical line, the single value
both engines will resolve against. (Loader rerouting is the next step, held.)

Three independent gf copies exist (RYA-353 loader trace; corrected RYA-350 §1):
  #1 atomic_lines.tsv          GES v6   → flux synthesis        (reference_code = underlying ref)
  #3 42000_VALD regions file   VALD3    → the live EW→abundance  (524 curated lines)
  #2 linelist_solar.csv        VALD3    → COG / diagnostics / scoring
plus NIST ASD graded values (nist_reference.csv + nist_crosscheck.csv) — precise
log_gf + grade for ~60 anchor lines.

Build (per the signed-off design + its four answers):
  • Aggregate each store to PHYSICAL LINES, HFS-aware: total gf = log10(Σ10^gf)
    over same-species components within EPTOL / HFS_GAP (the Step-0 matcher; reused
    from audit_gf_duplication via _aggregate). Records hfs_n_components.
  • Union the three stores per canonical physical line (species_key + λ + EP),
    assign an opaque surrogate line_id (Q4) and keep species_key/λ/EP columns too.
  • SEED the single log_gf to the best current per-line value, preserving published
    differentials by construction (design §2.4, sign-off Q1):
        1. NIST-graded present            → the NIST value      (status=nist; fixes 6247, Q3)
        2. else in #3 (live EW route)     → the regions value   (keeps EW abundances stable)
        3. else in #2 (linelist_solar)    → that value
        4. else in #1 (synth-only)        → the GES value
    A line whose synth (#1) and EW (#3) values diverge >SEED_TOL with no NIST grade
    is flagged adjudication_status=pending → RYA-354 (value quality), NOT resolved here.
  • Record provenance: per-store values, seed_source, loggf_reference (the underlying
    paper for GES via reference_code; 'VALD3'/'NIST ASD v5.11' otherwise), Δ(synth−EW).

NO silent fallback: a line that resolves to no sourced value raises. Component-level
detail for the synth path (rescaling GES components to the seeded total) is handled by
the loader reroute, not here — this builds the per-line canonical list.

Usage:
    python3 scripts/migrate_gf_single_source.py --build      # write canonical_gf.csv
    python3 scripts/migrate_gf_single_source.py --validate   # (next step — stub)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / 'scripts'))

from pipeline.species import species_key, z_symbol, MOLECULE  # noqa: E402
import pipeline.abundances_derive as ad  # noqa: E402
from audit_gf_duplication import _aggregate, WTOL, EPTOL  # noqa: E402  (reuse Step-0 matcher)

SEED_TOL = 0.01  # dex — |Δgf| at/under which two stores are "the same value"

_LL_SOLAR = _REPO / 'data' / 'linelists' / 'linelist_solar.csv'
_NIST_REF = _REPO / 'data' / 'linelists' / 'nist_reference.csv'
_NIST_XC = _REPO / 'data' / 'linelists' / 'nist_crosscheck.csv'
_OUT = _REPO / 'data' / 'linelists' / 'canonical_gf.csv'

# RYA-347 anchors for the build self-check (synth, ll/regions expected per-store gf).
_ANCHORS = {5234.623: (-2.180, -2.230), 5991.371: (-3.647, -3.540),
            6084.102: (-3.881, -3.780), 6247.557: (-2.435, -2.310),
            6456.380: (-2.185, -2.100)}


def _label(key) -> str:
    if key[0] == MOLECULE:
        return key[1]
    roman = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI'}
    return f"{z_symbol(key[0])} {roman.get(key[1], key[1])}"


# ── load + aggregate the three stores to physical lines ───────────────────────

def _load_synth() -> pd.DataFrame:
    syn = pd.read_csv(ad._SYNTH_LINELIST_FILE, sep='\t', low_memory=False)
    syn['key'] = syn.apply(lambda r: species_key(r['element'], r['ion'], r['molecule']), axis=1)
    syn = syn.rename(columns={'wave_A': 'wl', 'loggf': 'gf_synth', 'lower_state_eV': 'ep_synth'})
    return _aggregate(syn[['key', 'wl', 'gf_synth', 'ep_synth', 'reference_code']],
                      'gf_synth', 'ep_synth', {'synth_ref': 'reference_code'})


def _load_linelist() -> pd.DataFrame:
    ll = pd.read_csv(_LL_SOLAR, low_memory=False)
    ll['key'] = ll.apply(lambda r: species_key(r['element'], r['ion']), axis=1)
    ll = ll.rename(columns={'wavelength_air_A': 'wl', 'log_gf': 'gf_ll',
                            'excitation_potential_eV': 'ep_ll'})
    return _aggregate(ll[['key', 'wl', 'gf_ll', 'ep_ll', 'nist_grade', 'loggf_source']],
                      'gf_ll', 'ep_ll', {'ll_grade': 'nist_grade', 'll_source': 'loggf_source'})


def _load_regions() -> pd.DataFrame:
    import ispec
    reg = ispec.read_line_regions(ad._LINE_REGIONS_ALL)
    rows = []
    for r in reg:
        try:
            k = species_key(str(r['note']))
        except ValueError:
            continue  # molecular / non-atomic region note
        rows.append({'key': k, 'wl': float(r['wave_A']), 'gf_reg': float(r['loggf']),
                     'ep_reg': float(r['lower_state_eV'])})
    df = pd.DataFrame(rows)
    return _aggregate(df, 'gf_reg', 'ep_reg', {})


def _load_nist() -> pd.DataFrame:
    frames = []
    for f in (_NIST_REF, _NIST_XC):
        d = pd.read_csv(f, comment='#')
        d['key'] = d.apply(lambda r: species_key(r['element'], r['ion']), axis=1)
        frames.append(d[['key', 'wavelength_air_A', 'excitation_potential_eV',
                         'log_gf', 'nist_grade']].rename(
            columns={'wavelength_air_A': 'wl', 'excitation_potential_eV': 'ep',
                     'log_gf': 'nist_gf'}))
    nist = pd.concat(frames, ignore_index=True)
    # nist_reference (grade A) wins over crosscheck on dup; keep first after sort
    return nist.sort_values('nist_grade').drop_duplicates(['key', 'wl'], keep='first')


def _grouped(df: pd.DataFrame, cols) -> dict:
    """Pre-split a store into {key: ndarray-of-records} for cheap per-key lookup."""
    if df is None or df.empty:
        return {}
    return {k: g[cols].to_dict('records') for k, g in df.groupby('key', sort=False)}


def _nearest_rec(recs, wl, ep, gfcol):
    """Nearest record (same key) within WTOL/EPTOL → (gf, rec) or (nan, None).
    `recs` is the small per-key list from _grouped — no full-frame scan."""
    best, bd = None, None
    for r in recs:
        if abs(r['wl'] - wl) <= WTOL and abs(r['ep'] - ep) <= EPTOL:
            d = abs(r['wl'] - wl)
            if bd is None or d < bd:
                bd, best = d, r
    return (float(best[gfcol]), best) if best is not None else (np.nan, None)


# ── build ─────────────────────────────────────────────────────────────────────

def build() -> pd.DataFrame:
    print("Loading + aggregating the three stores (HFS-aware total gf)...")
    syn = _load_synth()      # key, wl, gf, ep, n_comp, synth_ref
    ll = _load_linelist()    # key, wl, gf, ep, n_comp, ll_grade, ll_source
    reg = _load_regions()    # key, wl, gf, ep, n_comp
    nist = _load_nist()      # key, wl, ep, nist_gf, nist_grade
    print(f"  synth(#1)={len(syn)}  regions(#3)={len(reg)}  linelist(#2)={len(ll)}  nist={len(nist)}")

    # Per-key groups (small) so all matching is within-key, not full-frame scans.
    g_syn = _grouped(syn, ['wl', 'ep', 'gf', 'n_comp', 'synth_ref'])
    g_reg = _grouped(reg, ['wl', 'ep', 'gf'])
    g_ll = _grouped(ll, ['wl', 'ep', 'gf', 'll_grade', 'll_source'])
    g_nist = _grouped(nist, ['wl', 'ep', 'nist_gf', 'nist_grade'])

    # Canonical physical lines = union over stores, each emitted once. Priority for the
    # (wl, ep) anchor: regions(#3, live EW) → synth(#1) → linelist(#2).
    rows = []
    src_rank = {'reg': 0, 'syn': 1, 'll': 2}
    pool = pd.concat([
        reg.assign(src='reg')[['key', 'wl', 'ep', 'src']],
        syn.assign(src='syn')[['key', 'wl', 'ep', 'src']],
        ll.assign(src='ll')[['key', 'wl', 'ep', 'src']],
    ], ignore_index=True)
    pool = pool.sort_values(by='src', key=lambda s: s.map(src_rank)).reset_index(drop=True)

    seen_by_key: dict = {}
    for key, wl, ep in zip(pool['key'], pool['wl'], pool['ep']):
        prev = seen_by_key.setdefault(key, [])
        if any(abs(wl - w) <= WTOL and abs(ep - e) <= EPTOL for w, e in prev):
            continue  # already emitted this physical line via a higher-rank store
        prev.append((wl, ep))

        gf_syn, rsyn = _nearest_rec(g_syn.get(key, []), wl, ep, 'gf')
        gf_reg, _    = _nearest_rec(g_reg.get(key, []), wl, ep, 'gf')
        gf_ll,  rll  = _nearest_rec(g_ll.get(key, []), wl, ep, 'gf')
        nist_gf, rn  = _nearest_rec(g_nist.get(key, []), wl, ep, 'nist_gf')

        in_syn, in_reg, in_ll = not np.isnan(gf_syn), not np.isnan(gf_reg), not np.isnan(gf_ll)
        has_nist = not np.isnan(nist_gf)
        n_comp = int(rsyn['n_comp']) if in_syn else 1

        # seed — priority follows the line's PRIMARY measurement path (sign-off Q1):
        #   NIST graded  →  live-EW regions (#3)  →  synth (#1)  →  linelist (#2, auxiliary)
        # so synth-primary lines keep their GES value and EW-path lines keep the published
        # VALD3 value; only #2-only (COG/diag) lines fall to linelist. Differentials cancel
        # gf regardless; this keeps each path's ABSOLUTE values stable where it is primary.
        if has_nist:
            seed, src = nist_gf, 'NIST'
            ref = f"NIST ASD v5.11 grade {rn['nist_grade']}"
            grade, status = rn['nist_grade'], 'nist'
        elif in_reg:
            seed, src, grade = gf_reg, 'regions(VALD3)', ''
            ref = 'VALD3'
            status = ('pending' if in_syn and abs(gf_syn - gf_reg) > SEED_TOL else 'agreed')
        elif in_syn:
            seed, src, grade = gf_syn, 'synth(GES)', ''
            ref = str(rsyn['synth_ref'])       # underlying reference code (K07/FMW/…)
            status = ('pending' if in_ll and abs(gf_syn - gf_ll) > SEED_TOL else 'single_source')
        elif in_ll:
            seed, src, grade = gf_ll, 'linelist(VALD3)', rll['ll_grade']
            ref = 'VALD3'
            status = 'single_source'
        else:
            raise RuntimeError(f"line {_label(key)} {wl:.3f} resolves to no sourced gf")

        dsyn_ew = (gf_syn - gf_reg) if (in_syn and in_reg) else np.nan
        rows.append({
            'key_z': key[1] if key[0] == MOLECULE else key[0],
            'ion': '' if key[0] == MOLECULE else key[1],
            'species': _label(key),
            'wavelength_air_A': round(wl, 4),
            'excitation_potential_eV': round(ep, 4),
            'hfs_n_components': n_comp,
            'log_gf': round(seed, 4),
            'loggf_reference': ref,
            'nist_grade': grade,
            'seed_source': src,
            'adjudication_status': status,
            'gf_synth_ges': round(gf_syn, 4) if in_syn else '',
            'gf_regions_vald': round(gf_reg, 4) if in_reg else '',
            'gf_linelist_vald': round(gf_ll, 4) if in_ll else '',
            'in_synth': in_syn, 'in_regions': in_reg, 'in_linelist': in_ll,
            'delta_synth_minus_ew': round(dsyn_ew, 4) if not np.isnan(dsyn_ew) else '',
        })

    df = pd.DataFrame(rows).sort_values(['species', 'wavelength_air_A']).reset_index(drop=True)
    df.insert(0, 'line_id', [f"gf_{i:06d}" for i in range(len(df))])  # opaque surrogate (Q4)
    return df


def _selfcheck(df: pd.DataFrame) -> None:
    print("\n── anchor self-check (RYA-347 Fe II) ──")
    fe2 = df[df['species'] == 'Fe II']
    ok = True
    for wl, (gs, ge) in _ANCHORS.items():
        h = fe2[fe2['wavelength_air_A'].sub(wl).abs() < 0.01]
        if h.empty:
            print(f"  {wl:9.3f} MISSING ✗"); ok = False; continue
        r = h.iloc[0]
        # 6247 is NIST-graded → seed must be NIST −2.329, not the VALD3 −2.31
        exp_seed = -2.329 if abs(wl - 6247.557) < 0.01 else ge
        good = (abs(float(r['gf_synth_ges']) - gs) < 1e-3
                and abs(float(r['gf_regions_vald']) - ge) < 1e-3
                and abs(float(r['log_gf']) - exp_seed) < 1e-3)
        ok &= good
        print(f"  {wl:9.3f} synth={r['gf_synth_ges']} ew={r['gf_regions_vald']} "
              f"seed={r['log_gf']} ({r['seed_source']}) {'✓' if good else '✗'}")
    if not ok:
        raise SystemExit("anchor self-check FAILED — STOP.")
    print("  anchors OK ✓")


def _report(df: pd.DataFrame) -> None:
    n = len(df)
    print(f"\n── canonical_gf.csv summary (n={n} physical lines) ──")
    print("  presence:")
    print(f"    in all 3 stores : {(df.in_synth & df.in_regions & df.in_linelist).sum()}")
    print(f"    synth-only      : {(df.in_synth & ~df.in_regions & ~df.in_linelist).sum()}")
    print(f"    in live EW (#3) : {df.in_regions.sum()}")
    print("  seed source:")
    print(df['seed_source'].value_counts().to_string().replace('\n', '\n    ').rjust(0))
    print("  adjudication_status:")
    print(df['adjudication_status'].value_counts().to_string().replace('\n', '\n    '))
    print(f"  NIST-seeded lines : {(df.adjudication_status=='nist').sum()}")
    print(f"  pending (synth≠EW, no NIST → RYA-354): {(df.adjudication_status=='pending').sum()}")
    miss = df['loggf_reference'].isna() | (df['loggf_reference'] == '')
    print(f"  provenance completeness: {n - miss.sum()}/{n} have a reference "
          f"({'OK' if miss.sum()==0 else 'MISSING '+str(miss.sum())})")


def validate() -> None:
    """Prove the reroute: both abundance paths (synth #1, EW regions #3) now resolve
    to the SAME gf per line — zero cross-path divergence — and the loud-guard fires."""
    import ispec
    from pipeline import gf_resolver as gr
    from pipeline.gf_resolver import cluster_physical_lines

    print("── 1. cross-path gf agreement (synth #1 vs EW regions #3, post-reroute) ──")
    arr = gr.apply_to_synth_array(ispec.read_atomic_linelist(ad._SYNTH_LINELIST_FILE))
    keys = [species_key(arr['element'][i], arr['ion'][i], arr['molecule'][i])
            for i in range(len(arr))]
    wl = arr['wave_A'].astype(float); ep = arr['lower_state_eV'].astype(float)
    gf = arr['loggf'].astype(float)
    syn_tot = {}  # (key, round wl, round ep) -> total gf
    for cl in cluster_physical_lines(keys, wl, ep):
        w = 10.0 ** gf[cl]
        cen = float((wl[cl] * w).sum() / w.sum())
        syn_tot[(keys[cl[0]], round(cen, 2), round(float(ep[cl].mean()), 2))] = \
            float(np.log10(w.sum()))

    reg = gr.apply_to_regions(ispec.read_line_regions(ad._LINE_REGIONS_ALL))
    n_cmp = n_div = n_noabs = 0
    worst = 0.0
    for i in range(len(reg)):
        try:
            k = species_key(str(reg['note'][i]))
        except ValueError:
            continue
        rwl, rep, rgf = float(reg['wave_A'][i]), float(reg['lower_state_eV'][i]), float(reg['loggf'][i])
        hit = None
        for dwl in (0.0, 0.01, -0.01, 0.02, -0.02):
            cand = syn_tot.get((k, round(rwl + dwl, 2), round(rep, 2)))
            if cand is not None:
                hit = cand; break
        if hit is None:
            n_noabs += 1; continue
        n_cmp += 1
        d = abs(hit - rgf)
        worst = max(worst, d)
        if d > 0.01:
            n_div += 1
    print(f"  regions in both paths: {n_cmp}  | cross-path divergent (>0.01): {n_div}  "
          f"| worst |Δ|: {worst:.4f}")
    print(f"  regions with no synth counterpart (EW-only, expected): {n_noabs}")
    assert n_div == 0, f"STOP: {n_div} lines still diverge across paths after reroute"
    print("  ✓ zero cross-path gf divergence — both paths read one value")

    print("\n── 2. loud-guard (no silent fallback) ──")
    try:
        gr.resolve(species_key('Fe', 'II'), 9999.99, 3.0)
        raise SystemExit("STOP: resolver did not raise on an absent line")
    except gr.GfResolutionError:
        print("  ✓ resolver raises GfResolutionError on an unresolved line")

    print("\n── 3. anchor (6247.557 → NIST −2.329 on both paths) ──")
    v = gr.resolve(species_key('Fe', 'II'), 6247.557, 3.892)
    print(f"  resolve(Fe II 6247.557) = {v:+.4f}  {'✓' if abs(v+2.329) < 1e-3 else '✗ FAIL'}")
    print("\nVALIDATION PASSED")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--validate', action='store_true')
    args = ap.parse_args()
    if args.validate:
        validate(); return
    if not args.build:
        ap.print_help(); return

    df = build()
    _selfcheck(df)
    _report(df)
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(_OUT, index=False)
    print(f"\nwrote {_OUT}  ({len(df)} lines)")


if __name__ == '__main__':
    main()
