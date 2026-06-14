# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         build_linelist_procyon.py
# Module:       scripts (linelist build tools)
# Description:  RYA-270 — build the Procyon F-star line pool
#               data/linelists/linelist_procyon.csv from the two RYA-223
#               VALD3 FTP extracts (3000–4800 Å blue + 3780–6910 Å optical).
#               Fe I/Fe II selection at Procyon's parameters via the VALD
#               central_depth computed for Teff=6530/logg=3.96/vmic=1.7.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-12
# Last modified: 2026-06-12
# Linear issue: RYA-270 — BUILD: Procyon F-star Fe I line pool
#
# -----------------------------------------------------------------------------
# KEY REFERENCES
# -----------------------------------------------------------------------------
# Ryabchikova et al. 2015 — Phys. Scr. 90, 054005 — VALD3 description
# Allende Prieto et al. 2002 — A&A 386, 1039 — Procyon parameters
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
# External: pandas, numpy, scipy (proximity flag)
# Data:     data/linelists/vald_procyon_raw/*.vald
#           data/linelists/linelist_solar.csv (GES/NIST carryover)
# =============================================================================

import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from config.constants import PATHS
from data.linelists.vald_parse import read_vald_header, parse_vald_long
from pipeline.build_linelist import compute_vald_proximity_flag

LINELISTS = PATHS['linelists']
OUT_PATH = LINELISTS / 'linelist_procyon.csv'

# RYA-223 extracts, labels per that issue's inspection report.
EXTRACTIONS = [
    ('optical_3780_6910', LINELISTS / 'vald_procyon_raw' / 'vald_procyon_optical.vald'),
    ('blue_3000_4800',    LINELISTS / 'vald_procyon_raw' / 'vald_procyon_blue.vald'),
]

# Fe selection window on VALD central_depth at the Procyon extraction params
# (RYA-270 spec): excludes noise-level lines (< 0.05) and saturated lines
# whose EWs are damping-dominated at F-star temperatures (> 0.80).
CD_MIN, CD_MAX = 0.05, 0.80

# Vetted Procyon exclusions established by RYA-249/263/268 (all Fe II).
# Matched by wavelength below; blend_flag=True with the citing note.
VETTED_PROCYON = [
    (5376.466, 'RYA-263'),
    (6149.258, 'RYA-249 vetting'),
    (6229.340, 'RYA-268'),
    (6247.557, 'RYA-249 vetting'),
    (6491.254, 'RYA-268'),
]
VETTING_TOL = 0.05  # Å

LOGGF_CONFLICT_TOL = 0.001
PARSE_FAIL_CRITICAL_FRAC = 0.001

CRITICALS = []


def critical(msg):
    CRITICALS.append(msg)
    print(f"\n  *** CRITICAL: {msg}")


def main():
    print("RYA-270 — Procyon F-star line pool → linelist_procyon.csv")

    # ── Parse extracts ───────────────────────────────────────────────────────
    frames = []
    fail_total, n_total = 0, 0
    for tag, path in EXTRACTIONS:
        hdr = read_vald_header(path)
        if hdr['truncated']:
            critical(f"{tag}: truncated VALD delivery")
            continue
        records, rep = parse_vald_long(path)
        fail_total += rep['n_failures']
        n_total += rep['n_parsed'] + rep['n_failures']
        df = pd.DataFrame(records)
        df['extraction_source'] = tag
        frames.append(df)
        print(f"  {tag:<18}: {hdr['wl_start']:.0f}–{hdr['wl_end']:.0f} Å, "
              f"parsed {rep['n_parsed']}/{hdr['n_selected']} selected, "
              f"failures {rep['n_failures']}, truncation: NO")
        if rep['examples']:
            for ex in rep['examples']:
                print(f"    parse failure: {ex}")
    if n_total and fail_total / n_total > PARSE_FAIL_CRITICAL_FRAC:
        critical(f"parse failure rate {fail_total/n_total:.2%} > 0.1%")
    if CRITICALS:
        return finish()

    merged = pd.concat(frames, ignore_index=True)

    # ── Cross-extraction overlap dedup (3780–4800 Å) ─────────────────────────
    # Same rules as RYA-269: dedup across extractions only (same-extraction
    # key-sharing = HFS/isotopic components, retained); keep the copy from the
    # more inclusive extraction. Both FTP requests used the same settings
    # (empirical min central_depth reported below); the optical extract is the
    # primary, so its copy wins on ties. log_gf disagreement → CRITICAL.
    for tag, _ in EXTRACTIONS:
        sub = merged[merged['extraction_source'] == tag]
        print(f"  {tag:<18}: min central_depth = {sub['central_depth'].min():.3f} "
              f"(empirical request-threshold proxy)")
    merged['dedup_key'] = list(zip(
        merged['element'], merged['ion'],
        merged['wavelength'].round(3), merged['e_low_eV'].round(3),
    ))
    pref = {tag: i for i, (tag, _) in enumerate(EXTRACTIONS)}  # optical first
    merged['src_pref'] = merged['extraction_source'].map(pref)
    key_best = merged.groupby('dedup_key')['src_pref'].transform('min')
    keep_mask = merged['src_pref'] == key_best
    kept, dropped = merged[keep_mask].copy(), merged[~keep_mask]

    conflicts = []
    kept_loggf = kept.groupby('dedup_key')['log_gf'].agg(list)
    for key, grp in dropped.groupby('dedup_key', sort=False):
        kept_vals = kept_loggf[key]
        for src, lg in zip(grp['extraction_source'], grp['log_gf']):
            if min(abs(lg - kv) for kv in kept_vals) > LOGGF_CONFLICT_TOL:
                conflicts.append((key, src, lg, kept_vals))
    if conflicts:
        critical(f"{len(conflicts)} cross-extraction duplicates disagree on "
                 f"log_gf > {LOGGF_CONFLICT_TOL} — refusing to pick silently")
        for c in conflicts[:5]:
            print(f"    {c}")
        return finish()
    print(f"  Overlap dedup: {len(merged)} → {len(kept)} "
          f"({len(dropped)} duplicates removed, 0 log_gf conflicts)")
    merged = kept

    # ── Proximity flag on the FULL all-species pool ──────────────────────────
    # A line's contamination comes from any nearby line, so the neighbour
    # density must be computed before subsetting to Fe (same formula as
    # linelist_solar — pipeline/build_linelist.py, RYA-209).
    prox_df = pd.DataFrame({'wavelength_air_A': merged['wavelength'],
                            'central_depth': merged['central_depth']})
    merged['vald_proximity_flag'] = compute_vald_proximity_flag(prox_df)

    # ── Fe I/Fe II selection ─────────────────────────────────────────────────
    is_fe = (merged['element'] == 'Fe') & merged['ion'].isin(['I', 'II'])
    in_window = merged['central_depth'].between(CD_MIN, CD_MAX)
    merged['priority'] = np.where(is_fe & in_window, 1, 0)

    fe_all = merged[is_fe]
    pool = merged[merged['priority'] > 0]
    fe1 = pool[pool['ion'] == 'I']
    fe2 = pool[pool['ion'] == 'II']
    print(f"\n  Fe candidates: {len(fe_all)} total "
          f"(Fe I {int((fe_all['ion']=='I').sum())}, Fe II {int((fe_all['ion']=='II').sum())})")
    print(f"  Selected (central_depth ∈ [{CD_MIN}, {CD_MAX}]): "
          f"{len(pool)} (Fe I {len(fe1)}, Fe II {len(fe2)})")
    print(f"  Rejected: {int((is_fe & ~in_window).sum())} "
          f"(too weak < {CD_MIN}: {int((is_fe & (merged['central_depth'] < CD_MIN)).sum())}, "
          f"saturated > {CD_MAX}: {int((is_fe & (merged['central_depth'] > CD_MAX)).sum())})")
    if len(fe2) < 5:
        critical(f"Fe II pool has {len(fe2)} lines after selection (< 5)")
        return finish()

    # EP coverage of the selected Fe I pool (excitation balance needs spread)
    ep = fe1['e_low_eV']
    bins = np.arange(0, 7.5, 1.0)
    hist, _ = np.histogram(ep, bins=bins)
    print(f"\n  Fe I EP histogram (selected pool):")
    for lo, n in zip(bins[:-1], hist):
        print(f"    {lo:.0f}–{lo+1:.0f} eV: {n:5d}  {'#' * min(n // 20, 40)}")
    print(f"  EP range: {ep.min():.3f}–{ep.max():.3f} eV, median {ep.median():.3f}")

    # Build-time quality components (full RYA-220 score needs measured EWs and
    # is computed at runtime by _compute_line_scores; distribution of the
    # linelist-side components documented here).
    for label, ser in [('vald_proximity_flag (selected Fe)', pool['vald_proximity_flag']),
                       ('central_depth (selected Fe)', pool['central_depth'])]:
        q = ser.quantile([0, .25, .5, .75, 1]).values
        print(f"  {label}: min={q[0]:.3f} q25={q[1]:.3f} med={q[2]:.3f} "
              f"q75={q[3]:.3f} max={q[4]:.3f}")

    # ── GES/NIST carryover from the solar list ───────────────────────────────
    sol = pd.read_csv(str(PATHS['linelist_solar']), low_memory=False,
                      usecols=['element', 'ion', 'wavelength_air_A',
                               'loggf_source', 'nist_grade'])
    sol_fe = sol[sol['element'] == 'Fe']
    merged['loggf_source'] = 'VALD3'
    merged['nist_grade'] = ''
    n_carried = 0
    fe_idx = merged.index[is_fe]
    for ion, sgrp in sol_fe.groupby('ion'):
        m_ion = merged.loc[fe_idx]
        m_ion = m_ion[m_ion['ion'] == ion]
        if m_ion.empty:
            continue
        sw = sgrp['wavelength_air_A'].values
        for idx, wl in zip(m_ion.index, m_ion['wavelength'].values):
            j = int(np.abs(sw - wl).argmin())
            if abs(sw[j] - wl) < VETTING_TOL:
                src = sgrp.iloc[j]['loggf_source']
                grd = sgrp.iloc[j]['nist_grade']
                if pd.notna(src):
                    merged.loc[idx, 'loggf_source'] = src
                if pd.notna(grd) and str(grd).strip():
                    merged.loc[idx, 'nist_grade'] = str(grd)
                    n_carried += 1
    n_solar_member = 0
    for ion, sgrp in sol_fe.groupby('ion'):
        m_ion = merged.loc[fe_idx]
        m_ion = m_ion[m_ion['ion'] == ion]
        sw = sgrp['wavelength_air_A'].values
        n_solar_member += sum(np.abs(sw - wl).min() < VETTING_TOL
                              for wl in m_ion['wavelength'].values)
    print(f"\n  Solar-list carryover: {n_solar_member} Fe lines also in "
          f"linelist_solar (loggf_source/nist_grade carried; "
          f"{n_carried} with a NIST grade)")

    # ── Vetted Procyon blend exclusions (RYA-249/263/268) ────────────────────
    merged['blend_flag'] = False
    merged['notes'] = ''
    for wl, issue in VETTED_PROCYON:
        m = (is_fe & (np.abs(merged['wavelength'] - wl) < VETTING_TOL))
        n_hit = int(m.sum())
        if n_hit:
            merged.loc[m, 'blend_flag'] = True
            merged.loc[m, 'notes'] = f'vetted Procyon exclusion ({issue})'
            print(f"  Vetted exclusion: {wl:.3f} Å — {n_hit} row(s) flagged ({issue})")
        else:
            print(f"  WARNING: vetted exclusion {wl:.3f} Å not matched in pool "
                  f"(not in cd window or extracts)")

    # ── Write output ─────────────────────────────────────────────────────────
    out = pd.DataFrame({
        'element'                 : merged['element'],
        'ion'                     : merged['ion'],
        'wavelength_air_A'        : merged['wavelength'],
        'excitation_potential_eV' : merged['e_low_eV'],
        'log_gf'                  : merged['log_gf'],
        'loggf_source'            : merged['loggf_source'],
        'nist_grade'              : merged['nist_grade'],
        'damping_rad'             : merged['damping_rad'],
        'damping_stark'           : merged['damping_stark'],
        'damping_vdW'             : merged['damping_vdW'],
        'central_depth'           : merged['central_depth'],
        'vald_proximity_flag'     : merged['vald_proximity_flag'].round(6),
        'blend_flag'              : merged['blend_flag'],
        'priority'                : merged['priority'],
        'notes'                   : merged['notes'],
        'extraction_source'       : merged['extraction_source'],
        'wavelength_convention'   : 'air',   # all lines ≥ 3000 Å
    }).sort_values('wavelength_air_A')

    header = f"""\
# linelist_procyon.csv — Procyon (HD 61421) F-star per-system linelist (RYA-270)
# Generated: {date.today().isoformat()}  by scripts/build_linelist_procyon.py
# Stellar params of the VALD extraction: Teff=6530 K, logg=3.96, [Fe/H]=-0.04,
#   vmic=1.7 km/s (Allende Prieto et al. 2002)
# Sources (RYA-223 FTP deliveries, both complete, HFS-unsplit convention):
#   optical_3780_6910: vald_procyon_raw/vald_procyon_optical.vald (primary)
#   blue_3000_4800   : vald_procyon_raw/vald_procyon_blue.vald
#   Overlap 3780-4800 A deduped cross-extraction, optical copy kept.
# priority=1 marks the Fe I/II working pool: central_depth in [{CD_MIN}, {CD_MAX}]
#   at Procyon params (excludes noise-level and damping-dominated lines).
#   All other species included at priority=0 for proximity/blend context.
# vald_proximity_flag computed on the full all-species pool (RYA-209 formula,
#   pipeline/build_linelist.py) — feeds the RYA-220 runtime quality scorer.
# blend_flag=True only for vetted Procyon exclusions (RYA-249/263/268).
# wavelength_convention: air throughout (no lines below 2000 A).
# nist_grade carried over from linelist_solar.csv where the line matches.
"""
    with open(OUT_PATH, 'w') as f:
        f.write(header)
        out.to_csv(f, index=False)
    rel = OUT_PATH.relative_to(REPO_ROOT)
    print(f"\n  Wrote {rel} ({len(out)} rows, "
          f"{len(fe1)} Fe I + {len(fe2)} Fe II at priority 1)")
    return finish()


def finish():
    if CRITICALS:
        print(f"\nRESULT: CRITICAL — {len(CRITICALS)} stop condition(s)")
        for c in CRITICALS:
            print(f"  - {c}")
        return 1
    print("\nRESULT: OK")
    return 0


if __name__ == '__main__':
    sys.exit(main())
