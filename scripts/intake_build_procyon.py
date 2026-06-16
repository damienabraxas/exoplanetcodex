#!/usr/bin/env python3
# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         intake_build_procyon.py
# Module:       scripts (VALD intake + per-system line-list build)
# Description:  RYA-309 Step 3.1 — intake gate + per-system rebuild of the
#               Procyon (HD 61421) VALD3 line list from the full HFS-ON
#               "Extract Stellar" delivery (all species).
#
#               The previous linelist_procyon.csv (RYA-270) was built from the
#               RYA-223 deliveries in the HFS-UNSPLIT convention — the exact HFS
#               gap RYA-270 flagged. This rebuilds it HFS-ON. Per the VALD gate,
#               HFS-ON is VERIFIED by the split-group signature, not trusted
#               from the submission form ("should have been on" is the
#               assumption that hid the original HFS bug).
#
#               Intake per delivery: truncation (line 1 = header) → parse
#               completeness → coverage → HFS split-group signature → verdict.
#               Full species retained (Fe I + Fe II usable for the RYA-309 run;
#               all species kept for the later RYA-308 all-element F-star row).
#
#               Mr. Code does NOT download from VALD — Ryan submits manually and
#               drops the raw .gz; this script verifies/processes on delivery.
#               The superseded HFS-OFF list is quarantined (git mv), never
#               deleted (DATA_SAFETY.md).
#
#               Wavelength convention: VALD delivers VACUUM below 2000 Å despite
#               the WL_air label. Recorded per line; NO vacuum->air conversion
#               here (done at the spectrum-matching boundary; cf. RYA-269/303).
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-14
# Last modified: 2026-06-14
# Linear issue: RYA-309 — RUN: Procyon Fe (I+II) dual-method vs literature
#
# -----------------------------------------------------------------------------
# PARAMS (Heiter+2015 GBS fundamental Teff/logg; extraction vmic from delivery)
#   Procyon HD 61421 : Teff 6530, logg 3.96, [Fe/H] -0.04, vmic 1.7 (extraction)
# -----------------------------------------------------------------------------
# DEPENDENCIES
#   External: pandas
#   Shared:   data/linelists/vald_parse.py (read_vald_header, parse_vald_long)
# =============================================================================

import gzip
import shutil
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from data.linelists.vald_parse import read_vald_header, parse_vald_long  # noqa: E402
from pipeline.build_linelist import compute_vald_proximity_flag           # noqa: E402

LINELISTS = REPO_ROOT / 'data' / 'linelists'
VALD_PROCYON = Path.home() / 'Documents' / 'Exoplanet Codex' / 'VALD' / 'Procyon'

VACUUM_BELOW_A = 2000.0
# RYA-270 Fe working-pool curation: priority=1 for Fe I/II whose central_depth
# falls in the working window (re-applied to the HFS-ON list so it is a drop-in
# for the old linelist consumed by lines_fit / abundances_derive).
CD_MIN, CD_MAX = 0.05, 0.80

# Procyon HFS-ON delivery (FTP, 2026-06-12), three disjoint bands 1150-11000 Å.
BANDS = ['uv', 'optical', 'nir']
BAND_RANGE = {
    'uv':      (1150.0, 3000.0),
    'optical': (3000.0, 6910.0),    # priority — HARPS coverage 3780-6910
    'nir':     (6910.0, 11000.0),
}
RANGE_TOL = 5.0

# HFS-capable species for the split-group signature (primary HFS-ON test).
HFS_SPECIES = {
    ('Mn', 'I'), ('Mn', 'II'), ('Co', 'I'), ('Co', 'II'), ('Cu', 'I'),
    ('V', 'I'), ('V', 'II'), ('Sc', 'I'), ('Sc', 'II'),
    ('Ba', 'II'), ('Eu', 'II'), ('Na', 'I'), ('Al', 'I'),
}

STAR = dict(name='Procyon', hd='HD 61421', teff=6530, logg=3.96, feh=-0.04,
            vmic=1.7)

OUT_PATH = LINELISTS / 'linelist_procyon.csv'
OLD_PATH = LINELISTS / 'linelist_procyon.csv'   # quarantined before overwrite
QUARANTINE = LINELISTS / 'linelist_procyon_rya270_hfsoff_quarantine.csv'


def read_gz_header(gz_path):
    with gzip.open(gz_path, 'rt', errors='replace') as f:
        line1 = f.readline().strip()
        truncated = line1.startswith('WARNING')
        meta = f.readline().strip() if truncated else line1
    fields = [p.strip() for p in meta.split(',')]
    return {'wl_start': float(fields[0]), 'wl_end': float(fields[1]),
            'n_selected': int(fields[2]), 'truncated': truncated}


def discover():
    """Map each .gz/plain VALD delivery in the Procyon folder to a band by its
    header wavelength range (robust to whatever request IDs VALD assigns)."""
    chosen = {}
    for gz in sorted(VALD_PROCYON.glob('*.gz')):
        try:
            h = read_gz_header(gz)
        except Exception:
            continue
        band = next((b for b, (lo, hi) in BAND_RANGE.items()
                     if abs(h['wl_start'] - lo) <= RANGE_TOL
                     and abs(h['wl_end'] - hi) <= RANGE_TOL), None)
        if band is None:
            continue
        score = (not h['truncated'], h['n_selected'], gz.stat().st_mtime)
        if band not in chosen or score > chosen[band][1]:
            chosen[band] = (gz, score)
    return {b: v[0] for b, v in chosen.items()}


def split_group_signature(records, window=0.2, min_cluster=3):
    hfs = [(r['element'], r['ion'], r['wavelength'])
           for r in records if (r['element'], r['ion']) in HFS_SPECIES]
    by_species = defaultdict(list)
    for elem, ion, wl in hfs:
        by_species[(elem, ion)].append(wl)
    in_groups = 0
    detail = {}
    for key, wls in sorted(by_species.items()):
        wls = sorted(wls)
        grp = 0
        i = 0
        while i < len(wls):
            j = i
            while j < len(wls) and wls[j] - wls[i] <= window:
                j += 1
            if j - i >= min_cluster:
                grp += (j - i)
            i = j
        if grp:
            detail[f'{key[0]} {key[1]}'] = grp
        in_groups += grp
    total = len(hfs)
    return total, in_groups, (in_groups / total if total else 0.0), detail


def intake_one(band, gz_path):
    raw_txt = LINELISTS / f'vald_procyon_{band}_hfson_raw.txt'
    req_lo, req_hi = BAND_RANGE[band]
    row = dict(band=band, file=gz_path.name,
               requested=f'{req_lo:.0f}-{req_hi:.0f}', truncated='', parsed='',
               coverage='', hfs_pct='', verdict='')
    with gzip.open(gz_path, 'rb') as fi, open(raw_txt, 'wb') as fo:
        shutil.copyfileobj(fi, fo)
    hdr = read_vald_header(raw_txt)
    if hdr['truncated']:
        row.update(truncated='YES', verdict='REJECT (truncated)')
        return None, row
    row['truncated'] = 'no'
    records, report = parse_vald_long(raw_txt)
    n, sel = report['n_parsed'], hdr['n_selected']
    row['parsed'] = f'{n}/{sel}' + ('' if n == sel else ' INCOMPLETE')
    wls = [r['wavelength'] for r in records]
    row['coverage'] = f'{min(wls):.1f}-{max(wls):.1f}' if wls else '—'
    total, in_grp, frac, detail = split_group_signature(records)
    row['hfs_pct'] = f'{frac*100:.0f}% ({in_grp}/{total})'
    row['_hfs_detail'] = detail
    if total == 0:
        row['verdict'] = 'ACCEPT (HFS inconclusive — no capable lines)'
    elif frac >= 0.40:
        row['verdict'] = 'ACCEPT'
    elif frac <= 0.05:
        row['verdict'] = 'REJECT (HFS-OFF)'
        return None, row
    else:
        row['verdict'] = f'AMBIGUOUS ({frac*100:.0f}%) — escalate'
        return None, row
    for r in records:
        r['_band'] = band
        r['_src'] = f'{band}_{gz_path.stem.split(".")[-1]}'
    return records, row


def build(all_records, registry):
    merged = pd.DataFrame(all_records)
    before = len(merged)
    merged = merged.drop_duplicates(
        subset=['element', 'ion', 'wavelength', 'log_gf'], keep='first')
    seam = before - len(merged)

    # RYA-270 curation re-applied to the HFS-ON list (drop-in for the old file):
    #   vald_proximity_flag — neighbour-density contamination on the FULL pool.
    #   priority=1          — Fe I/II working pool (central_depth ∈ [CD_MIN,CD_MAX]).
    prox = compute_vald_proximity_flag(pd.DataFrame({
        'wavelength_air_A': merged['wavelength'].to_numpy(),
        'central_depth': merged['central_depth'].to_numpy()}))
    is_fe = (merged['element'] == 'Fe') & merged['ion'].isin(['I', 'II'])
    in_win = merged['central_depth'].between(CD_MIN, CD_MAX)
    priority = np.where(is_fe & in_win, 1, 0)

    out = pd.DataFrame({
        'element': merged['element'],
        'ion': merged['ion'],
        'wavelength_air_A': merged['wavelength'],
        'excitation_potential_eV': merged['e_low_eV'],
        'log_gf': merged['log_gf'],
        'loggf_source': 'VALD3',
        'nist_grade': '',
        'damping_rad': merged['damping_rad'],
        'damping_stark': merged['damping_stark'],
        'damping_vdW': merged['damping_vdW'],
        'central_depth': merged['central_depth'],
        'vald_proximity_flag': np.asarray(prox),
        'blend_flag': False,
        'priority': priority,
        'notes': '',
        'extraction_source': merged['_src'],
        'wavelength_convention': merged['wavelength'].map(
            lambda w: 'vacuum' if w < VACUUM_BELOW_A else 'air'),
    }).sort_values('wavelength_air_A')
    src_lines = '\n'.join(
        f'#   {b:<8}: vald_procyon_{b}_hfson_raw.txt ({registry[b].name}), '
        f'{BAND_RANGE[b][0]:.0f}-{BAND_RANGE[b][1]:.0f} A'
        for b in BANDS if b in registry)
    header = f"""\
# linelist_procyon.csv — {STAR['name']} ({STAR['hd']}) per-system VALD3 linelist (RYA-309)
# Generated: {date.today().isoformat()}  by scripts/intake_build_procyon.py
# Stellar params (extraction): Teff={STAR['teff']} K, logg={STAR['logg']}, [Fe/H]={STAR['feh']:+.2f}, vmic={STAR['vmic']} km/s
#   (Heiter+2015 GBS fundamental Teff/logg; HFS-ON re-extraction)
# HFS splitting: ON (VERIFIED by split-group signature — supersedes the RYA-270
#   HFS-unsplit linelist_procyon.csv, quarantined as
#   linelist_procyon_rya270_hfsoff_quarantine.csv)
# Full species retained (Fe I + Fe II for the RYA-309 Fe run; all species kept
#   for the later RYA-308 all-element row).
# Extractions merged (all HFS-ON):
{src_lines}
# Wavelength convention: 'vacuum' for lambda < 2000 A (VALD delivers vacuum
#   below 2000 A despite the WL_air label), 'air' for lambda >= 2000 A.
#   NO vacuum->air conversion at build — done at the spectrum-matching boundary.
# Seam duplicates dropped at band endpoints: {seam}
# nist_grade empty pending NIST cross-validation; blend_flag=False everywhere.
"""
    with open(OUT_PATH, 'w') as f:
        f.write(header)
        out.to_csv(f, index=False)
    return out, seam


def main():
    print(f"{'='*72}\nRYA-309 Step 3.1 — Procyon HFS-ON VALD intake + linelist rebuild\n{'='*72}")
    registry = discover()
    rows, accepted = [], []
    for band in BANDS:
        gz = registry.get(band)
        if gz is None:
            rows.append(dict(band=band, file='—',
                             requested=f'{BAND_RANGE[band][0]:.0f}-{BAND_RANGE[band][1]:.0f}',
                             truncated='', parsed='', coverage='', hfs_pct='',
                             verdict='NOT SUBMITTED'))
            continue
        recs, row = intake_one(band, gz)
        rows.append(row)
        if recs is not None and row['verdict'].startswith('ACCEPT'):
            accepted.extend(recs)

    print("\nINTAKE GATE")
    cols = ['band', 'file', 'requested', 'truncated', 'parsed', 'coverage',
            'hfs_pct', 'verdict']
    print(pd.DataFrame(rows)[cols].to_string(index=False))

    print("\nHFS split-group detail (top species):")
    for r in rows:
        d = r.get('_hfs_detail')
        if d:
            top = ', '.join(f'{s}:{n}' for s, n in
                            sorted(d.items(), key=lambda x: -x[1])[:6])
            print(f"  {r['band']:<8} {r['hfs_pct']:<16} {top}")

    n_sub = sum(1 for r in rows if r['verdict'] != 'NOT SUBMITTED')
    n_acc = sum(r['verdict'].startswith('ACCEPT') for r in rows)
    print("\nBUILD")
    if n_acc == n_sub and n_acc > 0:
        # quarantine the superseded HFS-OFF list before overwrite (git-tracked)
        if OLD_PATH.exists() and not QUARANTINE.exists():
            import subprocess
            r = subprocess.run(['git', 'mv', str(OLD_PATH), str(QUARANTINE)],
                               cwd=REPO_ROOT, capture_output=True, text=True)
            if r.returncode != 0:           # not tracked or already moved
                shutil.move(str(OLD_PATH), str(QUARANTINE))
            print(f"  quarantined old HFS-OFF list -> {QUARANTINE.name}")
        out, seam = build(accepted, registry)
        fe1 = len(out[(out.element == 'Fe') & (out.ion == 'I')])
        fe2 = len(out[(out.element == 'Fe') & (out.ion == 'II')])
        nsp = out.groupby(['element', 'ion']).ngroups
        vac = (out.wavelength_convention == 'vacuum').sum()
        print(f"  {OUT_PATH.relative_to(REPO_ROOT)}: {len(out)} rows, "
              f"{out.wavelength_air_A.min():.1f}-{out.wavelength_air_A.max():.1f} A, "
              f"{nsp} species")
        print(f"    Fe I={fe1}  Fe II={fe2}  | vacuum<2000={vac} air>=2000={len(out)-vac}"
              f"  | seam dups dropped={seam}")
    else:
        print(f"  HELD — {n_acc}/{n_sub} bands ACCEPT; not rebuilding.")
    print(f"\n{'='*72}\nFiles: {n_acc} ACCEPT, "
          f"{sum(r['verdict'].startswith('REJECT') for r in rows)} REJECT.\n{'='*72}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
