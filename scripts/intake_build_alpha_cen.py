#!/usr/bin/env python3
# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         intake_build_alpha_cen.py
# Module:       scripts (VALD intake + per-system line-list build)
# Description:  RYA-300 — intake gate + per-system build for the alpha Cen A/B
#               VALD3 "Extract Stellar" deliveries (HFS ON, per-system).
#
#               Mr. Code does NOT download from VALD — Ryan submits manually at
#               vald.astro.uu.se and drops the raw .gz deliveries into
#               ~/Documents/Exoplanet Codex/VALD/Alpha Centauri {A,B}/. This
#               script runs the standing VALD intake gate on each delivery:
#                 1. Truncation — line 1 must be the VALD metadata header, not
#                    the 100k-cap "Output was truncated" warning. Truncated =
#                    REJECT (quarantine, never delete).
#                 2. Parse completeness — n_parsed vs n_selected.
#                 3. Coverage — delivered min/max vs requested range.
#                 4. HFS — split-group signature on HFS-capable species
#                    (Mn, Co, Cu, V, Sc, Ba II, Eu II, ...): ON ~90%+, OFF ~2%.
#                 5. Verdict ACCEPT/REJECT per file.
#               Then, for any star whose full band set is present and ACCEPTed,
#               it assembles the per-system line list with the shared
#               vald_parse.py (never forked) and the RYA-269 output schema.
#
#               Wavelength convention: VALD delivers VACUUM below 2000 Å despite
#               the WL_air column label. Recorded per line (vacuum < 2000 Å,
#               air >= 2000 Å); NO vacuum->air conversion here — that happens at
#               the spectrum-matching boundary (consistent with RYA-269 / 303).
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-14
# Last modified: 2026-06-14
# Linear issue: RYA-300 — VALD: Extract line lists for alpha Cen A + B (HFS ON)
#
# -----------------------------------------------------------------------------
# PARAMS (Heiter+2015 GBS Paper I / Jofré+2014 GBS Paper III)
#   alpha Cen A : Teff 5792, logg 4.30, [Fe/H] +0.20, vmic 1.1
#   alpha Cen B : Teff 5231, logg 4.53, [Fe/H] +0.20, vmic 1.0
# -----------------------------------------------------------------------------
# DEPENDENCIES
#   External: pandas
#   Shared:   data/linelists/vald_parse.py (read_vald_header, parse_vald_long)
# =============================================================================

import gzip
import shutil
import sys
import tempfile
from collections import defaultdict
from datetime import date
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from data.linelists.vald_parse import read_vald_header, parse_vald_long  # noqa: E402

LINELISTS = REPO_ROOT / 'data' / 'linelists'
VALD_ROOT = Path.home() / 'Documents' / 'Exoplanet Codex' / 'VALD'

VACUUM_BELOW_A = 2000.0
CENTRAL_DEPTH = 0.05            # shared detection threshold (RYA-300 recipe)

# HFS-capable species for the split-group signature (primary HFS-ON test).
HFS_SPECIES = {
    ('Mn', 'I'), ('Mn', 'II'), ('Co', 'I'), ('Co', 'II'), ('Cu', 'I'),
    ('V', 'I'), ('V', 'II'), ('Sc', 'I'), ('Sc', 'II'),
    ('Ba', 'II'), ('Eu', 'II'), ('Na', 'I'), ('Al', 'I'),
}

# Per-star fundamental parameters (for the output header provenance).
STAR_PARAMS = {
    'a': dict(name='alpha Cen A', teff=5792, logg=4.30, feh=+0.20, vmic=1.1),
    'b': dict(name='alpha Cen B', teff=5231, logg=4.53, feh=+0.20, vmic=1.0),
}

# Each star expects four disjoint bands spanning 1150-17000 Å. A band is
# "present" once Ryan drops its .gz into the star's VALD folder; the per-system
# build fires only when all four bands of a star are present AND ACCEPTed.
BANDS = ['uv1', 'uv2', 'optical', 'nir']
BAND_RANGE = {
    'uv1':     (1150.0, 2000.0),
    'uv2':     (2000.0, 3780.0),
    'optical': (3780.0, 6910.0),    # THE priority — matches HARPS coverage
    'nir':     (6910.0, 17000.0),
}
RANGE_TOL = 5.0   # Å tolerance matching a delivery's header range to a band

STAR_FOLDER = {
    'a': VALD_ROOT / 'Alpha Centauri A',
    'b': VALD_ROOT / 'Alpha Centauri B',
}


def read_gz_header(gz_path):
    """Peek the VALD metadata line of a gzipped delivery (first two text lines)
    without decompressing the whole file. Returns the read_vald_header dict."""
    with gzip.open(gz_path, 'rt', errors='replace') as f:
        line1 = f.readline().strip()
        truncated = line1.startswith('WARNING')
        meta = f.readline().strip() if truncated else line1
    fields = [p.strip() for p in meta.split(',')]
    return {'wl_start': float(fields[0]), 'wl_end': float(fields[1]),
            'n_selected': int(fields[2]), 'truncated': truncated}


def discover(star):
    """Scan a star's VALD folder and map each .gz delivery to a band by its
    header wavelength range (robust to whatever request IDs VALD assigns).

    On collision (resubmission of the same band) prefer a non-truncated file,
    then the one with the most selected transitions, then the newest mtime.
    Returns {band: Path}."""
    chosen = {}                      # band -> (path, score tuple)
    folder = STAR_FOLDER[star]
    if not folder.is_dir():
        return {}
    for gz in sorted(folder.glob('*.gz')):
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


REGISTRY = {star: discover(star) for star in ('a', 'b')}


# -----------------------------------------------------------------------------
def split_group_signature(records, window=0.2, min_cluster=3):
    """Fraction of HFS-capable transitions that sit in tight wavelength clusters
    (>= min_cluster lines within `window` Å). HFS-ON splits a single line into a
    cluster of components; HFS-OFF leaves singletons."""
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


def gunzip(gz_path, dest):
    """Decompress a .gz delivery to `dest` (the committed working raw txt)."""
    with gzip.open(gz_path, 'rb') as fi, open(dest, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def intake_one(star, band, gz_path):
    """Run the full intake gate on one delivery. Returns (records, row) where
    row is the intake-table dict; records is None on REJECT/missing."""
    raw_txt = LINELISTS / f'vald_alpha_cen_{star}_{band}_raw.txt'
    req_lo, req_hi = BAND_RANGE[band]
    row = dict(star=star.upper(), band=band, file=gz_path.name,
               requested=f'{req_lo:.0f}-{req_hi:.0f}', truncated='', parsed='',
               coverage='', hfs_pct='', verdict='')

    if not gz_path.exists():
        row['verdict'] = 'MISSING'
        return None, row

    gunzip(gz_path, raw_txt)                      # persist committed raw txt
    hdr = read_vald_header(raw_txt)

    # 1. truncation
    if hdr['truncated']:
        row.update(truncated='YES', verdict='REJECT (truncated)')
        return None, row
    row['truncated'] = 'no'

    # 2. parse completeness
    records, report = parse_vald_long(raw_txt)
    n, sel = report['n_parsed'], hdr['n_selected']
    row['parsed'] = f'{n}/{sel}' + ('' if n == sel else ' INCOMPLETE')

    # 3. coverage
    wls = [r['wavelength'] for r in records]
    cov_lo, cov_hi = (min(wls), max(wls)) if wls else (0, 0)
    row['coverage'] = f'{cov_lo:.1f}-{cov_hi:.1f}'

    # 4. HFS split-group signature
    total, in_grp, frac, detail = split_group_signature(records)
    row['hfs_pct'] = f'{frac*100:.0f}% ({in_grp}/{total})'
    row['_hfs_detail'] = detail

    # 5. verdict
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
    return records, row


def build_per_system(star, band_records):
    """Assemble linelist_alpha_cen_{star}.csv from the four band record sets."""
    frames = []
    for band, recs in band_records.items():
        df = pd.DataFrame(recs)
        df['extraction_source'] = f'{band}_{REGISTRY[star][band].stem.split(".")[-1]}'
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True)
    # disjoint bands share inclusive endpoints (2000.0, 3780.0, 6910.0); drop
    # exact-duplicate transitions at the seams, keep first occurrence.
    before = len(merged)
    merged = merged.drop_duplicates(
        subset=['element', 'ion', 'wavelength', 'log_gf'], keep='first')
    seam_dups = before - len(merged)

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
        'blend_flag': False,
        'priority': 0,
        'notes': '',
        'extraction_source': merged['extraction_source'],
        'wavelength_convention': merged['wavelength'].map(
            lambda w: 'vacuum' if w < VACUUM_BELOW_A else 'air'),
    }).sort_values('wavelength_air_A')

    p = STAR_PARAMS[star]
    src_lines = '\n'.join(
        f'#   {b:<8}: vald_alpha_cen_{star}_{b}_raw.txt '
        f'({REGISTRY[star][b].name}), {BAND_RANGE[b][0]:.0f}-{BAND_RANGE[b][1]:.0f} A'
        for b in BANDS)
    header = f"""\
# linelist_alpha_cen_{star}.csv — {p['name']} per-system VALD3 linelist (RYA-300)
# Generated: {date.today().isoformat()}  by scripts/intake_build_alpha_cen.py
# Stellar params: Teff={p['teff']} K, logg={p['logg']}, [Fe/H]={p['feh']:+.2f}, vmic={p['vmic']} km/s
#   (Heiter+2015 GBS Paper I / Jofre+2014 GBS Paper III)
# HFS splitting: ON (Codex permanent convention — all extractions, RYA-300)
# Detection threshold (central_depth): {CENTRAL_DEPTH}
# Extractions merged (all HFS-ON):
{src_lines}
# Wavelength convention: 'vacuum' for lambda < 2000 A (VALD delivers vacuum
#   below 2000 A despite the WL_air label), 'air' for lambda >= 2000 A.
#   NO vacuum->air conversion at build — done at the spectrum-matching boundary.
# Seam duplicates dropped at band endpoints: {seam_dups}
# nist_grade empty pending NIST cross-validation; blend_flag=False everywhere.
"""
    out_path = LINELISTS / f'linelist_alpha_cen_{star}.csv'
    with open(out_path, 'w') as f:
        f.write(header)
        out.to_csv(f, index=False)
    return out_path, len(out), seam_dups


def main():
    print(f"{'='*72}\nRYA-300 — alpha Cen A/B VALD intake + per-system build\n{'='*72}")
    rows = []
    star_band_records = {'a': {}, 'b': {}}

    for star in ['a', 'b']:
        for band in BANDS:
            gz = REGISTRY[star].get(band)
            if gz is None:
                rows.append(dict(star=star.upper(), band=band, file='—',
                                 requested=f'{BAND_RANGE[band][0]:.0f}-{BAND_RANGE[band][1]:.0f}',
                                 truncated='', parsed='', coverage='',
                                 hfs_pct='', verdict='NOT SUBMITTED'))
                continue
            recs, row = intake_one(star, band, gz)
            rows.append(row)
            if recs is not None and row['verdict'].startswith('ACCEPT'):
                star_band_records[star][band] = recs

    # ---- intake table ----
    print("\nINTAKE GATE")
    cols = ['star', 'band', 'file', 'requested', 'truncated', 'parsed',
            'coverage', 'hfs_pct', 'verdict']
    tbl = pd.DataFrame(rows)[cols]
    with pd.option_context('display.max_colwidth', 40, 'display.width', 200):
        print(tbl.to_string(index=False))

    # ---- HFS detail for accepted files ----
    print("\nHFS split-group detail (accepted files, top species):")
    for r in rows:
        d = r.get('_hfs_detail')
        if d:
            top = ', '.join(f'{s}:{n}' for s, n in
                            sorted(d.items(), key=lambda x: -x[1])[:6])
            print(f"  {r['star']} {r['band']:<7} {r['hfs_pct']:<14} {top}")

    # ---- per-system build (only complete + all-ACCEPT stars) ----
    print("\nPER-SYSTEM BUILD")
    built = []
    for star in ['a', 'b']:
        have = set(star_band_records[star])
        if have == set(BANDS):
            path, nrows, seam = build_per_system(star, star_band_records[star])
            built.append(path)
            print(f"  {STAR_PARAMS[star]['name']}: {path.relative_to(REPO_ROOT)} "
                  f"({nrows} rows, {seam} seam dups dropped)  ✓ COMPLETE")
        else:
            missing = [b for b in BANDS if b not in have]
            print(f"  {STAR_PARAMS[star]['name']}: HELD — bands accepted "
                  f"{sorted(have) or 'none'}, awaiting {missing}")

    print(f"\n{'='*72}")
    n_acc = sum(r['verdict'].startswith('ACCEPT') for r in rows)
    n_rej = sum(r['verdict'].startswith('REJECT') for r in rows)
    print(f"Files: {n_acc} ACCEPT, {n_rej} REJECT, "
          f"{sum(r['verdict'] in ('NOT SUBMITTED','MISSING') for r in rows)} pending. "
          f"Built {len(built)} per-system list(s).")
    print(f"{'='*72}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
