# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         merge_vald_55cnc.py
# Module:       scripts (linelist build tools)
# Description:  RYA-269 — unpack/verify/merge the three 55 Cnc A VALD3
#               extractions (optical + UV + NIR) into the per-system
#               linelist data/linelists/linelist_55cnc.csv.
#               NIST cross-validation is OUT of scope (RYA-64 follow-up).
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-12
# Last modified: 2026-06-12
# Linear issue: RYA-269 — BUILD: 55 Cnc A VALD UV + NIR unpack/verify/merge
#
# -----------------------------------------------------------------------------
# KEY REFERENCES
# -----------------------------------------------------------------------------
# Ryabchikova et al. 2015 — Phys. Scr. 90, 054005 — VALD3 description
# RYA-64 how-to comment — extraction settings and truncation post-mortem
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
# External: pandas
# Data:     data/linelists/vald_55cnc_raw.txt      (optical, complete 3780–6910 Å)
#           data/linelists/vald_55cnc_uv_raw.txt   (UV, request 019509)
#           data/linelists/vald_55cnc_nir_raw.txt  (NIR, 5000–30000 Å @ cd 0.05)
# =============================================================================

import sys
from datetime import date
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from config.constants import PATHS
from data.linelists.vald_parse import read_vald_header, parse_vald_long

LINELISTS = PATHS['linelists']
OUT_PATH = LINELISTS / 'linelist_55cnc.csv'

# Extraction registry. central_depth_threshold is the VALD detection threshold
# of the REQUEST (RYA-64 comments): optical + UV at 0.01, NIR re-extraction at
# 0.05 per the truncation post-mortem. Dedup keeps the copy from the LOWER
# threshold (more inclusive) extraction.
EXTRACTIONS = [
    # (source_tag, path, request_id, central_depth_threshold)
    ('optical_019385', LINELISTS / 'vald_55cnc_raw.txt',     'see header note', 0.01),
    ('uv',             LINELISTS / 'vald_55cnc_uv_raw.txt',  '019509',          0.01),
    ('nir',            LINELISTS / 'vald_55cnc_nir_raw.txt', 'not preserved',   0.05),
]
CD_THRESHOLD = {tag: cd for tag, _, _, cd in EXTRACTIONS}

# RYA-197 molecular species — VALD is non-authoritative for these (RYA-64).
MOLECULAR_RYA197 = ['CH', 'CN', 'C2', 'NH', 'OH', 'CO']
MOLECULAR_NOTE = ('MOLECULAR — VALD non-authoritative, '
                  'pending RYA-197 ExoMol/Brooke replacement')

# Vacuum/air boundary: VALD delivers vacuum wavelengths below 2000 Å despite
# the WL_air(A) column label. NO conversion at merge stage — conversion to
# match STIS vacuum frames happens at the spectrum-matching boundary (RYA-269).
VACUUM_BELOW_A = 2000.0

PARSE_FAIL_CRITICAL_FRAC = 0.001   # > 0.1 % unparseable lines = CRITICAL
LOGGF_CONFLICT_TOL = 0.001         # duplicate-record log_gf disagreement gate

CRITICALS = []


def critical(msg):
    CRITICALS.append(msg)
    print(f"\n  *** CRITICAL: {msg}")


def main():
    print("RYA-269 — 55 Cnc A VALD merge → linelist_55cnc.csv")

    # ── [1/4] Inventory ───────────────────────────────────────────────────────
    print(f"\n[1/4] Inventory: {len(EXTRACTIONS)} extractions found")
    headers = {}
    for tag, path, req_id, cd in EXTRACTIONS:
        if not path.exists():
            critical(f"{tag}: file not found — {path}")
            continue
        hdr = read_vald_header(path)
        headers[tag] = hdr
        trunc = 'YES' if hdr['truncated'] else 'NO'
        print(f"  {tag:<15}: {hdr['wl_start']:.0f}–{hdr['wl_end']:.0f} Å requested, "
              f"{hdr['n_selected']} transitions, truncation: {trunc}")
        if hdr['truncated']:
            critical(f"{tag}: VALD output truncated at 100000 lines — "
                     f"re-extract at higher central_depth before merging")
    if CRITICALS:
        return finish()

    # ── Parse all extractions ────────────────────────────────────────────────
    frames = []
    fail_total, fail_examples, n_lines_total = 0, [], 0
    for tag, path, req_id, cd in EXTRACTIONS:
        records, report = parse_vald_long(path)
        fail_total += report['n_failures']
        fail_examples += [(tag,) + ex for ex in report['examples']]
        n_lines_total += report['n_parsed'] + report['n_failures']
        delivered_min = min(r['wavelength'] for r in records)
        delivered_max = max(r['wavelength'] for r in records)
        hdr = headers[tag]
        complete = report['n_parsed'] == hdr['n_selected']
        print(f"  {tag:<15}: delivered {delivered_min:.3f}–{delivered_max:.3f} Å, "
              f"parsed {report['n_parsed']}/{hdr['n_selected']} selected "
              f"({'complete' if complete else 'INCOMPLETE'})")
        if not complete:
            critical(f"{tag}: parsed count {report['n_parsed']} != "
                     f"lines selected {hdr['n_selected']}")
        df = pd.DataFrame(records)
        df['extraction_source'] = tag
        frames.append(df)

    # UV must actually reach the 3780 Å request boundary (RYA-269 Step 2)
    uv = frames[[t for t, *_ in EXTRACTIONS].index('uv')]
    uv_max = uv['wavelength'].max()
    if uv_max < 3779.0:
        critical(f"UV extraction stops at {uv_max:.3f} Å — does not reach 3780 Å")
    else:
        print(f"  UV reach check: last UV line at {uv_max:.3f} Å vs 3780 Å request "
              f"boundary — OK")

    # ── [2/4] Coverage gaps ──────────────────────────────────────────────────
    spans = sorted((f['wavelength'].min(), f['wavelength'].max(),
                    f['extraction_source'].iloc[0]) for f in frames)
    gaps = []
    cur_lo, cur_hi, _ = spans[0]
    for lo, hi, tag in spans[1:]:
        if lo > cur_hi:
            gaps.append((cur_hi, lo))
        cur_hi = max(cur_hi, hi)
    if gaps:
        gap_str = '; '.join(f"{a:.3f}–{b:.3f} Å ({b-a:.3f} Å)" for a, b in gaps)
    else:
        gap_str = 'none — contiguous 1150–30000 Å coverage'
    print(f"[2/4] Coverage gaps: {gap_str}")
    for a, b in gaps:
        if b - a > 1.0:
            critical(f"coverage gap {a:.3f}–{b:.3f} Å exceeds 1 Å")

    if fail_total:
        print(f"  Parse failures: {fail_total} (first {min(5, len(fail_examples))}):")
        for ex in fail_examples[:5]:
            print(f"    {ex}")
    if n_lines_total and fail_total / n_lines_total > PARSE_FAIL_CRITICAL_FRAC:
        critical(f"parse failure rate {fail_total/n_lines_total:.2%} > 0.1%")
    if CRITICALS:
        return finish()

    # ── [3/4] Merge + dedup ──────────────────────────────────────────────────
    merged = pd.concat(frames, ignore_index=True)
    merged['dedup_key'] = list(zip(
        merged['element'], merged['ion'],
        merged['wavelength'].round(3), merged['e_low_eV'].round(3),
    ))

    # Overlap dedup is CROSS-EXTRACTION only ("extractions overlap" — RYA-269
    # spec). Records sharing the key WITHIN one extraction are hyperfine/
    # isotopic components of the same multiplet — physically distinct
    # transitions, never duplicates — and are all retained. For keys present
    # in 2+ extractions, keep ALL rows from the LOWER central_depth-threshold
    # (more inclusive) extraction and drop the other sources' copies; any
    # dropped copy whose log_gf has no kept counterpart within tolerance is a
    # genuine VALD inconsistency → CRITICAL.
    merged['cd_threshold'] = merged['extraction_source'].map(CD_THRESHOLD)
    key_min_cd = merged.groupby('dedup_key')['cd_threshold'].transform('min')
    keep_mask = merged['cd_threshold'] == key_min_cd
    kept = merged[keep_mask].copy()
    dropped = merged[~keep_mask]

    conflicts = []
    kept_loggf = kept.groupby('dedup_key')['log_gf'].agg(list)
    for key, grp in dropped.groupby('dedup_key', sort=False):
        kept_vals = kept_loggf[key]
        for src, lg in zip(grp['extraction_source'], grp['log_gf']):
            if min(abs(lg - kv) for kv in kept_vals) > LOGGF_CONFLICT_TOL:
                conflicts.append((key, src, lg, kept_vals))
    if conflicts:
        critical(f"{len(conflicts)} cross-extraction duplicate record(s) "
                 f"disagree on log_gf by > {LOGGF_CONFLICT_TOL} — "
                 f"refusing to pick silently")
        for key, src, lg, kept_vals in conflicts[:5]:
            print(f"    {key}: dropped {src} log_gf={lg} vs kept {kept_vals}")
        return finish()

    n_dups = len(dropped)
    n_intra = int(kept.duplicated('dedup_key').sum())
    merged = kept
    print(f"[3/4] Merge: {len(merged)} unique lines "
          f"({n_dups} duplicates removed, {len(conflicts)} log_gf conflicts)")
    print(f"       Intra-extraction key-sharing components retained "
          f"(HFS/isotopic, not duplicates): {n_intra}")

    # Molecular flagging — RYA-197 species only (per spec)
    merged['notes'] = ''
    mol_mask = merged['element'].isin(MOLECULAR_RYA197)
    merged.loc[mol_mask, 'notes'] = MOLECULAR_NOTE
    mol_counts = merged.loc[mol_mask, 'element'].value_counts()
    counts_str = ', '.join(f"{m}: {mol_counts.get(m, 0)}" for m in MOLECULAR_RYA197)
    print(f"       Molecular lines flagged for RYA-197: {int(mol_mask.sum())} "
          f"({counts_str})")

    # Other molecular species delivered by VALD but outside the RYA-197 six —
    # reported (not silently dropped, not RYA-197-flagged; Ryan to decide).
    _ATOMS = set("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V
        Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag
        Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
        Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Th U""".split())
    other_mol = merged.loc[~merged['element'].isin(_ATOMS) & ~mol_mask, 'element']
    if len(other_mol):
        oc = other_mol.value_counts()
        print(f"       Other molecular species (NOT RYA-197-flagged): "
              + ', '.join(f"{k}: {v}" for k, v in oc.items()))
    print(f"       Parse failures: {fail_total}")

    # CRITICAL science-line presence gates (RYA-269 smoke spec)
    def present(elem, ion, lo, hi):
        m = ((merged['element'] == elem) & (merged['ion'] == ion) &
             merged['wavelength'].between(lo, hi))
        return merged.loc[m, 'wavelength'].tolist()

    o_triplet = present('O', 'I', 7770.5, 7776.5)
    if len(o_triplet) < 3:
        critical(f"O I 7771/7774/7775 Å triplet absent/incomplete — found "
                 f"{o_triplet} — NIR extraction is wrong")
    else:
        print(f"       O I triplet check: {[f'{w:.3f}' for w in o_triplet]} ✓")
    li = present('Li', 'I', 6707.0, 6708.5)
    if not li:
        critical("Li I 6707 Å absent from merged output")
    else:
        li_src = merged.loc[(merged['element'] == 'Li') &
                            merged['wavelength'].between(6707.0, 6708.5),
                            'extraction_source'].unique().tolist()
        print(f"       Li I 6707 check: {[f'{w:.3f}' for w in li]} "
              f"(source: {li_src}) ✓")
    if CRITICALS:
        return finish()

    # ── [4/4] Write output ───────────────────────────────────────────────────
    out = pd.DataFrame({
        'element'                 : merged['element'],
        'ion'                     : merged['ion'],
        'wavelength_air_A'        : merged['wavelength'],
        'excitation_potential_eV' : merged['e_low_eV'],
        'log_gf'                  : merged['log_gf'],
        'loggf_source'            : 'VALD3',
        'nist_grade'              : '',          # populated in NIST follow-up phase
        'damping_rad'             : merged['damping_rad'],
        'damping_stark'           : merged['damping_stark'],
        'damping_vdW'             : merged['damping_vdW'],
        'central_depth'           : merged['central_depth'],
        'blend_flag'              : False,       # vetted exclusions only (RYA-209)
        'priority'                : 0,           # placeholder
        'notes'                   : merged['notes'],
        'extraction_source'       : merged['extraction_source'],
        'wavelength_convention'   : merged['wavelength'].map(
            lambda w: 'vacuum' if w < VACUUM_BELOW_A else 'air'),
    }).sort_values('wavelength_air_A')

    header = f"""\
# linelist_55cnc.csv — 55 Cnc A per-system VALD3 linelist (RYA-269)
# Generated: {date.today().isoformat()}  by scripts/merge_vald_55cnc.py
# Stellar params: Teff=5196 K, logg=4.41, [Fe/H]=+0.35, vmic=0.9 km/s
# Extractions merged:
#   optical_019385: vald_55cnc_raw.txt, 3780-6910 A, central_depth 0.01, complete
#       (NOTE: this is the complete 125615-line optical delivery; the web file
#        RyanSchmitt.019385 was truncated at 5603 A and is NOT merged here)
#   uv            : vald_55cnc_uv_raw.txt (RyanSchmitt.019509), 1150-3780 A,
#       central_depth 0.01
#   nir           : vald_55cnc_nir_raw.txt, 5000-30000 A, central_depth 0.05
# Wavelength convention: 'vacuum' for lambda < 2000 A (VALD delivers vacuum
#   below 2000 A despite the WL_air label), 'air' for lambda >= 2000 A.
#   NO vacuum->air conversion performed at merge stage — conversion to match
#   STIS vacuum frames happens at the spectrum-matching boundary.
# nist_grade empty pending NIST cross-validation (RYA-64 follow-up phase).
# blend_flag=False everywhere — vetted exclusions only (RYA-209 architecture).
"""
    with open(OUT_PATH, 'w') as f:
        f.write(header)
        out.to_csv(f, index=False)
    rel = OUT_PATH.relative_to(REPO_ROOT)
    print(f"[4/4] Wrote {rel} ({len(out)} rows)")
    return finish()


def finish():
    if CRITICALS:
        print(f"\nRESULT: CRITICAL — {len(CRITICALS)} stop condition(s); "
              f"no output written" if not OUT_PATH.exists()
              else f"\nRESULT: CRITICAL — {len(CRITICALS)} stop condition(s)")
        for c in CRITICALS:
            print(f"  - {c}")
        return 1
    print("\nRESULT: OK")
    return 0


if __name__ == '__main__':
    sys.exit(main())
