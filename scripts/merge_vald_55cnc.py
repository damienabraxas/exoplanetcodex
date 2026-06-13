# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         merge_vald_55cnc.py
# Module:       scripts (linelist build tools)
# Description:  RYA-269 — unpack/verify/merge the four 55 Cnc A VALD3
#               HFS-ON extractions (optical + UV-a + UV-b + NIR) into the
#               per-system linelist data/linelists/linelist_55cnc.csv.
#               NIST cross-validation is OUT of scope (RYA-64 follow-up).
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-12
# Last modified: 2026-06-13
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
# Data:     data/linelists/vald_55cnc_raw.txt       (optical, HFS-ON, 3780–6910 Å)
#           data/linelists/vald_55cnc_uv_a_raw.txt  (UV-a, HFS-ON, 1150–2000 Å, 019516)
#           data/linelists/vald_55cnc_uv_b_raw.txt  (UV-b, HFS-ON, 2000–3780 Å, 019517)
#           data/linelists/vald_55cnc_nir_raw.txt   (NIR, HFS-ON, 6910–17000 Å, 019518)
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

# Extraction registry. All four are HFS-ON (Codex convention, RYA-269 DECISION).
# central_depth_threshold is the VALD detection threshold of the REQUEST.
# Dedup keeps the copy from the LOWER threshold (more inclusive) extraction.
# In practice the four extractions do not overlap (the wavelength ranges are
# contiguous with small seam gaps), so n_dups = 0 is expected.
EXTRACTIONS = [
    # (source_tag, path, request_id, central_depth_threshold)
    ('uv_a',          LINELISTS / 'vald_55cnc_uv_a_raw.txt', '019516', 0.05),
    ('uv_b',          LINELISTS / 'vald_55cnc_uv_b_raw.txt', '019517', 0.05),
    ('optical_019385', LINELISTS / 'vald_55cnc_raw.txt',     'see header note', 0.01),
    ('nir',            LINELISTS / 'vald_55cnc_nir_raw.txt', '019518', 0.05),
]
CD_THRESHOLD = {tag: cd for tag, _, _, cd in EXTRACTIONS}

# RYA-197 molecular species — VALD is non-authoritative for these (RYA-64).
MOLECULAR_RYA197 = ['CH', 'CN', 'C2', 'NH', 'OH', 'CO']
MOLECULAR_NOTE = ('MOLECULAR — VALD non-authoritative, '
                  'pending RYA-197 ExoMol/Brooke replacement')

# Molecular hydrides — not RYA-197 CNO carriers; excluded from atomic EW only.
MOLECULAR_HYDRIDES = ['SiH', 'MgH']
HYDRIDE_NOTE = ('MOLECULAR HYDRIDE — excluded from atomic EW; '
                'NOT RYA-197 CNO scope')

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
        print(f"  {tag:<17}: {hdr['wl_start']:.0f}–{hdr['wl_end']:.0f} Å requested, "
              f"{hdr['n_selected']} transitions, truncation: {trunc}")
        if hdr['truncated']:
            critical(f"{tag}: VALD output truncated at 100000 lines — "
                     f"re-extract at higher central_depth before merging")
    if CRITICALS:
        return finish()

    # ── Parse all extractions ────────────────────────────────────────────────
    frames = []
    fail_total, fail_examples, n_lines_total = 0, [], 0
    tag_list = [tag for tag, *_ in EXTRACTIONS]
    for tag, path, req_id, cd in EXTRACTIONS:
        records, report = parse_vald_long(path)
        fail_total += report['n_failures']
        fail_examples += [(tag,) + ex for ex in report['examples']]
        n_lines_total += report['n_parsed'] + report['n_failures']
        delivered_min = min(r['wavelength'] for r in records)
        delivered_max = max(r['wavelength'] for r in records)
        hdr = headers[tag]
        complete = report['n_parsed'] == hdr['n_selected']
        print(f"  {tag:<17}: delivered {delivered_min:.3f}–{delivered_max:.3f} Å, "
              f"parsed {report['n_parsed']}/{hdr['n_selected']} selected "
              f"({'complete' if complete else 'INCOMPLETE'})")
        if not complete:
            critical(f"{tag}: parsed count {report['n_parsed']} != "
                     f"lines selected {hdr['n_selected']}")
        df = pd.DataFrame(records)
        df['extraction_source'] = tag
        frames.append(df)

    # UV-a must reach the 2000 Å boundary; UV-b must reach 3780 Å (RYA-269 spec)
    uv_a = frames[tag_list.index('uv_a')]
    uv_b = frames[tag_list.index('uv_b')]
    uv_a_max = uv_a['wavelength'].max()
    uv_b_max = uv_b['wavelength'].max()
    if uv_a_max < 1999.0:
        critical(f"UV-a extraction stops at {uv_a_max:.3f} Å — does not reach 2000 Å")
    else:
        print(f"  UV-a reach check: last line at {uv_a_max:.3f} Å vs 2000 Å boundary — OK")
    if uv_b_max < 3779.0:
        critical(f"UV-b extraction stops at {uv_b_max:.3f} Å — does not reach 3780 Å")
    else:
        print(f"  UV-b reach check: last line at {uv_b_max:.3f} Å vs 3780 Å boundary — OK")

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
        gap_str = 'none — contiguous coverage'
    print(f"[2/4] Coverage gaps: {gap_str}")
    for a, b in gaps:
        if b - a > 2.0:
            critical(f"coverage gap {a:.3f}–{b:.3f} Å exceeds 2 Å")

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

    # HFS-convention gate: all extractions are HFS-ON; log_gf conflicts == 0
    # is the operative confirmation (the 204-conflict failure from the HFS-OFF
    # run must not recur). Already ensured by the block above.
    n_hfs_conflicts = len(conflicts)
    if n_hfs_conflicts > 0:
        critical(f"HFS-convention conflict gate FAIL: {n_hfs_conflicts} conflicts — "
                 f"check that all extractions used HFS splitting ON")

    n_dups = len(dropped)
    n_intra = int(kept.duplicated('dedup_key').sum())
    merged = kept
    print(f"[3/4] Merge: {len(merged)} unique lines "
          f"({n_dups} duplicates removed, {len(conflicts)} log_gf conflicts)")
    print(f"       HFS-convention conflicts: {n_hfs_conflicts} ✓" if n_hfs_conflicts == 0
          else f"       HFS-convention conflicts: {n_hfs_conflicts} *** FAIL")
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

    # Molecular hydrides (SiH, MgH) — excluded from atomic EW, distinct from
    # the RYA-197 six (not CNO carriers, not ExoMol/Brooke scope). Reuses the
    # notes column, the same mechanism as the RYA-197 exclusion above.
    hydride_mask = merged['element'].isin(MOLECULAR_HYDRIDES)
    merged.loc[hydride_mask, 'notes'] = HYDRIDE_NOTE
    hyd_counts = merged.loc[hydride_mask, 'element'].value_counts()
    hyd_str = ', '.join(f"{m}: {hyd_counts.get(m, 0)}" for m in MOLECULAR_HYDRIDES)
    print(f"       Molecular hydrides flagged (excl. atomic EW, NOT RYA-197): "
          f"{int(hydride_mask.sum())} ({hyd_str})")

    # Residual: any non-atomic, non-RYA-197, non-hydride species — report only.
    _ATOMS = set("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V
        Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag
        Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
        Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Th U""".split())
    other_mol = merged.loc[
        ~merged['element'].isin(_ATOMS) & ~mol_mask & ~hydride_mask, 'element']
    if len(other_mol):
        oc = other_mol.value_counts()
        print(f"       Other unclassified molecular species: "
              + ', '.join(f"{k}: {v}" for k, v in oc.items()))
    print(f"       Parse failures: {fail_total}")

    # CRITICAL science-line presence gates (RYA-269 smoke spec)
    def present(elem, ion, lo, hi):
        m = ((merged['element'] == elem) & (merged['ion'] == ion) &
             merged['wavelength'].between(lo, hi))
        return merged.loc[m, 'wavelength'].tolist()

    # O I 7771/7774/7775 Å triplet — from NIR extraction (6910–17000 Å).
    # Gate requires the two STRONG members (7771.944, cd≈0.068; 7774.166, cd≈0.056)
    # which are both well above the cd=0.05 threshold. The weakest member
    # (7775.388, cd≈0.041) is genuinely below the 0.05 threshold at 55 Cnc params
    # and is correctly absent from any cd=0.05 extraction — not an extraction error.
    o_triplet = present('O', 'I', 7770.5, 7776.5)
    o_strong = [w for w in o_triplet if abs(w - 7771.944) < 0.5 or abs(w - 7774.166) < 0.5]
    if len(o_strong) < 2:
        critical(f"O I 7771/7774 Å (strong triplet members) absent — found "
                 f"{o_triplet} — NIR extraction is wrong")
    else:
        weak_note = ('; 7775.388 below cd=0.05 threshold (cd≈0.041) — expected absent'
                     if not any(abs(w - 7775.388) < 0.5 for w in o_triplet) else '')
        print(f"       O I triplet check: {[f'{w:.3f}' for w in o_triplet]} ✓{weak_note}")

    # Li I 6707 Å — from optical extraction (3780–6910 Å)
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
# HFS splitting: ON (Codex permanent convention — all four extractions, RYA-269 DECISION)
# Extractions merged (all HFS-ON):
#   uv_a          : vald_55cnc_uv_a_raw.txt (RyanSchmitt.019516), 1150-2000 A,
#       central_depth 0.05, vacuum wavelengths (VALD delivers vacuum below 2000 A)
#   uv_b          : vald_55cnc_uv_b_raw.txt (RyanSchmitt.019517), 2000-3780 A,
#       central_depth 0.05
#   optical_019385: vald_55cnc_raw.txt, 3780-6910 A, central_depth 0.01, complete
#       (NOTE: request ID not preserved in file; NOT the truncated web file
#        RyanSchmitt.019385 which was truncated at 5603 A and is quarantined)
#   nir           : vald_55cnc_nir_raw.txt (RyanSchmitt.019518), 6910-17000 A,
#       central_depth 0.05
# Wavelength convention: 'vacuum' for lambda < 2000 A (VALD delivers vacuum
#   below 2000 A despite the WL_air label), 'air' for lambda >= 2000 A.
#   NO vacuum->air conversion performed at merge stage — conversion to match
#   STIS vacuum frames happens at the spectrum-matching boundary.
# Superseded HFS-OFF files quarantined (not deleted) per RYA-269 spec.
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
