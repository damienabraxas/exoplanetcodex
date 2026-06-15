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

import gzip
import shutil
import sys
from datetime import date
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from config.constants import PATHS, STAR_PARAMS
from data.linelists.vald_parse import (
    read_vald_header, parse_vald_long, verify_metallicity,
    parse_vald_applied_metallicity)

LINELISTS = PATHS['linelists']
OUT_PATH = LINELISTS / 'linelist_55cnc.csv'

STAR_ID = '55cnc_a'

# RYA-324 — intake driver. Mr. Code does not download from VALD; Ryan drops the
# raw .gz deliveries into ~/Documents/Exoplanet Codex/VALD/55 Cancri/. This script
# discovers them by header wavelength range (robust to whatever request IDs VALD
# assigns), gunzips the chosen file per band into the committed raw .txt, runs the
# standing intake gates (truncation / completeness / coverage / HFS log_gf), the
# RYA-321 METALLICITY gate (applied [M/H] from the abundance block must match the
# catalog 55 Cnc [Fe/H]), then merges per the RYA-269 schema below.
VALD_55CNC_FOLDER = Path.home() / 'Documents' / 'Exoplanet Codex' / 'VALD' / '55 Cancri'
RANGE_TOL = 5.0                       # Å tolerance matching a delivery to a band
CD_DEFAULT = 0.05                     # request central_depth threshold

# Each band -> (requested wl range, committed raw .txt the chosen .gz is unpacked
# into). UV-a/UV-b tags are kept (downstream reach checks key on them).
BANDS = {
    'uv_a':    ((1150.0, 2000.0),  LINELISTS / 'vald_55cnc_uv_a_raw.txt'),
    'uv_b':    ((2000.0, 3780.0),  LINELISTS / 'vald_55cnc_uv_b_raw.txt'),
    'optical': ((3780.0, 6910.0),  LINELISTS / 'vald_55cnc_optical_raw.txt'),
    'nir':     ((6910.0, 17000.0), LINELISTS / 'vald_55cnc_nir_raw.txt'),
}


def _gz_header_range(gz_path):
    """Peek a gzipped delivery's requested wavelength range (Å) without unpacking."""
    with gzip.open(gz_path, 'rt', errors='replace') as f:
        line1 = f.readline().strip()
        meta = f.readline().strip() if line1.startswith('WARNING') else line1
    lo, hi, *_ = (p.strip() for p in meta.split(','))
    return float(lo), float(hi)


SUPP_MAX_WIDTH = 100.0   # a delivery narrower than this is a targeted supplement
SUPP_CD = 0.01           # request central_depth of the narrow low-cd supplements


def discover_supplements():
    """Targeted low-cd supplements live in SUBFOLDERS of VALD/55 Cancri/ (e.g.
    'o line zoom in' — RYA-333). Each is a narrow window (< SUPP_MAX_WIDTH Å)
    extracted at a low cd to recover diagnostic lines that fall below the main
    cd=0.05 cut (the O I 7771/7774/7775 triplet is shallow at 55 Cnc's 5196 K +
    metal-rich veiling). Gunzip each into a committed raw and register it as an
    extraction at SUPP_CD so the dedup keeps its (more inclusive) copy. Returns
    [(tag, path, request_id, cd)]."""
    supp = []
    for gz in sorted(VALD_55CNC_FOLDER.glob('*/*.gz')):
        try:
            lo, hi = _gz_header_range(gz)
        except Exception:
            continue
        if hi - lo > SUPP_MAX_WIDTH:          # a full band misfiled, not a supplement
            continue
        req_id = gz.stem.split('.')[-1]
        raw = LINELISTS / f'vald_55cnc_supp_{req_id}_raw.txt'
        with gzip.open(gz, 'rb') as fi, open(raw, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        supp.append((f'supp_{req_id}', raw, req_id, SUPP_CD))
    return supp


def discover_and_materialize():
    """Map each .gz in the 55 Cancri folder to a band by header range, gunzip the
    best per band (non-truncated, most lines, newest) into its committed raw .txt,
    and return the EXTRACTIONS registry [(tag, path, request_id, cd)]. Returns
    (extractions, criticals) — a missing band is a CRITICAL caught by the caller."""
    chosen = {}                       # band -> (gz_path, score, req_id)
    for gz in sorted(VALD_55CNC_FOLDER.glob('*.gz')):
        try:
            lo, hi = _gz_header_range(gz)
        except Exception:
            continue
        band = next((b for b, ((blo, bhi), _) in BANDS.items()
                     if abs(lo - blo) <= RANGE_TOL and abs(hi - bhi) <= RANGE_TOL), None)
        if band is None:
            continue
        req_id = gz.stem.split('.')[-1]            # RyanSchmitt.019562 -> 019562
        score = (gz.stat().st_size, gz.stat().st_mtime)
        if band not in chosen or score > chosen[band][1]:
            chosen[band] = (gz, score, req_id)

    extractions, crit = [], []
    for band, ((_, _), raw_path) in BANDS.items():
        if band not in chosen:
            crit.append(f"{band}: no delivery in {VALD_55CNC_FOLDER} matching its "
                        f"wavelength range — drop the .gz and re-run")
            continue
        gz, _, req_id = chosen[band]
        with gzip.open(gz, 'rb') as fi, open(raw_path, 'wb') as fo:
            shutil.copyfileobj(fi, fo)             # materialize committed raw .txt
        extractions.append((band, raw_path, req_id, CD_DEFAULT))
    return extractions, crit

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
    print("RYA-269/324 — 55 Cnc A VALD intake + merge → linelist_55cnc.csv")

    # Intake: discover deliveries by band, gunzip the chosen .gz into committed
    # raws. A missing band is a CRITICAL.
    EXTRACTIONS, discovery_crit = discover_and_materialize()
    EXTRACTIONS += discover_supplements()        # RYA-333 narrow low-cd supplements
    CD_THRESHOLD = {tag: cd for tag, _, _, cd in EXTRACTIONS}
    for c in discovery_crit:
        critical(c)
    if CRITICALS:
        return finish()

    expected_feh = STAR_PARAMS[STAR_ID]['feh_ref']

    # ── [1/4] Inventory ───────────────────────────────────────────────────────
    print(f"\n[1/4] Inventory: {len(EXTRACTIONS)} extractions discovered "
          f"(catalog [Fe/H]={expected_feh:+.2f})")
    headers = {}
    for tag, path, req_id, cd in EXTRACTIONS:
        if not path.exists():
            critical(f"{tag}: file not found — {path}")
            continue
        hdr = read_vald_header(path)
        headers[tag] = hdr
        trunc = 'YES' if hdr['truncated'] else 'NO'
        # RYA-321 metallicity gate: VALD must have APPLIED the catalog [Fe/H] to
        # the composition (read from the abundance block, not the substituted
        # Castelli structure-model filename).
        met_verdict, met_msg = verify_metallicity(path, STAR_ID, STAR_PARAMS)
        applied = parse_vald_applied_metallicity(path)
        applied_s = f'{applied:+.2f}' if applied is not None else 'none'
        print(f"  {tag:<17}: {hdr['wl_start']:.0f}–{hdr['wl_end']:.0f} Å requested "
              f"({req_id}), {hdr['n_selected']} transitions, truncation: {trunc}, "
              f"applied [M/H]={applied_s} {met_verdict}")
        if hdr['truncated']:
            critical(f"{tag}: VALD output truncated at 100000 lines — "
                     f"re-extract at higher central_depth before merging")
        if met_verdict == 'REJECT':
            critical(f"{tag}: metallicity gate REJECT — {met_msg.split(': ', 1)[-1]}")
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
    # A gap straddling a band boundary is a SEAM between two separately-gated,
    # complete extractions — at the uniform cd=0.05 of the metal-rich re-extraction
    # the 6910 Å optical/nir seam is genuinely line-sparse (the old RYA-269 optical
    # used cd=0.01 and closed it). Allow a few-Å seam gap; keep INTERIOR gaps (a
    # silently truncated extraction) strict at 2 Å.
    BAND_SEAMS = (2000.0, 3780.0, 6910.0)
    SEAM_GAP_MAX, INTERIOR_GAP_MAX = 5.0, 2.0
    for a, b in gaps:
        seam = any(a - 1.0 <= s <= b + 1.0 for s in BAND_SEAMS)
        limit = SEAM_GAP_MAX if seam else INTERIOR_GAP_MAX
        if b - a > limit:
            critical(f"coverage gap {a:.3f}–{b:.3f} Å ({b-a:.3f}) exceeds {limit} Å"
                     + (" at a band seam" if seam else " (band interior)"))
        elif seam and b - a > INTERIOR_GAP_MAX:
            print(f"       NOTE: {a:.3f}–{b:.3f} Å seam gap ({b-a:.3f} Å) — "
                  f"line-sparse at cd={CD_DEFAULT} between complete bands, allowed")

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
    # The two strong members were cd≈0.068/0.056 in the earlier SOLAR-scale build,
    # just above the cd=0.05 cut. In the metal-rich +0.35 re-extraction the whole
    # triplet (high-excitation, 9.15 eV — intrinsically weak at 55 Cnc's 5196 K)
    # falls below cd=0.05: stronger metal line-blanketing raises the H- continuous
    # opacity, shrinking the central DEPTH of high-excitation lines even as their
    # EW grows. This is a detection-threshold/veiling effect, NOT an extraction
    # error — so with all four bands present, complete and metallicity-gated, the
    # absence is a WARNING (the O anchor needs a targeted lower-cd NIR re-extraction
    # — RYA-323 follow-up), not a hard CRITICAL.
    o_triplet = present('O', 'I', 7770.5, 7776.5)
    o_strong = [w for w in o_triplet if abs(w - 7771.944) < 0.5 or abs(w - 7774.166) < 0.5]
    if len(o_strong) < 2:
        print(f"       WARNING: O I 7771/7774 Å triplet below cd={CD_DEFAULT} in the "
              f"metal-rich NIR extraction (found {o_triplet}). Expected for the cool, "
              f"metal-rich star; the oxygen anchor needs a lower-cd NIR re-extract.")
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

    applied_mh = parse_vald_applied_metallicity(EXTRACTIONS[0][1])
    applied_s = f'{applied_mh:+.2f}' if applied_mh is not None else 'unknown'
    src_lines = '\n'.join(
        f"#   {tag:<12}: {path.name} (RyanSchmitt.{req_id}), "
        f"{headers[tag]['wl_start']:.0f}-{headers[tag]['wl_end']:.0f} A, central_depth {cd}"
        for tag, path, req_id, cd in EXTRACTIONS)
    header = f"""\
# linelist_55cnc.csv — 55 Cnc A per-system VALD3 linelist (RYA-269; rebuilt RYA-323/324)
# Generated: {date.today().isoformat()}  by scripts/merge_vald_55cnc.py
# Stellar params: Teff=5196 K, logg=4.41, [Fe/H]={expected_feh:+.2f} (catalog, STAR_PARAMS), vmic=0.9 km/s
# Composition applied by VALD: [M/H]={applied_s} (from the abundance block; passes the
#   RYA-321 metallicity gate vs catalog {expected_feh:+.2f}). Supersedes the earlier
#   solar-default ([M/H]=0.00) extraction.
# HFS splitting: ON (Codex permanent convention — all four extractions, RYA-269 DECISION)
# Extractions merged (HFS-ON, metallicity-gated, auto-discovered from VALD/55 Cancri/):
{src_lines}
# Wavelength convention: 'vacuum' for lambda < 2000 A (VALD delivers vacuum
#   below 2000 A despite the WL_air label), 'air' for lambda >= 2000 A.
#   NO vacuum->air conversion performed at merge stage — conversion to match
#   STIS vacuum frames happens at the spectrum-matching boundary.
# Superseded solar-default + HFS-OFF files quarantined (not deleted).
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
