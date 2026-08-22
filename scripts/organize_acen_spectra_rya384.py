"""
scripts/organize_acen_spectra_rya384.py
=======================================
RYA-384 follow-up — build a CLEAN, vetted alpha Cen A/B folder tree for the Sirius transfer.

The current folders are NOT separated by star: UVES + the entire HST/STIS UV set are
DUPLICATED into both the "Alpha Centauri A" and "Alpha Centauri B" folders and internally MIX
A (HD128620) and B (HD128621). This routes every FITS by HEADER TRUTH (not folder) into a clean
tree, de-duplicating the A<->B folder copies. NON-DESTRUCTIVE: copies into a new vetted root;
the messy originals are left in place until you verify and delete them.

    python scripts/organize_acen_spectra_rya384.py --dry-run   # plan only
    python scripts/organize_acen_spectra_rya384.py             # copy

HD128620 = alpha Cen A (G2V); HD128621 = alpha Cen B (K1V).
"""
from __future__ import annotations
import argparse, warnings, glob, os, shutil, bz2, collections
from pathlib import Path
warnings.filterwarnings('ignore')
import pandas as pd
from astropy.io import fits as pf
# Standalone-script bootstrap (RYA-313): put the REPO ROOT on sys.path BEFORE
# importing config/pipeline, so this runs from any cwd. Derived from __file__.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

ROOT = codex_path('data.spectra_local')
SRC = [ROOT / "Alpha Centauri A", ROOT / "Alpha Centauri B"]
DST = ROOT / "Alpha Centauri (vetted)"

# instrument folder-name -> clean instrument label
INSTR = {
    'harps': 'HARPS', 'feros': 'FEROS', 'giraffe': 'GIRAFFE', 'uves': 'UVES',
    'alpha centauri b eso': 'UVES', 'chiron': 'CHIRON', 'phoenix': 'Phoenix-RAW',
    'uv-e140m': 'HST-STIS-FUV-E140M', 'uv-e140h': 'HST-STIS-FUV-E140H',
    'nuv-e230h': 'HST-STIS-NUV-E230H', 'nuv-e230m': 'HST-STIS-NUV-E230M',
}


from pipeline.star_id import resolve_star  # noqa: E402

def instr_label(folder_name: str) -> str:
    return INSTR.get(folder_name.lower(), folder_name)


def classify(h, folder_star: str, instr: str) -> tuple:
    """(star 'A'|'B'|'REVIEW', reason). Header truth wins; indeterminate -> provenance."""
    tgt = str(h.get('OBJECT') or h.get('TARGNAME') or '').upper().replace(' ', '')
    # RYA-964: the alias set lives in system_catalog.csv, not in this function. The inline
    # tests this replaces had to be kept in sync by hand in two files, and one of them
    # already needed a guard against `HD128620`'s test swallowing `HD128621`.
    raw = str(h.get('OBJECT') or h.get('TARGNAME') or '')
    sid = resolve_star(raw)
    if sid == 'alpha_cen_a':
        return 'A', f'target={tgt} (resolve_star -> alpha_cen_a)'
    if sid == 'alpha_cen_b':
        return 'B', f'target={tgt} (resolve_star -> alpha_cen_b)'
    if 'DARK' in tgt:
        return 'A', f'raw/dark ({tgt}) -> Phoenix-RAW'
    if 'GL-559' in tgt or 'GL559' in tgt:
        return 'REVIEW', f'Gl-559 system name (FEROS A/B disambig owed)'
    if instr == 'CHIRON':
        return 'B', 'CHIRON CTIO 2021, target in ext -> B by provenance (confirm owed)'
    # coord-named or empty: route by folder origin, flag
    return ('REVIEW', f'indeterminate target={tgt!r}, folder={folder_star} (verify)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    seen = {}            # filename -> (star, instr) to dedup A<->B duplicates
    rows, copied = [], 0
    for sd in SRC:
        folder_star = sd.name[-1]
        for inst in sorted(d for d in glob.glob(str(sd / '*')) if os.path.isdir(d)):
            label = instr_label(os.path.basename(inst))
            files = glob.glob(os.path.join(inst, '**', '*.fits'), recursive=True)
            files += glob.glob(os.path.join(inst, '**', '*.fits.bz2'), recursive=True)
            for f in files:
                base = os.path.basename(f)
                key = (label, base)
                if key in seen:                       # already routed from the other folder
                    rows.append(dict(src=f, base=base, instr=label, star='DUP', reason='dedup'))
                    continue
                try:
                    opener = bz2.BZ2File(f) if f.endswith('.bz2') else f
                    with pf.open(opener) as hd:
                        star, reason = classify(hd[0].header, folder_star, label)
                except Exception as e:
                    star, reason = 'REVIEW', f'open-err {type(e).__name__}'
                seen[key] = (star, label)
                sub = {'A': 'Alpha Cen A', 'B': 'Alpha Cen B', 'REVIEW': '_NEEDS-REVIEW'}[star]
                dest_dir = DST / sub / label
                rows.append(dict(src=f, base=base, instr=label, star=star, reason=reason,
                                 dest=str(dest_dir / base)))
                if not a.dry_run:
                    dest_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(f, dest_dir / base)
                    copied += 1

    df = pd.DataFrame(rows)
    # summary
    routed = df[df.star != 'DUP']
    print("=" * 92)
    print(f"  RYA-384 vetted alpha Cen tree {'(DRY RUN)' if a.dry_run else '-> ' + str(DST)}")
    print("=" * 92)
    for star in ['A', 'B', 'REVIEW']:
        sub = routed[routed.star == star]
        tag = {'A': 'Alpha Cen A', 'B': 'Alpha Cen B', 'REVIEW': '_NEEDS-REVIEW'}[star]
        print(f"\n  {tag}  ({len(sub)} files):")
        for instr, g in sub.groupby('instr'):
            print(f"    {instr:22} {len(g):3}   e.g. {g.iloc[0]['reason'][:48]}")
    print(f"\n  de-duplicated A<->B folder copies: {(df.star == 'DUP').sum()}")
    print(f"  total routed (unique): {len(routed)}   copied: {copied}")
    # manifest
    out = Path(__file__).resolve().parents[1] / 'data' / 'audit' / 'acen_holdings_rya384'
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / 'reorg_manifest.csv', index=False)
    print(f"  [manifest] {out / 'reorg_manifest.csv'}")


if __name__ == '__main__':
    main()
