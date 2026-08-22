"""
scripts/audit_acen_holdings_rya384.py
=====================================
RYA-384 — alpha Cen A & B spectral-holdings audit (Part B/C), codex-data-audit approach
(reuses the 370/377 logic). READ-ONLY. Walks the external spectra root, and for every
instrument folder reports: file count, the HEADER-TRUTH target split (HD128620=alpha Cen A,
HD128621=alpha Cen B), sampled wavelength range / SNR / date, and a binary-contamination
verdict (the alpha Cen folders are NOT cleanly separated by star -- UVES + HST STIS are
duplicated into both folders and internally MIX A and B; selection must be by header target,
never folder name).

Part A (VALD coverage) is scripts/audit_vald_inventory_rya376.py.

    python scripts/audit_acen_holdings_rya384.py
"""
from __future__ import annotations
import warnings, glob, os, collections, bz2
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

DATA_ROOT = codex_path('data.spectra_local')
OUT = Path(__file__).resolve().parents[1] / 'data' / 'audit' / 'acen_holdings_rya384'
# HD128620 = alpha Cen A (G2V); HD128621 = alpha Cen B (K1V)
A_IDS = ('HD128620', 'ALF_CEN_A', 'ALPHA-CEN-A', 'ALPHACENA')
B_IDS = ('HD128621', 'ALPHA-CEN-B', 'ALPHACENB', 'ALF_CEN_B')


from pipeline.star_id import resolve_star  # noqa: E402


def _target(h):
    o = (h.get('OBJECT') or h.get('TARGNAME') or '')
    return str(o)


def _classify(name: str) -> str:
    """RYA-964: one alias lookup, from `system_catalog.csv`.

    The substring tests this replaces needed an explicit guard because `HD128621` contains
    no `HD128620` but the reverse-ordered `x in u` scan made the two orderings matter. An
    exact match on a normalised label has no such failure mode.
    """
    sid = resolve_star(name)
    return {'alpha_cen_a': 'A', 'alpha_cen_b': 'B'}.get(sid, '?')


def _open_header(f):
    if f.endswith('.bz2'):
        return pf.open(bz2.BZ2File(f))
    return pf.open(f)


def scan_instrument(folder: Path) -> dict:
    files = glob.glob(os.path.join(folder, '**', '*.fits'), recursive=True)
    files += glob.glob(os.path.join(folder, '**', '*.fits.bz2'), recursive=True)
    split = collections.Counter()
    wmin = wmax = None
    snr_sample = date_sample = instr = None
    for i, f in enumerate(files):
        try:
            with _open_header(f) as hd:
                h = hd[0].header
                split[_classify(_target(h))] += 1
                if i == 0:
                    instr = h.get('INSTRUME')
                    date_sample = h.get('DATE-OBS') or h.get('MJD-OBS')
                    snr_sample = (h.get('SNR') or h.get('HIERARCH ESO DRS SPE EXT SN0')
                                  or h.get('HIERARCH ESO QC ORDER MIDDLE SNR'))
                    wmin = h.get('WAVELMIN') or h.get('HIERARCH ESO INS WLEN MIN')
                    wmax = h.get('WAVELMAX') or h.get('HIERARCH ESO INS WLEN MAX')
        except Exception:
            split['(err)'] += 1
    a, b, q = split.get('A', 0), split.get('B', 0), split.get('?', 0)
    if a and b:
        verdict = f'MIXED A/B ({a} A, {b} B) -- select by header target, NOT folder'
    elif a and not b:
        verdict = 'clean A'
    elif b and not a:
        verdict = 'clean B'
    else:
        verdict = 'target indeterminate from header (format-specific read needed)'
    return dict(instrument=os.path.basename(folder), n_files=len(files), instr_hdr=instr,
                n_A=a, n_B=b, n_unknown=q, wmin=wmin, wmax=wmax, snr_sample=snr_sample,
                date_sample=date_sample, verdict=verdict)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for star in ['Alpha Centauri A', 'Alpha Centauri B']:
        sd = DATA_ROOT / star
        for inst in sorted(d for d in glob.glob(str(sd / '*')) if os.path.isdir(d)):
            if not any(glob.glob(os.path.join(inst, '**', e), recursive=True)
                       for e in ('*.fits', '*.fits.bz2')):
                continue
            r = scan_instrument(Path(inst))
            r = {'folder_star': star[-1], **r}
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / 'spectral_inventory.csv', index=False)
    print("=" * 100)
    print("  RYA-384 — alpha Cen spectral holdings (header-truth target split per instrument)")
    print("=" * 100)
    for _, r in df.iterrows():
        print(f"  [{r['folder_star']}] {r['instrument']:22} {r['n_files']:3} files | "
              f"A={r['n_A']} B={r['n_B']} ?={r['n_unknown']} | {r['verdict']}")
    print(f"\n  [out] {OUT / 'spectral_inventory.csv'}")


if __name__ == '__main__':
    main()
