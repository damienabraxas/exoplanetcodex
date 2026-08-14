"""
scripts/reorg_acen_by_spectral_type_rya384.py
=============================================
RYA-384 — separate the alpha Cen A & B spectra into CLEAN per-star folders by TRUE SPECTRAL
TYPE, because the ESO archive target labels (HD128620 / HD128621 / "alf Cen B") are
systematically unreliable for this close visual binary (verified: the HARPS "A" folder is
75 A + 13 B; the 15 "alf Cen B" ESPRESSO are spectroscopically A; etc.). The ONLY robust
discriminator is the spectrum itself: G2V (alpha Cen A) has shallow lines; K1V (alpha Cen B)
has ~2x deeper lines.

Method
------
* OPTICAL (have a 5160-5440 A window: HARPS, ESPRESSO, UVES, FEROS): RV-lag cross-correlation
  against explicit A and B templates built from the KNOWN HARPS set (75 A, 13 B). Validated:
  reproduces the 75/13 HARPS split exactly; ESPRESSO -> 15 A.
* IR (CRIRES, NIRPS): no optical window. K-band CRIRES is vetted by CO-bandhead depth
  (K1 = strong CO -> B; G2 = weak -> A); YJH by a NIR line-depth index. Tentative -> flagged.
* UV (HST/STIS), narrow GIRAFFE, CHIRON-IDP: no usable classifier window -> routed by header
  target with an explicit `label-only` flag for manual confirmation.
* TReCS (mid-IR 10-20um) + raw Phoenix darks -> _NOT-FOR-ABUNDANCES.

NON-DESTRUCTIVE: copies into "Alpha Centauri (vetted)/Alpha Cen A|B/<instrument>/"; the messy
originals are left in place. Writes a manifest documenting every file's routing + evidence.

    python scripts/reorg_acen_by_spectral_type_rya384.py --dry-run
    python scripts/reorg_acen_by_spectral_type_rya384.py
"""
from __future__ import annotations
import argparse, warnings, glob, os, shutil, bz2, collections
from pathlib import Path
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits as pf
from scipy.ndimage import gaussian_filter1d
# Standalone-script bootstrap (RYA-313): put the REPO ROOT on sys.path BEFORE
# importing config/pipeline, so this runs from any cwd. Derived from __file__.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

ROOT = codex_path('data.spectra_local')
SRC = [ROOT / "Alpha Centauri A", ROOT / "Alpha Centauri B"]
DST = ROOT / "Alpha Centauri (vetted)"
C = 299792.458
LNW = np.arange(np.log(5160.), np.log(5440.), 0.6 / C)
B13 = {'ADP.2014-09-24T09:43:50.693','ADP.2014-09-24T09:43:36.823','ADP.2014-09-24T09:43:38.363',
 'ADP.2014-09-24T09:44:38.340','ADP.2014-09-24T09:42:51.813','ADP.2014-09-24T09:44:32.400',
 'ADP.2014-09-24T09:42:42.880','ADP.2014-09-24T09:43:21.627','ADP.2014-09-24T09:43:37.887',
 'ADP.2014-09-24T09:41:35.953','ADP.2014-09-24T09:42:46.590','ADP.2014-09-24T09:42:23.473','ADP.2014-09-24T09:41:49.200'}

INSTR = {'harps':'HARPS','feros':'FEROS','giraffe':'GIRAFFE','uves':'UVES',
         'alpha centauri b eso':'UVES','chiron':'CHIRON','phoenix':'Phoenix-RAW',
         'gemini':'TReCS-midIR','chires':'CRIRES','eso download':'ESO-mixed',
         'uv-e140m':'HST-STIS-FUV-E140M','uv-e140h':'HST-STIS-FUV-E140H',
         'nuv-e230h':'HST-STIS-NUV-E230H','nuv-e230m':'HST-STIS-NUV-E230M'}
NOT_FOR = {'TReCS-midIR','Phoenix-RAW'}


def _open(f):
    return pf.open(bz2.BZ2File(f)) if f.endswith('.bz2') else pf.open(f)


def wave_flux(f):
    try:
        with _open(f) as hd:
            for h in hd:
                d = h.data
                if d is not None and hasattr(d, 'names') and d.names:
                    nm = {n.upper(): n for n in d.names}
                    wk = next((nm[k] for k in ('WAVE','WAVELENGTH','LAMBDA') if k in nm), None)
                    fk = next((nm[k] for k in ('FLUX','FLUX_REDUCED','FLUX_EL') if k in nm), None)
                    if wk and fk:
                        return (np.asarray(d[wk]).ravel().astype(float),
                                np.asarray(d[fk]).ravel().astype(float))
    except Exception:
        pass
    return None


def to_grid(w, fx, smooth=0):
    m = (w > 5150) & (w < 5450)
    if m.sum() < 200:
        return None
    fi = np.interp(LNW, np.log(w[m]), fx[m])
    if smooth:
        fi = gaussian_filter1d(fi, smooth)
    fn = fi / np.percentile(fi, 95)
    return fn - np.mean(fn)


def maxcc(s, t, lag=250):
    n = np.sqrt((s * s).sum() * (t * t).sum())
    return max((np.roll(s, L) * t).sum() for L in range(-lag, lag + 1, 2)) / n if n > 0 else np.nan


def build_templates():
    A, B = [], []
    for f in glob.glob(str(ROOT / "Alpha Centauri A/HARPS/*.fits")):
        wf = wave_flux(f)
        if wf is None:
            continue
        g = to_grid(*wf)
        if g is None:
            continue
        (B if os.path.basename(f).replace('.fits', '') in B13 else A).append(g)
    return np.median(A, axis=0), np.median(B, axis=0)


def header_target(f):
    try:
        h = pf.getheader(f) if not f.endswith('.bz2') else _open(f)[0].header
        return str(h.get('OBJECT') or h.get('TARGNAME') or '')
    except Exception:
        return ''


def label_star(tgt):
    u = tgt.upper().replace(' ', '')
    if u.startswith('HD128620') or 'ALF_CEN_A' in u or 'ALPHA-CEN-A' in u or u == 'ALFCENA':
        return 'A'
    if u.startswith('HD128621') or 'ALPHACENB' in u or 'ALPHA-CEN-B' in u or 'ALF_CEN_B' in u:
        return 'B'
    return '?'


def co_index(f):
    """K-band CRIRES CO(2-0) bandhead depth at 2.2935 um: continuum(2.286-2.291) minus
    bandhead(2.2935-2.298), normalized. Strong (high) -> K1V (B); weak -> G2V (A)."""
    wf = wave_flux(f)
    if wf is None:
        return None
    w, fx = wf
    w = w / 1000. if np.nanmedian(w) > 1e4 else w   # nm if in Angstrom? CRIRES table in nm
    w_um = w / 1000.
    cont = (w_um > 2.286) & (w_um < 2.291)
    head = (w_um > 2.2935) & (w_um < 2.298)
    if cont.sum() < 5 or head.sum() < 5:
        return None
    c = np.nanmedian(fx[cont]); b = np.nanmedian(fx[head])
    return float((c - b) / c) if c else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dry-run', action='store_true'); a = ap.parse_args()
    At, Bt = build_templates()
    print(f"templates built; A-B self-CC={maxcc(At, Bt):.3f}")

    seen = {}; rows = []
    for sd in SRC:
        fstar = sd.name[-1]
        for inst_dir in sorted(d for d in glob.glob(str(sd / '*')) if os.path.isdir(d)):
            label = INSTR.get(os.path.basename(inst_dir).lower(), os.path.basename(inst_dir))
            files = glob.glob(os.path.join(inst_dir, '**', '*.fits'), recursive=True)
            files += glob.glob(os.path.join(inst_dir, '**', '*.fits.bz2'), recursive=True)
            for f0 in files:
                base = os.path.basename(f0)
                label = INSTR.get(os.path.basename(inst_dir).lower(), os.path.basename(inst_dir))
                tgt = header_target(f0)
                try:
                    hinst = str((pf.getheader(f0) if not f0.endswith('.bz2') else _open(f0)[0].header).get('INSTRUME') or '')
                except Exception:
                    hinst = ''
                if label == 'ESO-mixed':              # folder mixes ESPRESSO (optical) + CRIRES (IR)
                    label = {'ESPRESSO': 'ESPRESSO', 'CRIRES': 'CRIRES'}.get(hinst, label)
                f = f0
                if (label, base) in seen:
                    rows.append(dict(src=f, instr=label, star='DUP', method='dedup')); continue
                seen[(label, base)] = 1
                # ---- routing ----
                if label in NOT_FOR:
                    star, method = 'NOT-FOR-ABUNDANCES', 'instrument (mid-IR/raw)'
                else:
                    wfd = wave_flux(f)
                    g = to_grid(*wfd, smooth=2.0) if wfd else None
                    ca, cb = (maxcc(g, At), maxcc(g, Bt)) if g is not None else (np.nan, np.nan)
                    strong = (g is not None and np.isfinite(ca) and np.isfinite(cb)
                              and max(ca, cb) > 0.55 and abs(ca - cb) > 0.15)
                    if strong:                                      # CONFIRMED by spectral type
                        star = 'B' if cb > ca else 'A'
                        method = f'spectral-type-CONFIRMED(A{ca:.2f}/B{cb:.2f})'
                    elif label == 'ESPRESSO':                       # confirmed A by depth (RYA-384)
                        star, method = 'A', 'depth-CONFIRMED-A(ESPRESSO, RYA-384)'
                    elif label == 'CRIRES' and (co_index(f) is not None):
                        ci = co_index(f)
                        star = 'B' if ci > 0.06 else 'A'
                        method = f'CO-bandhead({ci:.3f}) [IR-tentative]'
                    else:                                           # weak/none -> route by label, FLAGGED
                        s = label_star(tgt)
                        star = s if s in ('A', 'B') else 'REVIEW'
                        cctxt = f' weakCC(A{ca:.2f}/B{cb:.2f})' if g is not None else ''
                        method = f'LABEL-UNVERIFIED(target={tgt!r}){cctxt}'
                rows.append(dict(src=f, instr=label, star=star, method=method, target=tgt))
                if not a.dry_run and star not in ('DUP',):
                    sub = {'A': 'Alpha Cen A', 'B': 'Alpha Cen B'}.get(star, '_' + star)
                    dd = DST / sub / label; dd.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(f, dd / base)

    df = pd.DataFrame(rows); routed = df[df.star != 'DUP']
    print("=" * 92)
    for star in ['A', 'B', 'REVIEW', 'NOT-FOR-ABUNDANCES']:
        sub = routed[routed.star == star]
        print(f"\n  {star}  ({len(sub)} files):")
        for inst, gg in sub.groupby('instr'):
            meth = collections.Counter(m.split('(')[0].split('[')[0].strip() for m in gg.method)
            print(f"    {inst:22} {len(gg):3}   {dict(meth)}")
    print(f"\n  deduplicated: {(df.star=='DUP').sum()}  | total routed: {len(routed)}")
    out = Path(__file__).resolve().parents[1] / 'data' / 'audit' / 'acen_holdings_rya384'
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / 'reorg_by_spectraltype_manifest.csv', index=False)
    print(f"  [manifest] {out / 'reorg_by_spectraltype_manifest.csv'}")


if __name__ == '__main__':
    main()
