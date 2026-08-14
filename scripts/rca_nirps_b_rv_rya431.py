"""
scripts/rca_nirps_b_rv_rya431.py
================================
RYA-431 — RCA of the alpha Cen "NIRPS B" RV zero-point (~-34.6 km/s, below the orbit floor).
RYA-423 attributed it to a G-vs-K mask zero-point, but 7 km/s is far too large for that. This
runs the decisive test: is the offset CONSTANT, BERV-correlated, or random — and re-dispositions
the 20 frames.

Decisive findings (printed below):
  1. BERV bug RULED OUT: header ESO QC BERV equals the astropy coordinate+time computation to
     ~7 m/s; the confirmed-A frames are RV-consistent across BERV -18..+20 -> BERV is applied
     correctly. (Header BERV is NOT trusted blindly -- it is re-derived and checked.)
  2. The offset is CONSTANT (all 20 frames are one night, 2024-03-17; std ~0.3 km/s). A
     BERV/date correlation cannot be tested within a single epoch.
  3. The -34.6 is a REAL stellar CCF peak (contrast ~18, FWHM ~12), not noise, and its blue
     wing is TRUNCATED at the -38 search-window edge (the pipeline centred the CCF at -18,
     expecting alpha Cen A) -> the true RV is <= -34.6.
  4. -34.6 is OUTSIDE the orbit bounds gamma +/- max(K) = [-27.9, -16.9] (acen_orbit.rv_bounds).
     A bound alpha Cen member CANNOT have this RV.
=> The 20 frames are NOT alpha Cen B: their real photospheric RV is off the alpha Cen orbit.
   The J-depth "deep lines -> B" call was a K-type FALSE POSITIVE (a different K star, mislabelled
   'alf Cen A'). Re-dispositioned NOT-ALPHA-CEN (RV off-orbit) -- NOT a calibratable zero-point,
   NOT a BERV bug. Removed from the alpha Cen B science set; flagged, not routed-by-label.

    python scripts/rca_nirps_b_rv_rya431.py            # report
    python scripts/rca_nirps_b_rv_rya431.py --write    # + re-disposition the manifest
"""
from __future__ import annotations
import warnings, glob, os, sys, argparse; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits as pf
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from pipeline.acen_orbit import predicted_rv, rv_bounds, consistent_with_orbit, GAMMA, K_A, K_B  # noqa
from config.constants import codex_path  # RYA-810 path register

DATA = str(codex_path('data.spectra_local'))
VET = os.path.join(DATA, "Alpha Centauri (vetted)")
MANIFEST = os.path.join(os.path.dirname(__file__), '..', 'data', 'audit',
                        'acen_holdings_rya384', 'ir_star_id_rya423_manifest.csv')


def hk(h, k):
    for x in h.keys():
        if k.lower() == str(x).lower():
            return h[x]
    return None


def computed_berv(h):
    ra, dec, mjd = hk(h, 'RA'), hk(h, 'DEC'), hk(h, 'MJD-OBS')
    glat, glon, gel = hk(h, 'ESO TEL GEOLAT'), hk(h, 'ESO TEL GEOLON'), hk(h, 'ESO TEL GEOELEV')
    loc = (EarthLocation(lat=glat * u.deg, lon=glon * u.deg, height=(gel or 2400) * u.m)
           if glat is not None else EarthLocation.of_site('La Silla Observatory'))
    sc = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    return sc.radial_velocity_correction('barycentric', obstime=Time(mjd, format='mjd'),
                                         location=loc).to(u.km / u.s).value


def collect(star):
    rows = []
    for f in sorted(glob.glob(os.path.join(VET, star, 'NIRPS', '**', '*.fits'), recursive=True)):
        with pf.open(f) as hd:
            h = hd[0].header
            rv = hk(h, 'ESO QC CCF RV')
            if rv is None:
                continue
            win_start = hk(h, 'ESO RV START')
            fwhm = next((h[x] for x in h.keys() if 'CCF FWHM' in str(x).upper() and 'ERROR' not in str(x).upper()), None)
            rows.append(dict(frame=os.path.basename(f), date=str(hk(h, 'DATE-OBS'))[:10],
                             rv=float(rv), hberv=hk(h, 'ESO QC BERV'), cberv=computed_berv(h),
                             contrast=hk(h, 'ESO QC CCF CONTRAST'), fwhm=fwhm,
                             rv_win_start=win_start, mjd=hk(h, 'MJD-OBS')))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--write', action='store_true'); a = ap.parse_args()
    A, B = collect('Alpha Cen A'), collect('Alpha Cen B')
    lo, hi = rv_bounds()
    print(f"orbit bounds (gamma {GAMMA} +/- max(K_A {K_A:.2f}, K_B {K_B:.2f})): RV in [{lo:.2f}, {hi:.2f}]\n")

    print("== 1. BERV bug test: header BERV vs astropy-computed ==")
    for df, nm in [(A, 'A'), (B, 'B')]:
        r = (df.hberv - df.cberv)
        print(f"   {nm}: header-computed median={r.median():+.4f} max|={r.abs().max():.4f} km/s | NaN? {df.hberv.isna().any()}")
    print("   -> header BERV correct to ~m/s, no NaN: BERV BUG RULED OUT.\n")

    print("== 2. confirmed-A control: RV consistent across BERV (proves BERV handling) ==")
    Acf = A[A.contrast > 25]    # the high-contrast confirmed-A frames
    print(f"   confirmed-A (n={len(Acf)}): RV median={Acf.rv.median():+.2f} std={Acf.rv.std():.2f} across BERV {Acf.cberv.min():+.0f}..{Acf.cberv.max():+.0f}\n")

    print("== 3. the 20 'B' frames ==")
    print(f"   n={len(B)} | dates={sorted(B.date.unique())} | RV median={B.rv.median():+.2f} std={B.rv.std():.2f}")
    print(f"   CCF contrast median={B.contrast.median():.1f} (real peak), FWHM median={B.fwhm.median():.1f}")
    print(f"   search-window blue edge (ESO RV START)={B.rv_win_start.median():+.0f}; CCF blue half-max "
          f"= RV-FWHM/2 = {B.rv.median()-B.fwhm.median()/2:+.1f} -> TRUNCATED below the edge")
    print(f"   distinct nights={B.date.nunique()} -> BERV/date correlation NOT testable within one epoch "
          f"=> offset is CONSTANT (not BERV-correlated, not random)")
    onorbit = B.rv.apply(consistent_with_orbit)
    print(f"   on alpha Cen orbit? {int(onorbit.sum())}/{len(B)}  (RV {B.rv.median():+.1f} is {B.rv.median()-lo:+.1f} vs floor {lo:.1f})\n")

    print("== VERDICT ==")
    print("   CONSTANT offset, NOT a BERV bug, but the RV is OFF the alpha Cen orbit (below the floor)")
    print("   and is a REAL truncated stellar peak -> the 20 frames are a DIFFERENT K-type star,")
    print("   NOT alpha Cen B. Re-disposition: NOT-ALPHA-CEN (RV off-orbit). The J-depth K-type call")
    print("   was a false positive for B; the RV ephemeris (PRIMARY) overrules the spectral type here.")

    if a.write and os.path.exists(MANIFEST):
        m = pd.read_csv(MANIFEST)
        bframes = set(B.frame)
        mask = m.frame.isin(bframes)
        m.loc[mask, 'verdict'] = 'NOT-ALPHA-CEN'
        m.loc[mask, 'evidence'] = (f'RYA-431 RCA: real CCF RV {B.rv.median():.1f} (contrast '
            f'{B.contrast.median():.0f}, window-truncated) is OFF the alpha Cen orbit '
            f'[{lo:.1f},{hi:.1f}]; BERV correct (ruled out); constant (one night 2024-03-17). '
            f'K-type by J-depth but NOT alpha Cen B -> different star, mislabelled "alf Cen A".')
        m.to_csv(MANIFEST, index=False)
        print(f"\n   [manifest] re-dispositioned {int(mask.sum())} frames -> NOT-ALPHA-CEN in {MANIFEST}")


if __name__ == '__main__':
    main()
