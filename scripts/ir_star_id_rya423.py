"""
scripts/ir_star_id_rya423.py
============================
RYA-423 — IR-native star-ID + blend vetting for the alpha Cen CRIRES/NIRPS frames, which the
RYA-384 optical (HARPS-template) classifier cannot reach.

PRIMARY  (model-independent): absolute RV vs the alpha Cen AB orbital ephemeris
         (pipeline.acen_orbit, Kervella 2016). NIRPS supplies a pipeline CCF RV per frame;
         CRIRES does not (and its reduced spectra are telluric-dominated -> RV not recoverable
         per-frame without telluric correction, the same wall RYA-384 hit).
SECONDARY (corroboration only -- NOT the arbiter, to avoid CO/abundance circularity): the
         RYA-384 spectral-type signal (J-band line depth + calibrated flux), and the CCF
         contrast (K-mask deep vs G-mask).
BLEND/INDETERMINATE: low CCF contrast (few lines -> possible hot standard, e.g. 'Star S5'),
         or an absolute RV inconsistent with BOTH predicted stars, is flagged loudly.

NIR caveat the data forced: NIRPS absolute RV carries a MASK-DEPENDENT zero-point (G vs K mask),
so the absolute RV cleanly confirms the G-type (A) branch but is offset for the K-type (B)
branch -- so RV PRIMARY confirms A + the binary separation; the spectral-type SECONDARY carries
B, with the RV offset flagged, never silently overridden.

    python scripts/ir_star_id_rya423.py
"""
from __future__ import annotations
import warnings, glob, os, sys; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from astropy.io import fits as pf
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from pipeline.acen_orbit import predicted_rv, consistent_with_orbit, rv_bounds, SOURCE, GAMMA, K_A, K_B   # noqa: E402
from config.constants import codex_path  # RYA-810 path register

DATA = str(codex_path('data.spectra_local'))
VET = os.path.join(DATA, "Alpha Centauri (vetted)")
OUT = os.path.join(os.path.dirname(__file__), '..', 'data', 'audit', 'acen_holdings_rya384')
RV_TOL = 2.5     # km/s; a frame "matches" a predicted star within this
CONTRAST_MIN = 8.0   # below this the CCF barely caught lines -> possible standard / blend


def hk(h, k):
    for x in h.keys():
        if k.lower() == str(x).lower():
            return h[x]
    return None


def jdepth(f):
    try:
        with pf.open(f) as hd:
            for e in hd[1:]:
                d = e.data
                if d is not None and hasattr(d, 'names') and d.names:
                    nm = {n.upper(): n for n in d.names}
                    wk = next((nm[k] for k in ('WAVE', 'WAVELENGTH') if k in nm), None)
                    fk = next((nm[k] for k in ('FLUX', 'FLUX_REDUCED') if k in nm), None)
                    if wk and fk:
                        w = np.asarray(d[wk]).ravel().astype(float); fx = np.asarray(d[fk]).ravel().astype(float)
                        wu = w / 1e4 if np.nanmedian(w) > 1e4 else w / 1000.
                        m = (wu > 1.20) & (wu < 1.235); g = m & np.isfinite(fx) & (fx > 0)
                        if g.sum() < 50:
                            return None
                        return float(np.mean(1 - np.clip(fx[g] / np.percentile(fx[g], 95), 0, 1)))
    except Exception:
        return None


def verdict(rv, pa, pb, jd, contrast):
    """Combine PRIMARY (RV) + SECONDARY (J-depth/flux/contrast) into A / B / NOT-ALPHA-CEN /
    INDETERMINATE.

    RYA-431 correction: a REAL absolute RV that lies OUTSIDE the hard orbit bounds gamma +/-
    max(K) (acen_orbit.rv_bounds) means the frame is NOT a bound alpha Cen member -- a NIRPS
    mask zero-point cannot move a precision-RV CCF by >5 km/s, so an off-orbit RV is a real
    velocity, and the spectral type CANNOT override an orbital impossibility. This supersedes
    the old "NIR mask zero-point, flagged not overridden" branch (the 20 2024-03-17 "B" frames
    sit at -34.6, 6 km/s below the floor -> a different K star mislabelled 'alf Cen A')."""
    st = 'A' if (jd is not None and jd < 0.45) else ('B' if (jd is not None and jd > 0.55) else '?')
    lo, hi = rv_bounds()
    rv_match = None
    if rv is not None and np.isfinite(rv):
        if abs(rv - pa) <= RV_TOL and abs(rv - pa) < abs(rv - pb):
            rv_match = 'A'
        elif abs(rv - pb) <= RV_TOL and abs(rv - pb) < abs(rv - pa):
            rv_match = 'B'
        elif not consistent_with_orbit(rv):
            rv_match = 'OFF-ORBIT'      # outside gamma +/- max(K) -> not a bound alpha Cen member
    # flags
    if contrast is not None and contrast < CONTRAST_MIN:
        return 'INDETERMINATE', f'low CCF contrast {contrast:.1f} (few lines: possible hot standard/blend); spec-type={st}; RV={rv}'
    if rv_match == 'OFF-ORBIT':
        return 'NOT-ALPHA-CEN', (f'PRIMARY RV {rv:+.1f} is OUTSIDE the alpha Cen orbit bounds '
                    f'[{lo:.1f},{hi:.1f}] (gamma {GAMMA} +/- max(K) {max(K_A,K_B):.1f}) -- a real '
                    f'velocity (NIRPS mask cannot offset a CCF by >5 km/s), so NOT a bound member; '
                    f'spec-type={st} (K-type but a DIFFERENT star). RV ephemeris overrules spec-type.')
    if rv_match == 'A' and st in ('A', '?'):
        return 'A', f'PRIMARY RV matches A ({rv:+.1f} vs {pa:+.1f}); spec-type={st}: AGREE'
    if rv_match == 'B' and st in ('B', '?'):
        return 'B', f'PRIMARY RV matches B ({rv:+.1f} vs {pb:+.1f}); spec-type={st}: AGREE'
    if rv_match and st in ('A', 'B') and rv_match != st:
        return 'INDETERMINATE', f'PRIMARY-vs-SECONDARY DISAGREE: RV->{rv_match}, spec-type->{st} [human review]'
    if st in ('A', 'B'):
        return st, f'spec-type={st} (RV unavailable/inconclusive)'
    return 'INDETERMINATE', f'no clean discriminator (RV={rv}, spec-type={st})'


def main():
    rows = []
    # ---- NIRPS (pipeline CCF RV) ----
    for star_dir, _ in [('Alpha Cen A', 'A'), ('Alpha Cen B', 'B')]:
        for f in sorted(glob.glob(os.path.join(VET, star_dir, 'NIRPS', '**', '*.fits'), recursive=True)):
            with pf.open(f) as hd:
                h = hd[0].header
                mjd = hk(h, 'MJD-OBS'); rv = hk(h, 'ESO QC CCF RV'); con = hk(h, 'ESO QC CCF CONTRAST')
                lab = str(hk(h, 'ESO OBS TARG NAME')); dt = str(hk(h, 'DATE-OBS'))[:10]
            if mjd is None:
                continue
            p = predicted_rv(mjd); jd = jdepth(f)
            v, ev = verdict(rv, p['rv_A'], p['rv_B'], jd, con)
            rows.append(dict(instr='NIRPS', frame=os.path.basename(f), date=dt, header_label=lab,
                             obs_rv=round(rv, 2) if rv is not None else None,
                             pred_A=round(p['rv_A'], 2), pred_B=round(p['rv_B'], 2),
                             contrast=round(con, 1) if con is not None else None,
                             jdepth=round(jd, 3) if jd is not None else None,
                             verdict=v, evidence=ev))
    # ---- CRIRES (no pipeline RV; reduced spectra telluric-dominated) ----
    #
    # 🔴 RYA-972 — THIS GLOB FAILED SILENTLY OFF-SIRIUS, AND THE DIAGNOSIS WAS WRONG.
    #
    # RYA-972 was filed as "globs a misspelt CHIRES dir -> never executed, silent empty
    # glob", and RYA-971 carried that forward as "the CRIRES star-ID branch is dead code".
    # Measured: `CHIRES` IS the directory's name on disk, the code matches it, and on
    # Sirius this glob returns 16 FITS files. The branch is not dead and the spelling is
    # not the fault — renaming it to CRIRES would BREAK the one machine where it works.
    #
    # The real fault is environmental and silent: `data.spectra_local` resolves to
    # /mnt/codex-data/spectra/_mac_import_20260816, which exists on Sirius and NOT on the
    # Mac. Same code, 16 matches there, 0 here — and nothing said so. A run off-Sirius
    # simply contributed no CRIRES rows and reported success, which is how "never executed"
    # became a believed fact about the code rather than about the machine.
    #
    # Loud now, on the same principle as the Kitt Peak atlas loader: name the missing
    # directory rather than let every frame report as absent (RYA-833 — an absence is a
    # hypothesis, never a conclusion).
    _crires_dir = os.path.join(DATA, 'Alpha Centauri A/CHIRES')
    if not os.path.isdir(_crires_dir):
        raise SystemExit(
            f"CRIRES frame directory not found: {_crires_dir}\n"
            f"`data.spectra_local` resolves to {DATA!r}, which does not exist here — the "
            f"alpha Cen frames live on Sirius. Run this there, or stage the data. "
            f"Refusing to emit a star-ID table with the CRIRES branch silently empty: "
            f"that is what made RYA-972 look like a spelling bug for a release.")
    for f in sorted(glob.glob(os.path.join(_crires_dir, '*.fits'))):
        with pf.open(f) as hd:
            h = hd[0].header
            lab = str(hk(h, 'ESO OBS TARG NAME')); dt = str(hk(h, 'DATE-OBS'))[:10]
            setid = next((str(h[x]) for x in h.keys() if 'WLEN ID' in str(x)), '')
        # CRIRES: no pipeline RV (PRIMARY unavailable) and the reduced spectra are telluric-
        # dominated -> spectral type not separable (RYA-384). CO is corroboration ONLY (the
        # ticket's no-circularity rule: CO is the science, not the ID arbiter), so even the
        # 'Star S5' K-band deep CO does NOT promote to B on its own -> INDETERMINATE.
        co_note = ' [K-band deep CO suggests a cool star, corroboration-only, NOT an ID]' if 'S5' in lab and setid.startswith('K') else ''
        v = 'INDETERMINATE'
        ev = (f'CRIRES {setid}: PRIMARY (RV) unavailable (no pipeline RV; telluric-dominated reduced '
              f'spectrum); spectral type not separable (RYA-384); CO is corroboration-only.{co_note} '
              f'label={lab} -> needs telluric correction + IR templates (downstream gate).')
        rows.append(dict(instr='CRIRES', frame=os.path.basename(f), date=dt, header_label=lab,
                         obs_rv=None, pred_A=None, pred_B=None, contrast=None, jdepth=None,
                         verdict=v, evidence=ev))

    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, 'ir_star_id_rya423_manifest.csv'), index=False)
    print(f"orbit: {SOURCE}")
    print(f"derived K_A={K_A:.2f} K_B={K_B:.2f} gamma={GAMMA}\n")
    print("=== per-frame verdict tally ===")
    for instr in ('NIRPS', 'CRIRES'):
        sub = df[df.instr == instr]
        print(f"  {instr}: " + "  ".join(f"{k}={v}" for k, v in sub['verdict'].value_counts().items()))
    print("\n=== NIRPS detail (one row per date/label group) ===")
    g = df[df.instr == 'NIRPS'].groupby(['date', 'header_label', 'verdict']).agg(
        n=('frame', 'size'), obs_rv=('obs_rv', 'median'), pred_A=('pred_A', 'median'),
        pred_B=('pred_B', 'median'), contrast=('contrast', 'median'), jdepth=('jdepth', 'median')).reset_index()
    print(g.to_string(index=False))
    print(f"\n[manifest] {os.path.join(OUT, 'ir_star_id_rya423_manifest.csv')}")
    ind = df[df.verdict == 'INDETERMINATE']
    print(f"\nINDETERMINATE frames: {len(ind)} (enumerated in manifest)")


if __name__ == '__main__':
    main()
