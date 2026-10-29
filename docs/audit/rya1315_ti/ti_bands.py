import sys, numpy as np, pandas as pd
sys.path.insert(0, '.'); sys.path.insert(0, 'scripts')
from pipeline.gf_grades import is_gf_graded
from measure_band_ew import load_window_ex
c = pd.read_csv('data/linelists/canonical_gf.csv', low_memory=False)
t = c[c.species == 'Ti I']; g = t[is_gf_graded(t)].copy()
# band -> list of (instrument, holding)
B = [('near-UV', 3000, 3780, [('kpno_solar_atlas', 'solar_kpno_molecfit_corrected')]),
     ('VIS', 3780, 6910, [('kpno_solar_atlas', 'solar_kpno_molecfit_corrected')]),
     ('red-opt', 6910, 9199, [('iag_fts_solar_atlas', 'solar_iag'), ('kpno_solar_atlas', 'solar_kpno_molecfit_corrected')]),
     ('NIR-Y', 9199, 11159.9, [('kpno_solar_atlas', 'solar_kpno_molecfit_corrected')]),
     ('J', 11159.9, 13000, [('crires_plus', 'solar_crires_plus_j_rya1219'), ('kpno_solar_atlas', 'solar_kpno_molecfit_corrected')]),
     ('H', 15007, 17494, [('crires_plus', 'solar_crires_plus_h_rya1232')]),
     ('K', 19452, 24846, [('crires_plus', 'solar_crires_plus_k_rya1219')])]
TH = 5040 / 5772
out = []
for band, lo, hi, hold in B:
    x = g[g.wavelength_air_A.between(lo, hi)]
    for _, r in x.iterrows():
        w = float(r.wavelength_air_A)
        row = dict(band=band, wave=round(w, 3), ep=r.excitation_potential_eV, loggf=r.log_gf,
                   src=r.lab_source_tag if isinstance(r.lab_source_tag, str) else r.gf_tier,
                   sigma=r.gf_sigma_dex, strength=round(r.log_gf - TH * r.excitation_potential_eV, 2))
        for inst, h in hold:
            try:
                win = load_window_ex(inst, w, 1.0, holding=h)
                wl, f = np.asarray(win[0] if isinstance(win, tuple) else win.wave), None
            except Exception as e:
                row[f'err_{h}'] = str(e)[:80]; continue
            wl = np.asarray(win.wave if hasattr(win, 'wave') else win[0], float)
            f = np.asarray(win.flux if hasattr(win, 'flux') else win[1], float)
            ok = np.isfinite(f)
            wl, f = wl[ok], f[ok]
            if len(wl) < 10: continue
            res = w / 3e5 * 3          # ~3 km/s core window
            core = np.abs(wl - w) <= max(res, 0.02)
            if not core.any(): continue
            cont = np.nanpercentile(f[np.abs(wl - w) <= 1.0], 95)
            depth = 1 - np.nanmin(f[core]) / cont
            near = np.abs(wl - w) <= 0.5
            crowd = 1 - np.nanmean(f[near]) / cont       # mean absorption in +-0.5 A
            # resolved neighbours: local minima deeper than 5% within 0.08-0.4 A
            fm = f / cont; nb = 0
            for i in range(1, len(fm) - 1):
                d = abs(wl[i] - w)
                if 0.08 < d <= 0.4 and fm[i] < fm[i-1] and fm[i] <= fm[i+1] and fm[i] < 0.95:
                    nb += 1
            row.update({'holding': h, 'depth': round(depth, 3), 'crowd': round(crowd, 3),
                        'nbrs': nb, 'cont': round(cont, 3)})
            break
        out.append(row)
df = pd.DataFrame(out)
df.to_csv(sys.argv[1], index=False)
print(df.groupby('band').size())
