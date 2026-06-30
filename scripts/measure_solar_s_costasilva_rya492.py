"""
scripts/measure_solar_s_costasilva_rya492.py
============================================
RYA-492 — settle solar S's +0.40-over-Asplund offset by adopting the Costa Silva+2020
(A&A 634, A136) atlas-tuned log gf for S I Multiplet 8, then re-measuring.

RYA-491 measured A(S) ~7.52 on the CORRECT lines (6743.5 + 6757.1, after dropping the
Costa-Silva-discarded 6748) but with the GES synthesis gf, and flagged the residual as
"it's the gf". This script adopts the published Costa Silva gf and re-measures, honestly,
against the same solar atlas.

PROVENANCE (NOT fabricated) — Costa Silva+2020 Table 1, pulled verbatim from the arXiv
LaTeX source (arXiv:1912.08659, AA2019_36523.tex, the table 'atomic:parameters'). All
S I, Multiplet 8, EP from the table:

  lambda(air)  EP     log gf   reference
  6743.4834    7.866  -1.27    emp. (Costa Silva adapted, from Kurucz2004/Wiese1969 base)
  6743.5400    7.866  -0.95    Kurucz2004
  6743.6401    7.866  -0.93    emp. (Costa Silva adapted)
  6756.7500    7.870  -1.67    Kurucz2004
  6757.0000    7.870  -0.83    Kurucz2004   (orig 6756.960, wavelength shifted)
  6757.1750    7.870  -0.24    Kurucz2004   (orig 6757.150, wavelength shifted)
The 6748 line was DISCARDED by Costa Silva (incorrect model fitting); we already dropped it.

KEY PRE-ANALYSIS (surfaced, not hidden): our GES synthesis list already carries the 6757
triplet at Costa-Silva-IDENTICAL gf (6756.75 -1.67, 6756.96 -0.83, 6757.15 -0.24), so that
window does not move. The only gf changes are the 6743 triplet (our -1.310/-1.066/-0.957 ->
CS -1.27/-0.95/-0.93). Those are SMALL (+0.04..+0.12 dex, lines get slightly stronger ->
A drops slightly). NOTE a grid offset: our third 6743 component sits at 6743.580, Costa
Silva's at 6743.6401 (0.06 A) — we adopt CS's gf onto the nearest component and flag it.
Mechanically this cannot move A(S) by 0.40 dex; the expectation per the ticket's own
escalation branch is a residual -> RYA-161/162 astrophysical-gf-floor, because Costa Silva's
7.12 is an ATLAS+CODE-tuned differential zero-point (MOOG + Kurucz ATLAS), not a transferable
absolute gf. We run it to confirm empirically.

Scope: analysis/scratch only — no canonical_gf / gold-reference edits unless it lands ~7.12.
STOP at the verdict.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.cno_synthesis as cs                                          # noqa: E402
from scripts.procyon_oi777_continuum_localize_rya478 import (               # noqa: E402
    _resources, SOLAR_VMAC, SOLAR_VSINI)
from config.constants import get_star_params                                # noqa: E402


def _params(star='solar'):
    rec = get_star_params(star)
    return {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
            'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}

RESULTS = ROOT / 'data' / 'results'
TMP = '/tmp/ispec_solar_s_rya492'
S_WINDOWS = [(6743.53, 0.35), (6757.15, 0.35)]        # same windows as RYA-491 track3_S
ASPLUND_S = 7.12                                       # Costa Silva calibration target / Asplund 2021

# Costa Silva+2020 Table 1 — {nearest-our-wavelength: (cs_lambda, cs_loggf, cs_ref)}.
# Matched onto the loaded GES S I components by nearest wavelength within each EP group.
COSTA_SILVA_GF = {
    'mult8_6743': [(6743.4834, -1.27, 'emp'), (6743.5400, -0.95, 'Kurucz2004'),
                   (6743.6401, -0.93, 'emp')],
    'mult8_6757': [(6756.7500, -1.67, 'Kurucz2004'), (6757.0000, -0.83, 'Kurucz2004'),
                   (6757.1750, -0.24, 'Kurucz2004')],
}


def _rule(c='='):
    print(c * 80)


def _s_rows(ll):
    """Indices of S I Multiplet-8 components (EP ~7.86-7.88) in the two windows."""
    el = np.array([str(x).strip() for x in ll['element']])
    wA = ll['wave_A'].astype(float)
    ep = ll['lower_state_eV'].astype(float)
    m = (el == 'S 1') & (ep > 7.85) & (ep < 7.90) & (
        ((wA > 6743.0) & (wA < 6744.0)) | ((wA > 6756.5) & (wA < 6757.6)))
    return np.where(m)[0]


def override_costa_silva(ll):
    """Transiently set the Mult-8 S I loggf to Costa Silva+2020 (nearest-wavelength map).
    Returns a list of (wl, old, new, cs_lambda, ref) for the audit trail."""
    wA = ll['wave_A'].astype(float)
    idx = list(_s_rows(ll))
    targets = COSTA_SILVA_GF['mult8_6743'] + COSTA_SILVA_GF['mult8_6757']
    changes = []
    used = set()
    for cs_wl, cs_gf, ref in targets:
        cand = [i for i in idx if i not in used]
        j = min(cand, key=lambda i: abs(wA[i] - cs_wl))
        used.add(j)
        old = float(ll['loggf'][j])
        ll['loggf'][j] = cs_gf
        changes.append(dict(our_wl=round(float(wA[j]), 4), cs_wl=cs_wl, old=round(old, 4),
                            new=cs_gf, d=round(cs_gf - old, 4), ref=ref,
                            grid_offset=round(float(wA[j]) - cs_wl, 4)))
    return changes


def fit_S(ll, iso, sab, params, atm, grid, tag):
    s_codes = cs._atom_codes(('S',), cs._load_synth_resources()[2], sab)
    w_nm, flx = cs._load_observed_spectrum('solar')
    per = []
    for wl, hw in S_WINDOWS:
        r = cs._fit_element(w_nm, flx, atm, params, 'S', {'S': 7.2}, s_codes,
                            ((wl - hw, wl + hw),), False, (cs.HARPS_R, SOLAR_VMAC, SOLAR_VSINI),
                            5.8, 8.2, ll, iso, sab, TMP)
        g = grid[(abs(grid['wave_A'] - 6757.15) < 0.1) & (grid['teff_K'] == 5772)]
        d = float(g['delta_nlte'].iloc[0]) if len(g) else -0.015
        per.append(dict(window=wl, a_lte=round(r['A_X'], 3), chi2r=round(r['red_chi2'], 1),
                        nlte=round(d, 4), a_nlte=round(r['A_X'] + d, 3), status=r['status']))
        print(f"    [{tag}] S I {wl:>8}: A_LTE {per[-1]['a_lte']}  chi2r {per[-1]['chi2r']}  "
              f"NLTE {d:+.4f}  A_NLTE {per[-1]['a_nlte']}  [{r['status']}]")
    vals = [p['a_nlte'] for p in per if np.isfinite(p['a_nlte'])]
    a = round(float(np.median(vals)), 3) if vals else None
    sig = round(float(np.std(vals, ddof=1)), 3) if len(vals) > 1 else 0.0
    print(f"    [{tag}] A(S)_NLTE median {a} (sigma {sig})  vs {ASPLUND_S} "
          f"({(a or 0) - ASPLUND_S:+.3f})")
    return dict(per_window=per, a_nlte=a, sigma=sig)


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, sab, codes = _resources()
    cs._CHEM = None
    params = _params('solar')
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    grid = pd.read_csv(ROOT / 'data' / 'nlte_grids' / 'S_Amarsi2025_PySME.csv')

    _rule(); print("  RYA-492 — solar S re-measure: GES gf (control) vs Costa Silva+2020 gf"); _rule()

    print("\n  CONTROL — current GES synthesis gf (reproduce RYA-491 ~7.52):")
    ctrl = fit_S(ll, iso, sab, params, atm, grid, 'GES')

    print("\n  Costa Silva+2020 Table 1 gf override (arXiv:1912.08659, verified):")
    changes = override_costa_silva(ll)
    for c in changes:
        flag = f"  [grid offset {c['grid_offset']:+.3f} A]" if abs(c['grid_offset']) > 0.02 else ""
        print(f"    our {c['our_wl']:>9} <- CS {c['cs_wl']:>9} ({c['ref']:>10}): "
              f"loggf {c['old']:+.3f} -> {c['new']:+.3f} (d {c['d']:+.3f}){flag}")

    print("\n  COSTA SILVA — re-fit with adopted gf:")
    cosi = fit_S(ll, iso, sab, params, atm, grid, 'CS')

    _rule(); print("  VERDICT (validate-don't-tune)"); _rule()
    a_ctrl, a_cs = ctrl['a_nlte'], cosi['a_nlte']
    moved = round((a_cs or 0) - (a_ctrl or 0), 3)
    resid = round((a_cs or 0) - ASPLUND_S, 3)
    settled = a_cs is not None and abs(resid) <= 0.10
    print(f"  A(S): GES {a_ctrl}  ->  Costa Silva gf {a_cs}  (moved {moved:+.3f})")
    print(f"  residual vs 7.12 target: {resid:+.3f}")
    if settled:
        print("  VERDICT: S SETTLED on Costa Silva gf -> bankable measured solar S (gf-floor resolved).")
    else:
        print(f"  VERDICT: S STILL HIGH on the published gf (residual {resid:+.3f}). The gf swap "
              f"accounts for only {moved:+.3f} dex — Costa Silva's 7.12 is an ATLAS+CODE-tuned "
              f"differential zero-point (MOOG + Kurucz ATLAS), NOT a transferable absolute gf. "
              f"ESCALATE to RYA-161/162 as a genuine astrophysical-gf-floor residual (like Cu). "
              f"Do NOT tune. S stays measured-not-banked.")

    out = {
        'ticket': 'RYA-492 solar S Costa Silva gf re-measure',
        'gf_provenance': 'Costa Silva+2020 A&A 634 A136 Table 1 (arXiv:1912.08659 LaTeX source, verified)',
        'control_ges': ctrl, 'costa_silva': cosi,
        'gf_changes': changes,
        'a_s_ges': a_ctrl, 'a_s_costa_silva': a_cs, 'moved': moved,
        'residual_vs_712': resid, 'settled': settled,
        'verdict': 'banked' if settled else 'escalate_rya161_gf_floor',
        'note_6757_window': 'our 6757 triplet already carries Costa-Silva-identical gf -> no change there; '
                            'the entire lever is the 6743 triplet',
        'caveat_grid_offset': 'our third 6743 component is at 6743.580 vs Costa Silva 6743.6401 (0.06 A); '
                              'CS gf adopted onto the nearest component, flagged',
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'solar_s_costasilva_rya492.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'solar_s_costasilva_rya492.json'}")
    return out


if __name__ == '__main__':
    main()
