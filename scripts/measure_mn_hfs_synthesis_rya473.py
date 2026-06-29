#!/usr/bin/env python3
"""
scripts/measure_mn_hfs_synthesis_rya473.py
==========================================
RYA-473 — MEASURE solar Mn via HFS-resolved synthesis on the Den Hartog e6S→z6P
triplet (6013.51 / 6016.67 / 6021.82). The REAL unlock: RYA-468 found — empirically,
on the solar spectrum — that the gf grade was never Mn's final blocker; SATURATION is.

THE BLOCKER (RYA-468, empirical): these three lines now carry a graded gf (Den Hartog
+2011, MED tier, on the canonical table) yet SAT-cull through the EW path — measured
EW ~90-100 mÅ, REW −4.78..−4.82, over the −4.90 saturation knee; each is HFS-split
(hfs_n=6). A single-profile EW fit over-saturates them → KEPT=0, Mn stays no-value. The
gf is correct; the lines needed the right *method*. Same root cause as Cu/Eu/Li: the
WRONG MEASUREMENT TOOL for an HFS element.

THE FIX (reuse, not rebuild): synthesise. The GES atomic line list the production
synthesis path already loads (`GESv6_atom_hfs_iso`, RYA-353) carries the ⁵⁵Mn hyperfine
components of the e6S→z6P triplet NATIVELY — 6 resolved components per feature, gf-summing
to exactly the Den Hartog values (−0.354 / −0.181 / −0.054). Those components ARE the Den
Hartog Lawler Sobeck Sneden Cowan 2011 (ApJS 194,35) hfs pattern (GES reference tag
"DLSSC" / "Den Hartog…"), i.e. Den Hartog Table 4 — cited, NOT invented here. A flux-space
synthesis fit (`cno_synthesis._fit_element`, the RYA-338 synthesis-v2 machinery, the same
path RYA-411 validated for Mn and RYA-466 used for Cu) handles the HFS the EW path could
not, varying only A(Mn) per line and reading χ²ᵣ against the observed HARPS solar flux.
No new synthesis engine; no new line data.

gf PROVENANCE (verify, don't re-edit): RYA-468 already adjudicated the canonical-table
6013/16/21 to Den Hartog (MED tier). The synthesis reads the GES list, which independently
already carries the DLSSC hfs gf — so the synthesis gf-sum is checked here against the
RYA-468 canonical values (must match within rounding); we report it, we do not transcribe
or tune it.

NLTE:
  * The MPIA / Bergemann Mn NLTE grid (data/nlte_grids/Mn_Bergemann_MPIA.csv, on disk) is
    the registered Mn NLTE source; its solar Mn correction on its reference lines is
    small + positive (~+0.107 dex).
  * RYA-411 finding (carried as a LOUD caveat): the grid's nodes are the HIGH-EP lines
    4998/6304/6306/6867 — the low-EP strongly-HFS triplet 6013/16/21 is NOT a grid node,
    and its HFS-RESOLVED Amarsi delta differs (HFS desaturation). That line-exact value
    needs the Amarsi GALAH PySME departure grid (nlte_Mn_scatt_pysme.grd), which is offline
    in this worktree.
  * So: TRY the live Amarsi HFS-resolved delta (RYA-411 path); if the .grd is offline,
    APPLY the MPIA grid solar Mn δ (+0.107) as a vendored correction, LOUDLY FLAGGED —
    never silently LTE, never re-fitted to Asplund. The primary reported value is
    A(Mn)_LTE; A(Mn)_NLTE carries the caveat.

VALIDATE-DON'T-TUNE: gf cited (Den Hartog, in the GES list), HFS cited (Den Hartog Table 4,
in the GES list), NLTE from the registered grid. No value pulled toward Asplund
(A(Mn)☉ = 5.42 ref only). A synthesised Mn that still sits high/low is a FINDING.

Out:
  data/audit/mn_hfs_synthesis/solar_mn_hfs_synthesis_rya473.json   (phase_c folds this)
  data/audit/mn_hfs_synthesis/solar_mn_hfs_synthesis_rya473.csv    (per-line table)

    python -m scripts.measure_mn_hfs_synthesis_rya473
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import warnings
from datetime import date
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from config.constants import ISPEC_DIR, HARPS_R, SOLAR_ASPLUND2021, get_star_params  # noqa: E402

sys.path.insert(0, str(ISPEC_DIR))
import ispec  # noqa: E402

from pipeline.abundances_derive import (  # noqa: E402
    _load_atmosphere, _load_synth_resources, _load_observed_spectrum,
    _ISPEC_SOLAR_ABUND_FILE,
)
from pipeline.cno_synthesis import _atom_codes, _solar_A, _fit_element  # noqa: E402

CANON = _REPO / 'data' / 'linelists' / 'canonical_gf.csv'
OUT_DIR = _REPO / 'data' / 'audit' / 'mn_hfs_synthesis'
TMP_DIR = '/tmp/ispec_mn473'

# ── The Den Hartog e6S→z6P triplet (air Å) + per-line fit half-window (Å) ──────
# Graded to Den Hartog+2011 in RYA-468 (canonical −0.354/−0.181/−0.054, MED tier);
# the GES synth list carries the matching 6-component HFS pattern (DLSSC).
MN_LINES = [
    (6013.51, 0.40),
    (6016.67, 0.40),
    (6021.82, 0.40),
]
# RYA-468 canonical (Den Hartog) gf-sum per feature — the value the GES HFS gf must match.
MN_DH_GF = {6013.51: -0.354, 6016.67: -0.181, 6021.82: -0.054}

# MPIA / Bergemann Mn NLTE grid reference lines (the grid's actual nodes — high EP).
MN_MPIA_REF_LINES = [4998.13, 6304.91, 6306.34, 6867.11]
MN_NLTE_REF = ('MPIA SpectrumTools / Bergemann Mn NLTE grid '
               '(data/nlte_grids/Mn_Bergemann_MPIA.csv); solar Mn δ ~+0.107 dex on its '
               'high-EP reference lines (4998/6304/6306/6867). RYA-411: the low-EP HFS '
               'triplet 6013/16/21 is NOT a grid node — its HFS-resolved Amarsi delta '
               'differs (HFS desaturation) and needs the offline nlte_Mn_scatt_pysme.grd.')

_CHEM = None


def _solar_params():
    rec = get_star_params('solar')
    return {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
            'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}, rec


# ── HFS / gf provenance (the GES feature is resolved + carries Den Hartog gf) ──

def verify_hfs_and_gf(ll) -> dict:
    """For each triplet feature: count the GES ⁵⁵Mn hyperfine components within ±0.25 Å
    and gf-sum them. The premise (RYA-468): >1 component (the EW path could not measure
    this) AND the gf-sum matches the RYA-468 Den Hartog canonical value (the synthesis
    independently already uses the DLSSC gf)."""
    notes = np.array([str(x) for x in ll['element']])
    w = np.asarray(ll['wave_A'], float)
    gflog = np.asarray(ll['loggf'], float)
    out = {}
    for wl, _ in MN_LINES:
        m = np.where((notes == 'Mn 1') & (np.abs(w - wl) < 0.25))[0]
        gf_sum = float(np.log10(np.sum(10.0 ** gflog[m]))) if len(m) else float('nan')
        out[round(wl, 3)] = dict(n_components=int(len(m)),
                                 gf_sum_loggf=round(gf_sum, 3),
                                 dh_canonical=MN_DH_GF[wl],
                                 gf_matches_dh=bool(abs(gf_sum - MN_DH_GF[wl]) < 0.02))
    return out


# ── flux-space synthesis fit (single free element) ────────────────────────────

def _measure(ow, of, atm, params, ll, iso, sab, broad, verbose=True):
    codes = _atom_codes(('Mn',), _CHEM, sab)
    A0 = _solar_A(('Mn',), _CHEM, sab)
    state = dict(A0)
    per_line = []
    for wl, hw in MN_LINES:
        t0 = time.time()
        r = _fit_element(ow, of, atm, params, 'Mn', state, codes,
                         ((wl - hw, wl + hw),), False, broad,
                         A0['Mn'] - 1.2, A0['Mn'] + 1.2,
                         ll, iso, sab, TMP_DIR)
        r['line'] = round(wl, 3)
        r['wall_s'] = round(time.time() - t0, 1)
        per_line.append(r)
        if verbose:
            print(f"    Mn {wl:9.3f}  A={r['A_X']}  χ²ᵣ={r['red_chi2']}  "
                  f"σfit={r['sigma_fit']}  npix={r['n_pix']}  [{r['status']}]  ({r['wall_s']}s)")
    return per_line


def _summarise(per_line):
    vals = [r['A_X'] for r in per_line if isinstance(r['A_X'], float) and np.isfinite(r['A_X'])]
    if not vals:
        return dict(n=0, mean=None, median=None, scatter=None)
    return dict(n=len(vals), mean=round(float(np.mean(vals)), 3),
                median=round(float(np.median(vals)), 3),
                scatter=round(float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0, 3))


def mn_nlte_delta(params) -> dict:
    """Try the live Amarsi HFS-resolved Mn delta on the triplet (RYA-411 path); fall back
    to the registered MPIA/Bergemann grid solar δ (on disk) if the Amarsi .grd is offline.
    Never silent, never re-fitted to Asplund."""
    # 1. live Amarsi HFS-resolved (the line-exact value RYA-411 would compute)
    live_fail = None
    try:
        from pipeline import pysme_nlte as P
        from scripts.resolve_camn_stops_rya411 import _build_lines
        P._A_SUN.setdefault('Mn', SOLAR_ASPLUND2021['Mn'])
        lines, _ = _build_lines('Mn', [wl for wl, _ in MN_LINES], hfs=True)
        P.NLTE_LINES['Mn'] = lines
        star = {'teff': params['teff_K'], 'logg': params['logg'],
                'feh': params['feh'], 'vmic': params['vturb_kms']}
        res = P.nlte_delta('Mn', star=star)
        return dict(delta=round(float(res['delta_median']), 3), live=True,
                    source='Amarsi GALAH PySME Mn HFS-resolved delta on the 6013/16/21 '
                           'triplet (live, RYA-411 path); line-exact, NOT the MPIA high-EP '
                           'value.')
    except Exception as exc:
        live_fail = type(exc).__name__         # bind before the except block clears `exc`
    # 2. vendored: MPIA/Bergemann grid solar Mn δ (on-disk CSV), high-EP reference lines
    try:
        from pipeline.nlte_corrections import _mpia_element_delta
        ds = [_mpia_element_delta('Mn', w, params['teff_K'], params['logg'], params['feh'])
              for w in MN_MPIA_REF_LINES]
        ds = [float(d) for d in ds if d is not None and np.isfinite(d)]
        med = round(float(np.median(ds)), 3)
        return dict(delta=med, live=False,
                    source=f'{MN_NLTE_REF} (live Amarsi .grd offline: {live_fail}; '
                           f'MPIA grid solar δ {med:+.3f} on high-EP reference lines applied '
                           f'as a VENDORED correction, FLAGGED — the triplet-exact HFS-resolved '
                           f'δ would differ, RYA-411).')
    except Exception as exc2:
        return dict(delta=None, live=False,
                    source=f'NLTE unavailable (Amarsi .grd offline: {live_fail}; MPIA grid '
                           f'unreadable: {type(exc2).__name__}) — A(Mn)_LTE only, FLAGGED.')


# ── main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    global _CHEM
    ap = argparse.ArgumentParser()
    ap.parse_args(argv)

    Path(TMP_DIR).mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    params, rec = _solar_params()
    print(f"\n{'='*78}\n  RYA-473 — solar Mn via HFS-resolved synthesis on the Den Hartog "
          f"e6S→z6P triplet\n  (Teff={params['teff_K']:.0f} logg={params['logg']:.2f} "
          f"[Fe/H]={params['feh']:+.2f} ξ={params['vturb_kms']:.2f})\n{'='*78}")

    ll, iso, _CHEM = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    ow, of = _load_observed_spectrum('solar')
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    broad = (float(HARPS_R), float(rec.get('vmac', 1.6) or 1.6), float(rec.get('vsini', 1.9) or 1.9))

    # ── 1. HFS + gf provenance (GES carries the Den Hartog Table-4 pattern) ────
    hfs = verify_hfs_and_gf(ll)
    print("\n  ── HFS resolution + gf provenance (GES ⁵⁵Mn components per feature) ──")
    for wl, d in hfs.items():
        print(f"    Mn {wl:9.3f}: {d['n_components']} components, gf-sum {d['gf_sum_loggf']:+.3f} "
              f"(Den Hartog canonical {d['dh_canonical']:+.3f}) "
              f"→ {'MATCH' if d['gf_matches_dh'] else 'MISMATCH!'}")
    assert all(d['n_components'] > 1 for d in hfs.values()), "a Mn feature collapsed to <=1 component"
    assert all(d['gf_matches_dh'] for d in hfs.values()), "GES gf-sum != Den Hartog canonical"
    print("    → all features HFS-resolved (6 comps); GES gf = Den Hartog Table 4 (DLSSC). "
          "EW path SAT-culls these (RYA-468); synthesis measures them.")

    # ── 2. Mn HFS-resolved LTE synthesis fit ──────────────────────────────────
    print("\n  ── Mn I HFS-resolved synthesis fit (flux space) ──")
    mn_pl = _measure(ow, of, atm, params, ll, iso, sab, broad)
    mn_sum = _summarise(mn_pl)
    a_lte = mn_sum['median']
    asp = SOLAR_ASPLUND2021['Mn']

    # ── 3. NLTE ───────────────────────────────────────────────────────────────
    nd = mn_nlte_delta(params)
    a_nlte = (round(a_lte + nd['delta'], 3)
              if (a_lte is not None and nd['delta'] is not None) else a_lte)
    print(f"\n  A(Mn)_LTE  median {a_lte} (mean {mn_sum['mean']}, σ {mn_sum['scatter']}, n={mn_sum['n']})")
    if nd['delta'] is not None:
        print(f"  NLTE Δ {nd['delta']:+.3f} ({'live Amarsi HFS' if nd['live'] else 'vendored MPIA grid'}) "
              f"→ A(Mn)_NLTE {a_nlte}  ({a_nlte - asp:+.3f} vs Asplund {asp:.2f})")
    else:
        print(f"  NLTE unavailable → A(Mn)_LTE only ({(a_lte or 0) - asp:+.3f} vs Asplund {asp:.2f})")
    print(f"    NLTE provenance: {nd['source']}")

    # ── 4. write products ─────────────────────────────────────────────────────
    payload = {
        'ticket': 'RYA-473', 'star': 'solar', 'generated': date.today().isoformat(),
        'engine': 'pipeline.cno_synthesis._fit_element (synthesis v2, RYA-338) on '
                  'GESv6_atom_hfs_iso (HFS-resolved); Turbospectrum 1D-LTE flux fit',
        'params': params,
        'Mn': {
            'method': 'HFS-resolved synthesis (flux-space χ²ᵣ fit) on the Den Hartog '
                      'e6S→z6P triplet 6013/16/21',
            'why_synthesis': 'EW path SAT-culls these (RYA-468: REW −4.78..−4.82 over the '
                             '−4.90 knee, HFS-split hfs_n=6); synthesis on the GES HFS line '
                             'list (6 components/feature) measures them.',
            'gf_provenance': 'Den Hartog Lawler Sobeck Sneden Cowan 2011 ApJS 194,35 '
                             '(MED tier, adjudicated RYA-468); the GES synth list carries the '
                             'matching 6-component hfs pattern (DLSSC = Den Hartog Table 4).',
            'hfs_provenance': 'Den Hartog+2011 Table 4 hfs components (native in GES '
                              'GESv6_atom_hfs_iso) — cited, not invented.',
            'hfs_components': {str(k): v for k, v in hfs.items()},
            'A_lte_median': a_lte, 'A_lte_mean': mn_sum['mean'],
            'scatter': mn_sum['scatter'], 'n_lines': mn_sum['n'],
            'nlte_delta': nd['delta'], 'nlte_live': nd['live'], 'nlte_source': nd['source'],
            'A_nlte': a_nlte, 'asplund2021': asp,
            'delta_vs_asplund': round((a_nlte if a_nlte is not None else a_lte) - asp, 3)
                                if (a_lte is not None) else None,
            'per_line': [{k: r[k] for k in ('line', 'A_X', 'red_chi2', 'sigma_fit',
                                            'n_pix', 'status')} for r in mn_pl],
        },
    }
    out_json = OUT_DIR / 'solar_mn_hfs_synthesis_rya473.json'
    out_json.write_text(json.dumps(payload, indent=2, default=float))

    out_csv = OUT_DIR / 'solar_mn_hfs_synthesis_rya473.csv'
    with open(out_csv, 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['element', 'line_A', 'A_X_lte', 'red_chi2', 'sigma_fit', 'n_pix', 'status'])
        for r in mn_pl:
            wr.writerow(['Mn', r['line'], r['A_X'], r['red_chi2'], r['sigma_fit'],
                         r['n_pix'], r['status']])

    print(f"\n  wrote {out_json.relative_to(_REPO)}\n        {out_csv.relative_to(_REPO)}")
    final = a_nlte if a_nlte is not None else a_lte
    print(f"\n  VERDICT MOVE: Mn no-value → MEASURED via HFS synthesis. "
          f"A(Mn)_{'NLTE' if nd['delta'] is not None else 'LTE'} {final} "
          f"({final - asp:+.3f} vs Asplund {asp:.2f}) n={mn_sum['n']} σ={mn_sum['scatter']}.")
    return payload


if __name__ == '__main__':
    main()
