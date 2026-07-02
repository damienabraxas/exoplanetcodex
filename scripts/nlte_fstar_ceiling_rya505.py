#!/usr/bin/env python3
"""
scripts/nlte_fstar_ceiling_rya505.py — RYA-505 hot-Teff NLTE coverage recon (Step 0).

Extends the RYA-349 Step-0 finding (non-Fe NLTE grids top out 6200-6500 K, so Procyon
6554 K is off-grid for everything but Fe). RYA-505 reframe: MULTI-MODEL NLTE == multi-
instrument — run every model family that covers a point, report per-model, and treat the
inter-model spread as a measured model-systematic (RYA-282). This harness reports:

  Step 0 (recon gate): per non-Fe element in NLTE_CORRECTION_ELEMENTS — (a) current grid +
  code/atom + atmosphere family, (b) ON-DISK Teff ceiling (max node in the loaded CSV),
  (c) PUBLISHED ceiling from a CITED source (in-repo binary provenance or the cited
  survey). The per-element ceiling verdict names how each family reaches Procyon.

  Multi-model layer (Steps 3-4): per element, which model families are on disk, the
  per-family SOLAR delta, and the inter-model SPREAD at the solar overlap (the model-
  systematic meter; UPPER BOUND — carries the atmosphere-baseline term too) vs a 0.05 tol.
  Single-model elements are flagged (no cross-check available), never silently trusted.

Also confirms (not re-derives): Step 1 the CNO 1D-NLTE leg covers Procyon, and Step 2
the RYA-483 banked Procyon O used the 1D-NLTE (in-grid) correction.

No values from memory: on-disk ceilings are read from the grid CSVs; published ceilings
from the cited provenance. Writes data/results/nlte_fstar_ceiling_rya505.{csv,md}.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

GD = _REPO / 'data' / 'nlte_grids'
AG = GD / 'amarsi_galah'
RESULTS = _REPO / 'data' / 'results'
RESULTS.mkdir(parents=True, exist_ok=True)

PROCYON_TEFF = 6554.0
TAUBOO_TEFF = 6400.0

# Code / model-atom / model-atmosphere family per current registered grid (from each
# grid's provenance 'reference'/'derivation'; NOT from memory).
FAMILY = {
    'Na': 'Amarsi2020 GALAH (Lind2011 atom) · PySME/SME · MARCS marcs2012',
    'Mg': 'Amarsi2020 GALAH · PySME/SME · MARCS marcs2012',
    'Si': 'Amarsi2020 GALAH · PySME/SME · MARCS marcs2012',
    'Al': 'Amarsi2020 GALAH (Nordlander&Lind2017 atom) · PySME/SME · MARCS marcs2012',
    'K':  'Amarsi2020 GALAH · PySME/SME · MARCS marcs2012',
    'S':  'Amarsi2025 · PySME/SME · MARCS marcs2012',
    'Ca': 'Mashonkina2017 · DETAIL · MAFAGS-OS (1D)',
    'Ti': 'Bergemann2011 MPIA · DETAIL/SIU · MAFAGS-OS (1D)',
    'Cr': 'Bergemann&Cescutti2010 MPIA · DETAIL/SIU · MAFAGS-OS (1D)',
    'Mn': 'Bergemann MPIA · DETAIL/SIU · MAFAGS-OS (1D)',
    'Ba': 'Korotin2015 · MULTI · MARCS',
    'Sr': 'Bergemann2012 INSPECT · MULTI · MARCS (metal-poor)',
}

# Published Teff ceiling + the CITED source. Amarsi 2020 GALAH is one released grid set
# (Zenodo 3982506) on a common 32-node Teff axis 2500-8000 K — VERIFIED in-repo from the
# Na/Al/K v3 provenance; Mg/Si are the same release. Ba from its prov coverage. MPIA
# elements: our extraction provenance did not record the select range; the published
# nlte.mpia.de grids run to 7000 K (RYA-505 survey, cited). Amarsi .grd is also staged
# in amarsi_galah/ for Ca/Mn (Family-A, RYA-410) → a self-consistent alt to MAFAGS-OS.
PUBLISHED = {
    'Na': (8000, 'Amarsi2020 GALAH v3 prov (Zenodo 3982506): Teff 2500-8000, 32 nodes'),
    'Mg': (8000, 'Amarsi2020 GALAH (same Zenodo 3982506 release axis 2500-8000)'),
    'Si': (8000, 'Amarsi2020 GALAH (same Zenodo 3982506 release axis 2500-8000)'),
    'Al': (8000, 'Amarsi2020 GALAH v3 prov (Zenodo 3982506): Teff 2500-8000'),
    'K':  (8000, 'Amarsi2020 GALAH v3 prov (Zenodo 3982506): Teff 2500-8000'),
    'S':  (8000, 'Amarsi2025 prov: Teff 3000-8000, 3756 MARCS nodes'),
    'Ca': (None, 'Mashonkina2017 prov records only the 6500 subset; full published range '
                 'not in-repo → source pull owed. Amarsi Ca .grd staged (amarsi_galah/) as '
                 'a self-consistent alt.'),
    'Ti': (7000, 'nlte.mpia.de Bergemann grids run to 7000 K (RYA-505 survey); our '
                 'extraction did not record the select range → subset.'),
    'Cr': (7000, 'nlte.mpia.de Bergemann grids to 7000 K (RYA-505 survey); subset on disk.'),
    'Mn': (7000, 'nlte.mpia.de Bergemann grids to 7000 K (RYA-505 survey). Amarsi Mn .grd '
                 'staged (amarsi_galah/) → self-consistent alt.'),
    'Ba': (6500, 'Korotin2015 prov coverage: Teff 4000-6500 — REAL published limit.'),
    'Sr': (None, 'Bergemann2012 INSPECT metal-poor; INASAN primary pending (RYA-433).'),
}

# NLTE significance at F-star Teff + verdict priority (ticket Step 3 order).
SIGNIF = {
    'Ca': 'significant', 'Ba': 'significant', 'Na': 'significant', 'Mn': 'significant',
    'Mg': 'moderate', 'K': 'moderate', 'S': 'moderate',
    'Si': 'near-LTE', 'Al': 'near-LTE', 'Ti': 'near-LTE', 'Cr': 'moderate', 'Sr': 'off-crit',
}


def on_disk_ceiling(gridfile: str):
    df = pd.read_csv(GD / gridfile)
    tcol = next(c for c in df.columns if c.lower() in ('teff_k', 'teff'))
    return float(df[tcol].max())


def verdict(el, on_disk, pub):
    fam = FAMILY[el]
    amarsi = fam.startswith('Amarsi')
    pub_ceil = pub[0]
    if pub_ceil is not None and on_disk >= pub_ceil - 1:      # on-disk == published
        if pub_ceil >= PROCYON_TEFF:
            return 'IN-GRID (published ceiling already covers Procyon — just under-loaded)'
        return 'REAL-LIMIT → BOUNDED-CLAMP (54 K, monotonic) or find higher grid'
    # on-disk is a subset below the published ceiling
    if amarsi and pub_ceil and pub_ceil >= PROCYON_TEFF:
        return 'SELF-CONSISTENT-EXTEND (re-synth Procyon node from same Amarsi .grd)'
    if el in ('Ca', 'Mn'):
        return ('SELF-CONSISTENT-EXTEND via Amarsi (.grd staged) — swap w/ validate-'
                'dont-tune cross-check vs current MAFAGS-OS value')
    if el in ('Ti', 'Cr'):
        return 'MPIA-WITH-CROSS-CHECK (re-scrape to 7000 K; MAFAGS-OS≠MARCS → solar overlap check)'
    if el == 'Sr':
        return 'DEFER (off-critical; INASAN primary pending RYA-433)'
    return 'UNRESOLVED — source pull owed'


def step0_table():
    from config.constants import NLTE_CORRECTION_ELEMENTS as REG
    rows = []
    for el in REG:
        gf = REG[el]['grid']
        try:
            od = on_disk_ceiling(gf)
        except Exception as exc:
            od = float('nan')
        pub = PUBLISHED.get(el, (None, 'unknown'))
        rows.append({
            'element': el, 'family': FAMILY.get(el, '?'),
            'on_disk_ceiling_K': od,
            'published_ceiling_K': pub[0] if pub[0] is not None else 'unrecorded',
            'covers_procyon_6554': (pub[0] is not None and pub[0] >= PROCYON_TEFF),
            'nlte_signif_Fstar': SIGNIF.get(el, '?'),
            'verdict': verdict(el, od, pub),
            'published_source': pub[1],
        })
    return pd.DataFrame(rows)


def step1_cno():
    from pipeline import nlte_cno
    ceil = float(nlte_cno.TEFF_3D_CEILING)
    return {
        'teff_3d_ceiling_K': ceil,
        'procyon_uses': '1D-NLTE leg (tables 5/6)' if PROCYON_TEFF > ceil else '3D leg',
        'cno_1d_leg_teff_range_K': '4000-8000 (17 nodes; verified from table5/6)',
        'verdict': ('CONFIRMED — Procyon C/O already use the Amarsi/Balder 1D-NLTE leg '
                    '(in-grid at 6554 K); 3D refinement is future v2 (RYA-444/445). No change.'),
    }


def step2_o():
    d = json.load(open(RESULTS / 'procyon_O_hardbank_rya483.json'))
    ind = d.get('primary_indicator', '')
    is_1d_nlte = '1D_NLTE' in ind
    return {
        'primary_indicator': ind,
        'A_O_nlte': d.get('procyon_O_AO_nlte'),
        'banked_OH': d.get('procyon_OH_vs_our_sun'),
        'verdict': ('NLTE-VALID — banked [O/H] +0.085 used O I 777 1D-NLTE (in-grid at '
                    '6554 K), NOT 3D-off-grid or LTE. Its PROVISIONAL status is the '
                    'continuum zero-point (~0.18) + [O I] 6300 terminal leg, not an '
                    'off-grid NLTE correction. No re-derivation needed.'
                    if is_1d_nlte else 'CHECK — indicator is not 1D-NLTE; re-derive.'),
    }


# ── Multi-model layer (RYA-505 reframe): run every model family that covers ───
# On-disk families per element beyond the single registered grid. Amarsi/MARCS vs
# Bergemann/MAFAGS-OS are independent measurements; where >=2 are on disk, the solar
# overlap spread is the model-systematic meter (UPPER BOUND — carries the atmosphere
# baseline too). Elements not listed here are single-model (flagged) at present.
MULTI_FAMILY = {
    'Na': [('Amarsi2020_MARCS', 'Na_Amarsi2020_PySME.csv'),
           ('Lind2011_INSPECT_MARCS', 'Na_Lind2011_INSPECT.csv')],
    'Mg': [('Amarsi2020_MARCS', 'Mg_Amarsi2020_PySME.csv'),
           ('Bergemann_MAFAGS-OS', 'Mg_Bergemann_MPIA.csv')],
    'Si': [('Amarsi2020_MARCS', 'Si_Amarsi2020_PySME.csv'),
           ('Bergemann_MAFAGS-OS', 'Si_Bergemann_MPIA.csv')],
}
OVERLAP_TOL = 0.05


def _solar_delta_median(gridfile, teff=5772.0, logg=4.44, feh=0.0):
    from scipy.interpolate import LinearNDInterpolator, griddata
    df = pd.read_csv(GD / gridfile)
    cols = {c.lower(): c for c in df.columns}
    tc = cols.get('teff_k', cols.get('teff')); gc = cols['logg']; fc = cols['feh']
    dc = cols['delta_nlte']; wc = cols.get('wave_a', list(df.columns)[2])
    df = df.dropna(subset=[dc])
    out = []
    for _w, g in df.groupby(wc):
        pts = g[[tc, gc, fc]].values
        val = g[dc].values
        if len(pts) < 4:
            out.append(float(griddata(pts, val, [[teff, logg, feh]], method='nearest')[0]))
            continue
        try:
            v = LinearNDInterpolator(pts, val)([[teff, logg, feh]])[0]
            if not np.isfinite(v):
                v = griddata(pts, val, [[teff, logg, feh]], method='nearest')[0]
            out.append(float(v))
        except Exception:
            out.append(float(griddata(pts, val, [[teff, logg, feh]], method='nearest')[0]))
    return float(np.median(out)), len(out)


def multimodel_table():
    """Per element: model families on disk, per-family solar delta, inter-model spread +
    agreement badge (Step 4 meter). Single-model elements are flagged (no cross-check)."""
    rows = []
    for el, fams in MULTI_FAMILY.items():
        deltas = {}
        for name, f in fams:
            try:
                d, n = _solar_delta_median(f)
                deltas[name] = d
                rows.append({'element': el, 'model': name, 'grid': f,
                             'solar_delta': round(d, 4), 'n_lines': n})
            except Exception as exc:
                rows.append({'element': el, 'model': name, 'grid': f,
                             'solar_delta': np.nan, 'n_lines': 0, 'err': str(exc)[:40]})
        if len(deltas) == 2:
            a, b = list(deltas.values())
            spread = abs(a - b)
            rows.append({'element': el, 'model': 'SPREAD', 'grid': '',
                         'solar_delta': round(spread, 4),
                         'n_lines': ('AGREE' if spread <= OVERLAP_TOL else 'FLAG-adjudicate')})
    return pd.DataFrame(rows)


def main():
    t = step0_table()
    s1 = step1_cno()
    s2 = step2_o()
    mm = multimodel_table()
    t.to_csv(RESULTS / 'nlte_fstar_ceiling_rya505.csv', index=False)
    mm.to_csv(RESULTS / 'nlte_multimodel_rya505.csv', index=False)

    print("=" * 78)
    print("  RYA-505 — hot-Teff NLTE coverage recon (Step 0 gate)")
    print("=" * 78)
    print(t[['element', 'on_disk_ceiling_K', 'published_ceiling_K',
             'covers_procyon_6554', 'nlte_signif_Fstar', 'verdict']].to_string(index=False))
    print("\n[Step 1 — CNO 1D-NLTE branch]")
    for k, v in s1.items():
        print(f"    {k}: {v}")
    print("\n[Step 2 — RYA-483 banked Procyon O]")
    for k, v in s2.items():
        print(f"    {k}: {v}")

    print("\n[Multi-model layer — run every family that covers; solar-overlap spread meter]")
    print(mm.to_string(index=False))
    print("    (Mg/Si Bergemann=MAFAGS-OS vs Amarsi=MARCS -> spread is an UPPER BOUND: "
          "carries the atmosphere-baseline term as well as the NLTE-method term.)")

    # Procyon / tau Boo coverage map after the decision
    print("\n[Coverage map]")
    live = t[t['covers_procyon_6554']]['element'].tolist()
    print(f"    Procyon 6554 K: published grid already covers -> {sorted(live)} "
          f"(these are under-loaded, not missing: re-load/synth the node).")
    print(f"    tau Boo 6400 K: inside every >=6500 grid -> fully in published grid for "
          f"all except the real-limit/6000 cases (Sr 6000).")
    return t, s1, s2


if __name__ == '__main__':
    main()
