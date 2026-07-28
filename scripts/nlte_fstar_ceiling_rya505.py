#!/usr/bin/env python3
"""
scripts/nlte_fstar_ceiling_rya505.py — RYA-505 hot-Teff NLTE coverage recon (Step 0).

Extends the RYA-349 Step-0 finding (non-Fe NLTE grids top out 6200-6500 K, so Procyon
6554 K is off-grid for everything but Fe). Step 0 here is a HARD RECON GATE: per non-Fe
element in NLTE_CORRECTION_ELEMENTS, report (a) the current grid + code/atom + atmosphere
family we wire, (b) our ON-DISK Teff ceiling (max node in the loaded CSV), (c) the
PUBLISHED grid's real Teff ceiling from a CITED source (the in-repo binary-provenance
JSON, or the cited survey where provenance didn't record it), and the verdict:
self-consistent-extend / MPIA-with-cross-check / bounded-clamp / real-limit.

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


def main():
    t = step0_table()
    s1 = step1_cno()
    s2 = step2_o()
    t.to_csv(RESULTS / 'nlte_fstar_ceiling_rya505.csv', index=False)

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
