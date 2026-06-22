"""
tests/test_nlte_na_mg_ba_rya165.py
==================================
RYA-165 (rescoped) — register Na I / Mg I / Ba II (+ the Mn I / Si I gap-fill) in
NLTE_CORRECTION_ELEMENTS and drop in the vendored grids (Na = Lind 2011 / INSPECT,
Mg/Mn/Si = MPIA MAFAGS-OS, Ba = Korotin 2015).
This ticket is REGISTRATION + grid wiring only: the generic per-line apply machinery is
RYA-235's. Guards here: the vendored grids load + interpolate to their documented
provenance solar anchors; the apply layer corrects Na/Mg (neutral) + Ba (ion II) and
tags each with its CORRECT per-source provenance flag (a non-MPIA grid is never
mislabelled 'MPIA'); and — per the curate-first contract — NLTE on a raw pool is NOT
asserted onto Asplund here (that acceptance is gated on a curated pool, separately).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import nlte_corrections as N            # noqa: E402
from config.constants import NLTE_CORRECTION_ELEMENTS as REG  # noqa: E402

GRIDS = ROOT / 'data' / 'nlte_grids'
SOLAR = dict(teff_K=5777, logg=4.4, feh=0.0)
# RYA-410 re-sourced Na/Mg/Si onto the Amarsi-2020 PySME grids (single source, [Fe/H]->+0.6,
# closing the 55-Cnc clamp); the PySME-Amarsi solar delta reproduced the prior MPIA/INSPECT
# value within tol before the swap. Ba/Mn stayed (Mn: PySME cross-check STOPPED on HFS).
NEW = {
    'Na': {'ion': 1, 'grid': 'Na_Amarsi2020_PySME.csv', 'flag': 'NLTE_Amarsi2020_PySME_1D'},
    'Mg': {'ion': 1, 'grid': 'Mg_Amarsi2020_PySME.csv', 'flag': 'NLTE_Amarsi2020_PySME_1D'},
    'Ba': {'ion': 2, 'grid': 'Ba_Korotin2015.csv',      'flag': 'NLTE_Korotin2015_1D'},
    'Mn': {'ion': 1, 'grid': 'Mn_Bergemann_MPIA.csv',   'flag': 'NLTE_MPIA_MAFAGS_1D'},  # gap-fill (stayed MPIA)
    'Si': {'ion': 1, 'grid': 'Si_Amarsi2020_PySME.csv', 'flag': 'NLTE_Amarsi2020_PySME_1D'},  # negligible, documented
}


# ── registration ──────────────────────────────────────────────────────────────
def test_na_mg_ba_registered_with_correct_ion_and_grid():
    for el, want in NEW.items():
        assert el in REG, f"{el} not registered in NLTE_CORRECTION_ELEMENTS"
        assert int(REG[el]['ion']) == want['ion']
        assert REG[el]['grid'] == want['grid']
        assert (GRIDS / REG[el]['grid']).exists()
        assert (GRIDS / (Path(REG[el]['grid']).stem + '.prov.json')).exists()  # provenance sidecar
        assert REG[el]['ref']


# ── grids load + reproduce their documented provenance solar anchors ──────────
@pytest.mark.parametrize('el', ['Na', 'Mg', 'Ba', 'Mn', 'Si'])
def test_grid_loads_and_matches_provenance_solar_anchor(el):
    prov = json.loads((GRIDS / (Path(REG[el]['grid']).stem + '.prov.json')).read_text())
    anchor = prov['solar_anchor']
    c = N._load_mpia_element_grid(el)
    assert len(c['waves']) >= 2 and c['ion'] == REG[el]['ion']
    # interpolate at the anchor's own stellar point
    ds = [N._mpia_element_delta(el, float(w), anchor['teff_K'], anchor['logg'], anchor['feh'])
          for w in c['waves']]
    ds = [d for d in ds if np.isfinite(d)]
    assert ds, f"{el}: no finite solar delta"
    if 'line_A' in anchor:                                # single-line anchor (Ba)
        d = N._mpia_element_delta(el, float(anchor['line_A']),
                                  anchor['teff_K'], anchor['logg'], anchor['feh'])
        assert abs(d - anchor['delta_nlte']) < 1e-3
    else:                                                 # median anchor (Na/Mg)
        assert abs(float(np.median(ds)) - anchor['delta_nlte']) < 0.02


# ── apply layer tags the correct per-source provenance flag ──────────────────
def _per_line(el, a_1dlte, n=4):
    waves = N._load_mpia_element_grid(el)['waves'][:n]
    ion = 'II' if int(REG[el]['ion']) == 2 else 'I'
    return [{'element': el, 'ion': ion, 'wavelength_air_A': float(w), 'a_1dlte': a_1dlte}
            for w in waves]


def test_apply_flags_each_source_correctly_not_mislabelled_mpia():
    res = pd.DataFrame([
        {'element': 'Na', 'ion': 1, 'A_X': 6.24, 'n_lines': 2},
        {'element': 'Mg', 'ion': 1, 'A_X': 7.55, 'n_lines': 3},
        {'element': 'Ba', 'ion': 2, 'A_X': 2.27, 'n_lines': 4},
    ])
    pl = pd.DataFrame(_per_line('Na', 6.24) + _per_line('Mg', 7.55) + _per_line('Ba', 2.27))
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pl)

    flag = lambda el, ion: out[(out.element == el) & (out.ion == ion)].iloc[0]['nlte_flag']
    assert flag('Na', 1) == 'NLTE_Amarsi2020_PySME_1D'    # RYA-410 re-source; NOT mislabelled MPIA
    assert flag('Ba', 2) == 'NLTE_Korotin2015_1D'         # NOT mislabelled MPIA
    assert flag('Mg', 1) == 'NLTE_Amarsi2020_PySME_1D'    # RYA-410: Mg now Amarsi, not MPIA
    # ref column carries the true citation (Na model atom is still Lind 2011)
    assert 'Lind' in out[out.element == 'Na'].iloc[0]['nlte_ref']
    assert 'Korotin' in out[out.element == 'Ba'].iloc[0]['nlte_ref']


def test_ba_ii_corrected_and_ba_i_left_alone():
    # Ba is the registry's first ion-II element — the ion gate must target II, not I.
    res = pd.DataFrame([
        {'element': 'Ba', 'ion': 2, 'A_X': 2.27, 'n_lines': 4},
        {'element': 'Ba', 'ion': 1, 'A_X': 2.20, 'n_lines': 1},   # not in registry-ion → untouched
    ])
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pd.DataFrame(_per_line('Ba', 2.27)))
    assert out[(out.element == 'Ba') & (out.ion == 2)].iloc[0]['nlte_flag'] == 'NLTE_Korotin2015_1D'
    assert out[(out.element == 'Ba') & (out.ion == 1)].iloc[0]['nlte_flag'] == '1D_LTE'


# ── curate-first contract carried over to the new elements ───────────────────
def test_mn_sizeable_si_clean_solar_nlte():
    # the gap-fill characterisation: Mn I NLTE is real (~+0.10), Si I is negligible
    # (~−0.004 → registered as 'NLTE-clean & documented', not silently 1D-LTE).
    mn = np.median([N._mpia_element_delta('Mn', float(w), 5772, 4.438, 0.0)
                    for w in N._load_mpia_element_grid('Mn')['waves']])
    si = np.median([N._mpia_element_delta('Si', float(w), 5772, 4.438, 0.0)
                    for w in N._load_mpia_element_grid('Si')['waves']])
    assert mn > 0.05                                   # Mn NLTE is sizeable
    assert abs(si) < 0.02                              # Si NLTE is negligible (clean)


def test_validation_still_deferred_not_asserted_onto_asplund():
    # registration does NOT claim Asplund acceptance — applying NLTE to a raw pool just
    # shifts the (possibly un-curated) number. We assert the MACHINERY ran, not science.
    res = pd.DataFrame([{'element': 'Na', 'ion': 1, 'A_X': 6.51, 'n_lines': 2}])  # raw-high Na (239 +0.27)
    out = N.apply_element_nlte_corrections(res, SOLAR, per_line_df=pd.DataFrame(_per_line('Na', 6.51)))
    row = out.iloc[0]
    assert row['nlte_flag'] == 'NLTE_Amarsi2020_PySME_1D'        # applied (RYA-410 re-source)
    assert np.isfinite(row['A_X_nlte']) and row['n_nlte_lines'] >= 1
