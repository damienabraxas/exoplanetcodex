"""Ion identity must survive discovery, level matching and synthesis (RYA-1250)."""
import numpy as np
import pandas as pd
import pytest
from pipeline import pysme_nlte as pn
from pipeline.element_prepare import nlte_has_line


@pytest.fixture
def overlapping_grid(monkeypatch):
    from pipeline import nlte_bfactor_synth as bs
    # Both ions have similar excitation energies; unrestricted matching picks I.
    grid = {'energy': np.array([1., 2., 1.001, 2.001]),
            'J': np.array([.5, 1.5, .5, 1.5]),
            'spec': [b'Mg 1', b'Mg 1', b'Mg 2', b'Mg 2'],
            'conf': [b'neutral_low', b'neutral_up', b'ion_low', b'ion_up'],
            'term': [b'2S', b'2P*', b'2S', b'2P*']}
    monkeypatch.setattr(bs, 'read_amarsi_grid', lambda element: grid)
    return grid


@pytest.mark.parametrize('matcher', [pn.auto_labels, pn.labels_by_selection_rules])
def test_levels_are_matched_inside_requested_ion(overlapping_grid, matcher):
    assert matcher('Mg', 1., 2., ion='II')[:2] == ('ion_low 2S', 'ion_up 2P*')
    assert matcher('Mg', 1., 2.)[:2] == ('neutral_low 2S', 'neutral_up 2P*')
    with pytest.raises(ValueError, match='no levels'):
        matcher('Mg', 1., 2., ion=3)


def line():
    return (7877., .389, 9.99554, .5, 11.5691, 1.5, '2p6.4p 2P*', '2p6.4d 2D', -7.)


def test_ionized_rows_preserve_identity_and_broadening():
    row = pn._linelist_rows('Mg', [line() + (None, {'gamrad': 8., 'gamqst': -5.})], ion='II')[0]
    assert (row['species'], row['ionization']) == ('Mg 2', 2)
    assert (row['gamrad'], row['gamqst'], row['gamvw']) == (8., -5., -7.)
    with pytest.raises(ValueError, match='explicit radiative'):
        pn._linelist_rows('Mg', [line()], ion=2)


def test_neutral_defaults_are_unchanged():
    row = pn._linelist_rows('Mg', [line()])[0]
    assert (row['species'], row['ionization'], row['gamrad'], row['gamqst']) == ('Mg 1', 1, 7.8, 0.)


def test_ionized_calls_cannot_reuse_neutral_diagnostics():
    with pytest.raises(ValueError, match='explicit diagnostic'):
        pn.nlte_delta('Na', ion=2)
    with pytest.raises(ValueError, match='explicit diagnostic'):
        pn.nlte_delta('Na', lines=[5682.633], ion=2)


def test_coverage_and_failed_derivations_are_ion_specific():
    table = pd.DataFrame({'ion': [1], 'wave_A': [5000.]})
    prov = {'not_derivable': {'5100.000': 'legacy neutral'},
            'not_derivable_by_ion': {'2': {'5200.000': 'ionized failure'}}}
    for w in (5000., 5100.):
        assert nlte_has_line(table, prov, 'Mg I', w)
        assert not nlte_has_line(table, prov, 'Mg II', w)
    assert nlte_has_line(table, prov, 'Mg II', 5200.)
    assert not nlte_has_line(table, prov, 'Mg I', 5200.)


def test_explicit_rows_need_no_neutral_registry():
    assert pn._select_lines('Mg', [line()]) == [line()]


@pytest.mark.parametrize('flags', [[], [False], [True, False]])
def test_partial_or_missing_departure_application_fails(flags):
    with pytest.raises(RuntimeError, match='not applied'):
        pn.assert_nlte_engaged(flags, [{}])


def test_departure_application_is_not_inferred_from_nonzero_delta():
    pn.assert_nlte_engaged([True, True], [{}, {}])


def test_extension_reads_the_selected_ion_and_preserves_stark(monkeypatch, tmp_path):
    from scripts import extend_nlte_grid as ext
    monkeypatch.setattr(ext, 'ELEMENT', 'Mg')
    monkeypatch.setattr(ext, 'ION', 2)
    monkeypatch.setattr(ext, 'ROOT', tmp_path)
    folder = tmp_path / 'data/linelists'
    folder.mkdir(parents=True)
    pd.DataFrame({'species': ['Mg I', 'Mg II'], 'wavelength_air_A': [5000., 5000.],
                  'excitation_potential_eV': [1., 2.], 'log_gf': [-1., -.2]}).to_csv(
                      folder / 'canonical_gf.csv', index=False)
    synth = pd.DataFrame({'wave_A': [5000.], 'lower_state_eV': [2.],
                          'upper_state_eV': [4.], 'waals': [-7.], 'rad': [8.],
                          'stark': [-5.], 'lower_j': [.5], 'upper_j': [1.5]})
    monkeypatch.setattr(ext, '_synth_rows', lambda wl: synth)
    result = ext.atomic(5000.)
    assert result[:4] == (-.2, 2., 4., -7.)
    assert result[4] == {'gamrad': 8., 'gamqst': -5.}
    monkeypatch.setattr(ext, '_synth_rows', lambda wl: pd.concat([synth, synth]))
    with pytest.raises(ValueError, match='identified synthesis component'):
        ext.atomic(5000.)


def test_extension_write_keeps_other_ion_and_keys_failures(monkeypatch, tmp_path):
    import json
    import concurrent.futures
    from scripts import extend_nlte_grid as ext
    folder = tmp_path / 'data/nlte_grids'
    folder.mkdir(parents=True)
    path = folder / 'Mg_Amarsi2020_PySME.csv'
    original = {'element': 'Mg', 'ion': 1, 'wave_A': 5000., 'teff_K': 5772,
                'logg': 4.44, 'feh': 0., 'delta_nlte': -.02}
    pd.DataFrame([original]).to_csv(path, index=False)
    report = {'nlte_missing': [{'species': 'Mg II', 'wavelength_A': w} for w in (5000., 5100.)]}
    (tmp_path / 'solar_Mg.json').write_text(json.dumps(report))
    monkeypatch.setattr(ext, 'ROOT', tmp_path)
    monkeypatch.setattr(ext, 'OUT', tmp_path)
    monkeypatch.setattr(ext, 'NODES', [(5772, 4.44, 0.)])
    monkeypatch.setattr(ext, 'build_lines', lambda waves: waves)
    monkeypatch.setattr(ext, '_node', lambda args: {5000.: -.1, 5100.: float('nan')})
    # Serial stand-in tests CSV/provenance, not process scheduling or radiative transfer.
    class Pool:
        def __init__(self, *args): pass
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def map(self, fn, args): return map(fn, args)
    monkeypatch.setattr(concurrent.futures, 'ProcessPoolExecutor', Pool)
    assert ext.main(['--element', 'Mg', '--ion', 'II', '--write']) == 0
    saved = pd.read_csv(path)
    assert saved.iloc[0].to_dict() == original
    assert saved.iloc[1]['ion'] == 2 and saved.iloc[1]['delta_nlte'] == -.1
    prov = json.loads(path.with_suffix('.prov.json').read_text())
    assert '5100.000' in prov['not_derivable_by_ion']['2']
    assert not nlte_has_line(saved, prov, 'Mg I', 5100.)
