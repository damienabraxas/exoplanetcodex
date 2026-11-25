import json
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pytest

from pipeline.multi_atom import read_multi_atom, to_lightweaver
from scripts.rya1282_extend_oxygen_atom import build


def test_radiative_extension_preserves_baseline_and_resolves_target_components(tmp_path):
    root = Path.home() / 'codex'
    base_path = root / 'bride_inputs/rya1281/atom.o41f'
    li = root / 'rya1395/data/reference/li2023_oi/raw'
    if not all(p.is_file() for p in (base_path, li / 'tablea1.tsv', li / 'tablea2.tsv')):
        pytest.skip('Held O41f or Li 2023 source tables are absent')

    path = build(base_path, li / 'tablea1.tsv', li / 'tablea2.tsv', tmp_path / 'trial')
    base, expanded = read_multi_atom(base_path), read_multi_atom(path)
    assert len(expanded.levels) == 53
    assert len(expanded.lines) == 399
    assert expanded.levels[:39] == base.levels[:39]
    assert expanded.levels[51:] == [replace(level, index=level.index + 12)
                                     for level in base.levels[39:]]
    assert expanded.lines[:256] == base.lines
    for block in ('continua', 'collisions'):
        originals, remapped = getattr(base, block), getattr(expanded, block)
        assert len(originals) == len(remapped)
        for before, after in zip(originals, remapped):
            for key, value in asdict(before).items():
                expected = value + 12 if key in ('i', 'j') and value >= 40 else value
                np.testing.assert_equal(getattr(after, key), expected)

    receipt = json.loads(path.with_name('receipt.json').read_text())
    assert receipt['status'] == 'RADIATIVE_ONLY_HOLD'
    model = to_lightweaver(expanded)
    windows = [(369.3, 369.4), (1589.1, 1589.8), (1687.0, 1688.0), (1745.7, 1746.0)]
    added = [line.lambda0 for line in model.lines[256:]]
    assert [sum(lo < w < hi for w in added) for lo, hi in windows] == [3, 9, 3, 6]


def test_norad_ground_target_covers_each_added_level(tmp_path):
    root = Path.home() / 'codex'
    base_path = root / 'bride_inputs/rya1281/atom.o41f'
    li = root / 'rya1395/data/reference/li2023_oi/raw'
    norad = root / 'bride_mac_checks/rya1282-o53-radiative-v2'
    inputs = (base_path, li / 'tablea1.tsv', li / 'tablea2.tsv',
              norad / 'o1.en.ls.txt', norad / 'o1.ptpx.txt')
    if not all(path.is_file() for path in inputs):
        pytest.skip('Held O41f, Li 2023 or NORAD source tables are absent')

    path = build(base_path, li / 'tablea1.tsv', li / 'tablea2.tsv', tmp_path / 'trial', norad=norad)
    expanded = read_multi_atom(path)
    old = expanded.continua[:22]
    new = expanded.continua[22:]
    wavelength_floor = min(c.wavelength[-1] for c in old)
    assert {c.i for c in new} == set(range(40, 52))
    assert len(new) == 12
    for continuum in new:
        assert continuum.j == 52
        threshold = expanded.level(continuum.j).E_cm - expanded.level(continuum.i).E_cm
        assert continuum.wavelength[0] == pytest.approx(1e8 / threshold)
        assert continuum.alpha0 == continuum.alpha[0]
        assert continuum.nlamb > 2
        assert continuum.wavelength[-1] >= wavelength_floor
    assert len(to_lightweaver(expanded).continua) == 34


def test_ls_hydrogen_excitation_round_trips_through_multi_ch(tmp_path):
    from scripts.rya1282_add_oxygen_hydrogen_rates import CM_TO_K, add_rates

    path = Path.home() / 'codex/bride_mac_checks/rya1282-o53-radiative-v4/atom.o53-radiative'
    if not path.is_file():
        pytest.skip('Radiative O53 artifact is absent')
    rates = tmp_path / 'rates.txt'
    rates.write_text('lower upper 4000 7000\n4p.5P 5d.5Do 1.2e-12 3.0e-12\n')
    previous = read_multi_atom(path)
    atom = read_multi_atom(add_rates(path, rates, tmp_path / 'trial', expected_pairs=1))
    added = atom.collisions[len(previous.collisions):]
    assert len(added) == 15
    for i in (24, 25, 26):
        low = atom.level(i)
        excitation = np.zeros(2)
        for coll in added:
            if coll.i != i:
                continue
            high = atom.level(coll.j)
            excitation += coll.rates * high.g / low.g * np.exp(
                -(high.E_cm - low.E_cm) * CM_TO_K / coll.temperature)
        np.testing.assert_allclose(excitation, [1.2e-12, 3.0e-12], rtol=1e-8, atol=0)
    assert len(to_lightweaver(atom).continua) == 34


def test_van_regemorter_neutral_rate_uses_multi_ce_units():
    from scripts.rya1282_add_oxygen_electron_rates import multi_ce, neutral_gaunt
    from scripts.rya1282_add_oxygen_hydrogen_rates import CM_TO_K

    np.testing.assert_allclose(neutral_gaunt([.01, .1, 1, 4, 10]), [1.160, .493, .100, .040, .023])
    delta = 6000.0
    temperature = delta * CM_TO_K  # y = 1, so P(y) = 0.100 for a neutral target.
    wavelength_a = 1e8 / delta
    f, g_low, g_high = .02, 3, 5
    spontaneous_a = 6.670e15 * (g_low / g_high) * f / wavelength_a**2
    expected_down = 20.60 * (wavelength_a * 1e-8)**3 * spontaneous_a * .100 / np.sqrt(temperature)
    multi_down = multi_ce(delta, f, temperature) * (g_low / g_high) * np.sqrt(temperature)
    assert multi_down == pytest.approx(expected_down, rel=1e-12, abs=0)


def test_forbidden_lines_reuse_ground_configuration_levels(tmp_path):
    from scripts.rya1282_add_forbidden_oxygen import build

    root = Path(__file__).resolve().parents[1]
    primary = root / 'data/linelists/primary_gf'
    source = Path.home() / 'codex/bride_inputs/rya1281/atom.o41f'
    if not source.is_file():
        pytest.skip('Held O41f atom is absent')
    nist = [primary / f'nist_asd_OI_{span}.tsv' for span in ('900_3000', '3000_25000')]
    out = build(source, nist, tmp_path / 'forbidden')
    old, new = read_multi_atom(source), read_multi_atom(out)
    assert new.levels == old.levels
    assert new.lines[:-3] == old.lines
    assert [(t.j,t.i) for t in new.lines[-3:]] == [(5,4),(4,2),(4,3)]
    assert new.lines[-3].gRad == pytest.approx(1.34311735)
    assert new.lines[-2].f == pytest.approx(1.84343e-11, abs=0)
    assert len(to_lightweaver(new).lines) == 259
