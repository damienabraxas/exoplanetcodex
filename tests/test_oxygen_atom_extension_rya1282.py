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
