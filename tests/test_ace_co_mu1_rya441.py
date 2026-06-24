"""
tests/test_ace_co_mu1_rya441.py
===============================
RYA-441 — guard the disposable mu=1 intensity diagnostic. Fast, deterministic
checks: the intensity-column remap (the correctness-critical piece -- it must NOT
silently substitute flux), and the bsyn-stdin injection / Popen restoration. The
full mu=1 A(C) fit (~45 s, needs the iSpec install) is the CLI smoke test
`python scripts/ace_co_mu1_bracket.py --validate`, not CI.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
import ace_co_mu1_bracket as mu1  # noqa: E402
import ispec.synth.turbospectrum as ts  # noqa: E402


def test_remap_picks_intensity_column(tmp_path):
    # bsyn intensity RESULTFILE: wave, FLUX(col1), 2x abs, mu-intensity(col4).
    # Remap must rewrite col1 := col4 (intensity), so iSpec's data[:,1] = intensity.
    rf = tmp_path / 'result.txt'
    rows = np.array([
        [22958.16, 0.70124, 2.289e5, 8.05e4, 0.71861],   # deep line: flux<intensity
        [22958.17, 0.99999, 3.26e5, 1.12e5, 0.99999],
    ])
    header = '# mu-points   1.000000E+00\n'
    rf.write_text(header + '\n'.join('  '.join(f'{x:.5f}' for x in r) for r in rows) + '\n')
    stats = {'remapped': 0}
    mu1._remap_intensity_column(str(rf), stats)
    d = np.loadtxt(rf)
    assert d.shape[1] == 2                       # collapsed to (wave, value)
    assert stats['remapped'] == 1
    np.testing.assert_allclose(d[:, 0], rows[:, 0], rtol=0, atol=1e-3)
    np.testing.assert_allclose(d[:, 1], rows[:, 4], rtol=0, atol=1e-4)   # = intensity, not flux
    assert d[0, 1] != pytest.approx(rows[0, 1])  # NOT the flux value


def test_remap_refuses_flux_format(tmp_path):
    # A flux-shaped RESULTFILE (<5 cols) must RAISE, never silently pass flux as mu=1.
    rf = tmp_path / 'flux.txt'
    np.savetxt(rf, np.array([[22958.16, 0.70124, 2.289e5], [22958.17, 0.99999, 3.26e5]]))
    with pytest.raises(RuntimeError, match='intensity remap FAILED'):
        mu1._remap_intensity_column(str(rf), {'remapped': 0})


def test_injection_writes_mupoints_and_restores_popen():
    real = ts.subprocess.Popen
    with mu1.turbospectrum_intensity_at_mu(1.0):
        assert ts.subprocess.Popen is not real          # patched inside
        assert mu1._MU_POINTS_FILE.read_text().split() == ['1', '1.0000']
    assert ts.subprocess.Popen is real                  # restored on exit


def test_injection_swaps_flux_token(monkeypatch):
    # Drive the patched Popen with a fake process; confirm the bsyn stdin flux token
    # is flipped to Intensity + MU-POINTS before reaching the real communicate.
    seen = {}

    class FakeProc:
        def communicate(self, input=None, timeout=None):
            seen['input'] = input
            return (b'', b'')

    monkeypatch.setattr(ts.subprocess, 'Popen', lambda *a, **k: FakeProc())
    with mu1.turbospectrum_intensity_at_mu(1.0):
        patched = ts.subprocess.Popen
        p = patched(['bsyn_lu'])
        p.communicate(input=b"'INTENSITY/FLUX:' 'Flux'\n'RESULTFILE :' '/tmp/none'")
    assert b"'INTENSITY/FLUX:' 'Intensity'" in seen['input']
    assert b"'MU-POINTS:'" in seen['input']
    assert b"'Flux'" not in seen['input']


def test_constants_inherited_from_440():
    # Single-source: reference + tolerance come from the 440 harness / constants.
    assert mu1.ace.A_C_REF == 8.46
    assert mu1.A_C_FLUX_440 == 8.497
    assert mu1.ace.NEAR_REF_DEX == 0.10
