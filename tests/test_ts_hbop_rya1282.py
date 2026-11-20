"""RYA-1282: the Python port of TS hbop.f's continuum (contonly) against the compiled Fortran.

The referee is Turbospectrum's own hbop.f + hlinprof.f, compiled here with gfortran and
called as absko.f calls it (LTE, contonly=.true.). Measured: median 2.9e-6, max 9.3e-5
(single-precision Fortran), identical zero pattern, 5 plasmas x 319 wavelengths 70-380 nm.
Pinned values guard the port where no compiler/TS source is held.
"""
import shutil, subprocess
import numpy as np
import pytest
from pipeline.ts_hbop import hbop_continuum, hbf, wcalc

DRIVER = '''      program drv
      implicit none
      integer ios
      integer nlo(1), nup(1)
      real*8 hlam(1), wave, srcf
      real nh, nhe, ne, t, dop, total, cont, npop(1), bd(1)
      logical contonly, lineonly
      nlo(1)=1; nup(1)=2; hlam(1)=1215.67d0
      dop=0.; npop(1)=0.; bd(1)=1.; srcf=0.d0
 10   read(*,*,iostat=ios) wave, t, nh, nhe, ne
      if (ios .ne. 0) stop
      contonly=.true.; lineonly=.false.
      call hbop(wave,0,nlo,nup,hlam,nh,nhe,ne,t,dop,npop,bd,0,
     &  1.,1.,srcf,1.,total,cont,contonly,lineonly,.false.,.false.)
      write(*,'(1p,e23.15)') cont
      goto 10
      end
'''
WAVES = np.r_[np.linspace(700, 3800, 311), 911.0, 913.0, 1215.0, 1216.0, 1300., 3645., 3648., 3700.]
PLASMAS = [(3600., 1e16, 1e15, 1e10), (4500., 3e17, 3e16, 3e12), (6300., 1e17, 1e16, 1e14),
           (9800., 1e17, 1e16, 3e15), (5000., 1e15, 1e14, 1e9)]


def _ts_source():
    from pathlib import Path
    for p in (Path.home() / 'codex/engines/Turbospectrum_NLTE/source',):
        if (p / 'hbop.f').exists():
            return p
    return None


def test_port_matches_the_compiled_fortran(tmp_path):
    src = _ts_source(); fc = shutil.which('gfortran')
    if src is None or fc is None:
        pytest.skip('TS source or gfortran not held here (pinned-value test still runs)')
    (tmp_path / 'drv.f').write_text(DRIVER)
    subprocess.run([fc, '-O0', '-std=legacy', '-w', '-o', str(tmp_path / 'drv'), str(tmp_path / 'drv.f'),
                    str(src / 'hbop.f'), str(src / 'hlinprof.f')], check=True)
    lines = [f'{float(w):.6f} {T} {nh} {nhe} {ne}' for (T, nh, nhe, ne) in PLASMAS for w in WAVES]
    out = subprocess.run([str(tmp_path / 'drv')], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout
    ref = np.array([float(x) for x in out.split() if 'E' in x]).reshape(len(PLASMAS), len(WAVES))
    mine = np.array([hbop_continuum(WAVES, [T], [nh], [nhe], [ne])[:, 0] for (T, nh, nhe, ne) in PLASMAS])
    assert np.array_equal(ref > 0, mine > 0)
    m = ref > 0
    r = np.abs(mine[m] / ref[m] - 1)
    assert np.median(r) < 1e-5 and r.max() < 3e-4                # measured 2.9e-6 / 9.3e-5


def test_pinned_values_and_structure():
    # values from the compiled Fortran (gfortran 16.1, TS v20 hbop.f)
    v = hbop_continuum([1500.0], [4500.], [3e17], [3e16], [3e12])[0, 0]
    assert v == pytest.approx(1.333980617346242e-05, rel=3e-4)
    # n=1 cross-section stops at the Balmer edge table limit (BPz 2025-08-25 fix)
    assert hbf(3646.9, 1) > 0 and hbf(3647.0, 1) == 0
    # occupation probability: 1 for an isolated atom, falls with density
    assert wcalc(0., 0., 0., 2., 5000.) == pytest.approx(1.)
    assert wcalc(1e18, 1e15, 1e17, 10., 5000.) < wcalc(1e15, 1e12, 1e14, 10., 5000.)
    # the dissolved Lyman continuum is present redward of Ly-alpha and tapers off
    a, b = hbop_continuum([1300., 1800.], [4500.], [3e17], [3e16], [3e12])[:, 0]
    assert a > b > 0
