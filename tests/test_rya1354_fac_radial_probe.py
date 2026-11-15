from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_radial_probe import input_text,wave,overlaps,OLD,NEW


def test_radial_probe_has_no_new_CI_optimization_or_rate():
    t=input_text(100)
    assert 'RestorePotential' in t and 'OptimizeRadial' not in t and 'Structure(' not in t and 'TransitionTable' not in t
    assert t.count('WaveFuncTable(')==2*len(OLD)+len(NEW)
    assert 'SetConfigEnergyMode(-1)' in t and t.count('RadialOverlaps(')==10


def test_wave_headers_samples_and_bound_identity_guard():
    t='# n = 4\n# kappa = -1\n# energy = -1.00000000E+00\n# ilast = 1\n# qrnorm = 1.0\n0 0.1 1 2 3 4 0 0 0 0 0 0\n1 0.2 1 2 2 3 0 0 0 0 0 0\n'
    r,values=wave(t,4,-1);assert r['samples']==2 and values[0]==(0,.1,1.,3.,4.)
    for bad in (t.replace('energy = -1','energy = 1'),t.replace('ilast = 1','ilast = 2'),t.replace('1 0.2','1 0.05')):
        with pytest.raises(ValueError):wave(bad,4,-1)


def test_overlap_unit_norm_and_pair_completeness_guard():
    t='4 -1 -1 -1 4 -1 -1 -1 1.00000E+00 -1 -1\n'
    assert overlaps(t,-1,((4,-1),))['pairs']==1
    with pytest.raises(ValueError):overlaps(t.replace('1.00000E+00','9.99999E-01'),-1,((4,-1),))
    with pytest.raises(ValueError):overlaps(t,-1,((3,-1),(4,-1)))
