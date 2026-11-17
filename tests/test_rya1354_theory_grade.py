import json
import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_theory_grade import STATUS, accuracy_class, cancellation_factor, dT, dT_tilde, group_of  # noqa: E402


@pytest.mark.parametrize('d,cf,cls',[(0.03,0.1,'A'),(0.0300001,0.1,'B'),(0.03,0.0999,'B'),(0.10,0.5,'B'),
    (0.10,0.05,'C'),(0.25,1,'C'),(0.2500001,1,'D'),(0.25,0.01,'D'),(0.50,1,'D'),(0.50,0.01,'E'),(0.51,1,'E'),(0.9,0,'E')])
def test_li2023_table5_boundaries(d,cf,cls):
    assert accuracy_class(d,cf)==cls


def test_dT_and_groups():
    assert dT(4,5)==pytest.approx(0.2) and dT(5,4)==dT(4,5) and dT(0,0)==0
    with pytest.raises(ValueError):dT(-1,1)
    assert [group_of(a) for a in (1e-3,1e-2,5,1e3,1e6,3.5e6)]==[0,1,1,3,4,5]
    rows=dT_tilde([dict(A_B=1e8,A_C=1.0e8),dict(A_B=2e8,A_C=1.6e8),dict(A_B=10,A_C=1)])
    assert rows[0]['dT_tilde']==pytest.approx(0.1) and rows[1]['dT_tilde']==pytest.approx(0.2)
    assert rows[2]['dT_tilde']==pytest.approx(0.9) and rows[2]['group_size']==1
    with pytest.raises(ValueError):accuracy_class(float('nan'),1)
    assert STATUS.startswith('CANDIDATE')


def test_cf_sums_per_csf_pair_before_magnitudes():
    ctx=dict(B=[dict(lower=0,upper=0,term=1.0,flags=0),dict(lower=0,upper=0,term=-0.5,flags=0),
                dict(lower=1,upper=0,term=-0.25,flags=0),dict(lower=1,upper=1,term=9.0,flags=1)])
    # pairs: (0,0)=0.5, (1,0)=-0.25 -> CF=(0.25/0.75)^2; the cut term is excluded
    assert cancellation_factor(ctx)==pytest.approx((0.25/0.75)**2)


def test_actual_v46_fac_channels_cf_bounds_cancellation_index():
    from rya1354_fac_e1_components_audit import decompose, parse
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_mixcut_control_complete_v46',repo_root=ROOT)
    _,ctxs=parse(files['runs/fac_mixcut_control_v46/ground_coherent100_B/components.tsv'].decode())
    seen=0
    for ctx in ctxs:
        cf=cancellation_factor(ctx);ci=decompose(ctx)['cancellation_index']
        assert 0<=cf<=1 and cf>=1/ci**2*(1-1e-12)
        seen+=1
    assert seen==5
