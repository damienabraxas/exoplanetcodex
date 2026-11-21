import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import rya1354_mcdhf_s0_invariance as I  # noqa: E402
import rya1354_mcdhf_eol as M  # noqa: E402

TABLE=""" Electric 2**( 1)-pole transitions
 =================================

 Upper       Lower 
 Lev  J P   Lev  J P       E (Kays)         A (s-1)          gf            S
 f2  1    1 -  f1  1    2 +       76794.98 C  5.12345D+08  3.00000D-01  1.20000D+00
                                           B  5.00000D+08  2.90000D-01  1.10000D+00
 f1  3    1 +  f2  2    2 -       10000.00 C  1.00000D+06  1.00000D-02  2.00000D-07
                                           B  9.00000D+05  9.00000D-03  1.00000D-07
"""


def test_line_strength_parser_reads_both_gauges():
    s=I.line_strengths(TABLE)
    assert s[('f2',1,1,'f1',1,2)]=={'E':76794.98,'S_C':1.2,'S_B':1.1}
    assert len(s)==2


def test_compare_tolerances_and_floor():
    s=I.line_strengths(TABLE)
    t={k:dict(v) for k,v in s.items()}
    t[('f1',3,1,'f2',2,2)]['S_B']*=2     # below the 1e-6 au floor: ignored
    e={'even:1:0':-74.8,'odd:1:0':-74.3}
    ok=I.compare(e,dict(e),s,t)
    assert ok['agree'] and ok['max_rel_dS']['S_B_pairs']==1
    t[('f2',1,1,'f1',1,2)]['S_C']*=1.002
    assert not I.compare(e,dict(e),s,t)['agree']
    e2=dict(e);e2['odd:1:0']+=1e-6
    assert not I.compare(e,e2,s,s)['agree'] and I.compare(e,{**e,'odd:1:0':e['odd:1:0']+4e-7},s,s)['agree']
    with pytest.raises(ValueError,match='pair sets'):I.compare(e,e,s,{})


def test_predicted_invariant_sets_are_the_complete_mr_manifolds():
    """Each predicted set is {1s,2s} or the full set of MR orbitals of one l above n=2."""
    for p in ('even','odd'):
        top=M.mr_max(p)
        manifolds=[{f'{n}{l}' for n in range({'s':3,'p':3,'d':3,'f':4}[l],top[l]+1)} for l in top if top[l]>2]
        assert sorted(map(sorted,I.INVARIANT_SETS[p]))==sorted(map(sorted,[{'1s','2s'}]+manifolds))
    assert I.in_invariant_set('6p-','even') and not I.in_invariant_set('2p-','even') and I.in_invariant_set('4s','odd')


def test_v56_strict_stalls_lie_in_predicted_sets():
    """The orbitals that stalled in v56 (strict) are all predicted-redundant."""
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_controls_v56',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_controls_v56/receipt.json'])
    for name in ('even_method3_all','odd_method3_all'):
        c=r['runs'][name];acc=c['scf']['native_function_ACC']
        stalled=[o for o,v in c['scf']['reported_orbital_self_consistency'].items() if v>acc]
        assert stalled and all(I.in_invariant_set(o,c['parity']) for o in stalled),(name,stalled)


def test_energies_parser():
    lbl=' Pos   J   Parity      Energy Total      Comp. of ASF\n  1    0     +           -74.848089984     100.000%\n   0.9 0.8 x\n  1    1     +           -74.8     100.000%\n'
    assert I.energies(lbl)=={(0,0):-74.848089984,(1,1):-74.8}
