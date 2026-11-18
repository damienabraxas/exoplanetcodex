import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import rya1354_mcdhf_s0_sequential as S  # noqa: E402
import rya1354_mcdhf_eol as M  # noqa: E402
import rya1354_mcdhf_s0_controls as C  # noqa: E402


def test_steps_cover_li_targets_and_every_mr_orbital_once():
    for p in ('even','odd'):
        assert S.STEPS[p][-1][0]=='ci' and S.STEPS[p][-1][2]==[]
        assert sorted(S.included(p,len(S.STEPS[p])-1))==sorted(M.TARGET_FAMILIES[p])
        assert S.expected(S.included(p,len(S.STEPS[p])-1))==M.expected_targets(p)
        varied=[o for _,_,v in S.STEPS[p] for o in v]
        assert len(varied)==len(set(varied))
        mr=set();kept,_=__import__('rya1354_li_scheme_sizing').references('mcdhf',p,6)
        import re
        for line in kept:mr|={n+l for n,l in re.findall(r'(\d+)([spdfgh])\(',line)}
        assert set(varied)==mr


def test_final_step_generator_equals_v53_s0_list_definition():
    for p in ('even','odd'):
        k=len(S.STEPS[p])-1
        ours=S.generator_stdin(p,k).split('\n');v53=M.generator_stdin(p,0).split('\n')
        assert ours[:2]==v53[:2] and sorted(ours[2:-5])==sorted(v53[2:-5]) and ours[-4]==v53[-4]==M.TWOJ[p]
        assert S.active([l for l in ours[2:-5] if l])==M.active_set(p,0)


def test_step_ranges_and_lists():
    assert S.generator_stdin('even',0)=='*\n0\n1s(2,i)2s(2,i)2p(4,*)\n\n2s,2p\n0,4\n0\nn\n'
    assert S.twoj_range(S.included('odd',0))=='2,6'
    assert S.twoj_range(S.included('even',1))=='0,6'
    assert S.config_lines('odd',S.included('odd',3))==['1s(2,i)2s(2,i)2p(3,*)3s(1,*)','1s(2,i)2s(2,i)2p(3,*)4s(1,*)',
        '1s(2,i)2s(2,i)2p(3,*)5s(1,*)','1s(2,i)2s(2,i)2p(3,*)6s(1,*)']
    assert S.config_lines('odd',S.included('odd',5))==['1s(2,i)2s(2,i)2p(3,*)3d(1,*)','1s(2,i)2s(2,i)2p(3,*)4d(1,*)']


def test_rydberg_steps_use_their_own_series_without_ground_states():
    """v59 FAILED_RETAINED: 4f collapsed with the 2p4 ground states in its EOL."""
    for p in ('even','odd'):
        for label,fams,varied in S.STEPS[p][1:-1]:
            l=varied[0][-1]
            assert ('ground',2) not in fams and all(f[0].endswith('.'+l) for f in fams if f!=('2D.s',3))
            ns=sorted(n for f,n in fams if f.startswith('4S'))
            assert ns==list(range(ns[0],int(varied[0][:-1])+1))   # every lower orbital of that l stays in the list


def test_twoj_range_refuses_a_gap(monkeypatch):
    from collections import Counter
    monkeypatch.setattr(S,'expected',lambda fams:{1:Counter(),3:Counter()})
    with pytest.raises(ValueError,match='contiguous'):S.twoj_range([])


def test_rmcdhf_stdin_is_v56_path_with_ladder_half_grid_and_named_varied_orbitals():
    import math
    s=S.rmcdhf_stdin('1,2\n1\n',['3p'])
    ref=C.rmcdhf_stdin('method3_all','even','1,2\n1\n')
    # lifetime ladder half-grid prompt answers (ran natively in lifetime_ladder_v3)
    ladder_grid=f'n\nn\ny\nn\ny\n{math.exp(-65/16)/8:.17g}\n0.03125\n0\n439\ny\n{.0625**6:.17g}\n'
    assert s==ref.replace('n\nn\nn\nn\n',ladder_grid,1).replace('5\n*\n*\n','5\n3p*\n3p*\n').replace('\n\n\n*\n\n','\n\n\n3p*\n\n')


def test_rwfnestimate_half_grid_matches_ladder_initialization():
    from rya1354_radial_controls import initialization_input
    lad=initialization_input(3,439,.03125)           # previous.w, then method 3, then n
    assert S.STEP_ESTIMATE==lad.replace('*\n3\n*\nn\n','*\n2\n*\n3\n*\nn\n')
    assert S.FIRST_ESTIMATE==lad.split('1\nprevious.w')[0]+'2\n*\n3\n*\nn\n'


def test_identity_for_partial_families():
    fams=S.included('odd',0)
    lab={'4S.s':'1s(2).2s(2).2p(3)4S.3s_{t}','2D.s':'1s(2).2s(2).2p(3)2D.3s_{t}'}
    out=[' Pos   J   Parity      Energy Total      Comp. of ASF'];i=0
    for fam,n in fams:
        for t,js in M.TERMS[fam]:
            for j in js:
                i+=1;out+=[f'  {i}    {j}     -           -74.8{i:08d}     100.000%',f'      0.99  0.98  {lab[fam].format(t=t)}']
    got=S.identity('\n'.join(out)+'\n','odd',fams)
    assert sum(sum(v.values()) for v in got.values())==6
    with pytest.raises(ValueError,match='identity'):S.identity('\n'.join(out[:-2])+'\n','odd',fams)


def test_v57_archive_records_the_divergent_start():
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_invariance_failed_v57',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_invariance_v57/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'even_B/rmcdhf' in r['error'] and r['verdict'] is None
    assert r['runs']['even_A']['outcome']=='CONVERGED_IDENTITY_DIFFERS'
    assert 'continuum state' in files['runs/mcdhf_s0_invariance_v57/even_B/rmcdhf.stdout'].decode()


def test_v58_archive_records_the_1s_grid_floor():
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_sequential_failed_v58',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_sequential_v58/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'even_0_core/rmcdhf' in r['error']
    last=files['runs/mcdhf_s0_sequential_v58/even_0_core/rmcdhf.stdout'].decode().rsplit('Iteration number',1)[1]
    res={l.split()[0]:float(l.split()[4].replace('D','E')) for l in last.splitlines() if l.split()[:1] in (['1s'],['2s'],['2p-'],['2p'])}
    assert res['1s']>.0625**6 and all(v<.0625**6 for k,v in res.items() if k!='1s')


def test_step_estimate_reads_every_earlier_step_then_tf_then_hydrogenic():
    e=S.step_estimate(3)
    assert e.startswith(S.FIRST_ESTIMATE.split('2\n*\n3')[0])
    assert e.endswith('1\nprev0.w\n*\n1\nprev1.w\n*\n1\nprev2.w\n*\n2\n*\n3\n*\nn\n')


def test_ci_step_stdin_varies_nothing():
    s=S.rmcdhf_stdin('1\n',[]).split('\n')
    i=s.index('5');assert s[i+1:i+3]==['',''] and s[i+3]=='500'


def test_v59_archive_records_core_and_np_success_and_4f_collapse():
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_sequential_failed_v59',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_sequential_v59/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'even_5_4f/rmcdhf' in r['error']
    ok=[k for k,c in r['steps'].items() if c.get('leading_families')]
    assert ok==['even_0_core','even_1_3p','even_2_4p','even_3_5p','even_4_6p']
    assert all(r['steps'][k]['scf']['all_reported_functions_below_ACC'] for k in ok)
    last=files['runs/mcdhf_s0_sequential_v59/even_5_4f/rmcdhf.stdout'].decode().rsplit('Iteration number',1)[1]
    e={l.split()[0]:float(l.split()[1].replace('D','E')) for l in last.splitlines() if l.split()[:1] in (['4f'],['4f-'])}
    assert e['4f']>0.5 and e['4f-']<0.05     # j=7/2 collapsed, j=5/2 Rydberg-like


def test_v60_archive_f72_has_no_weight():
    """v60 FAILED_RETAINED: np series qualified again; in the 4f-only step the targeted states carry ~no f7/2."""
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_sequential_failed_v60',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_sequential_v60/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'even_5_4f/rmcdhf' in r['error']
    assert [k for k,c in r['steps'].items() if c.get('leading_families')]==['even_0_core','even_1_3p','even_2_4p','even_3_5p','even_4_6p']
    out=files['runs/mcdhf_s0_sequential_v60/even_5_4f/rmcdhf.stdout'].decode()
    occ=[float(x.replace('D','E')) for x in out.rsplit('Generalised occupation numbers:',1)[1].split()[:6]]
    assert occ[4]>0.99 and occ[5]<1e-6       # 4f- carries the electron, 4f (j=7/2) none
    assert 'continuum state' in out
