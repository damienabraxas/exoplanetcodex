import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import rya1354_mcdhf_eol as M  # noqa: E402
from rya1354_li_scheme_sizing import stdin_for, TWOJ as V47_TWOJ  # noqa: E402


def test_targets_are_li_45_even_36_odd_lowest_states():
    assert sum(M.roots('even').values())==45 and sum(M.roots('odd').values())==36
    assert M.roots('even')=={0:6,1:11,2:14,3:8,4:4,5:2}
    assert M.roots('odd')=={0:3,1:11,2:12,3:7,4:3}
    # every target J lies inside the generated 2J range
    for p in ('even','odd'):
        lo,hi=map(int,M.TWOJ[p].split(','));assert all(lo<=2*j<=hi for j in M.roots(p))


def test_layer_lists_extend_v47_mcdhf_definition():
    """n6 differs from the v47 sized list only in the wider 2J range (needed for the J<=5 targets)."""
    for p in ('even','odd'):
        ours=M.generator_stdin(p,6).split('\n');v47=stdin_for('mcdhf',p,6).split('\n')
        diff=[i for i,(x,y) in enumerate(zip(ours,v47)) if x!=y]
        assert len(ours)==len(v47) and [ours[i] for i in diff]==[M.TWOJ[p]] and [v47[i] for i in diff]==[V47_TWOJ[p]]
    assert M.generator_stdin('even',0).split('\n')[-3]=='0'


def test_active_sets_keep_mr_orbitals_and_vary_only_new_shell():
    assert M.mr_max('even')=={'s':2,'p':6,'f':5} and M.mr_max('odd')=={'s':6,'p':2,'d':5}
    assert M.active_set('even',3)=='3s,6p,3d,5f' and M.active_set('odd',4)=='6s,4p,5d,4f'
    for p in ('even','odd'):
        assert M.active_set(p,6)=='6s,6p,6d,6f,6g,6h'
        seen=set()
        for n in M.LAYERS:
            new=M.new_orbitals(p,n)
            assert all(o.endswith('*') for o in new)   # GETRSL: '3d' alone is only j=5/2
            names={o[:-1] for o in new}
            assert not names&seen and all(int(o[:-1][:-1])==n for o in new);seen|=names
            assert not any(M.mr_max(p).get(o[-2],0)>=n for o in new)


def test_levels_stdin_follows_native_block_order_and_refuses_mismatch():
    s=M.levels_stdin([0,1,2,3,4],'odd')
    assert s.split('\n')[:5]==['1,2,3','1,2,3,4,5,6,7,8,9,10,11','1,2,3,4,5,6,7,8,9,10,11,12','1,2,3,4,5,6,7','1,2,3']
    assert M.levels_stdin([4,3,2,1,0],'odd').split('\n')[0]=='1,2,3'
    with pytest.raises(ValueError,match='differ'):M.levels_stdin([0,1,2,3],'odd')
    with pytest.raises(ValueError,match='differ'):M.levels_stdin([0,1,2,3,4,5],'odd')


def test_block_order_reads_native_list(tmp_path):
    p=tmp_path/'rcsf.out'
    rec=lambda j:f'  2p-( 2)  2p ( 1)  3s ( 1)\n               3/2      1/2\n                           {j}-\n'
    p.write_text('Core subshells:\n  1s\nPeel subshells:\n  2s   2p-  2p   3s\nCSF(s):\n'+rec(2)+rec(2)+' *\n'+rec(1)+' *\n'+rec(0))
    assert M.block_order(p)==[2,1,0]


def lbl(levels):
    out=[' Pos   J   Parity      Energy Total      Comp. of ASF']
    for i,(j,par,label) in enumerate(levels,1):
        out+=[f'  {i}    {j}     {par}           -74.8{i:08d}     100.000%',f'         0.99    0.98   {label}',
              '         0.05    0.003   1s(2).2s(2).2p(3)2D.3p_3P']
    return '\n'.join(out)+'\n'


def synthetic(parity):
    lab={'ground':'1s(2).2s(2).2p(4)_{t}','4S.p':'1s(2).2s(2).2p(3)4S.{n}p_{t}','4S.f':'1s(2).2s(2).2p(3)4S.{n}f_{t}',
         '4S.s':'1s(2).2s(2).2p(3)4S.{n}s_{t}','4S.d':'1s(2).2s(2).2p(3)4S.{n}d_{t}','2D.s':'1s(2).2s(2).2p(3)2D.{n}s_{t}'}
    levels=[]
    for fam,n in M.TARGET_FAMILIES[parity]:
        for term,js in M.TERMS[fam]:
            for j in js:levels.append((j,'+' if parity=='even' else '-',lab[fam].format(n=n,t=term)))
    return levels


@pytest.mark.parametrize('parity',['even','odd'])
def test_identity_accepts_targets_and_refuses_intruder(parity):
    levels=synthetic(parity)
    got=M.identity(lbl(levels),parity)
    assert sum(sum(c.values()) for c in got.values())==(45 if parity=='even' else 36)
    j,par,_=levels[-1]
    with pytest.raises(ValueError,match='identity'):
        M.identity(lbl(levels[:-1]+[(j,par,'1s(2).2s(2).2p(3)2P.3p_3D' if parity=='even' else '1s(2).2s(2).2p(3)2P.3s_3P')]),parity)
    with pytest.raises(ValueError,match='parity'):
        M.identity(lbl([(j,'-' if par=='+' else '+',x) for j,par,x in levels]),parity)


def test_family_labels():
    assert M.family('1s(2).2s(2).2p(4)_3P')==('ground',2)
    assert M.family('1s(2).2s(2).2p(3)4S.5f_5F')==('4S.f',5)
    assert M.family('1s(2).2s(2).2p(3)2D.3s_1D')==('2D.s',3)
    assert M.family('1s(2).2s(2).2p(2)3P.3s(2)_3P')==('other',0)


def test_rmcdhf_stdin_order_matches_getold():
    """GETOLD: default? -> level lists (LODSTATE) -> weights (GETOLDWT, 5 standard) -> varied -> spectroscopic -> iterations."""
    s=M.rmcdhf_stdin('1,2\n1\n','3s*,3d*','')
    assert s=='y\n1,2\n1\n5\n3s*,3d*\n\n200\n'


def test_v53_archive_records_the_5f_stall():
    """v53 FAILED_RETAINED: strict iteration limit with 5f (j=7/2) frozen at 1.66e-1."""
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_eol_failed_v53',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_eol_v53/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'MAXIMUM ITERATIONS' in r['error'] and list(r['stages'])==['s0_even']
    assert r['stages']['s0_even']['blocks']=={'0+':20,'1+':47,'2+':56,'3+':40,'4+':22,'5+':10}
    out=files['runs/mcdhf_eol_v53/s0_even/rmcdhf.stdout'].decode()
    last=out.rsplit('Iteration number',1)[1]
    row=[l.split() for l in last.splitlines() if l.split()[:1]==['5f']][0]
    assert float(row[4].replace('D','E'))>0.1 and 'could not be solved using method 1' in out
