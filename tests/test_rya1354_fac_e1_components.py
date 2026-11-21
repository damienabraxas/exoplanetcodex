from pathlib import Path
import math
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_e1_components_audit import parse,verify,decompose,grid_context,configuration_shift,gamma,U
from rya1354_fac_e1_components import input_text,FORBIDDEN,UNSTABLE,STATE_FILES
from rya1354_fac_radial_rates import rate_input

FIX=ROOT/'tests/fixtures/rya1354_e1_components_smoke'


def mapping():
    """Synthetic level inventory mirroring the sealed v41 J/parity layout (not its energies)."""
    rows=[(0,.0,4,'+'),(1,.02,2,'+'),(2,.03,0,'+'),(3,.08,4,'+'),(4,.2,0,'+'),(5,.36,4,'-'),(6,.38,2,'-')]
    return dict(targets={'3S':[6]},levels={str(i):dict(energy_au=e,two_J=j,parity=p) for i,e,j,p in rows})


def fixture(g):
    header,contexts=parse((FIX/f'components_gauge{g}.tsv').read_text())
    raw={}
    for line in (FIX/f'raw_gauge{g}.tsv').read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split();raw[(int(f[1]),int(f[0]))]=float(f[12])
    return header,contexts,raw


def test_component_calls_sit_inside_the_native_e1_photon_grid_context():
    for scale in (10,100):
        for gauge in (1,2):
            text,requests,e1=input_text(scale,gauge,mapping())
            original=rate_input(scale,gauge,mapping())[0]
            assert e1==[0,1,2,3,4] and text.count('RawE1Components(')==5
            assert text.replace(''.join(f"RawE1Components('components.tsv',{l},6,{gauge})\n" for l in e1),'')==original
            calls=text.index('RawE1Components(');last=text.index("RawTransition('raw.tsv',4,6,-1)")
            nxt=text.index('TransitionTable(',last)
            assert last<calls<nxt and text.count('TransitionTable(',0,calls)==1
            assert not any(x in text for x in FORBIDDEN)


def test_measured_null_keeps_state_files_under_byte_identity():
    assert not set(STATE_FILES)&UNSTABLE
    assert {'sfac.stdout','sfac.stderr','o.basis'}.isdisjoint(UNSTABLE)


def test_smoke_fixture_reproduces_native_amplitudes_in_both_branches():
    for g in (1,2):
        header,contexts,raw=fixture(g)
        assert float(header['ANGZCUT'])==1e-8 and [c['swapped'] for c in contexts]==[1,0]
        for c in contexts:
            v=verify(c,raw[(c['lower'],c['upper'])],float(header['ANGZCUT']))
            assert v['native_amplitude']==raw[(c['lower'],c['upper'])]
            assert abs(v['reconstruction_difference'])<=v['reconstruction_bound']
            d=decompose(c)
            assert math.isclose(math.fsum(r['amplitude'] for r in d['configuration_pairs']),d['amplitude'],rel_tol=0,abs_tol=4*U*d['abs_sum'])
            assert d['cancellation_index']>=1
            assert grid_context(c)['inside']
        assert contexts[0]['B'][6]['flags']==8 and decompose(contexts[0])['cut']['terms']==1


def edited(g,tag,index,column,value):
    lines=(FIX/f'components_gauge{g}.tsv').read_text().splitlines()
    hits=[i for i,l in enumerate(lines) if l.startswith(tag+'\t')]
    f=lines[hits[index]].split('\t');f[column]=value;lines[hits[index]]='\t'.join(f)
    return parse('\n'.join(lines))


@pytest.mark.parametrize('tag,index,column,value,match',[
    ('A',0,9,'1.0','running sum'),
    ('B',0,10,'0.5','angular x mixing'),
    ('B',0,11,'1','cut flags'),
    ('B',0,12,'1.0','product'),
    ('L',0,5,'-0.8','norm'),
])
def test_any_altered_record_is_rejected(tag,index,column,value,match):
    header,contexts=edited(2,tag,index,column,value)
    raw=fixture(2)[2]
    with pytest.raises(ValueError,match=match):
        verify(contexts[0],raw[(0,6)],float(header['ANGZCUT']))


def test_native_amplitude_must_match_exactly_not_approximately():
    header,contexts,raw=fixture(2);c=contexts[1]
    with pytest.raises(ValueError,match='native amplitude'):
        verify(c,math.nextafter(raw[(1,6)],0),float(header['ANGZCUT']))


def test_structure_errors_are_not_silently_parsed():
    text=(FIX/'components_gauge1.tsv').read_text()
    with pytest.raises(ValueError,match='header'):parse(text.replace('FAC_E1_COMPONENTS_V44','X'))
    with pytest.raises(ValueError,match='Unterminated|Missing'):parse(text.rsplit('\nE\t',1)[0])
    with pytest.raises(ValueError,match='outside'):parse(text+'Z\t1\n')
    with pytest.raises(ValueError,match='Unknown'):parse(text.replace('\nE\t0\t','\nZ\t1\nE\t0\t',1))


def test_summation_bound_and_configuration_shift_accounting():
    assert gamma(1)>U and gamma(10)<11*U
    with pytest.raises(ValueError):gamma(2**53)
    _,c1,_=fixture(1);_,c2,_=fixture(2)
    s=configuration_shift(decompose(c1[0]),decompose(c2[0]))
    assert math.isclose(s['accounted'],s['total_change'],rel_tol=1e-12,abs_tol=1e-15)


def test_actual_v44_archive_replays_every_component_and_the_machine_report():
    import json
    from rya1354_fac_e1_components_report import summarize
    base=ROOT/'data/audit/rya1354'
    result=summarize(base,ROOT)
    assert result==json.loads((base/'source_curation_v2/fac_e1_components_v44.json').read_text())
    assert result['archive']['native_paths']==882 and result['archive']['referenced_files']==744
    assert result['scientific_admission']=='HOLD' and not result['adoption_allowed'] and result['accepted_lifetime'] is None
    assert len(result['cases'])==8 and all(len(c['channels'])==5 for c in result['cases'].values())
    for name,case in result['cases'].items():
        for ch in case['channels'].values():
            v=ch['verification']
            assert v['native_amplitude']==ch['decomposition']['amplitude'] or abs(v['reconstruction_difference'])<=v['reconstruction_bound']
            assert abs(ch['decomposition']['cut']['amplitude'])<1e-5*abs(v['native_amplitude'])+1e-6
        # Reference 3P2: native photon grid equals exact own-energy radials to rounding.
        assert case['reference_radial_vs_v43']['max_abs_relative']<1e-15
    # Mechanism: the lower-level 2p3 3p admixture via the 3s-3p radial dominates the gauge split
    # in the equal-reference model and is minor in the ground model.
    def adm(name):
        d=result['cases'][name]['channels']['0']['decomposition']
        return {(r['key']['lower_config'],r['key']['upper_config']):r['amplitude'] for r in d['configuration_pairs']}[('2p3.3p1','2p3.3s1')]
    assert adm('equal_coherent100_B')<-0.9 and adm('equal_coherent100_C')>0.15
    assert abs(adm('ground_coherent100_B'))<0.08 and abs(adm('ground_coherent100_C'))<0.011


def test_mixcut_control_inserts_only_setmixcut_before_structure():
    from rya1354_fac_e1_components import MIXCUT_CHANGED,level_table
    for gauge in (1,2):
        v44=input_text(100,gauge,mapping())[0];v45=input_text(100,gauge,mapping(),0)[0]
        assert v45.replace('SetMixCut(0)\n','',1)==v44 and v45.count('SetMixCut(')==1
        assert v45.index('SetMixCut(0)\n')+len('SetMixCut(0)\n')==v45.index("Structure('o.lev.b', ")
    with pytest.raises(ValueError,match='declared full-eigenvector'):input_text(100,2,mapping(),1e-6)
    assert MIXCUT_CHANGED=={'o.basis','o_export_000.cm','raw.tsv'} and not MIXCUT_CHANGED&{'ci.pot','source.pot','o_export_000.c'}
    assert level_table('A\nTSess\t= 1\nB')==level_table('A\nTSess\t= 2\nB')!=level_table('A\nTSess\t= 2\nC')


def test_mixcut_review_view_passes_exact_v44_input_and_retains_strictness(tmp_path,monkeypatch):
    import rya1354_fac_e1_components as m
    v44=input_text(100,2,mapping())[0];v45=input_text(100,2,mapping(),0)[0]
    out=tmp_path/'out';out.mkdir();(out/'input.sf').write_text(v45);(out/'o.basis').write_text('basis')
    seen={}
    def fake(view,layout,scale,declared_input=None):
        seen.update(input=(view/'input.sf').read_text(),declared=declared_input,basis=(view/'o.basis').read_text(),
                    link=(view/'o.basis').is_symlink())
        return {'ok':True}
    monkeypatch.setattr(m,'ci_review',fake)
    assert m.mixcut_ci_review(out,None,100,v45,0)=={'ok':True}
    assert seen==dict(input=v44,declared=v44,basis='basis',link=True) and (out/'input.sf').read_text()==v45
    with pytest.raises(ValueError,match='single SetMixCut'):m.mixcut_ci_review(out,None,100,v44+'SetMixCut(0)\n',0)
    with pytest.raises(ValueError,match='single SetMixCut'):m.mixcut_ci_review(out,None,100,v45.replace('SetMixCut(0)\n','SetMixCut(0)\nSetMixCut(0)\n'),0)


def test_actual_v46_full_eigenvector_control_replays_and_bounds_the_mixing_cut():
    import json
    from rya1354_fac_e1_components_report import summarize_mixcut
    base=ROOT/'data/audit/rya1354'
    result=summarize_mixcut(base,ROOT)
    assert result==json.loads((base/'source_curation_v2/fac_mixcut_control_v46.json').read_text())
    assert result['retained_failure']['status']=='FAILED_RETAINED' and result['retained_failure']['canary_files']==56
    assert result['retained_failure']['completed_native_case_reproduced']==['ground_coherent100_B']
    assert result['scientific_admission']=='HOLD' and not result['adoption_allowed'] and result['accepted_lifetime'] is None
    assert set(result['cases'])=={f'{m}_coherent100_{g}' for m in ('ground','equal') for g in 'BC'}
    for case in result['cases'].values():
        assert case['byte_identical_v41_files']==53 and abs(case['total_rate_s1']['ratio']-1)<1e-5
        for low,ch in case['channels'].items():
            full=ch['stored_basis']['full_mixing']
            assert all(full[s]['dropped_by_mixing_cut']==0 for s in ('0','1'))
            if low in ('0','1','2'):  # 3P2,1,0: the stored-mixing truncation is not the gauge split
                assert abs(ch['relative_amplitude_change'])<1e-5
