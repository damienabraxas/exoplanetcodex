import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts  # noqa: E402
from rya1354_li_scheme_sizing import LAYERS, inventory_stream, stdin_for  # noqa: E402
import rya1354_rci_feasibility as F  # noqa: E402


@pytest.fixture(scope='module')
def native():
    files,_=load_parts(ROOT/'data/audit/rya1354/mr_inventory_v5')
    return files


def test_generator_reproduces_v47_stdin_and_mr_differs_only_in_excitations():
    for n in LAYERS:
        for stage in ('mcdhf','rci'):
            for parity in ('even','odd'):
                assert F.generator_stdin(stage,parity,n,2)==stdin_for(stage,parity,n)
                assert F.generator_stdin(stage,parity,n,0).splitlines()[:-2]==stdin_for(stage,parity,n).splitlines()[:-2]
                assert F.generator_stdin(stage,parity,n,0).splitlines()[-2]=='0'
    assert F.generator_stdin('rci','odd',4,2,F.COST_TWOJ).splitlines()[-3]=='2,2'


def test_interact_table_parses_native_format():
    rows=''.join('   %2d      %7d   %10d   %10d\n'%r for r in ((1,7,142302,9001),(2,9,205391,12345)))
    out=' Block    MR NCSF   Before NCSF   After NCSF\n'+rows+' RCSFinteract: Execution complete.\n'
    assert F.interact_table(out)==[dict(block=1,MR=7,before=142302,after=9001),dict(block=2,MR=9,before=205391,after=12345)]
    with pytest.raises(ValueError,match='absent'):F.interact_table('no table')


def _split(text):
    head=text.splitlines()[:5];body=text.splitlines()[5:];blocks=[[]]
    for line in body:
        if line.strip()=='*':blocks.append([])
        else:blocks[-1].append(line)
    return head,[[b[i:i+3] for i in range(0,len(b),3)] for b in blocks]


def _join(head,blocks,trim=False):
    out=list(head)
    for i,b in enumerate(blocks):
        if i:out.append(' *')
        for r in b:out+=[(x.rstrip() if trim else x) for x in r]
    return '\n'.join(out)+'\n'


def test_check_reduction_accepts_consistent_lists_and_refuses_each_defect(native,tmp_path):
    text=native[next(p for p in native if p.endswith('odd_CV_MR/rcsf.out'))].decode()
    head,blocks=_split(text)
    assert all(len(b)>=2 for b in blocks)
    full,red,mr=tmp_path/'full',tmp_path/'red',tmp_path/'mr'
    full.write_text(_join(head,blocks))
    red.write_text(_join(head,[b[:-1] for b in blocks],trim=True))      # rcsfinteract TRIMs its lines
    mr.write_text(_join(head,[b[:1] for b in blocks]))
    inv=lambda p:inventory_stream(p,'odd','rci')
    fi,ri,mi=inv(full),inv(red),inv(mr)
    labels=sorted(fi['blocks'],key=lambda s:(s[-1],int(s[:-1])))
    table=[dict(block=i+1,MR=mi['blocks'][k],before=fi['blocks'][k],after=ri['blocks'][k]) for i,k in enumerate(labels)]
    assert F.check_reduction(table,fi,ri,mi,full,red,mr)==labels
    bad=[dict(t,after=t['after']+1) if i==0 else t for i,t in enumerate(table)]
    with pytest.raises(ValueError,match='After'):F.check_reduction(bad,fi,ri,mi,full,red,mr)
    bad=[dict(t,before=t['before']-1) if i==0 else t for i,t in enumerate(table)]
    with pytest.raises(ValueError,match='Before'):F.check_reduction(bad,fi,ri,mi,full,red,mr)
    # a reduced list missing an MR CSF
    red2=tmp_path/'red2';red2.write_text(_join(head,[b[1:] for b in blocks]))
    ri2=inv(red2);t2=[dict(t,after=ri2['blocks'][k]) for t,k in zip(table,labels)]
    with pytest.raises(ValueError,match='omits an MR'):F.check_reduction(t2,fi,ri2,mi,full,red2,mr)
    # a reduced CSF absent from the full list
    full3=tmp_path/'full3';full3.write_text(_join(head,[b[1:] for b in blocks]))
    fi3=inv(full3);t3=[dict(t,before=fi3['blocks'][k]) for t,k in zip(table,labels)]
    with pytest.raises(ValueError,match='not a subset'):F.check_reduction(t3,fi3,ri,mi,full3,red,mr)


def test_timed_budget_kills_and_refusal_raises(tmp_path):
    import os,stat
    exe=tmp_path/'slow';exe.write_text('#!/bin/sh\nsleep 30\n');exe.chmod(exe.stat().st_mode|stat.S_IEXEC)
    if not Path('/usr/bin/time').exists() or os.system('/usr/bin/time -v true >/dev/null 2>&1'):
        pytest.skip('GNU time -v absent on this host (present on Sirius)')
    with pytest.raises(F.BudgetExceeded):F.timed(exe,'',tmp_path,1)


def test_timed_passes_absolute_time_output_from_relative_cwd(tmp_path,monkeypatch):
    """v48 FAILED_RETAINED: GNU time resolved a relative -o path against the child's cwd."""
    seen={}
    class Stop(Exception):pass
    def fake(args,**kw):
        seen.update(args=args,cwd=kw['cwd']);raise Stop
    monkeypatch.setattr(F.subprocess,'Popen',fake)
    work=tmp_path/'runs'/'case';work.mkdir(parents=True);monkeypatch.chdir(tmp_path)
    with pytest.raises(Stop):F.timed(Path('/bin/sh'),'',Path('runs/case'),5)
    out=Path(seen['args'][seen['args'].index('-o')+1])
    assert out.is_absolute() and out==work.resolve()/'sh.time' and Path(seen['cwd']).is_absolute()


def test_mr_input_takes_full_header_and_keeps_native_mr_body(native,tmp_path):
    """v49 FAILED_RETAINED: a 0-excitation MR header lists only occupied orbitals (not a prefix)."""
    text=native[next(p for p in native if p.endswith('odd_CV_MR/rcsf.out'))].decode()
    lines=text.split('\n');orbs=lines[3].split()
    extra='   5g-  5g   6s'
    full=tmp_path/'full';full.write_text('\n'.join(lines[:3]+[lines[3]+extra]+lines[4:]))
    mr=tmp_path/'mr';mr.write_text(text);dest=tmp_path/'rcsfmr.inp'
    d=F.mr_input(full,mr,dest)
    out=dest.read_text().split('\n')
    assert out[3]==lines[3]+extra and out[4:]==lines[4:] and d['native_mr_orbitals']==orbs
    assert F.peel_header(dest)==F.peel_header(full)
    # refusals: an MR orbital absent from the full list, and a reordered full header
    full.write_text('\n'.join(lines[:3]+['  '+'   '.join(orbs[:-1])]+lines[4:]))
    with pytest.raises(ValueError,match='ordered subset'):F.mr_input(full,mr,tmp_path/'x')
    full.write_text('\n'.join(lines[:3]+['  '+'   '.join(orbs[1:]+orbs[:1])]+lines[4:]))
    with pytest.raises(ValueError,match='ordered subset'):F.mr_input(full,mr,tmp_path/'y')


def _receipt_pack(tmp_path,*,listed_bytes=None,list_it=True):
    import hashlib,json
    from rya1354_pack_native_stream import pack
    root=tmp_path/'root';run=root/'runs/r';(run/'a').mkdir(parents=True)
    (run/'a/rcsfinteract.stdin').write_text('2\n');big=b'csf'*100;(run/'a/rcsf.out').write_bytes(big)
    rec=lambda b:dict(sha256=hashlib.sha256(b).hexdigest(),bytes=len(b))
    files={'a/rcsfinteract.stdin':rec(b'2\n')}
    if list_it:files['a/rcsf.out']=dict(rec(big),bytes=listed_bytes or len(big))
    (run/'receipt.json').write_text(json.dumps(dict(status='OK',files=files)))
    off={'runs/r/a/rcsf.out':dict(**rec(big),regenerate_stdin='runs/r/a/rcsfinteract.stdin')}
    out=tmp_path/'out';pack(root,[Path('runs/r')],out,off_archive_regenerable=off)
    return out


def test_native_receipt_vouches_for_off_archive_file(tmp_path):
    files,manifest=load_parts(_receipt_pack(tmp_path))
    assert 'runs/r/a/rcsf.out' not in files and set(manifest['off_archive_regenerable'])=={'runs/r/a/rcsf.out'}


@pytest.mark.parametrize('kw,msg',[(dict(listed_bytes=1),'differs from its native receipt'),(dict(list_it=False),'not covered')])
def test_native_receipt_off_archive_refusals(tmp_path,kw,msg):
    with pytest.raises(ValueError,match=msg):_receipt_pack(tmp_path,**kw)


def test_packer_refuses_off_archive_file_changed_on_disk(tmp_path):
    from rya1354_pack_native_stream import pack
    root=tmp_path/'root';(root/'runs/r').mkdir(parents=True);(root/'runs/r/x').write_bytes(b'abc')
    with pytest.raises(ValueError,match='differs on disk'):
        pack(root,[Path('runs/r')],tmp_path/'out',off_archive_regenerable={'runs/r/x':dict(sha256='0'*64,bytes=3,regenerate_stdin='-')})


def test_v50_archive_estimate_output_name():
    """v50 FAILED_RETAINED: rwfnestimate wrote rwfn.inp; the driver had linked rwfn.out."""
    import json
    files,m=load_parts(ROOT/'data/audit/rya1354/rci_feasibility_failed_v50',repo_root=ROOT)
    run='runs/rci_feasibility_v50/cost_odd_J1_n4/run/'
    r=json.loads(files['runs/rci_feasibility_v50/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'rwfn.out' in r['error']
    assert run+F.ESTIMATE_OUTPUT in files and run+'rwfn.out' not in files
    assert 'rwfnestimate' in files[run+'rwfnestimate.stdin'].decode() or files[run+'rwfnestimate.stdin']==b'y\n2\n*\n'


def test_retained_v48_v49_archives_record_their_failures():
    import json
    for v,msg in (('v48','rcsfgenerate.time'),('v49','peel orbital headers differ')):
        files,_=load_parts(ROOT/f'data/audit/rya1354/rci_feasibility_failed_{v}',repo_root=ROOT)
        r=json.loads(files[f'runs/rci_feasibility_{v}/receipt.json'])
        assert r['status']=='FAILED_RETAINED' and msg in r['error']


def test_section_b_steps_run_under_declared_address_space_cap(tmp_path,monkeypatch):
    seen={}
    class Stop(Exception):pass
    def fake(args,**kw):seen['args']=args;raise Stop
    monkeypatch.setattr(F.subprocess,'Popen',fake)
    with pytest.raises(Stop):F.timed(Path('/bin/sh'),'',tmp_path,5,F.MEM_CAP_BYTES)
    a=seen['args'];i=a.index('/usr/bin/prlimit')
    assert a[i+1]==f'--as={F.MEM_CAP_BYTES}' and a[i+2]=='--' and a[-1]=='/bin/sh' and i>a.index('-o')
    (tmp_path/'nocap').mkdir()
    with pytest.raises(Stop):F.timed(Path('/bin/sh'),'',tmp_path/'nocap',5)
    assert '/usr/bin/prlimit' not in seen['args']


def test_rwfnestimate_stdin_falls_back_to_screened_hydrogenic():
    """v51 FAILED_RETAINED: Thomas-Fermi could not build 5g and the remaining-subshell prompt hit EOF."""
    assert F.RWFN_STDIN.split('\n')[:5]==['y','2','*','3','*'] and issubclass(F.MemoryCapReached,F.BudgetExceeded)


def test_memory_marker_is_a_measured_stop_only_under_a_cap(tmp_path):
    import os,stat
    if os.system('/usr/bin/time -v true >/dev/null 2>&1') or not Path('/usr/bin/prlimit').exists():
        pytest.skip('GNU time/prlimit absent on this host (present on Sirius)')
    exe=tmp_path/'alloc';exe.write_text('#!/bin/sh\necho "         Unable to allocate pnt"\nexit 1\n')
    exe.chmod(exe.stat().st_mode|stat.S_IEXEC)
    (tmp_path/'a').mkdir();(tmp_path/'b').mkdir()
    with pytest.raises(F.MemoryCapReached):F.timed(exe,'',tmp_path/'a',5,F.MEM_CAP_BYTES)
    with pytest.raises(ValueError,match='Native refusal'):F.timed(exe,'',tmp_path/'b',5)


def test_v52_archive_completes_section_b_and_matches_section_a():
    """v52 COMPLETE: every rung's DCB table agrees with its own inventories, and the n6 J=1
    count equals the v51 section A odd n6 DCB 2J=2 count; costs grow with CSF count."""
    import json
    files,m=load_parts(ROOT/'data/audit/rya1354/rci_feasibility_complete_v52',repo_root=ROOT)
    r=json.loads(files['runs/rci_feasibility_v52/receipt.json'])
    assert r['status']=='FEASIBILITY_COMPLETED_NO_ATOMIC_PREDICTION' and r['budget_stops']==[]
    assert r['scientific_admission']=='HOLD'
    a,_=load_parts(ROOT/'data/audit/rya1354/rci_feasibility_failed_v51',repo_root=ROOT)
    sa=json.loads(a['runs/rci_feasibility_v51/receipt.json'])
    import hashlib
    assert hashlib.sha256(a['runs/rci_feasibility_v51/receipt.json']).hexdigest()==r['section_a_receipt_sha256']
    rungs=[r['cases'][f'cost_odd_J1_n{n}'] for n in (4,5,6)]
    for c in rungs:
        (row,)=c['DCB']['table']
        assert (row['before'],row['after'],row['MR'])==(c['full']['CSFs'],c['DCB']['CSFs'],c['MR']['CSFs'])
        assert all(s['returncode']==0 and not s['semantic_errors'] for s in c['steps'].values())
        assert c['steps']['rci']['mem_cap_bytes']==F.MEM_CAP_BYTES
    n6=[c for k,c in sa['cases'].items() if 'n6' in k and 'odd' in k]
    assert any(c['DCB']['blocks'].get('1-')==rungs[2]['DCB']['CSFs'] for c in n6)
    for lo,hi in zip(rungs,rungs[1:]):
        assert hi['DCB']['CSFs']>lo['DCB']['CSFs'] and hi['angular_bytes']>lo['angular_bytes']
        assert hi['steps']['rci']['wall_s']>lo['steps']['rci']['wall_s']
    off=m['off_archive_regenerable']
    assert off and all(v['regenerate_stdin'] in files for v in off.values())
