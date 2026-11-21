import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts  # noqa: E402
from rya1354_li_scheme_sizing import active, inventory_stream, references, stdin_for  # noqa: E402
from rya1354_mr_inventory import inventory  # noqa: E402


@pytest.fixture(scope='module')
def native():
    files,_=load_parts(ROOT/'data/audit/rya1354/mr_inventory_v5')
    return files


def test_stream_inventory_matches_tested_inventory_on_native_lists(native,tmp_path):
    cases=[p for p in native if p.endswith('/rcsf.out')]
    assert len(cases)==8
    for path in cases:
        parity='even' if '/even_' in path else 'odd';core='_CV_' in path
        old=inventory(native[path].decode(),parity,core)
        f=tmp_path/'rcsf.out';f.write_bytes(native[path])
        stage='rci' if core else 'mcdhf'
        new=inventory_stream(f,parity,stage)
        assert (new['CSFs'],new['blocks'])==(old['CSFs'],old['blocks'])


def test_stream_inventory_keeps_refusals(native,tmp_path):
    path=next(p for p in native if p.endswith('odd_CV_MR/rcsf.out'))
    text=native[path].decode();f=tmp_path/'rcsf.out'
    for bad,msg in ((text.replace('1s ( 2)','1s ( 1)',1),'occupation'),
                    (text+'\n'.join(text.splitlines()[5:8])+'\n','Duplicate')):
        f.write_text(bad)
        with pytest.raises(ValueError,match=msg):inventory_stream(f,'odd','rci')
    f.write_text(text)
    with pytest.raises(ValueError,match='parity'):inventory_stream(f,'even','rci')
    vv=next(p for p in native if p.endswith('odd_VV_MR/rcsf.out'))
    f.write_bytes(native[path])  # CV list has open 2s: refused as an mcdhf-stage list
    with pytest.raises(ValueError,match='occupation'):inventory_stream(f,'odd','mcdhf')
    f.write_bytes(native[vv]);assert inventory_stream(f,'odd','mcdhf')['CSFs']>0


def test_li_table1_families_and_layer_truncation():
    kept,omitted=references('rci','even',6)
    assert '1s(2,i)2s(2,*)2p(4,*)' in kept and '1s(2,i)2s(2,*)2p(3,*)8p(1,*)' in omitted
    assert '1s(2,i)2s(1,*)2p(4,*)7s(1,*)' in omitted and '1s(2,i)2s(1,*)2p(4,*)6d(1,*)' in kept
    kept,omitted=references('rci','even',9)
    assert not omitted and len(kept)==1+6+2+5+4
    kept,omitted=references('mcdhf','odd',6)
    assert not omitted and all('2s(2,i)' in c for c in kept) and len(kept)==4+3
    assert active(6)=='6s,6p,6d,6f,6g,6h' and active(9)=='9s,9p,9d,9f,9g,9h'
    assert stdin_for('mcdhf','odd',6).endswith('\n\n6s,6p,6d,6f,6g,6h\n2,4\n2\nn\n')


def _sidecar_pack(tmp_path,*,mutate=lambda inv,off:None,receipt=None):
    import hashlib,json
    from rya1354_pack_native_stream import pack
    root=tmp_path/'root';run=root/'runs/r';run.mkdir(parents=True);(root/'ev').mkdir()
    (run/'receipt.json').write_text(json.dumps(receipt or dict(status='OK')))
    (run/'a').mkdir();(run/'a/rcsfgenerate.stdin').write_text('*\n')
    big=b'csf'*100
    rec=lambda b:dict(sha256=hashlib.sha256(b).hexdigest(),bytes=len(b))
    inv={'receipt.json':rec((run/'receipt.json').read_bytes()),'a/rcsfgenerate.stdin':rec(b'*\n'),'a/rcsf.out':rec(big)}
    off={'runs/r/a/rcsf.out':dict(**rec(big),regenerate_stdin='runs/r/a/rcsfgenerate.stdin')}
    mutate(inv,off)
    (root/'ev/files.json').write_text(json.dumps(dict(run='runs/r',files=inv)))
    out=tmp_path/'out'
    pack(root,[Path('ev'),Path('runs/r')],out,receipt_inventories={'runs/r/receipt.json':'ev/files.json'},off_archive_regenerable=off)
    return out


def test_sidecar_inventory_accepts_declared_off_archive_file(tmp_path):
    files,manifest=load_parts(_sidecar_pack(tmp_path))
    assert 'runs/r/a/rcsf.out' not in files and set(manifest['off_archive_regenerable'])=={'runs/r/a/rcsf.out'}


@pytest.mark.parametrize('mutate,receipt,msg',[
    (lambda inv,off:inv['a/rcsfgenerate.stdin'].update(sha256='0'*64),None,'bytes/hash differ'),
    (lambda inv,off:off['runs/r/a/rcsf.out'].update(bytes=1),None,'missing or changed'),
    (lambda inv,off:inv.pop('a/rcsfgenerate.stdin'),None,'absent from sidecar'),
    (lambda inv,off:inv.pop('a/rcsf.out'),None,'not covered'),
    (lambda inv,off:off['runs/r/a/rcsf.out'].update(regenerate_stdin='runs/r/none'),None,'regeneration input'),
    (lambda inv,off:None,dict(status='OK',files={}),'already carries'),
])
def test_sidecar_inventory_refusals(tmp_path,mutate,receipt,msg):
    with pytest.raises(ValueError,match=msg):_sidecar_pack(tmp_path,mutate=mutate,receipt=receipt)


def test_receipt_without_files_and_without_sidecar_still_fails(tmp_path):
    from rya1354_pack_native_stream import pack
    root=tmp_path/'root';(root/'runs/r').mkdir(parents=True);(root/'runs/r/receipt.json').write_text('{}')
    with pytest.raises(KeyError):pack(root,[Path('runs/r')],tmp_path/'out')


def test_actual_v47_sizing_archive():
    import json
    files,manifest=load_parts(ROOT/'data/audit/rya1354/li_scheme_sizing_v47')
    r=json.loads(files['runs/li_scheme_sizing_v47/receipt.json'])
    assert r['status']=='SIZING_COMPLETED_NO_ATOMIC_PREDICTION' and len(r['cases'])==16
    assert len(manifest['off_archive_regenerable'])==8
    assert all(k.split('/')[2].startswith('rci_') for k in manifest['off_archive_regenerable'])
    for name,c in r['cases'].items():
        rcsf=f'runs/li_scheme_sizing_v47/{name}/rcsf.out'
        if name.startswith('mcdhf_'):
            assert inventory_stream_bytes(files[rcsf],name)==(c['CSFs'],c['blocks'])
        assert c['CSFs']==sum(c['blocks'].values())
        assert files[f'runs/li_scheme_sizing_v47/{name}/rcsfgenerate.stdin'].decode()==stdin_for(*name.split('_')[:2],int(name[-1]))
    assert r['cases']['mcdhf_even_n9']['blocks']['2+']==146166 and r['cases']['rci_even_n9']['blocks']['2+']==2107505


def inventory_stream_bytes(data,name):
    import tempfile
    with tempfile.NamedTemporaryFile() as f:
        f.write(data);f.flush()
        out=inventory_stream(f.name,name.split('_')[1],name.split('_')[0])
    return out['CSFs'],out['blocks']
