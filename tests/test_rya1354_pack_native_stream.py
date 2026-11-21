from pathlib import Path
import hashlib,json,os,sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_pack_native_stream import pack
from rya1354_native_parts import load_parts


def test_single_pass_bundle_reconstructs_fragmented_random_zero_small_and_empty_files(tmp_path):
    run=tmp_path/'run';run.mkdir();data={'random.bin':os.urandom(12000),'zero.bin':b'\0'*9000,'small':b'exact\n','empty':b''}
    for name,value in data.items():(run/name).write_bytes(value)
    receipt={'files':{k:{'bytes':len(v),'sha256':hashlib.sha256(v).hexdigest()} for k,v in data.items()}}
    (run/'receipt.json').write_text(json.dumps(receipt));pack(tmp_path,[Path('run')],tmp_path/'out',chunk=2048)
    restored,manifest=load_parts(tmp_path/'out')
    assert {k:restored['run/'+k] for k in data}==data
    assert manifest['fragments'] and all(v['bytes']<=5*1024*1024 for v in manifest['parts'].values())


def test_original_receipt_corruption_is_rejected_by_existing_decoder(tmp_path):
    run=tmp_path/'run';run.mkdir();(run/'data').write_bytes(b'actual')
    (run/'receipt.json').write_text(json.dumps({'files':{'data':{'bytes':6,'sha256':'0'*64}}}))
    with pytest.raises(ValueError,match='missing or changed'):pack(tmp_path,[Path('run')],tmp_path/'out')


def _earlier(tmp_path,payload):
    old=tmp_path/'repo/old';old.mkdir(parents=True);(old/'big.bin').write_bytes(payload)
    pack(tmp_path/'repo',[Path('old')],tmp_path/'repo/archive_old',chunk=2048)
    return hashlib.sha256(payload).hexdigest()


def test_v3_references_restore_identical_bytes_from_pinned_earlier_archive(tmp_path):
    payload=os.urandom(9000);digest=_earlier(tmp_path,payload)
    stage=tmp_path/'stage';run=stage/'new';run.mkdir(parents=True);(run/'own.txt').write_bytes(b'new\n')
    receipt={'files':{'own.txt':{'bytes':4,'sha256':hashlib.sha256(b'new\n').hexdigest()},'copy.bin':{'bytes':9000,'sha256':digest}}}
    (run/'receipt.json').write_text(json.dumps(receipt))
    refs={'new/copy.bin':dict(sha256=digest,bytes=9000,archive='archive_old',archive_path='old/big.bin')}
    pack(stage,[Path('new')],tmp_path/'repo/archive_new',references=refs,repo_root=tmp_path/'repo')
    files,manifest=load_parts(tmp_path/'repo/archive_new',repo_root=tmp_path/'repo')
    assert files['new/copy.bin']==payload and manifest['schema'].endswith('.v3')
    assert manifest['native_paths_preserved_byte_for_byte']==3 and 'new/copy.bin' not in manifest['files']
    with pytest.raises(ValueError,match='repository root'):load_parts(tmp_path/'repo/archive_new')


def test_v3_rejects_wrong_reference_and_changed_earlier_archive(tmp_path):
    payload=os.urandom(5000);digest=_earlier(tmp_path,payload)
    stage=tmp_path/'stage';run=stage/'new';run.mkdir(parents=True)
    (run/'receipt.json').write_text(json.dumps({'files':{'copy.bin':{'bytes':5000,'sha256':digest}}}))
    bad={'new/copy.bin':dict(sha256='0'*64,bytes=5000,archive='archive_old',archive_path='old/big.bin')}
    with pytest.raises(ValueError,match='bytes/hash'):
        pack(stage,[Path('new')],tmp_path/'repo/a1',references=bad,repo_root=tmp_path/'repo')
    good={'new/copy.bin':dict(sha256=digest,bytes=5000,archive='archive_old',archive_path='old/big.bin')}
    pack(stage,[Path('new')],tmp_path/'repo/a2',references=good,repo_root=tmp_path/'repo')
    m=tmp_path/'repo/archive_old/manifest.json';m.write_text(m.read_text()+' ')
    import rya1354_native_parts;rya1354_native_parts._REFERENCE_CACHE.clear()
    with pytest.raises(ValueError,match='manifest differs'):load_parts(tmp_path/'repo/a2',repo_root=tmp_path/'repo')


def test_v3_self_reference_restores_a_file_stored_once_in_the_same_archive(tmp_path):
    stage=tmp_path/'stage';(stage/'run/b').mkdir(parents=True);(stage/'run/c').mkdir()
    payload=os.urandom(7000);(stage/'run/b/o.basis').write_bytes(payload)
    d=hashlib.sha256(payload).hexdigest()
    (stage/'run/receipt.json').write_text(json.dumps({'files':{'b/o.basis':{'bytes':7000,'sha256':d},'c/o.basis':{'bytes':7000,'sha256':d}}}))
    refs={'run/c/o.basis':dict(sha256=d,bytes=7000,archive='.',archive_path='run/b/o.basis')}
    pack(stage,[Path('run')],tmp_path/'a',chunk=2048,references=refs)
    files,manifest=load_parts(tmp_path/'a')
    assert files['run/c/o.basis']==files['run/b/o.basis']==payload and manifest['referenced_archives']=={'.':'.'}
    bad={'run/c/o.basis':dict(sha256=d,bytes=7000,archive='.',archive_path='run/x')}
    with pytest.raises(ValueError,match='not stored here'):pack(stage,[Path('run')],tmp_path/'b',references=bad)
