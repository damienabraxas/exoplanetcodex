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
