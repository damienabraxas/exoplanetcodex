from pathlib import Path
import json
import struct
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_stream_transport import transport,record_ranges,inventory,configurations,EXTENDED

@pytest.fixture(scope='module')
def capsule():
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_full_valence_complete_v28')
    return files

@pytest.fixture
def case(tmp_path,capsule):
    prefix='runs/fac_full_valence_v28/ground_coherent10/'
    for fn in ('o.basis','o_export_000.c','o_export_000.cm'):(tmp_path/fn).write_bytes(capsule[prefix+fn])
    return tmp_path

def check(p):return transport(p/'o.basis',p/'o_export_000.c',p/'o_export_000.cm')

def test_existing_complete_audit_reproduced(case,capsule):
    receipt=json.loads(capsule['runs/fac_full_valence_v28/receipt.json'])
    assert check(case)==receipt['cases']['ground_coherent10']['transport']

@pytest.mark.parametrize('mutation',['sign','id','missing','duplicate','marker'])
def test_corrupt_native_or_binary_evidence_is_rejected(case,mutation):
    p=case/'o.basis';lines=p.read_text().splitlines();start=next(i for i,x in enumerate(lines) if 'Mixing Coefficients' in x)
    at=next(i for i in range(start+1,len(lines)) if lines[i].strip() and not lines[i].startswith('#') and float(lines[i].split()[9])!=0)
    f=lines[at].split()
    if mutation=='sign':f[9]=format(-float(f[9]),'.8E');lines[at]=' '.join(f)
    elif mutation=='id':f[0]='999999';lines[at]=' '.join(f)
    elif mutation=='missing':lines.pop(at)
    elif mutation=='duplicate':lines.insert(at,lines[at])
    else:
        q=case/'o_export_000.cm';data=bytearray(q.read_bytes());data[0]^=1;q.write_bytes(data)
    p.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError):check(case)

def test_prospective_extended_inventory():
    configs=configurations(EXTENDED);counts=inventory(EXTENDED,configs)
    assert len(configs)==114 and sum(counts.values())==6600 and max(counts.values())==916
    assert sum((j+1)*n for (p,j),n in counts.items())==35960


def test_extended_native_input_preserves_fixed_potential_protocol():
    from rya1354_fac_radial_ci import input_text
    s=input_text()
    assert s.count('Config(')==114 and s.count('Structure(')==1
    assert s.count("RestorePotential('source.pot')")==1 and 'OptimizeRadial(' not in s
    assert "BasisTable('o.basis',0)" in s and "BasisTable('o_export',10)" in s
    assert "SetConfigEnergyMode(-1)" in s and "SetOption('radial:config_energy',-1)" in s
    assert s.count('WaveFuncTable(')==33


def test_principal_complex_identity_preserves_previous_transport(case,capsule):
    from rya1354_fac_stream_transport_v34 import transport as corrected
    expected=json.loads(capsule['runs/fac_full_valence_v28/receipt.json'])['cases']['ground_coherent10']['transport']
    assert corrected(case/'o.basis',case/'o_export_000.c',case/'o_export_000.cm')==expected


def test_principal_complex_corruption_rejected(case):
    from rya1354_fac_stream_transport_v34 import transport as corrected
    p=case/'o.basis';lines=p.read_text().splitlines()
    at=next(i for i,x in enumerate(lines) if x.strip() and not x.startswith('#') and 'Basis Table' not in x)
    f=lines[at].split();f[7]='99*4';lines[at]=' '.join(f);p.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError,match='identity'):corrected(p,case/'o_export_000.c',case/'o_export_000.cm')


def test_complete_extended_capsule_replays_and_both_failures_remain_retained():
    from rya1354_fac_radial_ci_audit import summarize
    report=summarize(ROOT/'data/audit/rya1354')
    assert report['retained_v33_status']==report['retained_v34_status']=='FAILED_RETAINED'
    assert not report['native_recalculation'] and not report['adoption_allowed']
    case=report['cases']['ground_coherent10']
    assert case['transport']['levels']==case['transport']['CSFs']==6600
    assert case['transport']['signed_printed_coefficients_verified']==4118356
    assert case['post_CI_orbitals_unchanged']
    assert case['potential']['protected_header_grid_potential_byte_identical']
