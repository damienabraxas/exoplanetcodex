"""Offline, hash-checked native-table summary of the exposed continuation capsule."""
from pathlib import Path
import hashlib
import json
import re
import tarfile

ROOT=Path(__file__).resolve().parents[1]
NUMBER=r'[\d.D+\-]+'


def main(destination=None):
    capsule=ROOT/'data/audit/rya1354/continuation_20261108_compact.tar.gz'
    with tarfile.open(capsule) as archive:
        files={m.name:archive.extractfile(m).read() for m in archive.getmembers() if m.isfile()}
    def text(name):
        return files[name].decode()
    receipts={}
    for name,content in files.items():
        if not name.endswith('/receipt.json'):
            continue
        receipt=json.loads(content)
        parent=name.rsplit('/',1)[0]
        for relative,record in receipt['files'].items():
            key=parent+'/'+relative
            if key not in files:
                if not Path(relative).name.startswith('mcp.'):
                    raise ValueError(f'Unexpected omitted artifact: {key}')
                continue
            if hashlib.sha256(files[key]).hexdigest()!=record['sha256']:
                raise ValueError(f'Artifact changed: {key}')
        receipts[parent]=dict(status=receipt['status'],error=receipt.get('error'),receipt_sha256=hashlib.sha256(content).hexdigest())
    def e1(path):
        pattern=(rf'f2\s+(\d+)\s+(\d+)\s+[+-]\s+f1\s+(\d+)\s+(\d+)\s+[+-]\s+([\d.]+)\s+C\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})'
                 rf'\s+B\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})')
        result=[]
        for u,uj,l,lj,en,ac,gc,sc,ab,gb,sb in re.findall(pattern,text(path)):
            conv=lambda v:float(v.replace('D','E'))
            ab,ac,gb,gc=map(conv,(ab,ac,gb,gc))
            result.append(dict(upper_ordinal=int(u),upper_J=int(uj),lower_ordinal=int(l),lower_J=int(lj),
                               native_energy_cm1=float(en),A_B_s1=ab,A_C_s1=ac,gf_B=gb,gf_C=gc,
                               gauge_fraction_A=abs(ab-ac)/max(ab,ac)))
        if not result:
            raise ValueError(f'No native E1 rows in {path}')
        return result
    coverage='runs/coverage_v8'
    if receipts[coverage]['status']!='NATIVE_OPERATORS_COMPLETED_REVIEW_REQUIRED':
        raise ValueError('Coverage run did not complete')
    expected={'triplet_3s':('3s_3S',1),'triplet_3p':('3p_3P',3),
              'quintet_3d':('3d_5D',5),'quintet_4d':('4d_5D',5)}
    for case,(label,count) in expected.items():
        if text(f'{coverage}/{case}/state.lsj.lbl').count(label)!=count:
            raise ValueError(f'Native LS composition identity refused: {case}')
    coverage_tables={case:e1(f'{coverage}/{case}/lower.upper.t') for case in
                     ['UV_E1','IR_3d_E1','optical_4d_E1','triplet_8446_E1']}
    if [len(coverage_tables[k]) for k in ['IR_3d_E1','optical_4d_E1','triplet_8446_E1']] != [9,9,3]:
        raise ValueError('Reference multiplet component inventory changed')
    ladder={f'n{n}':e1(f'runs/valence_sd_v4/n{n}_E1/lower.upper.t') for n in (3,4)}
    result=dict(schema='rya1354.native-continuation-summary.v1',scientific_admission='HOLD',
                exposure='NONBLIND_EXPLORATORY',capsule_sha256=hashlib.sha256(capsule.read_bytes()).hexdigest(),
                omitted='mcp.* angular caches remain on Sirius; all other receipt-listed artifacts hash-checked',
                receipts=receipts,ladder=ladder,coverage=coverage_tables,
                forbidden=dict(status='NATIVE_M1_AND_E2_COMPUTED_NOT_ADOPTED',
                               table=f'{coverage}/ground_M1_E2/lower.upper.t',
                               table_sha256=hashlib.sha256(files[f'{coverage}/ground_M1_E2/lower.upper.t']).hexdigest(),
                               warning='Identical states in f1/f2 yield duplicate directional records; do not sum duplicates. Very large E2 gauge disagreements. No E1 conversion, uncertainty or lifetime supplied.'),
                warning='No active-space convergence, core-valence closure, Hamiltonian sensitivity or calibrated uncertainty established; frozen-core triplet 3p differs from other relaxed reference cases')
    (destination or ROOT/'data/audit/rya1354/nonblind_v1/native_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:len(v) for k,v in coverage_tables.items()})


if __name__=='__main__':
    main()
