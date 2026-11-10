"""Validate pinned author downloads and receipt corrected BC2016 chemistry data."""
import argparse,json,hashlib,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
PIN='cd2f54ed46687659017b343763e95478ed93653e'
FILES=['ReadMe','ReadMe_Bug_fix_Oct-Nov2022','table1.dat','table6.dat','table7.dat','table7_vNov2022.dat','table8_vNov2022.dat','LICENSE']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--metadata',type=Path,required=True);p.add_argument('--license-metadata',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_h2_tables import read_wide_table,h2_dissociation_eV
    directory=require_codex_path('experiments.bride_molecular_tables')
    if(directory/'receipt.json').exists():raise FileExistsError('Immutable intake exists')
    metadata={r['name']:r for r in json.loads(args.metadata.read_text())};metadata['LICENSE']=json.loads(args.license_metadata.read_text())
    for name in FILES:
        path=directory/(name+'.download');data=path.read_bytes();blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if blob!=metadata[name]['sha']:raise ValueError('Pinned Git blob differs: '+name)
        if(directory/name).exists():raise FileExistsError('Existing source file')
    for name,kind in [('table6.dat','partition'),('table7.dat','equilibrium'),('table7_vNov2022.dat','equilibrium'),('table8_vNov2022.dat','partition')]:read_wide_table(directory/(name+'.download'),kind)
    d0=h2_dissociation_eV(directory/'table1.dat.download')
    result=dict(schema='codex.bride_molecular_source/1',source_commit=PIN,scientific_admission='HOLD',generator=fingerprint(__file__),
        source='https://github.com/barklem/public-data/tree/'+PIN+'/partition-functions_and_equilibrium-constants',
        license='BSD-3-Clause; copyright notice retained',adopted_H2_D0_eV=d0,files={},
        scope='Full source tables held; only neutral H2 interface being qualified. Not adoption of production CNO chemistry. Original CDS partial holding unchanged.')
    for name in FILES:
        path=directory/name;(directory/(name+'.download')).rename(path);path.chmod(0o444)
        relative=name if name=='LICENSE'else'partition-functions_and_equilibrium-constants/'+name
        result['files'][name]=dict(git_blob_sha=metadata[name]['sha'],fingerprint=fingerprint(path),url=f'https://raw.githubusercontent.com/barklem/public-data/{PIN}/{relative}')
    (directory/'receipt.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Validated8 immutable author files; D0',d0)


if __name__=='__main__':main()
