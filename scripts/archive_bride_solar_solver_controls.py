"""Preserve exact compact solver-attempt evidence and receipt omitted raw arrays."""
import argparse,hashlib,json,shutil
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 a.output.mkdir(parents=True,exist_ok=False);rows=[]
 for name in ['rya1281-solar-matched-v2','rya1281-solar-matched-dry-v2','rya1281-solar-mali-v1','rya1281-mali-operator-audit-v1']:
  source=a.source_root/name
  for path in sorted(source.rglob('*')):
   if not path.is_file():continue
   relative=Path(name)/path.relative_to(source);data=path.read_bytes()
   archived=path.name!='rates.npz'
   rows.append(dict(path=str(relative),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),archived=archived))
   if archived:
    dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    if dest.read_bytes()!=data:raise ValueError('Archive bytes changed')
 (a.output/'receipt.json').write_text(json.dumps(dict(scientific_admission='HOLD',note='All failed runs remain incomplete; failure seeds are initialization only. Original full rates remain external and receipted; no solver or science admission.',files=rows),indent=2)+'\n')

if __name__=='__main__':main()
