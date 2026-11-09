"""Copy immutable numerical-control evidence; receipt raw files left external."""
import argparse,hashlib,json,shutil
from pathlib import Path


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--exclude-transition-run',action='store_true');p.add_argument('--transition-only',action='store_true');p.add_argument('--license-source',type=Path);a=p.parse_args()
 if a.license_source:
  target=a.output/'LICENSE.lightweaver'
  if not a.output.is_dir()or target.exists():raise ValueError('Require existing archive and a fresh license path')
  shutil.copyfile(a.license_source,target)
  if target.read_bytes()!=a.license_source.read_bytes():raise ValueError('License bytes changed')
  return
 baseline=a.source_root/'rya1281-solar-transition-v1/baseline'
 if not a.exclude_transition_run and not any((baseline/f).is_file()for f in ['report.json','failure.json']):raise ValueError('Transition baseline is still incomplete; do not snapshot a live run')
 groups=['rya1281-solar-anderson-dry-v1','rya1281-solar-anderson-v1','rya1281-solver-conditioning-v1','rya1281-continuum-overlap-v2','rya1281-continuum-overlap-v3','rya1281-continuum-overlap-v4','rya1281-continuum-overlap-v5','rya1281-mali-operator-audit-v2','rya1281-mali-operator-audit-v3','rya1281-overlap-plugin-v1','rya1281-transition-plugin-v1','rya1281-solar-transition-dry-v1','rya1281-solar-transition-v1']
 if a.exclude_transition_run:groups.remove('rya1281-solar-transition-v1')
 if a.transition_only:
  if a.exclude_transition_run:raise ValueError('Conflicting snapshot scopes')
  groups=['rya1281-solar-transition-dry-v1','rya1281-solar-transition-v1']
  for directory in baseline.parent.iterdir():
   if directory.is_dir()and not any((directory/f).is_file()for f in ['report.json','failure.json']):raise ValueError('Owned transition case is still live/incomplete')
 rows=[];a.output.mkdir(parents=True,exist_ok=False)
 for name in groups:
  source=a.source_root/name
  for path in sorted(source.rglob('*')):
   if not path.is_file():continue
   relative=Path(name)/path.relative_to(source);data=path.read_bytes()
   build=name in ['rya1281-overlap-plugin-v1','rya1281-transition-plugin-v1']
   archived=(path.name in ['receipt.json','overlap_guard.cpp','SimdFullIterationTemplates.hpp','build.log'])if build else path.name!='rates.npz'
   rows.append(dict(path=str(relative),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),archived=archived))
   if archived:
    dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    if dest.read_bytes()!=data:raise ValueError('Archive bytes changed')
 (a.output/'receipt.json').write_text(json.dumps(dict(schema='codex.bride_solver_acceleration_archive/1',scientific_admission='HOLD',files=rows,note='Exact compact archival copies. Installed backend untouched; compiled plugins and full J/rate arrays remain externally receipted. Failure seeds are initialization only.'),indent=2)+'\n')

if __name__=='__main__':main()
