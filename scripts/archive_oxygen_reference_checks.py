"""Byte-copy compact reference evidence; receipt all raw files without size bypass."""
import argparse,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--output-opacity-supplement',action='store_true');p.add_argument('--boundary-eos-supplement',action='store_true');p.add_argument('--molecular-preflight-supplement',action='store_true');p.add_argument('--hop-provider-supplement',action='store_true');p.add_argument('--hop-trace-supplement',action='store_true');a=p.parse_args()
    if sum([a.output_opacity_supplement,a.boundary_eos_supplement,a.molecular_preflight_supplement,a.hop_provider_supplement,a.hop_trace_supplement])>1:raise ValueError('Choose one supplement')
    extra_sources={}
    groups={
      'rya1281-independent-oxygen-populations-v1':'analysis',
      'rya1371-oxygen-reference-v1':'native',
      'rya1371-oxygen-spectrum-v1':'native',
      'rya1371-oxygen-dense-v1':'native',
      'rya1371-oxygen-frequency-v1':'native',
      'rya1371-oxygen-angle-v1':'failed-native',
      'rya1371-oxygen-angle6-v1':'native',
      'rya1371-bride-reference-transfer-v1':'invalid',
      'rya1371-bride-reference-transfer-v2':'analysis',
      'rya1371-bride-dense-transfer-v1':'analysis',
      'rya1371-oxygen-dense-analysis-v3':'analysis',
      'rya1371-oxygen-frequency-analysis-v1':'analysis',
      'rya1371-oxygen-angle6-analysis-v1':'analysis',
      'rya1371-dense-spectrum-comparison-v3':'analysis',
      'rya1371-reference-refinement-v2':'analysis',
    }
    if a.output_opacity_supplement:groups={'rya1371-oxygen-output-bf-v1':'native','rya1371-oxygen-output-bf-analysis-v1':'analysis','rya1371-output-opacity-audit-v1':'analysis'}
    if a.boundary_eos_supplement:groups={
      'rya1371-bride-zero-cap-transfer-v1':'analysis','rya1371-zero-cap-comparison-v1':'analysis',
      'rya1371-consumed-state-v1':'invalid','rya1371-consumed-state-cutoff1K-v1':'invalid',
      'rya1371-consumed-state-v2':'analysis','rya1371-consumed-state-cutoff1K-v2':'analysis',
      'rya1371-oxygen-species-v1':'native','rya1371-oxygen-species-analysis-v1':'analysis',
      'rya1371-oxygen-no-molecules-v1':'failed-native',
      'rya1371-oxygen-molecular-cutoff1K-v1':'native','rya1371-oxygen-molecular-cutoff1K-analysis-v1':'analysis',
      'rya1371-molecular-reservoir-v1':'analysis'}
    if a.molecular_preflight_supplement:groups={'rya1371-bride-molecular-preflight-v2':'analysis'}
    if a.hop_provider_supplement:groups={'rya1371-hop-state-provider-v2':'analysis'}
    if a.hop_trace_supplement:
        from config.constants import require_codex_path
        build=require_codex_path('engines.hop_trace_provider_build');groups={'rya1371-hop-trace-provider-v1':'analysis','rya1371-hop-opacity-provider-v2':'analysis','provider_build':'native'};extra_sources['provider_build']=build
        for n in [1,2,3]:
            name='prototype_build_'+str(n);groups[name]='receipt-only';extra_sources[name]=build.parent/('hop-provider-build-v'+str(n))
    a.output.mkdir(parents=True,exist_ok=False);records=[]
    for name,kind in groups.items():
        source=extra_sources.get(name,a.raw_root/name)
        if not source.is_dir():raise FileNotFoundError(source)
        for path in sorted(source.rglob('*')):
            rel=path.relative_to(source)
            if not path.is_file()or path.is_symlink()or 'input'in rel.parts:continue
            receipt=fingerprint(path);copy=kind=='analysis'or kind!='receipt-only'and len(rel.parts)==1 and path.suffix in ['.json','.log','.nml','.txt']
            if receipt['bytes']>5*1024*1024:copy=False
            target=a.output/name/rel
            if copy:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
                if fingerprint(target)!=receipt or fingerprint(path)!=receipt:raise ValueError('Archive byte identity failed')
            records.append(dict(raw=str(path),receipt=receipt,archived=str(target.relative_to(a.output))if copy else None,disposition=kind))
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',files=records,limitations=['Invalid single-ray v1 arrays remain external and inadmissible; receipts do not validate them.','Native binary/raw output remains external with byte receipts; compact analysis copied exactly.','Older exploratory readers mislabeled native per-Angstrom intensity as per-Hz; only dense-analysis-v3 and later qualified readers used here.'])
    (a.output/'archive_receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(raw_files=len(records),copied=sum(x['archived']is not None for x in records))))


if __name__=='__main__':main()
