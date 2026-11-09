"""Byte-copy compact reference evidence; receipt all raw files without size bypass."""
import argparse,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--output-opacity-supplement',action='store_true');a=p.parse_args()
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
    a.output.mkdir(parents=True,exist_ok=False);records=[]
    for name,kind in groups.items():
        source=a.raw_root/name
        if not source.is_dir():raise FileNotFoundError(source)
        for path in sorted(source.rglob('*')):
            rel=path.relative_to(source)
            if not path.is_file()or path.is_symlink()or 'input'in rel.parts:continue
            receipt=fingerprint(path);copy=kind=='analysis'or len(rel.parts)==1 and path.suffix in ['.json','.log','.nml','.txt']
            if receipt['bytes']>5*1024*1024:copy=False
            target=a.output/name/rel
            if copy:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
                if fingerprint(target)!=receipt or fingerprint(path)!=receipt:raise ValueError('Archive byte identity failed')
            records.append(dict(raw=str(path),receipt=receipt,archived=str(target.relative_to(a.output))if copy else None,disposition=kind))
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',files=records,limitations=['Invalid single-ray v1 arrays remain external and inadmissible; receipts do not validate them.','Native binary/raw output remains external with byte receipts; compact analysis copied exactly.','Older exploratory readers mislabeled native per-Angstrom intensity as per-Hz; only dense-analysis-v3 and later qualified readers used here.'])
    (a.output/'archive_receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(raw_files=len(records),copied=sum(x['archived']is not None for x in records))))


if __name__=='__main__':main()
