"""Byte-copy compact reference evidence; receipt all raw files without size bypass."""
import argparse,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--output-opacity-supplement',action='store_true');p.add_argument('--boundary-eos-supplement',action='store_true');p.add_argument('--molecular-preflight-supplement',action='store_true');p.add_argument('--hop-provider-supplement',action='store_true');p.add_argument('--hop-trace-supplement',action='store_true');p.add_argument('--hop-fixed-transfer-supplement',action='store_true');p.add_argument('--hop-formal-stark-supplement',action='store_true');p.add_argument('--hop-helium-supplement',action='store_true');p.add_argument('--hop-profile-source-supplement',action='store_true');p.add_argument('--hop-table-profile-supplement',action='store_true');p.add_argument('--hop-atmospheric-profile-supplement',action='store_true');p.add_argument('--accurate-native-profile-supplement',action='store_true');p.add_argument('--vdw-wavelength-supplement',action='store_true');p.add_argument('--native-vacuum-centres-supplement',action='store_true');p.add_argument('--multi-width-supplement',action='store_true');p.add_argument('--original-source-widths-supplement',action='store_true');a=p.parse_args()
    if sum([a.output_opacity_supplement,a.boundary_eos_supplement,a.molecular_preflight_supplement,a.hop_provider_supplement,a.hop_trace_supplement,a.hop_fixed_transfer_supplement,a.hop_formal_stark_supplement,a.hop_helium_supplement,a.hop_profile_source_supplement,a.hop_table_profile_supplement,a.hop_atmospheric_profile_supplement,a.accurate_native_profile_supplement,a.vdw_wavelength_supplement,a.native_vacuum_centres_supplement,a.multi_width_supplement])>1:raise ValueError('Choose one supplement')
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
    if a.hop_fixed_transfer_supplement:groups={**{'rya1371-hop-fixed-transfer-dry-v'+str(n):'failed-native'if n<4 else'analysis'for n in [1,2,3,4]},'rya1371-hop-fixed-transfer-v1':'analysis','rya1371-hop-fixed-transfer-analysis-v1':'analysis'}
    if a.hop_formal_stark_supplement:groups={
      'rya1371-oxygen-linear-output-v1':'native','rya1371-oxygen-linear-analysis-v1':'analysis',
      'rya1371-hop-linear-transfer-v1':'analysis','rya1371-hop-linear-comparison-v1':'analysis',
      'rya1371-oxygen-linear-no-stark-v1':'native','rya1371-oxygen-linear-no-stark-analysis-v1':'analysis',
      'rya1371-hop-no-stark-transfer-v1':'analysis','rya1371-hop-no-stark-comparison-v1':'analysis',
      'rya1371-hop-formal-stark-audit-v1':'analysis'}
    if a.hop_helium_supplement:groups={'rya1371-hop-no-he-transfer-v1':'invalid','rya1371-hop-no-he-transfer-v2':'analysis','rya1371-hop-no-he-comparison-v2':'analysis','rya1371-hop-helium-audit-v1':'analysis'}
    if a.hop_profile_source_supplement:groups={'rya1371-hop-no-he-transfer-v3':'analysis','rya1371-hop-he-guard-replay-v1':'analysis','rya1371-oxygen-native-table-audit-v1':'invalid','rya1371-oxygen-native-table-audit-v2':'analysis'}
    if a.hop_table_profile_supplement:groups={
      'rya1371-hop-table-centres-dry-v1':'failed-native','rya1371-hop-table-centres-v1':'failed-native','rya1371-hop-table-centres-dry-v2':'analysis',
      'rya1371-hop-table-centres-v2':'analysis','rya1371-hop-table-strengths-v1':'analysis','rya1371-hop-table-combined-v1':'analysis',
      'rya1371-hop-table-centres-comparison-v2':'analysis','rya1371-hop-table-strengths-comparison-v1':'analysis','rya1371-hop-table-combined-comparison-v1':'analysis',
      'rya1371-hop-table-audit-v1':'analysis','rya1371-qlg-profile-probe-v1':'receipt-only','rya1371-qlg-profile-probe-v2':'kernel'}
    if a.hop_atmospheric_profile_supplement:groups={'rya1371-optical-damping-samples-v1':'receipt-only','rya1371-optical-damping-samples-v2':'analysis','rya1371-qlg-atmospheric-probe-v1':'receipt-only','rya1371-qlg-atmospheric-probe-v2':'receipt-only','rya1371-qlg-atmospheric-probe-v3':'kernel'}
    if a.accurate_native_profile_supplement:
        from config.constants import require_codex_path
        build=require_codex_path('engines.openmulti3d_profile_reference_build');extra_sources['profile_build']=build;groups={'profile_build':'native'}
        for n in [1,2]:
            name='prototype_profile_build_'+str(n);groups[name]='receipt-only';extra_sources[name]=build.parent/('profile-reference-build-v'+str(n))
        groups.update({'rya1371-oxygen-accurate-profile-v1':'native','rya1371-oxygen-accurate-profile-analysis-v1':'analysis','rya1371-hop-accurate-profile-transfer-v1':'analysis','rya1371-hop-accurate-profile-comparison-v1':'analysis','rya1371-accurate-profile-audit-v1':'analysis','rya1371-original-profile-replay-v1':'native','rya1371-original-profile-replay-analysis-v1':'analysis','rya1371-original-dispatch-audit-v1':'analysis'})
    if a.vdw_wavelength_supplement:groups={'rya1371-oxygen-no-vdw-v1':'native','rya1371-oxygen-no-vdw-analysis-v1':'analysis','rya1371-hop-no-vdw-transfer-v1':'analysis','rya1371-hop-no-vdw-comparison-v1':'analysis','rya1371-vdw-profile-audit-v1':'analysis','rya1371-native-wavelength-roundtrip-v1':'receipt-only','rya1371-native-wavelength-roundtrip-v2':'kernel','rya1371-oxygen-transfer-residuals-v1':'analysis'}
    if a.native_vacuum_centres_supplement:groups={'rya1371-hop-native-vacuum-centres-v1':'analysis','rya1371-hop-native-vacuum-centres-comparison-v1':'analysis','rya1371-native-vacuum-centres-audit-v1':'analysis','rya1371-oxygen-transfer-residuals-v2':'analysis'}
    if a.multi_width_supplement:
        from config.constants import require_codex_path
        extra_sources['width_oracle_build']=require_codex_path('engines.multi_width_oracle');groups={'width_oracle_build':'native','rya1371-multi-width-comparison-v1':'analysis','rya1371-multi-width-laws-v1':'receipt-only','rya1371-multi-width-laws-v2':'kernel'}
    if a.original_source_widths_supplement:
        if any([a.output_opacity_supplement,a.boundary_eos_supplement,a.molecular_preflight_supplement,a.hop_provider_supplement,a.hop_trace_supplement,a.hop_fixed_transfer_supplement,a.hop_formal_stark_supplement,a.hop_helium_supplement,a.hop_profile_source_supplement,a.hop_table_profile_supplement,a.hop_atmospheric_profile_supplement,a.accurate_native_profile_supplement,a.vdw_wavelength_supplement,a.native_vacuum_centres_supplement,a.multi_width_supplement]):raise ValueError('Source widths archive cannot combine supplements')
        groups={name:'analysis'for name in ['rya1371-multi-width-law-export-v1','rya1371-hop-original-source-widths-dry-v1','rya1371-hop-original-source-widths-v1','rya1371-hop-original-source-widths-comparison-v1','rya1371-source-width-backend-audit-v1','rya1371-radiative-width-audit-v1']}
    a.output.mkdir(parents=True,exist_ok=False);records=[]
    for name,kind in groups.items():
        source=extra_sources.get(name,a.raw_root/name)
        if not source.is_dir():raise FileNotFoundError(source)
        for path in sorted(source.rglob('*')):
            rel=path.relative_to(source)
            if not path.is_file()or path.is_symlink()or 'input'in rel.parts:continue
            receipt=fingerprint(path);copy=kind=='analysis'or kind!='receipt-only'and len(rel.parts)==1 and path.suffix in ['.json','.log','.nml','.txt']
            if kind=='kernel':copy=len(rel.parts)==1 and path.suffix in ['.json','.log','.npz']
            if receipt['bytes']>5*1024*1024:copy=False
            target=a.output/name/rel
            if copy:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
                if fingerprint(target)!=receipt or fingerprint(path)!=receipt:raise ValueError('Archive byte identity failed')
            records.append(dict(raw=str(path),receipt=receipt,archived=str(target.relative_to(a.output))if copy else None,disposition=kind))
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',files=records,limitations=['Invalid single-ray v1 arrays remain external and inadmissible; receipts do not validate them.','Native binary/raw output remains external with byte receipts; compact analysis copied exactly.','Older exploratory readers mislabeled native per-Angstrom intensity as per-Hz; only dense-analysis-v3 and later qualified readers used here.'])
    (a.output/'archive_receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(raw_files=len(records),copied=sum(x['archived']is not None for x in records))))


if __name__=='__main__':main()
