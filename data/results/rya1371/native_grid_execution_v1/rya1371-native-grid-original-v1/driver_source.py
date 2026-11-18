"""Prepare receipt-bound original/enriched native frequency-grid control plans."""
import argparse,copy,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    from pipeline.native_frequency_grid import load_native_frequency_grid
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline-output',type=Path,required=True);a=p.parse_args()
    qualified=ROOT/'config/experiments/rya1371_openmulti3d_oxygen_shared16.json';base=json.loads(qualified.read_text())
    source_key='experiments.openmulti3d_oxygen_shared16_run';source=require_codex_path(source_key)
    src_report=source/'report.json';src_wave=source/'native/save/all_lam.bin'
    doc=json.loads(src_report.read_text())
    if not doc['computational_completed']or doc['plan']!=fingerprint(qualified):raise ValueError('Wrong native baseline producer')
    expected=next(r['receipt']for r in doc['output_files']if r['path']=='native/save/all_lam.bin')
    if fingerprint(src_wave)!=expected:raise ValueError('Changed native baseline grid')
    old=np.fromfile(src_wave,dtype='<f8')
    if len(old)!=39841:raise ValueError('Unexpected original grid')
    a.baseline_output.mkdir(parents=True,exist_ok=False)
    grid=a.baseline_output/'frequency_grid_vac_A.txt';np.savetxt(grid,old,fmt='%.17e',header='Original native shared16 grid; vacuum Angstrom')
    producer=dict(generator=fingerprint(__file__),computational_completed=True,scientific_admission='HOLD',all_original_nodes_exact=True,proposed_points=len(old),grid=fingerprint(grid),inputs={str(src_report):fingerprint(src_report),str(src_wave):fingerprint(src_wave),str(qualified):fingerprint(qualified)},scope='Original-grid reproduction prerequisite; no new RT or enriched nodes.')
    (a.baseline_output/'report.json').write_text(json.dumps(producer,indent=2)+'\n')
    for kind,key,name in [('baseline','experiments.native_grid_original_control','baseline_plan.json'),('continuum-refinement','experiments.native_continuum_grid_proposal','refinement_plan.json')]:
        directory=require_codex_path(key);plan=copy.deepcopy(base)
        plan['native_frequency_grid_control']=dict(kind=kind,resource=key,producer=fingerprint(directory/'report.json'),grid=fingerprint(directory/'frequency_grid_vac_A.txt'),source_resource=source_key,source_report=fingerprint(src_report),source_wave=expected)
        plan['scope']=base['scope']+' Explicit retained-node native grid '+kind+' control, same atom/atmosphere/abundance/physics/limits. Original-grid reproduction required before continuum enrichment. No scientific admission.'
        load_native_frequency_grid(plan,resolver=require_codex_path)
        (a.baseline_output/name).write_text(json.dumps(plan,indent=2)+'\n')
    print(json.dumps(dict(computational_completed=True,original_points=len(old),plans=['baseline_plan.json','refinement_plan.json'])))


if __name__=='__main__':main()
