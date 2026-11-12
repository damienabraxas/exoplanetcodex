"""Compare native frequency/radiation controls and prove output-only UV/IR identity."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
CASES={'NQ2':('rya1371-oxygen-accurate-profile-v1','rya1371-oxygen-accurate-profile-analysis-v1','rya1371-native-gross-analysis-v1'),'NQ4':('rya1371-oxygen-frequency4-accurate-v1','rya1371-oxygen-frequency4-accurate-analysis-v3','rya1371-native-gross-frequency4-analysis-v1'),'shared4':('rya1371-oxygen-shared4-v1','rya1371-oxygen-shared4-analysis-v1','rya1371-native-gross-shared4-analysis-v1')}


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.bride_observables import window_integral
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();checks={};cases={};metrics={};rows=[];areas=[]
    for key,(run_name,analysis_name,physical_name)in CASES.items():
        run=a.raw_root/run_name;analysis=a.raw_root/analysis_name;physical=a.raw_root/physical_name;doc=json.loads((run/'report.json').read_text());ad=json.loads((analysis/'report.json').read_text());pr=json.loads((physical/'report.json').read_text())
        if not doc['computational_completed']or not pr['computational_completed']or ad['native_report']!=fingerprint(run/'report.json')or pr['inputs'].get(str(run/'report.json'))!=fingerprint(run/'report.json'):raise ValueError('Wrong qualified native producer')
        for q in [run/'report.json',analysis/'report.json',analysis/'arrays.npz',analysis/'spectra.npz',physical/'report.json',physical/'rates.npz']:checks[str(q)]=fingerprint(q)
        for entry in doc['output_files']:
            q=run/entry['path']
            if fingerprint(q)!=entry['receipt']:raise ValueError('Changed native output')
            checks[str(q)]=entry['receipt']
        with np.load(analysis/'arrays.npz')as d:state={name:d[name].copy()for name in ['height_m','consumed_state_cgs','native_nstar_cm3','source_logtau500_on_native_height']}
        with np.load(physical/'rates.npz')as d:values={name:d[name].copy()for name in ['R','C','n','star','post']}
        with np.load(analysis/'spectra.npz')as d:spectra={name:d[name].copy()for name in ['wave_vac_A','wave_air_A','Inu_cgs','continuum_Inu_cgs']}
        if key!='NQ2':
            for name in state:np.testing.assert_array_equal(state[name],cases['NQ2']['state'][name])
            for name in ['star','C']:np.testing.assert_array_equal(values[name],cases['NQ2']['values'][name])
            np.testing.assert_array_equal(spectra['wave_vac_A'],cases['NQ2']['spectra']['wave_vac_A'])
        for centre in [6158.17,6300.304,6363.776,7771.944,7774.166,7775.388,8446.4]:areas.append(dict(case=key,centre_air_A=centre,scope='absent-transition control'if centre==6363.776 else'finite theoretical O-only disk-centre intensity',area_mA=1000*window_integral(spectra['wave_air_A'],spectra['Inu_cgs']/spectra['continuum_Inu_cgs'],centre-.5,centre+.5,max_gap=.0051)))
        cases[key]=dict(run=run,state=state,values=values,spectra=spectra)
        metrics[key]=dict(seconds=doc['seconds'],physical_pre_update_residual=pr['max_physical_residual'],observer_original_binaries_identical=pr['original_binary_files_bit_identical'],native_final_iteration=pr['final_stage'][0]['iteration'],nq_fac=doc['inputs']['controls']['nq_fac'],decouple_continuum=doc['inputs']['controls'].get('decouple_continuum',True))
    for new,old in [('NQ4','NQ2'),('shared4','NQ4')]:
        ratio=cases[new]['values']['post']/cases[old]['values']['post'];idx=np.unravel_index(np.argmax(abs(ratio-1)),ratio.shape);delta_areas=[x['area_mA']-next(y['area_mA']for y in areas if y['case']==old and y['centre_air_A']==x['centre_air_A'])for x in areas if x['case']==new and x['centre_air_A']!=6363.776]
        metrics[new]['relative_to_'+old]=dict(population_ratio_min=float(ratio.min()),population_ratio_max=float(ratio.max()),max_relative_change=float(abs(ratio[idx]-1)),worst_level_1based=int(idx[0]+1),worst_depth_bottom_to_top=int(idx[1]),max_absolute_six_window_area_difference_mA=float(np.max(abs(np.array(delta_areas)))))
        for k in range(41):
            for z in range(56):rows.append(dict(case=new,reference=old,level_1based=k+1,depth_bottom_to_top=z,population_ratio=float(ratio[k,z])))
    own=require_codex_path('experiments.hop_ghost_coupled_oxygen')
    with np.load(own/'populations.npz')as d:own_n=d['populations'].copy();own_star=d['nstar'].copy()
    checks[str(own/'populations.npz')]=fingerprint(own/'populations.npz');checks[str(own/'report.json')]=fingerprint(own/'report.json')
    for key in cases:
        native=cases[key]['values'];np.testing.assert_array_equal(own_star,native['star'][:,::-1]*1e6);ratio=own_n/(native['post'][:,::-1]*1e6);metrics[key]['previous_own_coupled_max_population_relative_difference']=float(np.max(abs(ratio-1)))
    output_proofs={}
    for folder in ['rya1371-oxygen-uv-ir-output-v1','rya1371-oxygen-catalog-edges-v1']:
        directory=a.raw_root/folder;doc=json.loads((directory/'report.json').read_text())
        if not doc['computational_completed']:raise ValueError('Incomplete multiband output')
        checks[str(directory/'report.json')]=fingerprint(directory/'report.json');count=0
        for q in sorted((cases['NQ4']['run']/'native/save').glob('*.bin')):
            if q.name=='out_lam.bin'or q.name.startswith(('ie_','flux_','Jnu_out_')):continue
            counterpart=directory/'native/save'/q.name
            if q.read_bytes()!=counterpart.read_bytes():raise ValueError('Output-only control changed SE '+q.name)
            checks[str(counterpart)]=fingerprint(counterpart);count+=1
        output_proofs[folder]=dict(original_internal_binary_files_identical=count,excluded_output_labels=['out_lam.bin','ie_*','flux_*','Jnu_out_*'],scope='All original non-output native binary files, including populations/rates/internal radiation, compared exactly. Flux-labelled files are output-window moments and are not certified disk-integrated flux under the disk-centre-only quadrature. Sparse output windows only; no full-band science validation.')
    a.output.mkdir(parents=True,exist_ok=False)
    for name,records in [('population_changes.csv',rows),('optical_areas.csv',areas)]:
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed refinement audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=metrics,output_only_proofs=output_proofs,scientific_admission='HOLD',computational_completed=True,scope='Native numerical/radiation convention sensitivity, not parameter fitting, full convergence or observed-solar validation. Rates evaluated at their own matched pre-update states, not forced coefficient parity. AIR convention only for inherited optical windows.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(metrics))
if __name__=='__main__':main()
