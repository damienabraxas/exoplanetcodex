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
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--shared-convergence',action='store_true');a=p.parse_args();checks={};cases={};metrics={};rows=[];areas=[]
    configurations=CASES if not a.shared_convergence else {f'shared{factor}':(f'rya1371-oxygen-shared{factor}-v1',f'rya1371-oxygen-shared{factor}-analysis-v1',f'rya1371-native-gross-shared{factor}-analysis-v1')for factor in [4,8,16]}
    base_key=next(iter(configurations))
    for key,(run_name,analysis_name,physical_name)in configurations.items():
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
        if key!=base_key:
            for name in state:np.testing.assert_array_equal(state[name],cases[base_key]['state'][name])
            for name in ['star','C']:np.testing.assert_array_equal(values[name],cases[base_key]['values'][name])
            np.testing.assert_array_equal(spectra['wave_vac_A'],cases[base_key]['spectra']['wave_vac_A'])
        for centre in [6158.17,6300.304,6363.776,7771.944,7774.166,7775.388,8446.4]:areas.append(dict(case=key,centre_air_A=centre,scope='absent-transition control'if centre==6363.776 else'finite theoretical O-only disk-centre intensity',area_mA=1000*window_integral(spectra['wave_air_A'],spectra['Inu_cgs']/spectra['continuum_Inu_cgs'],centre-.5,centre+.5,max_gap=.0051)))
        cases[key]=dict(run=run,state=state,values=values,spectra=spectra)
        metrics[key]=dict(seconds=doc['seconds'],physical_pre_update_residual=pr['max_physical_residual'],observer_original_binaries_identical=pr['original_binary_files_bit_identical'],native_final_iteration=pr['final_stage'][0]['iteration'],nq_fac=doc['inputs']['controls']['nq_fac'],decouple_continuum=doc['inputs']['controls'].get('decouple_continuum',True),pre_update_physical_residual_within_original_rtol=bool(pr['max_physical_residual']<=doc['inputs']['controls']['rtol']))
    comparisons=[('NQ4','NQ2'),('shared4','NQ4')]if not a.shared_convergence else [('shared8','shared4'),('shared16','shared8')]
    for new,old in comparisons:
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
    probes=[('NQ4','rya1371-oxygen-uv-ir-output-v1'),('NQ4','rya1371-oxygen-catalog-edges-v1')]if not a.shared_convergence else [(f'shared{factor}',f'rya1371-oxygen-shared{factor}-uv-ir-v1')for factor in [8,16]]
    for source_key,folder in probes:
        directory=a.raw_root/folder;doc=json.loads((directory/'report.json').read_text())
        if not doc['computational_completed']:raise ValueError('Incomplete multiband output')
        checks[str(directory/'report.json')]=fingerprint(directory/'report.json');count=0
        for q in sorted((cases[source_key]['run']/'native/save').glob('*.bin')):
            if q.name=='out_lam.bin'or q.name.startswith(('ie_','flux_','Jnu_out_')):continue
            counterpart=directory/'native/save'/q.name
            if q.read_bytes()!=counterpart.read_bytes():raise ValueError('Output-only control changed SE '+q.name)
            checks[str(counterpart)]=fingerprint(counterpart);count+=1
        output_proofs[folder]=dict(original_internal_binary_files_identical=count,excluded_output_labels=['out_lam.bin','ie_*','flux_*','Jnu_out_*'],scope='All original non-output native binary files, including populations/rates/internal radiation, compared exactly. Flux-labelled files are output-window moments and are not certified disk-integrated flux under the disk-centre-only quadrature. Sparse output windows only; no full-band science validation.')
    multiband=[]
    if a.shared_convergence:
        probe_values={}
        for factor in [8,16]:
            folder=a.raw_root/f'rya1371-oxygen-shared{factor}-uv-ir-analysis-v1'
            doc=json.loads((folder/'report.json').read_text())
            if doc['native_report']!=fingerprint(a.raw_root/f'rya1371-oxygen-shared{factor}-uv-ir-v1/report.json'):raise ValueError('Wrong multiband analysis producer')
            for name in ['report.json','spectra.npz']:checks[str(folder/name)]=fingerprint(folder/name)
            with np.load(folder/'spectra.npz')as d:probe_values[factor]={key:d[key].copy()for key in ['wave_vac_A','Inu_cgs','continuum_Inu_cgs','air_frame_supported']}
        coarse=probe_values[8];fine=probe_values[16]
        np.testing.assert_array_equal(coarse['wave_vac_A'],fine['wave_vac_A']);np.testing.assert_array_equal(coarse['air_frame_supported'],fine['air_frame_supported'])
        wave=coarse['wave_vac_A']/10.;coarse_norm=coarse['Inu_cgs']/coarse['continuum_Inu_cgs'];fine_norm=fine['Inu_cgs']/fine['continuum_Inu_cgs']
        from pipeline.reference_engines import output_grid_vacuum_A
        plan=ROOT/'config/experiments/rya1371_openmulti3d_oxygen_shared16_uv_ir.json';probe_plan=json.loads(plan.read_text());checks[str(plan)]=fingerprint(plan)
        expected,_=output_grid_vacuum_A(probe_plan['output_grid']);np.testing.assert_allclose(wave,expected/10.,rtol=0,atol=2e-12)
        for lo,hi in probe_plan['output_grid']['windows']:
            mask=(wave>=lo-1e-12)&(wave<=hi+1e-12)
            if np.count_nonzero(mask)<2 or not np.isfinite(coarse_norm[mask]).all()or not np.isfinite(fine_norm[mask]).all():raise ValueError('Incomplete numerical UVIR window')
            multiband.append(dict(vacuum_min_nm=lo,vacuum_max_nm=hi,points=int(np.count_nonzero(mask)),max_absolute_normalized_intensity_change=float(np.max(abs(fine_norm[mask]-coarse_norm[mask]))),finite_window_area_delta_mA=1000*(window_integral(wave*10.,fine_norm,lo*10.,hi*10.,max_gap=.0051)-window_integral(wave*10.,coarse_norm,lo*10.,hi*10.,max_gap=.0051)),scope='Finite theoretical disk-centre output window; UV vacuum axis only; no observed/full-band/abundance validation'))
    a.output.mkdir(parents=True,exist_ok=False)
    for name,records in [('population_changes.csv',rows),('optical_areas.csv',areas)]+([('multiband_changes.csv',multiband)]if multiband else[]):
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
    if a.shared_convergence:
        import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        fig,axes=plt.subplots(len(multiband),2,figsize=(13,12),constrained_layout=True)
        for row,entry in zip(axes,multiband):
            lo=entry['vacuum_min_nm'];hi=entry['vacuum_max_nm'];mask=(wave>=lo-1e-12)&(wave<=hi+1e-12)
            row[0].plot(wave[mask],coarse_norm[mask],label='NQFAC8');row[0].plot(wave[mask],fine_norm[mask],'--',label='NQFAC16');row[0].legend(fontsize=7);row[0].set_ylabel('I / continuum')
            row[1].plot(wave[mask],fine_norm[mask]-coarse_norm[mask]);row[1].axhline(0,color='gray',lw=.7);row[1].set_ylabel('NQ16 minus NQ8')
            for ax in row:ax.set_xlabel('VACUUM wavelength / nm');ax.set_title(f'{lo:g}--{hi:g} nm');ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        fig.suptitle('Held O41/MARCS5750g4.5A8.70 shared-field frequency refinement\nSparse theoretical disk-centre probes; not observed/full-band validation')
        fig.savefig(a.output/'multiband_refinement.png',dpi=150);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed refinement audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=metrics,multiband_changes=multiband,output_only_proofs=output_proofs,scientific_admission='HOLD',computational_completed=True,scope='Native numerical/radiation convention sensitivity, not parameter fitting, full convergence or observed-solar validation. Rates evaluated at their own matched pre-update states, not forced coefficient parity. AIR convention only for inherited optical windows.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(metrics))
if __name__=='__main__':main()
