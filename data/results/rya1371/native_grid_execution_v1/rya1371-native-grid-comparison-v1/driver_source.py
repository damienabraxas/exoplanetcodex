"""Compare native retained-node grid solves with exact physical-input checks."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from config.constants import require_codex_path
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();checks={str(Path(__file__)):fingerprint(__file__)}
    groups=[('baseline','experiments.openmulti3d_oxygen_shared16_physical_rates','experiments.openmulti3d_oxygen_shared16_analysis','experiments.openmulti3d_oxygen_shared16_run'),('refined','experiments.openmulti3d_oxygen_grid_refined_physical_rates','experiments.openmulti3d_oxygen_grid_refined_analysis','experiments.openmulti3d_oxygen_grid_refined_run')]
    cases={};spectra={};metrics={};grids={}
    for name,physical,analysis,run in groups:
        physical,analysis,run=[require_codex_path(key)for key in [physical,analysis,run]]
        for path in [physical/'report.json',physical/'rates.npz',analysis/'report.json',analysis/'spectra.npz',run/'report.json',run/'native/save/all_lam.bin']:checks[str(path)]=fingerprint(path)
        report=json.loads((physical/'report.json').read_text());run_doc=json.loads((run/'report.json').read_text());an=json.loads((analysis/'report.json').read_text())
        if not report['computational_completed']or not run_doc['computational_completed']or an['native_report']!=fingerprint(run/'report.json'):raise ValueError('Unqualified native producer')
        if report['inputs'].get(str(run/'report.json'))!=fingerprint(run/'report.json'):raise ValueError('Observer baseline mismatch')
        with np.load(physical/'rates.npz')as x:cases[name]={key:x[key].copy()for key in x.files}
        with np.load(analysis/'spectra.npz')as x:spectra[name]={key:x[key].copy()for key in x.files}
        grids[name]=np.fromfile(run/'native/save/all_lam.bin',dtype='<f8')
        q=cases[name]['R']+cases[name]['C'];n=cases[name]['n'];inc=np.sum(q*n[None,:,:],axis=1);out=q.sum(axis=0)*n
        residual=float(np.max(abs(inc-out)/(inc+out)))
        np.testing.assert_allclose(residual,report['max_physical_residual'],rtol=1e-9,atol=1e-15)
        metrics[name]=dict(points=len(grids[name]),max_pre_update_physical_residual=residual,max_update_relative=report['max_update_relative'],observer_unchanged_original_binaries=report['original_binary_files_bit_identical'])
    old,new=cases['baseline'],cases['refined'];idx=np.searchsorted(grids['refined'],grids['baseline'])
    np.testing.assert_array_equal(grids['refined'][idx],grids['baseline'])
    for key in ['C','star']:np.testing.assert_array_equal(new[key],old[key])
    ratios={key:new[key]/old[key]-1 for key in ['n','post']};rows=[]
    for key in ratios:
        loc=np.unravel_index(np.argmax(abs(ratios[key])),ratios[key].shape)
        metrics[key+'_mesh_change']=dict(max_relative=float(abs(ratios[key][loc])),worst_level_1based=int(loc[0]+1),worst_depth_top_to_bottom=int(55-loc[1]),triplet_levels6_8_9_10_max_relative=float(np.max(abs(ratios[key][[5,7,8,9]]))))
        for level in range(41):
            for depth in range(56):rows.append(dict(state=key,level_1based=level+1,depth_top_to_bottom=55-depth,relative_population_change=float(ratios[key][level,depth])))
    s0,s1=spectra['baseline'],spectra['refined']
    for key in ['wave_vac_A','wave_air_A','air_frame_supported']:np.testing.assert_array_equal(s0[key],s1[key])
    wave=s0['wave_air_A'];b=s0['Inu_cgs']/s0['continuum_Inu_cgs'];r=s1['Inu_cgs']/s1['continuum_Inu_cgs'];windows=[]
    plan=ROOT/'config/experiments/rya1371_openmulti3d_oxygen_shared16.json';checks[str(plan)]=fingerprint(plan);plan_doc=json.loads(plan.read_text())
    for lo,hi in plan_doc['output_grid']['windows']:
        mask=(wave>=lo-1e-8)&(wave<=hi+1e-8)
        if mask.sum()<2 or np.any(np.diff(wave[mask])>.0051):raise ValueError('Incomplete fixed optical window')
        windows.append(dict(air_min_A=lo,air_max_A=hi,points=int(mask.sum()),max_absolute_normalized_intensity_change=float(np.max(abs(r[mask]-b[mask]))),finite_window_area_change_mA=float(-1000*np.trapezoid(r[mask]-b[mask],wave[mask]))))
    metrics['spectra']=dict(max_absolute_normalized_intensity_change=float(np.max(abs(r-b))),max_absolute_finite_window_area_change_mA=max(abs(w['finite_window_area_change_mA'])for w in windows))
    a.output.mkdir(parents=True,exist_ok=False)
    for name,values in [('populations.csv',rows),('windows.csv',windows)]:
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(values[0]),lineterminator='\n');w.writeheader();w.writerows(values)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    im=axes[0].imshow(100*ratios['post'][:,::-1],aspect='auto',origin='lower',extent=[-.5,55.5,.5,41.5],cmap='RdBu_r');fig.colorbar(im,ax=axes[0],label='Stored population change (%)')
    axes[0].set(xlabel='Native depth index (top to bottom)',ylabel='Source level (1-based)',title='Original 39,841 to enriched 59,003 grid')
    mask=(wave>=7770)&(wave<=7777);axes[1].plot(wave[mask],1e6*(r[mask]-b[mask]));axes[1].set(xlabel='AIR wavelength (Angstrom)',ylabel='Normalized intensity change (1e-6)',title='Theoretical disk-centre O I triplet window')
    fig.suptitle('Held 5750K / logg4.5 / A(O)8.70 — numerical control; science HOLD');fig.savefig(a.output/'grid_comparison.png',dpi=160);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed grid comparison input')
    report=dict(generator=fingerprint(__file__),inputs=checks,computational_completed=True,scientific_admission='HOLD',metrics=metrics,windows=windows,collision_and_nstar_exact=True,all_original_nodes_retained=True,scope='Two self-consistent native numerical-grid solves with identical atom/atmosphere/physics and exact C/nStar. Physical rates are independently reconstructed at each exact PRE-update state; stored post-update populations and emergent optical windows compared separately. Changed R includes field/population/BB+BF sampling feedback, not a fixed-population quadrature-only attribution. Finite AIR-window areas are not full equivalent widths or observed abundance errors. No canonical-solar/flux/fullband/3D/model admission.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(metrics))


if __name__=='__main__':main()
