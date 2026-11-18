"""Audit native O rates/consumed mesh and compare independent departures by height."""
import argparse,csv,importlib.util,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.reference_engines import load_plan_engine,fingerprint
    from pipeline.openmulti3d_output import read_rate_pair,read_patch_set,interpolate_positive,per_angstrom_to_per_hz
    from pipeline.bride_rate_audit import process_flow
    from pipeline.wavelength_util import vac_to_air
    run=json.loads((a.source/'report.json').read_text());_,app,_,engine=load_plan_engine(dict(engine=run['engine'],profile_variant=run.get('profile_variant','original')));checks=[]
    if not run['computational_completed']or run['engine_receipt']!=engine:raise ValueError('Unqualified execution/pin')
    for entry in run['output_files']:
        p=a.source/entry['path']
        if fingerprint(p)!=entry['receipt']:raise ValueError('Changed native output')
        checks.append((p,entry['receipt']))
    native=a.source/'native';save=native/'save';parts=[];rates=[]
    for path in sorted(save.glob('rates_O_*.bin')):
        data=read_rate_pair(path,save/path.name.replace('rates_','colrates_'));n=data['populations_cm3'];s=data['nstar_cm3'];net,gross=process_flow(data['collision_rates_s_inv']+data['radiative_rates_s_inv'],n);res=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));cons=float(np.max(abs(n.sum(axis=0)/s.sum(axis=0)-1)));parts.append(data);rates.append(dict(file=path.name,iteration=data['iteration'],rate_residual=res,conservation=cons))
    if len(parts)!=8:raise ValueError('Expected all eight declared1D depth slabs')
    n=np.concatenate([d['populations_cm3']for d in parts],axis=1);s=np.concatenate([d['nstar_cm3']for d in parts],axis=1);b=n/s
    if n.shape!=(41,56):raise ValueError('Changed native mesh layout')
    ap=read_patch_set(save,'atmos');atmos=np.concatenate([ap[k].reshape(7,14,order='F')for k in sorted(ap)],axis=0)
    mesh=(native/'atmos_mesh.txt').read_text().splitlines();dz=np.concatenate([np.full(int(mesh[i].split()[2]),float(mesh[i+1].split()[2]))for i in range(0,len(mesh),2)])
    # read_atmos_mod%interpolate uses newz=lower+i*ds, not cell midpoints.
    height=np.cumsum(dz)*1e6
    reference_path=require_codex_path('experiments.bride_oxygen_reference');reference=json.loads(reference_path.read_text());source=reference['atmosphere'];source_h=np.array(source['height_m'])[::-1];source_h-=source_h[0]
    foreign=np.array(reference['reference']['departures'])[:,::-1];external,interp=interpolate_positive(source_h,foreign,height)
    bride_dir=ROOT/'data/results/rya1281/solar_full_cap_v1/rya1281-solar-full-cap-v1/baseline';br=json.loads((bride_dir/'report.json').read_text())
    if br['atom']!=run['inputs']['atom']or reference['atom']!=br['atom']:raise ValueError('Source atom identity mismatch')
    with np.load(bride_dir/'populations.npz')as d:bb=d['populations']/d['nstar'];bride,bi=interpolate_positive(source_h,bb[:,::-1],height)
    local_tau=np.interp(np.clip(height,source_h[0],source_h[-1]),source_h,np.log10(reference['reference']['tau'])[::-1]);diagnostic=(local_tau>=-3)&(local_tau<=0)
    rows=[]
    for i in range(41):
        for k,z in enumerate(height):rows.append(dict(level_1based=i+1,depth_bottom_to_top=k,height_above_source_bottom_m=z,source_interpolated_logtau500=local_tau[k],native_b=b[i,k],held_multi_b_log_height=external[i,k],bride_b_log_height=bride[i,k],native_minus_multi=b[i,k]-external[i,k],native_minus_bride=b[i,k]-bride[i,k]))
    a.output.mkdir(parents=True,exist_ok=False)
    for name,rr in [('rates.csv',rates),('populations.csv',rows)]:
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rr[0]),lineterminator='\n');w.writeheader();w.writerows(rr)
    np.savez_compressed(a.output/'arrays.npz',height_m=height,native_populations_cm3=n,native_nstar_cm3=s,native_b=b,held_multi_b_on_native_height=external,bride_b_on_native_height=bride,consumed_state_cgs=atmos,source_logtau500_on_native_height=local_tau)
    spectra=False
    if(save/'out_lam.bin').exists():
        # The current engine writes an explicit set manifest; require it rather
        # than inferring physics from a three-component array or importing its GUI stack.
        labels={}
        for line in (native/'sets.txt').read_text().splitlines():
            if not line.strip()or line.lstrip().startswith('#'):continue
            number,kind,element,offset=line.split()
            if element!='-'or offset!='-'or kind in labels:raise ValueError('Unexpected abundance-offset/output set')
            labels[kind]=int(number)-1
        if set(labels)!={'nlte','lte','continuum'}or set(labels.values())!={0,1,2}:raise ValueError('Incomplete spectrum set identity')
        ie=read_patch_set(save,'ie');positive=sorted(k for k in ie if k[1]>0);combined=np.concatenate([ie[k]for k in positive],axis=0)
        if combined.ndim!=5 or combined.shape[1:4]!=(1,1,1)or combined.shape[4]!=3:raise ValueError('Not the declared disk-centre intensity layout')
        I=combined[:,0,0,0,labels['nlte']];Ilte=combined[:,0,0,0,labels['lte']];Icont=combined[:,0,0,0,labels['continuum']]
        vac=np.fromfile(save/'out_lam.bin',dtype='<f8');air=vac_to_air(vac)
        if len(air)!=len(I)or np.any(Icont<=0):raise ValueError('Output wavelength/intensity layout mismatch')
        # rt_patch_mod.f90::save_intensity writes per vacuum Angstrom,
        # not per Hz. Preserve its native values and apply the Jacobian.
        np.savez_compressed(a.output/'spectra.npz',wave_vac_A=vac,wave_air_A=air,Ilambda_cgs_per_A=np.asarray(I).squeeze(),LTE_Ilambda_cgs_per_A=np.asarray(Ilte).squeeze(),continuum_Ilambda_cgs_per_A=np.asarray(Icont).squeeze(),Inu_cgs=per_angstrom_to_per_hz(vac,np.asarray(I).squeeze()),LTE_Inu_cgs=per_angstrom_to_per_hz(vac,np.asarray(Ilte).squeeze()),continuum_Inu_cgs=per_angstrom_to_per_hz(vac,np.asarray(Icont).squeeze()));spectra=True
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,i in zip(axes.flat,[5,7,8,9]):
        for values,label in [(b,'OpenMULTI3D'),(external,'Held MULTI / height interpolation'),(bride,'Bride / height interpolation')]:ax.plot(local_tau,values[i],label=label)
        ax.set_xlabel('Source-interpolated log tau500');ax.set_ylabel('b');ax.set_title('O level '+str(i+1));ax.legend(fontsize=7)
    fig.suptitle('Same source atom/atmosphere; consumed grid/EOS/opacity/collision differences remain');fig.tight_layout();fig.savefig(a.output/'departures.png',dpi=150);plt.close(fig)
    for p,fp in checks:
        if fingerprint(p)!=fp:raise ValueError('Native output changed during analysis')
    doc=dict(schema='codex.openmulti3d_oxygen_audit/1',scientific_admission='HOLD',profile_variant=run.get('profile_variant','original'),generator=fingerprint(__file__),native_report=fingerprint(a.source/'report.json'),engine_receipt=engine,rates=rates,stored_rate_max_residual=max(x['rate_residual']for x in rates),max_conservation=max(x['conservation']for x in rates),distinct_final_iterations=sorted(set(x['iteration']for x in rates)),reference=fingerprint(reference_path),bride=fingerprint(bride_dir/'report.json'),height_interpolation=interp,bride_height_interpolation=bi,triplet_diagnostic_native_vs_multi_max_abs=float(np.max(abs(b[[5,7,8,9]][:,diagnostic]-external[[5,7,8,9]][:,diagnostic]))),triplet_diagnostic_native_vs_bride_max_abs=float(np.max(abs(b[[5,7,8,9]][:,diagnostic]-bride[[5,7,8,9]][:,diagnostic]))),spectra=spectra,limitations=['Stored physical-rate balance is independently reconstructed, not a fresh fixed-population radiation replay.','Native grid is pressure-graded and PCHIP/log interpolated; reference b uses declared positive linear height interpolation. No row-index identity claimed.','Native EOS/background opacity, collider populations, collision interpolation and top boundary differ; no identical-physics certification.','Diagnostic source logtau interval is not a contribution function; no abundance/production/model admission.'])
    (a.output/'report.json').write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n');print(json.dumps(doc))

if __name__=='__main__':main()
