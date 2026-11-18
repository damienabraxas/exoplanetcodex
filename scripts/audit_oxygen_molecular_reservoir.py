"""Native molecular-equilibrium intervention; atomic reservoirs versus selected molecules."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.openmulti3d_output import read_patch_set
    from pipeline.bride_observables import window_integral
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline-run','control-run','baseline-analysis','control-analysis','original-analysis','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();inputs={};loaded=[];species=[]
    def load(path):
        inputs[str(path)]=fingerprint(path)
        with np.load(path)as z:return {k:z[k].copy()for k in z.files}
    runs=[json.loads((folder/'report.json').read_text())for folder in [a.baseline_run,a.control_run]]
    for r in runs:
        if not r['computational_completed']:raise ValueError('Incomplete native intervention')
    c=[dict(r['inputs']['controls'])for r in runs]
    if c[0].pop('molecular_equilibrium')!='native'or c[1].pop('molecular_equilibrium')!='suppressed-above-1K-control'or c[0]!=c[1]:raise ValueError('Not the declared one-axis chemistry control')
    for key in ['atom','atmosphere']:
        if runs[0]['inputs'][key]!=runs[1]['inputs'][key]:raise ValueError('Physical source identity differs')
    for folder,run in zip([a.baseline_analysis,a.control_analysis],[a.baseline_run,a.control_run]):
        values=load(folder/'arrays.npz');values.update(load(folder/'spectra.npz'));loaded.append(values);sp={}
        for key in ['O_I','O_II','CO','OH','H2O']:
            save=run/'native/save';parts=read_patch_set(save,key)
            sp[key]=np.concatenate([parts[k].ravel(order='F')for k in sorted(parts)]).astype(float)
            for path in save.glob(key+'*'):inputs[str(path)]=fingerprint(path)
        species.append(sp)
    b,d=loaded;original=load(a.original_analysis/'arrays.npz')
    np.testing.assert_array_equal(b['native_populations_cm3'],original['native_populations_cm3'])
    np.testing.assert_array_equal(b['consumed_state_cgs'][:,:3],d['consumed_state_cgs'][:,:3])
    np.testing.assert_array_equal(b['height_m'],d['height_m']);np.testing.assert_array_equal(b['wave_vac_A'],d['wave_vac_A'])
    s,t=species;atomic=s['O_I']+s['O_II'];selected=s['CO']+s['OH']+s['H2O'];fraction=selected/(atomic+selected)
    rows=[dict(depth_bottom_to_top=k,height_m=b['height_m'][k],temperature_K=b['consumed_state_cgs'][k,2],atomic_O_cm3=atomic[k],CO_cm3=s['CO'][k],OH_cm3=s['OH'][k],H2O_cm3=s['H2O'][k],selected_molecule_O_fraction=fraction[k],cutoff_atomic_O_cm3=t['O_I'][k]+t['O_II'][k],cutoff_CO_cm3=t['CO'][k],cutoff_OH_cm3=t['OH'][k],cutoff_H2O_cm3=t['H2O'][k])for k in range(len(atomic))]
    areas=[];wave=b['wave_air_A']
    for centre in [6158.17,6300.304,7771.944,7774.166,7775.388,8446.4]:
        vb=window_integral(wave,b['Inu_cgs']/b['continuum_Inu_cgs'],centre-.5,centre+.5,max_gap=.0051)*1000;vd=window_integral(wave,d['Inu_cgs']/d['continuum_Inu_cgs'],centre-.5,centre+.5,max_gap=.0051)*1000
        areas.append(dict(centre_air_A=centre,baseline_area_mA=vb,cutoff_area_mA=vd,change_mA=vd-vb))
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'reservoir.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    for path,receipt in inputs.items():
        if fingerprint(path)!=receipt:raise ValueError('Changed chemistry input')
    report=dict(generator=fingerprint(__file__),scientific_admission='HOLD',inputs=inputs,export_replay_populations_bit_identical=True,consumed_T_ne_rho_bit_identical=True,selected_molecule_O_fraction_range=[float(fraction.min()),float(fraction.max())],max_cutoff_selected_molecule_density_cm3=float(max(np.max(t[key])for key in ['CO','OH','H2O'])),windows=areas,limitations=['CO/OH/H2O are selected O reservoirs, not a closed budget over every molecule and ion.','Cutoff control changes native molecular equilibrium and the resulting colliders/background/SE together; not a unique mechanism attribution or preferred model.','No validated molecular NLTE, observation fit, full atmosphere/EOS parity or production admission.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['selected_molecule_O_fraction_range','max_cutoff_selected_molecule_density_cm3','windows']}))


if __name__=='__main__':main()
