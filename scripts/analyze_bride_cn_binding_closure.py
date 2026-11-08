"""Separate legacy native density scaling and H-normalization conventions."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint,read_marcs_structure
 import lightweaver.constants as c
 from lightweaver.atomic_table import PeriodicTable
 plan_path=ROOT/'config/experiments/rya1281_cn_binding_closure.json';plan=json.loads(plan_path.read_text());report=json.loads((a.source/'report.json').read_text())
 if not report['arithmetic_passed']or not report['molecular_switch_state_identity']:raise ValueError('Binding input not admitted as a numerical inventory')
 data=np.genfromtxt(a.source/'state_binding.csv',delimiter=',',names=True);composition=json.loads((a.source/'composition.json').read_text());species=np.load(a.source/'implemented_species.npz')
 native=require_codex_path('engines.ispec')/'synthesizer/turbospectrum/source/eqmol_pe_lu.f'
 if native.read_text().count('      rho=rho*1.2123e-8/tt\n')!=1:raise ValueError('Changed legacy density equation')
 cn_dir=require_codex_path('experiments.bride_cn_coefficients');archive=ROOT/'data/results/rya1281/cn_coefficients_v1'
 for name in ['report.json','CN/states.csv']:
  if(cn_dir/name).read_bytes()!=(archive/name).read_bytes():raise ValueError('Changed original CN report/state')
 cn=np.genfromtxt(cn_dir/'CN/states.csv',delimiter=',',names=True)
 for name in ['temperature_K','electron_pressure_dyn_cm2','gas_pressure_dyn_cm2','rho_g_cm3']:
  if not np.array_equal(cn[name].reshape(3,56),np.tile(data[name],(3,1))):raise ValueError('H inventory not bound to original CN state '+name)
 model=require_codex_path('grids.marcs_standard')/report['model']['name']
 if fingerprint(model)!=report['model']:raise ValueError('Changed original model')
 s=read_marcs_structure(model);H=next(x for x in composition if x['Z']==1);press=species['pressure_dyn_cm2'];mass=species['mass_amu']
 # Atomic masses are exported in amu, while the native density equation uses
 # the literal factor above. Evaluate the explicit LW mass-unit convention
 # separately rather than silently replacing the native density calculation.
 amu_g=c.Amu*1000;legacy=float(np.float32(1.2123e-8));electron_mass_amu=5.4858e-4
 # Read the actual native electron-mass literal and default-REAL representation.
 if native.read_text().count('        atmass(99)=5.4858e-4')!=1:raise ValueError('Changed native electron mass convention')
 electron_mass_amu=float(np.float32(electron_mass_amu));weighted=(press*mass).sum(axis=1)+data['electron_pressure_dyn_cm2']*electron_mass_amu
 rho_mass=weighted*amu_g/(data['boltz_erg_K']*data['temperature_K'])
 Xnetwork=data['H_nuclei_pressure_dyn_cm2']*H['mean_mass_amu']/weighted;Hmass_LW=PeriodicTable['H'].mass
 density_factor=rho_mass/(s['rho_kg_m3']/1000);Xfactor=Xnetwork/s['hydrogen_mass_fraction'];Hmass_factor=Hmass_LW/H['mean_mass_amu'];reconstructed=density_factor*Xfactor*Hmass_factor
 decomposition_error=float(np.max(np.abs(reconstructed/data['network_nH_over_MARCS_LW']-1)))
 if decomposition_error>64*np.finfo(float).eps:raise ValueError('Hydrogen unit-factor decomposition failed')
 rows=[]
 for j,row in enumerate(data):rows.append(dict(depth_index=j+1,legacy_density_factor_over_LW_amu_per_k=float(legacy*row['boltz_erg_K']/amu_g),rho_mass_unit_g_cm3=float(rho_mass[j]),rho_mass_unit_over_model=float(density_factor[j]),network_H_mass_fraction=float(Xnetwork[j]),network_over_model_H_mass_fraction=float(Xfactor[j]),LW_over_native_H_mean_mass=float(Hmass_factor),reconstructed_nH_ratio=float(reconstructed[j]),native_nH_ratio=float(row['network_nH_over_MARCS_LW'])))
 a.output.mkdir(parents=True,exist_ok=False)
 with(a.output/'closure.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 result=dict(schema='codex.bride_cn_binding_closure/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),source_report=fingerprint(a.source/'report.json'),source_states=fingerprint(a.source/'state_binding.csv'),source_species=fingerprint(a.source/'implemented_species.npz'),source_composition=fingerprint(a.source/'composition.json'),original_cn_states=fingerprint(cn_dir/'CN/states.csv'),native_EOS_source=fingerprint(native),LW_constants=fingerprint(c.__file__),LW_Amu_kg=c.Amu,LW_KBoltzmann_J_K=c.KBoltzmann,LW_H_mean_mass_amu=Hmass_LW,native_H_mean_mass_amu=H['mean_mass_amu'],native_electron_mass_amu=electron_mass_amu,legacy_density_factor=legacy,original_CN_atmosphere_fields_exact=True,H_ratio_factorization_max_relative_difference=decomposition_error,summary={key:[min(r[key]for r in rows),max(r[key]for r in rows)]for key in rows[0]if key!='depth_index'},limitations=['LW atomic-mass-unit convention is an explicit comparative unit reference, not an adopted replacement for native constants or an independent EOS.','Source mass-density factor differs from amu/k by about0.8percent; no upstream physical bug or correction is adjudicated.','Factorization is dimensional bookkeeping on the same exported state; not independent chemical validation.','CN populations per gram already use atmosphere density; multiplying their existing absolute opacity by EOS/atmosphere density would introduce a new inconsistent rescaling.','No native packet/provider/model admission, RT/SE or abundance change.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['summary']));print(f'Read accepted native inventory, original CN/model and explicit constants; wrote closure.csv/report.json to {a.output}')

if __name__=='__main__':main()
