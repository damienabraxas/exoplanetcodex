"""Source-based H-pair/Thomson closure; no fitted continuum coefficients."""
import argparse,csv,json,re,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def pair_tables(text):
    clean=re.sub(r'/\*.*?\*/','',text,flags=re.S);arrays=[]
    for name in ['lambdaFFPlus','tempFFPlus','kappaFFPlus']:
        found=re.findall(r'\b'+name+r'\s*\[[^\]]+\]\s*=\s*\{(.*?)\};',clean,flags=re.S)
        if len(found)!=1:raise ValueError('Ambiguous H2+ table '+name)
        arrays.append(np.array([float(x.strip())for x in found[0].split(',')if x.strip()]))
    wave,temp,values=arrays
    if min(len(wave),len(temp))<2 or len(values)!=len(wave)*len(temp)or not all(np.isfinite(x).all()for x in arrays)or np.any(np.diff(wave)<=0)or np.any(np.diff(temp)<=0)or np.any(values<0):raise ValueError('Invalid H-pair table dimensions/domain')
    return wave,temp,values.reshape(len(wave),len(temp))


def interpolate_pair(wave,temp,tables):
    L,T,Q=tables;wave=np.asarray(wave);temp=np.asarray(temp)
    if np.any(wave<L[0])or np.any(wave>=L[-1])or np.any(temp<T[0])or np.any(temp>T[-1]):raise ValueError('Pair table domain requires explicit no-extrapolation')
    i=np.clip(np.searchsorted(L,wave,side='right')-1,0,len(L)-2);j=np.clip(np.searchsorted(T,temp,side='right')-1,0,len(T)-2);f=(wave-L[i])/(L[i+1]-L[i]);g=(temp-T[j])/(T[j+1]-T[j])
    return (1-f[:,None])*(1-g[None,:])*Q[i[:,None],j[None,:]]+f[:,None]*(1-g[None,:])*Q[i[:,None]+1,j[None,:]]+(1-f[:,None])*g[None,:]*Q[i[:,None],j[None,:]+1]+f[:,None]*g[None,:]*Q[i[:,None]+1,j[None,:]+1]


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 import lightweaver as lw
 import lightweaver.constants as c
 from lightweaver import rh_atoms
 import lightweaver.LwCompiled as compiled
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.wavelength_util import vac_to_air
 plan_path=ROOT/'config/experiments/rya1281_H_background_components.json';plan=json.loads(plan_path.read_text());parent=json.loads((a.source/'report.json').read_text())
 if fingerprint(compiled.__file__)!=parent['backend']or fingerprint(c.__file__)!=parent['LW_constants']:raise ValueError('Changed paired-control backend/constants')
 arrays=[];inputs=[]
 for case in parent['cases']:
  folder=a.source/case['case']
  for name,key in [('input.json','input'),('background_window.npz','background_window')]:
   if fingerprint(folder/name)!=case[key]:raise ValueError('Changed paired background input')
  inputs.append(json.loads((folder/'input.json').read_text()));raw=np.load(folder/'background_window.npz');arrays.append({k:raw[k].copy()for k in raw.files});raw.close()
 first,second=arrays;A,B=inputs;wave=first['wavelength_vac_nm'];T=np.array(A['temperature_K']);ne=np.array(A['ne_m3']);H=np.array(A['nH_m3']);r=np.array(B['nH_m3'])/H
 if A['temperature_K']!=B['temperature_K']or A['ne_m3']!=B['ne_m3']or not np.array_equal(wave,second['wavelength_vac_nm']):raise ValueError('Component-control fixed fields differ')
 source=require_codex_path('engines.bride_eval')/'lw_src/Source/Background.cpp';original=source.read_bytes();tables=pair_tables(original.decode());K=interpolate_pair(wave,T,tables);pops=lw.lte_pops(rh_atoms.H_6_atom(),T,ne,H,debye=True)
 pair=(pops[0][None,:]*1e-29)*(pops[-1][None,:]*1e-20)*K
 sigma=8*np.pi/3*(c.QElectron/(np.sqrt(4*np.pi*c.Epsilon0)*(np.sqrt(c.MElectron)*c.CLight)))**4;thomson=ne*sigma;Bnu=lw.planck(T[None,:],wave[:,None])
 prediction={'chi':r[None,:]*first['chi']+(r*r-r)[None,:]*pair+(1-r)[None,:]*thomson[None,:],'eta':r[None,:]*first['eta']+(r*r-r)[None,:]*pair*Bnu,'sca':r[None,:]*first['sca']+(1-r)[None,:]*thomson[None,:]}
 bound=float(plan['arithmetic_eps_multiple']*np.finfo(float).eps);errors={k:float(np.max(abs(v/second[k]-1)))for k,v in prediction.items()};rows=[]
 for j in range(56):
  row=dict(depth_index=j+1,H_ratio=float(r[j]),pair_fraction_of_baseline_chi_min=float(np.min(pair[:,j]/first['chi'][:,j])),pair_fraction_of_baseline_chi_max=float(np.max(pair[:,j]/first['chi'][:,j])),Thomson_fraction_of_baseline_sca_min=float(np.min(thomson[j]/first['sca'][:,j])),Thomson_fraction_of_baseline_sca_max=float(np.max(thomson[j]/first['sca'][:,j])))
  for key,value in prediction.items():row[key+'_max_relative_closure_error']=float(np.max(abs(value[:,j]/second[key][:,j]-1)))
  rows.append(row)
 a.output.mkdir(parents=True,exist_ok=False)
 with(a.output/'summary.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 # Match the parent's predeclared original AIR sample coordinates by their
 # recorded vacuum indices, without a roundtrip conversion in the physics.
 samples=np.load(a.source/'sampled_backgrounds.npz');indices=np.searchsorted(wave,samples['wavelength_vac_nm']);np.testing.assert_array_equal(wave[indices],samples['wavelength_vac_nm']);np.savez(a.output/'sampled_components.npz',wavelength_vac_nm=wave[indices],pair_chi_m_inv=pair[indices],Thomson_chi_m_inv=thomson,**{k+'_predicted':v[indices]for k,v in prediction.items()});samples.close()
 if source.read_bytes()!=original:raise ValueError('Inspected source changed')
 result=dict(schema='codex.bride_H_background_components/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),parent=fingerprint(a.source/'report.json'),inspected_source=fingerprint(source),backend=parent['backend'],constants=parent['LW_constants'],closure_max_relative_errors=errors,prospective_arithmetic_bound=bound,closure_passed=all(x<=bound for x in errors.values()),summary={k:[min(row[k]for row in rows),max(row[k]for row in rows)]for k in ['pair_fraction_of_baseline_chi_min','pair_fraction_of_baseline_chi_max','Thomson_fraction_of_baseline_sca_min','Thomson_fraction_of_baseline_sca_max']},limitations=['Source-table/runtime closure on fixed state, not independent physical continuum validation or source-to-binary rebuild attestation.','H0*H+ pair continuum is present despite empty LW molecular population table; no additional molecular equilibrium inferred.','No fitted amplitude, background/default change, formal solution or O NLTE solve.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['closure_max_relative_errors','closure_passed','summary']}))
 if not result['closure_passed']:raise RuntimeError('H-background component closure failed; no bound adjustment')

if __name__=='__main__':main()
