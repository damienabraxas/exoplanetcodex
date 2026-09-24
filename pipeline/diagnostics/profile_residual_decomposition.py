from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--engine-root',type=Path,required=True); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--region',default='nir_cn_iag'); ap.add_argument('--a-n',type=float,required=True); ap.add_argument('--c',type=float,default=8.46); ap.add_argument('--o',type=float,default=8.73)
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=False); sys.path.insert(0,str(a.engine_root.resolve()))
    from pipeline import cno_synthesis as c
    from pipeline.uncertainty_stack import params_and_deltas
    region=c.REGIONS[a.region]; diag=[d for d in c.REGION_DIAGNOSTICS[a.region] if d.element=='N' and d.role=='primary'][0]
    broad=c.preflight(region,'solar',(diag,)); params,_=params_and_deltas('solar'); llpath,_,gf=c.region_atomic_linelist(region,(diag,)); ll,iso,chem=c._load_synth_resources(str(llpath),apply_canonical_gf=(gf or {}).get('apply_canonical_gf',False)) if llpath else c._load_synth_resources(); sab=c.ispec.read_solar_abundances(c._ISPEC_SOLAR_ABUND_FILE); codes=c._atom_codes(('C','N','O','Ni'),chem,sab); ow,of=c._load_region_spectrum('solar',region); mask=np.zeros(len(ow),bool)
    for lo,hi in diag.windows_A: mask|=(ow*10>=lo)&(ow*10<=hi)
    ow,of=ow[mask],of[mask]; atm=c._load_atmosphere(params['teff_K'],params['logg'],params['feh'],params['vturb_kms']); state=c._seed_abundances('solar',params,codes,c._solar_A(('C','N','O','Ni'),chem,sab),params['feh']); state.update({'C':a.c,'O':a.o,'N':a.a_n})
    rows=[]; tmp=str(a.out/'synth'); Path(tmp).mkdir()
    for lo,hi in diag.windows_A:
        m=(ow*10>=lo)&(ow*10<=hi); ow_i,of_i=ow[m],of[m]; sw=np.arange(lo/10,hi/10+c._WSTEP_NM*.5,c._WSTEP_NM)
        sf=c._synth_window(sw,atm,params,ll,iso,sab,c._fixed_ab(state,codes),broad,True,tmp); model=np.interp(ow_i,sw,sf,left=1.0,right=1.0); r=of_i-model; z=r/c._SIGMA_FLUX
        rows.append({'window_air_A':[lo,hi],'n_pix':int(len(r)),'residual_mean_flux':float(np.mean(r)),'residual_rms_flux':float(np.sqrt(np.mean(r*r))),'residual_std_flux':float(np.std(r,ddof=1)),'lag1_cov_flux':float(np.cov(r[:-1],r[1:],ddof=1)[0,1]) if len(r)>2 else None,'lag1_corr':float(np.corrcoef(r[:-1],r[1:])[0,1]) if len(r)>2 else None,'chi2':float(np.sum(z*z))})
    out={'schema':'rya1220.profile_residual_decomposition.v1','route':a.region,'fit_A_N':a.a_n,'C':a.c,'O':a.o,'indicator_id':'CN_AX_IR:12C14N:A-X(0-0):windows_sha256:'+hashlib.sha256(json.dumps([list(x) for x in diag.windows_A]).encode()).hexdigest(),'windows':rows,'aggregate':{'n_pix':sum(x['n_pix'] for x in rows),'rms_flux_weighted':float(np.sqrt(sum(x['n_pix']*x['residual_rms_flux']**2 for x in rows)/sum(x['n_pix'] for x in rows))),'mean_lag1_corr':float(np.mean([x['lag1_corr'] for x in rows]))}}
    (a.out/'profile_residuals.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__': main()
