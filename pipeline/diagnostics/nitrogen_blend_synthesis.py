from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--engine-root',type=Path,required=True); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--a-n',type=float,default=8.188); ap.add_argument('--c',type=float,default=8.46); ap.add_argument('--o',type=float,default=8.73); a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=False); sys.path.insert(0,str(a.engine_root.resolve()))
    from pipeline import cno_synthesis as c
    from pipeline.uncertainty_stack import params_and_deltas
    region=c.REGIONS['nir_cn_iag']; ds=[d for d in c.REGION_DIAGNOSTICS['nir_cn_iag'] if d.element=='N' and d.role=='primary']
    params,_=params_and_deltas('solar'); broad=c.preflight(region,'solar',tuple(ds)); llpath,_,gf=c.region_atomic_linelist(region,tuple(ds)); ll,iso,chem=c._load_synth_resources(str(llpath),apply_canonical_gf=(gf or {}).get('apply_canonical_gf',False)) if llpath else c._load_synth_resources(); sab=c.ispec.read_solar_abundances(c._ISPEC_SOLAR_ABUND_FILE); codes=c._atom_codes(('C','N','O','Ni'),chem,sab); ow,of=c._load_region_spectrum('solar',region); atm=c._load_atmosphere(params['teff_K'],params['logg'],params['feh'],params['vturb_kms']); state=c._seed_abundances('solar',params,codes,c._solar_A(('C','N','O','Ni'),chem,sab),params['feh']); state.update({'C':a.c,'N':a.a_n,'O':a.o}); fa=c._fixed_ab(state,codes); rows=[]; tmp=str(a.out/'synth'); Path(tmp).mkdir()
    elements=np.asarray(ll['element']).astype(str) if 'element' in ll.dtype.names else None
    n_i_mask=np.char.upper(elements)=='N 1' if elements is not None else None
    for d in ds:
      lo,hi=d.windows_A[0]; m=(ow*10>=lo)&(ow*10<=hi); ow_i,of_i=ow[m],of[m]; sw=np.arange(lo/10,hi/10+c._WSTEP_NM*.5,c._WSTEP_NM)
      full=c._synth_window(sw,atm,params,ll,iso,sab,fa,broad,True,tmp)
      # This CN window is a molecular-N diagnostic.  Record whether an atomic N I
      # component is actually present before constructing the counterfactual.
      no_n=c._synth_window(sw,atm,params,ll[~n_i_mask],iso,sab,fa,broad,True,tmp) if elements is not None else full
      no_mol=c._synth_window(sw,atm,params,ll,iso,sab,fa,broad,False,tmp)
      def stat(x):
       y=np.interp(ow_i,sw,x,left=1,right=1); return {'mean_depth':float(np.mean(1-y)),'rms_vs_observed':float(np.sqrt(np.mean((of_i-y)**2))),'mean_abs_delta_from_full':float(np.mean(np.abs(y-np.interp(ow_i,sw,full,left=1,right=1))))}
      rows.append({'key':d.key,'window_air_A':[lo,hi],'n_pix':int(len(ow_i)),'atomic_N_I_lines_in_linelist':int(np.sum(n_i_mask)) if n_i_mask is not None else None,'full_synthesis':stat(full),'without_N_I':stat(no_n),'without_molecules':stat(no_mol),'N_I_model_effect_depth':float(stat(full)['mean_depth']-stat(no_n)['mean_depth']),'molecular_model_effect_depth':float(stat(full)['mean_depth']-stat(no_mol)['mean_depth'])})
    (a.out/'nitrogen_blend_synthesis.json').write_text(json.dumps({'schema':'rya1220.nitrogen_blend_synthesis.v2','route':'solar IAG NIR CN','A_N':a.a_n,'A_C':a.c,'A_O':a.o,'method':'full Turbospectrum component synthesis; atomic-N-I and no-molecule counterfactuals','windows':rows},indent=2)+'\n')
if __name__=='__main__': main()
