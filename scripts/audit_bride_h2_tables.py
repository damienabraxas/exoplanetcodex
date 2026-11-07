"""Bounded source-table and fixed-ne H chemistry checks; no RT/EOS admission."""
import argparse,json,sys,platform
from pathlib import Path
import numpy as np
from scipy.constants import k,h,e,atomic_mass
from scipy.interpolate import CubicSpline,PchipInterpolator
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only science control')
    from config.constants import require_codex_path
    from pipeline.bride_h2_tables import H2Tables,read_wide_table
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_opacity_binding import binding_digest
    import lightweaver as lw
    from lightweaver import rh_atoms
    from lightweaver.atomic_table import PeriodicTable
    source=require_codex_path('experiments.bride_molecular_tables');table=H2Tables(source)
    path=ROOT/'data/results/rya1281/solar_bound_v1/runs/rya1281-solar-bound-v1/solar5772_bound_Fe_sampled/experimental_input.json'
    b=json.loads(path.read_text());binding=b['local_binding']
    if binding_digest(binding)!=b['binding_sha256']:raise ValueError('Changed binding')
    T=np.array(binding['fields']['temperature_K']);N=np.array(binding['nH_m3'])
    a=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=np.array(binding['fields']['height_m']),temperature=T.copy(),ne=np.array(binding['ne_m3']),nHTot=N.copy(),vturb=np.full(len(T),binding['scalars']['vturb_m_s']),vlos=np.zeros(len(T)))
    abundance=lw.AtomicAbundance({PeriodicTable[i+1]:v for i,v in enumerate(binding['native_abundances_dex'])})
    eq=lw.RadiativeSet([rh_atoms.H_6_atom(),rh_atoms.He_9_atom()],abundance=abundance).compute_eq_pops(a)
    state=eq.atomicPops[PeriodicTable['H']];neutral=state.n[[l.stage==0 for l in state.model.levels]].sum(axis=0)
    f0=neutral/state.n.sum(axis=0);minus=eq.HminPops/neutral
    fields=table.hydrogen_network(T,N,f0,minus)
    # Independent bounded root in fractional HI, rather than the quadratic formula.
    roots=np.array([brentq(lambda x:(1/f+m)*x+2*ph*n*x*x-1,0.,1.,xtol=1e-14) for f,m,ph,n in zip(f0,minus,table.phi(T),N)])
    root_error=float(np.max(abs(roots-fields['HI_m3']/N)))
    if root_error>1e-12 or np.max(abs(fields['nuclei_residual']))>1e-12:raise ValueError('Hydrogen accounting failed')
    oldT,old=read_wide_table(source/'table7.dat','equilibrium')
    _,corrected=read_wide_table(source/'table7_vNov2022.dat','equilibrium')
    np.testing.assert_array_equal(old['H2'],corrected['H2'])
    x=np.log10(table.temperature);z=table.trend;loo={}
    for method in ('not-a-knot','natural','pchip'):
        errors=[]
        for i in range(1,len(x)-1):
            mask=np.arange(len(x))!=i
            spline=PchipInterpolator(x[mask],z[mask],extrapolate=False) if method=='pchip' else CubicSpline(x[mask],z[mask],bc_type=method,extrapolate=False)
            errors.append(float(np.expm1(-np.log(10)*(spline(x[i])-z[i]))))
        loo[method]=dict(temperature_K=table.temperature[1:-1].tolist(),relative_phi_error=errors,max_absolute_relative_error=max(map(abs,errors)))
    natural=CubicSpline(x,z,bc_type='natural',extrapolate=False)
    boundary_delta=np.expm1(-np.log(10)*(natural(np.log10(T))-table.interpolator(np.log10(T))))
    # Eq 7 uses modern SI constants and LW's atomic H mass; published Kp stays authoritative.
    mass=PeriodicTable['H'].mass*atomic_mass
    recon=1.5*np.log(2*np.pi*(mass/2)*k*table.temperature/h**2)+np.log(k*table.temperature)+2*np.log(table.partition_HI)-np.log(table.partition_H2)-table.dissociation_eV*e/(k*table.temperature)
    equation_delta=recon/np.log(10)-table.log10_Kp_Pa
    args.output.mkdir(parents=True,exist_ok=False)
    csv=args.output/'hydrogen_column.csv'
    np.savetxt(csv,np.column_stack([np.arange(len(T)),T,N,table.phi(T),f0,minus,*fields.values()]),delimiter=',',header='depth_index,T_K,total_H_m3,phi_m3,atomic_neutral_fraction,Hminus_per_HI,'+','.join(fields),comments='')
    report=dict(schema='codex.bride_h2_tables_audit/1',scientific_admission='HOLD',generator=fingerprint(__file__),adapter=fingerprint(ROOT/'pipeline/bride_h2_tables.py'),source_receipt=table.receipt,input=fingerprint(path),source_domain_qualified_K=table.domain,solar_temperature_range_K=[float(min(T)),float(max(T))],covered_depth_count=len(T),interpolation='Paper Sec 3: cubic spline of log10(Kp)+D0*5040/T versus log10(T); not-a-knot boundary; no extrapolation',leave_one_knot_out=loo,max_solar_natural_boundary_relative_phi_delta=float(max(abs(boundary_delta))),H2_corrected_original_exactly_equal=True,max_nuclei_residual=float(max(abs(fields['nuclei_residual']))),max_independent_root_fraction_error=root_error,max_H2_over_total_H=float(max(fields['H2_m3']/N)),max_HI_depletion_relative_atomic_reservoir=float(max(1-fields['HI_m3']/(N*f0))),equation7_reconstruction_delta_dex_range=[float(min(equation_delta)),float(max(equation_delta))],equation7_residual_disposition='Unassigned finite source reconstruction difference; published Kp unchanged, no fitted constants or correction',column=fingerprint(csv),scope='Fixed-ne H I/H II/H minus/H2 nuclei bookkeeping only. No metal/CNO network, charge feedback, continuum, collisions, RT or O SE update. Knot holdout errors are interpolation diagnostics, not physical accuracy bounds.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:report[key] for key in ('covered_depth_count','max_nuclei_residual','max_independent_root_fraction_error','max_H2_over_total_H','max_solar_natural_boundary_relative_phi_delta','equation7_reconstruction_delta_dex_range')}))

if __name__=='__main__':main()
