"""Fresh fixed-population radiation/rate replay of one exact-bound solar control."""
import argparse
import json
from pathlib import Path
import platform
import sys
import time
import subprocess
from importlib.metadata import version
from dataclasses import replace

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 bounded replay only')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.bride_parameter_checks import interpolate_marcs_structure,marcs_filename
    from pipeline.bride_solar_control import density_column_structure,bound_fe_factory
    from pipeline.bride_opacity_binding import load_bound_coefficients
    from pipeline.bride_top_cap import FiniteColumnCap,converge_cap_radiation
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
    from pipeline.bride_charge_contract import install_author_mn_pair
    from pipeline.bride_rate_audit import process_flow,collision_classes
    import lightweaver.constants as const
    origin=args.case_directory;summary=json.loads((origin.parent/'summary.json').read_text())
    entry=next(x for x in summary['cases']if x['name']==origin.name)
    if fingerprint(origin/'manifest.json')!=entry['manifest']:raise ValueError('Changed run manifest')
    m=json.loads((origin/'manifest.json').read_text());case=m['case']
    if m['status']!='CONVERGED_EXPERIMENT' or case.get('experimental_column_contract')!='density_integral_log_native_top' or m['experimental_collision_convention']!='barklem2018_MN_LTE_pair':
        raise ValueError('Unsupported or incomplete bound experiment')
    if fingerprint(origin/'experimental_input.json')!=m['experimental_bound_input']['receipt']:raise ValueError('Changed experiment binding')
    if fingerprint(ROOT/'pipeline/bride_reference.py')!=summary['context_builder'] or fingerprint(ROOT/'pipeline/bride_solar_control.py')!=m['experimental_bound_input']['adapter']:
        raise ValueError('Changed experimental context/provider code')
    for name in ('pipeline/bride_fe_background.py','pipeline/bride_top_cap.py'):
        recorded=subprocess.check_output(['git','show',summary['code_commit']+':'+name],cwd=ROOT)
        if recorded!=(ROOT/name).read_bytes():raise ValueError('Changed inherited provider/cap equation')
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path)!=summary['atom']:raise ValueError('Changed atom')
    source=replace(read_multi_atom(atom_path),abundance=case['abundance'])
    store=require_codex_path('grids.marcs_standard');nodes=[]
    for (t,g),receipt in zip([(t,g)for t in case['bracket_teff']for g in case['bracket_logg']],m['atmosphere']['sources']):
        path=store/marcs_filename(t,g)
        if fingerprint(path)!=receipt:raise ValueError('Changed atmosphere source')
        nodes.append(((t,g),read_marcs_structure(path)))
    base,interpolation=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),thermodynamic_scale=case['thermodynamic_scale'])
    s,contract=density_column_structure(base,nodes,interpolation)
    certificate=require_codex_path('experiments.bride_bound_point_coefficients')/case['bound_point_coefficients']
    provider,bound=bound_fe_factory(certificate,base,s,source.abundance,column_contract=contract)
    if bound!=json.loads((origin/'experimental_input.json').read_text()):raise ValueError('Reconstructed binding differs')
    data,_=load_bound_coefficients(certificate,base,source.abundance)
    with np.load(origin/'populations.npz')as saved:pop=saved['populations'].copy();nstar=saved['nstar'].copy()
    args.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
    ctx=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],fe_background=provider,
        extra_wavelengths_vac_nm=data['wavelength_vac_nm'],experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))
    active=ctx.activeAtoms[0]
    np.testing.assert_allclose(active.nStar,nstar,rtol=64*np.finfo(float).eps,atol=0)
    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    spline=install_native_collision_control(ctx,source,oracle,source_directions=True,reversible_ch0=False)
    author=install_author_mn_pair(ctx,source,require_codex_path('experiments.bride_barklem_oxygen'),oracle)
    if spline!=m['native_spline_control'] or author!=m['author_charge_control']:raise ValueError('Changed collision contract')
    active.n[:]=pop
    radiation=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=start+180)
    np.testing.assert_array_equal(active.n,pop)
    groups=collision_classes(active.atomicModel,ctx.atmos.pyAtmos,ctx.kwargs['eqPops'])
    active.compute_collisions(fillDiagonal=False);collision=sum(groups.values())
    np.testing.assert_allclose(collision,active.C,rtol=256*np.finfo(float).eps,atol=0)
    radiative=np.zeros_like(collision)
    for transition in active.trans:
        radiative[transition.j,transition.i]+=transition.Rij
        radiative[transition.i,transition.j]+=transition.Rji
    net,gross=process_flow(collision+radiative,pop)
    residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
    if residual>=m['controls']['rtol']:raise ValueError('Fresh physical rate residual exceeds original solve tolerance')
    wave=np.array(ctx.spect.wavelength);nu=const.CLight/(wave[:,None]*1e-9)
    B=2*const.HPlanck*nu**3/const.CLight**2/np.expm1(const.HPlanck*nu/(const.KBoltzmann*s['temperature_K']))
    path=args.output/'local_nlte_rates.npz'
    np.savez_compressed(path,populations=pop,nstar=nstar,logtau500=s['logtau_500'],temperature_K=s['temperature_K'],
        wavelength_vac_nm=wave,Jnu_W_m2_sr_Hz=np.array(ctx.spect.J),Bnu_W_m2_sr_Hz=B,
        radiative_rates_s_inv=radiative,collision_rates_s_inv=collision,net_flow_m3_s=net,gross_flow_m3_s=gross)
    result=dict(schema='codex.bride_bound_rate_replay/1',status='DIAGNOSTIC_COMPLETE',scientific_admission='HOLD',
        generator=fingerprint(__file__),context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        source_manifest=fingerprint(origin/'manifest.json'),source_populations=fingerprint(origin/'populations.npz'),
        atom=summary['atom'],experimental_collision_convention=m['experimental_collision_convention'],
        experimental_bound_input=m['experimental_bound_input'],experimental_column_contract=contract,
        finite_column_top_cap=m['finite_column_top_cap'],modes=[dict(name='local_nlte',array=fingerprint(path),
            max_physical_rate_residual=residual,radiation=radiation)],elapsed_seconds=time.perf_counter()-start,
        interpretation='Fresh matched radiation with immutable saved solar populations, no SE update or external-reference comparison; original outer1e-6 tolerance retained.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(residual=residual,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
