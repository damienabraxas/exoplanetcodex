"""Actual LW synthetic background-array diagnostic; no RT or physical CN state."""
import argparse,copy,json,os,platform,resource,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def synthetic_problem(*,wrapper=False,enabled=False,zero=False,packet_mutator=None):
    import lightweaver as lw
    import lightweaver.constants as c
    from lightweaver import rh_atoms
    from pipeline.bride_cn_provider import SyntheticCnPacket,binding,synthetic_cn_factory
    from pipeline.wavelength_util import air_to_vac
    plan=json.loads((ROOT/'config/experiments/rya1281_cn_provider_arrays.json').read_text());n=plan['depths']
    atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=np.array(plan['height_m']),temperature=np.array(plan['temperature_K'],dtype=float),ne=np.array(plan['ne_m3']),nHTot=np.array(plan['nH_m3']),vturb=np.full(n,plan['xi_m_s']),vlos=np.zeros(n),logG=2.438);atmos.quadrature(2)
    rad=lw.RadiativeSet([rh_atoms.H_6_atom()],abundance=lw.AtomicAbundance());rad.set_active('H');eq=rad.compute_eq_pops(atmos)
    centre=float(air_to_vac([plan['line_air_A']])[0])*1e-10
    width=np.sqrt(2*c.KBoltzmann*atmos.temperature/(plan['mass_amu']*c.Amu)+atmos.vturb**2)/centre
    state,abundance=binding(atmos,rad.abundance)
    packet=SyntheticCnPacket(purpose=plan['purpose'],state_binding=state,abundance_binding=abundance,line_id='synthetic-CN-array-01',lambda0_air_A=plan['line_air_A'],lower_population_m3=tuple(np.zeros(n)if zero else plan['lower_population_m3']),width_Hz=tuple(width),damping_a=tuple(plan['gamma_s_inv']/(4*np.pi*width)),stimulated_emission=tuple(-np.expm1(-c.HPlanck*c.CLight/(c.KBoltzmann*centre*atmos.temperature))),integrated_strength_m2_Hz=tuple(np.full(n,c.QElectron**2/(4*c.Epsilon0*c.MElectron*c.CLight))),gf=plan['gf'])
    if packet_mutator is not None:packet=packet_mutator(packet)
    instances=[]
    factory=None
    if wrapper:
        base=synthetic_cn_factory(packet,enabled=enabled,cn_already_present=False,support_vac_nm=plan['support_vac_nm'],mapping=plan['mapping'])
        def factory(e,r,w):
            obj=base(e,r,w);instances.append(obj);return obj
    spect=rad.compute_wavelength_grid(extraWavelengths=np.linspace(*plan['support_vac_nm'],41))
    ctx=lw.Context(atmos,spect,eq,Nthreads=1,backgroundProvider=factory)
    return ctx,packet,instances


def checks():
    import lightweaver as lw
    from pipeline.bride_cn_provider import BIND_FIELDS
    bare,_,_=synthetic_problem();disabled,_,di=synthetic_problem(wrapper=True);enabled,packet,ei=synthetic_problem(wrapper=True,enabled=True);zero,_,zi=synthetic_problem(wrapper=True,enabled=True,zero=True)
    fields=('chi','eta','sca')
    for other in [disabled,zero]:
        for key in fields:np.testing.assert_array_equal(getattr(bare.background,key),getattr(other.background,key))
    provider=ei[0];mask=provider.mask;bound_multiple=16;eps=np.finfo(float).eps;errors={}
    for key,added in [('chi',provider.added_chi),('eta',provider.added_eta)]:
        total=np.asarray(getattr(enabled.background,key));base=np.asarray(getattr(bare.background,key));residual=np.abs((total-base)[mask]-added);bound=bound_multiple*eps*(np.abs(total[mask])+np.abs(base[mask])+np.abs(added))
        if np.any(residual>bound):raise ValueError('Synthetic provider additive array check failed')
        np.testing.assert_array_equal(total[~mask],base[~mask]);errors[key]=float(np.max(residual/np.maximum(bound,np.finfo(float).tiny)))
    np.testing.assert_array_equal(enabled.background.sca,bare.background.sca)
    expected=lw.planck(enabled.atmos.pyAtmos.temperature[None,:],provider.wave[mask,None])
    if not np.all(provider.added_chi>0)or not np.isfinite(expected).all():raise ValueError('Nonvacuous positive CN/Planck check required')
    kirchhoff=float(np.max(np.abs(provider.added_eta/provider.added_chi/expected-1)))
    if not np.isfinite(kirchhoff)or kirchhoff>32*eps:raise ValueError('Synthetic provider Kirchhoff check failed')
    snapshot=[np.array(getattr(enabled.background,k))for k in fields];enabled.background.update_background(enabled.atmos)
    for key,value in zip(fields,snapshot):np.testing.assert_array_equal(getattr(enabled.background,key),value)
    clone=copy.deepcopy(enabled)
    for key in fields:
        np.testing.assert_array_equal(getattr(clone.background,key),getattr(enabled.background,key))
        if np.shares_memory(np.asarray(getattr(clone.background,key)),np.asarray(getattr(enabled.background,key))):raise ValueError('Synthetic clone shares background storage')
    refusals=[]
    for key in BIND_FIELDS:
        ctx,_,_=synthetic_problem(wrapper=True,enabled=True);before=[np.array(getattr(ctx.background,k))for k in fields];array=getattr(ctx.atmos.pyAtmos,key);array[0]=np.nextafter(array[0],np.inf)
        try:ctx.background.update_background(ctx.atmos)
        except ValueError:refusals.append(key)
        else:raise ValueError('Mutated synthetic state accepted '+key)
        for name,value in zip(fields,before):np.testing.assert_array_equal(getattr(ctx.background,name),value)
    for kind in ('wavelength','runtime_abundance'):
        ctx,_,objects=synthetic_problem(wrapper=True,enabled=True);before=[np.array(getattr(ctx.background,k))for k in fields]
        if kind=='wavelength':ctx.spect.wavelength[0]=np.nextafter(ctx.spect.wavelength[0],np.inf)
        else:
            from lightweaver.atomic_table import PeriodicTable
            table=objects[0].radset.abundance.abundance;key=PeriodicTable['C'];table[key]=np.nextafter(table[key],np.inf)
        try:ctx.background.update_background(ctx.atmos)
        except ValueError:refusals.append(kind)
        else:raise ValueError('Mutated CN binding accepted '+kind)
        for key,value in zip(fields,before):np.testing.assert_array_equal(getattr(ctx.background,key),value)
    from dataclasses import replace
    try:synthetic_problem(wrapper=True,enabled=True,packet_mutator=lambda p:replace(p,abundance_binding=tuple([p.abundance_binding[0]*2]+list(p.abundance_binding[1:]))))
    except ValueError:refusals.append('abundance_binding')
    else:raise ValueError('Wrong abundance binding accepted')
    requests=di[0].audit.requests+ei[0].audit.requests+zi[0].audit.requests
    if any(key=='CN'for _,key in requests):raise ValueError('Baseline requested CN')
    return dict(disabled_identity=True,zero_population_identity=True,scattering_identity=True,outside_support_identity=True,repeat_identity=True,deepcopy_identity=True,additive_error_over_rounding_bound=errors,kirchhoff_max_relative_difference=kirchhoff,mutation_refusals=refusals,base_lookup_inventory=sorted(set(key for _,key in requests)),shape=list(np.asarray(enabled.background.chi).shape),support_rows=int(mask.sum()),packet=packet,contexts=[bare,disabled,enabled,zero])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();before=time.perf_counter()
    import lightweaver,lightweaver.LwCompiled as compiled
    from importlib.metadata import version
    from pipeline.bride_reference import fingerprint
    plan_path=ROOT/'config/experiments/rya1281_cn_provider_arrays.json';plan=json.loads(plan_path.read_text())
    if version('lightweaver')!=plan['backend_version']or platform.system()!='Darwin':raise ValueError('Wrong synthetic diagnostic backend/host')
    a.output.mkdir(parents=True,exist_ok=False);binary=fingerprint(compiled.__file__);result=checks();contexts=result.pop('contexts');packet=result.pop('packet')
    from dataclasses import asdict
    (a.output/'packet.json').write_text(json.dumps(asdict(packet),indent=2,allow_nan=False)+'\n')
    arrays={label+'_'+key:np.array(getattr(ctx.background,key))for label,ctx in zip(['bare','disabled','enabled','zero'],contexts)for key in ['chi','eta','sca']};arrays['wavelength_vac_nm']=np.array(contexts[0].spect.wavelength)
    np.savez(a.output/'arrays.npz',**arrays)
    if fingerprint(compiled.__file__)!=binary:raise ValueError('Installed backend changed')
    result.update(schema='codex.bride_cn_provider_arrays/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),provider=fingerprint(ROOT/'pipeline/bride_cn_provider.py'),kernel=fingerprint(ROOT/'pipeline/bride_cn_diagnostic.py'),plan=fingerprint(plan_path),backend_version=version('lightweaver'),backend_binary=binary,python=sys.version,numpy_version=np.__version__,seconds_including_array_checks=time.perf_counter()-before,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,limitations=['H-only disposable synthetic state and synthetic CN populations, not native56depth/solar chemistry or a new oxygen result.','Runtime lookup inventory plus separately inspected LW source supports this baseline path; no independent source-to-binary rebuild attestation.','No formal solution, RT/SE, abundance, production/default model or primary provenance admission.'])
    (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['shape','support_rows','kirchhoff_max_relative_difference','mutation_refusals','base_lookup_inventory','seconds_including_array_checks','process_peak_rss_bytes']}));print(f'Read disposable synthetic plan and installed backend; wrote packet, arrays and report to {a.output}')

if __name__=='__main__':main()
