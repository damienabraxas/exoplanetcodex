"""Explicit paired-H O contexts and common-fraction initialization, no defaults."""
import numpy as np


def scale_initial_fractions(populations,target_total):
    pop=np.asarray(populations,float);total=np.asarray(target_total,float)
    if pop.ndim!=2 or total.shape!=(pop.shape[1],)or not np.isfinite(pop).all()or np.any(pop<=0)or not np.isfinite(total).all()or np.any(total<=0):raise ValueError('Invalid common fraction initialization')
    return pop/pop.sum(axis=0)[None,:]*total[None,:]


def make_context(state,structure,atom,*,rays,qscale,extra_vac_nm,cap_length_m,ng,fixed=None,continuum_only=False):
    import lightweaver as lw
    from lightweaver import rh_atoms
    from lightweaver.atomic_table import PeriodicTable
    from pipeline.bride_binding_controls import BindingControl
    from pipeline.bride_reference import PASSIVE_ATOM_FACTORIES
    from pipeline.multi_atom import to_lightweaver
    from pipeline.bride_top_cap import FiniteColumnCap
    from pipeline.bride_cn_provider import EqLookupAudit
    if not isinstance(state,BindingControl)or state.case not in ('model_mass_H','native_network_H'):raise ValueError('Explicit paired-H state required')
    extra=np.asarray(extra_vac_nm,float)
    if extra.ndim!=1 or len(extra)<2 or not np.isfinite(extra).all()or np.any(extra<=0)or np.any(np.diff(extra)<=0):raise ValueError('Invalid explicit diagnostic endpoints')
    atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=np.array(state.height_m),temperature=np.array(state.temperature_K),ne=np.array(state.ne_m3),nHTot=np.array(state.nH_m3),vturb=np.array(state.vturb_m_s),vlos=np.zeros(len(state.nH_m3)),logG=structure['logg_cgs']-np.log10(100.),lowerBc=lw.ThermalisedRadiation(),upperBc=FiniteColumnCap(cap_length_m));atmos.quadrature(rays)
    model=to_lightweaver(atom)
    if continuum_only:model.lines=[]
    if qscale!=1:
        for line in model.lines:line.quadrature.Nlambda=qscale*(line.quadrature.Nlambda-1)+1
    abundance={PeriodicTable[z+1]:value for z,value in enumerate(structure['abundances_dex'])};abundance[PeriodicTable['O']]=atom.abundance
    rad=lw.RadiativeSet([getattr(rh_atoms,x)()for x in PASSIVE_ATOM_FACTORIES]+[model],abundance=lw.AtomicAbundance(abundance))
    if fixed is None:rad.set_active('O')
    else:rad.set_detailed_static('O')
    eq=rad.compute_eq_pops(atmos,mols=lw.MolecularTable([]))
    if fixed is not None:
        pop=eq.atomicPops[PeriodicTable['O']].n
        if fixed.shape!=pop.shape or not np.isfinite(fixed).all()or np.any(fixed<=0)or not np.allclose(fixed.sum(axis=0),pop.sum(axis=0),rtol=1e-10,atol=0):raise ValueError('Final-state continuum binding failed')
        pop[:]=fixed
    audits=[]
    def factory(e,r,w):
        proxy=EqLookupAudit(e);audits.append(proxy);return lw.LwCompiled.FastBackground(proxy,r,w,Nthreads=1)
    spect=rad.compute_wavelength_grid(extraWavelengths=extra)
    options=lw.NgOptions(Norder=3,Nperiod=5,Ndelay=10,threshold=.05,lowerThreshold=0.)if ng and fixed is None else None
    ctx=lw.Context(atmos,spect,eq,Nthreads=1,formalSolver='piecewise_besser_1d',conserveCharge=False,ngOptions=options,backgroundProvider=factory);ctx.depthData.fill=True
    for field,expected in [('temperature',state.temperature_K),('ne',state.ne_m3),('nHTot',state.nH_m3),('vturb',state.vturb_m_s),('height',state.height_m)]:
        if not np.array_equal(getattr(atmos,field),expected):raise ValueError('Context changed paired state '+field)
    if len(eq.molecularTable.indices):raise ValueError('Unexpected LW molecular populations')
    return ctx,audits
