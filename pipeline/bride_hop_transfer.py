"""Native-state fixed-population transfer only; unqualified coupled SE refuses."""
import re
import numpy as np
from lightweaver.LwCompiled import BackgroundProvider
from lightweaver.broadening import StandardLineBroadener
from pipeline.hop_packet import compute_packet
from pipeline.reference_engines import fingerprint,load_engine


def native_lobatto5():
    """Read the pinned native five-ray constants at their compiled default precision."""
    _,app,_,_=load_engine('openmulti3d');source=app/'SRC/ray_mod.f90';text=source.read_text().split('SUBROUTINE set_lobatto_quadratures',1)[1].split('END SUBROUTINE',1)[0]
    block=re.search(r'case\s*\(5\)(.*?)case\s*\(6\)',text,re.S|re.I)
    if block is None:raise ValueError('Changed native quadrature source')
    arrays=[]
    for name in ['mu','wmu']:
        match=re.search(r'\b'+name+r'\s*=\s*\(/(.*?)/\)',block[1],re.S|re.I)
        if match is None:raise ValueError('Missing native angular table')
        values=re.findall(r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[deDE][+-]?\d+)?',match[1])
        # The pinned native run reports default REAL kind4. These literals have
        # no D exponent: Fortran rounds them to real4 before assigning real8.
        if len(values)!=5 or any('d'in x.lower()for x in values):raise ValueError('Changed native literal precision/count')
        arrays.append(np.asarray(values,dtype=np.float32).astype(float))
    mu,weight=arrays
    if not np.all((mu>0)&(mu<=1))or not np.isclose(weight.sum(),1):raise ValueError('Invalid native angular table')
    return mu,weight,dict(source=fingerprint(source),literal_precision='default real4 promoted to real8; no weight renormalization',mu=mu.tolist(),weight=weight.tolist())


class Perturbers:
    """Broadening-only density interface; never an H level-population packet."""
    def __init__(self,packet,*,helium=True):
        if not isinstance(helium,bool):raise ValueError('Require explicit He broadening diagnostic boolean')
        self.packet=packet;self.helium=helium;self.abundance=self;self.reads=[]

    def __getitem__(self,key):
        if not isinstance(key,str)and getattr(key,'Z',None)==2:
            if not self.helium:
                self.reads.append('explicit He broadening-off diagnostic, EOS/continuum unchanged')
                return np.zeros_like(self.packet['He_I_m3'])
            self.reads.append('effective HeI/HI perturber ratio, not chemical abundance')
            return self.packet['He_I_m3']/self.packet['H_I_m3']
        if isinstance(key,str)and key=='H':
            owner=self
            class Rows:
                def __getitem__(self,index):
                    row,tail=(index[0],index[1:])if isinstance(index,tuple)else(index,())
                    if row not in [0,-1]:raise ValueError('Unqualified hydrogen excited-state perturber lookup')
                    name='H_I_m3'if row==0 else'H_II_m3';owner.reads.append(name+' as explicit broadening perturber, not H model populations')
                    return owner.packet[name][tail]if tail else owner.packet[name]
            return Rows()
        if isinstance(key,str)and key=='H-':
            self.reads.append('native Hminus density')
            return self.packet['H_-_m3']
        raise ValueError('Unqualified perturber lookup '+str(key))


def verify_helium_reads(reads,enabled):
    """Require the declared control to reach actual compiled-context broadeners."""
    on=any('effective HeI/HI' in x for x in reads)
    off=any('explicit He broadening-off' in x for x in reads)
    if (enabled and (not on or off))or(not enabled and (not off or on)):
        raise ValueError('Declared He broadening control did not reach actual perturber reads')


class NativePerturberBroadener(StandardLineBroadener):
    """Retain the original LW formula/coefficients with explicit native densities."""
    def __init__(self,inner,perturbers):self.inner=inner;self.perturbers=perturbers
    def setup(self,line):self.inner.setup(line)
    def broaden(self,atmos,eqPops):return self.inner.broaden(atmos,self.perturbers)


class NativeHopBackground(BackgroundProvider):
    def __init__(self,eqPops,radSet,wavelength,*,state,atom,abundances,records,coupled_binding=None):
        from lightweaver.atomic_table import PeriodicTable
        super().__init__(eqPops,radSet,wavelength)
        if coupled_binding is None:
            if radSet.activeAtoms or len(radSet.detailedAtoms)!=1 or radSet.detailedAtoms[0].element!=PeriodicTable['O']:raise ValueError('HOP background qualified for fixed native-O transfer only; no coupled SE')
        else:
            from pipeline.bride_collision_binding import FixedNativeCollisionBinding
            if not isinstance(coupled_binding,FixedNativeCollisionBinding)or not coupled_binding.receipt or len(radSet.activeAtoms)!=1 or radSet.detailedAtoms or radSet.activeAtoms[0]is not coupled_binding.atom or (len(coupled_binding.atom.collisions)!=1 or coupled_binding.atom.collisions[0]is not coupled_binding):raise ValueError('Unqualified explicit active-O collision contract')
        self.state={k:np.array(v,copy=True)for k,v in state.items()};self.wave=np.asarray(wavelength).copy();self.original_wave=wavelength
        for value in self.state.values():value.flags.writeable=False
        self.packet,self.meta=compute_packet(self.state['temperature_K'],self.state['ne_m3'],self.state['rho_kg_m3'],self.wave,abundances,atom,vmic_m_s=float(self.state['vturb_m_s'][0]))
        np.testing.assert_array_equal(self.packet['nHTot_m3'],self.state['nHTot_m3'])
        self.check(eqPops.atmosphere);records.append(dict(wavelengths=len(self.wave),depths=len(state['temperature_K']),packet=self.meta))

    def check(self,atmos):
        for name,key in [('temperature','temperature_K'),('ne','ne_m3'),('nHTot','nHTot_m3'),('vturb','vturb_m_s'),('height','height_m')]:
            if not np.array_equal(getattr(atmos,name),self.state[key]):raise ValueError('Changed bound HOP state '+name)
        if not np.array_equal(self.wave,self.original_wave):raise ValueError('Changed bound HOP wavelengths')

    def compute_background(self,atmos,chi,eta,sca):
        self.check(atmos.pyAtmos)
        for output,key in [(chi,'chi_m_inv'),(eta,'eta_SI'),(sca,'scat_m_inv')]:
            if output.shape!=self.packet[key].shape:raise ValueError('Background wavelength/depth axis mismatch')
            output[:]=self.packet[key]


def make_context(state,structure,atom,populations,nstar,*,extra_vac_nm,continuum_only=False,formal_solver='piecewise_besser_1d',helium_broadening=True,table_control=None,vdw_broadening=True,original_width_law=None,collision_binding=None,coupled_diagnostic=False,top_boundary=None):
    """Bind complete native O state and HOP background before Context compilation."""
    if formal_solver not in ['piecewise_besser_1d','piecewise_linear_1d']:raise ValueError('Unsupported fixed-transfer formal solver')
    if not isinstance(vdw_broadening,bool):raise ValueError('Require explicit VDW diagnostic boolean')
    if original_width_law is not None and (helium_broadening or not vdw_broadening):raise ValueError('Original source law owns full VDW/He prescription; no added He factor or VDW omission')
    if not isinstance(coupled_diagnostic,bool):raise ValueError('Require explicit coupled diagnostic boolean')
    if coupled_diagnostic and (collision_binding is None or original_width_law is None or table_control is not None or continuum_only):raise ValueError('Coupled diagnostic requires full source atom/width/collision contract; no transfer overrides')
    import lightweaver as lw
    from lightweaver import rh_atoms
    from lightweaver.atomic_table import PeriodicTable
    from pipeline.multi_atom import to_lightweaver
    from pipeline.bride_reference import PASSIVE_ATOM_FACTORIES
    from pipeline.bride_top_cap import FiniteColumnCap
    if not np.all(np.asarray(state['vturb_m_s'])==state['vturb_m_s'][0]):raise ValueError('Only depth-constant vmic is qualified')
    expected=(len(atom.levels),len(state['temperature_K']))
    if any(x.shape!=expected or not np.isfinite(x).all()or np.any(x<=0)for x in [populations,nstar]):raise ValueError('Invalid complete native atomic state')
    np.testing.assert_allclose(populations.sum(axis=0),nstar.sum(axis=0),rtol=1e-12,atol=0)
    seed,seed_meta=compute_packet(state['temperature_K'],state['ne_m3'],state['rho_kg_m3'],np.array([500.,845.]),structure['abundances_dex'],atom,vmic_m_s=float(state['vturb_m_s'][0]))
    state={key:np.array(value,copy=True)for key,value in state.items()};state['nHTot_m3']=seed['nHTot_m3'].copy()
    for value in state.values():value.flags.writeable=False
    atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=state['height_m'].copy(),temperature=state['temperature_K'].copy(),ne=state['ne_m3'].copy(),nHTot=state['nHTot_m3'].copy(),vturb=state['vturb_m_s'].copy(),vlos=np.zeros(expected[1]),logG=structure['logg_cgs']-2,lowerBc=lw.ThermalisedRadiation(),upperBc=FiniteColumnCap(0.)if top_boundary is None else top_boundary)
    mu,weight,angles=native_lobatto5();atmos.quadrature(Nrays=5,mu=mu,wmu=weight)
    model=to_lightweaver(atom);perturbers=Perturbers(seed,helium=helium_broadening)
    table_audit=None
    if table_control is not None:
        from pipeline.bride_native_table_control import apply_table_control
        table_audit=apply_table_control(model,*table_control)
    source_width_reads=[]
    for index,line in enumerate(model.lines):
        line.quadrature.Nlambda=2*(line.quadrature.Nlambda-1)+1
        if original_width_law is None:
            if not vdw_broadening:line.broadening.elastic=[]
            line.broadening.elastic=[NativePerturberBroadener(x,perturbers)for x in line.broadening.elastic]
        else:
            from pipeline.multi_width_law import OriginalMultiVdwLaw
            width=OriginalMultiVdwLaw(original_width_law.records[index],seed['H_I_m3'],original_width_law.temperature_interval);line.broadening.elastic=[width];line.setup(model);source_width_reads.append(width.audit)
    if continuum_only:model.lines=[]
    abundance={PeriodicTable[z+1]:value for z,value in enumerate(structure['abundances_dex'])};abundance[PeriodicTable['O']]=atom.abundance
    rad=lw.RadiativeSet([getattr(rh_atoms,name)()for name in PASSIVE_ATOM_FACTORIES]+[model],abundance=lw.AtomicAbundance(abundance))
    if coupled_diagnostic:rad.set_active('O')
    else:rad.set_detailed_static('O')
    eq=rad.compute_eq_pops(atmos,mols=lw.MolecularTable([]));oxygen=eq.atomicPops[PeriodicTable['O']]
    oxygen.nStar[:]=nstar;oxygen.nTotal[:]=nstar.sum(axis=0);oxygen.n[:]=populations
    if collision_binding is not None:
        collision_binding.setup(model);model.collisions=[collision_binding]
    spect=rad.compute_wavelength_grid(extraWavelengths=np.asarray(extra_vac_nm));records=[]
    def factory(e,r,w):return NativeHopBackground(e,r,w,state=state,atom=atom,abundances=structure['abundances_dex'],records=records,coupled_binding=collision_binding if coupled_diagnostic else None)
    ng=lw.NgOptions(Norder=3,Nperiod=5,Ndelay=10,threshold=.05,lowerThreshold=0.)if coupled_diagnostic else None
    ctx=lw.Context(atmos,spect,eq,Nthreads=1,formalSolver=formal_solver,conserveCharge=False,backgroundProvider=factory,ngOptions=ng);ctx.depthData.fill=True
    if not continuum_only and vdw_broadening and original_width_law is None:verify_helium_reads(perturbers.reads,helium_broadening)
    if not vdw_broadening and (perturbers.reads or any(line.broadening.elastic for line in model.lines)):raise ValueError('VDW-off diagnostic retained elastic broadening reads/terms')
    if original_width_law is not None and not continuum_only and (len(source_width_reads)!=256 or any(x['calls']==0 for x in source_width_reads)or perturbers.reads):raise ValueError('Original full-width law did not reach actual model without extra perturbation')
    return ctx,dict(coupled_native_contract_diagnostic=coupled_diagnostic,original_width_law=None if original_width_law is None else dict(source_mass_amu=original_width_law.source_mass_amu,temperature_numeric_interval=original_width_law.temperature_interval,inputs=original_width_law.receipt,actual_calls=source_width_reads),vdw_broadening=vdw_broadening,elastic_terms_remaining=sum(len(line.broadening.elastic)for line in model.lines),table_control=table_audit,helium_broadening=helium_broadening,formal_solver=formal_solver,seed=seed_meta,angles=angles,background_records=records,perturbers=perturbers.reads,scope='Fixed complete native O n/nstar/totals; native EOS/continuum and perturber densities; LW atomic formulas with explicit LW or original-source broadening prescription. Native source constants retained; lower ghost/formal-solver conventions remain distinct. No SE or production admission.')
