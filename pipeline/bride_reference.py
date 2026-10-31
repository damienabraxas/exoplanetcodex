"""Read exact held comparison inputs; no RT/SE or atmosphere interpolation."""
from pathlib import Path
import hashlib
import numpy as np

PASSIVE_ATOM_FACTORIES=('H_6_atom','He_9_atom','C_atom','N_atom','MgI_66_atom',
                      'Si_atom','Al_atom','Fe_atom','CaII_atom','Na_atom','S_atom')


def marcs_number_densities(structure):
    """Pinned-engine unit/constants boundary, shared by solve and replay."""
    import lightweaver.constants as constants
    from lightweaver.atomic_table import PeriodicTable
    ne=structure['electron_pressure_Pa']/(constants.KBoltzmann*structure['temperature_K'])
    nH=structure['rho_kg_m3']*structure['hydrogen_mass_fraction']/(constants.Amu*PeriodicTable['H'].mass)
    return ne,nH


def make_reference_context(structure,source_atom,*,rays,qscale,lte=False,ng=False,
                           fixed_populations=None,continuum_only=False,lte_h_background_lines=False,
                           lte_uv_background_lines=False,finite_column_top_cap=False,
                           extra_wavelengths_vac_nm=None,fe_background=None):
    """Exact benchmark physics for solve/replay; no RT, SE or abundance fitting.

    A fixed population replay is detailed-static, including its original NLTE
    bound-free opacity. continuum_only removes bound-bound transitions only;
    it never substitutes an LTE line profile for the continuum.
    """
    import lightweaver as lw
    from lightweaver import rh_atoms
    from lightweaver.atomic_table import PeriodicTable
    from pipeline.multi_atom import to_lightweaver
    ne,nH=marcs_number_densities(structure)
    n=structure['ndep']
    cap=None
    if finite_column_top_cap:
        from pipeline.bride_top_cap import cap_from_structure
        cap=cap_from_structure(structure)
    atmosphere=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,
        depthScale=structure['height_m'].copy(),temperature=structure['temperature_K'].copy(),
        ne=ne.copy(),nHTot=nH.copy(),vturb=np.full(n,structure['vturb_m_s']),
        vlos=np.zeros(n),logG=structure['logg_cgs']-np.log10(100.),
        lowerBc=lw.ThermalisedRadiation(),upperBc=cap if cap is not None else lw.ZeroRadiation())
    atmosphere.quadrature(rays)
    model=to_lightweaver(source_atom)
    if continuum_only:model.lines=[]
    if qscale!=1:
        for line in model.lines:line.quadrature.Nlambda=qscale*(line.quadrature.Nlambda-1)+1
    abundance_data={PeriodicTable[i+1]:a for i,a in enumerate(structure['abundances_dex'])}
    abundance_data[PeriodicTable['O']]=source_atom.abundance
    radset=lw.RadiativeSet([getattr(rh_atoms,name)() for name in PASSIVE_ATOM_FACTORIES]+[model],
                           abundance=lw.AtomicAbundance(abundance_data))
    if lte or fixed_populations is not None:radset.set_detailed_static('O')
    else:radset.set_active('O')
    if lte_h_background_lines:
        # Explicit LTE/CRD opacity control, not an H NLTE or PRD solution.
        from lightweaver.atomic_model import LineType
        hydrogen=radset.atoms[PeriodicTable['H']]
        for line in hydrogen.lines:line.type=LineType.CRD
        radset.set_detailed_static('H')
    if lte_uv_background_lines:
        # Bounded native-atom UV opacity diagnostic, not complete line blanketing.
        # Detailed-static species also contribute their native continua/sampling.
        from lightweaver.atomic_model import LineType
        for background in radset.atoms.values():
            if background.element == model.element:
                continue
            selected=[line for line in background.lines if 100. <= line.lambda0 <= 300.]
            if not selected:
                continue
            background.lines=selected
            for line in selected:
                line.type=LineType.CRD
            radset.set_detailed_static(background.element)
    eq=radset.compute_eq_pops(atmosphere)
    if fixed_populations is not None:
        actual=eq.atomicPops[model.element].n
        if fixed_populations.shape!=actual.shape or not np.isfinite(fixed_populations).all() or np.any(fixed_populations<=0):
            raise ValueError('Invalid fixed population replay')
        if not np.allclose(fixed_populations.sum(axis=0),actual.sum(axis=0),rtol=1e-10,atol=0):
            raise ValueError('Replay population totals differ from source atmosphere/abundance')
        actual[:]=fixed_populations
    if extra_wavelengths_vac_nm is not None:
        extra_wavelengths_vac_nm=np.asarray(extra_wavelengths_vac_nm,dtype=float)
        if (extra_wavelengths_vac_nm.ndim!=1 or not extra_wavelengths_vac_nm.size
                or not np.isfinite(extra_wavelengths_vac_nm).all()
                or np.any(extra_wavelengths_vac_nm<=0)):
            raise ValueError('Extra continuum diagnostic wavelengths must be finite positive vacuum nm')
    spect=radset.compute_wavelength_grid(extraWavelengths=extra_wavelengths_vac_nm)
    options=lw.NgOptions(Norder=3,Nperiod=5,Ndelay=10,threshold=.05,lowerThreshold=0.) if ng and not lte else None
    ctx=lw.Context(atmosphere,spect,eq,Nthreads=1,formalSolver='piecewise_besser_1d',
                   conserveCharge=False,ngOptions=options,backgroundProvider=fe_background)
    if cap is not None:ctx.depthData.fill=True
    return ctx


def fingerprint(path):
    path=Path(path)
    return dict(name=path.name,bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def read_marcs_structure(path):
    """Read the two standard MARCS structure tables in their printed precision.

    Source Depth is inward cm; output height is outward m. Pe is dyn/cm2;
    density is g/cm3. Hydrogen total is rho*X/mH, not neutral H or NLTE pops.
    Only data/source units are parsed here; no EOS or physics is substituted.
    """
    text=Path(path).read_text();lines=text.splitlines()
    model=lines[0].strip()
    count_line=next(i for i,l in enumerate(lines)if 'Number of depth points' in l)
    n=int(lines[count_line].split()[0])
    first=next(i for i,l in enumerate(lines)if 'lgTauR' in l and 'lgTau5' in l)
    second=next(i for i,l in enumerate(lines)if 'KappaRoss' in l and 'Density' in l)
    a=np.asarray([[float(x)for x in l.split()]for l in lines[first+1:first+1+n]])
    b=np.asarray([[float(x)for x in l.split()]for l in lines[second+1:second+1+n]])
    if a.shape!=(n,9)or b.shape!=(n,8)or not np.isfinite(a).all()or not np.isfinite(b).all():
        raise ValueError('Malformed/nonfinite MARCS structure tables')
    if not np.array_equal(a[:,0],np.arange(1,n+1))or not np.array_equal(a[:,0],b[:,0]):
        raise ValueError('MARCS structure row identities differ')
    if not np.array_equal(a[:,1],b[:,1]):
        raise ValueError('MARCS optical-depth identities differ across tables')
    abundance_start=next(i for i,l in enumerate(lines)if 'Logarithmic chemical number abundances' in l)
    abundances=[float(v)for line in lines[abundance_start+1:count_line]for v in line.split()]
    if len(abundances)!=92 or not np.isfinite(abundances).all()or abundances[0]!=12.:
        raise ValueError('Expected 92 finite MARCS elemental abundances with H=12')
    X=float(next(l for l in lines if 'are X, Y' in l).split()[0])
    gravity=float(next(l for l in lines if 'Surface gravity' in l).split()[0])
    vturb=float(next(l for l in lines if 'Microturbulence parameter' in l).split()[0])*1000.
    if not (0<X<1 and gravity>0 and vturb>=0 and np.all(a[:,4:7]>0) and np.all(b[:,3]>0)):
        raise ValueError('Invalid MARCS physical fields')
    if not np.all(np.diff(a[:,3])>0):
        raise ValueError('MARCS depth must increase inward')
    if np.any(b[:,7]<=0) or not np.all(np.diff(b[:,7])>0):
        raise ValueError('MARCS RHOX column mass must be positive and increase inward')
    return dict(model_id=model,ndep=n,height_m=-a[:,3]*.01,
        temperature_K=a[:,4].copy(),electron_pressure_Pa=a[:,5]*.1,
        rho_kg_m3=b[:,3]*1000.,hydrogen_mass_fraction=X,
        abundances_dex=abundances,logg_cgs=float(np.log10(gravity)),vturb_m_s=vturb,
        logtau_ross=a[:,1].copy(),logtau_500=a[:,2].copy(),
        column_mass_kg_m2=b[:,7].copy()*10.)


def read_exact_departure_record(binary,aux_text,*,model_id,abundance,nlevels,ndep):
    """Read one native MULTI grid record with exact identity, never nearest-node.

    The verified Gerber layout is 500 ASCII bytes, two little-endian i4,
    ndep linear tau doubles, then (level,depth) little-endian doubles.
    Return raw-record hash; do not copy/hash the multi-GB grid onto the Mac.
    """
    hits=[]
    for line in aux_text.splitlines():
        fields=line.split()
        if len(fields)>=9 and not line.lstrip().startswith('#'):
            if fields[0].strip("'")==model_id and float(fields[7])==abundance:
                hits.append(fields)
    if len(hits)!=1:
        raise ValueError(f'Expected one exact model/abundance record, found {len(hits)}')
    row=hits[0];offset=int(row[8])-1
    if offset<1000:raise ValueError('Invalid one-based record pointer')
    size=500+8+8*ndep*(1+nlevels)
    with Path(binary).open('rb')as f:f.seek(offset);raw=f.read(size)
    if len(raw)!=size:raise ValueError('Truncated native departure record')
    identity=raw[:500].decode('ascii').strip()
    dimensions=np.frombuffer(raw[500:508],dtype='<i4')
    if identity!=model_id or tuple(dimensions)!=(ndep,nlevels):
        raise ValueError('Native record identity or dimensions do not match requested input')
    tau=np.frombuffer(raw[508:508+8*ndep],dtype='<f8').copy()
    populations=np.frombuffer(raw[508+8*ndep:],dtype='<f8').reshape(nlevels,ndep).copy()
    if not np.isfinite(tau).all()or np.any(tau<=0)or not np.all(np.diff(tau)>0):
        raise ValueError('Invalid native linear optical-depth grid')
    if not np.isfinite(populations).all()or np.any(populations<0):
        raise ValueError('Invalid native departure coefficients')
    return dict(model_id=identity,abundance=abundance,pointer_1based=offset+1,
        record_sha256=hashlib.sha256(raw).hexdigest(),tau=tau,departures=populations)
