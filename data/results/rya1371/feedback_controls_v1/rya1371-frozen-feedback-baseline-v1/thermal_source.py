"""Opt-in stationary fixed-population thermal opacity/emissivity diagnostic.

Retains original absorption strengths and photon-unit rate normalization. This
is a Kirchhoff-matched sensitivity control, not a complete native operator.
"""
import numpy as np


def thermal_correction(Vij,g_old,bare_old,star_ratio,nupper,wave_nm,T,h,k,c):
    """Return delta chi [m^-1], delta eta [SI/m] and new inverse factors."""
    V=np.asarray(Vij,float);g=np.asarray(g_old,float);bare=np.asarray(bare_old,float)
    ratio=np.asarray(star_ratio,float);upper=np.asarray(nupper,float);wave=np.asarray(wave_nm,float);T=np.asarray(T,float)
    shape=(len(wave),len(T))
    if any(x.shape!=shape for x in [V,g,bare])or ratio.shape!=T.shape or upper.shape!=T.shape:
        raise ValueError('Thermal coefficient shapes differ')
    if any(not np.isfinite(x).all()or np.any(x<=0)for x in [g,bare,ratio,upper,wave,T])or not np.isfinite(V).all()or np.any(V<0)or any(not np.isfinite(x)or x<=0 for x in [h,k,c]):
        raise ValueError('Invalid thermal coefficients/state')
    nu=c/(wave[:,None]*1e-9);beta=np.exp(-h*nu/(k*T[None,:]));new_g=ratio[None,:]*beta;new_bare=2*h*nu**3/c**2
    return upper[None,:]*V*(g-new_g),upper[None,:]*V*(new_g*new_bare-g*bare),new_g,new_bare


def transition_fields(ctx,transition,h,k,c):
    import lightweaver.constants as lw
    a=ctx.activeAtoms[0];wave=np.asarray(transition.wavelength);T=np.asarray(ctx.atmos.pyAtmos.temperature)
    i=transition.i;j=transition.j;shape=(len(wave),len(T))
    if transition.type=='Line':
        phi=np.asarray(transition.phi);base=phi[:,0,0,:]
        np.testing.assert_array_equal(phi,np.broadcast_to(base[:,None,None,:],phi.shape))
        model=transition.transModel
        V=lw.HPlanck*lw.CLight/(4*np.pi)*(model.lambda0/wave[:,None])*model.Bij*base
        g=np.full(shape,model.Bji/model.Bij);bare=np.full(shape,model.Aji/model.Bji)
    elif transition.type=='Continuum':
        nu=lw.CLight/(wave[:,None]*1e-9)
        V=np.broadcast_to(np.asarray(transition.alpha)[:,None],shape)
        g=a.nStar[i][None,:]/a.nStar[j][None,:]*np.exp(-lw.HPlanck*nu/(lw.KBoltzmann*T[None,:]))
        bare=np.broadcast_to(2*lw.HPlanck*nu**3/lw.CLight**2,shape)
    else:raise ValueError('Unknown source transition')
    return thermal_correction(V,g,bare,a.nStar[i]/a.nStar[j],a.n[j],wave,T,h,k,c)


def apply_control(ctx,components,constants):
    import copy
    import lightweaver as lw
    from lightweaver.LwCompiled import BackgroundProvider
    if components not in [[],['BB'],['BF'],['BB','BF']]:raise ValueError('Unknown explicit thermal components')
    if np.any(ctx.atmos.pyAtmos.vlos!=0):raise ValueError('Thermal diagnostic requires stationary profiles')
    wave=np.asarray(ctx.spect.wavelength);delta_chi=np.zeros_like(ctx.background.chi);delta_eta=np.zeros_like(ctx.background.eta);base_chi=np.array(ctx.background.chi).copy();base_eta=np.array(ctx.background.eta).copy();counts={'BB':0,'BF':0}
    for t in ctx.activeAtoms[0].trans:
        kind='BB'if t.type=='Line'else'BF'
        local=np.asarray(t.wavelength);idx=np.searchsorted(wave,local);np.testing.assert_array_equal(wave[idx],local)
        a=ctx.activeAtoms[0];i=t.i;j=t.j;T=np.asarray(ctx.atmos.pyAtmos.temperature)
        import lightweaver.constants as c
        if kind=='BB':
            p=np.asarray(t.phi);profile=p[:,0,0,:];np.testing.assert_array_equal(p,np.broadcast_to(profile[:,None,None,:],p.shape));m=t.transModel
            V=c.HPlanck*c.CLight/(4*np.pi)*(m.lambda0/local[:,None])*m.Bij*profile
            g=np.full_like(V,m.Bji/m.Bij);bare=np.full_like(V,m.Aji/m.Bji)
        else:
            nu=c.CLight/(local[:,None]*1e-9);V=np.broadcast_to(np.asarray(t.alpha)[:,None],(len(local),len(T)))
            g=a.nStar[i][None,:]/a.nStar[j][None,:]*np.exp(-c.HPlanck*nu/(c.KBoltzmann*T[None,:]));bare=np.broadcast_to(2*c.HPlanck*nu**3/c.CLight**2,V.shape)
        base_chi[idx]+=a.n[i][None,:]*V-a.n[j][None,:]*g*V;base_eta[idx]+=a.n[j][None,:]*g*bare*V
        if kind in components:
            dc,de,_,_=transition_fields(ctx,t,**constants);delta_chi[idx]+=dc;delta_eta[idx]+=de
        counts[kind]+=1
    chi=base_chi+delta_chi;eta=base_eta+delta_eta;sca=np.array(ctx.background.sca).copy()
    if any(not np.isfinite(x).all()or np.any(x<0)for x in [chi,eta,sca]):raise ValueError('Invalid complete thermal coefficients; no clipping allowed')
    if counts!={'BB':256,'BF':22}:raise ValueError('Incomplete fixed-opacity atom')
    for x in [chi,eta,sca]:x.flags.writeable=False
    original_state={name:np.array(getattr(ctx.atmos.pyAtmos,name)).copy()for name in ['height','temperature','ne','nHTot','vturb','vlos']}
    class FrozenProvider(BackgroundProvider):
        def compute_background(self,atmos,out_chi,out_eta,out_sca):
            for name,values in original_state.items():np.testing.assert_array_equal(getattr(atmos.pyAtmos,name),values)
            for out,value in [(out_chi,chi),(out_eta,eta),(out_sca,sca)]:
                if out.shape!=value.shape:raise ValueError('Changed frozen coefficient grid')
                out[:]=value
    spect=copy.deepcopy(ctx.kwargs['spect']);spect.radSet.set_passive('O')
    formal=lw.Context(copy.deepcopy(ctx.atmos.pyAtmos),spect,ctx.eqPops,Nthreads=1,formalSolver='piecewise_linear_1d',conserveCharge=False,backgroundProvider=FrozenProvider)
    if formal.activeAtoms or formal.detailedAtoms:raise ValueError('Frozen total provider must not double-count any atom')
    np.testing.assert_array_equal(formal.spect.wavelength,wave);formal.depthData.fill=True
    return formal,dict(components=components,counts=counts,constants_SI=constants,max_abs_delta_chi_m_inv=float(np.max(abs(delta_chi))),max_abs_delta_eta=float(np.max(abs(delta_eta))),scope='Complete positive frozen total opacity/emissivity provider in a separate passive-only formal context; exact source populations/profiles retained separately for rates. Original source strengths, provider packet/scattering and photon-unit normalization unchanged. No SE/preconditioner or production qualification.')


def physical_rates(ctx,components,constants,Jnu=None):
    import lightweaver.constants as lw
    from pipeline.absolute_rate_integrals import bb_rates,bf_rates,trapezoid_measure
    a=ctx.activeAtoms[0];R=np.zeros((41,41,len(ctx.atmos.pyAtmos.temperature)));wave=np.asarray(ctx.spect.wavelength);W=float(np.sum(ctx.atmos.pyAtmos.wmu));T=np.asarray(ctx.atmos.pyAtmos.temperature)
    for t in a.trans:
        local=np.asarray(t.wavelength);idx=np.searchsorted(wave,local);J=np.asarray(ctx.spect.J if Jnu is None else Jnu)[idx];i=t.i;j=t.j;kind='BB'if t.type=='Line'else'BF'
        if kind=='BB':
            model=t.transModel;profile=np.asarray(t.phi)[:,0,0,:]
            up,down=bb_rates(local,profile,J,lambda0_nm=model.lambda0,Aji=model.Aji,Bij=model.Bij,Bji=model.Bji,angular_sum=W)
            if kind in components:
                _,_,g,bare=transition_fields(ctx,t,**constants)
                measure=trapezoid_measure(local)[:,None];norm=np.sum(measure*profile,axis=0)*W;weight=measure*profile/norm[None,:]*(model.lambda0/local[:,None])
                down=model.Bij*np.sum(weight*g*(bare*W+J),axis=0)
        else:
            values=dict(h=lw.HPlanck,k=lw.KBoltzmann,c=lw.CLight)
            if kind in components:values=constants
            up,down=bf_rates(local,np.asarray(t.alpha),J,nstar_ratio=a.nStar[i]/a.nStar[j],T_K=T,angular_sum=W,**values)
            if kind in components:
                # Preserve original LW photon-unit normalization; do not fold
                # an extra 1/h scale change into the field-feedback experiment.
                up*=constants['h']/lw.HPlanck;down*=constants['h']/lw.HPlanck
        R[j,i]+=up;R[i,j]+=down
    if not np.isfinite(R).all()or np.any(R<0):raise ValueError('Invalid modified physical rates')
    return R
