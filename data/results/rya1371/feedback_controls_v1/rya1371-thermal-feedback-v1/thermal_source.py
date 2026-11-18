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
    from pipeline.reference_engines import fingerprint
    if components not in [['BB'],['BF'],['BB','BF']]:raise ValueError('Unknown explicit thermal components')
    if np.any(ctx.atmos.pyAtmos.vlos!=0):raise ValueError('Thermal diagnostic requires stationary profiles')
    wave=np.asarray(ctx.spect.wavelength);delta_chi=np.zeros_like(ctx.background.chi);delta_eta=np.zeros_like(ctx.background.eta);counts={'BB':0,'BF':0}
    for t in ctx.activeAtoms[0].trans:
        kind='BB'if t.type=='Line'else'BF'
        if kind not in components:continue
        local=np.asarray(t.wavelength);idx=np.searchsorted(wave,local);np.testing.assert_array_equal(wave[idx],local)
        dc,de,_,_=transition_fields(ctx,t,**constants);delta_chi[idx]+=dc;delta_eta[idx]+=de;counts[kind]+=1
    chi=np.asarray(ctx.background.chi);eta=np.asarray(ctx.background.eta)
    if any(not np.isfinite(x).all()for x in [delta_chi,delta_eta])or np.any(chi+delta_chi<0)or np.any(eta+delta_eta<0):raise ValueError('Thermal correction gives invalid background storage; no clipping allowed')
    chi[:]+=delta_chi;eta[:]+=delta_eta
    return dict(components=components,counts=counts,constants_SI=constants,max_abs_delta_chi_m_inv=float(np.max(abs(delta_chi))),max_abs_delta_eta=float(np.max(abs(delta_eta))),scope='Fixed complete populations only; additive process-local storage represents changes to active O inverse/emission coefficients. Original absorption strengths, background provider packet, scattering and photon-unit normalization unchanged. No SE/preconditioner or production qualification.')


def physical_rates(ctx,components,constants):
    import lightweaver.constants as lw
    from pipeline.absolute_rate_integrals import bb_rates,bf_rates,trapezoid_measure
    a=ctx.activeAtoms[0];R=np.zeros((41,41,len(ctx.atmos.pyAtmos.temperature)));wave=np.asarray(ctx.spect.wavelength);W=float(np.sum(ctx.atmos.pyAtmos.wmu));T=np.asarray(ctx.atmos.pyAtmos.temperature)
    for t in a.trans:
        local=np.asarray(t.wavelength);idx=np.searchsorted(wave,local);J=np.asarray(ctx.spect.J)[idx];i=t.i;j=t.j;kind='BB'if t.type=='Line'else'BF'
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
