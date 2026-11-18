"""Independent stationary total chi/thermal-eta reconstruction; no formal solve."""
import numpy as np


def reconstruct_coefficients(ctx,oxygen):
    import lightweaver.constants as c
    wave=np.asarray(ctx.spect.wavelength);T=np.asarray(ctx.atmos.pyAtmos.temperature)
    chi=np.array(ctx.background.chi).copy();eta=np.array(ctx.background.eta).copy();counts={'BB':0,'BF':0}
    for t in oxygen.trans:
        local=np.asarray(t.wavelength);index=np.searchsorted(wave,local);np.testing.assert_array_equal(wave[index],local)
        if t.type=='Line':
            p=np.asarray(t.phi);profile=p[:,0,0,:]
            np.testing.assert_array_equal(p,np.broadcast_to(profile[:,None,None,:],p.shape))
            m=t.transModel;V=c.HPlanck*c.CLight/(4*np.pi)*(m.lambda0/local[:,None])*m.Bij*profile
            g=m.Bji/m.Bij;bare=m.Aji/m.Bji;kind='BB'
        elif t.type=='Continuum':
            nu=c.CLight/(local[:,None]*1e-9);V=np.asarray(t.alpha)[:,None]
            g=oxygen.nStar[t.i][None,:]/oxygen.nStar[t.j][None,:]*np.exp(-c.HPlanck*nu/(c.KBoltzmann*T[None,:]));bare=2*c.HPlanck*nu**3/c.CLight**2;kind='BF'
        else:raise ValueError('Unknown coefficient transition')
        chi[index]+=oxygen.n[t.i][None,:]*V-oxygen.n[t.j][None,:]*g*V
        eta[index]+=oxygen.n[t.j][None,:]*g*bare*V
        counts[kind]+=1
    if counts!={'BB':256,'BF':22}or any(not np.isfinite(x).all()or np.any(x<=0)for x in [chi,eta]):raise ValueError('Invalid/incomplete reconstructed whole-O coefficients')
    return chi,eta,counts
