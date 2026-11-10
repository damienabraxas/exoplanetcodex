"""Opt-in one-line natural-only LTE FeII intervention on an unchanged provider."""
from functools import partial
import numpy as np
from lightweaver.LwCompiled import BackgroundProvider
from pipeline.bride_trace_fe_lines import validate_lines


def feii_factory(base_factory,line,support):
    lo,hi=support
    if not np.isfinite(support).all() or not 0<lo<hi:raise ValueError('Invalid FeII support')
    if line is not None:
        validate_lines([line])
        if not lo<line['lambda0_vac_nm']<hi:raise ValueError('FeII centre outside support')
    return partial(FeIIControl,base_factory=base_factory,line=None if line is None else dict(line),support=tuple(support))


class FeIIControl(BackgroundProvider):
    def __init__(self,eqpops,radset,wavelength,*,base_factory,line,support):
        from lightweaver.atomic_table import PeriodicTable,KuruczPfTable
        super().__init__(eqpops,radset,wavelength)
        self.base=base_factory(eqpops,radset,wavelength);self.line=line;self.wave=np.asarray(wavelength)
        self.mask=(self.wave>=support[0])&(self.wave<=support[1]);self.cache=None
        atom=radset.atoms[PeriodicTable['Fe']]
        if atom not in radset.passiveAtoms:raise ValueError('FeII control requires passive Fe')
        if any(support[0]<=t.lambda0<=support[1] for t in atom.lines):raise ValueError('Native Fe line overlaps added support')
        self.partition=KuruczPfTable(atomicAbundance=radset.abundance)[PeriodicTable['Fe']]
        self.mass=PeriodicTable['Fe'].mass

    def line_width(self,atmosphere):
        return self.line['gamma_rad_s_inv']

    def compute_background(self,atmos,chi,eta,sca):
        self.base.compute_background(atmos,chi,eta,sca)
        if self.line is None or not self.mask.any():return
        a=atmos.pyAtmos;fields=[np.asarray(getattr(a,k)) for k in ('temperature','ne','nHTot','vturb')]
        if self.cache is None:
            from pipeline.bride_lte_blend import lte_line_coefficients
            p=self.partition
            if a.temperature.min()<p.Tpf[0] or a.temperature.max()>p.Tpf[-1]:raise ValueError('FeII partition domain')
            ions=p.lte_ionisation(a);U=np.exp(np.interp(a.temperature,p.Tpf,p.pf[1]))
            l=self.line
            x,e=lte_line_coefficients(self.wave[self.mask],a.temperature,ions[1]/U,a.vturb,
                lambda0_nm=l['lambda0_vac_nm'],mass_amu=self.mass,excitation_eV=l['excitation_eV'],
                loggf=l['loggf'],gamma_s_inv=self.line_width(a))
            self.cache=([v.copy() for v in fields],x,e)
        if any(not np.array_equal(a,b) for a,b in zip(fields,self.cache[0])):raise ValueError('FeII cached fields changed')
        chi[self.mask]+=self.cache[1];eta[self.mask]+=self.cache[2]
        if not np.isfinite(chi).all() or not np.isfinite(eta).all():raise ValueError('Invalid FeII coefficients')
