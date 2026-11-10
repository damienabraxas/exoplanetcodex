"""Opt-in natural-Voigt Fe I trace lines on an unchanged bound background.

Fixed-atmosphere LTE Kurucz partition/ionization diagnostic; not a complete
pressure-broadened line list, atmosphere/EOS update or Fe abundance measurement.
"""
from functools import partial
import numpy as np


def validate_lines(lines):
    required={'line_id','lambda0_vac_nm','excitation_eV','loggf','gamma_rad_s_inv'}
    seen=set()
    for row in lines:
        if set(row)!=required or row['line_id'] in seen:
            raise ValueError('Invalid or duplicate trace Fe line')
        seen.add(row['line_id'])
        values=[row[k] for k in required if k!='line_id']
        if (not np.isfinite(values).all() or row['lambda0_vac_nm']<=0
                or row['excitation_eV']<0 or row['gamma_rad_s_inv']<=0):
            raise ValueError('Missing/nonphysical trace Fe data')


def trace_fe_factory(base_factory, lines, support_vac_nm):
    validate_lines(lines)
    lo,hi=support_vac_nm
    if not np.isfinite([lo,hi]).all() or not 0<lo<hi:
        raise ValueError('Invalid explicit trace support')
    if any(not lo<r['lambda0_vac_nm']<hi for r in lines):
        raise ValueError('Line centre outside trace support')
    return partial(TraceFeIBackground,base_factory=base_factory,lines=tuple(dict(r)for r in lines),support=(lo,hi))


# Importing this module is intentionally confined to the installed LW environment.
from lightweaver.LwCompiled import BackgroundProvider


class TraceFeIBackground(BackgroundProvider):
    def __init__(self,eqpops,radset,wavelength,*,base_factory,lines,support):
        from lightweaver.atomic_table import PeriodicTable,KuruczPfTable
        super().__init__(eqpops,radset,wavelength)
        validate_lines(lines)
        self.base=base_factory(eqpops,radset,wavelength)
        self.wave=np.asarray(wavelength);self.lines=lines
        self.mask=(self.wave>=support[0])&(self.wave<=support[1])
        self.eq=eqpops;self.model=radset.atoms[PeriodicTable['Fe']]
        if self.model not in radset.passiveAtoms:
            raise ValueError('Trace Fe requires passive Fe; refuse double-installed lines')
        self.partition=KuruczPfTable(atomicAbundance=radset.abundance)[PeriodicTable['Fe']]
        self.mass=PeriodicTable['Fe'].mass
        self.cache=None;self.neutral_ratio_range=None

    def compute_background(self,atmos,chi,eta,sca):
        # Delegate every time: the complete bound-input/cap/composition checks remain live.
        self.base.compute_background(atmos,chi,eta,sca)
        if not self.lines or not self.mask.any():return
        a=atmos.pyAtmos
        actual=[np.asarray(getattr(a,k)) for k in ('temperature','ne','nHTot','vturb')]
        if self.cache is None:
            from pipeline.bride_lte_blend import lte_line_coefficients
            if a.temperature.min()<self.partition.Tpf[0] or a.temperature.max()>self.partition.Tpf[-1]:
                raise ValueError('Outside Fe partition temperature domain')
            ions=self.partition.lte_ionisation(a)
            U=np.exp(np.interp(a.temperature,self.partition.Tpf,self.partition.pf[0]))
            neutral=[l.stage==0 for l in self.model.levels]
            native=self.eq.atomicPops[self.model.element].n[neutral].sum(axis=0)
            ratio=ions[0]/native
            if not np.isfinite(ratio).all() or np.any(ratio<=0):raise ValueError('Invalid Fe LTE reservoir')
            self.neutral_ratio_range=[float(ratio.min()),float(ratio.max())]
            extinction=np.zeros((self.mask.sum(),len(a.temperature)));emission=np.zeros_like(extinction)
            for line in self.lines:
                x,e=lte_line_coefficients(self.wave[self.mask],a.temperature,ions[0]/U,a.vturb,
                    lambda0_nm=line['lambda0_vac_nm'],mass_amu=self.mass,
                    excitation_eV=line['excitation_eV'],loggf=line['loggf'],gamma_s_inv=self.line_width(line,a))
                extinction+=x;emission+=e
            self.cache=([v.copy() for v in actual],extinction,emission)
        if any(not np.array_equal(x,y) for x,y in zip(actual,self.cache[0])):
            raise ValueError('Cached trace coefficients refuse changed local fields')
        chi[self.mask]+=self.cache[1];eta[self.mask]+=self.cache[2]
        if not np.isfinite(chi).all() or not np.isfinite(eta).all():raise ValueError('Nonfinite mixed opacity')

    def line_width(self,line,atmosphere):
        return line["gamma_rad_s_inv"]
