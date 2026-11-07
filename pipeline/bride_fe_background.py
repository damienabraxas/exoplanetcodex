"""Opt-in bounded Fe I LTE-background intervention at one verified native node.

Replace, never add twice, passive Fe I bound-free extinction/emission only.
The packet is source-derived; no fitted opacity multiplier is accepted.
"""
from functools import partial
import json
from pathlib import Path

import numpy as np
from lightweaver import constants as c
from lightweaver.atomic_table import PeriodicTable
from lightweaver.LwCompiled import BackgroundProvider,FastBackground


def load_packet(directory):
    from pipeline.bride_reference import fingerprint
    directory=Path(directory);manifest=json.loads((directory/'manifest.json').read_text())
    if manifest['schema']!='codex.bride_fe_packet/1':raise ValueError('Wrong Fe packet schema')
    path=directory/manifest['array']['name']
    if fingerprint(path)!=manifest['array']:raise ValueError('Changed Fe packet')
    data=dict(np.load(path));data['receipt']=fingerprint(directory/'manifest.json')
    wave=data['wavelength_vac_nm'];n=len(data['expected_temperature_K'])
    if (wave.ndim!=1 or len(wave)<2 or np.any(np.diff(wave)<=0)
            or data['opacity_m2_kg'].shape!=(len(wave),n)
            or data['opacity_net_per_neutral_Fe_m2'].shape!=(len(wave),n)
            or any(not np.isfinite(v).all() for k,v in data.items()if isinstance(v,np.ndarray))
            or np.any(data['opacity_m2_kg']<0) or np.any(data['opacity_net_per_neutral_Fe_m2']<0)):
        raise ValueError('Invalid Fe packet arrays')
    return data


def fe_factory(directory,normalization):
    if normalization not in ('raw_specific','local_neutral'):
        raise ValueError('Declare raw_specific or local_neutral Fe normalization')
    return partial(FeIBackground,packet=load_packet(directory),normalization=normalization)


class FeIBackground(BackgroundProvider):
    def __init__(self,eqPops,radSet,wavelength,*,packet,normalization):
        super().__init__(eqPops,radSet,wavelength)
        self.base=FastBackground(eqPops,radSet,wavelength,Nthreads=1)
        self.eq=eqPops;self.radset=radSet;self.wave=np.asarray(wavelength)
        self.packet=packet;self.normalization=normalization
        self.model=radSet.atoms[PeriodicTable[26]]
        if self.model not in radSet.passiveAtoms:raise ValueError('Fe intervention requires passive LTE Fe')
        self.mask=(self.wave>=packet['wavelength_vac_nm'][0])&(self.wave<=packet['wavelength_vac_nm'][-1])

    def compute_background(self,atmos,chi,eta,sca):
        py=atmos.pyAtmos;p=self.packet
        for actual,key in ((py.temperature,'expected_temperature_K'),(py.ne,'expected_ne_m3'),(py.nHTot,'expected_nH_m3')):
            if not np.array_equal(actual,p[key]):raise ValueError('Fe packet refuses changed T/ne/nH or depth interpolation')
        self.base.compute_background(atmos,chi,eta,sca)
        if not self.mask.any():return
        waves=self.wave[self.mask];population=self.eq.atomicPops[self.model.element].n
        thermal=-np.expm1(-c.HC/(c.KBoltzmann*1e-9*waves[:,None]*py.temperature[None,:]))
        old=sum((t.alpha(waves)[:,None]*thermal*population[t.i][None,:]
                 for t in self.model.continua if self.model.levels[t.i].stage==0),np.zeros_like(thermal))
        field='opacity_m2_kg' if self.normalization=='raw_specific'else 'opacity_net_per_neutral_Fe_m2'
        new=np.column_stack([np.interp(waves,p['wavelength_vac_nm'],p[field][:,d]) for d in range(len(py.temperature))])
        if self.normalization=='raw_specific':new*=p['expected_rho_kg_m3'][None,:]
        elif self.normalization=='local_neutral':
            neutral=np.array([l.stage==0 for l in self.model.levels]);new*=population[neutral].sum(axis=0)[None,:]
        else:raise ValueError('Unsupported Fe normalization')
        planck=2*c.HPlanck*c.CLight/(waves[:,None]*1e-9)**3/np.expm1(c.HC/(waves[:,None]*1e-9*c.KBoltzmann*py.temperature[None,:]))
        chi[self.mask]+=new-old;eta[self.mask]+=(new-old)*planck
        if (not np.isfinite(chi).all() or not np.isfinite(eta).all()
                or np.any(chi<sca) or np.any(eta<0)):
            raise ValueError('Invalid Fe replacement coefficients')
