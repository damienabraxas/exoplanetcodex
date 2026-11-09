"""Explicit legacy-GENCOL diagnostic conventions; never automatic/adopted fixes.

Inspected GENCOL source at FORWARD commit15145239: CH is downward, CH0 input
order is source->destination, CH0 interpolates log(rate) linearly in T.
RH collision.c and modelatom.c use other conventions. Exact grid-build parity
is not established. Keep the existing converter as a named negative control.
"""
import numpy as np


def install_gencol_control(ctx, source, *, change_ch):
    """Replace only CH0 and optionally CH in a fresh context, before solving."""
    from lightweaver.collisional_rates import TemperatureInterpolationRates
    from weno4 import weno4

    class DirectedCH0(TemperatureInterpolationRates):
        def setup(self, atom):
            # Deliberately do NOT call sorting setup: j is source, i destination.
            self.atom=atom
            self.temperature=np.asarray(self.temperature)
            self.rates=np.asarray(self.rates)
            if np.any(self.rates<=0):raise ValueError('Log CH0 requires positive coefficients')

        def compute_rates(self, atmos, eqpops, matrix):
            if atmos.temperature.min()<self.temperature[0] or atmos.temperature.max()>self.temperature[-1]:
                raise ValueError('CH0 outside supplied temperature domain')
            coefficient=np.exp(np.interp(atmos.temperature,self.temperature,np.log(self.rates)))
            matrix[self.i,self.j]+=coefficient*eqpops['H'][0]

    class DownwardCH(TemperatureInterpolationRates):
        def compute_rates(self, atmos, eqpops, matrix):
            coefficient=weno4(atmos.temperature,self.temperature,self.rates)
            if np.any(coefficient<0):raise ValueError('Negative CH interpolation')
            down=coefficient*eqpops['H'][0]
            nstar=eqpops.atomicPops[self.atom.element].nStar
            matrix[self.i,self.j]+=down
            matrix[self.j,self.i]+=down*nstar[self.j]/nstar[self.i]

    if len(ctx.activeAtoms)!=1:raise ValueError('Expected one active atom')
    model=ctx.activeAtoms[0].atomicModel
    if len(model.collisions)!=len(source.collisions):raise ValueError('Collision layout differs')
    replacement=[]
    for raw,current in zip(source.collisions,model.collisions):
        if raw.kind=='CH0':
            # Raw first/second indices preserved by parser as j/i respectively.
            rate=DirectedCH0(j=raw.j-1,i=raw.i-1,
                temperature=raw.temperature,rates=raw.rates*1e-6)
            rate.setup(model);replacement.append(rate)
        elif raw.kind=='CH' and change_ch:
            rate=DownwardCH(j=current.j,i=current.i,temperature=current.temperature,rates=current.rates.copy())
            rate.setup(model);replacement.append(rate)
        else:replacement.append(current)
    model.collisions=replacement
