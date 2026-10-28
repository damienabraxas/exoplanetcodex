"""Direct calls to a pinned local MULTI TAUTSP/PPVALU diagnostic oracle."""
import ctypes as ct
import json
from pathlib import Path
import threading
import numpy as np

_LOCK=threading.Lock()  # Native INTERV has SAVE state; serialize oracle calls.


class NativeTautsp:
    def __init__(self, directory):
        from pipeline.bride_reference import fingerprint
        directory=Path(directory)
        manifest=json.loads((directory/'manifest.json').read_text())
        if manifest['source_commit']!='4d6aa3c39d634b38da772e469a5692435eba1354':raise ValueError('Wrong native source')
        library=directory/'libmulti_spline.dylib'
        if fingerprint(library)!=manifest['library']:raise ValueError('Changed spline oracle library')
        self.receipt=fingerprint(directory/'manifest.json')
        self.library=ct.CDLL(str(library))
        d=ct.POINTER(ct.c_double);i=ct.POINTER(ct.c_int)
        self.library.tautsp_.argtypes=[d,d,i,d,d,d,d,i,i,i]
        self.library.ppvalu_.argtypes=[d,d,i,i,d,i,i]
        self.library.ppvalu_.restype=ct.c_double

    def interpolate(self, grid, values, temperatures, *, logarithmic, allow_negative=False):
        grid=np.ascontiguousarray(grid,dtype=np.float64)
        values=np.ascontiguousarray(values,dtype=np.float64)
        query=np.ascontiguousarray(temperatures,dtype=np.float64)
        if (grid.ndim!=1 or query.ndim!=1 or not query.size or values.shape!=grid.shape or len(grid)<4
                or not np.isfinite(grid).all() or not np.isfinite(values).all()
                or not np.isfinite(query).all() or np.any(np.diff(grid)<=0)
                or query.min()<grid[0] or query.max()>grid[-1]):
            raise ValueError('Native spline requires finite in-domain data and >=4 ordered nodes')
        if logarithmic:
            if np.any(values<=0):raise ValueError('Native log interpolation requires positive coefficients')
            values=np.log(values)
        n=ct.c_int(len(grid));gamma=ct.c_double(5.5);pieces=ct.c_int();order=ct.c_int();flag=ct.c_int()
        work=np.zeros((len(grid),6),order='F');breaks=np.zeros(2*len(grid)+1)
        coef=np.zeros((4,2*len(grid)+1),order='F')
        ptr=lambda a:a.ctypes.data_as(ct.POINTER(ct.c_double))
        with _LOCK:
            self.library.tautsp_(ptr(grid),ptr(values),ct.byref(n),ct.byref(gamma),ptr(work),ptr(breaks),ptr(coef),ct.byref(pieces),ct.byref(order),ct.byref(flag))
            if flag.value!=1 or order.value!=4:raise ValueError('Native TAUTSP construction failed')
            output=np.empty_like(query)
            for index,T in np.ndenumerate(query):
                x=ct.c_double(T);derivative=ct.c_int(0);status=ct.c_int()
                output[index]=self.library.ppvalu_(ptr(breaks),ptr(coef),ct.byref(pieces),ct.byref(order),ct.byref(x),ct.byref(derivative),ct.byref(status))
                if status.value!=0:raise ValueError('Native PPVALU rejected query')
        if logarithmic:output=np.exp(output)
        if not np.isfinite(output).all() or (not allow_negative and np.any(output<0)):raise ValueError('Invalid native interpolated coefficient')
        return output


def install_native_collision_control(ctx, source, oracle, *, source_directions, reversible_ch0=False):
    """Precompute exact held TAUTSP coefficients, retaining local LW LTE/H fields.

    Source semantics follow mul23_opacu.f GENCL; alternate directions isolate
    interpolation from RH/LW direction conventions. No native EOS parity implied.
    """
    from lightweaver.collisional_rates import CollisionalRates
    import lightweaver.constants as c
    if len(ctx.activeAtoms)!=1:raise ValueError('One active atom required')
    if reversible_ch0 and not source_directions:raise ValueError('Reverse CH0 hypothesis requires explicit source direction')
    atom=ctx.activeAtoms[0].atomicModel
    T=ctx.atmos.pyAtmos.temperature.copy()
    report=dict(oracle=oracle.receipt,source_directions=source_directions,reversible_ch0_LTE_reservoir_hypothesis=reversible_ch0,
        temperature_unchanged_required=True,CH_negative_spline_values_clipped_as_native=0)
    replacement=[]
    for raw in source.collisions:
        if raw.kind not in ('CE','CI','CH','CH0'):raise ValueError('Unsupported native collision kind')
        grid=raw.temperature
        coefficient=oracle.interpolate(grid,raw.rates,T,logarithmic=raw.kind in ('CE','CI','CH0'),allow_negative=raw.kind=='CH')
        if raw.kind=='CH':
            report['CH_negative_spline_values_clipped_as_native']+=int(np.count_nonzero(coefficient<0))
            coefficient=np.maximum(coefficient,0.)
        coefficient=coefficient*1e-6

        class NativeRate(CollisionalRates):
            def setup(self,model):self.atom=model
            def compute_rates(self,atmos,eqpops,matrix):
                if not np.array_equal(atmos.temperature,T):raise ValueError('Native rate control atmosphere changed')
                lower,upper=self.i,self.j;q=self.coefficient
                nstar=eqpops.atomicPops[self.atom.element].nStar
                if self.source_kind=='CE':
                    down=q*atmos.ne*np.sqrt(T)*self.atom.levels[lower].g/self.atom.levels[upper].g
                    up=down*nstar[upper]/nstar[lower]
                elif self.source_kind=='CI':
                    gap=(self.atom.levels[upper].E-self.atom.levels[lower].E)*c.HC_CM
                    up=q*atmos.ne*np.sqrt(T)*np.exp(-gap/(c.KBoltzmann*T))
                    down=up*nstar[lower]/nstar[upper]
                elif self.source_kind=='CH':
                    if source_directions:
                        down=q*eqpops['H'][0];up=down*nstar[upper]/nstar[lower]
                    else:
                        up=q*eqpops['H'][0];down=up*nstar[lower]/nstar[upper]
                else:
                    if source_directions:
                        matrix[self.raw_destination,self.raw_source]+=q*eqpops['H'][0]
                        if reversible_ch0:
                            matrix[self.raw_source,self.raw_destination]+=q*eqpops['H'][0]*nstar[self.raw_source]/nstar[self.raw_destination]
                    else:matrix[lower,upper]+=q*eqpops['H'][0]
                    return
                matrix[lower,upper]+=down;matrix[upper,lower]+=up

        rate=NativeRate(i=min(raw.i,raw.j)-1,j=max(raw.i,raw.j)-1)
        rate.source_kind=raw.kind;rate.coefficient=coefficient
        rate.raw_source=raw.j-1;rate.raw_destination=raw.i-1;rate.setup(atom)
        if raw.kind=='CH0' and reversible_ch0:
            if atom.levels[rate.raw_destination].stage!=atom.levels[rate.raw_source].stage+1:
                raise ValueError('Reverse CH0 hypothesis requires explicit adjacent-ion-stage forward process')
        replacement.append(rate)
    atom.collisions=replacement
    return report
