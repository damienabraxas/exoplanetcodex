"""Pinned neutral-H2 table adapter and fixed-ne hydrogen-nuclei accounting.

Pressure equilibrium convention Kp=p(H I)^2/p(H2), in Pa. No extrapolation,
charged-molecule network, metal/CNO chemistry, or complete EOS is supplied.
"""
import json
from pathlib import Path
import numpy as np
from scipy.constants import k as KB, e as EV
from scipy.interpolate import CubicSpline


def read_wide_table(path,kind):
    lines=Path(path).read_text().splitlines()
    expected={'partition':'Partition functions Q','equilibrium':'Equilibrium constants  log10(pK)'}[kind]
    if lines[0].strip()!=expected:raise ValueError('Wrong published table kind')
    n=int(lines[1]);header=lines[2].split()
    if header[:2]!=['T','[K]']:raise ValueError('Missing explicit K temperature header')
    T=np.array(list(map(float,header[2:])))
    if T.shape!=(42,)or not np.isfinite(T).all()or np.any(np.diff(T)<=0)or np.any(T<=0):raise ValueError('Invalid published grid')
    rows={}
    for line in lines[3:]:
        if not line.strip():continue
        fields=line.split();name=fields[0]
        if name in rows or len(fields)!=43:raise ValueError('Ambiguous/incomplete species row')
        values=np.array(list(map(float,fields[1:])))
        if not np.isfinite(values).all()or(kind=='partition'and np.any(values<=0)):raise ValueError('Invalid tabulated values')
        rows[name]=values
    if len(rows)!=n:raise ValueError('Published species count mismatch')
    return T,rows


def h2_dissociation_eV(path):
    rows=[l.split()for l in Path(path).read_text().splitlines()if l.split()and l.split()[0]=='H2']
    if len(rows)!=1 or rows[0][1:3]!=['H','H']:raise ValueError('Missing unique neutral-H2 reaction')
    value=float(rows[0][-2])
    if not np.isfinite(value)or value<=0:raise ValueError('Missing adopted dissociation energy')
    return value


class H2Tables:
    """Experimental 1000-10000 K adapter; exact source knots, no observer fit."""
    def __init__(self,directory):
        from pipeline.bride_reference import fingerprint
        directory=Path(directory);receipt=json.loads((directory/'receipt.json').read_text())
        if receipt.get('schema')!='codex.bride_molecular_source/1':raise ValueError('Wrong source receipt schema')
        if receipt['source_commit']!='cd2f54ed46687659017b343763e95478ed93653e':raise ValueError('Wrong source revision')
        required={'ReadMe','ReadMe_Bug_fix_Oct-Nov2022','table1.dat','table6.dat','table7.dat','table7_vNov2022.dat','table8_vNov2022.dat','LICENSE'}
        if set(receipt['files'])!=required:raise ValueError('Incomplete source receipt')
        for name,entry in receipt['files'].items():
            if fingerprint(directory/name)!=entry['fingerprint']:raise ValueError('Changed chemistry source '+name)
        T,K=read_wide_table(directory/'table7_vNov2022.dat','equilibrium')
        tq,Q=read_wide_table(directory/'table6.dat','partition');ta,A=read_wide_table(directory/'table8_vNov2022.dat','partition')
        np.testing.assert_array_equal(T,tq);np.testing.assert_array_equal(T,ta)
        self.receipt=fingerprint(directory/'receipt.json');self.dissociation_eV=h2_dissociation_eV(directory/'table1.dat')
        keep=T>=1000.;self.temperature=T[keep];self.log10_Kp_Pa=K['H2'][keep]
        self.partition_H2=Q['H2'][keep];self.partition_HI=A['H_I'][keep]
        self.domain=(float(self.temperature[0]),float(self.temperature[-1]))
        self.log_phi_knots=np.log(KB*self.temperature)-np.log(10)*self.log10_Kp_Pa
        # Paper section 3 prescription; boundary convention is an explicit operator choice.
        self.trend=self.log10_Kp_Pa+self.dissociation_eV*5040/self.temperature
        self.interpolator=CubicSpline(np.log10(self.temperature),self.trend,bc_type='not-a-knot',extrapolate=False)
        for a in (self.temperature,self.log10_Kp_Pa,self.partition_H2,self.partition_HI,self.log_phi_knots,self.trend):a.setflags(write=False)

    def log_phi(self,temperature):
        T=np.asarray(temperature,dtype=float)
        if not np.isfinite(T).all()or np.any(T<self.domain[0])or np.any(T>self.domain[1]):raise ValueError('Outside qualified H2 table domain; no extrapolation/zero fallback')
        log_kp=self.interpolator(np.log10(T))-self.dissociation_eV*5040/T
        return np.log(KB*T)-np.log(10)*log_kp

    def phi(self,temperature):
        with np.errstate(over='raise',under='raise',invalid='raise'):result=np.exp(self.log_phi(temperature))
        if not np.isfinite(result).all()or np.any(result<=0):raise ValueError('H2 coefficient not representable')
        return result

    def hydrogen_network(self,temperature,total_H_m3,neutral_atomic_fraction,hminus_per_neutral):
        T,N,f0,minus=np.broadcast_arrays(temperature,total_H_m3,neutral_atomic_fraction,hminus_per_neutral)
        if (not all(np.isfinite(x).all()for x in (N,f0,minus))or np.any(N<=0)or np.any(f0<=0)or np.any(f0>1)or np.any(minus<0)):
            raise ValueError('Invalid fixed-ne hydrogen network inputs')
        phi=self.phi(T);a=1/f0+minus
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            neutral=2*N/(a+np.sqrt(a*a+8*phi*N))
            h2=phi*neutral**2;positive=neutral*(1/f0-1);negative=minus*neutral
        fields=dict(HI_m3=neutral,Hplus_m3=positive,Hminus_m3=negative,H2_m3=h2)
        if not all(np.isfinite(v).all()and np.all(v>=0)for v in fields.values()):raise ValueError('Invalid hydrogen solution')
        fields['nuclei_residual']=(neutral+positive+negative+2*h2)/N-1
        return fields
