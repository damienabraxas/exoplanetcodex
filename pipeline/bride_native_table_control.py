"""Fixed-population optical transfer diagnostics for a consumed native atom table."""
from dataclasses import fields
import copy
import numpy as np
from lightweaver.atomic_model import VoigtLine
from pipeline.wavelength_util import air_to_vac

WINDOWS=((6154.,6162.),(6298.,6303.),(6361.,6366.),(7770.,7777.),(8443.,8449.))

class NativeTableLine(VoigtLine):
    """Explicit transfer wavelength override; forbidden as a coupled-SE atom."""
    @property
    def lambda0_m(self):
        return self.native_center_m


def apply_table_control(model,path,mode):
    if mode not in ['centres','strengths','centres_and_strengths']:raise ValueError('Unknown native table diagnostic')
    text=path.read_text().splitlines();starts=[i for i,x in enumerate(text)if x.split()[:2]==['i','j']and'lam_air'in x]
    if len(starts)!=1:raise ValueError('Ambiguous native line table')
    start=starts[0]+1;changes=[];lines=[]
    for index,line in enumerate(model.lines):
        row=text[start+index].split()
        if len(row)!=13 or tuple(map(int,row[:3]))!=(index,line.i,line.j):raise ValueError('Native line table identity mismatch')
        centre=float(row[3]);loggf=float(row[5])
        if not np.isfinite([centre,loggf]).all():raise ValueError('Nonfinite native line input')
        if not any(lo<=centre<=hi for lo,hi in WINDOWS):lines.append(line);continue
        old_lambda=line.lambda0;old_f=line.f;new=copy.deepcopy(line)
        if mode in ['centres','centres_and_strengths']:
            new=NativeTableLine(**{f.name:getattr(new,f.name)for f in fields(new)if f.init})
            new.native_center_m=float(air_to_vac(centre))*1e-10
        if mode in ['strengths','centres_and_strengths']:new.f=10**loggf/model.levels[line.i].g
        new.setup(model);lines.append(new)
        changes.append(dict(line_index=index,i=line.i,j=line.j,source_vac_nm=old_lambda,transfer_vac_nm=new.lambda0,source_f=old_f,transfer_f=new.f,native_table_air_A=centre,native_loggf=loggf))
    if len(changes)!=10:raise ValueError('Changed optical control support; expected ten source transitions')
    model.lines=lines
    return dict(mode=mode,changes=changes,scope='Ten transitions in predeclared six-window optical support only; canonical AIR-to-vacuum conversion once. All other transitions, atom levels, continua and collisions unchanged. Transfer override is not a self-consistent energy-level/SE atom; fixed populations only. No empirical shift/strength fit or physical/default adoption.')
