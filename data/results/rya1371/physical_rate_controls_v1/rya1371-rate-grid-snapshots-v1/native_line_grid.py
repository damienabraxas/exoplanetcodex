"""Opt-in source-derived frequency offsets; preserves atomic photon identities."""
import json
import numpy as np
from lightweaver.atomic_model import LineQuadrature
from pipeline.reference_engines import fingerprint


class NativeFractionGrid(LineQuadrature):
    def __init__(self,offsets):
        offsets=np.array(offsets,dtype=float,copy=True)
        if offsets.ndim!=1 or len(offsets)<5 or len(offsets)%2!=1 or not np.isfinite(offsets).all()or np.any(np.diff(offsets)<=0)or offsets[len(offsets)//2]!=0 or np.max(abs(offsets))>=1:raise ValueError('Invalid native fractional frequency offsets')
        offsets.flags.writeable=False;self.offsets=offsets;self.Nlambda=len(offsets)
    def wavelength(self,line,vMicroChar=None):return line.lambda0*(1+self.offsets)
    def __repr__(self):return f'NativeFractionGrid(Nlambda={self.Nlambda})'
    def __hash__(self):return hash(self.offsets.tobytes())


class NativeLineGridControl:
    def __init__(self,packet,receipt):self.packet=packet;self.receipt=receipt
    def apply(self,model):
        if len(model.lines)!=256 or len(self.packet['grids'])!=256:raise ValueError('Require complete native O41 grid control')
        changes=[]
        for index,(line,row)in enumerate(zip(model.lines,self.packet['grids'])):
            if (index,line.i,line.j,line.f)!=(row['index'],row['i'],row['j'],row['f']):raise ValueError('Source line identity mismatch in grid control')
            old=line.wavelength();centre=line.lambda0;strength=line.f
            line.quadrature=NativeFractionGrid(row['offsets_fraction']);line.quadrature.setup(line);new=line.wavelength()
            if line.lambda0!=centre or line.f!=strength:raise ValueError('Grid control changed photon centre/strength')
            changes.append(dict(index=index,i=line.i,j=line.j,original_endpoint_vac_nm=[float(old[0]),float(old[-1])],control_endpoint_vac_nm=[float(new[0]),float(new[-1])],control_points=len(new)))
        return dict(receipt=self.receipt,changes=changes,scope=self.packet['scope'])


def load_line_grid_control(spec,atom_file,*,resolver):
    folder=resolver(spec['resource']);checks={}
    for name in ['grids.json','report.json']:
        path=folder/name
        if fingerprint(path)!=spec[name]:raise ValueError('Changed native grid-control receipt')
        checks[str(path)]=fingerprint(path)
    packet=json.loads((folder/'grids.json').read_text());report=json.loads((folder/'report.json').read_text())
    if packet['schema']!='codex.native_line_grid_control/1'or packet['scientific_admission']!='HOLD'or packet['atom']!=fingerprint(atom_file)or not report['computational_completed']or report['packet']!=fingerprint(folder/'grids.json')or packet['nq_fac']!=2 or packet['qnorm_km_s']!=10 or packet['legacy']is not False:raise ValueError('Unqualified native grid control')
    from pipeline.reference_engines import load_plan_engine
    _,_,_,engine=load_plan_engine(dict(engine='openmulti3d',profile_variant='accurate-kernel-control'))
    if packet['engine']!=engine:raise ValueError('Changed native grid engine identity')
    for path,fp in report['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed native grid oracle input')
        checks[path]=fp
    return NativeLineGridControl(packet,checks),checks
