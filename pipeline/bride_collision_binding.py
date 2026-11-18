"""Opt-in fixed-background native-O collision binding; no active-SE admission."""
import json
from pathlib import Path
import numpy as np
from lightweaver.collisional_rates import CollisionalRates
from pipeline.reference_engines import fingerprint
from pipeline.native_oxygen_collisions import collision_matrix

class FixedNativeCollisionBinding(CollisionalRates):
    """Whole-atom collision matrix with explicit immutable physical-state binding."""
    def __init__(self,matrix,state,nstar,nHTot,receipt):
        super().__init__(j=1,i=0)
        self.matrix=np.array(matrix,copy=True);self.star=np.array(nstar,copy=True)
        self.fields={name:np.array(state[key],copy=True)for name,key in [('temperature','temperature_K'),('ne','ne_m3'),('height','height_m'),('vturb','vturb_m_s')]}
        self.fields['nHTot']=np.array(nHTot,copy=True);self.receipt=receipt;self.calls=0
        n=len(self.fields['temperature'])
        if self.matrix.shape!=(41,41,n)or self.star.shape!=(41,n)or any(x.shape!=(n,)for x in self.fields.values()):raise ValueError('Invalid whole-O collision binding dimensions')
        if not np.isfinite(self.matrix).all()or np.any(self.matrix<0)or np.any(self.matrix[np.arange(41),np.arange(41)]!=0)or not np.isfinite(self.star).all()or np.any(self.star<=0):raise ValueError('Invalid whole-O collision binding rates')
        for x in [self.matrix,self.star,*self.fields.values()]:x.flags.writeable=False
    def setup(self,atom):
        if atom.element.Z!=8 or len(atom.levels)!=41:raise ValueError('Native collision binding requires exact O41 model')
        self.atom=atom;self.level_signature=tuple((x.E,x.g,x.stage)for x in atom.levels)
    def compute_rates(self,atmos,eqPops,Cmat):
        for name,value in self.fields.items():
            if not np.array_equal(getattr(atmos,name),value):raise ValueError('Changed collision-bound state '+name)
        oxygen=eqPops.atomicPops[self.atom.element]
        if tuple((x.E,x.g,x.stage)for x in self.atom.levels)!=self.level_signature or len(self.atom.collisions)!=1 or self.atom.collisions[0]is not self:raise ValueError('Changed collision-bound level identity or rate stacking')
        if oxygen.model is not self.atom or not np.array_equal(oxygen.nStar,self.star)or not np.array_equal(oxygen.nTotal,self.star.sum(axis=0)):raise ValueError('Changed collision-bound O identity/nStar/total')
        if Cmat.shape!=self.matrix.shape or np.any(Cmat):raise ValueError('Collision binding requires empty whole-atom matrix; no rate stacking')
        Cmat[:]=self.matrix;self.calls+=1

def load_collision_binding(item,atom_file,source,state,nstar):
    from config.constants import require_codex_path
    root=require_codex_path(item['resource']);report_file=root/'report.json';matrix_file=root/'matrices.npz'
    if fingerprint(report_file)!=item['report_receipt']or fingerprint(matrix_file)!=item['matrix_receipt']or fingerprint(atom_file)!=item['atom']:raise ValueError('Changed source collision packet/atom')
    report=json.loads(report_file.read_text())
    if not report['numerical_contract_pass']or not report['stored_upward_float32_bit_identical']or report['stored_upward_cells']!=45920 or report['counts']!={'CE':741,'CI':39,'CH':302,'CH0':24}or source.symbol!='O'or len(source.levels)!=41 or len(source.lines)!=256:raise ValueError('Unqualified native collision packet')
    checks={str(report_file):item['report_receipt'],str(matrix_file):item['matrix_receipt'],str(atom_file):item['atom']}
    for path,fp in report['inputs'].items():
        if fingerprint(path)!=fp:raise ValueError('Changed native collision provenance '+path)
        checks[path]=fp
    table=Path(next(p for p in report['inputs']if p.endswith('/native/atoms/O.txt')));lines=table.read_text().splitlines();h=next(i for i,x in enumerate(lines)if x.split()[:3]==['key','i','j']);records=[x.split()for x in lines[h+1:]];lh=next(i for i,x in enumerate(lines)if x.split()==['ev','g','ion','label']);levels=[list(map(float,x.split()[:4]))for x in lines[lh+1:lh+42]]
    packet=Path(next(p for p in report['inputs']if p.endswith('/packet.npz')))
    with np.load(packet)as d:
        for key in ['temperature_K','ne_m3','rho_kg_m3']:np.testing.assert_array_equal(d[key],state[key])
        ground=d['H_ground_m3'].copy()*1e-6;neutral=d['H_I_m3'].copy()*1e-6;nH=d['nHTot_m3'].copy()
    with np.load(matrix_file)as d:
        np.testing.assert_array_equal(d['nstar_cm3']*1e6,nstar);q=d['coefficients_cm3'].copy();expected=d['native_contract'].copy()
    matrix=collision_matrix(records,q,state['temperature_K'],state['ne_m3']*1e-6,ground,neutral,nstar*1e-6,np.array([x[1]for x in levels]),np.array([x[2]for x in levels]),report['native_e_k'])
    # Multiplication/division across SI and CGS can round; source output precision
    # is much coarser. Require a fixed float64 arithmetic bound, no fitted scale.
    np.testing.assert_allclose(matrix,expected,rtol=256*np.finfo(float).eps,atol=0)
    return FixedNativeCollisionBinding(matrix,state,nstar,nH,checks),checks
