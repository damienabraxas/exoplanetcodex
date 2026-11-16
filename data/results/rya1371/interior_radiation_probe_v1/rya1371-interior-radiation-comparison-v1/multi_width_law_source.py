"""Explicit original-MULTI O width component; source owns He treatment, no fitted scale."""
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from lightweaver.broadening import StandardLineBroadener
from pipeline.reference_engines import fingerprint

@dataclass(frozen=True)
class WidthRecord:
    index:int
    i_1based:int
    j_1based:int
    source_gvdw:float
    natural_gamma_s_inv:float
    coefficient_m3_s_inv_K_minus_beta:float
    beta:float
    family:str

@dataclass(frozen=True)
class WidthLaw:
    records:tuple
    temperature_interval:tuple
    source_mass_amu:float
    receipt:dict


def load_width_law(item,atom_file,atom):
    from config.constants import require_codex_path
    root=require_codex_path(item['resource']);path=root/'law.json';report_file=root/'report.json'
    if fingerprint(path)!=item['receipt']or fingerprint(report_file)!=item['report_receipt']:raise ValueError('Changed source-width law packet')
    d=json.loads(path.read_text());report=json.loads(report_file.read_text())
    if report['output']!=item['receipt']or d['schema']!='codex.original_multi_vdw_law/1'or d['purpose']!='experimental-fixed-O-width-comparison'or d['scientific_admission']!='HOLD'or d['source_pin']!='4d6aa3c39d634b38da772e469a5692435eba1354':raise ValueError('Unqualified source width purpose/schema')
    if atom.symbol!='O'or len(atom.levels)!=41 or len(atom.lines)!=256 or d['atom']!=fingerprint(atom_file)or d['source_atomic_mass_amu']!=atom.weight:raise ValueError('Source width atom mismatch')
    expected_cases={'half_neutral_H','plus_10K_frozen_EOS','xi_0p25kms','xi_2kms','xi_4kms'}
    if d['roundoff_bound']!=float(128*np.finfo(np.float64).eps)or {x['case']for x in d['control_checks']}!=expected_cases:raise ValueError('Changed source-law validation contract')
    if d['normalization_sample_index']!=0 or d['baseline_max_relative_error']>d['roundoff_bound']or len(d['control_checks'])!=5 or any(x['max_relative_error']>d['roundoff_bound']for x in d['control_checks']):raise ValueError('Unqualified native source-law reconstruction')
    records=tuple(WidthRecord(**r)for r in d['records'])
    if len(records)!=256:raise ValueError('Incomplete source width records')
    for i,(record,line)in enumerate(zip(records,atom.lines)):
        expected_beta=(1-(line.gVdw-np.floor(line.gVdw)))/2 if line.gVdw>=20 else float(np.float32(.3))
        if (record.index,record.i_1based,record.j_1based,record.source_gvdw,record.natural_gamma_s_inv)!=(i,line.i,line.j,line.gVdw,line.gRad)or line.gStark!=0 or record.beta!=expected_beta or not np.isfinite(record.coefficient_m3_s_inv_K_minus_beta)or record.coefficient_m3_s_inv_K_minus_beta<=0:raise ValueError('Source width transition/recipe mismatch')
    interval=tuple(d['temperature_numeric_interval_K'])
    if len(interval)!=2 or not np.isfinite(interval).all()or not 0<interval[0]<interval[1]:raise ValueError('Invalid source width numeric domain')
    checks={str(path):item['receipt'],str(report_file):item['report_receipt']}
    return WidthLaw(records,interval,atom.weight,checks),checks


class OriginalMultiVdwLaw(StandardLineBroadener):
    """Gamma in s^-1, total neutral H in m^-3; no H level-population lookup."""
    def __init__(self,record,neutral_H_m3,temperature_interval):
        self.record=record;self.neutral_H_m3=np.array(neutral_H_m3,copy=True,dtype=float);self.neutral_H_m3.flags.writeable=False;self.temperature_interval=temperature_interval
        if self.neutral_H_m3.ndim!=1 or not np.isfinite(self.neutral_H_m3).all()or np.any(self.neutral_H_m3<0):raise ValueError('Invalid total-neutral-H width density')
        self.audit=dict(line_index=record.index,calls=0,density_role='HOP total neutral H, not H ground/level populations',prescription='original-MULTI2.3-DAMP',extra_helium_factor=False)
    def setup(self,line):
        if (line.i+1,line.j+1)!=(self.record.i_1based,self.record.j_1based)or line.atom.element.Z!=8:raise ValueError('Source width model-line mismatch')
    def broaden(self,atmos,eqPops):
        T=np.asarray(atmos.temperature)
        if T.shape!=self.neutral_H_m3.shape or not np.isfinite(T).all()or np.any(T<self.temperature_interval[0])or np.any(T>self.temperature_interval[1]):raise ValueError('Outside verified source width numeric T/shape support; refresh source oracle')
        with np.errstate(over='raise',invalid='raise'):gamma=self.record.coefficient_m3_s_inv_K_minus_beta*self.neutral_H_m3*T**self.record.beta
        if not np.isfinite(gamma).all()or np.any(gamma<0):raise ValueError('Invalid original source width')
        self.audit['calls']+=1
        return gamma
