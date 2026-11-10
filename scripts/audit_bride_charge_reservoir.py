"""Inspect LTE collider reservoirs without solving RT/SE or assigning rate provenance."""
import sys,json,argparse
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from config.constants import require_codex_path
from pipeline.bride_reference import read_marcs_structure,make_reference_context,fingerprint
from pipeline.multi_atom import read_multi_atom
from lightweaver.atomic_table import PeriodicTable
import lightweaver.constants as c
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
refpath=args.reference
r=json.loads(refpath.read_text()); ap=require_codex_path('experiments.bride_oxygen_atom'); mp=require_codex_path('grids.marcs_standard')/r['atmosphere_file']['name']
if fingerprint(ap)!=r['atom'] or fingerprint(mp)!=r['atmosphere_file']:
 raise ValueError('Changed reference atom or atmosphere')
s=read_marcs_structure(mp);a=read_multi_atom(ap);ctx=make_reference_context(s,a,rays=3,qscale=1)
eq=ctx.kwargs['eqPops']; h=eq.atomicPops[PeriodicTable['H']];o=eq.atomicPops[PeriodicTable['O']]
T=s['temperature_K'];ne=ctx.atmos.pyAtmos.ne
CI=c.HPlanck**2/(2*np.pi*c.MElectron*c.KBoltzmann)
hm_ground=ne*h.nStar[0]*.25*(CI/T)**1.5*np.exp(c.E_ION_HMIN/(c.KBoltzmann*T))
rows=[]
for target in [-4.9177,-3.,-1.2621,0.]:
 d=int(np.argmin(abs(s['logtau_500']-target)));rows.append(dict(depth=d,logtau=float(s['logtau_500'][d]),T_K=float(T[d]),H_ground_m3=float(h.nStar[0,d]),Hplus_m3=float(h.nStar[-1,d]),Hminus_LW_m3=float(eq.HminPops[d]),Hminus_ground_Saha_m3=float(hm_ground[d]),LW_Hminus_over_ground_Saha=float(eq.HminPops[d]/hm_ground[d]),O_lower_over_ion_LTE=float(o.nStar[5,d]/o.nStar[39,d]),qreverse_over_qforward_if_Hminus=float(h.nStar[0,d]/hm_ground[d]*o.nStar[5,d]/o.nStar[39,d])))
ch0=[x for x in a.collisions if x.kind=='CH0']
out=dict(schema='codex.bride_charge_reservoir_probe/1',scientific_admission='HOLD',source_atom=fingerprint(ap),source_atmosphere=fingerprint(mp),reference=fingerprint(refpath),probe=fingerprint(Path(__file__)),CH0_count=len(ch0),raw_first_indices=sorted(set(x.j for x in ch0)),raw_second_indices=sorted(set(x.i for x in ch0)),rows=rows,interpretation='No rate/source provenance assignment. Hminus candidate uses ground-state Saha reservoir, separately from existing LW Hminus implementation. No RT/SE solve, coefficient change or model admission.')
p=args.output;p.mkdir(parents=True,exist_ok=False);(p/'report.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
