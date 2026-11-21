"""Transport and raw-query qualification only; unresolved native refusals block totals."""
import json
import math
from pathlib import Path
from rya1354_native_parts import load_parts
from rya1354_fac_recoupling_audit import summarize
from rya1354_fac_oxygen_levels import inspect_inventory
from rya1354_fac_controls import diagnostics
from rya1354_fac_rates import rate_input


def inspect_raw(text, requests, gauge):
    rows=[line.split() for line in text.splitlines() if line.strip() and not line.startswith('#')]
    expected={(r['upper'],r['lower'],r['multipole']) for r in requests}
    seen=set();unresolved=[];finite=0
    for row in rows:
        if len(row)!=15:raise ValueError('Raw column count differs')
        u,l,ju,jl,nele,m,g,mode,status,nz,nmk=map(int,row[:11])
        key=(u,l,m)
        if key not in expected or key in seen:raise ValueError('Raw channel missing/duplicate/unexpected')
        seen.add(key)
        if (ju,nele,g,mode)!=(2,8,gauge,0):raise ValueError('Raw identity/gauge/mode differs')
        energy,amplitude,gf,rate=map(float,row[11:])
        if not math.isfinite(energy) or energy<=0:raise ValueError('Invalid raw energy')
        if status:
            if not all(math.isnan(x) for x in (amplitude,gf,rate)):
                raise ValueError('Native refusal must retain NaN, never zero')
            unresolved.append(dict(upper=u,lower=l,multipole=m,status=status,angular_nz=nz,nmk=nmk))
        else:
            if not all(math.isfinite(x) for x in (amplitude,gf,rate)) or min(gf,rate)<0:
                raise ValueError('Invalid native rate')
            finite+=1
    if seen!=expected:raise ValueError('Raw channel missing/duplicate/unexpected')
    return dict(queries=len(rows),finite_native_rates=finite,unresolved=unresolved,
                total_rate=None,lifetime_ns=None)


def summarize_inventory(base):
    files,_=load_parts(base/'fac_rates_complete_v21')
    control,_=load_parts(base/'fac_coherent_complete_v18')
    mapped=summarize(base);root='runs/fac_model_rates_v21/'
    receipt=json.loads(files[root+'receipt.json'])
    if receipt['status']!='FAC_MODEL_RATES_COMPLETED_REVIEW_REQUIRED':raise ValueError('Run incomplete')
    report=dict(status='RAW_INVENTORY_VERIFIED_RATE_AUDIT_PENDING',scientific_admission='HOLD',
                exposure='EXPOSED_NONBLIND',accepted_lifetime=None,uncertainty=None,covariance=None,
                calibrated_sigma=None,cases={})
    for case,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
        source='runs/fac_coherent_controls_v18/'+case+'/'
        levels=inspect_inventory(control[source+'o.lev'].decode())
        for suffix,gauge in [('B',2),('C',1)]:
            name=case+'_'+suffix;prefix=root+name+'/'
            text,requests=rate_input(mult,scale,gauge,mapped['cases'][case]['recoupling'],levels)
            if files[prefix+'input.sf'].decode()!=text:raise ValueError('Rate input differs')
            diagnostics(files[prefix+'sfac.stdout'].decode(),mult,256)
            for fn in ('o_export_000.c','o_export_000.cm'):
                if files[prefix+fn]!=control[source+fn]:raise ValueError('Mapped state differs')
            report['cases'][name]=inspect_raw(files[prefix+'raw.tsv'].decode(),requests,gauge)
    return report

if __name__=='__main__':
    print(json.dumps(summarize_inventory(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
