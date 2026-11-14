"""Complete expanded raw/native table replay, dynamic writer omissions; metadata."""
from collections import defaultdict
import json
import hashlib
import math
from pathlib import Path
import struct
from rya1354_native_parts import load_parts
from rya1354_fac_expanded_rates import rate_input
from rya1354_fac_rate_audit import constants,formula
from rya1354_fac_controls import diagnostics


def review(raw,table,requests,gauge,mapping,c):
    expected={(r['upper'],r['lower'],r['multipole']) for r in requests};rows={};states=mapping['levels']
    for line in raw.splitlines():
        if not line.strip() or line.startswith('#'):continue
        v=line.split()
        if len(v)!=15:raise ValueError('Expanded raw column count differs')
        u,l,ju,jl,n,m,g,mode,rc,nz,nmk=map(int,v[:11]);key=(u,l,m)
        if key not in expected or key in rows:raise ValueError('Expanded raw duplicate/unexpected channel')
        su,sl=states[str(u)],states[str(l)]
        if (ju,jl,n,g,mode)!=(su['two_J'],sl['two_J'],8,gauge,0) or ju!=2:raise ValueError('Expanded raw identity/gauge/mode differs')
        if rc or nmk:raise ValueError('Unresolved expanded native status/MBPT; never zero-fill')
        e,s,gf,a=map(float,v[11:])
        if not all(math.isfinite(x) for x in (e,s,gf,a)) or e<=0 or min(gf,a)<0:raise ValueError('Invalid expanded raw rate')
        # Binary C/CM stores a common energy offset; allow only reconstruction ulps.
        ep=su['energy_au']-sl['energy_au'];bound=4*math.ulp(max(abs(su['energy_au']),abs(sl['energy_au'])))
        if abs(e-ep)>bound:raise ValueError('Expanded transition energy differs from mapped binary energies')
        if not all(math.isclose(x,y,rel_tol=2e-15,abs_tol=0) for x,y in zip(formula(m,e,s,c),(gf,a))):raise ValueError('Expanded raw formula differs')
        rows[key]=(ju,jl,e,s,gf,a)
    if set(rows)!=expected:raise ValueError('Missing expanded raw channel')
    group=defaultdict(float)
    for (u,l,m),v in rows.items():group[m]+=v[5]/c['RATE_AU']
    required={key for key,v in rows.items() if v[3]!=0 and group[key[2]]>=1e-30}
    found=set();counts=defaultdict(int);blocks=[];header={}
    for line in table.splitlines():
        if '=' in line:
            key,val=map(str.strip,line.split('=',1))
            if key=='NELE':header={};blocks.append(header)
            if blocks:header[key]=val
            continue
        v=line.split()
        if len(v)!=8 or not v[0].isdigit():continue
        u,ju,l,jl=map(int,v[:4]);m=int(header['MULTIP']);key=(u,l,m)
        if (header['NELE'],int(header['GAUGE']),header['MODE'])!=('8',gauge,'0'):raise ValueError('Expanded table header differs')
        if key not in required or key in found:raise ValueError('Expanded table unexpected/duplicate row')
        found.add(key);counts[m]+=1;rju,rjl,e,s,gf,a=rows[key];sf=struct.unpack('f',struct.pack('f',s))[0];pg,pa=formula(m,e,sf,c)
        if (ju,jl)!=(rju,rjl) or v[4:]!=[format(x,'.6E') for x in (e*c['HARTREE_EV'],pg,pa,sf)]:raise ValueError('Expanded printed float32 formula differs')
    if found!=required:raise ValueError('Expanded printed channel missing')
    if len({int(h['MULTIP']) for h in blocks})!=len(blocks):raise ValueError('Duplicate expanded operator block')
    if any(counts[int(h['MULTIP'])]!=int(h['NTRANS']) for h in blocks):raise ValueError('Expanded printed block count differs')
    total=math.fsum(v[5] for v in rows.values())
    return dict(channels=len(rows),printed_rows_verified=len(found),writer_suppressed=[dict(upper=u,lower=l,multipole=m,A_s1=rows[(u,l,m)][5]) for u,l,m in sorted(set(rows)-found)],
        E3_A_s1={str(l):v[5] for (u,l,m),v in rows.items() if m==-3},exploratory_model_total_A_s1=total,exploratory_model_lifetime_ns=1e9/total)


def summarize(base):
    rates,_=load_parts(base/'fac_expanded_rates_complete_v25');angular,_=load_parts(base/'fac_expanded_angular_complete_v24');source,_=load_parts(base/'fac_space_complete_v23')
    mapping=json.loads(angular['evidence/fac_expanded_replay_v24/report.json']);root='runs/fac_expanded_rates_v25/'
    receipt=json.loads(rates[root+'receipt.json'])
    if receipt['status']!='FAC_EXPANDED_RATES_COMPLETED_REVIEW_REQUIRED':raise ValueError('Expanded rates incomplete')
    parent,_=load_parts(base/'fac_rates_complete_v21')
    for key,data in [('source_receipt_sha256',source['runs/fac_space_v23/receipt.json']),('mapping_receipt_sha256',angular['runs/fac_expanded_recoupling_v24/receipt.json']),('rate_parent_receipt_sha256',parent['runs/fac_model_rates_v21/receipt.json'])]:
        if receipt[key]!=hashlib.sha256(data).hexdigest():raise ValueError('Expanded parent receipt provenance differs')
    if receipt['probe']!=json.loads(parent['runs/fac_model_rates_v21/receipt.json'])['probe']:raise ValueError('Expanded diagnostic frontend provenance differs')
    c=constants(parent['runs/fac_model_rates_v21/source_documents/faclib/consts.h'].decode())
    report=dict(schema='rya1354.fac-expanded-rate-audit.v25',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,cases={})
    for name,m in mapping['cases'].items():
        rank=1 if '_r1_' in name else 2;mult,scale=(.03,100) if name.endswith('100') else (.3,10)
        for suffix,gauge in [('B',2),('C',1)]:
            case=name+'_'+suffix;p=root+case+'/';text,requests,lower=rate_input(rank,mult,scale,gauge,m)
            if rates[p+'input.sf'].decode()!=text or not diagnostics(rates[p+'sfac.stdout'].decode(),mult,256)['printed_rule_verified']:raise ValueError('Expanded rate input/numerical guard differs')
            for fn in ('o_export_000.c','o_export_000.cm'):
                if rates[p+fn]!=source['runs/fac_space_v23/'+name+'/'+fn]:raise ValueError('Expanded rate mapped state differs')
            report['cases'][case]=review(rates[p+'raw.tsv'].decode(),rates[p+'rates.tr'].decode(),requests,gauge,m,c)
    return report

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
