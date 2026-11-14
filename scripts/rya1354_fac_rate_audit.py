"""Independent pinned FAC formula/table audit; minimal-model totals only."""
import json
import math
import re
import struct
from pathlib import Path
from rya1354_native_parts import load_parts
from rya1354_fac_rate_inventory import summarize_inventory


def constants(text):
    return {k:float(re.search(r'^#define\s+_'+k+r'\s+([0-9.Ee+-]+)',text,re.M)[1])
            for k in ('FINE_STRUCTURE_CONST','RATE_AU','HARTREE_EV')}


def formula(m,e,s,c):
    aw=c['FINE_STRUCTURE_CONST']*e
    x=s*s/(2*abs(m)+1.0);x*=e
    if abs(m)>1:x*=pow(aw,2*abs(m)-2)
    a=x*2.0*pow(aw,2)*c['FINE_STRUCTURE_CONST']
    return x,a/3.0*c['RATE_AU']


def review(raw,table,basis,c,gauge):
    # Exact qualified minimal basis: every odd/even pair differs only by 3s -> 2p.
    peel=basis.split('Peel subshells:')[1].split('CSF(s):')[0].split()
    if peel!=['1s','2s','2p-','2p','3s']:raise ValueError('E3 basis proof no longer applies')
    triples=[x for x in basis.splitlines()[5:] if not x.startswith(' *') and not x.startswith('*')]
    if len(triples)!=45:raise ValueError('Minimal basis CSF inventory differs')
    for i in range(0,len(triples),3):
        occupations={}
        for n,l,q in re.findall(r'(\d+)([spdfg])[ -]\(\s*(\d+)\)',triples[i]):
            key=(int(n),l);occupations[key]=occupations.get(key,0)+int(q)
        ground={(1,'s'):2,(2,'s'):2,(2,'p'):4}
        odd={(1,'s'):2,(2,'s'):2,(2,'p'):3,(3,'s'):1}
        expected=ground if triples[i+2].rstrip().endswith('+') else odd
        if occupations!=expected:raise ValueError('Minimal basis occupation proof differs')
    records={};zeros=[]
    for line in raw.splitlines():
        if not line.strip() or line.startswith('#'):continue
        v=line.split();u,l,ju,jl,n,m,g,mode,rc,nz,nmk=map(int,v[:11])
        e,s,gf,a=map(float,v[11:]);key=(u,l,m)
        if rc:
            # Rank3 fails single-electron triangle: j(3s)=1/2, j(2p)<=3/2.
            if (u,l,m,rc,nz,nmk) not in ((6,0,-3,-1,0,0),(6,3,-3,-1,0,0)):
                raise ValueError('Unexplained native refusal')
            zeros.append(dict(upper=u,lower=l,multipole=m,reason='Minimal-basis one-electron j triangle forbids rank3; native angular_nz0/nmk0 corroborates.'))
            continue
        predicted=formula(m,e,s,c)
        if not all(math.isclose(x,y,rel_tol=2e-15,abs_tol=0) for x,y in zip(predicted,(gf,a))):
            raise ValueError('Raw native formula mismatch')
        records[key]=(ju,jl,e,s,gf,a)
    found=set();headers={};counts={};blocks=[]
    for line in table.splitlines():
        if '=' in line:
            k,val=map(str.strip,line.split('=',1))
            if k=='NELE':headers={};blocks.append(headers)
            if blocks:headers[k]=val
            continue
        v=line.split()
        if len(v)!=8 or not v[0].isdigit():continue
        u,ju,l,jl=map(int,v[:4]);m=int(headers['MULTIP']);key=(u,l,m)
        if headers['NELE']!='8' or int(headers['GAUGE'])!=gauge or headers['MODE']!='0':
            raise ValueError('Printed table identity differs')
        if key in found or key not in records:raise ValueError('Unexpected/duplicate table row')
        found.add(key);counts[m]=counts.get(m,0)+1
        rju,rjl,e,s,gf,a=records[key]
        sf=struct.unpack('f',struct.pack('f',s))[0]
        pg,pa=formula(m,e,sf,c)
        expected=[format(x,'.6E') for x in (e*c['HARTREE_EV'],pg,pa,sf)]
        if (ju,jl)!=(rju,rjl) or v[4:]!=expected:raise ValueError('Printed float32 table mismatch')
    for h in blocks:
        if counts.get(int(h['MULTIP']),0)!=int(h['NTRANS']):raise ValueError('Block count differs')
    omitted=set(records)-found
    # Only M3 is suppressed by native trd <1E-30 (in atomic rate units).
    if omitted!={(6,5,3)} or records[(6,5,3)][5]/c['RATE_AU']>=1e-30:
        raise ValueError('Unexplained writer omission')
    if len(zeros)!=2 or len(records)!=11 or len(found)!=10:raise ValueError('Incomplete channel accounting')
    total=math.fsum(r[5] for r in records.values())
    return dict(printed_rows_verified=len(found),raw_nonzero_channels=len(records),
        explicit_model_angular_zeros=zeros,writer_suppressed_nonzero=[dict(upper=6,lower=5,multipole=3,A_s1=records[(6,5,3)][5])],
        exploratory_model_total_A_s1=total,exploratory_model_lifetime_ns=1e9/total)


def summarize(base):
    inventory=summarize_inventory(base)
    files,_=load_parts(base/'fac_rates_complete_v21');root='runs/fac_model_rates_v21/'
    c=constants(files[root+'source_documents/faclib/consts.h'].decode())
    report=dict(status='MINIMAL_MODEL_RATE_ACCOUNTING_VERIFIED',scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Exploratory one-photon totals in the declared15-level minimal FAC model only; not a physical accuracy or calibrated-grade qualification.',cases={})
    for name in inventory['cases']:
        p=root+name+'/'
        report['cases'][name]=review(files[p+'raw.tsv'].decode(),files[p+'rates.tr'].decode(),
            files[p+'o_export_000.c'].decode(),c,2 if name.endswith('_B') else 1)
    return report

if __name__=='__main__':
    print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
