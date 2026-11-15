"""Sirius-only full-closure identity/seniority/sign negative controls."""
import argparse
import json
from pathlib import Path
import platform
import re
import struct
from rya1354_fac_rates import read_files
from rya1354_fac_full_ls_audit import review
from rya1354_fac_recoupling_audit import records


def main():
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Full angular canaries belong on Sirius')
    ap=argparse.ArgumentParser()
    for n in ('source','converted','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args();s=read_files(a.source/'ground_coherent100');c=read_files(a.converted/'ground_coherent100')
    positive=review(s,c);result={}
    for mode in ('d4_seniority','vector_sign','matrix_sign','explicit_stdout_ID'):
        changed=dict(c)
        if mode=='d4_seniority':
            lines=c['state.lsj.c'].decode().splitlines();changed_one=False
            for i,line in enumerate(lines[:-1]):
                if re.search(r'3d\(\s*4\)',line) and '3P2' in lines[i+1]:
                    lines[i+1]=lines[i+1].replace('3P2','3P9',1);changed_one=True;break
            if not changed_one:raise ValueError('No actual d4 seniority canary')
            changed['state.lsj.c']='\n'.join(lines).encode()
        elif mode=='vector_sign':
            lines=c['state.lsj.j'].decode().splitlines();i=next(i for i,l in enumerate(lines) if 'Ssms=' in l)+2
            lines[i]=re.sub(r'(-?\d+\.\d{8})',lambda m:m[0][1:] if m[0].startswith('-') else '-'+m[0],lines[i],count=1)
            changed['state.lsj.j']='\n'.join(lines).encode()
        elif mode=='matrix_sign':
            rr=records(c['state.lsj.T']);_,_,_,nb=struct.unpack('<4i',rr[1])
            for i in range(3+nb,len(rr)):
                if len(rr[i])==16:
                    jj,l,v=struct.unpack('<iid',rr[i])
                    if abs(v)>.01:rr[i]=struct.pack('<iid',jj,l,-v);break
            changed['state.lsj.T']=b''.join(struct.pack('<i',len(r))+r+struct.pack('<i',len(r)) for r in rr)
        else:
            lines=c['jj2lsj.stdout'].decode().splitlines();waiting=False;changed_one=False
            for i,line in enumerate(lines):
                if 'Weights of major contributors to ASF in LS-coupling:' in line:waiting=True;continue
                if waiting and re.match(r'^\s*\d+\s+\d+\s+[+-]\s+.*\bof\b',line):
                    lines[i]=re.sub(r'^(\s*)\d+',r'\g<1>9999',line,count=1);changed_one=True;break
            if not changed_one:raise ValueError('No actual stdout ID canary')
            changed['jj2lsj.stdout']='\n'.join(lines).encode()
        try:review(s,changed)
        except ValueError as e:result[mode]=str(e)
        else:raise ValueError('Full LS mutation admitted: '+mode)
    a.out.write_text(json.dumps(dict(status='POSITIVE_AND_FOUR_NEGATIVE_CONTROLS_PASS',positive_full_LS_coefficients=positive['full_LS_coefficients_reconstructed'],positive_label_ID_overflows=positive['label_ID_overflows'],results=result),indent=2)+'\n')

if __name__=='__main__':main()
