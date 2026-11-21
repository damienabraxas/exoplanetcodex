"""Sirius-only negative tests of full expanded angular replay; files unchanged."""
import argparse
import json
from pathlib import Path
import platform
import re
import struct
from rya1354_fac_rates import read_files
from rya1354_fac_expanded_ls_audit import review
from rya1354_fac_recoupling_audit import records


def main():
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Angular canaries belong on Sirius')
    ap=argparse.ArgumentParser()
    for n in ('source','converted','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args();s=read_files(a.source/'n3_r1_coherent10');c=read_files(a.converted/'n3_r1_coherent10');results={}
    for mode in ('vector_sign','matrix_sign','original_ID','occupation'):
        changed=dict(c)
        if mode=='vector_sign':
            lines=c['state.lsj.j'].decode().splitlines()
            i=next(i for i,l in enumerate(lines) if 'Ssms=' in l)+2
            lines[i]=re.sub(r'(-?\d+\.\d{8})',lambda m:m[0][1:] if m[0].startswith('-') else '-'+m[0],lines[i],count=1)
            changed['state.lsj.j']='\n'.join(lines).encode()
        elif mode=='matrix_sign':
            rr=records(c['state.lsj.T']);_,_,_,nb=struct.unpack('<4i',rr[1])
            for i in range(3+nb,len(rr)):
                if len(rr[i])==16:
                    jj,l,v=struct.unpack('<iid',rr[i])
                    if abs(v)>.01:rr[i]=struct.pack('<iid',jj,l,-v);break
            changed['state.lsj.T']=b''.join(struct.pack('<i',len(r))+r+struct.pack('<i',len(r)) for r in rr)
        elif mode=='original_ID':changed['state.lsj.lbl']=re.sub(rb'(?m)^\s*2\s+0\s+\+',b'999 0 +',c['state.lsj.lbl'],count=1)
        else:changed['state.lsj.c']=c['state.lsj.c'].replace(b'2p( 4)',b'2p( 3)',1)
        try:review(s,changed,1)
        except ValueError as e:results[mode]=str(e)
        else:raise ValueError('Canary silently admitted: '+mode)
    a.out.write_text(json.dumps(dict(status='ALL_FOUR_MUTATIONS_REJECTED',results=results),indent=2)+'\n')

if __name__=='__main__':main()
