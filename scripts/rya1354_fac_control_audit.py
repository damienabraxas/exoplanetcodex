"""Replay capped/coherent FAC controls and exported mixing records; no LS/rates."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import struct

from rya1354_native_parts import load_parts
from rya1354_fac_controls import diagnostics,input_text as outer_input
from rya1354_fac_coherent import input_text as coherent_input
from rya1354_fac_oxygen_levels import inspect_inventory

ROOT=Path(__file__).resolve().parents[1]
FAILED='runs/fac_numerical_controls_v18/'
COHERENT='runs/fac_coherent_controls_v18/'


def mixing(data):
    """Validate pinned FAC G92MIX layout and report algebraic normalization."""
    records=[];offset=0
    while offset<len(data):
        if offset+8>len(data):
            raise ValueError('Truncated FAC mixing record')
        size=struct.unpack_from('<i',data,offset)[0]
        if size<0 or offset+size+8>len(data) or struct.unpack_from('<i',data,offset+size+4)[0]!=size:
            raise ValueError('FAC mixing record markers differ')
        records.append(data[offset+4:offset+size+4]);offset+=size+8
    if len(records)<2 or records[0]!=b'G92MIX' or len(records[1])!=24:
        raise ValueError('FAC mixing header differs')
    nele,ncsf,norb,nlev,nvalues,nblocks=struct.unpack('<6i',records[1])
    if nele!=8 or nlev!=15 or len(records)!=2+4*nblocks:
        raise ValueError('FAC mixing declared inventory differs')
    blocks=[];indices=[];total_basis=0;total_values=0;max_residual=0.0
    for pos in range(2,len(records),4):
        if len(records[pos])!=20:
            raise ValueError('FAC mixing block header size differs')
        number,nb,nv,jplus,parity=struct.unpack('<5i',records[pos])
        if nb<1 or nv!=nb or parity not in (-1,1) or jplus<1 or number!=len(blocks)+1:
            raise ValueError('FAC minimal-model mixing block differs')
        if (len(records[pos+1])!=4*nv or len(records[pos+2])!=8*(nv+1)
                or len(records[pos+3])!=8*nb*nv):
            raise ValueError('FAC mixing block array dimensions differ')
        ids=struct.unpack('<'+'i'*nv,records[pos+1]);indices.extend(ids)
        energies=struct.unpack('<'+'d'*(nv+1),records[pos+2])
        values=struct.unpack('<'+'d'*(nb*nv),records[pos+3])
        if not all(math.isfinite(x) for x in (*energies,*values)):
            raise ValueError('Nonfinite FAC mixing entry')
        for i in range(nv):
            for j in range(nv):
                dot=sum(values[i*nb+k]*values[j*nb+k] for k in range(nb))
                max_residual=max(max_residual,abs(dot-(1 if i==j else 0)))
        total_basis+=nb;total_values+=nb*nv
        blocks.append(dict(number=number,basis=nb,levels=nv,two_J=jplus-1,parity=parity,indices=list(ids)))
    if set(indices)!=set(range(15)) or len(indices)!=15 or total_basis!=ncsf or total_values!=nvalues:
        raise ValueError('FAC mixing aggregate dimensions/indices differ')
    return dict(electrons=nele,CSFs=ncsf,orbitals=norb,levels=nlev,blocks=blocks,
                max_matrix_orthonormality_residual=max_residual,
                scope='Binary dimensions/finite entries/normalization diagnostic only; no LS recoupling or physical state admission.')


def summarize(base):
    failed,fm=load_parts(base/'fac_controls_failed_v18')
    coherent,cm=load_parts(base/'fac_coherent_complete_v18')
    fr=json.loads(failed[FAILED+'receipt.json']);cr=json.loads(coherent[COHERENT+'receipt.json'])
    if (fr['status']!='FAILED_RETAINED' or 'tight10' not in fr['error']
            or cr['status']!='FAC_COHERENT_CONTROLS_COMPLETED_REVIEW_REQUIRED'
            or fr['pin']!=cr['pin'] or fr['binary_sha256']!=cr['binary_sha256']):
        raise ValueError('FAC numerical campaign dispositions/provenance differ')
    cases={}
    spec=[('cap_canary',3.0,1,1),('explicit_default',3.0,256,1),('tight10',.3,256,1),
          ('coherent10',.3,256,10),('coherent100',.03,256,100)]
    for name,multiplier,limit,scale in spec:
        files=failed if scale==1 else coherent;root=(FAILED if scale==1 else COHERENT)+name+'/'
        expected=outer_input(multiplier,limit) if scale==1 else coherent_input(multiplier,scale)
        if files[root+'input.sf'].decode()!=expected:
            raise ValueError('FAC declared numerical controls differ')
        d=diagnostics(files[root+'sfac.stdout'].decode(),multiplier,limit)
        should_pass=name not in ('cap_canary','tight10')
        if d['printed_rule_verified']!=should_pass:
            raise ValueError('FAC numerical disposition differs from native guards')
        item=dict(diagnostics=d,levels_written=root+'o.lev' in files)
        if should_pass:
            levels=inspect_inventory(files[root+'o.lev'].decode())
            item['levels']=levels
            item['mixing_export']=mixing(files[root+'o_export_000.cm'])
            text=files[root+'o_export_000.c'].decode()
            lines=[l for l in text.splitlines()[5:] if not l.startswith(' *') and not l.startswith('*')]
            if len(lines)!=45 or any(sum(int(x) for x in re.findall(r'\(\s*(\d+)\)',lines[i]))!=8 for i in range(0,45,3)):
                raise ValueError('FAC exported CSF electron/record inventory differs')
        cases[name]=item
    if any(k.startswith(FAILED+'tight100/') for k in failed):
        raise ValueError('Failed outer-only ladder silently continued')
    defaults=failed[FAILED+'source_documents/faclib/consts.h'].decode()
    manual=failed[FAILED+'source_documents/doc/manual.tex'].decode()
    if (re.search(r'#define OPTTOL\s+3\.0',defaults) is None or re.search(r'#define OPTNITER\s+2048',defaults) is None
            or '$10^{-6}$' not in manual):
        raise ValueError('Pinned default discrepancy evidence differs')
    return dict(schema='rya1354.fac-numerical-audit.v18',cases=cases,
        failed_receipt_sha256=hashlib.sha256(failed[FAILED+'receipt.json']).hexdigest(),
        coherent_receipt_sha256=hashlib.sha256(coherent[COHERENT+'receipt.json']).hexdigest(),
        default_discrepancy='Pinned source3.0/2048 versus manual prose1e-6/100. Runtime source governs; explicit settings retained.',
        numerical_result='Tightening outer multiplier alone capped without qualifying; coherently tightened inner/outer settings satisfy pinned printed/native-marker conditions. No physical accuracy estimate.',
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        physical_LS_terms_admitted=False,failed_parts=fm['parts'],coherent_parts=cm['parts'],
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
