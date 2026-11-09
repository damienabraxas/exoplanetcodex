"""Read-only CN source/partition/isotope contract; no chemistry or RT execution."""
import argparse,csv,json,shlex,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def cn_partition_record(text):
    records=[]
    for number,line in enumerate(text.splitlines(),1):
        if not line.strip() or line.lstrip().startswith('*'):continue
        fields=shlex.split(line)
        if len(fields)<5 or fields[2].strip()!='C N':continue
        nq=int(fields[4]);coeff=np.asarray([float(x.replace('D','E'))for x in fields[5:5+nq]])
        tail=list(map(int,fields[5+nq:]))
        if nq<1 or nq>10 or len(coeff)!=nq or tail!=[2,6,1,7,1]or not np.isfinite(coeff).all():raise ValueError('Malformed neutral CN partition entry')
        records.append(dict(line=number,raw=line,code=fields[1].strip(),file_dissociation=float(fields[3]),coefficients=coeff.tolist()))
    if len(records)!=1:raise ValueError('Expected exactly one neutral CN partition record')
    return records[0]


def project_partition(temperature,coefficients):
    t=np.asarray(temperature,dtype=np.float32);c=np.asarray(coefficients,dtype=float)
    if not np.isfinite(t).all()or np.any(t<=0)or c.ndim!=1 or not len(c)or not np.isfinite(c).all():raise ValueError('Invalid partition projection input')
    # Native TEMP is default REAL: LOG(TEMP) is evaluated before assignment to
    # DOUBLE LOGT; Horner coefficients/accumulator are double. Not a libm parity claim.
    x=np.log(t).astype(float);q=np.full(t.shape,c[-1])
    for a in c[-2::-1]:q=q*x+a
    q=np.exp(q)
    if not np.isfinite(q).all()or np.any(q<=0):raise ValueError('Invalid projected partition')
    return q.astype(np.float32)


def isotope_weights(text):
    rows={}
    for line in text.splitlines():
        f=line.split()
        if len(f)!=4:raise ValueError('Malformed isotope source row')
        z=int(float(f[0]));mass=int(f[1]);value=float(f[3])
        if (z,mass)in rows or not np.isfinite(value)or not 0<=value<=1:raise ValueError('Invalid isotope source')
        rows[z,mass]=value
    pairs=[((6,12),(7,14)),((6,13),(7,14)),((6,12),(7,15)),((6,13),(7,15))];out=[]
    for carbon,nitrogen in pairs:
        raw=rows[carbon]*rows[nitrogen];rounded=float(f'{rows[carbon]:.3f}')*float(f'{rows[nitrogen]:.3f}')
        out.append(dict(isotopologue=f'{carbon[1]}C{nitrogen[1]}N',source_product=raw,serialized_product=rounded,relative_change=rounded/raw-1 if raw else None))
    return rows,out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure
    plan_path=ROOT/'config/experiments/rya1281_cn_interface.json';plan=json.loads(plan_path.read_text());baseline=ROOT/plan['baseline'];run=json.loads(baseline.read_text())
    engine=require_codex_path('engines.ispec');source=engine/'synthesizer/turbospectrum/source';data=engine/'synthesizer/turbospectrum/DATA';pin=plan['source_commit'];checks=[]
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=engine,text=True).strip()!=pin:raise ValueError('Wrong source lineage')
    def protected(path,expected=None):
        relative=path.relative_to(engine)
        if expected is None:
            if path.read_bytes()!=subprocess.check_output(['git','show',pin+':'+str(relative)],cwd=engine):raise ValueError('Changed native source '+str(relative))
        elif fingerprint(path)!=expected:raise ValueError('Changed separately held input '+str(relative))
        receipt=fingerprint(path);checks.append(dict(path=str(relative),receipt=receipt,authority='pinned native git'if expected is None else 'accepted original input receipt'));return receipt
    snippets=[]
    owners={'bsyn.f':['ntt(k)=partryck(k,molindex)/boltz/t(k)/ro(k)*symmfactor','call partffordepth(ntau,t,lele,fpartition)','ntt(k)=ntt(k)*isotopfrac(atom(i),isotope(i))','f=10**(gfelog)'],
      'partffordepth.f':['fpartition(k)=sngl(qmol(molindex))'],
      'molecpartf.f':['logt=log(temp)','qq=qq*logt + qcoeff(n,j)','qmol(j)=exp(qq)'],
      'depth.f':['N(I) = G*NTT(i)*EX10(-CHIE*THETA(i))/Q(I)','STIM(I) = 1. - EXP(-CONSTA/T(i))'],
      'die_pe_lu.f':['qprod=1./qmol(j)']}
    for filename,needles in owners.items():
        path=source/filename;protected(path);lines=path.read_text().splitlines()
        for needle in needles:
            hits=[dict(line=i,text=line)for i,line in enumerate(lines,1)if needle in line]
            if len(hits)!=1:raise ValueError('Ambiguous source contract '+needle)
            snippets.append(dict(file=filename,**hits[0]))
    partition_path=data/'IRWIN_molecules_v15.1.dat';protected(partition_path);protected(data/'LISTE_molecules_all_v19.1.dat')
    record=cn_partition_record(partition_path.read_text());model=require_codex_path('grids.marcs_standard')/run['model']['name']
    if fingerprint(model)!=run['model']:raise ValueError('Changed reference model')
    structure=read_marcs_structure(model);temperature=structure['temperature_K'];Q=project_partition(temperature,record['coefficients'])
    isotope_path=engine/'input/isotopes/SPECTRUM.lst';protected(isotope_path,run['isotopes'])
    if fingerprint(isotope_path)!=run['isotopes']:raise ValueError('Changed original isotope source')
    rows,weights=isotope_weights(isotope_path.read_text())
    old=require_codex_path('experiments.bride_native_molecules')/'CN/input.txt';case=next(c for c in run['cases']if c['name']=='CN')
    if fingerprint(old)!=case['input']:raise ValueError('Changed native CN input')
    for z,m in [(6,12),(6,13),(7,14),(7,15)]:
        expected=f'{z}.{m:03d} {rows[z,m]:.3f}'
        if old.read_text().splitlines().count(expected)!=1:raise ValueError('Native input isotope serialization mismatch')
    inventory=json.loads((ROOT/'data/results/rya1281/missing_opacity_v1/report.json').read_text());families=[]
    for name in plan['families']:
        path=engine/'input/linelists/turbospectrum/molecules'/name;protected(path,next(item['installed']for item in inventory['files']if item['filename']==name));header=path.read_text().splitlines()[:2]
        families.append(dict(filename=name,header=header,declared_count=int(shlex.split(header[0])[-1])))
    for item in checks:
        if fingerprint(engine/item['path'])!=item['receipt']:raise ValueError('Source mutated during audit')
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'partition_projection.csv').open('w',newline='')as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['depth_index','temperature_K','source_polynomial_Q_default_real'])
        w.writerows((i,float(t),float(q))for i,(t,q)in enumerate(zip(temperature,Q)))
    report=dict(schema='codex.bride_cn_interface/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),baseline=fingerprint(baseline),source_pin=pin,scope=plan['scope'],protected_sources=checks,source_equations=snippets,model=run['model'],ndep=len(Q),temperature_min_K=float(min(temperature)),temperature_max_K=float(max(temperature)),partition_min=float(min(Q)),partition_max=float(max(Q)),partition_entry=record,isotope_weights=weights,serialized_cn_input=case['input'],line_lists=families,
      limitations=['Polynomial projection is not native compiled-runtime or independent physical validation. No temperature validity bounds are established by the held coefficient row.','Q is shared for elemental CN; isotope masses/fractions enter separately in inspected native source. No independent nuclear-spin/statistical convention is qualified.','Full-network CN partial pressure and native EOS state are not exported by this audit; no chemistry parity or numerical CN opacity coefficient comparison.','Three-decimal isotope products describe original diagnostic serialization only; preserve the original result and do not silently repair it. Missing 13C15N line file is not a completeness or physical error bound.','Primary release/frame/probability/partition/isotope provenance remains blocked; no Bride insertion or abundance claim.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['ndep','temperature_min_K','temperature_max_K','partition_min','partition_max','isotope_weights']}));print(f'Read {len(checks)} protected sources and unchanged model/native input; wrote report.json and partition_projection.csv to {a.output}')

if __name__=='__main__':main()
