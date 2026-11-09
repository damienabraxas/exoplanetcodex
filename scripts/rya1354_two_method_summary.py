"""Hash-checked n3/n4/n5 comparison; shared omitted physics is not uncertainty."""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile

ROOT=Path(__file__).resolve().parents[1]
NUM=r'[\d.D+\-]+'


def scf_diagnostics(stdout,stderr,*,function_accy=.0625**6):
    """Report actual native termination/residuals; energy stop is not all-function convergence."""
    iterations=re.findall(r'Iteration number\s+(\d+)',stdout)
    if not iterations:
        raise ValueError('Native SCF iteration inventory absent')
    final=stdout.rsplit('Iteration number',1)[1]
    residuals={}
    for line in final.splitlines():
        fields=line.split()
        if len(fields)>=11 and re.fullmatch(r'\d+[spdfgh]-?',fields[0]) and fields[2] in ('1','2','3','4'):
            residuals[fields[0]]=float(fields[4].replace('D','E'))
    if not residuals:
        raise ValueError('Native final orbital residual table absent')
    if not math.isfinite(function_accy) or function_accy<=0:
        raise ValueError('Invalid explicitly declared native function criterion')
    accy=function_accy
    declared=re.findall(r'(?:Revise the default|Revised) ACCY\s*=\s*([\d.EeDd+\-]+)',stderr)
    if declared:
        native=float(declared[-1].replace('D','E').replace('d','e'))
        if not math.isclose(native,accy,rel_tol=1e-12,abs_tol=0):
            raise ValueError('Declared function criterion differs from native ACCY prompt')
    limited='MAXIMUM ITERATIONS IN SCF EXCEEDED' in stderr.upper()
    return dict(iterations=int(iterations[-1]),iteration_limit_hit=limited,
                reported_orbital_self_consistency=residuals,max_reported_self_consistency=max(residuals.values()),
                native_function_ACC=accy,all_reported_functions_below_ACC=(max(residuals.values())<=accy),
                stopping_rule='Pinned scf.f90 accepts either SCNSTY<=ACCY or relative EOL-energy change <0.001*ACCY. Native termination alone does not prove all orbital residuals meet ACCY.')


def capsule_files(path,*,verify=True):
    with tarfile.open(path) as archive:
        # Standard tar hardlinks deduplicate byte-identical native artifacts.
        # extractfile resolves only archive members; no filesystem extraction.
        files={m.name:archive.extractfile(m).read() for m in archive.getmembers() if m.isfile() or m.islnk()}
    if verify:
        verify_artifacts(files)
    return files


def verify_artifacts(files):
    """Validate receipts after explicitly joining all parts of a split capsule."""
    for name,content in files.items():
        if not name.endswith('/receipt.json'):
            continue
        receipt=json.loads(content)
        parent=name.rsplit('/',1)[0]
        for relative,record in receipt['files'].items():
            key=parent+'/'+relative
            if key not in files and Path(relative).name.startswith('mcp.'):
                continue
            if key not in files or hashlib.sha256(files[key]).hexdigest()!=record['sha256']:
                raise ValueError(f'Native artifact missing or changed: {key}')


def occupations(text):
    return Counter(tuple(sorted((nl,int(q)) for nl,q in re.findall(r'([1-9][spdfg])\(\s*(\d+)\)',line)))
                   for line in text.splitlines() if '(' in line)


def parse_mchf(text):
    number=lambda v:float(v.replace('D','E'))
    wave=re.search(r'([\d.]+) CM-1\s+([\d.]+) ANGS\(VAC\)\s+([\d.]+) ANGS\(AIR\)',text)
    length=re.search(rf'E1\s+length:\s+S =\s+({NUM})\s+GF =\s+({NUM})\s+AKI =\s+({NUM})',text)
    velocity=re.search(rf'velocity:\s+S =\s+({NUM})\s+GF =\s+({NUM})\s+AKI =\s+({NUM})',text)
    if not all((wave,length,velocity)) or '2p(3)4S3.3s_5S' not in text or '2p(3)4S3.3p_5P' not in text:
        raise ValueError('Incomplete or wrong-identity native MCHF multiplet')
    en,vac,air=map(number,wave.groups())
    sl,gl,al=map(number,length.groups())
    sv,gv,av=map(number,velocity.groups())
    if not all(math.isfinite(v) and v>0 for v in (en,vac,air,sl,gl,al,sv,gv,av)):
        raise ValueError('Invalid native value')
    return dict(native_energy_cm1=en,native_vacuum_A=vac,printed_air_A_rejected=air,
                air_warning='Native probab.f uses 206030 in its dispersion numerator; disagrees with project wavelength SSOT. Air column is not used.',
                S_length_au=sl,S_velocity_au=sv,gf_length=gl,gf_velocity=gv,A_length_s1=al,A_velocity_s1=av,
                gauge_fraction_A=abs(al-av)/max(al,av),comparison_scope='LS_MULTIPLET_TOTAL_ONLY')


def summarize(destination=None):
    base=ROOT/'data/audit/rya1354'
    gpath=base/'n5_and_mchf_control_v1.tar.gz'
    mpath=base/'mchf_ladder_v3_all_attempts.tar.gz'
    grasp=capsule_files(gpath)
    mchf=capsule_files(mpath)
    original=capsule_files(base/'continuation_20261108_compact.tar.gz')
    # Regeneration in canonical orbital order changes recoupling presentation,
    # not the physical occupation inventory or number of CSFs per occupation.
    basis_checks=[]
    for n in (3,4):
        for tag in ('lower','upper'):
            old=mchf[f'runs/mchf_ladder_v1/n{n}_{tag}/cfg.inp'].decode()
            new=mchf[f'runs/mchf_ladder_v2/n{n}_{tag}/cfg.inp'].decode()
            if occupations(old)!=occupations(new):
                raise ValueError('Regenerated physical CSF occupation inventory changed')
            basis_checks.append(dict(n=n,state=tag,same_occupation_multiset=True,csfs=sum(occupations(new).values())))
    control=json.loads(grasp['validation/mchf_native_control_v1/receipt.json'])
    if not control['byte_identical_to_reference']:
        raise ValueError('Modern MCHF upstream control did not reproduce native reference bytes')
    old=json.loads((base/'nonblind_v1/native_summary.json').read_text())['ladder']
    g5=grasp['runs/valence_sd_v6/n5_E1/lower.upper.t'].decode()
    pattern=(rf'f2\s+\d+\s+([123])\s+\+\s+f1\s+1\s+2\s+-\s+([\d.]+)\s+C\s+({NUM})\s+({NUM})\s+({NUM})'
             rf'\s+B\s+({NUM})\s+({NUM})\s+({NUM})')
    g5rows=[]
    for j,en,ac,gc,sc,ab,gb,sb in re.findall(pattern,g5):
        conv=lambda v:float(v.replace('D','E'))
        ab,ac,gb,gc=map(conv,(ab,ac,gb,gc))
        g5rows.append(dict(upper_J=int(j),native_energy_cm1=float(en),A_B_s1=ab,A_C_s1=ac,
                           gf_B=gb,gf_C=gc,gauge_fraction_A=abs(ab-ac)/max(ab,ac)))
    if {r['upper_J'] for r in g5rows}!={1,2,3} or len(g5rows)!=3:
        raise ValueError('Missing or duplicate GRASP n5 component')
    old['n5']=g5rows
    rows=[]
    for n in (3,4,5):
        native=mchf[f'runs/mchf_ladder_v3/n{n}_E1/lower.upper.ls'].decode()
        m=parse_mchf(native)
        gs=sum(r['gf_B'] for r in old[f'n{n}'])
        gfiles=grasp if n==5 else original
        groot='runs/valence_sd_v6' if n==5 else 'runs/valence_sd_v4'
        diagnostics={tag:scf_diagnostics(gfiles[f'{groot}/n{n}_{tag}/rmcdhf.stdout'].decode(),
                                        gfiles[f'{groot}/n{n}_{tag}/rmcdhf.stderr'].decode()) for tag in ('lower','upper')}
        rows.append(dict(n=n,grasp_components=old[f'n{n}'],grasp_summed_gf_B=gs,mchf=m,
                         grasp_SCF_diagnostics=diagnostics,
                         relative_method_difference_gf=(gs-m['gf_length'])/m['gf_length']))
    result=dict(schema='rya1354.two-method-ladder.v2',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,
                capsules={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (gpath,mpath,base/'continuation_20261108_compact.tar.gz')},
                input_summary_sha256=hashlib.sha256((base/'nonblind_v1/native_summary.json').read_bytes()).hexdigest(),
                native_control=dict(energy_difference_hartree=control['energy_difference_hartree'],byte_identical=True),
                basis_checks=basis_checks,rows=rows,
                rejected_attempt=dict(run='mchf_ladder_v1',reason='Biotr reports inconsistent shared orbital order, exits zero and writes no .ls; raw receipt completion status is invalid for transition data'),
                warning='Close method agreement shares omitted correlation and lineage; it is not an experimental error estimate or controlled isolation of relativity. Independent radial/native Rydberg conventions remain explicit.',
                uncertainty=None,covariance=None,adoption_allowed=False)
    (destination or base/'source_curation_v2/two_method_ladder.json').write_text(json.dumps(result,indent=2)+'\n')
    print([(r['n'],r['grasp_summed_gf_B'],r['mchf']['gf_length'],r['relative_method_difference_gf']) for r in rows])


if __name__=='__main__':
    summarize()
