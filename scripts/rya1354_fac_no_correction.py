"""Fixed-SD potential x diagonal-correction missing cells; Sirius-only."""
import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import platform
import re
import shutil
import subprocess
from rya1354_fac_space import input_text as ground_input,inspect
from rya1354_fac_potential import input_text as equal_input
from rya1354_fac_expanded_rates import rate_input as ground_rate_input
from rya1354_fac_space_audit import transport
from rya1354_fac_controls import diagnostics
from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import read_files,verify_receipt
from rya1354_fac_recouple import CONTROLS,GRASP_PIN
from rya1354_solver_smoke import run
from rya1354_fac_expanded_ls_audit import review as angular_review
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants
from rya1354_fac_config_probe import build


def instrument(text):
    marker="Structure('o.lev.b', "
    if text.count(marker)!=1:raise ValueError('Configuration diagnostic insertion differs')
    return text.replace(marker,"ConfigDiagnostics('configs.before.tsv')\n"+marker,1)+"ConfigDiagnostics('configs.after.tsv')\n"


def input_text(policy,mult,scale):
    text=ground_input(2,mult,scale) if policy=='ground' else equal_input(mult,scale) if policy=='equal' else None
    if text is None:raise ValueError('Undeclared potential policy')
    text=text.replace('ConfigEnergy(0)',"SetConfigEnergyMode(-1)\nSetOption('radial:config_energy',-1)\nConfigEnergy(0)",1)
    return instrument(text)


def rate_input(policy,mult,scale,gauge,mapping):
    text,requests,lower=ground_rate_input(2,mult,scale,gauge,mapping);prefix=ground_input(2,mult,scale)
    if not text.startswith(prefix):raise ValueError('No-correction rate prefix differs')
    return input_text(policy,mult,scale)+text[len(prefix):],requests,lower


def config_state(text,disabled):
    header=re.match(r'# FAC_CONFIG_V1 mode=(-?\d+) groups=(\d+)\n',text)
    if not header:raise ValueError('Configuration state header differs')
    mode,ng=map(int,header.groups());rows=[];keys=set();groups=set()
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):continue
        f=line.split()
        if len(f)!=7:raise ValueError('Configuration state columns differ')
        g,i,nele=map(int,f[:3]);energy,delta,shift=map(float,f[3:6]);key=(g,i)
        if key in keys or nele!=8 or not all(math.isfinite(x) for x in (energy,delta,shift)):raise ValueError('Invalid configuration state')
        keys.add(key);groups.add(g);rows.append(dict(group=g,index=i,name=f[6],energy=energy,delta=delta,shift=shift))
    if not rows or groups!=set(range(ng)):raise ValueError('Incomplete configuration group state')
    for g in groups:
        ids=[r['index'] for r in rows if r['group']==g]
        if set(ids)!=set(range(len(ids))):raise ValueError('Incomplete configuration indices')
    if disabled and (mode!=-1 or any(r['delta']!=0 or r['shift']!=0 for r in rows)):raise ValueError('Configuration correction/shift is not disabled')
    return dict(mode=mode,groups=ng,configurations=len(rows),max_abs_delta=max(abs(r['delta']) for r in rows),max_abs_shift=max(abs(r['shift']) for r in rows))


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','grasp','isodata','baseline','enabled','minimal','probe','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic correction controls belong on Sirius')
    for n,v in vars(a).items():setattr(a,n,v.resolve())
    baseline=verify_receipt(a.baseline);enabled=verify_receipt(a.enabled);minimal=verify_receipt(a.minimal)
    if (baseline['status']!='FAC_EXPANDED_LEVELS_COMPLETED_REVIEW_REQUIRED' or enabled['status']!='FAC_POTENTIAL_CONTROL_COMPLETED_REVIEW_REQUIRED'
        or minimal['status']!='FAC_COHERENT_CONTROLS_COMPLETED_REVIEW_REQUIRED' or baseline['pin']!=PIN
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.grasp,text=True).strip()!=GRASP_PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.grasp,text=True).strip()
        or sha(a.source/'sfac/sfac')!=baseline['binary_sha256']):raise ValueError('No-correction source provenance differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*'):
        if p.is_file():shutil.copyfile(p,drivers/p.name)
    for name in ('doc/manual.tex','sfac/sfac.c','faclib/config.c','faclib/radial.c','faclib/structure.c','faclib/consts.h','faclib/dbase.c','faclib/transition.c'):
        dest=a.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/name,dest)
    receipt=dict(schema='rya1354.fac-no-correction.v27',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),fac_pin=PIN,grasp_pin=GRASP_PIN,
        baseline_receipt_sha256=sha(a.baseline/'receipt.json'),enabled_receipt_sha256=sha(a.enabled/'receipt.json'),minimal_receipt_sha256=sha(a.minimal/'receipt.json'),
        cases={},steps=[],scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    def sfac(binary,out,text):
        (out/'input.sf').write_text(text)
        with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
            p=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=900,check=False)
        if p.returncode:raise ValueError(f'No-correction native exit{p.returncode}: '+str(out))
    try:
        probe=build(a.source,a.probe,a.out/'frontend_evidence');receipt['probe']=probe;binary=Path(probe['path'])
        canary=a.out/'enabled_minimal_canary';canary.mkdir();original=a.minimal/'coherent100'
        sfac(binary,canary,instrument((original/'input.sf').read_text()))
        if not diagnostics((canary/'sfac.stdout').read_text(),.03,256)['printed_rule_verified']:raise ValueError('Diagnostic canary numerical failure')
        for fn in ('o_export_000.c','o_export_000.cm'):
            if (canary/fn).read_bytes()!=(original/fn).read_bytes():raise ValueError('Read-only diagnostic changed enabled minimal state')
        before=config_state((canary/'configs.before.tsv').read_text(),False);after=config_state((canary/'configs.after.tsv').read_text(),False)
        if before['mode']!=0 or before['max_abs_delta']==0 or after['max_abs_delta']==0:raise ValueError('Enabled diagnostic positive control failed')
        receipt['enabled_diagnostic_canary']=dict(before=before,after=after,original_CM_byte_identical=True)
        for policy in ('ground','equal'):
            for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
                name=policy+'_'+label;case=a.out/name;case.mkdir();levels=case/'levels';levels.mkdir();sfac(binary,levels,input_text(policy,mult,scale))
                d=diagnostics((levels/'sfac.stdout').read_text(),mult,256)
                if not d['printed_rule_verified']:raise ValueError('No-correction stopping failed')
                inspect((levels/'o.lev').read_text(),2);s=read_files(levels)
                states={n:config_state(s['configs.'+n+'.tsv'].decode(),True) for n in ('before','after')}
                if any(r['groups']!=15 for r in states.values()):raise ValueError('No-correction configuration group inventory differs')
                tr=transport(s['o.basis'].decode(),s['o_export_000.c'].decode(),s['o_export_000.cm'],2)
                angular=case/'angular';angular.mkdir()
                for old,new in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:shutil.copyfile(levels/old,angular/new)
                shutil.copyfile(a.isodata,angular/'isodata');step=run(a.grasp/'build/bin/jj2lsj',CONTROLS,angular,900);receipt['steps'].append(step)
                if step['semantic_errors']:raise ValueError('No-correction angular semantic refusal')
                c=read_files(angular)
                if c['state.c']!=s['o_export_000.c'] or c['state.cm']!=s['o_export_000.cm']:raise ValueError('No-correction angular input changed')
                mapping=angular_review(s,c,2);(case/'mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
                item=dict(diagnostics=d,transport=tr,configuration_state=states,targets=mapping['targets'],rates={});receipt['cases'][name]=item
                for suffix,gauge in [('B',2),('C',1)]:
                    out=case/suffix;out.mkdir();text,requests,lower=rate_input(policy,mult,scale,gauge,mapping);sfac(binary,out,text)
                    rd=diagnostics((out/'sfac.stdout').read_text(),mult,256)
                    if not rd['printed_rule_verified']:raise ValueError('No-correction rate stopping failure')
                    configs={n:config_state((out/('configs.'+n+'.tsv')).read_text(),True) for n in ('before','after')}
                    for fn in ('o_export_000.c','o_export_000.cm'):
                        if (out/fn).read_bytes()!=s[fn]:raise ValueError('No-correction regenerated rate state differs')
                    reviewed=rate_review((out/'raw.tsv').read_text(),(out/'rates.tr').read_text(),requests,gauge,mapping,constants((a.source/'faclib/consts.h').read_text()))
                    item['rates'][suffix]=dict(gauge=gauge,requests=requests,lower_states=lower,configuration_state=configs,diagnostics=rd,audit=reviewed)
                (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_NO_CORRECTION_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
