"""Fixed-SD-space explicit equal-reference potential control; Sirius-only."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
from rya1354_fac_space import input_text as ground_input,inspect
from rya1354_fac_space_audit import transport
from rya1354_fac_controls import diagnostics
from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import read_files,verify_receipt
from rya1354_fac_recouple import CONTROLS,GRASP_PIN
from rya1354_solver_smoke import run
from rya1354_fac_expanded_ls_audit import review as angular_review
from rya1354_fac_expanded_rates import rate_input as ground_rate_input
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants


def input_text(mult,scale):
    text=ground_input(2,mult,scale)
    if text.count("OptimizeRadial(['ground'])")!=1:raise ValueError('Ground potential insertion differs')
    return text.replace("OptimizeRadial(['ground'])","OptimizeRadial(['ground','odd3s'],[0.5,0.5])",1)


def rate_input(mult,scale,gauge,mapping):
    text,requests,lower=ground_rate_input(2,mult,scale,gauge,mapping)
    original=ground_input(2,mult,scale)
    if not text.startswith(original):raise ValueError('Potential rate prefix differs')
    return input_text(mult,scale)+text[len(original):],requests,lower


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','grasp','isodata','baseline','rate_parent','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic potential controls belong on Sirius')
    for n,v in vars(a).items():setattr(a,n,v.resolve())
    baseline=verify_receipt(a.baseline);parent=verify_receipt(a.rate_parent)
    if (baseline['status']!='FAC_EXPANDED_LEVELS_COMPLETED_REVIEW_REQUIRED' or parent['status']!='FAC_MODEL_RATES_COMPLETED_REVIEW_REQUIRED'
        or baseline['pin']!=PIN or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.grasp,text=True).strip()!=GRASP_PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.grasp,text=True).strip()):raise ValueError('Potential source provenance differs')
    original=a.source/'sfac/sfac';probe=parent['probe'];binary=Path(probe['path'])
    if sha(original)!=baseline['binary_sha256'] or sha(original)!=probe['upstream_binary_sha256'] or sha(binary)!=probe['sha256'] or sha(a.source/'libfac.a')!=probe['libfac_sha256']:raise ValueError('Potential executable/library differs')
    for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
        p=a.baseline/('n3_r2_'+label)
        if (p/'input.sf').read_text()!=ground_input(2,mult,scale):raise ValueError('Fixed-space baseline input differs')
        if not diagnostics((p/'sfac.stdout').read_text(),mult,256)['printed_rule_verified']:raise ValueError('Fixed-space baseline numerical failure')
        inspect((p/'o.lev').read_text(),2)
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    for name in ('doc/manual.tex','sfac/sfac.c','faclib/config.c','faclib/radial.c','faclib/consts.h','faclib/dbase.c','faclib/transition.c'):
        dest=a.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/name,dest)
    for name in ('getmixblock.f90','jj2lsj_code.f90','jj2lsj2K.f90','packLS.f90'):
        dest=a.out/'converter_source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.grasp/'src/appl/jj2lsj90'/name,dest)
    receipt=dict(schema='rya1354.fac-potential.v26',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        fac_pin=PIN,grasp_pin=GRASP_PIN,baseline_receipt_sha256=sha(a.baseline/'receipt.json'),rate_parent_receipt_sha256=sha(a.rate_parent/'receipt.json'),probe=probe,
        cases={},steps=[],scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Fixed788-level n3SD space, explicit0.5/0.5 ground/odd3s optimization weights only changed; cold two numerical cases, full angular replay, bothgauges and fresh complete channel inventory. Shared angular implementation; no accuracy claim.')
    def sfac(exe,out,text):
        (out/'input.sf').write_text(text)
        with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
            process=subprocess.run([str(exe),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=900,check=False)
        if process.returncode:raise ValueError(f'Potential native exit{process.returncode}: '+str(out))
    try:
        for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            case=a.out/label;case.mkdir();levels=case/'levels';levels.mkdir();sfac(original,levels,input_text(mult,scale))
            d=diagnostics((levels/'sfac.stdout').read_text(),mult,256)
            if not d['printed_rule_verified']:raise ValueError('Potential numerical stopping failed')
            inspect((levels/'o.lev').read_text(),2);s=read_files(levels)
            tr=transport(s['o.basis'].decode(),s['o_export_000.c'].decode(),s['o_export_000.cm'],2)
            angular=case/'angular';angular.mkdir()
            for old,new in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:shutil.copyfile(levels/old,angular/new)
            shutil.copyfile(a.isodata,angular/'isodata');step=run(a.grasp/'build/bin/jj2lsj',CONTROLS,angular,900);receipt['steps'].append(step)
            if step['semantic_errors']:raise ValueError('Potential converter semantic refusal')
            c=read_files(angular)
            if c['state.c']!=s['o_export_000.c'] or c['state.cm']!=s['o_export_000.cm']:raise ValueError('Potential angular input changed')
            mapping=angular_review(s,c,2);(case/'mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
            item=dict(diagnostics=d,transport=tr,targets=mapping['targets'],rates={});receipt['cases'][label]=item
            for suffix,gauge in [('B',2),('C',1)]:
                out=case/suffix;out.mkdir();text,requests,lower=rate_input(mult,scale,gauge,mapping);sfac(binary,out,text)
                rd=diagnostics((out/'sfac.stdout').read_text(),mult,256)
                if not rd['printed_rule_verified']:raise ValueError('Potential rate stopping failed')
                for fn in ('o_export_000.c','o_export_000.cm'):
                    if (out/fn).read_bytes()!=s[fn]:raise ValueError('Potential fresh rate state differs')
                reviewed=rate_review((out/'raw.tsv').read_text(),(out/'rates.tr').read_text(),requests,gauge,mapping,constants((a.source/'faclib/consts.h').read_text()))
                item['rates'][suffix]=dict(gauge=gauge,requests=requests,lower_states=lower,diagnostics=rd,audit=reviewed)
            (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_POTENTIAL_CONTROL_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
