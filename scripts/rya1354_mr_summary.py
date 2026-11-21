"""Replay native MR inventory and optionally completed lifetime experiments."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_lifetime_summary import load_parts, decay, scf_state
from rya1354_lifetime_ladder import STATES, identity
from rya1354_mr_inventory import inventory, references
from rya1354_triplet_correlation import require_triplet_quintet_order

ROOT=Path(__file__).resolve().parents[1]


def csf_records(text):
    lines=[line for line in text.splitlines()[5:] if not line.lstrip().startswith('*')]
    if len(lines)%3:
        raise ValueError('Incomplete CSF triples')
    return {tuple(lines[i:i+3]) for i in range(0,len(lines),3)}


def summarize_inventory(files):
    root='runs/mr_inventory_n4_v1'
    receipt=json.loads(files[root+'/receipt.json'])
    if receipt['status']!='INVENTORY_COMPLETED_NO_ATOMIC_PREDICTION':
        raise ValueError('Inventory not completed')
    cases={}
    records={}
    for name,case in receipt['cases'].items():
        parity,scope,reference=name.split('_')
        core=scope=='CV';mr=reference=='MR'
        text=files[root+'/'+name+'/rcsf.out'].decode()
        actual=inventory(text,parity,core)
        if any(case[k]!=v for k,v in actual.items()) or case['references']!=references(parity,core,mr):
            raise ValueError('Declared MR inventory differs from actual native records')
        expected={'0+','1+','2+'} if parity=='even' else {'1-','2-'}
        if set(actual['blocks'])!=expected:
            raise ValueError('Wrong J-block inventory')
        cases[name]=case
        records[name]=csf_records(text)
    expansions={}
    for parity in ('even','odd'):
        for scope in ('VV','CV'):
            prefix=parity+'_'+scope
            a,b=records[prefix+'_SR'],records[prefix+'_MR']
            if not a < b:
                raise ValueError('MR is not a strict CSF superset of its matched SR')
            expansions[prefix]=dict(SR_CSFs=len(a),MR_CSFs=len(b),additional_CSFs=len(b-a),SR_subset_MR=True)
    return dict(cases=cases,expansions=expansions,status=receipt['status'])


def fixed_cases(files,case_specs,root):
    ladder,_=load_parts(ROOT/'data/audit/rya1354/lifetime_ladder_v4')
    old,_=load_parts(ROOT/'data/audit/rya1354/lifetime_aux_v4b')
    generated,_=load_parts(ROOT/'data/audit/rya1354/mr_inventory_v5')
    cases={}
    for spec in case_specs:
        case=spec['name']; selected={}
        for state in STATES:
            path=root+'/'+case+'_'+state
            seed='runs/lifetime_ladder_v3/n4_'+state+'_verify'
            for name in ('state.w','isodata'):
                if files[path+'/'+name]!=ladder[seed+'/'+name]:
                    raise ValueError('Fixed radial/nucleus identity differs')
            if case=='DC' and files[path+'/state.c']!=ladder[seed+'/state.c']:
                raise ValueError('SR control CSFs differ')
            if case=='CV_2s_MR_Breit' and files[path+'/state.c']!=files[root+'/CV_2s_MR_DC_'+state+'/state.c']:
                raise ValueError('Core MR DC/Breit spaces differ')
            text=files[path+'/state.lsj.lbl'].decode();identity(text,state);selected[state]=text
            inventory(files[path+'/state.c'].decode(),'even' if state=='ground' else 'odd',spec['core'])
            if 'MR' in case:
                parity='even' if state=='ground' else 'odd'
                scope='CV' if spec['core'] else 'VV'
                expected=csf_records(generated[f'runs/mr_inventory_n4_v1/{parity}_{scope}_MR/rcsf.out'].decode())
                if state!='ground':
                    expected={r for r in expected if r[2].split()[-1]==('2-' if state=='quintet' else '1-')}
                if csf_records(files[path+'/state.c'].decode())!=expected:
                    raise ValueError('MR state CSFs differ from independently audited generator inventory')
        require_triplet_quintet_order(selected['triplet'],selected['quintet'])
        cases[case]=decay(files,root+'/'+case+'_',ci=True)
        if len(cases[case]['channels'])!=13:
            raise ValueError('Incomplete lifetime channel inventory')
    previous=decay(old,'runs/lifetime_core_n4_v1/DC_',ci=True)
    if cases['DC']!=previous:
        raise ValueError('Identical SR/DC model does not reproduce prior native output')
    return cases


def summarize(mr_bundle=None,relax_bundle=None):
    inventory_files,manifest=load_parts(ROOT/'data/audit/rya1354/mr_inventory_v5')
    result=dict(schema='rya1354.mr-lifetime-continuation.v5',exposure='EXPOSED_NONBLIND',
                scientific_admission='HOLD',adoption_allowed=False,uncertainty=None,covariance=None,
                calibrated_sigma=None,chi_square=None,candidate_doi=None,
                prediction_origin='THIS_UNPUBLISHED_EXPLORATORY_CALCULATION',
                generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                inventory=summarize_inventory(inventory_files),
                inventory_parts=manifest['parts'],fixed_orbital_MR=None,relaxed_n4=None,
                scope='Inventory alone is not an atomic prediction. Only separately verified completed subcases may supply exploratory lifetimes; pending cases remain null. No published-protocol reproduction or independent experimental review claimed.')
    partial,partial_manifest=load_parts(ROOT/'data/audit/rya1354/mr_partial_v5')
    checkpoint=json.loads(partial['evidence/mr_partial_checkpoint_v5.json'])
    if (checkpoint['cases']!=['DC','VV_MR_DC'] or checkpoint['scope']!=
            'POSTHOC_SNAPSHOT_OF_COMPLETED_SUBCASES_ONLY_NOT_FINAL_CAMPAIGN_RECEIPT'):
        raise ValueError('Unexpected partial-campaign checkpoint scope')
    for name,record in checkpoint['files'].items():
        if name in partial:
            raw=partial[name]
            if hashlib.sha256(raw).hexdigest()!=record['sha256'] or len(raw)!=record['bytes']:
                raise ValueError('Completed subcase differs from checkpoint')
        elif name in partial_manifest['omitted_angular_caches'] and Path(name).name.startswith('mcp.'):
            continue
        elif name in partial_manifest['external_restart_caches'] and Path(name).name=='rci.res':
            external=partial_manifest['external_restart_caches'][name]
            if any(external[k]!=record[k] for k in ('sha256','bytes')):
                raise ValueError('Completed subcase restart-cache declaration differs')
        else:
            raise ValueError('Completed subcase artifact absent from checkpoint archive')
    cases=fixed_cases(partial,[dict(name='DC',core=False),dict(name='VV_MR_DC',core=False)],
                      'runs/lifetime_mr_n4_v1')
    result['completed_MR_subcases']=dict(scope=checkpoint['scope'],
        captured_utc=checkpoint['captured_utc'],campaign_status_at_capture=checkpoint['campaign_status_at_capture'],
        cases=cases,parts=partial_manifest['parts'],same_space_DC_reproduced=True,
        pending='Full2s MR DC/Breit and separate n4 relaxation are not supplied by this partial checkpoint')
    if mr_bundle:
        files,parts=load_parts(mr_bundle)
        root='runs/lifetime_mr_n4_v1'
        receipt=json.loads(files[root+'/receipt.json'])
        if receipt['status']!='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED' or not receipt['multireference'] or receipt['n']!=4:
            raise ValueError('MR lifetime campaign not completed/declared')
        cases=fixed_cases(files,receipt['case_specs'],root)
        result['fixed_orbital_MR']=dict(cases=cases,parts=parts['parts'],
                                      same_space_DC_reproduced=True,
                                      native_numerical_caveat='Stationary-nucleus RCI warnings retained; no isotope precision or calibrated uncertainty claimed')
    if relax_bundle:
        files,parts=load_parts(relax_bundle);root='runs/lifetime_relax_n4_v1'
        receipt=json.loads(files[root+'/receipt.json'])
        if receipt['status']!='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED' or receipt['n']!=4:
            raise ValueError('n4 relaxation not completed/qualified')
        states={}
        for state in STATES:
            for suffix in ('','_verify'):
                path=root+'/'+state+suffix
                states[state+suffix]=scf_state(files,path,state)
            if files[root+'/'+state+'/state.w']!=files[root+'/'+state+'_verify/rwfn.inp']:
                raise ValueError('Fresh restart did not use the optimized wavefunction')
        require_triplet_quintet_order(files[root+'/triplet_verify/state.lsj.lbl'].decode(),files[root+'/quintet_verify/state.lsj.lbl'].decode())
        result['relaxed_n4']=dict(states=states,lifetime=decay(files,root+'/'),parts=parts['parts'])
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mr-bundle',type=Path)
    ap.add_argument('--relax-bundle',type=Path)
    ap.add_argument('--out',type=Path,default=ROOT/'data/audit/rya1354/source_curation_v2/mr_summary_v5.json')
    args=ap.parse_args()
    args.out.write_text(json.dumps(summarize(args.mr_bundle,args.relax_bundle),indent=2)+'\n')
