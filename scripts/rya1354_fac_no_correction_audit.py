"""Metadata and complete raw-rate audit of zero-correction controls."""
import hashlib
import json
from pathlib import Path
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_no_correction import input_text,rate_input,config_state,instrument
from rya1354_fac_space import inspect
from rya1354_fac_space_audit import transport
from rya1354_fac_controls import diagnostics
from rya1354_fac_expanded_ls_audit import ls_blocks,label_rows
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants
from rya1354_fac_recouple import CONTROLS


def configuration_keys(basis):
    keys=set();mode=False
    for line in basis.splitlines():
        if 'Basis Table' in line:mode=True;continue
        if 'Mixing Coefficients' in line:break
        if mode and line.strip() and not line.startswith('#'):
            f=line.split()
            if len(f)!=10:raise ValueError('Configuration basis layout differs')
            keys.add((int(f[4]),int(f[5])))
    if not keys:raise ValueError('No native configuration keys')
    return keys


def verify_config_inventory(text,keys,disabled=True):
    summary=config_state(text,disabled)
    actual={(int(v[0]),int(v[1])) for line in text.splitlines() if line.strip() and not line.startswith('#') for v in [line.split()]}
    if actual!=keys:raise ValueError('Configuration diagnostic inventory differs from native basis')
    return summary


def summarize(base):
    directory=base/'fac_no_correction_complete_v27';verify_directory(directory);files,_=load_parts(directory);root='runs/fac_no_correction_v27/'
    r=json.loads(files[root+'receipt.json']);baseline,_=load_parts(base/'fac_space_complete_v23');enabled,_=load_parts(base/'fac_potential_complete_v26');minimal,_=load_parts(base/'fac_coherent_complete_v18')
    if r['status']!='FAC_NO_CORRECTION_COMPLETED_REVIEW_REQUIRED':raise ValueError('No-correction study incomplete')
    for field,data in [('baseline_receipt_sha256',baseline['runs/fac_space_v23/receipt.json']),('enabled_receipt_sha256',enabled['runs/fac_potential_v26/receipt.json']),('minimal_receipt_sha256',minimal['runs/fac_coherent_controls_v18/receipt.json'])]:
        if r[field]!=hashlib.sha256(data).hexdigest():raise ValueError('No-correction parent provenance differs')
    original_probe=json.loads(enabled['runs/fac_potential_v26/receipt.json'])['probe']
    if any(r['probe'][k]!=original_probe[k] for k in ('libfac_sha256','upstream_binary_sha256')):raise ValueError('No-correction diagnostic library/original frontend provenance differs')
    source=files[root+'source_documents/faclib/radial.c'].decode();cfg=files[root+'source_documents/faclib/config.c'].decode();structure=files[root+'source_documents/faclib/structure.c'].decode()
    if 'if (optimize_control.mce < 0) return 0;' not in source or 'cfg->delta = 0.0;' not in cfg or '*x1 += ci->delta+ci->shift;' not in structure:raise ValueError('No-correction source proof differs')
    result=dict(schema='rya1354.fac-no-correction-audit.v27',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,cases={})
    p=root+'enabled_minimal_canary/';old='runs/fac_coherent_controls_v18/coherent100/'
    if files[p+'input.sf'].decode()!=instrument(minimal[old+'input.sf'].decode()):raise ValueError('Diagnostic canary input differs')
    if not diagnostics(files[p+'sfac.stdout'].decode(),.03,256)['printed_rule_verified']:raise ValueError('Diagnostic canary stopping differs')
    for fn in ('o_export_000.c','o_export_000.cm'):
        if files[p+fn]!=minimal[old+fn]:raise ValueError('Read-only diagnostic changed enabled state')
    positive={n:verify_config_inventory(files[p+'configs.'+n+'.tsv'].decode(),configuration_keys(files[p+'o.basis'].decode()),False) for n in ('before','after')}
    if any(v['mode']!=0 or v['max_abs_delta']==0 for v in positive.values()):raise ValueError('Config diagnostic positive control failed')
    result['enabled_minimal_diagnostic_canary']=positive
    c=constants(files[root+'source_documents/faclib/consts.h'].decode())
    for policy in ('ground','equal'):
        for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            name=policy+'_'+label;p=root+name+'/';lev=p+'levels/';angular=p+'angular/';mapping=json.loads(files[p+'mapping.json'])
            if files[lev+'input.sf'].decode()!=input_text(policy,mult,scale) or files[angular+'jj2lsj.stdin'].decode()!=CONTROLS:raise ValueError('No-correction input differs')
            d=diagnostics(files[lev+'sfac.stdout'].decode(),mult,256)
            if not d['printed_rule_verified']:raise ValueError('No-correction stopping failed')
            inspect(files[lev+'o.lev'].decode(),2)
            tr=transport(files[lev+'o.basis'].decode(),files[lev+'o_export_000.c'].decode(),files[lev+'o_export_000.cm'],2)
            keys=configuration_keys(files[lev+'o.basis'].decode());states={n:verify_config_inventory(files[lev+'configs.'+n+'.tsv'].decode(),keys) for n in ('before','after')}
            if tr!=r['cases'][name]['transport']:raise ValueError('No-correction transport receipt differs')
            for dest,original in [('state.c','o_export_000.c'),('state.cm','o_export_000.cm')]:
                if files[angular+dest]!=files[lev+original]:raise ValueError('No-correction converter input differs')
            blocks=ls_blocks(files[angular+'state.lsj.c'].decode(),2);names={x['name'] for b in blocks for x in b};labels=label_rows(files[angular+'state.lsj.lbl'].decode(),788)
            if len(mapping['levels'])!=788 or mapping['full_LS_coefficients_reconstructed']!=63316:raise ValueError('No-correction full replay count differs')
            for ident,row in labels.items():
                m=mapping['levels'][str(ident)]
                if row['components'][0][2]!=m['dominant'] or any(name not in names for _,_,name in row['components']) or row['components']!=[(x['coefficient'],x['weight'],x['name']) for x in m['components']]:raise ValueError('No-correction full LS metadata differs')
            targets={term:[int(i) for i,x in mapping['levels'].items() if x['dominant']=='1s(2).2s(2).2p(3)4S.3s_'+term] for term in ('3S','5S')}
            if targets!=mapping['targets'] or any(len(v)!=1 for v in targets.values()):raise ValueError('No-correction target assignment ambiguous')
            item=dict(diagnostics=d,transport=tr,configuration_state=states,targets=targets,target_weights={t:mapping['levels'][str(v[0])]['dominant_printed_weight'] for t,v in targets.items()},rates={})
            for suffix,gauge in [('B',2),('C',1)]:
                q=p+suffix+'/';text,requests,lower=rate_input(policy,mult,scale,gauge,mapping)
                if files[q+'input.sf'].decode()!=text or not diagnostics(files[q+'sfac.stdout'].decode(),mult,256)['printed_rule_verified']:raise ValueError('No-correction rate input/stopping differs')
                configs={n:verify_config_inventory(files[q+'configs.'+n+'.tsv'].decode(),keys) for n in ('before','after')}
                for fn in ('o_export_000.c','o_export_000.cm'):
                    if files[q+fn]!=files[lev+fn]:raise ValueError('No-correction regenerated rate state differs')
                audit=rate_review(files[q+'raw.tsv'].decode(),files[q+'rates.tr'].decode(),requests,gauge,mapping,c)
                if audit!=r['cases'][name]['rates'][suffix]['audit']:raise ValueError('No-correction raw/table audit differs')
                item['rates'][suffix]=dict(configuration_state=configs,audit=audit)
            result['cases'][name]=item
    from rya1354_fac_expanded_rate_audit import summarize as enabled_ground_audit
    from rya1354_fac_potential_audit import summarize as enabled_equal_audit
    eg=enabled_ground_audit(base);ee=enabled_equal_audit(base)
    result['fixed_space_comparison_coherent100']={}
    for policy in ('ground','equal'):
        on={g:(eg['cases']['n3_r2_coherent100_'+g]['exploratory_model_lifetime_ns'] if policy=='ground' else ee['cases']['coherent100']['rates'][g]['exploratory_model_lifetime_ns']) for g in ('B','C')}
        off={g:result['cases'][policy+'_coherent100']['rates'][g]['audit']['exploratory_model_lifetime_ns'] for g in ('B','C')}
        result['fixed_space_comparison_coherent100'][policy]=dict(correction_workflow_enabled_ns=on,correction_disabled_verified_zero_ns=off)
    result['enabled_auxiliary_optimization_convergence_qualified']=False
    result['qualification']='Enabled cases are inherited exposed workflow diagnostics; main-potential stopping qualified, but auxiliary ConfigEnergy0 group optimizations occur before coherent settings/printing and were not individually qualified. Disabled cases bypass that path; exact zero delta/shift directly measured. No physical-accuracy or calibrated-uncertainty claim.'
    return result

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
