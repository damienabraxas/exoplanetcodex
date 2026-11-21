"""Complete full-space native rate audit and paired fixed-potential comparison."""
import hashlib
import json
from pathlib import Path
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_full_rates import rate_input,AUDIT_SHA
from rya1354_fac_no_correction_audit import configuration_keys,verify_config_inventory
from rya1354_fac_controls import diagnostics
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants


def summarize(base):
    for name in ('fac_full_valence_complete_v28','fac_full_angular_complete_v29','fac_full_rates_complete_v30'):verify_directory(base/name)
    levels,_=load_parts(base/'fac_full_valence_complete_v28');angular,_=load_parts(base/'fac_full_angular_complete_v29');rates,_=load_parts(base/'fac_full_rates_complete_v30')
    source='runs/fac_full_valence_v28/';converted='runs/fac_full_recoupling_v29/';root='runs/fac_full_valence_rates_v30/'
    raw_mapping=angular['evidence/fac_full_replay_v29b/report.json'];mapping=json.loads(raw_mapping);receipt=json.loads(rates[root+'receipt.json'])
    if receipt['status']!='FAC_FULL_VALENCE_RATES_COMPLETED_REVIEW_REQUIRED' or hashlib.sha256(raw_mapping).hexdigest()!=AUDIT_SHA or receipt['audit_sha256']!=AUDIT_SHA:raise ValueError('Full rates/mapping not qualified')
    for field,data in [('source_receipt_sha256',levels[source+'receipt.json']),('mapping_receipt_sha256',angular[converted+'receipt.json'])]:
        if receipt[field]!=hashlib.sha256(data).hexdigest():raise ValueError('Full rate parent provenance differs')
    if receipt['probe']!=json.loads(levels[source+'receipt.json'])['probe']:raise ValueError('Full rate frontend differs')
    result=dict(schema='rya1354.fac-full-valence-rates-audit.v30',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,cases={})
    # Same pinned primary constants were retained in the level/source snapshot.
    c=constants(levels[source+'source_documents/faclib/consts.h'].decode())
    for name,m in mapping['cases'].items():
        policy=name.split('_')[0];mult,scale=(.03,100) if name.endswith('100') else (.3,10)
        keys=configuration_keys(levels[source+name+'/o.basis'].decode())
        for suffix,gauge in [('B',2),('C',1)]:
            case=name+'_'+suffix;p=root+case+'/';text,requests,lower=rate_input(policy,mult,scale,gauge,m)
            if rates[p+'input.sf'].decode()!=text or not diagnostics(rates[p+'sfac.stdout'].decode(),mult,256)['printed_rule_verified']:raise ValueError('Full rate declared input/stopping differs')
            for fn in ('o_export_000.c','o_export_000.cm','o.pot'):
                if rates[p+fn]!=levels[source+name+'/'+fn]:raise ValueError('Full rate mapped state/potential differs')
            configs={n:verify_config_inventory(rates[p+'configs.'+n+'.tsv'].decode(),keys) for n in ('before','after')}
            audit=rate_review(rates[p+'raw.tsv'].decode(),rates[p+'rates.tr'].decode(),requests,gauge,m,c)
            if audit!=receipt['cases'][case]['audit']:raise ValueError('Full raw/formatted audit differs')
            result['cases'][case]=dict(configuration_state=configs,audit=audit)
    from rya1354_fac_no_correction_audit import summarize as sd_audit
    sd=sd_audit(base);result['paired_SD_full_comparison_coherent100']={}
    for policy in ('ground','equal'):
        result['paired_SD_full_comparison_coherent100'][policy]={g:dict(SD_lifetime_ns=sd['cases'][policy+'_coherent100']['rates'][g]['audit']['exploratory_model_lifetime_ns'],full_lifetime_ns=result['cases'][policy+'_coherent100_'+g]['audit']['exploratory_model_lifetime_ns']) for g in ('B','C')}
    return result

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
