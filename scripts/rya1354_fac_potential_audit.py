"""Independent metadata/raw-rate replay of fixed-space potential control."""
import hashlib
import json
from pathlib import Path
from rya1354_native_parts import load_parts
from rya1354_fac_potential import input_text,rate_input
from rya1354_fac_space import inspect
from rya1354_fac_space_audit import transport
from rya1354_fac_controls import diagnostics
from rya1354_fac_expanded_ls_audit import ls_blocks,label_rows
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants
from rya1354_fac_recouple import CONTROLS


def verify_directory(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    expected={'manifest.json',*manifest['parts']}
    if {p.name for p in directory.iterdir()}!=expected:
        raise ValueError('Unexpected/unlisted potential evidence directory member')


def summarize(base):
    verify_directory(base/'fac_potential_complete_v26')
    files,_=load_parts(base/'fac_potential_complete_v26');root='runs/fac_potential_v26/'
    receipt=json.loads(files[root+'receipt.json']);baseline,_=load_parts(base/'fac_space_complete_v23');parent,_=load_parts(base/'fac_rates_complete_v21')
    if receipt['status']!='FAC_POTENTIAL_CONTROL_COMPLETED_REVIEW_REQUIRED':raise ValueError('Potential control incomplete')
    if receipt['baseline_receipt_sha256']!=hashlib.sha256(baseline['runs/fac_space_v23/receipt.json']).hexdigest() or receipt['rate_parent_receipt_sha256']!=hashlib.sha256(parent['runs/fac_model_rates_v21/receipt.json']).hexdigest():raise ValueError('Potential parent receipt differs')
    if receipt['probe']!=json.loads(parent['runs/fac_model_rates_v21/receipt.json'])['probe']:raise ValueError('Potential frontend provenance differs')
    c=constants(files[root+'source_documents/faclib/consts.h'].decode());result=dict(schema='rya1354.fac-potential-audit.v26',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='Fixed SD space; explicit equal-reference potential policy with retained theory-only ConfigEnergy diagonal-correction workflow. Not pure orbital-only isolation or physical convergence.',cases={})
    manual=files[root+'source_documents/doc/manual.tex'].decode();radial=files[root+'source_documents/faclib/radial.c'].decode();config=files[root+'source_documents/faclib/config.c'].decode()
    if 'configuration group is given an equal weight' not in manual or 'static int _acfg_wmode = 2;' not in radial or 'acfg->weight[i] += fabs(cfg->sweight);' not in config:raise ValueError('Weight-default discrepancy source differs')
    result['weight_default_discrepancy_preserved']=True
    for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
        p=root+label+'/';lev=p+'levels/';angular=p+'angular/';mapping=json.loads(files[p+'mapping.json'])
        if files[lev+'input.sf'].decode()!=input_text(mult,scale) or files[angular+'jj2lsj.stdin'].decode()!=CONTROLS:raise ValueError('Potential declared input differs')
        d=diagnostics(files[lev+'sfac.stdout'].decode(),mult,256)
        if not d['printed_rule_verified']:raise ValueError('Potential source stopping guard failed')
        inspect(files[lev+'o.lev'].decode(),2)
        tr=transport(files[lev+'o.basis'].decode(),files[lev+'o_export_000.c'].decode(),files[lev+'o_export_000.cm'],2)
        if tr!=receipt['cases'][label]['transport']:raise ValueError('Potential transport receipt differs')
        if files[angular+'state.c']!=files[lev+'o_export_000.c'] or files[angular+'state.cm']!=files[lev+'o_export_000.cm']:raise ValueError('Potential angular input differs')
        blocks=ls_blocks(files[angular+'state.lsj.c'].decode(),2);names={r['name'] for b in blocks for r in b};labels=label_rows(files[angular+'state.lsj.lbl'].decode(),788)
        if len(mapping['levels'])!=788 or mapping['full_LS_coefficients_reconstructed']!=63316:raise ValueError('Potential full-replay inventory differs')
        for ident,row in labels.items():
            m=mapping['levels'][str(ident)]
            if row['components'][0][2]!=m['dominant'] or any(name not in names for _,_,name in row['components']):raise ValueError('Potential LS label differs')
            if row['components']!=[(r['coefficient'],r['weight'],r['name']) for r in m['components']]:raise ValueError('Potential archived LS composition differs')
        targets={term:[int(i) for i,r in mapping['levels'].items() if r['dominant']=='1s(2).2s(2).2p(3)4S.3s_'+term] for term in ('3S','5S')}
        if targets!=mapping['targets'] or any(len(v)!=1 for v in targets.values()):raise ValueError('Potential target assignment ambiguous')
        item=dict(diagnostics=d,transport=tr,full_LS_coefficients_reconstructed=mapping['full_LS_coefficients_reconstructed'],max_transformation_orthonormality_residual=mapping['max_transformation_orthonormality_residual'],targets=targets,target_weights={term:mapping['levels'][str(ids[0])]['dominant_printed_weight'] for term,ids in targets.items()},rates={})
        for suffix,gauge in [('B',2),('C',1)]:
            r=p+suffix+'/';text,requests,lower=rate_input(mult,scale,gauge,mapping)
            if files[r+'input.sf'].decode()!=text or not diagnostics(files[r+'sfac.stdout'].decode(),mult,256)['printed_rule_verified']:raise ValueError('Potential rate input/stopping guard differs')
            for fn in ('o_export_000.c','o_export_000.cm'):
                if files[r+fn]!=files[lev+fn]:raise ValueError('Potential regenerated rate state differs')
            audit=rate_review(files[r+'raw.tsv'].decode(),files[r+'rates.tr'].decode(),requests,gauge,mapping,c)
            if audit!=receipt['cases'][label]['rates'][suffix]['audit']:raise ValueError('Potential native/replayed rate audit differs')
            item['rates'][suffix]=audit
        result['cases'][label]=item
    return result

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
