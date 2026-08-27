#!/usr/bin/env python3
"""Read-only Elgueta Fe versus CRIRES+ implementation audit (RYA-1060)."""
import csv,json,sys
from collections import defaultdict,Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];SRC=ROOT/'data/reference/elgueta2026_vizier';OUT=ROOT/'data/audit/rya1060_elgueta_fe_gap'
sys.path.insert(0,str(ROOT))
from pipeline.line_match import match
BANDS={'Y':'atomicy.dat','J':'atomicj.dat','H':'atomich.dat'}; GD={'depth':249,'sat':251,'purity':262,'gof':293,'robust':295}
RYA1054={9913.180,9944.207,10142.844,10216.313,10435.355}
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def fl(s,a,b):
 try:return float(s[a:b])
 except:return ''
def near(w,vals,t=.05):return any(abs(w-x)<=t for x in vals)
def write(n,rows,fields=None):
 with (OUT/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 can=read(ROOT/'data/linelists/canonical_gf.csv');idx=defaultdict(list)
 for c in can:
  idx[c['species'].replace(' ','')].append(c)
 rya904={float(r['wavelength_air_A']) for r in read(ROOT/'data/results/rya904/FeI_10280_10680_crires_plus_SYNTH_1D-LTE_lines.csv')}
 rya794={float(r['wave_A']) for r in read(ROOT/'data/results/rya794/vesta_crires_plus_Y_FeI_lines.csv')}
 rows=[]
 for band,name in BANDS.items():
  for raw in (SRC/name).read_text(errors='replace').splitlines():
   if len(raw)!=805:continue
   sp=raw[13:17].strip()
   if sp not in ('FeI','FeII'):continue
   w=fl(raw,0,12);ep=fl(raw,18,27);flags={k:raw[v].strip() for k,v in GD.items()}
   source=[c for c in idx[sp] if c.get('wavelength_air_A') and c.get('excitation_potential_eV')]
   result=match([w],[float(c['wavelength_air_A']) for c in source],want_ep=[ep],src_ep=[float(c['excitation_potential_eV']) for c in source],tol_A=.05,ep_tol_eV=.01)
   unique=result.n_resolved==1 and not result.ambiguous
   c=source[int(result.index[0])] if unique else {}; grade=c.get('gf_tier',''); gradeable=grade=='LAB' or bool(c.get('nist_grade',''))
   main_prod=band=='Y' and near(w,rya904); branch_prod=band=='Y' and near(w,RYA1054); candidate=band=='Y' and near(w,rya794) and not branch_prod
   assessed=flags['robust'] in ('Y','N'); robust=flags['robust']=='Y'; unassessed=not assessed
   if branch_prod:gap='IMPLEMENTED_PRODUCTION'
   elif main_prod or candidate:gap='IMPLEMENTED_CANDIDATE_OR_REJECTED'
   elif band in ('J','H'):gap='NOT_REACHED_OR_UNCONDITIONED'
   elif gradeable:gap='GRADEABLE_NOT_IMPLEMENTED'
   elif robust and c:gap='NOT_GRADEABLE_CURRENTLY'
   elif robust:gap='ELGUETA_ROBUST_NOT_IMPLEMENTED'
   elif assessed:gap='ELGUETA_NONROBUST_NOT_IMPLEMENTED'
   elif unassessed:gap='ELGUETA_UNASSESSED'
   else:gap='OTHER'
   failures=';'.join(k for k,v in flags.items() if k!='robust' and v!='Y') if assessed and not robust else ''
   rows.append({'band':band,'species':sp,'wavelength_A':w,'lower_ep_eV':ep,'elgueta_loggf':fl(raw,28,37),'gd_depth':flags['depth'],'gd_saturation':flags['sat'],'gd_purity':flags['purity'],'gd_gof':flags['gof'],'gd_robust':flags['robust'],'elgueta_status':'robust' if robust else ('assessed_nonrobust' if assessed else 'unassessed'),'failure_reasons':failures,'canonical_match':unique,'canonical_line_id':c.get('line_id',''),'canonical_grade':grade,'canonical_gf_source':c.get('loggf_reference',''),'gf_sigma_dex':c.get('gf_sigma_dex',''),'gradeable_currently':gradeable,'in_current_main_rya904':main_prod,'in_recent_rya1054_commit':branch_prod,'implementation_artifact':('622c954:data/products/solar/Fe.json' if branch_prod else ('data/results/rya904/FeI_10280_10680_crires_plus_SYNTH_1D-LTE_lines.csv' if main_prod else ('data/results/rya794/vesta_crires_plus_Y_FeI_lines.csv' if candidate else ''))),'conditioned_spectrum_reaches':band=='Y','coverage_note':'Y wide derivative exists only on RYA-1054 commit 622c954' if band=='Y' else f'no conditioned {band} derivative in current main','gap_class':gap,'recommended_next_action':'measurement_required' if gradeable and not branch_prod else ('condition_and_measure' if band in ('J','H') else 'preserve_existing')})
 write('elgueta_fe_transition_audit.csv',rows)
 summ=[]
 for key,g in sorted({(r['band'],r['species']):[x for x in rows if (x['band'],x['species'])==(r['band'],r['species'])] for r in rows}.items()):
  co=Counter(x['gap_class'] for x in g)
  summ.append({'band':key[0],'species':key[1],'elgueta_total':len(g),'gd_robust':sum(x['elgueta_status']=='robust' for x in g),'assessed_nonrobust':sum(x['elgueta_status']=='assessed_nonrobust' for x in g),'unassessed':sum(x['elgueta_status']=='unassessed' for x in g),'implemented_production':co['IMPLEMENTED_PRODUCTION'],'implemented_candidate_rejected':co['IMPLEMENTED_CANDIDATE_OR_REJECTED'],'gradeable_not_implemented':sum(x['gradeable_currently'] and not x['in_recent_rya1054_commit'] for x in g),'weak_gf_only':sum(not x['gradeable_currently'] for x in g),'not_reached_unconditioned':co['NOT_REACHED_OR_UNCONDITIONED'],'unexplained_delta':len(g)-sum(co.values())})
 write('band_species_reconciliation.csv',summ)
 write('gradeable_not_implemented.csv',[r for r in rows if r['gradeable_currently'] and not r['in_recent_rya1054_commit']])
 write('robust_blocked_by_conditioning.csv',[r for r in rows if r['gd_robust']=='Y' and r['band'] in ('J','H')])
 write('y_unassessed.csv',[r for r in rows if r['band']=='Y' and r['elgueta_status']=='unassessed'])
 write('attractive_weak_gf.csv',[r for r in rows if r['gd_robust']=='Y' and not r['gradeable_currently']])
 opportunity=[]
 for band in 'YJH':
  g=[r for r in rows if r['band']==band];add=sum(r['gradeable_currently'] and not r['in_recent_rya1054_commit'] for r in g)
  sig=[float(r['gf_sigma_dex']) for r in g if r['gradeable_currently'] and r['gf_sigma_dex'] not in ('',None)]
  opportunity.append({'band':band,'current_production_n':sum(r['in_recent_rya1054_commit'] for r in g),'additional_gradeable_candidates':add,'requiring_new_conditioning_measurement':sum(r['gradeable_currently'] and r['band'] in ('J','H') for r in g),'likely_gf_systematic_floor_dex':min(sig) if sig else 'unknown','principal_risk':'unassessed/telluric/continuum' if band=='Y' else 'missing conditioned derivative plus tellurics','abundance_estimate':'measurement_required'})
 write('opportunity_estimate.csv',opportunity)
 verdict='material gap' if sum(int(x['additional_gradeable_candidates']) for x in opportunity)>10 else 'small gap'
 (OUT/'recommendation.md').write_text(f'# RYA-1060 recommendation: {verdict}\n\nAudit only; no production inputs were changed. Current main and the unmerged RYA-1054 implementation commit are reported separately. Every Fe row has exactly one primary gap class and reconciliation deltas are zero. Candidate abundances remain `measurement_required`.\n')
 (OUT/'summary.json').write_text(json.dumps({'transitions':len(rows),'gap_classes':Counter(r['gap_class'] for r in rows),'verdict':verdict,'rya1054_commit':'622c954','current_main_has_rya1054':False},indent=2)+'\n')
if __name__=='__main__':main()
