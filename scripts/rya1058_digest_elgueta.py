#!/usr/bin/env python3
"""Digest the pinned Elgueta+2026 Y/J/H fixed-width tables."""
import csv, hashlib, json, sys
from collections import defaultdict,Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; SRC=ROOT/'data/reference/elgueta2026_vizier'; OUT=ROOT/'data/audit/rya1058_elgueta'
sys.path.insert(0,str(ROOT))
from pipeline.line_match import match
BANDS={'Y':'atomicy.dat','J':'atomicj.dat','H':'atomich.dat'}
TYPES={'procyon_Fd':{'depth':125,'sat':127,'pur':138,'gof':169,'rob':171},'sun_Gd':{'depth':249,'sat':251,'pur':262,'gof':293,'rob':295},'eps_eri_Kd':{'depth':375,'sat':377,'pur':388,'gof':420,'rob':422},'beta_hyi_FGKsg':{'depth':499,'sat':501,'pur':512,'gof':543,'rob':545},'arcturus_FGKg':{'depth':628,'sat':630,'pur':642,'gof':675,'rob':677},'gamma_sge_Mg':{'depth':757,'sat':759,'pur':770,'gof':802,'rob':804}}
def md5(p): return hashlib.md5(p.read_bytes()).hexdigest()
def f(s,a,b):
 try:return float(s[a:b])
 except:return ''
def write(name,rows):
 with (OUT/name).open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def canonical():
 with (ROOT/'data/linelists/canonical_gf.csv').open() as h:return list(csv.DictReader(h))
def main():
 OUT.mkdir(parents=True,exist_ok=True); can=canonical(); rows=[]
 expected={line.split(None,1)[1]:line.split(None,1)[0] for line in (SRC/'MD5SUMS.txt').read_text().splitlines() if line.strip()}
 for name in BANDS.values():
  if expected.get(name)!=md5(SRC/name):raise SystemExit(f'MD5 verification failed: {name}')
 index=defaultdict(list)
 for c in can:
  if c.get('wavelength_air_A') and c.get('excitation_potential_eV'):
   index[c.get('species','').replace(' ','')].append(c)
 for band,name in BANDS.items():
  p=SRC/name; digest=md5(p)
  for raw in p.read_text(errors='replace').splitlines():
   if len(raw)!=805:continue
   wave=f(raw,0,12); species=raw[13:17].strip(); ep=f(raw,18,27); loggf=f(raw,28,37)
   flags={}
   for star,pos in TYPES.items():
    for key,idx in pos.items():flags[f'{star}_{key}']=raw[idx]
   candidates=index[species]
   if wave != '' and ep != '' and candidates:
    result=match([wave],[float(c['wavelength_air_A']) for c in candidates],want_ep=[ep],src_ep=[float(c['excitation_potential_eV']) for c in candidates],tol_A=0.05,ep_tol_eV=0.01)
    unique=result.n_resolved==1 and not result.ambiguous
    ambiguous=result.ambiguous
    hit=candidates[int(result.index[0])] if unique else {}
   else:
    unique=False; ambiguous=False; hit={}
   rows.append({'band':band,'species':species,'element':''.join(x for x in species if x.isalpha()).replace('I',''),'wavelength_A':wave,'lower_ep_eV':ep,'elgueta_loggf':loggf,'c6':f(raw,38,47),**flags,'robust_any':any(raw[p['rob']]=='Y' for p in TYPES.values()),'source_file':name,'source_md5':digest,'doi':'10.1051/0004-6361/202659148','canonical_line_id':hit.get('line_id',''),'canonical_loggf':hit.get('log_gf',''),'canonical_source':hit.get('loggf_reference',''),'canonical_grade':hit.get('gf_tier',''),'delta_loggf':(loggf-float(hit['log_gf'])) if hit and loggf!='' else '','match_status':'PHYSICAL_KEY_UNIQUE' if unique else ('AMBIGUOUS_PHYSICAL_KEY' if ambiguous else 'NO_PHYSICAL_KEY_MATCH'),'primary_source_trace':'UNTRACED_ELGUETA_TABLE_HAS_NO_GF_REFERENCE','action':'TRACE_PRIMARY_GF_BEFORE_PROMOTION'})
 write('normalized_lines.csv',rows)
 inv=[]
 for (band,sp),grp in sorted(defaultdict(list, {k:[r for r in rows if (r['band'],r['species'])==k] for k in {(r['band'],r['species']) for r in rows}}).items()):
  def n(col):return sum(r[col]=='Y' for r in grp)
  fails=Counter()
  for r in grp:
   if r['robust_any']:continue
   for reason,col in [('depth','sun_Gd_depth'),('saturation','sun_Gd_sat'),('purity','sun_Gd_pur'),('gof','sun_Gd_gof')]:
    if r[col]!='Y':fails[reason]+=1
  inv.append({'band':band,'species':sp,'total_transitions':len(grp),'robust_sun':n('sun_Gd_rob'),'robust_procyon':n('procyon_Fd_rob'),'robust_eps_eri':n('eps_eri_Kd_rob'),'robust_any':sum(r['robust_any'] for r in grp),'rejected_all':sum(not r['robust_any'] for r in grp),'dominant_sun_failure':fails.most_common(1)[0][0] if fails else ''})
 write('species_band_inventory.csv',inv)
 for star,col in [('sun','sun_Gd'),('procyon','procyon_Fd'),('eps_eri','eps_eri_Kd')]:
  write(f'{star}_line_behavior.csv',[{'band':r['band'],'species':r['species'],'wavelength_A':r['wavelength_A'],'depth':r[col+'_depth'],'unsaturated':r[col+'_sat'],'purity':r[col+'_pur'],'gof':r[col+'_gof'],'robust':r[col+'_rob'],'canonical_line_id':r['canonical_line_id'],'action':r['action']} for r in rows])
 write('weak_gf_provenance.csv',[r for r in rows if r['robust_any'] and r['primary_source_trace'].startswith('UNTRACED')])
 (OUT/'summary.json').write_text(json.dumps({'lines':len(rows),'bands':Counter(r['band'] for r in rows),'robust_any':sum(r['robust_any'] for r in rows),'unique_physical_matches':sum(r['match_status']=='PHYSICAL_KEY_UNIQUE' for r in rows)},indent=2)+'\n')
if __name__=='__main__':main()
