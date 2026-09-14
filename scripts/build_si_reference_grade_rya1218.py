#!/usr/bin/env python3
"""Build the Si external-reference ledger and keep Reference Grade separate from gf tiers.

Reference Grade is membership in a registered external published line set.  It is orthogonal
 to primary-laboratory gf quality and to Codex/Deep measurement tiers.
"""
from pathlib import Path
import csv, json
ROOT=Path(__file__).resolve().parents[1]
CAN=ROOT/'data/audit/rya1218_si_protocol/canonical_si_census.csv'
AG=ROOT/'data/reference/si_agss21/si_agss21_lines.csv'
BJ=ROOT/'data/reference/si_bergemann2013/si_bergemann2013_j_lines.csv'
LAB=ROOT/'data/audit/rya1169_si_intake/si_primary_lab_gf_census.csv'
OUT=ROOT/'data/audit/rya1218_si_protocol/si_reference_grade_ledger.csv'

def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))

def key(w,ep): return (round(float(w),3),round(float(ep),3))
can=rows(CAN)
by_id={r['line_id']:r for r in can}
# AGSS rows already carry canonical IDs and source inclusion/exclusion.
for r in rows(AG):
    cid=r.get('canonical_line_id','').strip()
    if not cid or cid not in by_id: continue
    c=by_id[cid]; s=r['selection_status']
    c.setdefault('reference_sets','')
    if s=='USED_BY_SOURCE_ANALYSIS':
        c['reference_grade']='Reference Grade'; c['reference_sets']='si-agss21'
        c['reference_grade_basis']='AGSS21 Si lineage; source-selected row'
    else:
        c['reference_grade']='Reference Excluded'; c['reference_sets']='si-agss21'
        c['reference_grade_basis']='AGSS21 lineage row explicitly excluded by source analysis'
# Bergemann identity crosswalk is the authoritative canonical ID match for J.
for r in rows(BJ):
    hits=[c for c in can if c['species'].strip()=='Si I' and abs(float(c['wavelength_air_A'])-float(r['wavelength_air_A']))<=0.01 and abs(float(c['excitation_potential_eV'])-float(r['ep_eV']))<=0.01]
    if len(hits)!=1: continue
    c=hits[0]; c['reference_grade']='Reference Grade'
    sets=[x for x in c.get('reference_sets','').split('|') if x]
    if 'si-bergemann-j' not in sets: sets.append('si-bergemann-j')
    c['reference_sets']='|'.join(sets)
    c['reference_grade_basis']='Bergemann2013 Table 1 published J-band line'
# all remaining rows are explicitly pending, never silently downgraded.
for c in can:
    c.setdefault('reference_sets','')
    c.setdefault('reference_grade_basis','')
    if not c.get('reference_grade') or c['reference_grade'].strip() in {'','HOLD_best_available_per_line_adjudication'}:
        c['reference_grade']='Reference Pending'
        c['reference_grade_basis']='No registered external reference membership yet'
# Literature readiness and measurement status are separate from reference membership.
for name in ('primary_lab_status','primary_lab_source','primary_lab_loggf','primary_lab_sigma_dex',
             'measurement_status','measurement_reason'):
    for c in can: c.setdefault(name,'')
for r in rows(LAB):
    cid=(r.get('canonical_line_id') or '').strip()
    if not cid or cid not in by_id: continue
    c=by_id[cid]
    c['primary_lab_status']='MATCHED_PRIMARY_LAB' if r.get('join_status')=='matched' else 'LITERATURE_ONLY'
    c['primary_lab_source']=r.get('primary_lab_source','')
    c['primary_lab_loggf']=r.get('experimental_loggf','')
    c['primary_lab_sigma_dex']=r.get('lab_sigma_dex','')
    c['measurement_status']=r.get('measurement_status','NOT_MEASURED')
    c['measurement_reason']=r.get('measurement_reason','')
for c in can:
    if not c['measurement_status']: c['measurement_status']='NOT_MEASURED'
    if not c['measurement_reason']: c['measurement_reason']='No valid completed measurement linked yet'
# write ledger in canonical order with explicit axis names
fields=list(can[0])
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(can)
# summary
from collections import Counter
summ={'ticket':'RYA-1218','reference_sets':{'si-agss21':{'used':sum('si-agss21' in c['reference_sets'] and c['reference_grade']=='Reference Grade' for c in can),'excluded':sum('si-agss21' in c['reference_sets'] and c['reference_grade']=='Reference Excluded' for c in can)},'si-bergemann-j':{'used':sum('si-bergemann-j' in c['reference_sets'] for c in can)}},'reference_grade_counts':dict(Counter(c['reference_grade'] for c in can)),'definition':'Reference Grade is external published line-set membership; Codex/Deep remain independent measurement/gf axes.'}
(ROOT/'data/audit/rya1218_si_protocol/si_reference_grade_summary.json').write_text(json.dumps(summ,indent=2)+'\n')
print(json.dumps(summ,indent=2))
