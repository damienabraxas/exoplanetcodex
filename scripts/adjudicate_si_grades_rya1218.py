#!/usr/bin/env python3
"""Apply the Fe VALD/gf/depth recipe to the complete Si canonical census."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
C=ROOT/'data/audit/rya1218_si_protocol/canonical_si_census.csv'
LL=ROOT/'data/linelists/linelist_full.csv'
OUT=ROOT/'data/audit/rya1218_si_protocol/si_full_grading_matrix.csv'
c=pd.read_csv(C,low_memory=False); l=pd.read_csv(LL,low_memory=False)
c['element']=c.species.astype(str).str.split().str[0]; c['ion_stage']=c.species.astype(str).str.split().str[1]
for d in (c,l):
 d['wavelength_air_A']=pd.to_numeric(d['wavelength_air_A'],errors='coerce'); d['excitation_potential_eV']=pd.to_numeric(d['excitation_potential_eV'],errors='coerce')
ll=l[['element','ion','wavelength_air_A','excitation_potential_eV','central_depth','blend_flag']].copy()
keys=['element','ion','wavelength_air_A','excitation_potential_eV']
def collapse(g):
 vals=sorted(set(float(x) for x in g.central_depth.dropna()))
 return pd.Series({'central_depth': (vals[0] if len(vals)==1 else float('nan')), 'depth_ambiguous': len(vals)>1, 'blend_flag': bool(g.blend_flag.fillna(False).any())})
ll=ll.groupby(keys,as_index=False).apply(collapse,include_groups=False).reset_index()
ll=ll[keys+['central_depth','depth_ambiguous','blend_flag']]
ll=ll[ll.element.eq('Si')].copy()

# Canonical wavelengths are rounded differently from the synthesis linelist
# (typically by 1e-5--4e-4 A).  An exact four-column join silently converted
# most known `in_linelist` rows into DEPTH_UNMEASURED.  Match the physical
# transition with explicit tolerances, then preserve ambiguity when more than
# one candidate survives.  This is the same validate-don't-tune principle as
# the Fe recipe; the tolerances only absorb printed precision.
matched=[]
for _, r in c.iterrows():
    q=ll[(ll.element==r.element)&(ll.ion==r.ion_stage)]
    q=q[(q.wavelength_air_A-r.wavelength_air_A).abs()<=0.005]
    q=q[(q.excitation_potential_eV-r.excitation_potential_eV).abs()<=0.01]
    if q.empty:
        matched.append({'central_depth':float('nan'),'depth_ambiguous':False,
                        'blend_flag':False})
        continue
    # Select the nearest physical key; retain a depth ambiguity if equivalent
    # candidates carry different depths.
    q=q.assign(_distance=(q.wavelength_air_A-r.wavelength_air_A).abs() +
               0.1*(q.excitation_potential_eV-r.excitation_potential_eV).abs())
    best=q[q._distance <= q._distance.min()+1e-9]
    depths=sorted(set(float(x) for x in best.central_depth.dropna()))
    matched.append({'central_depth': depths[0] if len(depths)==1 else float('nan'),
                    'depth_ambiguous': len(depths)>1,
                    'blend_flag': bool(best.blend_flag.fillna(False).any())})
j=c.reset_index(drop=True).copy()
j=pd.concat([j, pd.DataFrame(matched)], axis=1)
# DH23 Table 3 supplies upper/lower energies and J for the two matched
# transitions.  The canonical census predates those columns, so its inherited
# HOLD marker must not veto a source-level identity that is now resolved.
dh23 = j.DH23_matches.astype(str).str.strip().isin(['1']) & j.wavelength_air_A.round(3).isin([3905.523, 4102.936])
j.loc[dh23, 'physical_identity_status'] = 'RESOLVED_DH23_TABLE3_levels_and_J'
def depth_route(d):
 if isinstance(d, tuple): return 'DEPTH_UNMEASURED'
 if pd.isna(d): return 'DEPTH_UNMEASURED'
 if d < .05: return 'BELOW_OBSERVABILITY_FLOOR'
 if d <= .60: return 'CODEX_DEPTH_WINDOW'
 return 'DEEP_DEPTH_ROUTE'
def gf_class(r):
 if str(r['DH23_matches']).strip() not in ('','0','nan'):
  return 'PRIMARY_LAB_CANDIDATE_DH23'
 return 'systematic:K07' if str(r['gf_tier']).upper() in ('VALD3','KURUCZ','') else 'UNRESOLVED_GF_SOURCE'
def grade(r):
 dr=r['depth_route']; gc=r['gf_quality_class']; ident=str(r['physical_identity_status'])
 if gc.startswith('PRIMARY_LAB') and 'HOLD' in ident: return 'LAB_CANDIDATE_TRANSITION_HOLD'
 if gc.startswith('PRIMARY_LAB') and dr=='CODEX_DEPTH_WINDOW': return 'CODEX_GRADED_CANDIDATE'
 if gc.startswith('PRIMARY_LAB') and dr=='DEEP_DEPTH_ROUTE': return 'DEEP_GRADED_CANDIDATE'
 if dr=='CODEX_DEPTH_WINDOW': return 'UNGRADED_CODEX_DEPTH_DIAGNOSTIC'
 if dr=='DEEP_DEPTH_ROUTE': return 'UNGRADED_DEEP_DEPTH_DIAGNOSTIC'
 return dr
j['gf_quality_class']=j.apply(gf_class,axis=1); j['solar_feature_depth']=j['central_depth']; j.loc[j.depth_ambiguous.fillna(False),'depth_route']='DEPTH_AMBIGUOUS'; j.loc[~j.depth_ambiguous.fillna(False),'depth_route']=j.loc[~j.depth_ambiguous.fillna(False),'solar_feature_depth'].map(depth_route); j['recipe_grade']=j.apply(grade,axis=1)
j['reference_role']=j.external_membership.fillna('not_in_reference_census'); j['transition_match_required']='yes'; j['gf_action']=j.gf_quality_class.map(lambda x:'charge cited laboratory sigma when transition identity is confirmed' if x.startswith('PRIMARY_LAB') else 'charge systematic:K07; retain line in ungraded pool')
cols=['line_id','band','species','wavelength_air_A','excitation_potential_eV','log_gf','loggf_reference','gf_tier','gf_quality_class','gf_sigma_dex','gf_source_doi','DH23_matches','physical_identity_status','solar_feature_depth','depth_ambiguous','depth_route','recipe_grade','reference_role','blend_flag','transition_match_required','gf_action']
j[cols].to_csv(OUT,index=False)
print('rows',len(j)); print(j.recipe_grade.value_counts().to_string()); print('\nby band'); print(j.assign(band=j.band).groupby(['band','recipe_grade'],dropna=False).size().to_string())
