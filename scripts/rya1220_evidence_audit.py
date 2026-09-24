#!/usr/bin/env python3
"""Read-only audit of the exact five N I indicators and existing CN products."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Transcribed from the user-supplied publisher PDF, Tables 1, 2, 3.
# pm -> mA = 10; air nm -> A = 10. Values are benchmarks, never fit inputs.
NI = [
    (7442.29,10.330,-.403,.310,.075,.235,7.767,7.779,7.808,7.818,7.798,7.806),
    (8216.33,10.336,.138,.770,0.,.770,7.767,7.782,7.810,7.824,7.796,7.806),
    (8629.23,10.690,.077,.620,.210,.410,7.759,7.772,7.815,7.826,7.803,7.812),
    (8683.40,10.330,.106,.865,.115,.750,7.768,7.784,7.808,7.823,7.794,7.805),
    (10108.9,11.753,.444,.275,.075,.200,7.772,7.783,7.838,7.848,7.827,7.835),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--prior', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    n = pd.DataFrame(NI, columns=['wavelength_air_A','lower_EP_eV','published_loggf',
        'feature_ew_pm','cn_blend_ew_pm','ni_ew_pm','A_3d_nlte','A_3d_lte',
        'A_mean3d_nlte','A_mean3d_lte','A_1d_nlte','A_1d_lte'])
    n['geometry'] = 'disk_center_intensity'
    n['source'] = 'amarsi2020_solar_n: Tables 1-3'
    n['ni_ew_mA'] = n.ni_ew_pm*10
    n['delta_3d_minus_1d_nlte'] = n.A_3d_nlte-n.A_1d_nlte
    n['delta_nlte_minus_lte_1d'] = n.A_1d_nlte-n.A_1d_lte
    canon_path=ROOT/'data/linelists/canonical_gf.csv'
    canon=pd.read_csv(canon_path)
    canon=canon[canon.species=='N I']
    rows=[]
    for r in n.to_dict('records'):
        m=canon[(abs(canon.wavelength_air_A-r['wavelength_air_A'])<=.05)
                & (abs(canon.excitation_potential_eV-r['lower_EP_eV'])<=.005)]
        if len(m)!=1:
            r.update(canonical_match_status=f'REVIEW_{len(m)}_matches')
        else:
            q=m.iloc[0]
            r.update(canonical_match_status='unique_wavelength_EP',physical_id=q.physical_id,
                     canonical_loggf=q.log_gf,gf_tier=q.gf_tier,nist_grade=q.nist_grade,
                     gf_reference=q.loggf_reference,gf_source_doi=q.gf_source_doi,
                     gf_difference=float(q.log_gf-r['published_loggf']))
        rows.append(r)
    pd.DataFrame(rows).to_csv(a.out/'ni_five_line_reference.csv',index=False)
    previous=[]
    for p in sorted(a.prior.glob('NI_*_lines.csv')):
        d=pd.read_csv(p)
        for r in n.to_dict('records'):
            m=d[(abs(d.wavelength_air_A-r['wavelength_air_A'])<=.05)
                & (abs(d.ep_eV-r['lower_EP_eV'])<=.005)]
            for q in m.to_dict('records'):
                previous.append({**q,'source_file':str(p),
                    'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
                    'literature_1d_lte_intensity':r['A_1d_lte'],
                    'geometry_comparison':'flux versus disk-center intensity; NOT matched',
                    'disposition':'diagnostic_only_pending_reproduction_and_full_budget'})
    pd.DataFrame(previous).to_csv(a.out/'ni_prior_line_audit.csv',index=False)
    mol_path=ROOT/'data/reference/amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv'
    mol=pd.read_csv(mol_path)
    mol=mol[mol.element_parameter=='logepsN'].copy()
    from pipeline.wavelength_util import vac_to_air
    mol['wavelength_air_A']=vac_to_air(mol.wavelength_vac_nm.to_numpy()*10)
    mol['delta_3d_minus_marcs']=mol.abundance_3d-mol.abundance_marcs
    mol['geometry']='disk_center_intensity'
    mol.to_csv(a.out/'molecular_reference_lines.csv',index=False)
    groups=mol.groupby(['species','system','delta_nu']).agg(
        n=('source_row','size'),vac_nm_min=('wavelength_vac_nm','min'),
        vac_nm_max=('wavelength_vac_nm','max'),A_3d_mean=('abundance_3d','mean'),
        A_marcs_mean=('abundance_marcs','mean'),
        delta_3d_minus_marcs_mean=('delta_3d_minus_marcs','mean'))
    groups.to_csv(a.out/'molecular_groups.csv')
    budget=[]
    from pipeline.nitrogen_uncertainty import REQUIRED_TERMS
    for p in sorted((ROOT/'data/audit/cno_synthesis').glob('solar*cno_provenance.json')):
        q=json.loads(p.read_text());u=q.get('uncertainty',{}).get('N',{})
        for term in REQUIRED_TERMS:
            budget.append(dict(product=p.name,term=term,
                state='reported_partial' if term=='statistical' and u.get('stat') is not None else 'not_demonstrated',
                reported_stat=u.get('stat'),reported_sys=u.get('sys'),reported_total=u.get('tot'),
                publishable=False))
    pd.DataFrame(budget).to_csv(a.out/'prior_budget_adjudication.csv',index=False)
    manifest={'ticket':'RYA-1220','calibration_applied':False,
        'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [canon_path,mol_path]},
        'ni_reference_count':len(n),'ni_prior_rows':len(previous),
        'molecular_n_reference_count':len(mol),'publication_status':'HELD',
        'limitations':['Existing products are not new reproduction runs.',
            'N I 1D comparison retains intensity/flux mismatch.',
            'Molecular table lacks rotational quantum IDs; identities need source-list crossmatch.',
            'No aggregate absolute-error quadrature is used as a differential error.']}
    (a.out/'audit_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))


if __name__=='__main__': main()
