"""Cross-check the published Bergemann J-band Si I lines against canonical Si rows."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
POOL=ROOT/'data/audit/rya1218_si_protocol/si_nir_line_pool.csv'
CENSUS=ROOT/'data/audit/rya1218_si_protocol/canonical_si_census.csv'
OUT=ROOT/'data/audit/rya1218_si_protocol/si_j_line_identity.csv'

def main():
 canon=list(csv.DictReader(CENSUS.open())); pool=list(csv.DictReader(POOL.open())); rows=[]
 for p in pool:
  w=float(p['wavelength_air_A']); c=min((r for r in canon if r['species']=='Si I'),key=lambda r:abs(float(r['wavelength_air_A'])-w))
  rows.append({'external_source':'Bergemann2013 Table 1','external_line':p['lower_level']+' -> '+p['upper_level'],'external_wavelength_air_A':p['wavelength_air_A'],'external_ep_eV':p['ep_eV'],'external_loggf':p['loggf'],'canonical_line_id':c['line_id'],'canonical_wavelength_air_A':c['wavelength_air_A'],'canonical_ep_eV':c['excitation_potential_eV'],'canonical_loggf':c['log_gf'],'delta_wavelength_mA':round((float(c['wavelength_air_A'])-w)*1000,3),'delta_loggf_dex':round(float(c['log_gf'])-float(p['loggf']),6),'canonical_gf_reference':c['loggf_reference'],'canonical_source':c['seed_source'],'canonical_physical_identity_status':c['physical_identity_status'],'reference_grade':c['reference_grade'],'codex_grade':c['codex_grade'],'deep_grade':c['deep_grade'],'identity_disposition':'EXTERNAL_LEVEL_IDENTITY_RECORDED_CANONICAL_LEVEL_FIELDS_MISSING','gf_disposition':'HOLD_VALD_ONLY_NO_PRIMARY_LAB_PROVENANCE'})
 with OUT.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
 print(f'wrote {OUT} ({len(rows)} rows)')
if __name__=='__main__': main()
