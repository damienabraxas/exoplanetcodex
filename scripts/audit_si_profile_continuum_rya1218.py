#!/usr/bin/env python3
"""Audit profile-fit and continuum evidence for the exact production Si pool."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'data/results/rya1218/abundance_run_20260916/harps/SiI_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_FROMEW_1D-LTE_lines.csv'
d=pd.read_csv(P); d=d[d.in_aggregate].copy()
cols=['wavelength_air_A','abundance','sigma_A','frac_rise_weaker','edge_distance_dex','red_chi2','continuum_level','continuum_method','profile_sigma_A','profile_sigma_floor_A','profile_gamma_A']
rows=d[cols].where(pd.notna(d[cols]),None).to_dict(orient='records')
out={'ticket':'RYA-1218','contract':'RYA-587','pool_species':'Si I','pool_n':len(d),'component_states':{'continuum':{'state':'HOLD','sigma_dex':None,'reason':'No continuum level/sideband perturbation or correlated continuum covariance is recorded on the production lines.'},'profile_ew':{'state':'HOLD','sigma_dex':None,'reason':'Recorded sigma_A and red_chi2 are fit diagnostics, not a calibrated profile likelihood; profile_sigma_A is absent.'}},'diagnostics':rows,'summary':{'median_sigma_A':float(d.sigma_A.median()),'median_red_chi2':float(d.red_chi2.median()),'max_red_chi2':float(d.red_chi2.max()),'n_with_profile_sigma':int(d.profile_sigma_A.notna().sum()),'n_with_continuum_metadata':int(d.continuum_level.notna().sum())},'next_action':'Run matched continuum and profile-likelihood perturbations on Sirius for these seven line windows.'}
(ROOT/'data/results/rya1218/si_profile_continuum_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out['summary'],indent=2))
