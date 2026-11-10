"""Explicit source-precision column ownership for the matched solar experiment."""
import numpy as np


def departure_seed(populations,source_nstar,target_nstar,target_total):
 """Map one positive b-pattern to each local LTE state; initialization only."""
 from pipeline.bride_H_nlte_control import scale_initial_fractions
 p,s,t=[np.asarray(x,float)for x in [populations,source_nstar,target_nstar]]
 if p.shape!=s.shape or s.shape!=t.shape or p.ndim!=2 or any(not np.isfinite(x).all()or np.any(x<=0)for x in [p,s,t]):raise ValueError('Invalid departure-pattern seed')
 return scale_initial_fractions((p/s)*t,target_total)


def column_state(structure):
 from pipeline.bride_binding_controls import BindingControl
 from pipeline.bride_reference import marcs_number_densities
 names=['height_m','temperature_K','electron_pressure_Pa','rho_kg_m3'];n=structure['ndep']
 if n!=56 or any(np.asarray(structure[k]).shape!=(n,)or not np.isfinite(structure[k]).all()for k in names)or any(np.any(np.asarray(structure[k])<=0)for k in names[1:])or not np.all(np.diff(structure['height_m'])<0)or not np.isfinite(structure['vturb_m_s'])or structure['vturb_m_s']<0:raise ValueError('Invalid matched-column fields')
 ne,nH=marcs_number_densities(structure)
 if not np.isfinite(ne).all()or np.any(ne<=0)or not np.isfinite(nH).all()or np.any(nH<=0):raise ValueError('Invalid explicit source-owned densities')
 return BindingControl('model_mass_H',tuple(structure['height_m']),tuple(structure['temperature_K']),tuple(structure['electron_pressure_Pa']),tuple(ne),tuple(nH),tuple(np.full(n,structure['vturb_m_s'])),'source-precision MARCS/interpolated Pe/(LW k*T)','source-precision MARCS/interpolated rho*X/(LW Amu*H mass)')


def prepare_columns(plan,root):
 from pipeline.bride_reference import read_marcs_structure,fingerprint
 from pipeline.bride_parameter_checks import marcs_filename,interpolate_marcs_structure,controlled_structure
 nodes=[];paths=[]
 for node in plan['nodes']:
  path=root/marcs_filename(node['teff'],node['logg'])
  if fingerprint(path)!=node['source']:raise ValueError('Changed predeclared source node')
  nodes.append(((node['teff'],node['logg']),read_marcs_structure(path)));paths.append(path)
 baseline=next(s for p,s in nodes if p==(5750,4.5));probe,probe_receipt=interpolate_marcs_structure(nodes,(5750,4.5),thermodynamic_scale=plan['interpolation']);identity={}
 for key in ['height_m','temperature_K','electron_pressure_Pa','rho_kg_m3','logtau_ross','logtau_500']:
  np.testing.assert_allclose(probe[key],baseline[key],rtol=512*np.finfo(float).eps,atol=1e-7 if key=='height_m'else 1e-12);identity[key]=float(np.max(abs(probe[key]-baseline[key])))
 columns={};receipts={}
 for case in plan['cases']:
  if case['native']:
   if(case['teff'],case['logg'])!=(5750.,4.5):raise ValueError('Unsupported native baseline')
   structure=baseline;receipt=dict(method='native source precision, no interpolation',target=[case['teff'],case['logg']])
  else:structure,receipt=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),thermodynamic_scale=plan['interpolation'])
  structure=controlled_structure(structure,xi=case.get('xi_km_s',plan['xi_km_s']),drop_bottom=0);columns[case['name']]=structure;receipts[case['name']]=receipt
 return columns,dict(sources=[fingerprint(p)for p in paths],case_interpolation=receipts,native_node_max_absolute_identity=identity,native_node_interpolation=probe_receipt)
