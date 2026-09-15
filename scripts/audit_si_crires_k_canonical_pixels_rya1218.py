"""Audit exact pixel reach for canonical Si I/Si II lines in corrected CRIRES+ K."""
from pathlib import Path
import csv,glob,json
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
CENSUS=ROOT/'data/audit/rya1218_si_protocol/canonical_si_census.csv'
OUT=ROOT/'data/audit/rya1218_si_protocol/si_k_canonical_pixel_validation.csv'

def main():
 lines=[]
 for r in csv.DictReader(CENSUS.open()):
  if r['species'] not in ('Si I','Si II'): continue
  w=float(r['wavelength_air_A'])
  if 19450<=w<=24860: lines.append(r)
 rows=[]
 for fn in sorted(glob.glob(str(ROOT/'data/results/rya1219_crires_products/K/*.fits'))):
  with fits.open(fn,memmap=False) as h:
   d=h['SPECTRUM'].data; wave=np.asarray(d['WAVE'],float); flux=np.asarray(d['FLUX'],float); trans=np.asarray(d['MTRANS'],float); det=np.asarray(d['DETEC'],int); applied=bool(h[0].header.get('TELLAPP',False)) and 'MTRANS' in [x.name for x in h]
  for r in lines:
   w=float(r['wavelength_air_A']); m=np.isfinite(wave)&(abs(wave-w)<=1); n=int(m.sum()); ds=sorted(set(det[m].tolist())); finite=bool(n and np.isfinite(flux[m]).all()); nonunity=bool(n and np.any(abs(trans[m]-1)>1e-8)); step=float(np.median(np.diff(np.sort(wave[m])))) if n>2 else 0; tol=max(.08,1.5*step)
   if not applied: st='HOLD_TELLURIC_EVIDENCE'
   elif n==0: st='HOLD_NO_PIXELS'
   elif len(ds)!=1: st='HOLD_DETECTOR_SPLIT'
   elif wave[m].min()>w-1+tol or wave[m].max()<w+1-tol: st='HOLD_TRUNCATED_WINDOW'
   elif not finite: st='HOLD_NONFINITE_FLUX'
   elif not nonunity: st='HOLD_UNITY_MTRANS'
   else: st='REACHED_EXACT_WINDOW'
   rows.append({'product':Path(fn).as_posix(),'species':r['species'],'canonical_line_id':r['line_id'],'wavelength_air_A':r['wavelength_air_A'],'ep_eV':r['excitation_potential_eV'],'loggf':r['log_gf'],'gf_tier':r['gf_tier'],'n_pixels':n,'detectors':';'.join(map(str,ds)),'status':st})
 with OUT.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
 counts={}
 for r in rows: counts[r['status']]=counts.get(r['status'],0)+1
 summary={'ticket':'RYA-1218','scope':'canonical Si I/Si II lines in K corrected products','n_lines':len(lines),'n_products':len(glob.glob(str(ROOT/'data/results/rya1219_crires_products/K/*.fits'))),'n_cells':len(rows),'counts':counts,'measurement_status':'HOLD_NO_ABUNDANCE'}
 (OUT.with_suffix('.json')).write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
