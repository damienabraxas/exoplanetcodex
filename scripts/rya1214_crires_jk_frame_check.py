#!/usr/bin/env python3
"""RYA-1214 — independent frame check for the J/K rest-frame products.

The conditioner's held-out test cannot catch a frame CONVENTION difference, so this measures
the median core velocity of isolated neutral lines (catalogued depth 0.3-0.9, no catalogued
line deeper than 0.05 within 0.35 A) against catalogue positions on every CRIRES+ holding
and two FTS atlases. The Elgueta Y/H products and the new J/K products should agree with
each other; the atlases keep the observed solar shift (gravitational redshift + convective).
Must run where the holdings load (Sirius).
"""
import sys; from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/"scripts"))
import numpy as np, pandas as pd
from measure_band_ew import load_window_ex
from pipeline import reflected_solar_rv as rrv
LL=pd.read_csv(ROOT/"data/linelists/linelist_solar.csv",low_memory=False,usecols=["element","ion","wavelength_air_A","central_depth"])
def cores(lo,hi,getter):
    s=LL[LL.wavelength_air_A.between(lo,hi)]
    c=s[(s.ion.astype(str).str.strip()=="I")&s.central_depth.between(0.3,0.9)]
    iso=[w for w in c.wavelength_air_A if ((s.wavelength_air_A-w).abs().between(1e-6,0.35)&(s.central_depth>0.05)).sum()==0]
    v=[]
    for w in iso:
        try: x,y=getter(w)
        except Exception: continue
        vv=rrv._core_velocity(x,y,w,win=0.25)
        if np.isfinite(vv): v.append(vv)
    v=np.array(v); return len(iso),len(v),np.median(v) if len(v) else np.nan,np.median(np.abs(v-np.median(v)))*1.4826 if len(v) else np.nan
def hold(inst,h):
    def g(w):
        win=load_window_ex(inst,w,1.0,holding=h); return np.asarray(win.wave,float),np.asarray(win.flux,float)
    return g
def csv(p):
    d=pd.read_csv(p); X,Y=d.wavelength_air_A.to_numpy(),d.flux_normalized.to_numpy()
    return lambda w:(X,Y)
rows=[]
for name,lo,hi,g in [("CRIRES+ Y-wide (Elgueta)",9820,10780,hold("crires_plus","solar_crires_plus_y_wide_rya1054")),
                     ("CRIRES+ H (Elgueta)",15050,17450,hold("crires_plus","solar_crires_plus_h_rya1094")),
                     ("KP-K05 (9820-10000)",9820,10000,hold("kpno_solar_atlas","solar_kpno_kurucz2005_corrected")),
                     ("IAG (9820-11080)",9820,11080,hold("iag_fts_solar_atlas","solar_iag")),
                     ("J rest (this ticket)",11200,13450,csv(ROOT/"data/results/rya1214_crires_jk/solar_crires_plus_j_rya1219_rest.csv")),
                     ("K rest (this ticket)",19500,24800,csv(ROOT/"data/results/rya1214_crires_jk/solar_crires_plus_k_rya1219_rest.csv"))]:
    n_iso,n_m,v,mad=cores(lo,hi,g)
    rows.append(dict(holding=name,lo_A=lo,hi_A=hi,isolated_lines=n_iso,measured=n_m,median_core_v_kms=round(float(v),3),mad_kms=round(float(mad),3)))
    print("%-26s isolated %3d measured %3d  median core v %+.3f  MAD %.3f km/s"%(name,n_iso,n_m,v,mad))
out=ROOT/"data/audit/rya1214_crires_jk/frame_convention_check.csv"
pd.DataFrame(rows).to_csv(out,index=False); print("wrote",out.relative_to(ROOT))
