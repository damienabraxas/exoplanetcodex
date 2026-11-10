"""Explicit fixed masks and no-clipping continuum controls; diagnostic only."""
import numpy as np
from pipeline.bride_continuum_audit import support_geometry


def shared_masks(wave,centre,synthesis,iag,*,inner,outer,percentile):
 w,s,i=[np.asarray(x,float)for x in (wave,synthesis,iag)]
 if w.ndim!=1 or s.shape!=w.shape or i.shape!=w.shape or not np.all(np.diff(w)>0)or not 0<inner<outer or not 0<percentile<100:raise ValueError('Invalid shared support layout')
 eligible=np.isfinite(w)&np.isfinite(s)&np.isfinite(i)&(s>0)&(i>0)&(abs(w-centre)>inner)&(abs(w-centre)<=outer)
 if not eligible.any():raise ValueError('No eligible shared support')
 intersection=eligible&(s>=np.percentile(s[eligible],percentile))&(i>=np.percentile(i[eligible],percentile));balanced=np.zeros(len(w),bool)
 for side in (w<centre,w>centre):
  side=eligible&side
  if side.any():balanced|=side&(i>=np.percentile(i[side],percentile))
 return dict(shared_top_decile_intersection=intersection,balanced_iag_per_side=balanced)


def fixed_fit(wave,observed,centre,reference,mask,*,maximum_shift,minimum_pixels,minimum_each_side,require_bracket):
 w,f,m=[np.asarray(x,float)for x in (wave,observed,reference)];mask=np.asarray(mask)
 if w.ndim!=1 or f.shape!=w.shape or m.shape!=w.shape or mask.shape!=w.shape or mask.dtype!=bool or not np.all(np.diff(w)>0)or not np.isfinite(maximum_shift)or maximum_shift<=0 or minimum_pixels<2 or minimum_each_side<1:raise ValueError('Invalid fixed-support contract')
 indices=np.flatnonzero(mask);blue=int(np.sum(w[indices]<centre));red=int(np.sum(w[indices]>centre));record=dict(applied=False,reason='',indices=indices.tolist(),count=len(indices),blue_count=blue,red_count=red,require_bracket=require_bracket,clip=False,geometry=None,level=None,slope_per_A=None)
 if len(indices)<minimum_pixels:record['reason']='REFUSED_INSUFFICIENT_PIXELS';return record,None
 if require_bracket and(blue<minimum_each_side or red<minimum_each_side):record['reason']='REFUSED_BRACKETING';return record,None
 if not np.isfinite(f[indices]).all()or not np.isfinite(m[indices]).all()or np.any(m[indices]<=0):raise ValueError('Invalid fixed-support values')
 geometry=support_geometry(w,centre,indices);record['geometry']={k:v for k,v in geometry.items()if k!='centre_weights'}
 if require_bracket and geometry['centre_outside_hull_A']>0:raise ValueError('Bracketing geometry contradiction')
 c=np.polyfit(w[indices]-centre,f[indices]/m[indices],1);level=float(c[-1]);slope=float(c[0]);record.update(level=level,slope_per_A=slope)
 estimate=float(np.dot(geometry['centre_weights'],f[indices]/m[indices]));bound=64*np.finfo(float).eps*geometry['design_condition_2']*max(1.,float(np.max(abs(f[indices]/m[indices]))))
 if not np.isfinite(c).all()or abs(estimate-level)>bound:raise ValueError('Fixed-support coefficient reconstruction failed')
 if abs(level-1)>maximum_shift:record['reason']='REFUSED_PINNED_LEVEL_BOUND';return record,None
 cont=np.polyval(c,w-centre)
 if not np.isfinite(cont).all()or np.any(cont<=0):record['reason']='REFUSED_NONPOSITIVE_CONTINUUM';return record,None
 record.update(applied=True,reason='applied');return record,cont
