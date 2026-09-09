"""RYA-1204 Lever B, measured the way production measures: run the RYA-759 near-UV Fe
product with molecular opacity OFF (today) and ON, and report the shift in A(Fe).

Only two things are patched, both narrowly: `use_molecules=True` on the synthesis call,
and the molecules directory iSpec globs. Everything else -- atlas, model, atomic list,
window rule, fitter -- is the production script, unmodified.
"""
import os, sys, numpy as np
SP=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,'/Users/ryanschmitt/codex/ispec'); sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195')
sys.path.insert(0,'/Users/ryanschmitt/codex/rya1195/scripts')
USE_MOL = os.environ.get('USE_MOL','0')=='1'
MOL=os.path.join(SP,'molecules')
_orig_symlink=os.symlink
def symlink(src,dst,*a,**k):
    if str(dst).endswith('/molecules'): src=MOL
    return _orig_symlink(src,dst,*a,**k)
os.symlink=symlink
import ispec
_gen=ispec.generate_spectrum
def gen(*a,**k):
    if USE_MOL: k['use_molecules']=True
    return _gen(*a,**k)
ispec.generate_spectrum=gen
import pipeline.abundances_derive as ad
ad.ispec.generate_spectrum=gen

sys.argv=['rya759_nearuv_fe_product.py','--limit',os.environ.get('LIMIT','8'),
          '--tag', 'rya1204_mol' if USE_MOL else 'rya1204_nomol',
          '--out', os.path.join(SP,'nearuv_'+('mol' if USE_MOL else 'nomol')+'.json')]
import runpy
runpy.run_path('/Users/ryanschmitt/codex/rya1195/scripts/rya759_nearuv_fe_product.py',
               run_name='__main__')
