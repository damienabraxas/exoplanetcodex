"""Explicit source-matched solar Teff/logg experiments; no physics admission."""
import copy,re
import numpy as np
from pipeline.bride_xi_sweep import reconstruct as reconstruct_base


def bracket(value,nodes):
    nodes=np.asarray(nodes,dtype=float)
    if nodes.ndim!=1 or len(nodes)<2 or not np.isfinite(nodes).all()or np.any(np.diff(nodes)<=0)or not nodes[0]<=value<=nodes[-1]:raise ValueError('Outside declared source grid')
    i=min(max(int(np.searchsorted(nodes,value,side='right'))-1,0),len(nodes)-2)
    return nodes[i:i+2].tolist()


def validate_plan(plan,reference):
    if plan['schema']!='codex.bride_mac_parameter_plan/1' or plan['sweep']['axis']!='solar_Teff_logg':raise ValueError('Unsupported thermal sweep')
    if plan['controls']!=reference['controls']:raise ValueError('Changed original convergence/budget controls')
    base=reference['cases'][0];mutable={'name','teff','logg','bracket_teff','bracket_logg','bound_point_coefficients'}
    if not plan['cases']:raise ValueError('Empty thermal sweep')
    names=set();points=set();nominal=0
    for c in plan['cases']:
        if {k:v for k,v in c.items()if k not in mutable}!={k:v for k,v in base.items()if k not in mutable}:raise ValueError('Changed fixed thermal-sweep inputs')
        if not re.fullmatch(r'[A-Za-z0-9_-]+',c['name'])or c['name']in names or c['bound_point_coefficients']!=c['name']:raise ValueError('Invalid case/certificate identity')
        t,g=c['teff'],c['logg'];axis=plan['sweep']['case_axes'][c['name']]
        if not np.isfinite([t,g]).all()or (t,g)in points:raise ValueError('Invalid/duplicate thermal point')
        if axis=='temperature'and(t==base['teff']or g!=base['logg']):raise ValueError('Temperature case must vary Teff only')
        elif axis=='gravity'and(g==base['logg']or t!=base['teff']):raise ValueError('Gravity case must vary logg only')
        elif axis=='nominal':
            if (t,g)!=(base['teff'],base['logg']):raise ValueError('Changed nominal thermal point')
            nominal+=1
        elif axis not in ('temperature','gravity'):raise ValueError('Unknown thermal axis')
        if c['bracket_teff']!=bracket(t,plan['source_grid']['teff_K'])or c['bracket_logg']!=bracket(g,plan['source_grid']['logg_cgs']):raise ValueError('Changed nearest-cell source brackets')
        names.add(c['name']);points.add((t,g))
    if nominal!=1 or set(plan['sweep']['case_axes'])!=names:raise ValueError('Require one nominal and exact axis accounting')
    return plan


def reconstruct(case):
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_parameter_checks import assert_collision_coverage
    values=reconstruct_base(case)
    atom=read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    assert_collision_coverage(atom,values[1]['temperature_K'])
    return values


def warm_start(reference_pop,reference_nstar,target_nstar,target_total):
    """Preserve source departures, then enforce the target atmosphere's O total."""
    pop,old,new=[np.asarray(v)for v in (reference_pop,reference_nstar,target_nstar)];total=np.asarray(target_total)
    if pop.ndim!=2 or pop.shape!=old.shape or old.shape!=new.shape or total.shape!=(pop.shape[1],):raise ValueError('Warm-start shape mismatch')
    if not all(np.isfinite(v).all()and np.all(v>0)for v in (pop,old,new,total)):raise ValueError('Invalid warm-start populations')
    if np.array_equal(old,new)and np.allclose(pop.sum(axis=0),total,rtol=1e-10,atol=0):return pop.copy()
    with np.errstate(over='raise',invalid='raise',divide='raise'):
        result=(pop/old)*new;result*=total/result.sum(axis=0)
    if not np.isfinite(result).all()or np.any(result<=0):raise ValueError('Unrepresentable warm start')
    return result


def validate_spline_receipt(actual,reference,atom,ndep):
    """Only native CH clipping count can vary with the new local temperature."""
    a=copy.deepcopy(actual);r=copy.deepcopy(reference);count=a.pop('CH_negative_spline_values_clipped_as_native');r.pop('CH_negative_spline_values_clipped_as_native')
    if a!=r or not isinstance(count,int)or isinstance(count,bool)or not 0<=count<=sum(c.kind=='CH'for c in atom.collisions)*ndep:raise ValueError('Changed collision source/prescription')


def effective_controls(original,refine=False):
    result=copy.deepcopy(original)
    if refine:result.update(rays=5,qscale=2,flux_rays=7)
    return result
