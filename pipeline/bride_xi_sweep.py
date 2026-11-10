"""Explicit first-axis sweep contracts; wider amplitudes need no runner rewrite."""
import copy,re
import numpy as np


def validate_plan(plan,reference):
    if plan['schema']!='codex.bride_mac_parameter_plan/1' or plan['sweep']['axis']!='microturbulence_km_s':
        raise ValueError('Unsupported sweep schema/axis')
    if plan['controls']!=reference['controls']:raise ValueError('Changed original convergence/budget controls')
    baseline=reference['cases'][0];mutable={'name','xi','bound_point_coefficients'}
    names=set();values=set()
    if not plan['cases']:raise ValueError('Empty sweep')
    for case in plan['cases']:
        if {k:v for k,v in case.items()if k not in mutable}!={k:v for k,v in baseline.items()if k not in mutable}:
            raise ValueError('This phase varies microturbulence only; changed fixed inputs')
        name=case['name'];xi=case['xi']
        if not re.fullmatch(r'[A-Za-z0-9_-]+',name)or name in names:raise ValueError('Unsafe/duplicate case name')
        if not np.isfinite(xi)or xi<0 or xi in values:raise ValueError('Invalid/duplicate microturbulence')
        if case['bound_point_coefficients']!=name:raise ValueError('Each sweep case requires its own certificate identity')
        names.add(name);values.add(xi)
    if baseline['xi']not in values:raise ValueError('Sweep requires its declared baseline')
    return plan


def source_input_xi(text,xi_km_s):
    """Replace the pinned explicit XIFIX value, refusing ambiguous source layout."""
    if not np.isfinite(xi_km_s)or xi_km_s<0:raise ValueError('Invalid producer xi')
    lines=text.splitlines();indices=[i for i,line in enumerate(lines)if line.startswith("'XIFIX:'")]
    if len(indices)!=1:raise ValueError('Missing/ambiguous XIFIX source')
    i=indices[0]
    if lines[i]!="'XIFIX:' 'T'"or i!=len(lines)-2 or float(lines[i+1])!=1.:raise ValueError('Changed pinned producer xi input')
    lines[i+1]=f'{xi_km_s:.17e}'
    return '\n'.join(lines)+'\n'


def reconstruct(case):
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure,fingerprint
    from pipeline.bride_parameter_checks import marcs_filename,interpolate_marcs_structure,controlled_structure
    from pipeline.bride_solar_control import density_column_structure
    store=require_codex_path('grids.marcs_standard');nodes=[];sources=[]
    for t in case['bracket_teff']:
        for g in case['bracket_logg']:
            path=store/marcs_filename(t,g);nodes.append(((t,g),read_marcs_structure(path)));sources.append(fingerprint(path))
    raw,interpolation=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),thermodynamic_scale=case['thermodynamic_scale'])
    base=controlled_structure(raw,xi=case['xi'],drop_bottom=case['drop_bottom'])
    structure,contract=density_column_structure(base,nodes,interpolation)
    return base,structure,contract,sources
