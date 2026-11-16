"""Explicit fixed-population LTE/inverse diagnostic, never coupled admission."""
from pathlib import Path
import json
import numpy as np
from pipeline.native_oxygen_collisions import collision_matrix


def held_lw_reference_star(atom,state,original_star):
    from lightweaver import atomic_set,constants as c
    from pipeline.multi_atom import to_lightweaver
    from pipeline.oxygen_lte_contract import held_lw_weights,normalized_weights
    if atom.symbol!='O'or len(atom.levels)!=41:raise ValueError('Joint LTE diagnostic requires complete O41')
    model=to_lightweaver(atom);total=np.asarray(original_star).sum(axis=0)
    candidate=atomic_set.lte_pops(model,state['temperature_K'],state['ne_m3'],total,debye=True)
    constants=dict(h=c.HPlanck,k=c.KBoltzmann,m=c.MElectron,q=c.QElectron,epsilon=c.Epsilon0)
    independent=normalized_weights(held_lw_weights([l.E_SI for l in model.levels],[l.g for l in model.levels],np.array([l.stage+1 for l in model.levels]),state['temperature_K'],state['ne_m3'],constants),total)
    np.testing.assert_allclose(candidate,independent,rtol=1e-12,atol=0)
    np.testing.assert_allclose(candidate.sum(axis=0),total,rtol=1e-12,atol=0)
    return candidate


def source_direction_inverse_control(records,q,T,ne,H,HI,original_star,candidate_star,E,g,e_k):
    """Rebuild each source channel with its ORIGINAL forward direction.

    CE/CH source forward is downward; CI/CH0 is upward. The earlier all-upward
    algebraic canary is deliberately not reused. Returns aggregate C plus exact
    per-kind source-forward identity checks. No source coefficient is fitted.
    """
    original=np.asarray(original_star);candidate=np.asarray(candidate_star)
    if original.shape!=candidate.shape or not np.allclose(original.sum(axis=0),candidate.sum(axis=0),rtol=1e-12,atol=0):raise ValueError('LTE inverse control changed total/shape')
    args=(T,ne,H,HI);tail=(E,g,e_k)
    matrix=collision_matrix(records,q,*args,candidate,*tail)
    summary={};lower,upper=np.triu_indices(41,1)
    for kind in ['CE','CI','CH','CH0']:
        indices=[i for i,r in enumerate(records)if r[1]==kind]
        if not indices:raise ValueError('Missing source collision channel '+kind)
        rows=[]
        for new_index,index in enumerate(indices):
            row=list(records[index]);row[0]=str(new_index);rows.append(row)
        before=collision_matrix(rows,q[indices],*args,original,*tail)
        after=collision_matrix(rows,q[indices],*args,candidate,*tail)
        target,source=(lower,upper)if kind in ['CE','CH']else(upper,lower)
        np.testing.assert_array_equal(after[target,source],before[target,source])
        summary[kind]=dict(records=len(rows),preserved_direction='downward'if kind in ['CE','CH']else'upward',source_forward_exact=True,source_positive_pair_depths=int(np.count_nonzero(before[target,source])))
    return matrix,summary


def rebind_source_inverse(binding,atom,state,candidate_star):
    """Use the already verified native binding's exact source/collider receipts."""
    from pipeline.bride_collision_binding import FixedNativeCollisionBinding
    from pipeline.reference_engines import fingerprint
    if not isinstance(binding,FixedNativeCollisionBinding)or not binding.receipt:raise ValueError('Require verified native collision binding')
    for path,expected in binding.receipt.items():
        if fingerprint(path)!=expected:raise ValueError('Changed collision source before inverse control')
    table=Path(next(p for p in binding.receipt if p.endswith('/native/atoms/O.txt')))
    packet=Path(next(p for p in binding.receipt if p.endswith('/packet.npz')))
    matrix_file=Path(next(p for p in binding.receipt if p.endswith('/matrices.npz')))
    report_path=matrix_file.parent/'report.json';report=json.loads(report_path.read_text())
    text=table.read_text().splitlines();header=next(i for i,x in enumerate(text)if x.split()[:3]==['key','i','j']);records=[x.split()for x in text[header+1:]]
    level_header=next(i for i,x in enumerate(text)if x.split()==['ev','g','ion','label']);levels=[list(map(float,x.split()[:4]))for x in text[level_header+1:level_header+42]]
    with np.load(packet)as d:
        for key in ['temperature_K','ne_m3','rho_kg_m3']:np.testing.assert_array_equal(d[key],state[key])
        H=d['H_ground_m3'].copy()*1e-6;HI=d['H_I_m3'].copy()*1e-6
    with np.load(matrix_file)as d:q=d['coefficients_cm3'].copy()
    matrix,summary=source_direction_inverse_control(records,q,state['temperature_K'],state['ne_m3']*1e-6,H,HI,binding.star*1e-6,candidate_star*1e-6,np.array([x[1]for x in levels]),np.array([x[2]for x in levels]),report['native_e_k'])
    return FixedNativeCollisionBinding(matrix,state,candidate_star,binding.fields['nHTot'],binding.receipt),dict(method='held-LW-finite-level-source-forward',channels=summary,max_relative_nstar_change=float(np.max(abs(candidate_star/binding.star-1))),scope='Prescribed native populations/total/colliders/source-forward channels retained; finite-level held-LW LTE and source-direction inverses only. Fixed-population diagnostic, not absolute EOS or coupled admission.')


def validate_joint_plan(plan):
    control=plan.get('joint_lte_control')
    if control!={'method':'held-LW-finite-level-source-forward'}or plan.get('emission_profile_control')!={'profile':'held-LW'}or 'fixed_population_thermal_control'in plan or plan.get('boundary_control')!='native-zero-ghost-linear'or 'native_physical_control'not in plan or 'continuum_interpolation_control'not in plan:
        raise ValueError('Require explicit matched-state held-LW joint LTE/alpha/emission diagnostic')
