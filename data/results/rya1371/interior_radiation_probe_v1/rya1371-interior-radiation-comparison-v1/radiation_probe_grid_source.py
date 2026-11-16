"""Bounded source-indicator interior radiation probe; no node removal or tiling."""
import json
from pathlib import Path
import numpy as np
from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS
from pipeline.saved_field_quadrature import gaussian_panel_batches

METHOD='vacuum-nm-radiation-probe-gauss4'


def registered_probe_destination(output,resource,*,resolver):
    expected=Path(resolver(resource)).resolve();output=Path(output).resolve()
    if output!=expected or output.exists():raise ValueError('Require fresh registered probe destination')
    return output


def gaussian_probe_grid(inherited,scores,budget=32):
    wave=np.asarray(inherited,float);scores=np.asarray(scores,float)
    if wave.ndim!=1 or len(wave)<2 or not np.isfinite(wave).all()or np.any(wave<=0)or np.any(np.diff(wave)<=0)or scores.shape!=(len(wave)-1,)or not np.isfinite(scores).all()or np.any(scores<0):raise ValueError('Invalid radiation probe grid/indicators')
    if isinstance(budget,bool)or not isinstance(budget,int)or budget<=0 or budget%4 or len(wave)+budget>MAX_PACKET_WAVELENGTHS:raise ValueError('Probe exceeds original20000point or Gaussian budget')
    positive=np.flatnonzero(scores>0)
    selected=positive[np.argsort(-scores[positive],kind='stable')[:budget//4]]
    if len(selected)!=budget//4:raise ValueError('Insufficient positive source-indicator panels')
    points=[];weights=[]
    for panel in selected:
        query,measure=next(gaussian_panel_batches(wave[panel:panel+2],4))
        points.extend(query);weights.extend(measure)
    points=np.asarray(points);weights=np.asarray(weights);expanded=np.sort(np.r_[wave,points])
    if np.any(np.diff(expanded)<=0)or len(expanded)!=len(wave)+budget:raise ValueError('Probe nodes duplicated/unrepresentable')
    np.testing.assert_array_equal(expanded[np.searchsorted(expanded,wave)],wave)
    return expanded,points,weights,selected


def physics_plan(plan):
    return {key:value for key,value in plan.items()if key not in ['scope','wavelength_grid_control']}


def load_probe_control(spec,plan,native,state,atom,pop,source_star,*,resolver):
    from pipeline.reference_engines import fingerprint
    from pipeline.joint_oxygen_lte_control import validate_joint_plan,held_lw_reference_star
    validate_joint_plan(plan)
    if set(spec)!={'resource','packet','producer','method'}or spec['method']!=METHOD:raise ValueError('Unknown explicit radiation probe control')
    directory=Path(resolver(spec['resource']));packet=directory/'probe_grid.npz';producer=directory/'report.json'
    if fingerprint(packet)!=spec['packet']or fingerprint(producer)!=spec['producer']:raise ValueError('Changed radiation probe packet/producer')
    doc=json.loads(producer.read_text())
    if not doc['computational_completed']or doc['method']!=METHOD or doc['budget']!=32 or doc['physics_plan']!=physics_plan(plan)or doc['output_resource']!=spec['resource']:raise ValueError('Wrong radiation probe physics/producer')
    checks={str(packet):fingerprint(packet),str(producer):fingerprint(producer)}
    for path,expected in doc['inputs'].items():
        if fingerprint(path)!=expected:raise ValueError('Changed source indicator provenance '+path)
        checks[path]=expected
    for name in ['arrays.npz','report.json']:
        if doc['inputs'].get(str(native/name))!=fingerprint(native/name):raise ValueError('Wrong probe native-state owner')
    with np.load(packet)as d:
        wave=d['wave_vac_nm'].copy();base=d['inherited_wave_vac_nm'].copy();scores=d['panel_scores'].copy();probe=d['probe_wave_vac_nm'].copy();weights=d['probe_weights_nm'].copy();panels=d['selected_panels'].copy()
        np.testing.assert_array_equal(d['populations'],pop);np.testing.assert_array_equal(d['native_nstar'],source_star)
        np.testing.assert_array_equal(d['declared_nstar'],held_lw_reference_star(atom,state,source_star))
    fine_path=Path(doc['fine_rates'])
    if str(fine_path)not in checks:raise ValueError('Missing parent radiation receipt')
    with np.load(fine_path)as d:
        np.testing.assert_array_equal(base,d['wave_vac_nm']);np.testing.assert_array_equal(pop,d['populations'])
    expected,query,measure,selected=gaussian_probe_grid(base,scores,32)
    for actual,wanted in [(wave,expected),(probe,query),(weights,measure),(panels,selected)]:np.testing.assert_array_equal(actual,wanted)
    checks[str(Path(__file__))]=fingerprint(__file__)
    return wave,dict(method=METHOD,coarse_points=len(base),refined_points=len(wave),probe_points=len(probe),selected_panels=len(panels),scope='All inherited nodes retained;32Gauss4interior samples selected by all-BF source-normalized represented-field indicators, not native residual. No physical/quadrature adoption.'),checks
