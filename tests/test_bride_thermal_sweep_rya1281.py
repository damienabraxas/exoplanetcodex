import json,copy
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]


def plans():return [json.loads((ROOT/'config/experiments'/name).read_text())for name in ('rya1281_mac_thermal_slice.json','rya1281_mac_solar_bound_control.json')]


def test_plan_axes_brackets_and_fixed_physics():
    from pipeline.bride_thermal_sweep import validate_plan,bracket
    p,r=plans();validate_plan(p,r)
    for key,value in (('xi',1.1),('abundance',8.8),('drop_bottom',1)):
        bad=copy.deepcopy(p);bad['cases'][0][key]=value
        with pytest.raises(ValueError,match='fixed'):validate_plan(bad,r)
    bad=copy.deepcopy(p);bad['cases'][0]['logg']=4.338
    with pytest.raises(ValueError,match='Teff only'):validate_plan(bad,r)
    with pytest.raises(ValueError):bracket(6100,[5500,5750,6000])


def test_departure_warm_start_uses_new_nuclei_total_and_preserves_nominal():
    from pipeline.bride_thermal_sweep import warm_start
    pop=np.array([[2.],[2.]]);old=np.array([[1.],[3.]]);new=np.array([[1.],[5.]])
    np.testing.assert_allclose(warm_start(pop,old,new,np.array([6.])),[[2.25],[3.75]],rtol=1e-15)
    np.testing.assert_array_equal(warm_start(pop,old,old,np.array([4.])),pop)
    with pytest.raises(ValueError):warm_start(pop,old,new,np.array([0.]))


def test_covered_reconstructions_and_hotter_refusal():
    from pipeline.bride_thermal_sweep import reconstruct
    p,_=plans()
    for case in p['cases']:
        base,s,contract,receipts=reconstruct(case)
        assert s['ndep']==56 and max(s['temperature_K'])<=10000
    hot=copy.deepcopy(p['cases'][3]);hot['teff']=5872
    with pytest.raises(ValueError,match='collision tables'):reconstruct(hot)


def test_refinement_changes_sampling_only_and_clipping_count_not_source():
    from pipeline.bride_thermal_sweep import effective_controls,validate_spline_receipt
    p,_=plans();r=effective_controls(p['controls'],True)
    assert {k for k in r if r[k]!=p['controls'][k]}=={'rays','qscale','flux_rays'}
    from types import SimpleNamespace
    atom=SimpleNamespace(collisions=[SimpleNamespace(kind='CH')])
    ref=dict(source_directions=True,CH_negative_spline_values_clipped_as_native=0)
    actual=dict(ref,CH_negative_spline_values_clipped_as_native=2)
    validate_spline_receipt(actual,ref,atom,56)
    actual['source_directions']=False
    with pytest.raises(ValueError):validate_spline_receipt(actual,ref,atom,56)


def test_analysis_rejects_changed_fields_and_diagnostic_column_even_with_new_digest():
    from config.constants import require_codex_path
    from scripts.analyze_bride_thermal_sweep import verify_case_binding
    from pipeline.bride_opacity_binding import binding_digest
    r=json.loads((require_codex_path('experiments.bride_mac_grid')/'rya1281-thermal-coupled-v1/thermal_nominal/report.json').read_text())
    c=json.loads((require_codex_path('experiments.bride_thermal_coefficients')/'thermal_nominal/manifest.json').read_text())
    verify_case_binding(r,c)
    bad=copy.deepcopy(r);b=bad['local_binding']['local_binding'];b['fields']['temperature_K'][0]+=1.;bad['local_binding']['binding_sha256']=binding_digest(b)
    with pytest.raises(ValueError,match='fields differ'):verify_case_binding(bad,c)
    bad=copy.deepcopy(r);b=bad['local_binding']['local_binding'];b['column_mass_kg_m2'][0]*=1.01;bad['local_binding']['binding_sha256']=binding_digest(b)
    with pytest.raises(ValueError,match='constructed thermal column'):verify_case_binding(bad,c)
