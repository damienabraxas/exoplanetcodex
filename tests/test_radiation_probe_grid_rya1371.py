import numpy as np
import pytest
from pipeline.radiation_probe_grid import gaussian_probe_grid


def test_probe_preserves_every_node_and_ranks_only_positive_indicators():
    wave=np.array([1.,1.5,2.,2.5,3.]);scores=np.array([0.,9.,2.,1.])
    expanded,points,weights,panels=gaussian_probe_grid(wave,scores,8)
    np.testing.assert_array_equal(panels,[1,2]);np.testing.assert_array_equal(expanded[np.searchsorted(expanded,wave)],wave)
    assert len(expanded)==len(wave)+8 and np.all(weights>0)
    assert np.all((points[:4]>1.5)&(points[:4]<2.))and np.all((points[4:]>2.)&(points[4:]<2.5))
    np.testing.assert_allclose(weights.sum(),1.,rtol=2e-16)
    # Shared score ties resolve by wavelength order, without a reference target.
    np.testing.assert_array_equal(gaussian_probe_grid(wave,np.ones(4),4)[3],[0])


def test_probe_refuses_cap_partial_groups_missing_indicators_and_bad_grid():
    with pytest.raises(ValueError,match='20000'):gaussian_probe_grid(np.arange(1.,19971.),np.ones(19969),32)
    for scores,budget in [([1.,1.],3),([0.,0.],4),([1.,np.nan],4),([-1.,1.],4)]:
        with pytest.raises(ValueError):gaussian_probe_grid([1.,2.,3.],scores,budget)
    with pytest.raises(ValueError):gaussian_probe_grid([1.,1.,2.],[1.,1.],4)


def test_loader_binds_parent_physics_and_rejects_removed_node_even_with_new_packet_hash(tmp_path,monkeypatch):
    import json
    from pipeline.reference_engines import fingerprint
    from pipeline.radiation_probe_grid import load_probe_control,physics_plan,METHOD
    monkeypatch.setattr('pipeline.joint_oxygen_lte_control.held_lw_reference_star',lambda atom,state,star:star.copy())
    plan=dict(joint_lte_control=dict(method='held-LW-finite-level-source-forward'),emission_profile_control=dict(profile='held-LW'),boundary_control='native-zero-ghost-linear',native_physical_control={},continuum_interpolation_control={})
    wave=np.arange(1.,11.);scores=np.ones(9);expanded,probe,weights,selected=gaussian_probe_grid(wave,scores)
    pop=np.ones((41,2));native=tmp_path/'native';native.mkdir();(native/'arrays.npz').write_bytes(b'synthetic receipt only');(native/'report.json').write_text('{}')
    parent=tmp_path/'parent.npz';np.savez(parent,wave_vac_nm=wave,populations=pop)
    packet=tmp_path/'probe_grid.npz';arrays=dict(wave_vac_nm=expanded,inherited_wave_vac_nm=wave,panel_scores=scores,probe_wave_vac_nm=probe,probe_weights_nm=weights,selected_panels=selected,populations=pop,native_nstar=pop,declared_nstar=pop);np.savez(packet,**arrays)
    report=tmp_path/'report.json';doc=dict(computational_completed=True,method=METHOD,budget=32,physics_plan=physics_plan(plan),output_resource='synthetic.probe',inputs={str(p):fingerprint(p)for p in [parent,native/'arrays.npz',native/'report.json']},fine_rates=str(parent));report.write_text(json.dumps(doc))
    spec=dict(resource='synthetic.probe',packet=fingerprint(packet),producer=fingerprint(report),method=METHOD)
    result,_,_=load_probe_control(spec,plan,native,{},None,pop,pop,resolver=lambda key:tmp_path);np.testing.assert_array_equal(result,expanded)
    with pytest.raises(ValueError,match='physics'):load_probe_control(spec,dict(plan,seconds=1),native,{},None,pop,pop,resolver=lambda key:tmp_path)
    # Refreshing a checksum cannot make a removed inherited node admissible.
    removed=np.flatnonzero(expanded==wave[3])[0];arrays['wave_vac_nm']=np.delete(expanded,removed);np.savez(packet,**arrays);bad=dict(spec,packet=fingerprint(packet))
    with pytest.raises(AssertionError):load_probe_control(bad,plan,native,{},None,pop,pop,resolver=lambda key:tmp_path)
    parent.write_bytes(b'changed producer')
    with pytest.raises(ValueError,match='provenance'):load_probe_control(bad,plan,native,{},None,pop,pop,resolver=lambda key:tmp_path)


def test_registered_output_may_be_absent_but_cannot_overwrite_or_change_location(tmp_path):
    from pipeline.radiation_probe_grid import registered_probe_destination
    destination=tmp_path/'fresh'
    assert registered_probe_destination(destination,'probe',resolver=lambda key:destination)==destination
    assert not destination.exists()
    with pytest.raises(ValueError):registered_probe_destination(tmp_path/'other','probe',resolver=lambda key:destination)
    destination.mkdir();protected=destination/'retained.txt';protected.write_text('preserve')
    with pytest.raises(ValueError):registered_probe_destination(destination,'probe',resolver=lambda key:destination)
    assert protected.read_text()=='preserve'
