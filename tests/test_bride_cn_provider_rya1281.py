import pytest
pytest.importorskip('lightweaver')


def test_actual_lw_synthetic_provider_array_invariants():
    from scripts.check_bride_cn_provider import checks
    result=checks()
    assert result['disabled_identity']and result['zero_population_identity']and result['scattering_identity']
    assert 'CN'not in result['base_lookup_inventory']
    assert len(result['mutation_refusals'])==9


def test_inventory_mapping_and_non_synthetic_packets_refused():
    from dataclasses import replace
    from scripts.check_bride_cn_provider import synthetic_problem
    from pipeline.bride_cn_provider import synthetic_cn_factory
    _,packet,_=synthetic_problem()
    kwargs=dict(enabled=True,cn_already_present=False,support_vac_nm=(777.,779.),mapping='exact_vacuum_frequency')
    for override in [dict(cn_already_present=True),dict(cn_already_present=None),dict(mapping='native_linear_air'),dict(support_vac_nm=(779.,780.))]:
        with pytest.raises(ValueError):synthetic_cn_factory(packet,**dict(kwargs,**override))
    with pytest.raises(ValueError):replace(packet,purpose='native_external_diagnostic')
