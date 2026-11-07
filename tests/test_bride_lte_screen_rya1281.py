import numpy as np
import pytest


def test_independent_saha_matches_native_and_conserves_species():
    from types import SimpleNamespace
    from lightweaver.atomic_table import KuruczPfTable,PeriodicTable
    from pipeline.bride_lte_screen import independent_ion_populations
    a=SimpleNamespace(temperature=np.array([4100.,6000.,9900.]),ne=np.array([1e16,1e19,1e21]),nHTot=np.array([1e19,1e22,1e24]),Nspace=3)
    for symbol in ('Fe','Ti','Si','Ni','Ca'):
        p=KuruczPfTable()[PeriodicTable[symbol]]
        actual=independent_ion_populations(p,a);native=p.lte_ionisation(a)
        np.testing.assert_allclose(actual,native,rtol=1e-12,atol=0)
        np.testing.assert_allclose(actual.sum(axis=0),p.abundance*a.nHTot,rtol=1e-14)


def test_source_refusals_do_not_choose_new_gf_or_components():
    from copy import deepcopy
    from pipeline.bride_lte_screen import source_eligibility
    item=dict(treatment='not_in_trace_subset',exact_identity_matches=1,canonical_gf_matches_any_source=True,
              canonical_row=dict(hfs_n_components='1'),source_rows=[dict(turbospectrum_rad='1e8',rad='8')])
    assert source_eligibility(item)=='SOURCE_ELIGIBLE_FOR_NATURAL_ONLY_SCREEN'
    assert source_eligibility(item,changed_on_main=True)=='MAIN_CATALOG_CHANGED'
    for field,value,reason in [('exact_identity_matches',2,'SOURCE_IDENTITY_NOT_UNIQUE'),('canonical_gf_matches_any_source',False,'CANONICAL_SOURCE_GF_DIFFER')]:
        bad=deepcopy(item);bad[field]=value;assert source_eligibility(bad)==reason
    item['source_rows'][0]['turbospectrum_rad']='0'
    assert source_eligibility(item)=='NATURAL_WIDTH_UNQUALIFIED'


def test_partition_domain_is_not_extrapolated():
    from types import SimpleNamespace
    from lightweaver.atomic_table import KuruczPfTable,PeriodicTable
    from pipeline.bride_lte_screen import independent_ion_populations
    p=KuruczPfTable()[PeriodicTable['Fe']]
    a=SimpleNamespace(temperature=np.array([p.Tpf[0]-1]),ne=np.ones(1),nHTot=np.ones(1),Nspace=1)
    with pytest.raises(ValueError,match='domain'):independent_ion_populations(p,a)
