import copy
import numpy as np
import pytest


def contexts():
    from tests.test_bride_solar_control_rya1281 import inputs
    from tests.test_bride_trace_fe_lines_rya1281 import line
    from tests.test_bride_fe_pressure_rya1281 import recipe
    from tests.test_bride_fe_stark_rya1281 import row
    from config.constants import require_codex_path
    from pipeline.bride_h2_tables import H2Tables
    from pipeline.bride_h2_colliders import h2_collider_factory
    from pipeline.bride_fe_stark import stark_fe_factory
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import FiniteColumnCap
    folder,base,s,atom,contract=inputs('fine_T5772')
    provider,_=bound_fe_factory(folder,base,s,atom.abundance,column_contract=contract)
    args=(provider,[line()],(777.,779.),[recipe()],[row()],'cowley')
    t=H2Tables(require_codex_path('experiments.bride_molecular_tables'))
    def ctx(factory):return make_reference_context(s,atom,rays=2,qscale=1,fe_background=factory,
        extra_wavelengths_vac_nm=np.array([777.5,778.,778.5]),experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))
    return ctx(stark_fe_factory(*args)),ctx(h2_collider_factory(*args,t,'off')),ctx(h2_collider_factory(*args,t,'molecular'))


def test_actual_h2_off_identity_and_partial_intervention():
    from lightweaver.atomic_table import PeriodicTable
    old,off,on=contexts()
    for field in ('chi','eta','sca'):np.testing.assert_array_equal(getattr(old.background,field),getattr(off.background,field))
    np.testing.assert_array_equal(on.background.sca,off.background.sca)
    for element in ('H','He','O','Fe'):
        np.testing.assert_array_equal(on.kwargs['eqPops'].atomicPops[PeriodicTable[element]].n,off.kwargs['eqPops'].atomicPops[PeriodicTable[element]].n)
    np.testing.assert_array_equal(on.atmos.pyAtmos.ne,off.atmos.pyAtmos.ne)
    assert np.any(np.asarray(on.background.chi)!=np.asarray(off.background.chi))
    clone=copy.deepcopy(on);np.testing.assert_array_equal(clone.background.chi,on.background.chi)


def test_h2_cache_refuses_changed_hminus_even_if_atomic_populations_unchanged():
    _,_,on=contexts()
    minus=on.kwargs['eqPops'].HminPops;minus[0]=np.nextafter(minus[0],np.inf)
    with pytest.raises(ValueError,match='hydrogen inputs'):on.background.update_background(on.atmos)
