from tests.bride_test_requirements import requires_lightweaver
import copy
import json
import numpy as np
import pytest


def inputs(name):
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_parameter_checks import interpolate_marcs_structure, marcs_filename
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_solar_control import density_column_structure
    certificate = require_codex_path('experiments.bride_bound_point_coefficients')/name
    store = require_codex_path('grids.marcs_standard')
    atom = read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    if name == 'native_identity':
        base = read_marcs_structure(store/marcs_filename(5750, 4.5))
        return certificate, base, copy.deepcopy(base), atom, None
    nodes = [((t, g), read_marcs_structure(store/marcs_filename(t, g))) for t in (5750, 6000) for g in (4., 4.5)]
    base, receipt = interpolate_marcs_structure(nodes, (5772, 4.438), thermodynamic_scale='log_positive_layers')
    structure, contract = density_column_structure(base, nodes, receipt)
    return certificate, base, structure, atom, contract


@requires_lightweaver
def test_actual_native_bound_provider_identity_and_deepcopy():
    from config.constants import require_codex_path
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_fe_background import fe_factory
    from pipeline.bride_solar_control import bound_fe_factory
    folder, base, s, atom, contract = inputs('native_identity')
    provider, record = bound_fe_factory(folder, base, s, atom.abundance)
    old = make_reference_context(s, atom, rays=2, qscale=1,
        fe_background=fe_factory(require_codex_path('experiments.bride_fe_packet'), 'local_neutral'))
    new = make_reference_context(s, atom, rays=2, qscale=1, fe_background=provider)
    for field in ('chi', 'eta', 'sca'):
        np.testing.assert_array_equal(getattr(old.background, field), getattr(new.background, field))
    cloned = copy.deepcopy(new)
    np.testing.assert_array_equal(cloned.background.chi, new.background.chi)


@requires_lightweaver
def test_actual_solar_explicit_cap_and_runtime_mismatch_refusal():
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import FiniteColumnCap
    from pipeline.bride_solar_control import bound_fe_factory
    folder, base, s, atom, contract = inputs('fine_T5772')
    provider, record = bound_fe_factory(folder, base, s, atom.abundance, column_contract=contract)
    length = s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]
    ctx = make_reference_context(s, atom, rays=2, qscale=1, fe_background=provider,
        experimental_column_cap=FiniteColumnCap(length))
    assert ctx.atmos.pyAtmos.zUpperBc.length_m == pytest.approx(97762.93901851328)
    for field in ('height', 'vturb', 'vlos', 'temperature', 'ne', 'nHTot'):
        values = getattr(ctx.atmos.pyAtmos, field); old = values.copy()
        values[0] = np.nextafter(values[0], np.inf)
        with pytest.raises(ValueError):
            ctx.background.update_background(ctx.atmos)
        values[:] = old
    ctx.background.update_background(ctx.atmos)
    ctx.atmos.pyAtmos.zUpperBc.length_m *= 1.001
    with pytest.raises(ValueError, match='cap length'):
        ctx.background.update_background(ctx.atmos)
    ctx.atmos.pyAtmos.zUpperBc.length_m = length
    with pytest.raises(ValueError, match='explicit experimental cap'):
        make_reference_context(s, atom, rays=2, qscale=1, fe_background=provider)
    changed = copy.deepcopy(s); changed['column_mass_kg_m2'][-1] *= 1.001
    with pytest.raises(ValueError, match='density-integral'):
        bound_fe_factory(folder, base, changed, atom.abundance, column_contract=contract)


@requires_lightweaver
def test_actual_composition_and_native_interpolated_cap_guards():
    from dataclasses import replace
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import cap_from_structure
    from pipeline.bride_solar_control import bound_fe_factory
    folder, base, s, atom, contract = inputs('fine_T5772')
    provider, record = bound_fe_factory(folder, base, s, atom.abundance, column_contract=contract)
    with pytest.raises(ValueError, match='composition/active'):
        make_reference_context(s, replace(atom, abundance=atom.abundance+.01), rays=2, qscale=1, fe_background=provider)
    with pytest.raises(ValueError, match='Native column mass'):
        cap_from_structure(base)
