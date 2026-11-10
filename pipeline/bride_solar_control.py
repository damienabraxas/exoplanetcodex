"""Explicit experimental density-column and exact-input Fe provider; no adoption."""
import copy
from functools import partial
from pathlib import Path
import numpy as np

from pipeline.bride_fe_background import FeIBackground
from pipeline.bride_opacity_binding import local_binding, validate_binding, load_bound_coefficients, binding_digest
from pipeline.bride_column_checks import candidate_column, integrated_column


def density_column_structure(base, nodes, interpolation):
    if 'column_mass_kg_m2' in base:
        raise ValueError('Experimental density construction requires a column-withheld interpolant')
    indices = interpolation['simplex_source_indices']
    top = candidate_column([nodes[i][1]['column_mass_kg_m2'] for i in indices], interpolation['weights'], scale='log')[0]
    result = copy.deepcopy(base)
    result['column_mass_kg_m2'] = integrated_column(base['height_m'], base['rho_kg_m3'], top)
    return result, dict(schema='codex.bride_experimental_column/1', method='density_integral_log_native_top',
        top_anchor_kg_m2=float(top), weights=interpolation['weights'], source_indices=indices,
        source_nodes=[dict(point=list(p), top_column_kg_m2=float(s['column_mass_kg_m2'][0])) for p, s in nodes],
        interpretation='Explicit diagnostic interpolant, not native/equilibrium column. Native mass-definition basis; top anchor interpolated, no fitted scale. Zero mass closure is construction, not validation.')


def bound_fe_factory(directory, base, structure, oxygen_abundance_dex, *, column_contract=None):
    from pipeline.bride_reference import fingerprint
    data, certificate = load_bound_coefficients(directory, base, oxygen_abundance_dex)
    # A new experiment may add only its explicitly declared column. Existing
    # certificates remain immutable and cannot authorize this by flipping flags.
    stripped = copy.deepcopy(structure)
    if 'column_mass_kg_m2' in base:
        stripped['column_mass_kg_m2'] = base['column_mass_kg_m2'].copy()
    else:
        stripped.pop('column_mass_kg_m2', None)
    validate_binding(certificate['local_binding'], stripped, oxygen_abundance_dex)
    length = None
    if column_contract is None:
        validate_binding(certificate['local_binding'], structure, oxygen_abundance_dex)
    else:
        if column_contract['method'] != 'density_integral_log_native_top' or 'column_mass_kg_m2' in base:
            raise ValueError('Unknown experimental column route')
        expected = integrated_column(structure['height_m'], structure['rho_kg_m3'], column_contract['top_anchor_kg_m2'])
        if not np.array_equal(expected, structure['column_mass_kg_m2']):
            raise ValueError('Changed experimental density-integral column')
        length = float(expected[0]/structure['rho_kg_m3'][0])
    binding = local_binding(structure, oxygen_abundance_dex)
    packet = dict(wavelength_vac_nm=data['wavelength_vac_nm'], opacity_m2_kg=data['FeI_opacity_m2_kg'],
        opacity_net_per_neutral_Fe_m2=data['FeI_opacity_per_neutral_m2'],
        expected_temperature_K=np.asarray(binding['fields']['temperature_K']), expected_ne_m3=np.asarray(binding['ne_m3']),
        expected_nH_m3=np.asarray(binding['nH_m3']), expected_rho_kg_m3=np.asarray(binding['fields']['rho_kg_m3']))
    factory = partial(BoundFeIBackground, packet=packet, normalization='local_neutral', binding=binding, cap_length=length)
    return factory, dict(schema='codex.bride_experimental_bound_input/1', scientific_admission='HOLD',
        source_certificate=certificate['binding_sha256'], source_array=certificate['array'],
        source_certificate_file=fingerprint(Path(directory)/'manifest.json'),
        local_binding=binding, binding_sha256=binding_digest(binding), column_contract=column_contract,
        normalization='local_neutral', provider_scope='Explicit experimental opt-in only; legacy native-only loader unchanged')


class BoundFeIBackground(FeIBackground):
    def __init__(self, eqPops, radSet, wavelength, *, packet, normalization, binding, cap_length):
        if normalization != 'local_neutral':
            raise ValueError('Bound experimental provider permits local-neutral normalization only')
        self.bound = copy.deepcopy(binding)
        self.cap_length = cap_length
        super().__init__(eqPops, radSet, wavelength, packet=packet, normalization=normalization)
        import lightweaver as lw
        from lightweaver.atomic_table import PeriodicTable
        abundances = {PeriodicTable[i+1]: a for i, a in enumerate(binding['native_abundances_dex'])}
        abundances[PeriodicTable['O']] = binding['scalars']['active_oxygen_abundance_dex']
        expected = lw.AtomicAbundance(abundances)
        if any(radSet.abundance[PeriodicTable[i]] != expected[PeriodicTable[i]] for i in range(1, 93)):
            raise ValueError('Bound provider refuses changed actual composition/active abundance')

    def compute_background(self, atmos, chi, eta, sca):
        py = atmos.pyAtmos
        b = self.bound
        for actual, expected in ((py.height, b['fields']['height_m']), (py.vturb, [b['scalars']['vturb_m_s']]*b['ndep']),
                                 (py.vlos, b['static_vlos_m_s'])):
            if not np.array_equal(actual, expected):
                raise ValueError('Bound provider refuses changed geometry/velocity')
        if self.cap_length is not None:
            from pipeline.bride_top_cap import FiniteColumnCap
            if not isinstance(py.zUpperBc, FiniteColumnCap) or py.zUpperBc.length_m != self.cap_length:
                raise ValueError('Bound provider requires its explicit experimental cap length')
        super().compute_background(atmos, chi, eta, sca)
