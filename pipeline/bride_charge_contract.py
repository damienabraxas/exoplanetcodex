"""Pinned O+H source comparison and explicitly experimental LS-to-J rate pair.

The data identify a coefficient relationship, not the original MULTI producer.
Statistical fine-structure branching and an LTE H-minus reservoir are assumptions.
"""
import json
from pathlib import Path

import numpy as np


# O41 one-based levels, author one-based LS state, sum of LS statistical weights.
TERM_MAP = (
    ((6,), 4, 5), ((7,), 5, 3), ((8, 9, 10), 7, 15),
    ((11, 12, 13), 8, 9), ((14,), 9, 5), ((15,), 10, 3),
    ((16, 17, 18, 19, 20), 11, 25), ((21, 22, 23), 12, 15),
    ((24, 25, 26), 13, 15), ((27, 28, 29), 14, 9),
)
SOURCE_COMMIT = 'cd2f54ed46687659017b343763e95478ed93653e'
RECEIPT_SHA256 = 'bcf3d889a969309a9c083aceaee4b129b1892de73cd9375db1b56b1969f1bdc0'


def read_author_rates(directory):
    """Verify acquired Git/SHA256 receipts before reading final/initial matrices."""
    import hashlib
    from pipeline.bride_reference import fingerprint
    directory = Path(directory)
    if fingerprint(directory / 'receipt.json')['sha256'] != RECEIPT_SHA256:
        raise ValueError('Changed author acquisition receipt')
    receipt = json.loads((directory / 'receipt.json').read_text())
    if receipt['commit'] != SOURCE_COMMIT:
        raise ValueError('Wrong author data revision')
    expected = {'states.dat', 'ReadMe.simple', 'ReadMe.cds'} | {
        f'{t}_K.rates' for t in range(1000, 10001, 1000)}
    if {entry['name'] for entry in receipt['files']} != expected:
        raise ValueError('Incomplete author source receipt')
    for entry in receipt['files']:
        data = (directory / entry['name']).read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if (blob != entry['git_blob'] or hashlib.sha256(data).hexdigest() != entry['sha256']
                or len(data) != entry['bytes']):
            raise ValueError('Changed author source: ' + entry['name'])
    tables = np.array([np.loadtxt(directory / f'{t}_K.rates')
                       for t in range(1000, 10001, 1000)])
    if tables.shape != (10, 21, 21) or not np.isfinite(tables).all() or np.any(tables < 0):
        raise ValueError('Invalid author rate matrices')
    return np.arange(1000., 10001., 1000.), tables, fingerprint(directory / 'receipt.json')


def install_author_mn_pair(ctx, source, directory, oracle):
    """Replace CH0 only with sourced MN, statistical J branching and LTE reverse.

    MN coefficients are cm^3/s in the author tables; backend rates are s^-1.
    The inverse rate uses the same local LTE reservoirs as the prior closure test.
    No arbitrary scale, abundance adjustment or outside-domain interpolation.
    """
    from lightweaver.collisional_rates import CollisionalRates
    import lightweaver.constants as c
    from pipeline.bride_reference import fingerprint
    if len(ctx.activeAtoms) != 1:
        raise ValueError('One active oxygen atom required')
    hydrogen = ctx.kwargs['eqPops'].atomicPops['H']
    if not np.array_equal(hydrogen.n, hydrogen.nStar):
        raise ValueError('Author MN pair requires LTE hydrogen')
    model = ctx.activeAtoms[0].atomicModel
    grid, tables, receipt = read_author_rates(directory)
    temperature = ctx.atmos.pyAtmos.temperature.copy()
    if len(model.collisions) != len(source.collisions):
        raise ValueError('Collision layout changed')
    mapping = {j: (state, weight) for levels, state, weight in TERM_MAP for j in levels}
    if {r.j for r in source.collisions if r.kind == 'CH0'} != set(mapping):
        raise ValueError('Unexpected CH0 level set')

    class QuantumPair(CollisionalRates):
        source_kind = 'CH0'

        def compute_rates(self, atmos, eqpops, matrix):
            if not np.array_equal(atmos.temperature, temperature):
                raise ValueError('Author control atmosphere changed')
            h = eqpops.atomicPops['H']
            if not np.array_equal(h.n, h.nStar):
                raise ValueError('Author MN pair requires LTE hydrogen')
            # Ground-H Saha attachment, separately labelled from LW HminPops.
            ci = c.HPlanck**2 / (2*np.pi*c.MElectron*c.KBoltzmann)
            hminus = atmos.ne*h.nStar[0]*.25*(ci/temperature)**1.5 * np.exp(
                c.E_ION_HMIN/(c.KBoltzmann*temperature))
            down = self.coefficient*hminus
            lte = eqpops.atomicPops[model.element].nStar
            up = down*lte[self.j]/lte[self.i]
            matrix[self.i, self.j] += down
            matrix[self.j, self.i] += up

    count = 0
    for index, raw in enumerate(source.collisions):
        if raw.kind != 'CH0':
            continue
        if raw.i != 40:
            raise ValueError('Unexpected CH0 ionic endpoint')
        state, weight = mapping[raw.j]
        lower, upper = raw.j-1, raw.i-1
        if model.levels[lower].stage != 0 or model.levels[upper].stage != 1:
            raise ValueError('Wrong oxygen stages')
        coefficient = oracle.interpolate(grid, tables[:, state-1, 18], temperature,
                                         logarithmic=True)
        rate = QuantumPair(i=lower, j=upper)
        rate.coefficient = coefficient*1e-6*model.levels[lower].g/weight
        model.collisions[index] = rate
        count += 1
    return dict(author_data=receipt,adapter=fingerprint(__file__),CH0_replaced=count,
                branching='MN term to final J: g_J / sum_g_term (assumed statistical branching)',
                reservoir='LTE ground-H Saha H-minus using local engine constants',
                inverse='Local O LTE detailed balance; not an independent forward-table reproduction',
                scientific_admission='HOLD; author MN control, not original producer reconstruction')
