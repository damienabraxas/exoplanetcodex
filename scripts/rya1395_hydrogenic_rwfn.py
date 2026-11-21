"""RYA-1395 v72: Z_eff-hydrogenic radial orbitals written as a GRASP G92RWF file, plus a reader.

Why (v71 FAILED_RETAINED): GRASP rwfnestimate's "screened hydrogenic" option takes its screening from the CSF
list's CORE subshells; our s0 lists declare none, so it built bare Z=8 orbitals (4f eps 2.000 Ha, <r> 2.25).
Here P(r) is the nonrelativistic hydrogenic radial function for Z_eff (default 1, the Rydberg limit of a neutral
atom) and Q = (dP/dr + kappa P / r) / (2c) (Pauli small component, GRASP sign convention), normalised so that
int (P^2 + Q^2) dr = 1, on GRASP's log grid R_i = RNT (exp((i-1) H) - 1).

G92RWF layout (Fortran unformatted, 4-byte record markers): 'G92RWF'; per orbital three records:
(NP int32, NAK int32, E float64 [binding, Hartree], MY int32); (P0, P(1:MY), Q(1:MY)) float64; R(1:MY) float64.
"""
import math
import struct
from math import factorial

import numpy as np
from scipy.special import eval_genlaguerre

C_AU = 137.035999084  # CODATA 2018 speed of light in atomic units
L_OF = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}


def kappa(l, j_minus):
    """GRASP NAK: j = l - 1/2 ('-' label) -> +l ; j = l + 1/2 -> -(l+1)."""
    return l if j_minus else -(l + 1)


def grasp_grid(rnt, h, n):
    i = np.arange(n)
    return rnt * np.expm1(i * h)


def hydrogenic(nn, l, r, zeff=1.0):
    rho = 2.0 * zeff * r / nn
    norm = math.sqrt((2.0 * zeff / nn) ** 3 * factorial(nn - l - 1) / (2 * nn * factorial(nn + l)))
    return r * norm * np.exp(-rho / 2) * rho ** l * eval_genlaguerre(nn - l - 1, 2 * l + 1, rho)


def orbital(nn, l, j_minus, r, zeff=1.0, tail=1e-14):
    """(nak, eps, my, P0, P, Q, R) truncated where |P| and |Q| stay below tail * max|P| to the end."""
    k = kappa(l, j_minus)
    P = hydrogenic(nn, l, r, zeff)
    dP = np.gradient(P, r, edge_order=2)
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = np.where(r > 0, (dP + k * P / np.where(r > 0, r, 1.0)) / (2 * C_AU), 0.0)
    dr = np.gradient(r)
    s = math.sqrt(np.sum((P ** 2 + Q ** 2) * dr))
    P, Q = P / s, Q / s
    big = np.abs(P) > tail * np.abs(P).max()
    my = int(np.nonzero(big)[0].max()) + 2
    my = min(my, len(r))
    i0 = 5  # leading coefficient from a small-r point (P ~ P0 r^(l+1))
    p0 = float(P[i0] / r[i0] ** (l + 1))
    eps = zeff ** 2 / (2.0 * nn ** 2)
    return k, eps, my, p0, P[:my], Q[:my], r[:my]


def _rec(b):
    return struct.pack('<i', len(b)) + b + struct.pack('<i', len(b))


def write_rwfn(path, orbitals):
    """orbitals: list of (n, nak, eps, my, p0, P, Q, R)."""
    out = _rec(b'G92RWF')
    for nn, nak, eps, my, p0, P, Q, R in orbitals:
        out += _rec(struct.pack('<iidi', nn, nak, eps, my))
        out += _rec(np.concatenate([[p0], P, Q]).astype('<f8').tobytes())
        out += _rec(np.asarray(R, dtype='<f8').tobytes())
    open(path, 'wb').write(out)


def read_rwfn(path):
    b = open(path, 'rb').read()
    recs, i = [], 0
    while i < len(b):
        n = struct.unpack('<i', b[i:i + 4])[0]
        recs.append(b[i + 4:i + 4 + n])
        i += 8 + n
    if recs[0] != b'G92RWF':
        raise ValueError(f'{path}: not a G92RWF file')
    orbs = []
    for k in range(1, len(recs) - 2, 3):
        nn, nak, eps, my = struct.unpack('<iidi', recs[k][:20])
        pq = np.frombuffer(recs[k + 1], dtype='<f8')
        R = np.frombuffer(recs[k + 2], dtype='<f8')[:my]
        orbs.append(dict(n=nn, nak=nak, eps=eps, my=my, p0=pq[0], P=pq[1:my + 1], Q=pq[my + 1:2 * my + 1], R=R))
    return orbs


def label(o):
    l = o['nak'] if o['nak'] > 0 else -o['nak'] - 1
    return f"{o['n']}{'spdfg'[l]}{'-' if o['nak'] > 0 else ''}"


def mean_r(o):
    w = np.gradient(o['R'])
    rho = o['P'] ** 2 + o['Q'] ** 2
    return float(np.sum(rho * o['R'] * w) / np.sum(rho * w))


def overlap(a, b):
    n = min(a['my'], b['my'])
    if not np.allclose(a['R'][:n], b['R'][:n], rtol=1e-10, atol=0):
        raise ValueError('overlap needs a common grid')
    w = np.gradient(a['R'][:n])
    return float(np.sum((a['P'][:n] * b['P'][:n] + a['Q'][:n] * b['Q'][:n]) * w))


def build(nl_list, rnt, h, n, zeff=1.0):
    """nl_list like ['4f', '5f']: both j components of each -> list for write_rwfn."""
    r = grasp_grid(rnt, h, n)
    out = []
    for nl in nl_list:
        nn, l = int(nl[:-1]), L_OF[nl[-1]]
        for j_minus in ((True, False) if l > 0 else (False,)):
            k, eps, my, p0, P, Q, R = orbital(nn, l, j_minus, r, zeff)
            out.append((nn, k, eps, my, p0, P, Q, R))
    return out
