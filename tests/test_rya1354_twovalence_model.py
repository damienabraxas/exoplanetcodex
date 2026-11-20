"""RYA-1354 v65-v68: invariants of the two-valence-electron model (fast; no full CI)."""
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import rya1354_onevalence_model as M  # noqa: E402
import rya1354_twovalence_model as T  # noqa: E402


def test_spin_orbit_reproduces_lande_interval_rule():
    """3s3p 3P with xi_p = 1: E_J = (A/2)[J(J+1) - L(L+1) - S(S+1)], A = xi/2 -> +1/2, -1/2, -1 for J = 2, 1, 0."""
    cf = [(0, 0, 1, 0)]
    W = T.so_matrix(None, {1: np.array([[1.0]])}, cf, 1, 1, cf, 1, 1)[0, 0]
    for J, exact in ((2, 0.5), (1, -0.5), (0, -1.0)):
        e = (-1) ** (1 + 1 + J) * T.sixj(1, 1, J, 1, 1, 1) * W
        assert e == pytest.approx(exact, abs=1e-12)


def test_spin_orbit_one_electron_l_dot_s():
    """Racah form of <l s j|l.s|l s j> equals [j(j+1) - l(l+1) - 3/4]/2."""
    for l in (1, 2, 3):
        s = Fr(1, 2)
        for j in (l - s, l + s):
            v = (-1) ** int(l + s + j) * M.sixj(Fr(l), s, j, s, Fr(l), Fr(1)) * math.sqrt(l * (l + 1) * (2 * l + 1)) * math.sqrt(1.5)
            assert v == pytest.approx(float((j * (j + 1) - l * (l + 1) - s * (s + 1)) / 2), abs=1e-12)


@pytest.mark.parametrize('L1,L2,S', [(0, 1, 0), (0, 1, 1), (1, 2, 1), (2, 3, 1), (1, 1, 1)])
def test_fine_structure_fractions_sum_to_multiplicity(L1, L2, S):
    """sum_JJ' (2J+1)(2J'+1){L1 J S; J' L2 1}^2 = 2S+1: J-resolved gf add up to the LS multiplet."""
    tot = sum((2 * J1 + 1) * (2 * J2 + 1) * T.sixj(L1, J1, S, J2, L2, 1) ** 2
              for J1 in range(abs(L1 - S), L1 + S + 1) for J2 in range(abs(L2 - S), L2 + S + 1))
    assert tot == pytest.approx(2 * S + 1, abs=1e-12)


def _ion_states(l, p, n):
    M.use('Mg II')
    M.X0, M.X1, M.NPTS = math.log(1e-5), math.log(1200.0), 8000
    e, z, _ = M.solve(l, p, n)
    return e, z


def test_velocity_equals_length_times_energy_for_a_local_potential():
    """Gauge identity <a|nabla|c> = (E_c - E_a)<a|r|c> holds to 1e-5 when s and p share one (local) potential."""
    fit = json.load(open(ROOT / 'data/audit/rya1354/mg2_model_v63/summary.json'))['potential_fit']
    p0 = tuple(fit['s']['params'])
    es, zs = _ion_states(0, p0, 4)
    ep, zp = _ion_states(1, p0, 4)
    x, r, dx = M.grid()
    B = SimpleNamespace(orb={0: (es, zs), 1: (ep, zp)}, x=x, r=r, dx=dx)
    V = T.velocity_radial(B, 0, 1)
    for i, j in ((1, 1), (2, 1), (1, 2)):
        RL = (zs[:, i] * zp[:, j] * r * dx).sum()
        assert V[i, j] == pytest.approx((ep[j] - es[i]) * RL, rel=1e-5)


def test_spin_orbit_is_linear_in_each_l():
    so_l = {1: {'a': 1.0, 'b': 2.0}, 2: {'a': 0.5}}
    assert T.combine_so(so_l, {1: 2.0, 2: 4.0}) == {'a': 4.0, 'b': 4.0}
    assert T.combine_so(so_l, {}) == {'a': 1.5, 'b': 2.0}


def test_mg1_ls_broken_configurations_are_nf_ng_not_nd():
    b = T.ls_broken(dict(T.SPECIES['Mg I']))
    assert {'3s.4f', '3s.5f', '3s.5g'} <= b
    assert not {'3s.3d', '3s.4d', '3s.3p', '3s.4p'} & b


def test_lab_labels_carry_odd_parity_star():
    lab = T.lab_rows(dict(T.SPECIES['Mg I']), '3s')
    keys = {(x['lo'], x['up']) for x in lab}
    assert (('3s.3p', '3P*', 0), ('3s.4s', '3S', 1)) in keys
    assert (('3s.3d', '1D', 2), ('3s.4f', '3F*', 3)) in keys


def test_parse_conf_maps_nist_configurations_to_ci_keys():
    nmin = {0: 3, 1: 3, 2: 3, 3: 4}
    assert T.parse_conf('3p2', nmin, '2p6') == (1, 0, 1, 0)
    assert T.parse_conf('3p.3d', nmin, '2p6') == (1, 0, 2, 0)
    assert T.parse_conf('2p6.3s2', nmin, '2p6') == (0, 0, 0, 0)
    assert T.parse_conf('3s.4f', nmin, '2p6') == (0, 0, 3, 0)
    assert T.parse_conf('3s.10d', nmin, '2p6') == (0, 0, 2, 7)


def test_equivalent_electron_terms_obey_pauli():
    """l^2 allows only L + S even: d^2 gives 1S 3P 1D 3F 1G, so a '3d2 3S' component cannot exist (NIST lists one)."""
    allowed = {(L, S) for L in range(5) for S in (0, 1) if (L + S) % 2 == 0}
    assert {(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)} == allowed
    assert (0, 1) not in allowed


def test_v69_convergence_rule_excludes_basis_degrading_variants():
    assert T.CONV_ENERGY_SLACK == pytest.approx(0.001)
    assert T.VARIANTS['rcomp20']['group'] == 'conv'


def test_sister_ion_lab_labels_accept_full_configurations():
    lab = T.lab_rows(dict(T.SPECIES['Al II']), '3s')
    keys = {(x['lo'], x['up']) for x in lab}
    assert (('3s2', '1S', 0), ('3s.3p', '1P*', 1)) in keys          # 1670 A resonance line
    assert any(x['lo'][0] == '3s.3p' and x['up'][0] == '3p2' for x in lab)


# ---------------------------------------------------------------- v70: broadening + lifetimes + d-scale rule
def test_abo_tables_reproduce_ges_packed_values():
    """Implementation control (declared v70): the b-triplet and 5528 ABO values carried by GES within 2 %."""
    import math
    import rya1354_vdw as V
    T_ = V.load_tables()
    ip = 61671.05
    ns = lambda E: math.sqrt(V.R_MG / (ip - E))
    s, a = V.abo(T_, 'p', ns(21877.0), 's', ns(41197.403))       # 3s3p 3P - 3s4s 3S
    assert abs(s / 729 - 1) < 0.02 and abs(a - 0.238) < 0.01
    s, a = V.abo(T_, 'p', ns(35051.264), 'd', ns(53134.642))      # 5528: 3s3p 1P - 3s4d 1D
    assert abs(s / 1461 - 1) < 0.02 and abs(a - 0.312) < 0.01
    assert V.abo(T_, 'p', ns(35051.264), 'd', ns(56308.381)) == (None, None)   # 4703: 5d beyond the table edge
    # sigma, alpha -> log gamma6 at 1e4 K agrees with the GES single-gamma column (-7.267)
    assert abs(V.gamma6_abo(729, 0.238, 24.305) - (-7.267)) < 0.01


def test_hydrogenic_r2_inverse():
    import rya1354_vdw as V
    for n, l in ((2.3, 0), (3.7, 1), (4.2, 2), (6.0, 3)):
        assert abs(V.nstar_from_r2(V.r2_hyd(n, l), l) - n) < 1e-8


def test_unsold_monotone_in_dr2():
    import rya1354_vdw as V
    assert V.gamma6_unsold(10, 24.3) < V.gamma6_unsold(100, 24.3)
    # gamma6 ~ C6^(2/5): a factor 10 in d<r2> is 0.4 dex
    assert abs(V.gamma6_unsold(100, 24.3) - V.gamma6_unsold(10, 24.3) - 0.4) < 1e-9


def test_inverted_core_fs_rule_is_on():
    assert T.INVERTED_FS_RULE is True
