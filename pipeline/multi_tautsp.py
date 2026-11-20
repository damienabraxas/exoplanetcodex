"""Pure-Python port of MULTI 2.3's collision-table interpolation (RYA-1282).

MULTI 2.3 gencol (mul23_opacu.f, pymulti 4d6aa3c) interpolates every collision table with
de Boor's taut spline TAUTSP (Carlsson 89-08-30 variant, GAMMA = 5.5) evaluated by PPVALU:
  * on log(rate) for CE, CI, CH0, CH+ and OHM; on the linear rate for CH (and the rest);
  * the temperature is CLAMPED to the table range before evaluation (no extrapolation);
  * tables with fewer than 4 points use SPLIN at the UNCLAMPED temperature (splin_small).
A port, not a re-derivation: statement order and the single-precision literals of the
REAL(8)-expanded source (1./6. and 1.E-11 are REAL*4 constants promoted to double) are kept,
so results match the compiled Fortran oracle to rounding. Verified against the pinned oracle
(bride_native_spline.NativeTautsp) in tests/test_multi_tautsp_rya1282.py.
"""
import numpy as np

_SIXTH = float(np.float32(1.0) / np.float32(6.0))     # Fortran 1./6. (REAL*4)
_TINY_ONEMZT = float(np.float32(1e-11))                # Fortran 1.E-11 (REAL*4)
LOG_KINDS = frozenset({'CE', 'CI', 'CH0', 'CH+', 'OHM'})


def tautsp(tau, gtau, gamma=5.5):
    """Taut cubic spline. Returns (breaks, coef[4, L], L) with L = number of breakpoints,
    exactly as the Carlsson TAUTSP returns it (PPVALU searches breaks[0:L])."""
    tau = np.asarray(tau, dtype=float); gtau = np.asarray(gtau, dtype=float)
    ntau = len(tau)
    if ntau < 4 or gtau.shape != tau.shape:
        raise ValueError('TAUTSP needs >= 4 points')
    if not (np.isfinite(tau).all() and np.isfinite(gtau).all()):
        raise ValueError('TAUTSP needs finite data')
    # 1-based work arrays, as in the Fortran
    S = [[0.0] * 7 for _ in range(ntau + 2)]
    T = [0.0] + list(map(float, tau)); G = [0.0] + list(map(float, gtau))
    brk = [0.0] * (3 * ntau + 3); cf = [[0.0] * (3 * ntau + 3) for _ in range(5)]
    ntaum1 = ntau - 1
    for i in range(1, ntaum1 + 1):
        S[i][1] = T[i + 1] - T[i]
        if not S[i][1] > 0.:
            raise ValueError('TAUTSP points disordered')
        S[i + 1][4] = (G[i + 1] - G[i]) / S[i][1]
    for i in range(2, ntaum1 + 1):
        S[i][4] = S[i + 1][4] - S[i][4]
    i = 2
    S[2][2] = S[1][1] / 3.
    method = 2; gam = float(gamma)
    if gam <= 0.:
        method = 1
    if gam > 3.:
        method = 3; gam = gam - 3.
    onemg3 = 1. - gam / 3.
    alph = lambda x: min(1.0, onemg3 / x)
    zeta = onemzt = zt2 = alpha = factor = ratio = entry3 = 0.0
    while True:                                                   # label 10
        z = .5
        if method == 1:
            pass
        else:
            skip = method == 2 and S[i][4] * S[i + 1][4] < 0.
            if not skip:                                          # label 12
                temp = abs(S[i + 1][4]); denom = abs(S[i][4]) + temp
                if denom != 0.:
                    z = temp / denom
                    if abs(z - .5) <= _SIXTH:
                        z = .5
        S[i][5] = z                                               # label 19
        if z < .5:
            zeta = gam * z; onemzt = 1. - zeta; zt2 = zeta ** 2
            alpha = alph(onemzt); factor = zeta / (alpha * (zt2 - 1.) + 1.)
            S[i][6] = zeta * factor / 6.
            S[i][2] = S[i][2] + S[i][1] * ((1. - alpha * onemzt) * factor / 2. - S[i][6])
            if S[i][2] <= 0.:
                S[i][2] = 1.
            S[i][3] = S[i][1] / 6.
        elif z == .5:
            S[i][2] = S[i][2] + S[i][1] / 3.
            S[i][3] = S[i][1] / 6.
        else:
            onemzt = gam * (1. - z); zeta = 1. - onemzt
            alpha = alph(zeta); factor = onemzt / (1. - alpha * zeta * (1. + onemzt))
            S[i][6] = onemzt * factor / 6.
            S[i][2] = S[i][2] + S[i][1] / 3.
            S[i][3] = S[i][6] * S[i][1]
        if i <= 2:                                                # label 25
            S[1][5] = .5; S[1][2] = S[1][1] / 6.; S[1][3] = S[2][2]
            entry3 = S[2][3]
            if z < .5:
                factr2 = zeta * (alpha * (zt2 - 1.) + 1.) / (alpha * (zeta * zt2 - 1.) + 1.)
                ratio = factr2 * S[2][1] / S[1][2]
                S[2][2] = factr2 * S[2][1] + S[1][1]
                S[2][3] = -factr2 * S[1][1]
            elif z == .5:
                ratio = S[2][1] / S[1][2]
                S[2][2] = S[2][1] + S[1][1]
                S[2][3] = -S[1][1]
            else:
                ratio = S[2][1] / S[1][2]
                S[2][2] = S[2][1] + S[1][1]
                S[2][3] = -S[1][1] * 6. * alpha * S[2][6]
            S[2][2] = ratio * S[1][3] + S[2][2]                   # label 29
            S[2][3] = ratio * entry3 + S[2][3]
            S[1][4] = S[2][4]
            S[2][4] = ratio * S[1][4]
        else:                                                     # label 30
            S[i][2] = ratio * S[i - 1][3] + S[i][2]
            S[i][4] = ratio * S[i - 1][4] + S[i][4]
        if z < .5:                                                # label 35
            ratio = -S[i][6] * S[i][1] / S[i][2]
            S[i + 1][2] = S[i][1] / 3.
        elif z == .5:
            ratio = -(S[i][1] / 6.) / S[i][2]
            S[i + 1][2] = S[i][1] / 3.
        else:
            ratio = -(S[i][1] / 6.) / S[i][2]
            S[i + 1][2] = S[i][1] * ((1. - zeta * alpha) * factor / 2. - S[i][6])
        i = i + 1                                                 # label 40
        if not i < ntaum1:
            break
    S[i][5] = .5
    entry = ratio * S[i - 1][3] + S[i][2] + S[i][1] / 3.
    S[i + 1][2] = S[i][1] / 6.
    S[i + 1][4] = ratio * S[i - 1][4] + S[i][4]
    if z < .5:
        ratio = S[i][1] * 6. * S[i - 1][6] * alpha / S[i - 1][2]
        S[i][2] = ratio * S[i - 1][3] + S[i][1] + S[i - 1][1]
        S[i][3] = -S[i - 1][1]
    elif z == .5:
        ratio = S[i][1] / S[i - 1][2]
        S[i][2] = ratio * S[i - 1][3] + S[i][1] + S[i - 1][1]
        S[i][3] = -S[i - 1][1]
    else:
        factr2 = onemzt * (alpha * (onemzt ** 2 - 1.) + 1.) / (alpha * (onemzt ** 3 - 1.) + 1.)
        ratio = factr2 * S[i][1] / S[i - 1][2]
        S[i][2] = ratio * S[i - 1][3] + factr2 * S[i - 1][1] + S[i][1]
        S[i][3] = -factr2 * S[i - 1][1]
    S[i][4] = ratio * S[i - 1][4]                                 # label 45
    ratio = -entry / S[i][2]
    S[i + 1][2] = ratio * S[i][3] + S[i + 1][2]
    S[i + 1][4] = ratio * S[i][4] + S[i + 1][4]
    S[ntau][4] = S[ntau][4] / S[ntau][2]
    while True:                                                   # label 50
        S[i][4] = (S[i][4] - S[i][3] * S[i + 1][4]) / S[i][2]
        i = i - 1
        if not i > 1:
            break
    S[1][4] = (S[1][4] - S[1][3] * S[2][4] - entry3 * S[3][4]) / S[1][2]
    brk[1] = T[1]; L = 1
    for i in range(1, ntaum1 + 1):                                # DO 70
        cf[1][L] = G[i]; cf[3][L] = S[i][4]
        divdif = (G[i + 1] - G[i]) / S[i][1]
        z = S[i][5]
        if z < .5 and z != 0.:
            zeta = gam * z; onemzt = 1. - zeta
            c = S[i + 1][4] / 6.; d = S[i][4] * S[i][6]
            L = L + 1
            dl = zeta * S[i][1]
            brk[L] = T[i] + dl
            zt2 = zeta ** 2; alpha = alph(onemzt); factor = onemzt ** 2 * alpha
            cf[1][L] = G[i] + divdif * dl + S[i][1] ** 2 * (d * onemzt * (factor - 1.) + c * zeta * (zt2 - 1.))
            cf[2][L] = divdif + S[i][1] * (d * (1. - 3. * factor) + c * (3. * zt2 - 1.))
            cf[3][L] = 6. * (d * alpha * onemzt + c * zeta)
            cf[4][L] = 6. * (c - d * alpha) / S[i][1]
            cf[4][L - 1] = cf[4][L] - 6. * d * (1. - alpha) / (dl * zt2)
            cf[2][L - 1] = cf[2][L] - dl * (cf[3][L] - (dl / 2.) * cf[4][L - 1])
        elif z == .5:
            cf[2][L] = divdif - S[i][1] * (2. * S[i][4] + S[i + 1][4]) / 6.
            cf[4][L] = (S[i + 1][4] - S[i][4]) / S[i][1]
        elif z > .5 and abs(gam * (1. - z)) >= _TINY_ONEMZT:
            onemzt = gam * (1. - z); zeta = 1. - onemzt
            alpha = alph(zeta)
            c = S[i + 1][4] * S[i][6]; d = S[i][4] / 6.
            dl = zeta * S[i][1]
            brk[L + 1] = T[i] + dl
            cf[2][L] = divdif - S[i][1] * (2. * d + c)
            cf[4][L] = 6. * (c * alpha - d) / S[i][1]
            L = L + 1
            cf[4][L] = cf[4][L - 1] + 6. * (1. - alpha) * c / (S[i][1] * onemzt ** 3)
            cf[3][L] = cf[3][L - 1] + dl * cf[4][L - 1]
            cf[2][L] = cf[2][L - 1] + dl * (cf[3][L - 1] + (dl / 2.) * cf[4][L - 1])
            cf[1][L] = cf[1][L - 1] + dl * (cf[2][L - 1] + (dl / 2.) * (cf[3][L - 1] + (dl / 3.) * cf[4][L - 1]))
        else:                                                     # label 65 (z == 0, or tiny onemzt)
            cf[2][L] = divdif; cf[3][L] = 0.; cf[4][L] = 0.
        L = L + 1                                                 # label 68
        brk[L] = T[i + 1]                                         # label 70
    breaks = np.array(brk[1:L + 1])
    coef = np.array([[cf[m][j] for j in range(1, L + 1)] for m in range(1, 5)])
    return breaks, coef, L


def ppvalu(breaks, coef, L, x):
    """PPVALU(JDERIV=0) with INTERV's search over breaks[0:L]. Vectorised. Below the first
    break MULTI returns 0 (MFLAG=-1); above the last, 0 (MFLAG=1); exactly at it, the last
    piece. gencol never queries outside (it clamps), so outside is refused here."""
    x = np.asarray(x, dtype=float)
    if np.any(x < breaks[0]) or np.any(x > breaks[L - 1]):
        raise ValueError('PPVALU query outside the spline (MULTI would return 0); clamp first')
    left = np.searchsorted(breaks[:L], x, side='right') - 1
    left = np.minimum(left, L - 2)
    h = x - breaks[left]
    out = np.zeros_like(x)
    fmmjdr = 4.0
    for m in (3, 2, 1, 0):                                        # Horner, as label 9
        out = (out / fmmjdr) * h + coef[m, left]
        fmmjdr -= 1.
    return out


def gencol_coefficient(kind, temperature_grid, table, temperatures):
    """MULTI 2.3 gencol: the collision coefficient of one table row at `temperatures`."""
    kind = kind.upper().strip()
    grid = np.asarray(temperature_grid, dtype=float); table = np.asarray(table, dtype=float)
    log = kind in LOG_KINDS or kind.startswith('OHM')
    if log:
        if np.any(table <= 0):
            raise ValueError(f'{kind}: log interpolation needs positive coefficients')
        table = np.log(table)
    if len(grid) < 4:
        out = splin_small(temperatures, grid, table)          # gencol: unclamped T
    else:
        breaks, coef, L = tautsp(grid, table, 5.5)
        t = np.clip(np.asarray(temperatures, dtype=float), grid[0], grid[-1])
        out = ppvalu(breaks, coef, L, t)
    if log:
        out = np.exp(out)
    if kind == 'CH':
        out = np.maximum(out, 0.)        # gencol: 'ct set to zero' for CH/CHI
    return out


_F32 = lambda x: float(np.float32(x))
_SPLIN_EPS = _F32(1.0e-8)
_SPLIN_OMEGA = _F32(1.0717968)


def splin_small(t, x, y):
    """MULTI SPLIN for the tables gencol sends it (N = 1, 2, 3), evaluated at UNCLAMPED t
    exactly as gencol calls it:
      N = 1: constant.  N = 2: straight line through both points, EXTRAPOLATED outside.
      N = 3: natural cubic spline (S2 = 0 at both ends) by SOR (omega 1.0717968, stop when
             the normalised update < 1e-8); outside the table, the end value is held.
    No native oracle (the pinned oracle holds TAUTSP/PPVALU/INTERV only): checked against
    the analytic line and scipy's natural cubic in tests/test_multi_tautsp_rya1282.py."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float); t = np.asarray(t, dtype=float)
    n = len(x)
    if n == 1:
        return np.full_like(t, y[0])
    if n == 2:
        return (y[1] - y[0]) * (t - x[0]) / (x[1] - x[0]) + y[0]
    if n != 3:
        raise ValueError('gencol uses SPLIN only for N < 4')
    xnorm = 1. / (x[2] - x[0])
    ynorm = y.max() - y.min()
    if ynorm == 0:
        ynorm = y.max()
    if ynorm == 0:
        ynorm = 1.
    ynorm = 1. / ynorm
    h = [(x[i + 1] - x[i]) * xnorm for i in range(2)]
    dely = [((y[i + 1] - y[i]) / h[i]) * ynorm for i in range(2)]
    h2 = h[0] + h[1]; b = 0.5 * h[0] / h2
    delsqy = (dely[1] - dely[0]) / h2
    s2 = [0., 2. * delsqy, 0.]; c = 3. * delsqy
    while True:                                   # SOR, one interior point
        w = (c - b * s2[0] - (0.5 - b) * s2[2] - s2[1]) * _SPLIN_OMEGA
        eta = abs(w)
        s2[1] = s2[1] + w
        if eta < _SPLIN_EPS:
            break
    s3 = [(s2[i + 1] - s2[i]) / h[i] for i in range(2)]
    out = np.empty_like(t)
    for k, tk in np.ndenumerate(t):
        if (tk - x[0]) * xnorm <= 0:
            out[k] = y[0]; continue
        if (tk - x[2]) * xnorm > 0:
            out[k] = y[2]; continue
        i = 0 if tk < x[1] else 1
        if tk == x[1]:
            out[k] = y[1]; continue
        ht1 = (tk - x[i]) * xnorm; ht2 = (tk - x[i + 1]) * xnorm
        ss2 = s2[i] + ht1 * s3[i]; prod = ht1 * ht2
        delsqs = (s2[i] + s2[i + 1] + ss2) / 6.
        out[k] = y[i] + (ht1 * dely[i] + prod * delsqs) / ynorm
    return out
