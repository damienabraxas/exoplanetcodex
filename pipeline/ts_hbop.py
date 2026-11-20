"""Turbospectrum's hydrogen bound-free continuum (hbop.f, contonly=.true.), ported.

RYA-1282, 2026-11-20. TS absko adds, AFTER detabs, Barklem & Piskunov's HBOP continuum:
H I bound-free from the 6 lowest levels with the occupation-probability formalism
(Hummer & Mihalas 1988; Hubeny, Hummer & Lanz 1994 ionic fit). Redward of each edge the
cross-section is extrapolated (Karzas-Latter Gaunt factors, log-quadratic) and multiplied
by the dissolved fraction D = 1 - w(n*)/w(n), with an exponential taper redward of the
level's alpha line. Since BPz 2025-08-25 the n=1 extrapolation runs to 3647 A, so in cool,
metal-poor photospheres the dissolved Lyman continuum is most of the 150-200 nm absorption.
LTE only (NL = 0), as absko calls it. Lines are NOT included (TS adds them separately).

Units as hbop: wavelength in A (TS convention: vacuum below 2000 A, air above), number
densities in cm^-3, result in cm^-1 (stimulated emission included).
Validated against the compiled Fortran (tests/test_ts_hbop_rya1282.py).
"""
import numpy as np

NLEVELS = 100
NBF = 6
_H, _C, _K = 6.62607e-27, 2.9979e10, 1.38065e-16
_HC = _H * _C
_IONH = 2.17991e-11
_A0 = 5.29177e-9

EHYD = np.empty(NLEVELS + 1)            # 1-based, cm^-1
EHYD[1:9] = [0.0, 82259.105, 97492.302, 102823.893, 105291.651, 106632.160, 107440.444, 107965.051]
for _i in range(9, NLEVELS + 1):
    EHYD[_i] = 109678.764 - 109677.576 / float(_i * _i)

# Karzas & Latter Gaunt factors from MARCS (hbop.f HBF), levels 1-7 (7 only for the n=6 limit)
_L = {1: [182., 241., 300., 356., 408., 456., 538., 631., 729., 821., 877., 912.],
      2: [215., 301., 398., 503., 614., 729., 965., 1313., 1824., 2525., 3144., 3647.],
      3: [222., 316., 424., 545., 677., 821., 1132., 1641., 2525., 4103., 6034., 6933., 7529., 8207.],
      4: [224., 321., 434., 561., 703., 859., 1205., 1799., 2918., 5252., 8895., 10998., 12576., 14588.],
      5: [226., 324., 438., 569., 715., 877., 1241., 1882., 3144., 6034., 11397., 15095., 18235., 22794.],
      6: [226., 325., 441., 573., 722., 887., 1262., 1931., 3282., 6564., 13448., 18916., 24120., 32915.],
      7: [227., 326., 442., 576., 726., 894., 1275., 1961., 3372., 6933., 15095., 22347., 29992., 36616., 44694.]}
_G = {1: [.9939, .9948, .9856, .9721, .9571, .9423, .9157, .8850, .8531, .8246, .8076, .7973],
      2: [1.043, 1.058, 1.063, 1.061, 1.056, 1.049, 1.033, 1.009, .9755, .9342, .9011, .8761],
      3: [1.053, 1.072, 1.081, 1.083, 1.081, 1.078, 1.068, 1.051, 1.024, .9854, .9458, .9293, .9189, .9075],
      4: [1.056, 1.077, 1.087, 1.091, 1.091, 1.089, 1.082, 1.069, 1.047, 1.014, .9736, .9542, .9408, .9247],
      5: [1.058, 1.080, 1.090, 1.095, 1.095, 1.094, 1.089, 1.078, 1.060, 1.030, .9925, .9719, .9563, .9358],
      6: [1.059, 1.081, 1.092, 1.097, 1.098, 1.097, 1.092, 1.083, 1.067, 1.041, 1.006, .9851, .9682, .9436]}


def wcalc(nh, ne, nhe, ns, temp):
    """Occupation probability (hbop.f WCALC); arrays broadcast."""
    nh, ne, nhe, ns, temp = (np.asarray(a, dtype=float) for a in (nh, ne, nhe, ns, temp))
    ns2 = ns * ns; ns4 = ns2 * ns2
    rih = np.sqrt(2.5 * ns4 + 0.5 * ns2) * _A0
    x1 = rih + 1.73 * _A0; x2 = rih + 1.02 * _A0
    wneutr = np.exp(-4.18879 * (nh * x1 ** 3 + nhe * x2 ** 3))
    k = np.where(ns > 3., 5.33333333 * ns / (ns + 1.) / (ns + 1.), 1.)
    ok = (ne > 10.) & (temp > 10.)
    nes = np.where(ok, ne, 1e3); ts = np.where(ok, temp, 1e3)
    a = 0.09 * np.exp(0.16667 * np.log(nes)) / np.sqrt(ts)
    x = np.exp(3.15 * np.log(1. + a))
    betac = 8.3e14 * np.exp(-0.66667 * np.log(nes)) * k / ns4
    f = 0.1402 * x * betac ** 3 / (1. + 0.1285 * x * betac * np.sqrt(betac))
    wion = np.where(ok, f / (1. + f), 1.)
    return wion * wneutr


def _logint(xs, ys, xi):
    """hbop.f LOGINT: 3-point polynomial in log10(y), points chosen as the Fortran does."""
    xs = np.asarray(xs); yl = np.log10(np.asarray(ys)); n = len(xs)
    xi = np.asarray(xi, dtype=float)
    j = np.searchsorted(xs, xi, side='left')       # count of X(I) < XI  == Fortran J (1-based)
    j = np.where(j == 0, 1, j)
    j = np.where(j > n - 3, n - 2, j)
    k = j - 1                                       # 0-based first point
    x0, x1, x2 = xs[k], xs[k + 1], xs[k + 2]
    y0, y1, y2 = yl[k], yl[k + 1], yl[k + 2]
    yi = (y0 * (xi - x1) * (xi - x2) / ((x0 - x1) * (x0 - x2))
          + y1 * (xi - x0) * (xi - x2) / ((x1 - x0) * (x1 - x2))
          + y2 * (xi - x0) * (xi - x1) / ((x2 - x0) * (x2 - x1)))
    return 10. ** yi


def hbf(wave_A, nlo):
    """Bound-free cross-section of level nlo (cm^2), hbop.f HBF (with the 2025 edge limit)."""
    w = np.asarray(wave_A, dtype=float)
    g = np.where(w < _L[nlo + 1][-1], _logint(_L[nlo], _G[nlo], w), 0.)
    return 1.044e-26 * g * w ** 3 / float(nlo) ** 5


def hbop_continuum(wave_A, T, nh, nhe, ne):
    """TS hbop CONTIN, LTE: (nwave, ndep) in cm^-1. nh = n(H I), nhe = n(He I), ne (cm^-3)."""
    w = np.asarray(wave_A, dtype=float)[:, None]
    T = np.asarray(T, dtype=float)[None, :]
    nh, nhe, ne = (np.asarray(a, dtype=float)[None, :] for a in (nh, nhe, ne))
    kt = _K * T
    lev = np.arange(1, NLEVELS + 1, dtype=float)[:, None, None]
    W = wcalc(nh, ne, nhe, lev, T)                                      # (NLEVELS, 1, ndep)
    G = 2. * lev * lev
    wge = W * G * np.exp(-_HC * EHYD[1:, None, None] / kt)
    nlte = nh * wge / wge.sum(axis=0)                                   # LTE populations
    out = np.zeros(np.broadcast_shapes(w.shape, T.shape))
    stim = np.exp(-_HC / (w * 1e-8) / kt)
    for i in range(1, NBF + 1):
        x = 1. / (i * i) - _HC / (w * 1e-8 * _IONH)                     # > 0 redward of the edge
        red = x > 0
        ns = np.where(red, 1. / np.sqrt(np.where(red, x, 1.)), 1.)
        wstar = wcalc(nh, ne, nhe, ns, T)
        d = np.where(red, 1. - wstar / W[i - 1], 1.)
        redcut = 1e8 / (EHYD[i + 1] - EHYD[i])
        taper = np.where(red & (w > redcut), np.exp(np.minimum(redcut - w, 0.) * 1e-2), 1.)
        n_i = nlte[i - 1]
        out += (n_i - n_i * stim) * d * hbf(w, i) * taper
    return out
