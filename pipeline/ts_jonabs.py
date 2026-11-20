"""Read and evaluate one component of Turbospectrum's continuous-opacity table (jonabs).

TS (and MARCS) take bound-free metal opacities from jonabs_vac_*.dat: per component, a
table XKAP(lambda, T) read by inabs.f and interpolated by inabs.f (lambda) then tabs.f (T);
detabs.f multiplies it by FAKT, e.g. Fe I (component 17): P(Fe I) / (rho k T U(Fe I)), i.e.
the table is SUM_levels g_i exp(-E_i/kT) sigma_i per (n / U) [cm^2], without stimulated
emission (detabs multiplies the total by 1 - exp(-hc/lambda kT) at the end).

Port of the interpolation (iSpec TS source, inabs.f / tabs.f):
  lambda: table nodes above 2000 A converted VACUUM -> AIR (vacair), queries in air.
          KVADL=0 linear; ILOGL=0 on the values; outside the table: MINEX/MAXEX = 0 -> zero.
  T:      ITETA=0 -> T itself (else theta = 5040/T). Equidistant SEGMENTS (a new segment
          wherever the step changes by >= 1e-4). ILOGT=1 -> log values (zeros -> 1e-37
          first). KVADT=1 -> 3-point Lagrange in the segment (last point shifts back);
          below the table MINET>0 linear extrapolation, else zero; above, MAXET likewise.
Used by pipeline.bride_fe_background.JonabsFeIBackground (RYA-1282).
"""
from pathlib import Path
import numpy as np


def _air(wave_vac_A):
    """TS vacair (Morton/Edlen as in iSpec TS): only applied above 2000 A."""
    from pipeline.wavelength_util import vac_to_air
    return vac_to_air(np.asarray(wave_vac_A, dtype=float))


def read_component(path, label):
    """Return the raw block (header flags, lambda nodes, T nodes, XKAP[lambda, T])."""
    lines = Path(path).read_text(errors='replace').splitlines()
    start = [i for i, l in enumerate(lines) if l.startswith(label)]
    if len(start) != 1:
        raise ValueError(f'{label!r}: expected exactly one jonabs component, found {len(start)}')
    i = start[0] + 1

    def flags(line):
        """'ILOGL=0   KVADL=0 ... NLATB=38' -> dict (values may follow '=' after spaces)."""
        parts = line.replace('=', ' = ').split()
        return {parts[j - 1]: int(parts[j + 1]) for j, tok in enumerate(parts) if tok == '='}

    def field(text):
        """One Fortran 10-character field, read as TS reads it (FORMAT 6F10.0 / 6E10.3 with
        blanks ignored): '  7.51 E08' = 7.51e8; a field with no digits ('   .   E  ') = 0."""
        t = text.replace(' ', '').replace('D', 'E')
        if not any(ch.isdigit() for ch in t):
            return 0.0
        if t.endswith('E') or t.endswith('E+') or t.endswith('E-'):
            t = t.rstrip('E+-')
        return float(t)

    def numbers(n):
        nonlocal i
        vals = []
        while len(vals) < n:
            take = min(6, n - len(vals))
            vals += [field(lines[i][10 * k:10 * (k + 1)]) for k in range(take)]
            i += 1
        return np.array(vals)

    lf = flags(lines[i]); i += 1
    lam = numbers(lf['NLATB'])
    tf = flags(lines[i]); i += 1
    ntet = tf['NTETB']
    tet = numbers(ntet) if ntet > 1 else np.array([np.nan])
    xkap = np.column_stack([numbers(lf['NLATB']) for _ in range(ntet)])
    return dict(label=label, lflags=lf, tflags=tf, lam_vac_A=lam, tet=tet, xkap=xkap)


class JonabsComponent:
    """One jonabs component, evaluated exactly as inabs.f + tabs.f do (see module doc)."""

    def __init__(self, path, label):
        b = read_component(path, label)
        lf, tf = b['lflags'], b['tflags']
        self.label = label
        self.minex, self.maxex = lf['MINEX'], lf['MAXEX']
        self.ilogl, self.kvadl = lf['ILOGL'], lf['KVADL']
        lam = b['lam_vac_A'].copy()
        # inabs: a discontinuity = two table points < 2 A apart (on the VACUUM table)
        close = np.diff(lam) < 2.
        self.idiskv = np.zeros(len(lam), bool)
        self.idiskv[:-1] |= close; self.idiskv[1:] |= close
        self.idisk = bool(close.any())
        hi = lam > 2000.
        lam[hi] = _air(lam[hi])
        self.lam_air_A = lam                               # TS works in air
        self.ilogt, self.kvadt = tf['ILOGT'], tf['KVADT']
        self.minet, self.maxet, self.iteta = tf['MINET'], tf['MAXET'], tf['ITETA']
        self.tet = b['tet']; self.xkap = b['xkap']
        if self.ilogl >= 1:                      # inabs: log BEFORE the lambda interpolation
            self.xkap = np.log(np.where(self.xkap > 0, self.xkap, 1e-37))
        if len(self.tet) < 2:
            raise ValueError(f'{label}: temperature-independent components not ported')
        # segments (inabs.f): TBOT, DELT, NTM
        tet = self.tet
        tbot, delt, ntm = [tet[0]], [tet[1] - tet[0]], [1]
        lf_ = 1
        for L in range(2, len(tet)):
            diff = tet[L] - tet[L - 1]
            if abs(1. - diff / delt[-1]) >= 1e-4:
                tbot.append(tet[L - 1]); delt.append(diff)
                ntm[-1] = lf_; ntm.append(0); lf_ = 0
            lf_ += 1
        ntm[-1] = lf_
        self.tbot, self.delt, self.ntm = tbot, delt, ntm

    # -- lambda (inabs.f) -------------------------------------------------------------
    def _lambda_values(self, wave_air_A):
        """XKAP interpolated to each wavelength, (nwave, nT), exactly as inabs.f: linear, or
        3-point Lagrange (KVADL >= 1) falling back to linear beside a table discontinuity;
        in log if ILOGL (exponentiated here unless ILOGT keeps it in log for tabs.f), or
        logged here if ILOGT only. Rows TS drops (outside, no extrapolation) are flagged."""
        x = np.asarray(wave_air_A, dtype=float); lt = self.lam_air_A; n = len(lt)
        out = np.zeros((len(x), self.xkap.shape[1]))
        live = np.ones(len(x), bool)
        dv = self.idiskv
        for j, w in enumerate(x):
            # 1-based lambi/ihelp as inabs
            if w < lt[0]:
                ihelp, lambi = 1, 1
            elif w >= lt[-1]:
                ihelp, lambi = n, n
            else:
                lambi = int(np.searchsorted(lt, w, side='right'))        # lt[lambi-1] <= w < lt[lambi]
                ihelp = lambi + 1
            mode = None
            if ihelp <= 1:
                if self.minex > 0:
                    mode = 'linear'
            elif lambi == n:
                if self.maxex > 0:
                    lambi -= 1; mode = 'linear'
            else:
                if self.kvadl < 1:
                    mode = 'linear'
                else:
                    lambi = min(n - 2, lambi)
                    if (not self.idisk) or not dv[lambi]:              # idiskv(lambi+1), 0-based lambi
                        mode = 'quadratic'
                    elif w > lt[lambi]:                                # xla > xlatb(lambi+1)
                        lambi += 1
                        mode = 'linear' if (dv[lambi] or lambi + 1 == n) else 'quadratic'
                    elif dv[lambi - 1] or lambi == 1:
                        mode = 'linear'
                    else:
                        lambi -= 1; mode = 'quadratic'
            if mode is None:
                live[j] = False
                continue
            k = lambi - 1                                              # 0-based first point
            if mode == 'linear':
                a2 = (w - lt[k]) / (lt[k + 1] - lt[k])
                out[j] = (1. - a2) * self.xkap[k] + a2 * self.xkap[k + 1]
            else:
                d1, d2, d3 = w - lt[k], w - lt[k + 1], w - lt[k + 2]
                x21, x32, x31 = lt[k + 1] - lt[k], lt[k + 2] - lt[k + 1], lt[k + 2] - lt[k]
                a1 = d2 * d3 / (x21 * x31); a2 = d1 * d3 / (x21 * x32); a3 = d1 * d2 / (x31 * x32)
                out[j] = a1 * self.xkap[k] - a2 * self.xkap[k + 1] + a3 * self.xkap[k + 2]
        if self.ilogl >= 1:
            if self.ilogt <= 0:
                out = np.exp(out)
        elif self.ilogt > 0:                     # inabs: zero -> 1e-37, then ALOG
            out = np.log(np.where(out > 0, out, 1e-37))
        return out, live

    # -- temperature (tabs.f) ---------------------------------------------------------
    def _t_weights(self, T):
        """(index0, weights) per temperature, or None where TS returns zero."""
        ts = 5040. / T if self.iteta > 0 else T
        if ts < self.tbot[0]:
            if self.minet <= 0:
                return None
            ap = (ts - self.tbot[0]) / self.delt[0]
            return 0, (1. - ap, ap)
        intap = 1
        for seg in range(len(self.tbot)):
            ap = (ts - self.tbot[seg]) / self.delt[seg]
            ip = int(ap)
            inta = ip + intap
            inap = self.ntm[seg] - 1 + intap
            if inta <= inap:
                if self.kvadt <= 0 or self.ntm[seg] == 1:
                    a2 = ap - ip
                    return inta - 1, (1. - a2, a2)
                if inta >= inap:
                    inta -= 1; ip -= 1
                d1 = ap - ip; d2 = d1 - 1.; d3 = d1 - 2.
                return inta - 1, (d2 * d3 * .5, -d1 * d3, d1 * d2 * .5)
            intap = inap + 1
        if self.maxet <= 0:
            return None
        ip = self.ntm[-1] - 1
        a2 = ap - ip
        return inap - 1, (1. - a2, a2)

    def __call__(self, wave_vac_nm, temperature_K):
        """Table value [cm^2 per (n/U)] at vacuum nm (queried in air, as TS) and T."""
        wave_air_A = _air(np.asarray(wave_vac_nm, dtype=float) * 10.)
        vals, live = self._lambda_values(wave_air_A)
        T = np.atleast_1d(np.asarray(temperature_K, dtype=float))
        out = np.zeros((len(wave_air_A), len(T)))
        for d, t in enumerate(T):
            w = self._t_weights(t)
            if w is None:
                continue
            k, wt = w
            v = sum(wi * vals[:, k + m] for m, wi in enumerate(wt))
            out[:, d] = np.exp(v) if self.ilogt > 0 else v
        out[~live] = 0.
        return out


def _irwin_row(path, Z, ion):
    code = f'{Z}.{ion:02d}'
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        p = line.split()
        if p[1] == code:
            rest = line.split("'")[2].split()
            npart = int(rest[1])
            return float(rest[0]), [float(x.replace('D', 'E').replace('d', 'E')) for x in rest[2:2 + npart]]
    raise ValueError(f'No Irwin partition function for {code}')


def irwin_ionisation_eV(path, Z, ion):
    """Ionisation energy of stage `ion` (0 = neutral) as TS holds it: the IP column of the
    NEXT stage's row (e.g. 'fe+' 7.902 eV is Fe I -> Fe II)."""
    return _irwin_row(path, Z, ion + 1)[0]


def irwin_partition(path, Z, ion, temperature_K):
    """TS partf.f atomic partition function from IRWIN_atoms_*.dat: ln U = sum a_i (ln T)^i
    for 1000 <= T <= 16000 K; T < 1000 -> U(1000); T > 16000 -> ln U linear in T through
    U(10000) and U(16000). `ion` 0 = neutral. Matches TS's own Fe I U to 2.6e-6 (packet)."""
    row = _irwin_row(path, Z, ion)[1]
    lnU = lambda T: np.polyval(row[::-1], np.log(T))
    T = np.asarray(temperature_K, dtype=float)
    lo = np.maximum(T, 1000.)
    out = lnU(np.minimum(lo, 16000.))
    hot = T > 16000.
    if np.any(hot):
        u1, u2 = lnU(10000.), lnU(16000.)
        out = np.where(hot, (u2 - u1) * (T - 16000.) / 6000. + u2, out)
    return np.exp(out)
