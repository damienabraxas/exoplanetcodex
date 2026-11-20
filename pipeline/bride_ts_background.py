"""Bride passive metal bound-free background from Turbospectrum/MARCS's own OP tables (RYA-1282).

Lightweaver's passive RH atoms carry a handful of metal continua each. Against TS (which
reproduces MARCS absorption to 0.12% at 277 nm, RYA-1281 ts_fe_continuum_review) they leave
the 200-350 nm background 1.5-15x too weak (Fe I, Al I, Mg I dominant) - exactly the range
that photoionises O I 3s 5S (edge 277 nm) and 3p 5P (431 nm).

For each passive metal stage TS tabulates (C, Mg, Al, Si, Ca, Fe; I and II) this provider
REPLACES Lightweaver's bound-free extinction of that stage, never adds twice:
    chi_TS = XKAP(lambda, T) * n(stage) / U_Irwin(T) * (1 - exp(-hc / lambda k T))
XKAP from jonabs (pipeline.ts_jonabs: inabs/tabs ported; Fe I matches native TS to 1.3e-4);
U from TS's Irwin fits (Fe I matches TS to 2.6e-6); n(stage) is the Bride's own LTE EOS.
A neutral stage the passive model lacks (Ca I: Lightweaver carries Ca II/III only) is ADDED,
its n(I)/U(I) by Saha from the model's own next stage with TS's Irwin U and ionisation energy;
the model's own populations are not altered (Ca I/Ca II < 1e-2 at the solar top).
Emissivity changes by delta-chi * B_nu (LTE background); scattering is unchanged. Active
atoms (O) are never touched. Where a TS table is zero/undefined, TS has no opacity, so the
replacement gives zero there too (TS = MARCS producer convention).

Extra TS components (RYA-1282, 2026-11-20; detabs.f multipliers, cgs):
    H+H  Doyle 1968 (cut at 1750 A): XKAP n(H I)^2 1e-50 x stim
    H2+  Mihalas:                    XKAP n(H I) n(H II) 1e-40 x stim
    CH, OH Kurucz (cm^2/g, stim already in the table): XKAP n(X) m(X)
Lightweaver's background carries none of them. In metal-poor photospheres they are most of
the 175-277 nm absorption (H+H 39-66%, OH 20-55%, H2+ 4-16% at [Fe/H] -2). n(CH), n(OH)
come from a background-only LTE EOS WITH molecules (Lightweaver's molecular equilibrium,
DefaultMolecules); that EOS also feeds the passive-atom continua (molecule-depleted, as TS).
The ACTIVE atom's populations are never taken from it.
Not yet included: TS absko adds Barklem's hbop hydrogen (Lyman-alpha wing with its 160 nm
H2 quasi-molecular satellite; up to 90% of TS's 150-200 nm total at [Fe/H] -2) and the NH
continuum after detabs.
"""
from functools import partial
from pathlib import Path
import numpy as np
from lightweaver import constants as c
from lightweaver.atomic_table import PeriodicTable
from lightweaver.LwCompiled import BackgroundProvider, FastBackground

# (element, ion stage 0-based) -> jonabs component label
TS_METAL_COMPONENTS = {
    ('C', 0): 'C_I ', ('C', 1): 'C_II ',
    ('Mg', 0): 'MgI ', ('Mg', 1): 'MgII ',
    ('Al', 0): 'AlI ', ('Al', 1): 'AlII ',
    ('Si', 0): 'SiI ', ('Si', 1): 'SiII ',
    ('Ca', 0): 'CaI ', ('Ca', 1): 'CaII ',
    ('Fe', 0): 'FeI ', ('Fe', 1): 'FeII ',
}


def ts_data_dir():
    from config.constants import require_codex_path
    return require_codex_path('engines.turbospectrum') / 'DATA'


TS_EXTRA_COMPONENTS = {'H+H': '(H+H)', 'H2+': 'H2+ ', 'CH': 'CHcont', 'OH': 'OHcont'}
_MOLECULE_MASS_AMU = {'CH': 13.02, 'OH': 17.01}          # detabs.f values
_AVOGADRO_TS = 6.023e23                                  # detabs.f value


def default_molecules():
    """Lightweaver's shipped molecule set (CO_NLTE variant excluded)."""
    import lightweaver as lw
    d = Path(lw.__file__).parent / 'Data' / 'DefaultMolecules'
    return lw.MolecularTable([str(p) for p in sorted(d.glob('*.molecule')) if 'NLTE' not in p.name])


def background_eos(radset, atmosphere):
    """LTE populations WITH molecules, for the passive background only."""
    return radset.compute_eq_pops(atmosphere, mols=default_molecules())


def ts_metal_factory(species=None, data_dir=None, background_eq=None, extras=True):
    """BackgroundProvider factory for lw.Context(backgroundProvider=...)."""
    return partial(TsMetalBackground, species=species, data_dir=data_dir,
                   background_eq=background_eq, extras=extras)


class TsMetalBackground(BackgroundProvider):
    JONABS = 'jonabs_vac_v19.2.dat'
    IRWIN = 'IRWIN_atoms_v07.3.dat'

    def __init__(self, eqPops, radSet, wavelength, *, species=None, data_dir=None,
                 background_eq=None, extras=True):
        from pipeline.ts_jonabs import JonabsComponent
        super().__init__(eqPops, radSet, wavelength)
        self._args = (eqPops, radSet, wavelength,
                      dict(species=species, data_dir=data_dir, background_eq=background_eq, extras=extras))
        eq = eqPops if background_eq is None else background_eq
        self.base = FastBackground(eq, radSet, wavelength, Nthreads=1)
        self.eq = eq
        self.wave = np.asarray(wavelength, dtype=float)
        d = ts_data_dir() if data_dir is None else data_dir
        self.irwin_path = d / self.IRWIN
        want = TS_METAL_COMPONENTS if species is None else {k: TS_METAL_COMPONENTS[k] for k in species}
        passive = {a.element for a in radSet.passiveAtoms}
        self.replaced, self.added, self.skipped = [], [], []
        self.terms = []
        for (el, ion), label in want.items():
            element = PeriodicTable[el]
            if element not in passive:
                self.skipped.append((el, ion, 'not a passive atom'))
                continue
            model = radSet.atoms[element]
            stages = {l.stage for l in model.levels}
            if ion in stages:
                mode = 'replace'
            elif ion == 0 and 1 in stages:
                mode = 'saha'                                  # e.g. Ca I from CaII_atom
            else:
                self.skipped.append((el, ion, 'stage absent from the model'))
                continue
            comp = JonabsComponent(d / self.JONABS, label)
            self.terms.append((el, ion, element, model, comp, mode))
            (self.replaced if mode == 'replace' else self.added).append((el, ion))
        self.extras = {}
        if extras:
            H = PeriodicTable['H']
            self._H_stage = np.array([l.stage for l in radSet.atoms[H].levels])
            for key, label in TS_EXTRA_COMPONENTS.items():
                if key in _MOLECULE_MASS_AMU and key not in eq.molecularTable:
                    self.skipped.append((key, None, 'molecule not in the background EOS'))
                    continue
                self.extras[key] = JonabsComponent(d / self.JONABS, label)
        self._cache_T = None

    def __deepcopy__(self, memo):
        """Lightweaver's compute_rays deep-copies the Context. The compiled FastBackground built
        on a separate (molecular) EOS cannot be copied, so rebuild from copied inputs."""
        import copy
        e, r, w, kw = copy.deepcopy(self._args, memo)
        out = type(self)(e, r, w, **kw)
        memo[id(self)] = out
        return out

    def _extra(self, T):
        """The TS extra components (H+H, H2+, CH, OH) per (wavelength, depth), 1/m."""
        out = np.zeros((len(self.wave), len(T)))
        if not self.extras:
            return out
        nH = np.asarray(self.eq.atomicPops[PeriodicTable['H']].n)
        nHI = nH[self._H_stage == 0].sum(axis=0) * 1e-6               # cm^-3
        nHII = nH[self._H_stage == 1].sum(axis=0) * 1e-6
        lamA = _air_A(self.wave)
        stim = -np.expm1(-28556. * (5040. / T)[None, :] / lamA[:, None])   # detabs.f
        for key, comp in self.extras.items():
            tab = comp(self.wave, T)
            if key == 'H+H':
                k = tab * (nHI ** 2 * 1e-50)[None, :] * stim
            elif key == 'H2+':
                k = tab * (nHI * nHII * 1e-40)[None, :] * stim
            else:
                n = np.asarray(self.eq[key]) * 1e-6
                k = tab * (n * _MOLECULE_MASS_AMU[key] / _AVOGADRO_TS)[None, :]
            out += k * 100.                                           # 1/cm -> 1/m
        return out

    def _ts_and_old(self, T):
        """Per replaced stage: (TS chi, Lightweaver chi) on the wavelength grid, (nwave, ndep)."""
        from pipeline.ts_jonabs import irwin_ionisation_eV, irwin_partition
        stim = -np.expm1(-c.HC / (c.KBoltzmann * 1e-9 * self.wave[:, None] * T[None, :]))
        new = np.zeros_like(stim); old = np.zeros_like(stim)
        for el, ion, element, model, comp, mode in self.terms:
            pops = self.eq.atomicPops[element]
            n = np.asarray(pops.n)
            stage = np.array([l.stage for l in model.levels])
            if mode == 'saha':
                n_up = n[stage == ion + 1].sum(axis=0)
                U_up = irwin_partition(self.irwin_path, element.Z, ion + 1, T)
                chi_J = irwin_ionisation_eV(self.irwin_path, element.Z, ion) * c.EV
                phi = (c.HPlanck ** 2 / (2 * np.pi * c.MElectron * c.KBoltzmann * T)) ** 1.5
                n_over_U = n_up / (2 * U_up) * self.ne * phi * np.exp(chi_J / (c.KBoltzmann * T))
                new += comp(self.wave, T) * 1e-4 * n_over_U[None, :] * stim
                continue
            n_stage = n[stage == ion].sum(axis=0)
            U = irwin_partition(self.irwin_path, element.Z, ion, T)
            new += comp(self.wave, T) * 1e-4 * (n_stage / U)[None, :] * stim
            for t in model.continua:
                if model.levels[t.i].stage == ion:
                    old += t.alpha(self.wave)[:, None] * stim * n[t.i][None, :]
        return new, old

    def compute_background(self, atmos, chi, eta, sca):
        self.base.compute_background(atmos, chi, eta, sca)
        if not self.terms and not self.extras:
            return
        T = np.asarray(atmos.pyAtmos.temperature, dtype=float)
        self.ne = np.asarray(atmos.pyAtmos.ne, dtype=float)
        new, old = self._ts_and_old(T)
        new = new + self._extra(T)
        planck = (2 * c.HPlanck * c.CLight / (self.wave[:, None] * 1e-9) ** 3
                  / np.expm1(c.HC / (self.wave[:, None] * 1e-9 * c.KBoltzmann * T[None, :])))
        delta = new - old
        chi += delta
        eta += delta * planck
        if not (np.isfinite(chi).all() and np.isfinite(eta).all()) or np.any(chi < sca) or np.any(eta < 0):
            raise ValueError('Invalid TS metal background replacement')


def _air_A(wave_vac_nm):
    """TS works in air above 2000 A (vacuum below), as its tables and STIM do."""
    from pipeline.wavelength_util import vac_to_air
    w = np.asarray(wave_vac_nm, dtype=float) * 10.
    return np.where(w > 2000., vac_to_air(w), w)
