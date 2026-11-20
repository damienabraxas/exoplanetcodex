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
"""
from functools import partial
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


def ts_metal_factory(species=None, data_dir=None):
    """BackgroundProvider factory for lw.Context(backgroundProvider=...)."""
    return partial(TsMetalBackground, species=species, data_dir=data_dir)


class TsMetalBackground(BackgroundProvider):
    JONABS = 'jonabs_vac_v19.2.dat'
    IRWIN = 'IRWIN_atoms_v07.3.dat'

    def __init__(self, eqPops, radSet, wavelength, *, species=None, data_dir=None):
        from pipeline.ts_jonabs import JonabsComponent
        super().__init__(eqPops, radSet, wavelength)
        self.base = FastBackground(eqPops, radSet, wavelength, Nthreads=1)
        self.eq = eqPops
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
        self._cache_T = None

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
        if not self.terms:
            return
        T = np.asarray(atmos.pyAtmos.temperature, dtype=float)
        self.ne = np.asarray(atmos.pyAtmos.ne, dtype=float)
        new, old = self._ts_and_old(T)
        planck = (2 * c.HPlanck * c.CLight / (self.wave[:, None] * 1e-9) ** 3
                  / np.expm1(c.HC / (self.wave[:, None] * 1e-9 * c.KBoltzmann * T[None, :])))
        delta = new - old
        chi += delta
        eta += delta * planck
        if not (np.isfinite(chi).all() and np.isfinite(eta).all()) or np.any(chi < sca) or np.any(eta < 0):
            raise ValueError('Invalid TS metal background replacement')
