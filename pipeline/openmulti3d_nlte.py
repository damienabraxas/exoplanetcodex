"""OpenMULTI3D 1D NLTE departures as an ORCH Engine-B deck — RYA-1371.

Ryan, 2026-11-18: OpenMULTI3D (Bergemann group, pinned in config/reference_engines.yaml) is a
publication-pipeline engine wherever it has a model atom; the Bride fills only the gaps.

This module is the OpenMULTI3D twin of `pipeline.gerber_nlte`. It exposes the SAME surface
the Engine-B fit loop calls (`has_abundance_axis`, `abundance_axis`, `for_node`,
`as_ispec_tuple`, `assert_linelist_supports_nlte`, `assert_depth_match`), so the departures
enter Turbospectrum through the one existing path. Everything else — line list, continuum
rule, fitter, window — is the ORCH's. OpenMULTI3D supplies ONLY b(level, depth).

DECK = departures solved by OpenMULTI3D on the SAME MARCS solar node the Gerber 1D decks use
(`gerber_nlte.MARCS_SOLAR`, 5750/4.5), over an abundance axis, mapped onto that model's 56
layers by column mass (see scripts/build_openmulti3d_deck.py for why not by tau). Paired by
layer index with the route's MARCS.GES atmosphere exactly as the Gerber deck is. Solar node
only, like the Gerber path.

The model atoms are the MPG Keeper files (grids.keeper_atoms_20261118 on Sirius; Mac working
copies in engines.openmulti3d_atoms), pinned by SHA-256 here. Where a Keeper atom is byte-
identical to the Gerber deck's (15 of 17 elements, docs/audit/rya1371/keeper_atom_inventory.md),
the two decks differ ONLY in the statistical-equilibrium solver, which is the comparison
worth having.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class OpenMULTI3DDeckError(RuntimeError):
    pass


def _axis(lo, hi, step=0.1):
    n = int(round((hi - lo) / step))
    return tuple(round(lo + step * i, 3) for i in range(n + 1))


#: key -> deck spec. `Z` and `nk` are checked against the atom header when the deck is built.
DECKS = {
    # Keeper o41f, updated 2026-11-17: "Barklem + Kaulakys data (co-added), closed by Drawin".
    # Identical to o41f_ref (= held atom) for every triplet departure (RYA-1371 COG).
    "O": dict(Z=8, nk=41, atom="atom.o41f",
              atom_sha256="5cc8a84d7ebc6ba23b22d882e4cf0e7a38c4a707bc441e323350236bfed21af5",
              abundance_axis=_axis(8.30, 9.10), convlim=1e-4, seconds=600, threads=4,
              h_collisions="Barklem + Kaulakys (co-added), closed by Drawin"),
    # Keeper o41j (new 2026-11-17): the H-collision systematic, a named variant, never blended.
    "O@o41j": dict(Z=8, nk=41, atom="atom.o41j",
                   atom_sha256="9dc06ba050bfa0db2af144bdddd5d796b18add30aff075bf9030be0eefaf749d",
                   abundance_axis=_axis(8.30, 9.10), convlim=1e-4, seconds=600, threads=4,
                   h_collisions="Belyaev + Kaulakys, closed by Drawin"),
    # RYA-731: Keeper ni538qm, byte-identical to the Gerber deck's (sha256 checked on Sirius
    # 2026-11-21). Ni I only: no Ni II bound-bound transition (nlte_atom_reach/atom.ni538qm.json).
    # Axis covers Gerber's solar-node range 5.1-6.8 and solar 6.2 +/- 1.0.
    "Ni": dict(Z=28, nk=538, atom="atom.ni538qm",
               atom_sha256="dcbebdf513f9a39643604ee2610b99fbfd0a942fdbe051f64d21f8714d3a0872",
               abundance_axis=_axis(5.10, 7.20), convlim=1e-4, seconds=3600, threads=4,
               h_collisions="as in the Keeper atom header (ni538qm)"),
}

_CACHE: dict = {}


def _codex_path(key):
    from config.constants import require_codex_path
    return require_codex_path(key)


def element_of(key: str) -> str:
    return key.split('@')[0]


def deck_spec(key: str) -> dict:
    if key not in DECKS:
        raise OpenMULTI3DDeckError(f'no OpenMULTI3D deck registered for {key!r} (registered: {sorted(DECKS)})')
    from pipeline.gerber_nlte import MARCS_SOLAR
    s = dict(DECKS[key])
    s['abundance_axis'] = list(s['abundance_axis'])
    s['atom_path'] = str(_codex_path('engines.openmulti3d_atoms') / s['atom'])
    s['atmosphere'] = str(MARCS_SOLAR)
    return s


def engine():
    """(app dir, executable, receipt) of the pinned, clean, receipted build."""
    from pipeline.reference_engines import load_engine
    _, app, exe, receipt = load_engine('openmulti3d')
    return app, exe, receipt


def deck_dir(key: str) -> Path:
    s = DECKS[key]
    return _codex_path('grids.openmulti3d_decks') / f"{key.replace('@', '_')}_{s['atom_sha256'][:8]}"


def ensure_deck(key: str) -> Path:
    """The verified deck directory for `key`, building it (minutes) when absent."""
    d = deck_dir(key)
    if not (d / 'deck.json').exists():
        work = Path(os.environ.get('CODEX_WORK_ROOT', Path.home() / 'codex'))
        py = os.environ.get('CODEX_OPENMULTI3D_PY', str(_codex_path('engines.openmulti3d_python') / 'bin' / 'python'))
        app, _, _ = engine()
        env = dict(os.environ, PYTHONPATH=f"{app}:{app / 'HOP'}:{ROOT}", CODEX_WORK_ROOT=str(work))
        r = subprocess.run([py, str(ROOT / 'scripts' / 'build_openmulti3d_deck.py'), '--key', key, '--out', str(d)],
                           env=env, capture_output=True, text=True)
        if r.returncode:
            raise OpenMULTI3DDeckError(f'building the {key} deck failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}')
    verify_deck(key, d)
    return d


def verify_deck(key: str, d: Path) -> dict:
    import hashlib
    meta = json.loads((d / 'deck.json').read_text())
    spec = deck_spec(key)
    _, _, receipt = engine()
    problems = []
    if meta['key'] != key:
        problems.append(f"key {meta['key']}")
    if meta['spec']['atom_sha256'] != spec['atom_sha256']:
        problems.append('atom hash')
    if meta['spec']['abundance_axis'] != spec['abundance_axis']:
        problems.append('abundance axis')
    if meta['engine']['source_commit'] != receipt['source_commit'] or meta['engine']['executable'] != receipt['executable']:
        problems.append('engine build')
    if meta['atmosphere_sha256'] != hashlib.sha256(Path(spec['atmosphere']).read_bytes()).hexdigest():
        problems.append('atmosphere')
    if meta['deck_npz_sha256'] != hashlib.sha256((d / 'deck.npz').read_bytes()).hexdigest():
        problems.append('deck.npz')
    if problems:
        raise OpenMULTI3DDeckError(f'{d}: deck no longer matches its inputs ({", ".join(problems)}); '
                                   f'move it aside and rebuild — never edit a deck in place')
    return meta


def has_abundance_axis(key: str) -> bool:
    deck_spec(key)
    return True


def abundance_axis(key: str) -> tuple:
    return tuple(DECKS[key]['abundance_axis']) if key in DECKS else deck_spec(key)


def _load(key: str):
    if key not in _CACHE:
        d = ensure_deck(key)
        with np.load(d / 'deck.npz') as z:
            _CACHE[key] = dict(axis=z['abundance_axis'].copy(), dep=z['departures'].copy(),
                               lgTau5=z['lgTau5'].copy(), deck_dir=str(d))
    return _CACHE[key]


def for_node(key: str, teff: float, logg: float, feh: float, node: str = 'solar',
             workdir: str | None = None, abundance: float | None = None) -> dict:
    """Departures at `abundance`, linear between the deck's axis nodes. Solar node only."""
    spec = deck_spec(key)
    if not (abs(teff - 5750.0) <= 60.0 and abs(logg - 4.50) <= 0.10 and abs(feh) <= 0.10):
        raise OpenMULTI3DDeckError(
            f'the OpenMULTI3D deck is solved on the solar MARCS node only; asked for '
            f'teff={teff:.0f} logg={logg:.2f} feh={feh:+.2f}. Another star needs its own solve.')
    if abundance is None:
        raise OpenMULTI3DDeckError(f'{key}: the deck resolves an abundance axis; an abundance is required')
    deck = _load(key)
    ax = deck['axis']
    if not (ax[0] - 1e-9 <= abundance <= ax[-1] + 1e-9):
        raise OpenMULTI3DDeckError(f'A = {abundance:.3f} outside the {key} deck axis [{ax[0]:.2f}, {ax[-1]:.2f}]')
    i = int(np.clip(np.searchsorted(ax, abundance) - 1, 0, len(ax) - 2))
    t = (abundance - ax[i]) / (ax[i + 1] - ax[i])
    dep = (1 - t) * deck['dep'][i] + t * deck['dep'][i + 1]
    return dict(abundance=float(abundance), ndep=int(dep.shape[0]), nk=int(dep.shape[1]), tau=deck['lgTau5'],
                departures=dep, corners=[], atom_path=spec['atom_path'], Z=spec['Z'],
                deck_abundance=float(abundance), deck_key=key, provider='openmulti3d', deck_dir=deck['deck_dir'])


class AbundanceDependentDepartures(tuple):
    """iSpec's departure tuple, unchanged in content, that SAYS its coefficients depend on
    the trial abundance. `synth_cache` reads this rather than asking `gerber_nlte` about
    the element: an OpenMULTI3D deck for an element whose Gerber deck is node-fixed (Fe)
    would otherwise be cached as node-fixed, and one for an element with no Gerber deck
    (O) would raise there. iSpec unpacks it as an ordinary tuple."""
    abundance_dependent = True


def as_ispec_tuple(parsed: dict, abundance: float) -> tuple:
    from pipeline import gerber_nlte
    return AbundanceDependentDepartures(gerber_nlte.as_ispec_tuple(parsed, abundance))


def assert_depth_match(parsed: dict, atmosphere) -> None:
    from pipeline import gerber_nlte
    return gerber_nlte.assert_depth_match(parsed, atmosphere)


def assert_linelist_supports_nlte(linelist, Z: int, element: str, wave_lo_A=None, wave_hi_A=None, ion=None) -> int:
    """The Gerber guard, unchanged: it reads the TS line list's NLTE labels (window- and
    stage-local) and refuses rather than let iSpec fall back to LTE silently."""
    from pipeline import gerber_nlte
    return gerber_nlte.assert_linelist_supports_nlte(linelist, Z, element_of(element), wave_lo_A, wave_hi_A, ion=ion)


def provenance(key: str) -> str:
    s = DECKS[key]
    return (f"OpenMULTI3D 1D NLTE departures (Bergemann group, arXiv 2611.09370, pinned build) for "
            f"MPG Keeper {s['atom']} [H collisions: {s['h_collisions']}], solved on the solar MARCS "
            f"node and applied through Turbospectrum on MARCS.GES (RYA-1371)")
