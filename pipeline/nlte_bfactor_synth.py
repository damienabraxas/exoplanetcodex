"""
pipeline/nlte_bfactor_synth.py
==============================
RYA-402 — departure-coefficient (b-factor) NLTE-synthesis capability.

The generalisation of the RYA-396 machinery from "apply a pre-computed delta" to
"COMPUTE the delta from departure coefficients": feed an Amarsi-2020 GALAH
departure-coefficient grid (b per depth per level) into Turbospectrum NLTE
synthesis, synthesise NLTE vs LTE for each line, and bisect to the line-by-line
abundance correction  delta = A(NLTE) - A(LTE).  Unlocks Al/K (RYA-401 follow-up),
complements the RYA-399 3D *delta*-application path (this is the b-factor-only path).

ENGINE FINDING (RYA-402, verified on disk — corrects the ticket premise)
------------------------------------------------------------------------
We do NOT need to build a Turbospectrum NLTE engine from scratch:
  * the compiled `bsyn_lu` ALREADY has the NLTE path (read_departure.f /
    read_nlteinfofile.f; control keys NLTEINFOFILE / MODELATOMFILE / DEPARTUREFILE;
    "no departure file -> coefficients = 1.0", i.e. unity departures reproduce LTE);
  * iSpec ALREADY interpolates departure grids and feeds bsyn
    (ispec.atmospheres.interpolate_nlte_departure_coefficients; departure grids live
    in $ISPEC_DIR/input/dep-grid/{El}_nlte_grid_data.h5).
So the remaining work is (1) a READER for the Amarsi grid, (2) feeding the
departures to Turbospectrum (the Gerber-2022 adaptation, or iSpec's dep-grid), and
(3) the b-factor -> delta EXTRACTION. (1) the PySMEGrid reader is BUILT + verified
against the real vendored Na grid (deep-layer b->1, the format check); (3) is
built + unit-tested; (2)/the live synthesis is wired and FAIL-LOUD until wired to
the Gerber-2022 TS departure machinery.

GRID FORMAT (intake-verified, RYA-402): the grids are PySME "DirectAccess file
Version 1.10" (.grd) binaries — a 4D departure grid over (Teff, logg, [Fe/H],
abund-offset) with a (ndepth, nlevel) b = n_NLTE/n_LTE block per node, plus the
model-atom level data. Read by `PySMEGrid`. Vendored (md5-verified, gitignored;
provenance JSON committed): Cu (v6), S (v7), Na/Al/K (v3).

VALIDATE-DON'T-TUNE (the critical guarantee, ticket Step 3)
----------------------------------------------------------
Before trusting the path for Al/K we reproduce an ALREADY-KNOWN correction: run the
vendored Na grid (INSPECT, delta = -0.107) through this synthesis path and confirm
it lands within tolerance. If it cannot, we STOP (loud) and do RCA — we NEVER fit
to Asplund. `--validate-against Na` is that gate.

DATA (Ryan fetches — see data/nlte_grids/amarsi2020_galah/SOURCE_rya402.json)
The Amarsi-2020 grids are a multi-GB Zenodo deposit; the source is pinned in
SOURCE_rya402.json and posted to the ticket. Na is needed FIRST (the validation).

Linear: RYA-402  (follow-up to RYA-401; complements RYA-399)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

_AMARSI_DIR = _REPO / 'data' / 'nlte_grids' / 'amarsi_galah'

# RYA-402 intake correction (Cu/S, 2026-06-21): the real Amarsi PySME layout is
#   atmos_{El}.txt  = (Teff, logg, [Fe/H]) MARCS NODE table
#   label_{El}.txt  = model-atom ENERGY LEVELS (not the nodes the Step-1 scaffold
#                     assumed) ; nlte_{El}_*.grd = the PySME BINARY departure array.
# Step 2 reworks the reader/driver for this layout via the Gerber-2022 TS adaptation.

# Known-delta anchors for the Step-3 machinery validation (reproduce, never fit).
# Na from the vendored INSPECT grid (RYA-396): the solar median is -0.107, from the
# Na I 3p->4d subordinate doublet 5682.633 (delta -0.103) and 5688.205 (-0.115) at
# (Teff 5772, logg 4.44, [Fe/H] 0). Lower level = Na I 3p 2P (E~2.10 eV), upper = 4d
# 2D (E~4.28 eV) — used to map the line to the grid's model-atom levels.
_KNOWN_DELTA = {
    'Na': {'delta': -0.107, 'tol': 0.03,
           'lines': [(5682.633, 2.102, 4.283), (5688.205, 2.104, 4.284)],  # (wave_A, Elow_eV, Eup_eV)
           'ref': 'Lind et al. 2011 INSPECT (vendored Na_Lind2011_INSPECT.csv, RYA-396)'},
}

# Elements this path is being built to unlock (Amarsi-2020 b-factor only).
TARGET_ELEMENTS = ['Al', 'K']


# ── 1. PySME departure-grid reader (the real .grd format) ────────────────────
# RYA-402 intake (2026-06-21) established the actual on-disk format: the Amarsi
# departure grids ship as a PySME "DirectAccess file Version 1.10" (.grd) — the
# legacy IDL/SME save format. Layout:
#   64-byte version string;
#   header (nblocks u8, dir_length i2, ndir u8);
#   ndir directory entries (key S256, size[23] i4, pointer i8);
#   arrays are memmapped at `pointer`, dtype from the IDL typecode, shape REVERSED
#   (IDL column-major -> NumPy row-major).
# Keys: axis vectors teff/grav/feh/abund; model-atom data conf/term/spec/J/energy;
# and one departure block PER MODEL keyed 't{teff}_g{logg}_m{feh}_a{abund}', each an
# (ndepth, nlevel) float64 array of b = n_NLTE / n_LTE. (Format per the PySME source,
# pysme.nlte.DirectAccessFile; cross-checked live against the vendored Na grid.)

# IDL save-format typecodes -> NumPy dtype.
_IDL_TYPECODE = {1: 'u1', 2: '<i2', 3: '<i4', 4: '<f4', 5: '<f8', 6: '<c8',
                 7: 'S', 9: '<c16', 12: '<u2', 13: '<u4', 14: '<i8', 15: '<u8'}

_MODEL_KEY_RE = re.compile(r't(-?\d+)_g([+-]\d+\.\d+)_m([+-]\d+\.\d+)_a([+-]\d+\.\d+)$')


class PySMEGrid:
    """A vendored Amarsi PySME departure-coefficient grid (.grd DirectAccessFile).

    Read-only, memmapped — never loads the multi-GB file into RAM. `departures()`
    snaps to the nearest (Teff, logg, [Fe/H], abund) node; quadrilinear interpolation
    is layered on top for the live synthesis. Loud on a missing file / unknown key —
    no fake."""

    def __init__(self, path):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(
                f"PySME departure grid not found: {self.path}. Fetch + intake per the "
                f"provenance JSON in {self.path.parent} (RYA-402 Step 1).")
        with open(self.path, 'rb') as f:
            self.version = np.fromfile(f, 'S64', 1)[0].decode('latin1').strip()
            major = int(self.version[26]) if len(self.version) > 26 else 1
            minor = int(self.version[28:30]) if len(self.version) > 29 else 0
            if (major, minor) < (1, 10):
                # v1.00 header uses u2 for nblocks/ndir
                hdt = np.dtype([('nblocks', '<u2'), ('dir_length', '<u2'), ('ndir', '<u2')])
            else:
                hdt = np.dtype([('nblocks', '<u8'), ('dir_length', '<i2'), ('ndir', '<u8')])
            hdr = np.fromfile(f, hdt, 1)[0]
            ddt = np.dtype([('key', 'S256'), ('size', '<i4', 23), ('pointer', '<i8')])
            self._dir = np.fromfile(f, ddt, int(hdr['ndir']))
        self.keys = [k.decode('latin1').strip() for k in self._dir['key']]
        self._index = {k: i for i, k in enumerate(self.keys)}
        nodes, node_keys = [], []
        for k in self.keys:
            m = _MODEL_KEY_RE.match(k)
            if m:
                node_keys.append(k)
                nodes.append([float(m.group(1)), float(m.group(2)),
                              float(m.group(3)), float(m.group(4))])
        self.node_keys = node_keys
        self.nodes = np.asarray(nodes, dtype=float)   # (nmodel, 4): teff,logg,feh,abund

    def get(self, key: str) -> np.ndarray:
        """The array stored under `key` (a real ndarray copy of the memmap)."""
        if key not in self._index:
            raise KeyError(f"key {key!r} not in {self.path.name}")
        i = self._index[key]
        s = self._dir['size'][i]
        ndim = int(s[0])
        shape = tuple(int(x) for x in s[1:1 + ndim])
        dt = _IDL_TYPECODE[int(s[1 + ndim])]
        if dt == 'S':                       # byte-string array: first dim = itemsize
            dt = f'S{shape[0]}'
            shape = shape[1:] if ndim > 1 else (1,)
        mm = np.memmap(self.path, mode='r', offset=int(self._dir['pointer'][i]),
                       dtype=dt, shape=shape[::-1])
        return np.array(mm)

    # axis vectors
    @property
    def teff(self):  return self.get('teff')
    @property
    def logg(self):  return self.get('grav')
    @property
    def feh(self):   return self.get('feh')
    @property
    def abund(self): return self.get('abund')

    @property
    def n_models(self) -> int:
        return len(self.node_keys)

    @property
    def n_levels(self) -> int:
        return int(self.get('energy').shape[-1])

    def level_energy(self):  return self.get('energy')   # eV, (nlevel,)
    def level_J(self):       return self.get('J')

    def nearest_node(self, teff, logg, feh, abund):
        """Index of the nearest model node (Teff scaled so 1000 K ~ 1 dex)."""
        scale = np.array([1000.0, 1.0, 1.0, 1.0])
        q = (np.array([teff, logg, feh, abund]) ) / scale
        d2 = (((self.nodes / scale) - q) ** 2).sum(axis=1)
        return int(np.argmin(d2))

    def departures(self, teff, logg, feh, abund) -> np.ndarray:
        """Departure block b at the nearest node, shape (nlevel, ndepth)."""
        return self.get(self.node_keys[self.nearest_node(teff, logg, feh, abund)])

    def level_for_energy(self, e_eV: float, tol: float = 0.05) -> int:
        """Index of the model-atom level nearest `e_eV` (raises if none within tol)."""
        E = self.get('energy')
        j = int(np.argmin(np.abs(E - e_eV)))
        if abs(float(E[j]) - e_eV) > tol:
            raise ValueError(f"no {self.path.stem} level within {tol} eV of {e_eV} "
                             f"(nearest {float(E[j]):.3f} eV)")
        return j

    def node_tau(self, teff, logg, feh):
        """The departure-grid optical-depth scale at the nearest (Teff,logg,[Fe/H])
        node — `tau` is stored (feh, logg, teff, ndepth) and is abund-independent."""
        tau = self.get('tau')
        it = int(np.argmin(np.abs(self.teff - teff)))
        ig = int(np.argmin(np.abs(self.logg - logg)))
        im = int(np.argmin(np.abs(self.feh - feh)))
        return np.asarray(tau[im, ig, it, :], dtype=float)


def read_amarsi_grid(element: str, base_dir: Path = _AMARSI_DIR) -> PySMEGrid:
    """Load the vendored PySME departure grid for `element` (the nlte_{El}*.grd in
    `base_dir`). Loud FileNotFoundError if not intaken."""
    matches = sorted(base_dir.glob(f'nlte_{element}_*.grd'))
    if not matches:
        raise FileNotFoundError(
            f"No PySME departure grid (nlte_{element}_*.grd) in {base_dir}. Intake it "
            f"first (verify md5, place, provenance) per the RYA-402 intake gate.")
    return PySMEGrid(matches[0])


# ── 3. b-factor -> abundance delta (the extraction; pure, unit-tested) ────────

def delta_from_ews(ew_lte_at: callable, ew_nlte: float,
                   a_lte0: float, a_lo: float = None, a_hi: float = None,
                   tol: float = 1e-3, max_iter: int = 40) -> float:
    """The NLTE abundance correction delta = A(NLTE) - A(LTE).

    Given the NLTE equivalent width `ew_nlte` (synthesised WITH departures) and an
    LTE EW-vs-abundance function `ew_lte_at(A)` (synthesised WITHOUT departures, a
    monotonic increasing function), bisect the LTE abundance A* that reproduces the
    NLTE EW; delta = A* - a_lte0, where a_lte0 is the input (LTE) abundance the NLTE
    synthesis used. This is the standard "what abundance shift mimics NLTE" — the
    same quantity INSPECT/Amarsi tabulate. Pure function: the caller supplies the
    synthesis (real Turbospectrum live, a monotonic stub in tests)."""
    lo = a_lte0 - 1.5 if a_lo is None else a_lo
    hi = a_lte0 + 1.5 if a_hi is None else a_hi
    f_lo, f_hi = ew_lte_at(lo) - ew_nlte, ew_lte_at(hi) - ew_nlte
    if f_lo == 0:
        return lo - a_lte0
    if f_hi == 0:
        return hi - a_lte0
    if np.sign(f_lo) == np.sign(f_hi):
        raise ValueError(
            f"NLTE EW {ew_nlte:.4g} not bracketed by LTE EW in [A={lo},{hi}] "
            f"(EW {ew_lte_at(lo):.4g}..{ew_lte_at(hi):.4g}) — widen the bracket.")
    a = a_lte0
    for _ in range(max_iter):
        a = 0.5 * (lo + hi)
        f = ew_lte_at(a) - ew_nlte
        if abs(f) < tol or (hi - lo) < tol:
            break
        if np.sign(f) == np.sign(f_lo):
            lo, f_lo = a, f
        else:
            hi, f_hi = a, f
    return a - a_lte0


# ── 2/4. Live synthesis driver (wired to the verified engine; fail-loud) ──────

def _ispec_nlte_available() -> bool:
    """True iff iSpec's NLTE departure machinery has grids to work with."""
    try:
        from config.constants import ISPEC_DIR
    except Exception:
        return False
    return (Path(ISPEC_DIR) / 'input' / 'dep-grid').exists()


def _weakline_delta_estimate(grid: PySMEGrid, line, teff, logg, feh) -> float:
    """A FIRST-ORDER departure-only abundance correction for one weak line:
    delta ~= -<log10 b_lower> averaged over the canonical weak-line formation region
    (tau in [0.05, 1.5]). This is NOT the real correction — it omits the line source
    function and the proper contribution-function weighting — and is used ONLY by the
    STOP-gate to demonstrate that the departure shortcut cannot reproduce the anchor."""
    wave, elow, eup = line
    b = grid.departures(teff, logg, feh, 0.0)        # (nlevel, ndepth)
    tau = grid.node_tau(teff, logg, feh)
    jl = grid.level_for_energy(elow)
    m = (tau > 0.05) & (tau < 1.5)
    return float(-np.mean(np.log10(b[jl][m])))


def validate_against(element: str = 'Na') -> dict:
    """Step 3 — the STOP-gate. Try to reproduce a KNOWN delta and refuse to proceed
    (loud) if the available method cannot. NEVER fits to the anchor.

    Current state: the only method available standalone is the departure-only
    shortcut, which DOES NOT reproduce the anchor (Na: ~-0.01 vs -0.107) — exactly as
    expected, because the NLTE EW correction is an RT-weighted source-function
    integral, not a departure average. So the gate FIRES: do not trust the shortcut,
    do not derive Al/K from it. The rigorous path is bsyn NLTE synthesis (see
    synth_ew_nlte_vs_lte), which is the remaining build."""
    if element not in _KNOWN_DELTA:
        raise KeyError(f"No known-delta anchor for {element}; anchors: {list(_KNOWN_DELTA)}")
    known = _KNOWN_DELTA[element]
    grid = read_amarsi_grid(element)                 # raises if not intaken
    ests = [_weakline_delta_estimate(grid, ln, 5772, 4.44, 0.0) for ln in known['lines']]
    shortcut = float(np.median(ests))
    ok = abs(shortcut - known['delta']) <= known['tol']
    if not ok:
        raise RuntimeError(
            f"STOP-gate (RYA-402 Step 3): the departure-only shortcut gives "
            f"delta={shortcut:+.3f} for {element}, but the anchor is "
            f"{known['delta']:+.3f} +/- {known['tol']} [{known['ref']}]. The shortcut "
            f"CANNOT reproduce the known value (as expected — the NLTE EW correction is "
            f"an RT-weighted source-function integral, not a departure average). Do NOT "
            f"trust it and do NOT derive {TARGET_ELEMENTS} from it. RIGOROUS PATH: bsyn "
            f"NLTE synthesis (synth_ew_nlte_vs_lte) with a TS model atom + per-line NLTE "
            f"labels. BLOCKER: the SME grids supply departures + level energies (J, conf, "
            f"term, energy, tau) but NOT the transitions/collisions a TS model atom needs "
            f"(SRC points to an IDL .sav not in the package). Resolve via the Gerber-2022 "
            f"TS-native model atoms or the original model-atom files, then re-run this gate.")
    return {'element': element, 'shortcut_delta': shortcut, 'anchor': known['delta'],
            'within_tol': ok}


def synth_ew_nlte_vs_lte(element: str, line_wave_A: float, stellar_params: dict,
                         grid: 'PySMEGrid' = None) -> tuple:
    """(EW_lte, EW_nlte) for one line via Turbospectrum bsyn NLTE vs LTE. FAIL-LOUD
    until the bsyn NLTE deck is built — no silent LTE fallback (RYA-289 anti-pattern).

    The reader (PySMEGrid) + engine (bsyn_lu has NLTE compiled; iSpec interpolates
    dep grids) are ready. The remaining wiring is the bsyn NLTE deck:
      1. interpolate the PySMEGrid departures b(level, depth) to the target params;
      2. write a TS departure file (on the grid `tau` scale) + a TS MODEL ATOM +
         NLTEINFOFILE (the Gerber-2022 adaptation), or iSpec's dep-grid HDF5;
      3. tag the line's lower/upper levels with NLTE labels (level_for_energy);
      4. run babsma_lu + bsyn_lu twice (with / without departures) -> two EWs.
    BLOCKER (RYA-402): the SME grids carry departures + levels but NOT a TS model
    atom (transitions/collisions). It must come from the Gerber-2022 TS-native atoms
    or the original model-atom files before this can run."""
    raise NotImplementedError(
        f"bsyn NLTE deck for {element} {line_wave_A:.3f} A is not built yet: the SME "
        f"departure grid lacks a TS model atom (transitions/collisions). Source the "
        f"Gerber-2022 TS-native atom (or the original model atom), then wire the bsyn "
        f"NLTE-vs-LTE run here. No silent LTE fallback.")


def _cli(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="RYA-402 b-factor NLTE synthesis")
    p.add_argument('--validate-against', metavar='EL',
                   help="reproduce a known delta (e.g. Na -> INSPECT -0.107)")
    p.add_argument('--elements', help="comma list to derive (e.g. Al,K)")
    p.add_argument('--status', action='store_true', help="report engine + data readiness")
    a = p.parse_args(argv)

    if a.status or (not a.validate_against and not a.elements):
        print("── RYA-402 b-factor NLTE-synthesis capability — readiness ──")
        print(f"  engine: bsyn_lu NLTE compiled = VERIFIED; iSpec dep-grid present = "
              f"{_ispec_nlte_available()}")
        print(f"  Amarsi grids dir: {_AMARSI_DIR}")
        for el in ['Na', 'Cu', 'S'] + TARGET_ELEMENTS:
            try:
                g = read_amarsi_grid(el)
                print(f"    {el}: grid PRESENT ({g.n_models} models, {g.n_levels} levels)")
            except FileNotFoundError:
                print(f"    {el}: grid NOT intaken")
        print("  Reader READY. Na validation + Cu/S/Al/K derivation: pending the "
              "Gerber-2022 TS departure adaptation (Step 2 wiring).")
        return
    if a.validate_against:
        try:
            print(validate_against(a.validate_against))
        except (RuntimeError, FileNotFoundError, NotImplementedError) as exc:
            print(f"BLOCKED (no fake): {exc}")
            raise SystemExit(2)
    if a.elements:
        for el in a.elements.split(','):
            try:
                print(f"  {el}: {synth_ew_nlte_vs_lte(el, 0.0, {})}")
            except (RuntimeError, FileNotFoundError, NotImplementedError) as exc:
                print(f"  {el}: BLOCKED (no fake): {exc}")
                raise SystemExit(2)


if __name__ == '__main__':
    _cli()
