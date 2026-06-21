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
So the remaining work is (1) a READER for the Amarsi ASCII b-factor format, (2) a
CONVERTER to iSpec's dep-grid HDF5 (or a direct bsyn driver), and (3) the
b-factor -> delta EXTRACTION. (1) and (3) are built + unit-tested here; (2)/the live
synthesis run is wired to the verified engine and FAIL-LOUD until the grids are in.

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
# Na from the vendored INSPECT grid (RYA-396): solar delta = -0.107.
_KNOWN_DELTA = {
    'Na': {'delta': -0.107, 'tol': 0.03,
           'ref': 'Lind et al. 2011 INSPECT (vendored Na_Lind2011_INSPECT.csv, RYA-396)'},
}

# Elements this path is being built to unlock (Amarsi-2020 b-factor only).
TARGET_ELEMENTS = ['Al', 'K']


# ── 1. Amarsi b-factor grid reader ───────────────────────────────────────────

class AmarsiBFactorGrid:
    """A parsed Amarsi-2020 GALAH departure-coefficient grid for one element.

    Amarsi ASCII layout (per the pinned SOURCE_rya402.json):
      label_{El}.txt : grid-node table — one row per model, columns
                       (teff, logg, feh, a_x, ...) defining the node order.
      nlte_{El}.txt  : the departure block — for each node, an (nlevel x ndepth)
                       array of b = n_NLTE / n_LTE.
    This reader is format-tolerant (whitespace ASCII) and validates shapes; it does
    NOT invent data — a missing/short file raises.
    """

    def __init__(self, element: str, nodes: np.ndarray, node_cols: list,
                 b: np.ndarray, ndepth: int, nlevel: int):
        self.element = element
        self.nodes = nodes            # (nnode, ncol) float
        self.node_cols = node_cols    # column names
        self.b = b                    # (nnode, nlevel, ndepth) float
        self.ndepth = ndepth
        self.nlevel = nlevel

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    def node_index(self, teff, logg, feh, a_x=None, tol=(30.0, 0.05, 0.05, 0.05)):
        """Nearest grid node to (teff, logg, feh[, a_x]); raises if outside tol."""
        cols = {c: i for i, c in enumerate(self.node_cols)}
        q = [teff, logg, feh] + ([a_x] if a_x is not None else [])
        keys = ['teff', 'logg', 'feh'] + (['a_x'] if a_x is not None else [])
        d2 = np.zeros(self.n_nodes)
        for k, qq, t in zip(keys, q, tol):
            d2 += ((self.nodes[:, cols[k]] - qq) / t) ** 2
        j = int(np.argmin(d2))
        for k, qq, t in zip(keys, q, tol):
            if abs(self.nodes[j, cols[k]] - qq) > t:
                raise ValueError(
                    f"{self.element}: requested {k}={qq} outside grid tol {t} "
                    f"(nearest node {k}={self.nodes[j, cols[k]]}).")
        return j

    def departures(self, teff, logg, feh, a_x=None) -> np.ndarray:
        """b[level, depth] at the nearest node (no interpolation — node-snap; the
        live path will interpolate via iSpec's Delaunay machinery)."""
        return self.b[self.node_index(teff, logg, feh, a_x)]


def _read_label(path: Path):
    """Parse label_{El}.txt -> (nodes ndarray, column names). The Amarsi label file
    is whitespace ASCII; the first non-comment line may name columns. We accept a
    header line beginning '#' with column names, else assume teff logg feh a_x."""
    cols, rows = None, []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            toks = s.lstrip('#').split()
            if any(t.lower() in ('teff', 'logg', 'feh') for t in toks):
                cols = [t.lower() for t in toks]
            continue
        rows.append([float(x) for x in s.split()])
    arr = np.array(rows, dtype=float)
    if cols is None:
        cols = (['teff', 'logg', 'feh', 'a_x'] + [f'c{i}' for i in range(arr.shape[1] - 4)])[:arr.shape[1]]
    return arr, cols


def read_amarsi_grid(element: str, base_dir: Path = _AMARSI_DIR) -> AmarsiBFactorGrid:
    """Load the Amarsi b-factor grid for `element` from `base_dir`. Raises (loud, no
    fake) if the files are absent — they are fetched by Ryan per SOURCE_rya402.json."""
    label = base_dir / f'label_{element}.txt'
    matches = sorted(base_dir.glob(f'nlte_{element}_*.txt')) + \
        ([base_dir / f'nlte_{element}.txt'] if (base_dir / f'nlte_{element}.txt').exists() else [])
    if not label.exists() or not matches:
        raise FileNotFoundError(
            f"Amarsi-2020 b-factor grid for {element} not found in {base_dir} "
            f"(need label_{element}.txt + nlte_{element}*.txt). Fetch from the Zenodo "
            f"deposit pinned in {base_dir / 'SOURCE_rya402.json'} (Step 1).")
    nodes, cols = _read_label(label)
    b, ndepth, nlevel = _read_departures(matches[0], n_nodes=len(nodes))
    return AmarsiBFactorGrid(element, nodes, cols, b, ndepth, nlevel)


def _read_departures(path: Path, n_nodes: int):
    """Parse nlte_{El}.txt -> b[node, level, depth]. The Amarsi block stores, per
    node, an (nlevel x ndepth) departure array; the header line of each node block
    gives 'ndepth nlevel'. Format-tolerant; validates the node count."""
    toks = path.read_text().split()
    i = 0
    blocks, ndepth = [], None
    while i < len(toks):
        nd, nl = int(float(toks[i])), int(float(toks[i + 1]))
        i += 2
        n = nd * nl
        arr = np.array(toks[i:i + n], dtype=float).reshape(nl, nd)
        i += n
        blocks.append(arr)
        ndepth = nd
    b = np.stack(blocks)
    if len(b) != n_nodes:
        raise ValueError(f"departure node count {len(b)} != label node count {n_nodes} "
                         f"in {path.name}")
    return b, ndepth, b.shape[1]


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


def synth_ew_nlte_vs_lte(element: str, line_wave_A: float, stellar_params: dict,
                         grid: 'AmarsiBFactorGrid' = None) -> tuple:
    """(EW_lte, EW_nlte) for one line via Turbospectrum, with and without the
    element's departures. Wired to the verified engine (bsyn_lu NLTE + iSpec
    interpolate_nlte_departure_coefficients). FAIL-LOUD until the departure grid is
    installed — no silent LTE fallback (that is exactly the RYA-289 anti-pattern).

    Resume point: convert the Amarsi grid (read_amarsi_grid) to iSpec's
    input/dep-grid/{El}_nlte_grid_data.h5 (same MARCS nodes), then call the live
    Turbospectrum NLTE synthesis here. The engine is confirmed present; only the
    converted grid is owed (gated on Ryan's fetch)."""
    if not _ispec_nlte_available():
        raise RuntimeError(
            f"NLTE synthesis for {element} {line_wave_A:.3f} A requires the departure "
            f"grid installed at $ISPEC_DIR/input/dep-grid/{element}_nlte_grid_data.h5. "
            f"It is not present. Fetch the Amarsi-2020 {element} grid (SOURCE_rya402.json), "
            f"convert via the Amarsi->iSpec dep-grid step, then re-run. Refusing to fake "
            f"an LTE result (no silent fallback).")
    raise NotImplementedError(
        "live Turbospectrum NLTE EW extraction — wire to "
        "ispec.interpolate_nlte_departure_coefficients + bsyn_lu once a dep-grid is installed")


def validate_against(element: str = 'Na') -> dict:
    """Step 3 — reproduce a KNOWN delta through the synthesis path. Loud STOP on
    mismatch; never fits to Asplund. Returns the verdict dict."""
    if element not in _KNOWN_DELTA:
        raise KeyError(f"No known-delta anchor for {element}; anchors: {list(_KNOWN_DELTA)}")
    known = _KNOWN_DELTA[element]
    grid = read_amarsi_grid(element)             # raises if not fetched
    # ... live synth path here once the dep-grid is installed ...
    raise RuntimeError(
        f"Cannot run the {element} machinery validation yet: the Amarsi-2020 {element} "
        f"departure grid is not installed (fetch + convert per SOURCE_rya402.json). "
        f"Target on success: reproduce delta = {known['delta']:+.3f} +/- {known['tol']} "
        f"[{known['ref']}]. STOP-gate: if the synthesis path cannot reproduce it, do RCA "
        f"before going near {TARGET_ELEMENTS} (never fit the anchor).")


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
        for el in ['Na'] + TARGET_ELEMENTS:
            try:
                g = read_amarsi_grid(el)
                print(f"    {el}: grid PRESENT ({g.n_nodes} nodes)")
            except FileNotFoundError:
                print(f"    {el}: grid NOT fetched yet (Ryan — SOURCE_rya402.json)")
        print("  Na validation + Al/K derivation are BLOCKED on the Amarsi fetch.")
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
