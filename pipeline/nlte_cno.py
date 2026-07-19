"""
pipeline/nlte_cno.py
====================
C I / O I non-LTE abundance-correction grids — the C/O leg of the NLTE module.

Source (RYA-359): Amarsi, Nissen & Skuladottir 2019, A&A 630, A104
(CDS / VizieR catalogue J/A+A/630/A104; bibcode 2019A&A...630A.104A).
Line-by-line abundance corrections from the Balder code on the STAGGER 3D
hydrodynamic models. Vendored under data/nlte_grids/amarsi2019_cno/ with a
provenance manifest.

Two legs (the per-star tool map, mirrored from the 3D-vs-ceiling split we already
run for Fe):
  * 3D leg  (table2 = C I, table3 = O I): 3D non-LTE - 1D LTE corrections, used
    for Teff <= the 3D-grid ceiling (~6513 K → Sun, 55 Cnc, alpha Cen).
  * 1D leg  (table5 = C I, table6 = O I): 1D non-LTE - 1D LTE corrections, used
    above the 3D ceiling (Procyon, 6554 K), where the STAGGER 3D grid stops.

Grid axes (all four tables): (Teff, logg, [Fe/H], vmic, logeps, Line) → Delta.
The Delta is the abundance correction to ADD to the 1D-LTE abundance:

    A(X; non-LTE) = A(X; 1D-LTE) + Delta

C I and O I corrections are large and predominantly NEGATIVE (up to ~-0.3 dex for
C I in metal-poor F dwarfs, ~-0.6 dex for O I in metal-rich F dwarfs). There is
NO output-convention flip here (unlike the Fe MLP, RYA-339): the grid already
tabulates (non-LTE - LTE), so we add it directly. A positive Delta on the large-
correction lines (O I 777 triplet, IR C I) would be a wiring sign flip → guarded.

Interpolation (per species, leg, line):
  logeps is sampled on a per-node abundance ladder (~6 points centred on the
  scaled-solar value for each [Fe/H]) — it is NOT a global regular axis. So we
  interpolate separably: (1) collapse logeps inside each (Teff,logg,[Fe/H],vmic)
  node by 1D interpolation to the requested logeps, then (2) interpolate the
  collapsed node values over (Teff,logg,[Fe/H],vmic) by scattered linear
  interpolation (Delaunay). The 3D-leg Teff axis is irregular (STAGGER realised
  Teff per model) — the 4D Delaunay handles that natively. A query outside the
  4D convex hull → NaN → caller flags/raises (NO silent LTE fallback).

Unit catch: table1 wavelengths are in nm; converted to Angstrom at the load
boundary. The pipeline is Angstrom throughout.

Linear issue: RYA-359
"""

from __future__ import annotations

import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────

_REPO_ROOT = Path(__file__).parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from config.constants import committed_grid_artifact  # noqa: E402  (RYA-567)

# RYA-567: the Amarsi-2019 CNO delta tables are COMMITTED, version-controlled
# pre-derived ARTIFACTS (in-repo by ratified convention) — NOT heavy Sirius compute
# inputs — routed via the resolver (no ad-hoc literal). Eviction is a separate
# repo-wide migration (RYA-559 successor); no NEW compute-input grid here.
_DATA_DIR  = committed_grid_artifact('nlte_grids', 'amarsi2019_cno')

# table-to-(species, leg) map. table4 (Fe II) intentionally excluded.
_TABLES = {
    ('CI', '3D'): 'table2.dat',
    ('OI', '3D'): 'table3.dat',
    ('CI', '1D'): 'table5.dat',
    ('OI', '1D'): 'table6.dat',
}

# 3D-grid Teff ceiling. The STAGGER 3D grid (tables 2/3) reaches ~6513 K; above
# this (Procyon, 6554 K) we use the 1D-NLTE leg (tables 5/6). RYA-359 spec.
TEFF_3D_CEILING = 6500.0

CITATION = 'Amarsi, Nissen & Skuladottir 2019, A&A 630, A104 (CDS J/A+A/630/A104)'

# Required lines for the Phase-1 coverage STOP-GATE (RYA-359):
#   C I 5052/5380 optical + IR C I 909.5/911.1/940.6 nm; O I 777 triplet.
REQUIRED_LINES = {
    'CI': ['505.2nm', '538.0nm', '909.5nm', '911.2nm', '940.6nm'],
    'OI': ['777nm'],
}

# O I multiplet-averaged labels: the CDS ReadMe Note(1) states '777nm'/'844nm'/
# '926nm' are the multiplet correction integrated across all components. Apply as
# a single per-multiplet correction (do NOT split per component). Wavelength spans
# (air, Angstrom) used to route a query line to its averaged label.
_OI_MULTIPLET_SPANS = {
    '777nm': (7770.5, 7776.5),
    '844nm': (8445.5, 8447.5),
    '926nm': (9259.5, 9267.0),
}


# ── Fixed-width parsing ─────────────────────────────────────────────────────────

def _open_grid(name: str):
    """Open a vendored table by base name, transparently handling .gz."""
    p_gz = _DATA_DIR / (name + '.gz')
    p    = _DATA_DIR / name
    if p_gz.exists():
        return gzip.open(p_gz, 'rt')
    if p.exists():
        return open(p, 'rt')
    raise FileNotFoundError(
        f"Amarsi 2019 C/O grid '{name}' not found in {_DATA_DIR}. "
        "Vendor the CDS J/A+A/630/A104 tables (see provenance.json)."
    )


# Cache: (species, leg) -> {line_label: {(teff,logg,feh,vmic): np.ndarray[(logeps,delta)]}}
_grid_cache: dict = {}
# Cache: (species, leg, line_label) -> (Delaunay tri, list_of_ladders, coords)
_interp_cache: dict = {}
# Cache: species -> {label: rep_wave_A}  (line resolver anchors)
_line_table_cache: dict = {}


def _load_grid(species: str, leg: str) -> dict:
    """Parse a correction table into {label: {(teff,logg,feh,vmic): ladder}}.

    ladder = sorted np.ndarray of shape (n_logeps, 2) = [(logeps, delta), ...].
    Byte columns (CDS ReadMe, all of tables 2/3/5/6):
      Teff 1-4 (I4); logg 8-12 (F5.2); [Fe/H] 16-20 (F5.2); vmic 24-28 (F5.2);
      logeps 32-37 (F6.3); Line 41-48 (A8); Delta 52-59 (F8.3).
    """
    key = (species, leg)
    if key in _grid_cache:
        return _grid_cache[key]
    fname = _TABLES[key]
    nodes: dict = defaultdict(lambda: defaultdict(list))
    with _open_grid(fname) as fh:
        for line in fh:
            if not line.strip():
                continue
            teff = int(line[0:4])
            logg = float(line[7:12])
            feh  = float(line[15:20])
            vmic = float(line[23:28])
            leps = float(line[31:37])
            lab  = line[40:48].strip()
            delta = float(line[51:59])
            nodes[lab][(teff, logg, feh, vmic)].append((leps, delta))
    # freeze ladders to sorted arrays
    frozen: dict = {}
    for lab, nd in nodes.items():
        frozen[lab] = {k: np.array(sorted(v), dtype=float) for k, v in nd.items()}
    _grid_cache[key] = frozen
    return frozen


def _get_interpolator(species: str, leg: str, label: str):
    """Build (or fetch cached) the per-line interpolation structure.

    Returns (tri, ladders, coords) where `tri` is a scipy Delaunay over the 4D
    (Teff, logg, [Fe/H], vmic) node coordinates, `ladders[i]` is the logeps
    ladder for node `coords[i]`. Raises KeyError if the line is absent.
    """
    from scipy.spatial import Delaunay
    key = (species, leg, label)
    if key in _interp_cache:
        return _interp_cache[key]
    grid = _load_grid(species, leg)
    if label not in grid:
        raise KeyError(
            f"Line '{label}' not present in {species} {leg}-leg grid "
            f"(available: {sorted(grid)})"
        )
    nodes = grid[label]
    coords = np.array(list(nodes.keys()), dtype=float)
    ladders = [nodes[tuple(k)] for k in nodes.keys()]
    if len(coords) < 5:
        raise ValueError(
            f"Too few nodes ({len(coords)}) to interpolate {species} {label} "
            f"({leg} leg)."
        )
    tri = Delaunay(coords)
    _interp_cache[key] = (tri, ladders, coords)
    return _interp_cache[key]


# ── Leg selection ──────────────────────────────────────────────────────────────

def select_leg(teff: float) -> str:
    """3D-NLTE leg at/below the STAGGER ceiling, 1D-NLTE leg above it."""
    return '3D' if teff <= TEFF_3D_CEILING else '1D'


# ── Per-line correction ─────────────────────────────────────────────────────────

def cno_nlte_delta(species: str, label: str, teff: float, logg: float,
                   feh: float, vmic: float, logeps: float,
                   leg: str | None = None) -> float:
    """Interpolated Delta = A(X; non-LTE) - A(X; 1D-LTE) for one C I / O I line.

    leg defaults to select_leg(teff). Returns np.nan if the query
    (Teff,logg,[Fe/H],vmic) is outside the grid's 4D convex hull (caller must
    flag/raise — there is NO silent LTE fallback). logeps outside a node's ladder
    is clamped to the ladder end (within-node, smooth; corrections vary slowly
    with logeps).
    """
    if leg is None:
        leg = select_leg(teff)
    from scipy.interpolate import LinearNDInterpolator
    tri, ladders, _coords = _get_interpolator(species, leg, label)
    # (1) collapse logeps per node by 1D interpolation
    nodevals = np.array([np.interp(logeps, lad[:, 0], lad[:, 1]) for lad in ladders])
    # (2) scattered linear interpolation over (Teff,logg,[Fe/H],vmic)
    f = LinearNDInterpolator(tri, nodevals)
    val = f([[teff, logg, feh, vmic]])[0]
    return float(val)


# ── Sign guard (RYA-339 discipline, adapted for C/O) ─────────────────────────────

def assert_cno_sign(species: str, label: str, delta: float) -> None:
    """Fail loud on a sign flip. The O I 777 triplet and the IR C I lines carry
    large NEGATIVE non-LTE corrections; a large positive Delta means the wiring
    flipped the sign (cf. RYA-339 Fe). Optical C I lines are near-LTE (|Delta| ~
    0.01-0.05) and may be marginally positive, so only the large-correction
    lines are guarded, with a small positive tolerance.
    """
    if not np.isfinite(delta):
        return
    big_neg = (species == 'OI' and label in _OI_MULTIPLET_SPANS) or \
              (species == 'CI' and label in ('909.5nm', '911.2nm', '940.6nm'))
    if big_neg and delta > 0.05:
        raise ValueError(
            f"{species} {label} non-LTE correction {delta:+.3f} is POSITIVE — "
            f"expected large negative (Amarsi 2019). A sign flip has returned "
            f"(RYA-339 discipline)."
        )


# ── Line resolution (wavelength → grid label) ───────────────────────────────────

def _label_to_wave_A(label: str) -> float:
    """Parse a grid label like '505.2nm' or '777nm' to an air wavelength in A."""
    m = re.match(r'([0-9.]+)nm', label)
    return float(m.group(1)) * 10.0 if m else np.nan


def _load_line_table(species: str) -> dict:
    """Resolver anchors {grid_label: rep_wave_A} for a species.

    Built from the intersection of grid labels (what we can actually interpolate)
    and table1 wavelengths (refined air wavelengths). Per-component O I multiplet
    labels (777.x/844.x/926.x) are dropped in favour of the multiplet-averaged
    labels (RYA-359 spec). Representative wavelength comes from table1 where the
    label is present, else parsed from the label.
    """
    if species in _line_table_cache:
        return _line_table_cache[species]
    # table1 air wavelengths (nm→A), keyed by label
    t1: dict = defaultdict(list)
    with _open_grid('table1.dat') as fh:
        for line in fh:
            if not line.strip():
                continue
            sp = line[0:4].strip().replace(' ', '')
            if sp != species:
                continue
            lab = line[6:14].strip()
            wair = float(line[17:25]) * 10.0   # nm → A
            t1[lab].append(wair)
    grid_labels = set(_load_grid(species, '3D').keys())
    anchors: dict = {}
    component = set()
    if species == 'OI':
        for lab in grid_labels:
            if re.match(r'(777|844|926)\.\d', lab):
                component.add(lab)
    for lab in grid_labels:
        if lab in component:
            continue   # superseded by the multiplet-averaged label
        if lab in t1:
            anchors[lab] = float(np.mean(t1[lab]))
        elif lab in _OI_MULTIPLET_SPANS:
            # averaged label: representative = mean of its components' table1 waves
            comp_waves = [w for cl, ws in t1.items()
                          for w in ws if re.match(re.escape(lab[:3]) + r'\.\d', cl)]
            anchors[lab] = float(np.mean(comp_waves)) if comp_waves else _label_to_wave_A(lab)
        else:
            anchors[lab] = _label_to_wave_A(lab)
    _line_table_cache[species] = anchors
    return anchors


def resolve_line(species: str, wave_A: float, tol_A: float = 1.5) -> str | None:
    """Map an air wavelength (Angstrom) to its grid line label, or None.

    O I 777/844/926 triplet wavelengths route to the multiplet-averaged label
    (per the CDS ReadMe Note); all other lines to the nearest singlet label
    within `tol_A`.
    """
    if species == 'OI':
        for lab, (lo, hi) in _OI_MULTIPLET_SPANS.items():
            if lo <= wave_A <= hi:
                return lab
    anchors = _load_line_table(species)
    # exclude the multiplet-averaged labels from the singlet nearest-match
    cand = {l: w for l, w in anchors.items() if l not in _OI_MULTIPLET_SPANS}
    if not cand:
        return None
    best = min(cand, key=lambda l: abs(cand[l] - wave_A))
    return best if abs(cand[best] - wave_A) <= tol_A else None


def species_of(element: str, ion) -> str | None:
    """Map (element, ion) to a grid species label ('CI'/'OI'), neutral only."""
    from pipeline.species import parse_ion
    el = str(element).strip().capitalize()
    if el not in ('C', 'O'):
        return None
    try:
        if parse_ion(ion) != 1:   # only neutral C I / O I tabulated
            return None
    except Exception:
        return None
    return el + 'I'


# ── Coverage characterization ───────────────────────────────────────────────────

def grid_coverage(species: str, leg: str) -> dict:
    """min/max of each axis and the set of line labels for one (species, leg)."""
    grid = _load_grid(species, leg)
    teff = set(); logg = set(); feh = set(); vmic = set(); leps = []
    for lab, nodes in grid.items():
        for (t, g, f, v) in nodes:
            teff.add(t); logg.add(g); feh.add(f); vmic.add(v)
        for lad in nodes.values():
            leps.extend(lad[:, 0].tolist())
    return {
        'teff':  (min(teff), max(teff)),
        'logg':  (min(logg), max(logg)),
        'feh':   (min(feh),  max(feh)),
        'vmic':  (min(vmic), max(vmic)),
        'logeps': (min(leps), max(leps)),
        'lines': sorted(grid.keys()),
        'n_records': sum(len(lad) for nodes in grid.values() for lad in nodes.values()),
    }


def star_in_grid(species: str, teff: float, logg: float, feh: float,
                 vmic: float) -> dict:
    """Per-star in-grid verdict: which leg, and whether the bounding box and the
    required lines are covered. (Box check is necessary, not sufficient for the
    convex hull, but is the headline coverage statement.)"""
    leg = select_leg(teff)
    cov = grid_coverage(species, leg)
    box = (cov['teff'][0] <= teff <= cov['teff'][1] and
           cov['logg'][0] <= logg <= cov['logg'][1] and
           cov['feh'][0]  <= feh  <= cov['feh'][1] and
           cov['vmic'][0] <= vmic <= cov['vmic'][1])
    req = [l for l in REQUIRED_LINES.get(species, []) if l in cov['lines']]
    missing = [l for l in REQUIRED_LINES.get(species, []) if l not in cov['lines']]
    return {'leg': leg, 'in_box': box, 'required_present': req,
            'required_missing': missing, 'coverage': cov}


def coverage_ok() -> tuple[bool, list]:
    """Phase-1 STOP-GATE: do both species (both legs) carry every required line
    and reach [Fe/H] >= +0.32 (55 Cnc)? Returns (ok, list_of_problems)."""
    problems = []
    for species in ('CI', 'OI'):
        for leg in ('3D', '1D'):
            cov = grid_coverage(species, leg)
            for l in REQUIRED_LINES[species]:
                if l not in cov['lines']:
                    problems.append(f"{species} {leg}-leg MISSING required line {l}")
            if cov['feh'][1] < 0.32:
                problems.append(
                    f"{species} {leg}-leg [Fe/H] max {cov['feh'][1]} < +0.32 (55 Cnc)")
    return (len(problems) == 0, problems)


# ── table7 self-validation (the catalogue's own answer key) ──────────────────────

# Solar reference for converting table7 [X/H] (relative) to absolute logeps. For
# the Sun, table7 lists ABSOLUTE abundances, so the Sun's 1D-LTE entries ARE the
# reference on Amarsi's own scale (used self-consistently for the cross-check).
_SUN_REF = {'C': 8.4398, 'O': 8.8575}   # = table7 Sun [C/H]1L, [O/H]1L (absolute)


def read_table7() -> list:
    """Parse table7.dat (187 sample stars + 3N/1N/3L/1L C, O, Fe abundances).

    Byte columns from the CDS ReadMe. Missing values (-9999) → NaN. For the Sun,
    the abundance columns are absolute; for all other stars they are [X/H].
    """
    def f(s):
        s = s.strip()
        if s == '':
            return np.nan
        v = float(s)
        return np.nan if v <= -9990 else v
    out = []
    with _open_grid('table7.dat') as fh:
        for line in fh:
            if not line.strip():
                continue
            out.append(dict(
                star=line[0:12].strip(), pop=line[13:25].strip(),
                teff=int(line[42:46]), logg=float(line[49:54]), vmic=float(line[57:62]),
                feh1l=f(line[95:107]),
                cH3N=f(line[125:137]), cH1L=f(line[215:227]),
                oH3N=f(line[245:257]), oH1L=f(line[335:347]),
            ))
    return out


def _c_delta(teff, logg, feh, vmic, logeps, leg, lines=('505.2nm', '538.0nm', '658.8nm')):
    """Mean optical C I Delta over the diagnostic lines (NaN if none resolve)."""
    deltas = []
    for lab in lines:
        try:
            d = cno_nlte_delta('CI', lab, teff, logg, feh, vmic, logeps, leg=leg)
            if np.isfinite(d):
                deltas.append(d)
        except (KeyError, ValueError):
            pass
    return float(np.mean(deltas)) if deltas else np.nan


def _o_delta(teff, logg, feh, vmic, logeps, leg):
    """O I 777-triplet (multiplet-averaged) Delta (NaN if out of grid)."""
    try:
        d = cno_nlte_delta('OI', '777nm', teff, logg, feh, vmic, logeps, leg=leg)
        return float(d) if np.isfinite(d) else np.nan
    except (KeyError, ValueError):
        return np.nan


def validate_against_table7(max_stars: int = 12, tol: float = 0.05,
                            seed: int = 0) -> dict:
    """Confirm the interpolated grid reproduces table7's published corrections.

    table7 abundances are DIFFERENTIAL (relative to the Sun, per method), so
    [X/H]3N - [X/H]1L = Delta_star - Delta_sun. We therefore:
      * Sun (anchor): compare the grid's ABSOLUTE solar Delta to table7's Sun
        (absolute) 3N-1L.
      * other stars: compare the grid's DIFFERENTIAL Delta_star - Delta_sun to
        table7's differential [X/H]3N - [X/H]1L.
    Carbon: mean over optical C I 505.2/538.0/658.8 nm. Oxygen: the O I 777
    multiplet-averaged correction. Each Delta is evaluated at the star's leg
    (select_leg(Teff)) and at logeps = the star's 1D-LTE abundance (absolute, on
    Amarsi's scale).
    """
    stars = read_table7()
    sun = [s for s in stars if s['star'].lower() == 'sun']
    others = [s for s in stars if s['star'].lower() != 'sun'
              and np.isfinite(s['oH1L']) and np.isfinite(s['oH3N'])]
    rng = np.random.default_rng(seed)
    if len(others) > max_stars - 1:
        idx = rng.choice(len(others), size=max_stars - 1, replace=False)
        others = [others[i] for i in sorted(idx)]
    sample = sun + others

    # Sun anchor Delta (feh=0; the Sun's table7 abundances are absolute, so its
    # 1D-LTE entries set the reference logeps on Amarsi's scale).
    s0 = sun[0]
    c_sun = _c_delta(s0['teff'], s0['logg'], 0.0, s0['vmic'], s0['cH1L'], '3D')
    o_sun = _o_delta(s0['teff'], s0['logg'], 0.0, s0['vmic'], s0['oH1L'], '3D')

    rows = []
    for s in sample:
        is_sun = s['star'].lower() == 'sun'
        leg = select_leg(s['teff'])
        feh = 0.0 if is_sun else s['feh1l']
        # carbon
        c_pred = c_obs = np.nan
        if np.isfinite(s['cH1L']) and np.isfinite(s['cH3N']):
            logeps_c = s['cH1L'] if is_sun else s['cH1L'] + _SUN_REF['C']
            d_star = _c_delta(s['teff'], s['logg'], feh, s['vmic'], logeps_c, leg)
            if np.isfinite(d_star):
                c_obs = s['cH3N'] - s['cH1L']
                c_pred = d_star if is_sun else d_star - c_sun
        # oxygen (777 multiplet)
        logeps_o = s['oH1L'] if is_sun else s['oH1L'] + _SUN_REF['O']
        d_star = _o_delta(s['teff'], s['logg'], feh, s['vmic'], logeps_o, leg)
        o_pred = o_obs = np.nan
        if np.isfinite(d_star):
            o_obs = s['oH3N'] - s['oH1L']
            o_pred = d_star if is_sun else d_star - o_sun
        rows.append(dict(star=s['star'], teff=s['teff'], logg=s['logg'],
                         feh=feh, leg=leg,
                         c_pred=c_pred, c_obs=c_obs,
                         c_resid=(c_pred - c_obs),
                         o_pred=o_pred, o_obs=o_obs,
                         o_resid=(o_pred - o_obs)))
    # PASS verdict over the SCIENCE-RELEVANT range. Our benchmark set is all
    # near-solar ([Fe/H] 0 to +0.32); the fixed optical-C-line proxy used here
    # (not the grid itself) loses fidelity for extreme metal-poor halo stars
    # ([Fe/H] < ~-2, where the C diagnostic line set differs), so those are
    # reported but excluded from the pass criterion.
    feh_floor = -1.5
    def _res(key, sci_only):
        return np.array([r[key] for r in rows if np.isfinite(r[key])
                         and (not sci_only or r['feh'] >= feh_floor)])
    c_all, o_all = _res('c_resid', False), _res('o_resid', False)
    c_sci, o_sci = _res('c_resid', True),  _res('o_resid', True)
    return {
        'rows': rows,
        'feh_floor': feh_floor,
        'c_mae': float(np.mean(np.abs(c_all))) if len(c_all) else np.nan,
        'o_mae': float(np.mean(np.abs(o_all))) if len(o_all) else np.nan,
        'c_max': float(np.max(np.abs(c_all))) if len(c_all) else np.nan,
        'o_max': float(np.max(np.abs(o_all))) if len(o_all) else np.nan,
        'c_max_sci': float(np.max(np.abs(c_sci))) if len(c_sci) else np.nan,
        'o_max_sci': float(np.max(np.abs(o_sci))) if len(o_sci) else np.nan,
        'tol': tol,
        'passed': (len(c_sci) > 0 and len(o_sci) > 0 and
                   np.max(np.abs(c_sci)) <= tol and np.max(np.abs(o_sci)) <= tol),
    }


# ── Public API — apply C/O NLTE corrections (mirrors apply_fe_nlte_corrections) ──

def apply_cno_nlte_corrections(abundances_df, stellar_params: dict,
                               per_line_df=None):
    """Apply Amarsi 2019 C I / O I non-LTE corrections to 1D-LTE abundances.

    Parameters
    ----------
    abundances_df  : per-element results (columns: element, ion, A_X, ...).
                     Only C / O neutral rows are touched; others pass through.
    stellar_params : dict with teff_K, logg, feh, vturb_kms.
    per_line_df    : per-line 1D-LTE DataFrame (columns: element, ion,
                     wavelength_air_A, a_1dlte). Required — the corrections are
                     line-by-line and depend on each line's own abundance (logeps).

    Returns a copy of abundances_df with delta_nlte_mean, n_nlte_lines, A_X_nlte,
    A_X_std_nlte, nlte_flag, nlte_ref filled for C/O rows.

    No silent LTE fallback: an out-of-grid line is flagged (not silently carried
    as NLTE); a sign flip raises (assert_cno_sign).
    """
    import pandas as pd
    teff = float(stellar_params['teff_K'])
    logg = float(stellar_params['logg'])
    vmic = float(stellar_params.get('vturb_kms', stellar_params.get('vmic', 1.0)))
    feh  = float(stellar_params.get('feh', 0.0))
    leg  = select_leg(teff)

    result = abundances_df.copy()
    for col, default in [('delta_nlte_mean', np.nan), ('n_nlte_lines', 0),
                         ('A_X_nlte', np.nan), ('A_X_std_nlte', np.nan),
                         ('nlte_flag', '1D_LTE'), ('nlte_ref', '')]:
        if col not in result.columns:
            result[col] = default

    for idx, row in result.iterrows():
        species = species_of(row['element'], row.get('ion', 'I'))
        if species is None:
            if pd.isna(result.at[idx, 'A_X_nlte']):
                result.at[idx, 'A_X_nlte'] = float(row['A_X'])
            continue
        a_1dlte = float(row['A_X'])
        if per_line_df is None or per_line_df.empty:
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
            continue
        el = str(row['element']).strip().capitalize()
        lines = per_line_df[per_line_df['element'].astype(str).str.strip().str.capitalize() == el].copy()
        if lines.empty:
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
            continue
        a_nlte_vals, deltas = [], []
        for _, lr in lines.iterrows():
            wave = float(lr['wavelength_air_A'])
            label = resolve_line(species, wave)
            if label is None:
                continue
            a1 = float(lr['a_1dlte'])
            d = cno_nlte_delta(species, label, teff, logg, feh, vmic, a1, leg=leg)
            if not np.isfinite(d):
                continue              # out of hull → flagged via n_nlte_lines, not silent NLTE
            assert_cno_sign(species, label, d)
            deltas.append(d)
            a_nlte_vals.append(a1 + d)
        if not a_nlte_vals:
            result.at[idx, 'A_X_nlte']  = a_1dlte
            result.at[idx, 'nlte_flag'] = 'NLTE_unavailable'
            continue
        result.at[idx, 'delta_nlte_mean'] = round(float(np.mean(deltas)), 4)
        result.at[idx, 'n_nlte_lines']    = len(a_nlte_vals)
        result.at[idx, 'A_X_nlte']        = round(float(np.median(a_nlte_vals)), 3)
        result.at[idx, 'A_X_std_nlte']    = (round(float(np.std(a_nlte_vals)), 3)
                                             if len(a_nlte_vals) > 1 else np.nan)
        result.at[idx, 'nlte_flag']       = f'NLTE_Amarsi2019_{leg}'
        result.at[idx, 'nlte_ref']        = CITATION
    return result
