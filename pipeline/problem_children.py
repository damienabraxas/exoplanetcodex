#!/usr/bin/env python3
"""
pipeline/problem_children.py
============================
RYA-463 — the MASTER problem-children line registry. One persistent, growing,
cross-star catalog of every line/element/scope that needs special handling
(synthesis / 3D / NLTE / HFS-summing / upper-limit / exclusion / gf-differential),
seeded from solar + the architectural findings, that PREDICTS where each problem
recurs by spectral type / metallicity so every future star's pre-run audit
(RYA-205) walks in knowing the landmines.

Two population sources (RYA-463 §2):
  (a) CURATED — the architectural findings that aren't single-line RYA-458 flags
      (CURATED_SEED below). Settled + merged; the high-value rows.
  (b) AUTO-AGGREGATE — harvest each completed star's RYA-458 *_ew_integrity.csv:
      genuine per-line outliers (BAD_FIT / ABUND_OUTLIER / LIT_DEVIATION) kept
      per-line; routine COG saturation collapsed per (species, ion). Idempotent
      upsert (dedup by key; append the star to observed_in).

This is a CATALOG, not a pipeline stage — it never mutates a measured value and
never re-derives an abundance. The prediction layer is the payoff: predict(Teff,
[Fe/H]) emits the expected problem children for a star before it runs.

Outputs (RYA-463): data/registry/problem_children.csv + docs/science/problem_children_registry.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import data_namespace as ns  # noqa: E402  RYA-469 per-star namespacing

REGISTRY_CSV = ROOT / 'data' / 'registry' / 'problem_children.csv'
REGISTRY_DOC = ROOT / 'docs' / 'science' / 'problem_children_registry.md'

# ── closed vocabularies (RYA-463 §2) ─────────────────────────────────────────
# problem_class for the CURATED layer; the AUTO layer may also use the raw RYA-458
# flag names (BAD_FIT / ABUND_OUTLIER / LIT_DEVIATION) which map in here.
CURATED_CLASSES = {
    'CONTINUUM_LIMITED', 'ATOMIC_BLEND', 'MOLECULAR_BLEND', 'MOLECULAR_SYNTH_ONLY',
    'HFS_SUMMING', 'NLTE_OWED', 'NLTE_VOID', 'BAD_GF', 'SATURATION_COG', 'DATA_GAP',
    'CITED_NOT_MEASURED', 'UPPER_LIMIT',
}
AUTO_CLASSES = {'SATURATION_COG', 'BAD_FIT', 'ABUND_OUTLIER', 'LIT_DEVIATION'}
TREATMENTS = {
    'synthesis', '3D', 'NLTE_grid', 'HFS_sum', 'deblend', 'upper_limit', 'exclude',
    'cited_substitution', 'per_region_source', 'astrophysical_gf_differential', 'none',
}
STATUSES = {'active', 'owed', 'resolved', 'predicted', 'monitoring'}
SEVERITIES = {'critical', 'high', 'medium', 'low'}

SCHEMA_COLUMNS = ['species', 'lambda_or_scope', 'problem_class', 'required_treatment',
                  'observed_in', 'amplifies_with', 'severity', 'governing_tickets',
                  'status', 'population_source', 'notes']

# Solar anchor (the reference point amplification axes are measured against).
SUN_TEFF, SUN_FEH = 5772.0, 0.0
# Curve-of-growth saturation ceiling (RYA-458 single source), and a documented
# rough Teff scaling for the heads-up (NOT a fit — an order-of-magnitude estimate).
REW_LINEAR_CEILING = -4.90
SOLAR_SAT_CEILING_MA = 100.0


# ── CURATED layer (RYA-463 §2(b) seed + the recent settled findings) ─────────
CURATED_SEED = [
    dict(species='O I', lambda_or_scope='[O I] 6300', problem_class='CONTINUUM_LIMITED',
         required_treatment='cited_substitution', observed_in='Sun', amplifies_with='[Fe/H]↑↑',
         severity='high', governing_tickets='447-455,460', status='active',
         notes='Continuum-limited + Ni I blend; O I 777 is PRIMARY, [O I] 6300 a cross-check '
               '(Caffau cited). Gets worse at high [Fe/H] (target risk). RYA-455.'),
    dict(species='C I', lambda_or_scope='5380.34', problem_class='SATURATION_COG',
         required_treatment='exclude', observed_in='Sun', amplifies_with='Teff↑',
         severity='high', governing_tickets='454,458', status='active',
         notes='Saturated (COG) + BAD_FIT (chi2r~103, RYA-454 2.8% continuum strike); EXCLUDED '
               'from the C verdict, other C indicators carry it (RYA-458).'),
    dict(species='Li I', lambda_or_scope='6707.84', problem_class='MOLECULAR_BLEND',
         required_treatment='upper_limit', observed_in='Sun', amplifies_with='Teff↓/[Fe/H]↑',
         severity='medium', governing_tickets='103,458', status='active',
         notes='CN-blended; a clean low value is a RED FLAG. Carried as UPPER_LIMIT, never a '
               'point value (RYA-103/458).'),
    dict(species='Eu II', lambda_or_scope='6645.13', problem_class='HFS_SUMMING',
         required_treatment='HFS_sum', observed_in='Sun', amplifies_with='weak-line',
         severity='medium', governing_tickets='102,458', status='active',
         notes='Weak HFS line; EW must be the summed-component total (RYA-102). Disposition '
               'RECOVERED at 6.80 mA.'),
    dict(species='NH', lambda_or_scope='~3360 (band head)', problem_class='CONTINUUM_LIMITED',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='blue/UV',
         severity='high', governing_tickets='451,460', status='active',
         notes='Blue-edge no-true-continuum (KP cont-SNR ~28) + needs UV molecular synthesis '
               '(TS molecular linelist absent). UNMEASURABLE today, flagged not forced (RYA-460).'),
    dict(species='CN', lambda_or_scope='3883 (violet B-X)', problem_class='MOLECULAR_SYNTH_ONLY',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='Teff↓',
         severity='medium', governing_tickets='369,460', status='active',
         notes='Molecular band, synthesis-only; demoted to N cross-check (RYA-369). Blue-edge on '
               'KP; C/O-coupled.'),
    dict(species='N I', lambda_or_scope='7442-8718 (red multiplets)', problem_class='DATA_GAP',
         required_treatment='per_region_source', observed_in='Sun', amplifies_with='Teff↑',
         severity='high', governing_tickets='369,459,460,491', status='owed',
         notes='NOT physics-owed: the N I NLTE grid is WIRED (RYA-369, GALAH). The gap is DATA-'
               'CHANNEL - the good N I red multiplets (7442/7468, 8216, 8680-8718) are red of the '
               'HARPS 6910 cutoff; HARPS only sees weak C/O-coupled CN red (6000-6200). Solar N '
               'needs the FTS atlas (KP/Delbouille, RYA-459 acquire / 460 wire) or UVES; NH 3360 '
               'is blue of 3780 -> HST/STIS UV (RYA-119). Pending the atlas-channel wired-check '
               '(RYA-491).'),
    dict(species='Cr I', lambda_or_scope='gf pool', problem_class='BAD_GF',
         required_treatment='astrophysical_gf_differential', observed_in='Sun', amplifies_with='all',
         severity='high', governing_tickets='398,399,161', status='active',
         notes='Residual gf-scale floor after the blind cull (Kurucz/ungraded gf); '
               'target-only astrophysical-gf differential owed (RYA-161/162). Canary +0.402.'),
    dict(species='Ti I', lambda_or_scope='gf pool', problem_class='BAD_GF',
         required_treatment='astrophysical_gf_differential', observed_in='Sun', amplifies_with='all',
         severity='high', governing_tickets='398,161', status='active',
         notes='gf-scale residual floor; differential-survey curation owed (RYA-161/162).'),
    dict(species='Ni I', lambda_or_scope='gf pool', problem_class='BAD_GF',
         required_treatment='astrophysical_gf_differential', observed_in='Sun', amplifies_with='all',
         severity='high', governing_tickets='398,161', status='active',
         notes='gf-scale residual floor; differential-survey curation owed (RYA-161/162).'),
    dict(species='Si I', lambda_or_scope='gf pool', problem_class='BAD_GF',
         required_treatment='astrophysical_gf_differential', observed_in='Sun', amplifies_with='all',
         severity='high', governing_tickets='398,399,161', status='active',
         notes='gf-scale residual floor; 3D is not the lever (RYA-399). Curation owed (RYA-161/162).'),
    dict(species='V I', lambda_or_scope='all', problem_class='NLTE_VOID',
         required_treatment='none', observed_in='Sun', amplifies_with='all',
         severity='medium', governing_tickets='404', status='active',
         notes='No NLTE grid exists for V (NLTE_VOID); carried LTE-flagged, no correction '
               'available (RYA-404).'),
    dict(species='Fe I', lambda_or_scope='Procyon gf outliers', problem_class='BAD_GF',
         required_treatment='astrophysical_gf_differential', observed_in='Procyon(pending)',
         amplifies_with='Teff↑', severity='medium', governing_tickets='281', status='predicted',
         notes='F-star bad-gf Fe I outliers + COG-knee shift; predicted from the Teff axis, '
               'confirmed when RYA-281/348/349 land. correct->GES or exclude.'),
    dict(species='CH/CN/C2/NH/OH/CO', lambda_or_scope='bands', problem_class='MOLECULAR_SYNTH_ONLY',
         required_treatment='synthesis', observed_in='all', amplifies_with='Teff↓',
         severity='high', governing_tickets='methodology', status='active',
         notes='Hard carve-out: molecular bands are synthesis-only (no EW). Strengthen toward '
               'cool stars (Teff↓).'),
    dict(species='P I', lambda_or_scope='FUV / near-IR 10581/10596', problem_class='DATA_GAP',
         required_treatment='per_region_source', observed_in='Sun', amplifies_with='—',
         severity='medium', governing_tickets='119,460', status='active',
         notes='FUV lines need HST/STIS (RYA-119); the near-IR multiplet is ground-reachable '
               '(KP, RYA-460) -> DATA_GAP superseded for the Sun, gf-limited (+1.2).'),
    # ── recent settled findings (RYA-460/462) ──
    dict(species='K I', lambda_or_scope='7665/7699', problem_class='NLTE_OWED',
         required_treatment='NLTE_grid', observed_in='Sun', amplifies_with='all',
         severity='low', governing_tickets='460,462', status='resolved',
         notes='RESOLVED: K_Amarsi2020 grid wired (RYA-462), delta -0.312 -> A(K) 5.099 PASS. '
               '7665 sits in the telluric O2 A-band; use 7699.'),
    dict(species='Sc II', lambda_or_scope='4246.82', problem_class='HFS_SUMMING',
         required_treatment='HFS_sum', observed_in='Sun', amplifies_with='blue/UV',
         severity='medium', governing_tickets='460', status='active',
         notes='Blue-edge HFS single line (KP SNR ~180); value close (+0.06) but LOW_CONFIDENCE '
               '-> HFS-resolved synthesis + a cleaner line owed (RYA-460).'),
    dict(species='Co I', lambda_or_scope='3845.46', problem_class='CONTINUUM_LIMITED',
         required_treatment='per_region_source', observed_in='Sun', amplifies_with='blue/UV',
         severity='medium', governing_tickets='460', status='active',
         notes='Blue-edge SNR-limited (KP ~24, chi2r~3100) -> value not trusted; extract cleaner '
               'red Co I lines (within KP 1300nm reach) + HFS (RYA-460).'),
    dict(species='Sr II', lambda_or_scope='4077/4161/4215/4305', problem_class='SATURATION_COG',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='[Fe/H]↑',
         severity='medium', governing_tickets='421,428,429,430,433,520,551', status='resolved',
         notes='Strong near-UV resonance doublet (4077/4215) SAT-culls the naive EW fitter and the '
               'congested near-UV inflates blended Gaussian widths (RYA-429 read 370-633 mA artifacts). '
               'MEASURED via Turbospectrum synthesis with in-window blend handling (RYA-551): '
               'A(Sr II)_NLTE ~2.74 (4077 primary, dEW/dA 203 mA/dex, INSPECT delta -0.004). '
               'required_treatment=synthesis so the RYA-520 raw-EW suppressor routes Sr II to the '
               'synthesis channel (Sr I 6617 EW leg preserved). NLTE grid Bergemann2012/INSPECT covers '
               'ONLY 4077/4215 (delta -0.004/-0.006); 4161/4305 LTE-flagged NLTE_unavailable (corr ~0 '
               'at solar). 4305 EXCLUDED (zero abundance sensitivity, 0.2 mA/dex). Out-of-hull above '
               '[Fe/H]+0.0 -> target risk. Production two-engine freeze pending RYA-527.'),
    # ── RYA-519 Part-B: element-level dispositions for the EW-thin metals (were
    #    silently falling back to stray auto per-line saturation flags) ──
    dict(species='Mn I', lambda_or_scope='6013/16/21 (e6S->z6P triplet)', problem_class='HFS_SUMMING',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='all',
         severity='medium', governing_tickets='411,468,473,476', status='resolved',
         notes='EW SAT-culls the HFS-split triplet (NOT an EW-pool element, RYA-306). MEASURED via '
               'HFS-resolved synthesis (RYA-473: Den Hartog gf+HFS) + triplet-exact Amarsi Mn NLTE '
               'grid (RYA-476) -> A(Mn)_NLTE 5.554 PASS. Grid on Sirius amarsi_galah.'),
    dict(species='Cu I', lambda_or_scope='HFS-split lines', problem_class='HFS_SUMMING',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='all',
         severity='medium', governing_tickets='402,466', status='active',
         notes='EW cannot measure the HFS-split Cu lines (over-saturated core). Route = HFS-'
               'resolved synthesis + RYA-402 b-factor NLTE (RYA-466); not an EW-pool element.'),
    dict(species='Ba II', lambda_or_scope='5853 (+HFS, strong)', problem_class='SATURATION_COG',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='[Fe/H]↑',
         severity='medium', governing_tickets='279,165', status='owed',
         notes='EW SAT-culls the strong Ba II line (5853 REW -4.90, at the knee, +HFS/blend). '
               'ROUTE EXISTS: synthesis of Ba II + the wired Korotin2015 NLTE grid (Ba is the '
               'majority ion). Measurable-owed, NOT a void.'),
    dict(species='Zr II', lambda_or_scope='4629/5350/5372 (strong)', problem_class='SATURATION_COG',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='all',
         severity='low', governing_tickets='279,458', status='owed',
         notes='Zr II is the majority ion (LTE-robust, the Sr II/V II logic); the strong Zr II '
               'lines EW SAT-cull (4629 REW -4.56) -> route = synthesis, no NLTE grid needed. '
               'Measurable-owed, NOT a void.'),
    dict(species='Y II', lambda_or_scope='dominant ion absent from pool', problem_class='DATA_GAP',
         required_treatment='per_region_source', observed_in='Sun', amplifies_with='blue/UV',
         severity='medium', governing_tickets='458', status='active',
         notes='Only 3 Y I lines in HARPS coverage, all saturated/high-err (Y I is trace; Y is '
               'dominantly Y II). No Y II lines in the current pool + no NLTE grid -> UN-MEASURABLE '
               'now; route = acquire blue Y II lines (UV/blue arm).'),
    dict(species='Mg I', lambda_or_scope='b-triplet 5167/72/83 + 5528 (strong)', problem_class='SATURATION_COG',
         required_treatment='synthesis', observed_in='Sun', amplifies_with='all',
         severity='medium', governing_tickets='410,235', status='resolved',
         notes='The strong optical Mg I lines (b-triplet + 5528) EW SAT-cull, so the EW-curation '
               'pool is thin. Value comes from the NLTE-registered path: Amarsi-2020 Mg grid, '
               'RYA-410 Family-A PASS (Na/Mg/Si) -> synthesis + wired grid, not the EW pool. '
               '(Surfaced by the RYA-519 Part-C guard as a routed non-emitter.)'),
]


# ── amplification-axis parsing + star classification ─────────────────────────
def _parse_amplifies(s):
    """Map the human amplifies_with string to a set of machine triggers."""
    s = (s or '').lower()
    trig = set()
    if 'all' in s:
        trig.add('all')
    if '[fe/h]↑↑' in s or 'fe/h]↑↑' in s:
        trig.add('feh_high_high')
    if '[fe/h]↑' in s or 'fe/h]↑' in s:
        trig.add('feh_high')
    if 'teff↑' in s:
        trig.add('teff_high')
    if 'teff↓' in s:
        trig.add('teff_low')
    if 'blue' in s or 'uv' in s:
        trig.add('blue_uv')
    if 'weak' in s:
        trig.add('weak_line')
    return trig


def star_axes(teff, feh):
    """Classify a star's (Teff, [Fe/H]) into the amplification axes it triggers."""
    ax = set(['all'])
    if teff >= 6200:
        ax.add('teff_high')
    if teff <= 5500:
        ax.add('teff_low')
    if feh >= 0.30:
        ax.update(['feh_high', 'feh_high_high'])
    elif feh >= 0.15:
        ax.add('feh_high')
    ax.add('blue_uv')        # any star has the blue region; relevance, severity from SNR
    ax.add('weak_line')
    return ax


def estimated_cog_ceiling_mA(teff):
    """Rough Teff-scaled saturation ceiling for the heads-up (documented estimate,
    NOT a fit): hotter stars have weaker lines, so the saturation onset shifts to
    higher EW. Linear in Teff from the solar single source."""
    return round(SOLAR_SAT_CEILING_MA * (teff / SUN_TEFF), 0)


# These classes are universal special-handling needs — they recur on EVERY star
# regardless of type (the handling never goes away), so the prediction always lists
# them (amplified or not).
_UNIVERSAL_CLASSES = {'MOLECULAR_SYNTH_ONLY', 'NLTE_OWED', 'NLTE_VOID', 'HFS_SUMMING',
                      'CITED_NOT_MEASURED', 'BAD_GF'}


# ── auto-aggregate the RYA-458 EW-integrity flags ────────────────────────────
def _flag_to_class(flag_str):
    """Map a RYA-458 ew_integrity flag-combo to a registry problem_class string."""
    flags = [f for f in str(flag_str).split(',') if f]
    parts = []
    for f in flags:
        if f == 'COG_FLAG':
            parts.append('SATURATION_COG')
        elif f in AUTO_CLASSES:
            parts.append(f)
    return '+'.join(dict.fromkeys(parts)) if parts else 'SATURATION_COG'


def aggregate_ew_integrity(star_id, path=None):
    """Harvest a star's RYA-458 *_ew_integrity.csv into registry rows. Genuine
    per-line outliers (anything beyond a lone COG_FLAG) are kept per-line; routine
    COG saturation is collapsed per (species, ion). Returns a list of row dicts."""
    if path is None:
        path = ns.output_path(star_id, 'ew_integrity.csv', create=False)  # RYA-469 namespaced
    path = Path(path)
    if not path.exists():
        return []
    df = pd.read_csv(path)
    df['ew_integrity'] = df['ew_integrity'].fillna('').astype(str)
    flagged = df[df['ew_integrity'].str.len() > 0].copy()
    rows = []

    # genuine per-line problem children: any flagged line that is NOT a lone COG_FLAG
    genuine = flagged[flagged['ew_integrity'] != 'COG_FLAG']
    for _, r in genuine.iterrows():
        cls = _flag_to_class(r['ew_integrity'])
        treat = ('exclude' if 'BAD_FIT' in cls else
                 'astrophysical_gf_differential' if 'ABUND_OUTLIER' in cls else
                 'deblend' if 'LIT_DEVIATION' in cls else 'exclude')
        rows.append(dict(
            species=f"{r['element']} {r['ion']}",
            lambda_or_scope=f"{float(r['wavelength_air_A']):.3f}",
            problem_class=cls, required_treatment=treat, observed_in=star_id,
            amplifies_with='Teff↑' if 'COG' in cls else 'all',
            severity='high' if 'BAD_FIT' in cls else 'medium',
            governing_tickets='458', status='active', population_source='auto_ew_integrity',
            notes=f"RYA-458 auto: flags=[{r['ew_integrity']}] disp={r.get('ew_disposition','-')} "
                  f"ew={float(r['ew_mA']):.1f} mA"))

    # routine COG saturation: collapse per (species, ion) into one summary row
    cog = flagged[flagged['ew_integrity'] == 'COG_FLAG']
    for (el, ion), grp in cog.groupby(['element', 'ion']):
        rows.append(dict(
            species=f"{el} {ion}", lambda_or_scope=f"strong-line pool ({len(grp)} lines)",
            problem_class='SATURATION_COG', required_treatment='exclude', observed_in=star_id,
            amplifies_with='Teff↓', severity='low', governing_tickets='279,458', status='active',
            population_source='auto_ew_integrity',
            notes=f"RYA-458 auto: {len(grp)} {el} {ion} lines above the REW saturation knee "
                  f"{REW_LINEAR_CEILING} (COG cull); routine, weaker lines used instead."))
    return rows


# ── idempotent upsert ─────────────────────────────────────────────────────────
def _key(row):
    return (str(row['species']).strip(), str(row['lambda_or_scope']).strip(),
            str(row['problem_class']).strip())


def upsert(existing, new_rows):
    """Merge new_rows into the existing registry. Dedup by (species, scope, class);
    on a hit, append the star to observed_in (no dup) and keep the existing row.
    Idempotent: re-running with the same rows is a no-op."""
    rows = (existing.to_dict('records') if isinstance(existing, pd.DataFrame)
            else list(existing or []))
    index = {_key(r): r for r in rows}
    for nr in new_rows:
        nr = {**{c: '' for c in SCHEMA_COLUMNS}, **nr}
        k = _key(nr)
        if k in index:
            cur = index[k]
            seen = [s for s in str(cur.get('observed_in', '')).split(';') if s]
            add = [s for s in str(nr.get('observed_in', '')).split(';') if s]
            for s in add:
                if s not in seen:
                    seen.append(s)
            cur['observed_in'] = ';'.join(seen)
        else:
            index[k] = nr
            rows.append(nr)
    df = pd.DataFrame(rows)
    for c in SCHEMA_COLUMNS:
        if c not in df.columns:
            df[c] = ''
    return df[SCHEMA_COLUMNS]


def build_registry(star_files=None):
    """Assemble the registry: curated seed + auto-aggregated star files.
    star_files: dict {star_id: path-or-None}. Default ingests solar."""
    curated = [{**r, 'population_source': 'curated'} for r in CURATED_SEED]
    df = upsert([], curated)
    if star_files is None:
        star_files = {'solar': None}
    for star_id, path in star_files.items():
        df = upsert(df, aggregate_ew_integrity(star_id, path))
    return df


def load_registry():
    """The persisted registry (built from the curated seed if the CSV is absent)."""
    return (pd.read_csv(REGISTRY_CSV) if REGISTRY_CSV.exists() else build_registry())


# ── per-element disposition accessor (RYA-519: the verdict's disposition source) ─
# The verdict layer asks the registry "what is X's treatment?" so a non-emitting
# expected metal carries a CITED disposition (NLTE_VOID / HFS_SUMMING / BAD_GF /
# DATA_GAP …) instead of silently vanishing. Element-level: aggregates across ions;
# the CURATED architectural row wins over an auto_ew_integrity per-line flag; among
# ties the highest severity wins. Returns None when the registry has no entry
# (Part C then treats an expected-but-unregistered non-emitter as a LOUD failure).
_SEV_RANK = {'critical': 0, 'high': 1, 'medium': 2, 'low': 3}


def disposition_for(species, registry=None):
    """Registry disposition for an element/species, or None if unregistered.

    `species` may be an element ('Mn', 'V') or an ion string ('Mn I', 'V I') —
    matched at the element level (first token). Returns a dict:
    {species, problem_class, required_treatment, status, governing_tickets,
     population_source, notes}."""
    if registry is None:
        registry = load_registry()
    el = str(species).split()[0]
    m = registry[registry['species'].astype(str).str.split().str[0] == el]
    if m.empty:
        return None
    curated = m[m['population_source'] == 'curated'] if 'population_source' in m.columns else m
    pick = curated if not curated.empty else m
    pick = pick.assign(_sev=pick['severity'].map(lambda x: _SEV_RANK.get(x, 9))).sort_values('_sev')
    r = pick.iloc[0]
    return {
        'species': r['species'],
        'problem_class': r['problem_class'],
        'required_treatment': r['required_treatment'],
        'status': r['status'],
        'governing_tickets': r['governing_tickets'],
        'population_source': r.get('population_source', ''),
        'notes': r.get('notes', ''),
    }


# ── prediction layer (the payoff) ─────────────────────────────────────────────
def predict(teff, feh, registry=None, star_name=None):
    """Given a star's (Teff, [Fe/H]), emit the expected problem children — the
    pre-run heads-up. Each entry: {species, scope, class, treatment, status, why}.
    status ∈ {expected, amplified, watch, resolved}."""
    if registry is None:
        registry = (pd.read_csv(REGISTRY_CSV) if REGISTRY_CSV.exists()
                    else build_registry())
    axes = star_axes(teff, feh)
    out = []
    for _, r in registry.iterrows():
        trig = _parse_amplifies(r['amplifies_with'])
        observed = [s for s in str(r['observed_in']).split(';') if s]
        seen_here = star_name is not None and any(star_name.lower() in s.lower() for s in observed)
        amplified = bool(trig & axes - {'all', 'blue_uv', 'weak_line'}) or (
            'all' in trig and (teff >= 6200 or teff <= 5500 or abs(feh) >= 0.15))
        universal = (r['problem_class'] in _UNIVERSAL_CLASSES) or ('all' in trig)
        # metallicity-axis indicators are always worth listing (the indicator choice
        # is in play on every star) — 'watch' when below threshold, 'amplified' above.
        feh_axis = bool(trig & {'feh_high', 'feh_high_high'})
        relevant = seen_here or universal or amplified or feh_axis or bool(trig & axes)
        if not relevant:
            continue
        if str(r['status']) == 'resolved':
            status = 'resolved'
        elif amplified:
            status = 'amplified'
        elif seen_here or universal:
            status = 'expected'
        else:
            status = 'watch'
        out.append({
            'species': r['species'], 'scope': r['lambda_or_scope'],
            'problem_class': r['problem_class'], 'treatment': r['required_treatment'],
            'severity': r['severity'], 'status': status,
            'tickets': r['governing_tickets'],
            'why': _why(r, axes, status)})
    # severity-then-status ordering for the heads-up
    sev_rank = {'critical': 0, 'high': 1, 'medium': 2, 'low': 3}
    out.sort(key=lambda d: (sev_rank.get(d['severity'], 9), d['status']))
    return out


def _why(r, axes, status):
    a = str(r['amplifies_with'])
    if status == 'amplified':
        return f"amplifies with {a} (this star's axis matches)"
    if status == 'resolved':
        return f"resolved ({r['governing_tickets']}) — listed for provenance"
    if status == 'watch':
        return f"could recur (amplifies {a})"
    return f"special handling always needed ({r['required_treatment']})"


def predict_headsup(teff, feh, registry=None, star_name=None):
    """A short human heads-up string for a star (the pre-run audit summary)."""
    preds = predict(teff, feh, registry, star_name)
    axes = star_axes(teff, feh)
    lines = [f"Problem-children heads-up for ({teff:.0f} K, [Fe/H]={feh:+.2f}"
             + (f", {star_name}" if star_name else "") + "):"]
    if 'teff_high' in axes:
        lines.append(f"  - F-star (Teff↑): COG saturation knee shifts up to ~"
                     f"{estimated_cog_ceiling_mA(teff):.0f} mA (est); watch strong-line saturation "
                     f"+ bad-gf Fe outliers (RYA-281).")
    if 'teff_low' in axes:
        lines.append("  - cool (Teff↓): molecular bands (CH/CN/C2/NH/OH/CO) strengthen — "
                     "synthesis-only handling amplified.")
    if 'feh_high_high' in axes:
        lines.append("  - very metal-rich ([Fe/H]↑↑): [O I] 6300 unreliable -> O I 777 primary; "
                     "blends + saturation amplified.")
    elif 'feh_high' in axes:
        lines.append("  - metal-rich ([Fe/H]↑): watch [O I] 6300, Sr II out-of-hull NLTE, "
                     "stronger blends.")
    # flagship O-indicator note (always called out — the RYA-455 decision)
    o_row = next((p for p in preds if p['species'].startswith('O I')), None)
    if o_row:
        if o_row['status'] == 'amplified':
            lines.append("  - [O I] 6300 UNRELIABLE at this metallicity -> use O I 777 (primary).")
        else:
            lines.append("  - [O I] 6300 OK at this [Fe/H] but WATCH it for metal-rich targets "
                         "(amplifies [Fe/H]↑↑) -> O I 777 is the safe primary.")
    n_amp = sum(p['status'] == 'amplified' for p in preds)
    n_exp = sum(p['status'] == 'expected' for p in preds)
    n_watch = sum(p['status'] == 'watch' for p in preds)
    lines.append(f"  -> {len(preds)} entries: {n_amp} amplified, {n_exp} always-on, "
                 f"{n_watch} watch. Top:")
    for p in preds[:6]:
        lines.append(f"     [{p['severity']:6s}|{p['status']:9s}] {p['species']} {p['scope']} "
                     f"-> {p['problem_class']} ({p['treatment']})")
    return '\n'.join(lines)


# ── doc rendering + build ─────────────────────────────────────────────────────
def render_doc(df):
    from datetime import date
    L = ['# RYA-463 — Master problem-children line registry\n',
         f'_Generated {date.today().isoformat()} by pipeline/problem_children.py. A CATALOG '
         '(never mutates a measured value); the cross-star feed for every per-star pre-run '
         'audit (RYA-205) and the line-level companion to the RYA-277 acceptance profiles._\n',
         f'**Rows:** {len(df)} '
         f"({int((df['population_source']=='curated').sum())} curated, "
         f"{int((df['population_source']=='auto_ew_integrity').sum())} auto-aggregated from "
         'RYA-458 EW-integrity).\n',
         '## Schema\n',
         '`' + ' | '.join(SCHEMA_COLUMNS) + '`\n',
         '- **problem_class** (curated): ' + ', '.join(sorted(CURATED_CLASSES)),
         '- **required_treatment**: ' + ', '.join(sorted(TREATMENTS)),
         '- **amplifies_with**: the Teff / [Fe/H] axis where it worsens — the prediction key.\n',
         '## Registry (curated layer)\n',
         '| ' + ' | '.join(SCHEMA_COLUMNS[:-1]) + ' |',
         '|' + '|'.join('---' for _ in SCHEMA_COLUMNS[:-1]) + '|']
    for _, r in df[df['population_source'] == 'curated'].iterrows():
        L.append('| ' + ' | '.join(str(r[c]) for c in SCHEMA_COLUMNS[:-1]) + ' |')
    L.append('\n## Auto-aggregated layer (RYA-458 EW-integrity)\n')
    L.append('| ' + ' | '.join(['species', 'lambda_or_scope', 'problem_class', 'observed_in', 'notes']) + ' |')
    L.append('|---|---|---|---|---|')
    for _, r in df[df['population_source'] == 'auto_ew_integrity'].iterrows():
        L.append(f"| {r['species']} | {r['lambda_or_scope']} | {r['problem_class']} | "
                 f"{r['observed_in']} | {r['notes']} |")
    L.append('\n## Prediction layer — sample heads-ups\n')
    for name, teff, feh in [('Sun', 5777, 0.0), ('Procyon', 6554, 0.01), ('55 Cnc', 5196, 0.32)]:
        L.append('```')
        L.append(predict_headsup(teff, feh, df, star_name=name))
        L.append('```\n')
    L.append('_Predictions for Procyon / 55 Cnc come straight from the curated `amplifies_with` '
             'axes — they work BEFORE those stars run (walk in knowing the landmines). The '
             'auto-layer folds in `procyon_ew_integrity.csv` the moment RYA-281/348/349 land._')
    return '\n'.join(L)


def build_and_write(star_files=None):
    df = build_registry(star_files)
    REGISTRY_CSV.parent.mkdir(parents=True, exist_ok=True)
    REGISTRY_DOC.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(REGISTRY_CSV, index=False)
    REGISTRY_DOC.write_text(render_doc(df))
    return df


if __name__ == '__main__':
    df = build_and_write()
    print(f"Wrote {REGISTRY_CSV.relative_to(ROOT)}  ({len(df)} rows: "
          f"{int((df['population_source']=='curated').sum())} curated + "
          f"{int((df['population_source']=='auto_ew_integrity').sum())} auto)")
    print(f"Wrote {REGISTRY_DOC.relative_to(ROOT)}\n")
    for name, teff, feh in [('Sun', 5777, 0.0), ('Procyon', 6554, 0.01), ('55 Cnc', 5196, 0.32)]:
        print(predict_headsup(teff, feh, df, star_name=name))
        print()
