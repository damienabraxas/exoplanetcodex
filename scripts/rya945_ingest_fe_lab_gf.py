#!/usr/bin/env python3
"""
RYA-945 — ingest the primary-laboratory Fe I / Fe II log gf backbone into canonical_gf.csv.

WHAT THIS IS FOR. The 148-line VIS Fe I pool reads 7.586 while the 9-line lab-gf pool reads
7.445 (RYA-819/831). That 0.14 dex is ATOMIC DATA, not measurement — RYA-855 moved 0 of 36
bars because every Fe cell is a mixed pool sitting on the Kurucz floor. RYA-824 proved the
lever is real and that it is bounded by OUR COVERAGE, not by lab availability: the lab
tables hold 250 Fe I VIS lines and our measured pool reached 9 of them. This ticket closes
that gap in the LINE LIST, so the value the deriver reports lands in the literature band by
construction rather than by offset cancellation. It does NOT re-derive anything.

RYA-834 EXTENDED the table redward by APPENDING rows. This is the other operation: it
UPDATES rows that already exist, across every band, under a stated precedence. So the
governing invariant is different and it is enforced explicitly —

    A ROW IS ONLY EVER REWRITTEN BY A STRICTLY BETTER TIER.

    LAB (DH14 / RU14 / BEL17 / DH19)  >  NIST grade C or better  >  Kurucz K##  >  VALD3

`_TIER_RANK` is the single place that ordering is written, `plan_updates` refuses any
downgrade, and `test_fe_lab_gf_ingest_rya945.py` re-derives the invariant from the written
file rather than trusting this script's own report.

🔴 THE TICKET'S PREMISE "Ruffoni 2014 (`RU`) IS PARTIALLY WIRED" IS A REFERENCE-CODE
COLLISION, AND IT MATTERS. `canonical_gf` carries 5,113 rows tagged `RU` — all of them
Fe II, 4200-9199 A. Ruffoni 2014 is an Fe I paper. VALD's own reference block, read from
`VALD/Procyon/RyanSchmitt.019514`, resolves the code: `Fe 2: Raassen & Uylings ... gf:RU`.
`RU` is Raassen & Uylings, an orthogonal-operator CALCULATION for Fe II — not a laboratory
measurement, and not Ruffoni. Ruffoni's Fe I values reach canonical_gf under the GES codes
`GESHRL14` and the truncated bibcode `2014MNRAS.` instead. Treating `RU` as "lab already
wired" would have credited 5,113 theoretical Fe II rows as primary laboratory data.

⚠️ MATCHING IS ON WAVELENGTH *AND* EXCITATION POTENTIAL, NEVER WAVELENGTH ALONE (RYA-780
manufactured a 2.85 dex discrepancy from one 0.06 A coincidence). Every matcher here is run
against a RANDOMISED-WAVELENGTH null so its rate is interpretable, and an AMBIGUOUS line —
more than one canonical row inside the window — is reported and left alone, never guessed.
That is the provenance-trust failure class the ticket names.

⚠️ HFS. `canonical_gf` rows are PHYSICAL lines built by `gf_resolver.cluster_physical_lines`
and carry `hfs_n_components`; the lab tables publish per-transition values. A lab value is
therefore only written onto a row with `hfs_n_components == 1`. Writing a single-component
lab gf onto a summed multi-component row would silently drop the other components' strength.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
LAB_FE1 = ROOT / 'data' / 'reference' / 'fe_gf_lab' / 'fe1_lab_loggf.csv'
LAB_FE2 = ROOT / 'data' / 'reference' / 'fe_gf_lab' / 'fe2_lab_loggf_dh19.csv'
NIST_DIR = ROOT / 'data' / 'linelists' / 'primary_gf'
PROBLEM = ROOT / 'data' / 'registry' / 'problem_children.csv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'rya945_fe_lab_gf'

#: The ticket's short tags, mapped from the vendored tables' own `source` strings. The tag
#: goes in the NEW `lab_source_tag` column; `loggf_reference` keeps the repo's existing
#: `PRIMARY LAB <Source>` wording, which RYA-834's test asserts on and which
#: `scripts/generate_sources_page.py` renders. Two vocabularies, one of them machine-short,
#: is better than renaming a string three other places already parse.
LAB_TAG = {
    'DenHartog2014': 'DH14',
    'Ruffoni2014': 'RU14',
    'Belmonte2017': 'BEL17',
    'DenHartog2019': 'DH19',
    # RYA-1052: Ruffoni's H-band paper. RU14 is their 2014 optical set, so RU13 by the
    # same rule. 🔴 It was missing here, and RYA-1047 consequently wrote the LONG name
    # 'Ruffoni2013' into `lab_source_tag` — a vocabulary this column does not use.
    'Ruffoni2013': 'RU13',
}
LAB_DOI = {
    'DenHartog2014': '10.1088/0067-0049/215/2/23',
    'Ruffoni2014': '10.1093/mnras/stu780',
    'Belmonte2017': '10.3847/1538-4357/aa8cd3',
    'DenHartog2019': '10.3847/1538-4365/ab322e',
    'Ruffoni2013': '10.1088/0004-637X/779/1/17',
}
LAB_VIZIER = {
    'DenHartog2014': 'VizieR J/ApJS/215/23',
    'Ruffoni2014': 'VizieR J/MNRAS/441/3127',
    'Belmonte2017': 'arXiv:1710.07571 source (J/ApJ/848/125 is a 404)',
    'DenHartog2019': 'published PDF (not on VizieR — 0 tables, CDS ftp 404)',
}

#: NIST's accuracy ladder, worst-case percent on A_ki. ⚠️ ENUMERATED IN FULL, '+' TIERS
#: INCLUDED — omitting them once put a B+ line (<=7%) BELOW a B line (<=10%), an inverted
#: ladder that silently demoted the better measurement (RYA-592).
NIST_ACC_PCT = {
    'AAA': 0.3, 'AA': 1.0, 'A+': 2.0, 'A': 3.0, 'B+': 7.0, 'B': 10.0,
    'C+': 18.0, 'C': 25.0, 'D+': 40.0, 'D': 50.0, 'E': 100.0,
}
#: "grade C or better" — the field's cut for solar Fe, <=0.12 dex. C is 25% = 0.097 dex,
#: so the cut is the ladder position AND the dex bound agreeing, not one standing in for
#: the other.
NIST_ADOPT_GRADES = ('AAA', 'AA', 'A+', 'A', 'B+', 'B', 'C+', 'C')
NIST_MAX_SIGMA_DEX = 0.12

TIER_LAB, TIER_NIST, TIER_KURUCZ, TIER_VALD, TIER_OTHER = (
    'LAB', 'NIST-C+', 'KURUCZ', 'VALD3', 'OTHER')
#: The precedence, written ONCE. A lab value is never overwritten by Kurucz or VALD, and
#: `plan_updates` asserts the rank strictly increases on every row it rewrites.
_TIER_RANK = {TIER_VALD: 0, TIER_OTHER: 1, TIER_KURUCZ: 2, TIER_NIST: 3, TIER_LAB: 4}

#: Tolerances. Both keys required, and the wavelength window is MEASURED rather than
#: chosen. Binning the 2,023 graded NIST matches by their wavelength residual and asking
#: how often the two sources then DISAGREE by more than 0.2 dex:
#:
#:     residual (A)   n      frac |delta| > 0.2 dex   median |delta|
#:     0     - 0.001  1783            1.9%              0.0003
#:     0.001 - 0.003   184            3.8%              0.0009
#:     0.003 - 0.006    29            3.4%              0.0006
#:     0.006 - 0.012    17           11.8%              0.0025
#:     0.012 - 0.020    10           30.0%              0.1201
#:
#: The last two bins are a different population: the disagreement rate rises six-fold and
#: the median residual disagreement rises 400-fold. That is the chance-pairing signature
#: RYA-822 diagnosed — coincidences pile up at the tolerance WALL while real matches sit
#: at the fourth decimal. The window is therefore set INSIDE the separation, at the
#: ticket's 5 mA, which the data independently vindicates. It costs 27 of 2,023 graded
#: matches (1.3%) and removes the contaminated tail. `--scan-tol` reprints the scan.
WTOL_A = 0.005
EPTOL_EV = 0.005

#: RYA-945 diagnostic-line criteria that a LINE LIST can carry. EP >= 1.2 eV: lines below
#: it read high (the ticket's cut, and RYA-714's).
DIAG_EP_MIN_EV = 1.2

#: Above this, a rewrite is not a refinement — it is one source contradicting another, and
#: it goes on the outlier ledger by name. 0.20 dex is the blanket Kurucz systematic
#: (RYA-161's K07_SYSTEMATIC_DEX): a move larger than the uncertainty being carried is not
#: absorbed by it.
OUTLIER_DEX = 0.20

#: How close a source's number must be to the row's number for that source's accuracy to
#: DESCRIBE it. 0.02 dex is below the 2-decimal quantisation of the VALD delivery, so it
#: admits "same number, different rounding" and nothing else (pipeline.gf_grades uses the
#: same bound for the same reason).
LOGGF_DESCRIBES_TOL = 0.02

#: 🔴 A BAND WHERE THE TICKET'S PRECEDENCE IS OVERRULED BY A MEASUREMENT.
#: RYA-834 matched NIST ASD Fe I against this band's own line list on wavelength AND EP:
#: 34 of 56 matched, median |delta| 0.0003 dex, max 0.0016. Agreement at the fourth decimal
#: is not corroboration — it is the SAME NUMBER arriving twice, which is RYA-760 (FMW *is*
#: a NIST compilation and VALD copies it). It ratified that here NIST is RECORDED and never
#: ADOPTED, and `test_redward_gf_extension_rya834` enforces it.
#:
#: That ruling is band-specific and evidence-backed, so it outranks this ticket's blanket
#: "NIST A/B/C beats Kurucz". Nothing is lost by honouring it: because the two numbers are
#: already equal, the ROW still collects the NIST grade and its cited sigma through the
#: lateral-metadata path in `apply_plan`. What it does not do is rewrite a VALUE on the
#: authority of a source that is echoing that same value back.
NIST_NO_ADOPT_BANDS = ((9199.9, 12935.0),)

BANDS = (('near-UV', 3000.0, 3780.0), ('VIS', 3780.0, 5500.0),
         ('red', 5500.0, 9199.9), ('NIR', 9199.9, 13000.0))


# ── tiering ──────────────────────────────────────────────────────────────────────────
def classify_tier(row: pd.Series) -> str:
    """The tier a canonical row ALREADY sits at, read from what it actually carries."""
    ref = str(row.get('loggf_reference') or '').strip()
    if ref.startswith('PRIMARY LAB '):
        return TIER_LAB
    grade = str(row.get('nist_grade') or '').strip()
    if grade in NIST_ADOPT_GRADES and ref.startswith('NIST'):
        return TIER_NIST
    if re.fullmatch(r'K\d{2}', ref):
        return TIER_KURUCZ
    if ref.startswith('VALD') or ref == '':
        return TIER_VALD
    return TIER_OTHER


# ── matching ─────────────────────────────────────────────────────────────────────────
def match_unique(src_w: np.ndarray, src_e: np.ndarray,
                 dst_w: np.ndarray, dst_e: np.ndarray,
                 wtol: float, eptol: float) -> np.ndarray:
    """For each SOURCE row, the index of its UNIQUE destination match; -1 none, -2 many."""
    out = np.full(len(src_w), -1, dtype=int)
    order = np.argsort(dst_w)
    dw, de = dst_w[order], dst_e[order]
    for i in range(len(src_w)):
        lo, hi = np.searchsorted(dw, (src_w[i] - wtol, src_w[i] + wtol))
        if hi <= lo:
            continue
        c = np.where(np.abs(de[lo:hi] - src_e[i]) <= eptol)[0]
        out[i] = int(order[lo + c[0]]) if c.size == 1 else (-2 if c.size > 1 else -1)
    return out


def null_rate(src_w, src_e, dst_w, dst_e, wtol, eptol, trials, seed) -> float:
    """The SAME matcher on randomised source wavelengths, EP kept real."""
    rng = np.random.default_rng(seed)
    lo, hi = float(dst_w.min()), float(dst_w.max())
    rates = []
    for _ in range(trials):
        w = rng.uniform(lo, hi, len(src_w))
        idx = match_unique(w, src_e, dst_w, dst_e, wtol, eptol)
        rates.append(float((idx >= 0).mean()))
    return float(np.mean(rates))


# ── the lab pass ─────────────────────────────────────────────────────────────────────
def match_lab(canon: pd.DataFrame, lab: pd.DataFrame, species: str,
              wtol: float, eptol: float) -> pd.DataFrame:
    """One row per LAB line, carrying the canonical row it identifies (or why not).

    Runs source-by-source so that a line measured by two labs is visible as such: the
    ticket asks for the inter-source spread as a diagnostic, and a single pooled match
    would hide it.
    """
    sub = canon[canon['_species'] == species]
    dst_w = sub.wavelength_air_A.to_numpy(float)
    dst_e = sub.excitation_potential_eV.to_numpy(float)
    dst_i = sub.index.to_numpy()

    recs = []
    for src, g in lab.groupby('source'):
        g = g.reset_index(drop=True)
        idx = match_unique(g.wavelength_air_A.to_numpy(float), g.elo_eV.to_numpy(float),
                           dst_w, dst_e, wtol, eptol)
        for k, j in enumerate(idx):
            recs.append({
                'lab_source': src, 'lab_tag': LAB_TAG[src], 'species': species,
                'lab_wavelength_A': float(g.wavelength_air_A[k]),
                'lab_elo_eV': float(g.elo_eV[k]),
                'lab_loggf': float(g.loggf[k]),
                'lab_sigma_dex': float(g.e_loggf_dex[k]),
                'canon_index': int(dst_i[j]) if j >= 0 else -1,
                'status': 'matched' if j >= 0 else ('ambiguous' if j == -2 else 'absent'),
            })
    m = pd.DataFrame(recs)
    hit = m.canon_index >= 0
    m['canon_wavelength_A'] = np.nan
    m['canon_ep_eV'] = np.nan
    m['canon_log_gf'] = np.nan
    m['canon_reference'] = None
    m['canon_hfs_n'] = np.nan
    if hit.any():
        ci = m.loc[hit, 'canon_index'].to_numpy()
        m.loc[hit, 'canon_wavelength_A'] = canon.loc[ci, 'wavelength_air_A'].to_numpy()
        m.loc[hit, 'canon_ep_eV'] = canon.loc[ci, 'excitation_potential_eV'].to_numpy()
        m.loc[hit, 'canon_log_gf'] = canon.loc[ci, 'log_gf'].to_numpy()
        m.loc[hit, 'canon_reference'] = canon.loc[ci, 'loggf_reference'].to_numpy()
        m.loc[hit, 'canon_hfs_n'] = pd.to_numeric(
            canon.loc[ci, 'hfs_n_components'], errors='coerce').to_numpy()
    m['delta_lab_minus_canon'] = m.lab_loggf - m.canon_log_gf
    return m


def resolve_lab_precedence(m: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One winning lab value per canonical row, plus the inter-source spread.

    ⚠️ HFS: a multi-component canonical row carries the SUM of its components' gf. A lab
    table publishes one component. Writing one onto the other loses the rest, so those
    rows are refused here and reported rather than silently written.
    """
    hit = m[(m.canon_index >= 0)].copy()
    hfs = hit[hit.canon_hfs_n > 1].copy()
    hit = hit[~(hit.canon_hfs_n > 1)]

    spread = (hit.groupby('canon_index')
                 .agg(n_sources=('lab_source', 'nunique'),
                      loggf_min=('lab_loggf', 'min'), loggf_max=('lab_loggf', 'max'),
                      sources=('lab_tag', lambda s: '+'.join(sorted(set(s)))))
                 .reset_index())
    spread['inter_source_spread_dex'] = spread.loggf_max - spread.loggf_min
    spread = spread[spread.n_sources > 1]
    resolve_lab_precedence.last_spread = spread

    # Precedence WITHIN the lab tier is the smallest published sigma, then the most recent
    # measurement. Both are properties of the measurement, not of our preference.
    year = hit.lab_source.str.extract(r'(\d{4})$')[0].astype(int)
    hit = hit.assign(_year=year).sort_values(
        ['canon_index', 'lab_sigma_dex', '_year'], ascending=[True, True, False])
    win = hit.drop_duplicates('canon_index', keep='first').drop(columns='_year')
    return win, hfs


# ── the NIST pass ────────────────────────────────────────────────────────────────────
def load_nist(species: str) -> pd.DataFrame:
    slug = species.replace(' ', '')
    files = sorted(NIST_DIR.glob(f'nist_asd_{slug}_*.tsv'))
    if not files:
        raise SystemExit(f"no NIST ASD pull for {species} in {NIST_DIR} — run "
                         f"scripts/rya822_pull_nist_nearuv.py --species '{species}'")
    d = pd.concat([pd.read_csv(f, sep='\t') for f in files], ignore_index=True)
    d['nist_grade'] = d.nist_grade.astype(str).str.strip().replace(
        {'nan': None, 'None': None, '': None})
    d['nist_sigma_dex'] = d.nist_grade.map(
        lambda g: np.log10(1 + NIST_ACC_PCT[g] / 100.0) if g in NIST_ACC_PCT else np.nan)
    d = d[d.wavelength_A.notna() & d.ei_eV.notna()].copy()
    # The pulls overlap at their shared edges; the same ASD line must not enter twice.
    d = d.drop_duplicates(subset=['wavelength_A', 'ei_eV']).reset_index(drop=True)
    d.attrs['files'] = [f.name for f in files]
    return d


def match_nist(canon: pd.DataFrame, nist: pd.DataFrame, species: str,
               wtol: float, eptol: float) -> pd.DataFrame:
    sub = canon[canon['_species'] == species]
    idx = match_unique(nist.wavelength_A.to_numpy(float), nist.ei_eV.to_numpy(float),
                       sub.wavelength_air_A.to_numpy(float),
                       sub.excitation_potential_eV.to_numpy(float), wtol, eptol)
    dst_i = sub.index.to_numpy()
    out = nist.copy()
    out['canon_index'] = [int(dst_i[j]) if j >= 0 else j for j in idx]
    out['species'] = species
    hit = out.canon_index >= 0
    out['canon_log_gf'] = np.nan
    out['canon_hfs_n'] = np.nan
    if hit.any():
        ci = out.loc[hit, 'canon_index'].to_numpy()
        out.loc[hit, 'canon_log_gf'] = canon.loc[ci, 'log_gf'].to_numpy()
        out.loc[hit, 'canon_hfs_n'] = pd.to_numeric(
            canon.loc[ci, 'hfs_n_components'], errors='coerce').to_numpy()
    out['delta_nist_minus_canon'] = out.log_gf - out.canon_log_gf
    # A canonical row reached by TWO distinct ASD lines is not identified by either.
    dup = out[hit].canon_index.value_counts()
    out['canon_shared_with'] = out.canon_index.map(dup).fillna(0).astype(int)
    return out


# ── the plan ─────────────────────────────────────────────────────────────────────────
def plan_updates(canon: pd.DataFrame, lab_win: pd.DataFrame,
                 nist: pd.DataFrame) -> pd.DataFrame:
    """Every rewrite this run would make, with the tier it moves FROM and TO."""
    rows = []
    for r in lab_win.itertuples():
        rows.append({
            'canon_index': r.canon_index, 'to_tier': TIER_LAB,
            'log_gf': r.lab_loggf, 'sigma_dex': r.lab_sigma_dex,
            'loggf_reference': f'PRIMARY LAB {r.lab_source}',
            'lab_source_tag': r.lab_tag, 'nist_grade': None,
            'gf_source_doi': LAB_DOI[r.lab_source],
            'adjudication_status': 'lab_rya945',
            'delta': r.lab_loggf - r.canon_log_gf,
        })
    claimed = {r['canon_index'] for r in rows}
    ok = nist[(nist.canon_index >= 0) & (nist.canon_shared_with == 1)
              & nist.log_gf.notna() & nist.nist_grade.isin(NIST_ADOPT_GRADES)
              & (nist.nist_sigma_dex <= NIST_MAX_SIGMA_DEX)
              & ~(nist.canon_hfs_n > 1)]
    # Best (smallest sigma) ASD row per canonical line, then never over a lab value.
    ok = ok.sort_values('nist_sigma_dex').drop_duplicates('canon_index', keep='first')
    no_adopt = 0
    for r in ok.itertuples():
        if r.canon_index in claimed:
            continue
        w = float(canon.at[r.canon_index, 'wavelength_air_A'])
        if any(lo < w <= hi for lo, hi in NIST_NO_ADOPT_BANDS):
            no_adopt += 1
            rows.append({
                'canon_index': r.canon_index,
                'to_tier': canon.at[r.canon_index, '_tier'],   # lateral by construction
                'log_gf': r.log_gf, 'sigma_dex': r.nist_sigma_dex,
                'loggf_reference': canon.at[r.canon_index, 'loggf_reference'],
                'lab_source_tag': 'NIST-C+ (recorded, not adopted: RYA-834)',
                'nist_grade': r.nist_grade,
                'gf_source_doi': 'https://physics.nist.gov/asd',
                'adjudication_status': canon.at[r.canon_index, 'adjudication_status'],
                'delta': r.log_gf - r.canon_log_gf,
            })
            continue
        rows.append({
            'canon_index': r.canon_index, 'to_tier': TIER_NIST,
            'log_gf': r.log_gf, 'sigma_dex': r.nist_sigma_dex,
            'loggf_reference': f'NIST-C+ {r.ref_transition_probability}'.strip(),
            'lab_source_tag': 'NIST-C+', 'nist_grade': r.nist_grade,
            'gf_source_doi': 'https://physics.nist.gov/asd',
            'adjudication_status': 'nist_rya945',
            'delta': r.log_gf - r.canon_log_gf,
        })
    plan_updates.n_nist_recorded_not_adopted = no_adopt
    plan = pd.DataFrame(rows)
    if plan.empty:
        return plan
    plan['from_tier'] = canon.loc[plan.canon_index, '_tier'].to_numpy()

    # 🔴 THE INVARIANT. A rewrite must strictly IMPROVE the tier. Anything else is a
    # downgrade — the failure the ticket names as "VALD/Kurucz never override lab" — and
    # it is dropped here rather than written and apologised for later.
    rank_from = plan.from_tier.map(_TIER_RANK)
    rank_to = plan.to_tier.map(_TIER_RANK)
    plan['kept'] = rank_to > rank_from
    return plan


def apply_plan(canon: pd.DataFrame, plan: pd.DataFrame) -> pd.DataFrame:
    out = canon.copy()
    keep = plan[plan.kept]
    ci = keep.canon_index.to_numpy()
    out.loc[ci, 'log_gf'] = keep.log_gf.to_numpy().round(4)
    out.loc[ci, 'loggf_reference'] = keep.loggf_reference.to_numpy()
    out.loc[ci, 'adjudication_status'] = keep.adjudication_status.to_numpy()
    ng = keep.nist_grade.to_numpy()
    out.loc[ci[pd.notna(ng)], 'nist_grade'] = ng[pd.notna(ng)]
    out.loc[ci, 'lab_source_tag'] = keep.lab_source_tag.to_numpy()
    out.loc[ci, 'gf_sigma_dex'] = keep.sigma_dex.to_numpy().round(4)
    out.loc[ci, 'gf_source_doi'] = keep.gf_source_doi.to_numpy()
    out.loc[ci, 'gf_tier'] = keep.to_tier.to_numpy()

    # A rewrite REFUSED as a lateral move still leaves the row without a cited sigma, and
    # a row that already carried the right VALUE with no accuracy attached is exactly the
    # gap RYA-824 measured: "the Kurucz floor is in the sigma we assign, not uniformly in
    # the values". So the refused laterals get their metadata even though their number
    # does not move -- and ONLY where the source's number IS the row's number, because a
    # grade must describe the value actually used (RYA-850).
    lat = plan[~plan.kept & (plan.delta.abs() <= LOGGF_DESCRIBES_TOL)]
    if len(lat):
        li = lat.canon_index.to_numpy()
        out.loc[li, 'gf_sigma_dex'] = lat.sigma_dex.to_numpy().round(4)
        out.loc[li, 'lab_source_tag'] = lat.lab_source_tag.to_numpy()
        out.loc[li, 'gf_source_doi'] = lat.gf_source_doi.to_numpy()
    return out


# ── the diagnostic flag ──────────────────────────────────────────────────────────────
def excluded_lines() -> dict[str, set[float]]:
    """Line-level exclusions from the problem-children registry (RYA-463).

    ⚠️ THIS IS NOT THE WHOLE OF THE TICKET'S "not blend-flagged". `blend_flag` is a column
    of the MEASURED EW table, set per observation by the fitter, and applied downstream by
    `abundances_derive._prefilter`. A line list cannot carry it. What a line list CAN carry
    is the registry's standing per-line exclusions, so those are what this flag uses, and
    the measurement-time filter still runs after it. Stated rather than quietly narrowed.
    """
    d = pd.read_csv(PROBLEM)
    out: dict[str, set[float]] = {}
    for r in d.itertuples():
        if str(r.required_treatment).strip() != 'exclude':
            continue
        sp = str(r.species).strip()
        for tok in re.findall(r'\d{4,5}\.\d+', str(r.lambda_or_scope)):
            out.setdefault(sp, set()).add(float(tok))
    return out


def flag_diagnostic(canon: pd.DataFrame) -> pd.DataFrame:
    ex = excluded_lines()
    tier_ok = canon['gf_tier'].isin([TIER_LAB, TIER_NIST])
    ep_ok = pd.to_numeric(canon.excitation_potential_eV, errors='coerce') >= DIAG_EP_MIN_EV
    is_fe = canon['_species'].isin(['Fe I', 'Fe II'])
    excl = pd.Series(False, index=canon.index)
    for sp, waves in ex.items():
        m = canon['_species'] == sp
        for w in waves:
            excl |= m & (canon.wavelength_air_A.sub(w).abs() <= 0.05)
    canon['is_diagnostic'] = (is_fe & tier_ok & ep_ok & ~excl)
    canon['diagnostic_excluded_by_registry'] = (is_fe & tier_ok & ep_ok & excl)
    return canon


def demote_unverified_nist(canon: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Strip the NIST-C+ tier from rows this run could not stand behind.

    🔴 THIS TEST FOUND TWO PRE-EXISTING DEFECTS AND THE FIX IS TO REFUSE THE CLAIM, NOT TO
    INVENT THE MISSING NUMBER.

    (a) FIVE near-UV rows carry a NIST value on a canonical row with
        `hfs_n_components == 2` (RYA-822). A single-transition ASD gf cannot describe a
        gf-summed two-component cluster — it is the same HFS error this script refuses to
        make on the way in, arriving from a previous ticket.
    (b) FOUR rows tagged `NIST ASD v5.11 grade A/B` (RYA-354) carry a grade with no
        verified sigma. RYA-852 already showed two of them are STALE: 6149.246 and
        6247.557 are stored as B while live ASD reports E (0.301) and D (0.176) —
        understated 7.3x and 4.3x. Re-deriving a sigma from the stored letter would
        launder that staleness into a number.

    Both are demoted to OTHER, which keeps the value and the note and only withdraws the
    claim that the value is graded. They are written out by name so the withdrawal is
    auditable rather than a silent count change.
    """
    is_fe = canon['_species'].isin(['Fe I', 'Fe II'])
    hfs = pd.to_numeric(canon.hfs_n_components, errors='coerce')
    bad = (is_fe & (canon.gf_tier == TIER_NIST)
           & (canon.gf_sigma_dex.isna() | (hfs > 1)))
    out = canon.loc[bad, ['species', 'wavelength_air_A', 'excitation_potential_eV',
                          'log_gf', 'loggf_reference', 'nist_grade',
                          'hfs_n_components', 'adjudication_status']].copy()
    out['demoted_because'] = np.where(
        hfs[bad] > 1, 'multi-component cluster: a per-transition ASD gf cannot describe '
                      'a gf-summed row', 'grade recorded with no sigma this run verified')
    canon.loc[bad, 'gf_tier'] = TIER_OTHER
    return canon, out


def band_of(w: float) -> str:
    for name, lo, hi in BANDS:
        if lo <= w < hi:
            return name
    return 'out-of-band'


# ── report ───────────────────────────────────────────────────────────────────────────
def coverage(canon: pd.DataFrame) -> pd.DataFrame:
    d = canon[canon['_species'].isin(['Fe I', 'Fe II'])].copy()
    d['band'] = d.wavelength_air_A.map(band_of)
    tab = (d.groupby(['_species', 'band', 'gf_tier']).size()
             .unstack('gf_tier').fillna(0).astype(int))
    for c in (TIER_LAB, TIER_NIST, TIER_KURUCZ, TIER_VALD, TIER_OTHER):
        if c not in tab.columns:
            tab[c] = 0
    tab = tab[[TIER_LAB, TIER_NIST, TIER_KURUCZ, TIER_VALD, TIER_OTHER]]
    tab['total'] = tab.sum(axis=1)
    diag = d[d.is_diagnostic].groupby(['_species', 'band']).size()
    tab['diagnostic'] = diag.reindex(tab.index).fillna(0).astype(int)
    return tab.reindex([(s, b) for s in ('Fe I', 'Fe II') for _, b, *_ in
                        [(0, n, 0, 0) for n, _, _ in BANDS]]).dropna(how='all').astype(int)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='write canonical_gf.csv')
    ap.add_argument('--scan-tol', action='store_true')
    ap.add_argument('--null-trials', type=int, default=10)
    ap.add_argument('--wtol-A', type=float, default=WTOL_A)
    ap.add_argument('--eptol-eV', type=float, default=EPTOL_EV)
    a = ap.parse_args()

    canon = pd.read_csv(CANON, low_memory=False)
    canon['_species'] = canon.species.astype(str).str.strip()
    canon['_tier'] = canon.apply(classify_tier, axis=1)
    for col, default in (('lab_source_tag', None), ('gf_sigma_dex', np.nan),
                         ('gf_source_doi', None), ('is_diagnostic', False)):
        if col not in canon.columns:
            canon[col] = default
    canon['gf_tier'] = canon['_tier']

    lab1 = pd.read_csv(LAB_FE1)
    lab2 = pd.read_csv(LAB_FE2)[['source', 'wavelength_air_A', 'elo_eV', 'loggf',
                                 'e_loggf_dex']]
    print(f"canonical_gf : {len(canon)} rows  "
          f"(Fe I {int((canon._species == 'Fe I').sum())}, "
          f"Fe II {int((canon._species == 'Fe II').sum())})")
    print(f"lab Fe I     : {len(lab1)}  " + ', '.join(
        f"{LAB_TAG[k]} {v}" for k, v in lab1.source.value_counts().items()))
    print(f"lab Fe II    : {len(lab2)}  DH19")

    if a.scan_tol:
        for species, lab in (('Fe I', lab1), ('Fe II', lab2)):
            sub = canon[canon._species == species]
            dw = sub.wavelength_air_A.to_numpy(float)
            de = sub.excitation_potential_eV.to_numpy(float)
            print(f"\n  {species} tolerance scan (lab line -> canonical row):")
            for wt in (0.005, 0.01, 0.02, 0.05):
                for et in (0.002, 0.005, 0.05):
                    idx = match_unique(lab.wavelength_air_A.to_numpy(float),
                                       lab.elo_eV.to_numpy(float), dw, de, wt, et)
                    nul = null_rate(lab.wavelength_air_A.to_numpy(float),
                                    lab.elo_eV.to_numpy(float), dw, de, wt, et,
                                    a.null_trials, 945)
                    print(f"    wtol={wt:<6} eptol={et:<6} unique={int((idx >= 0).sum()):4d} "
                          f" ambiguous={int((idx == -2).sum()):3d}  null={100 * nul:5.1f}%")
        return

    OUT_AUDIT.mkdir(parents=True, exist_ok=True)
    report: dict = {'ticket': 'RYA-945', 'tolerances': {
        'wavelength_A': a.wtol_A, 'excitation_eV': a.eptol_eV}, 'sources': {}}

    lab_m = pd.concat([match_lab(canon, lab1, 'Fe I', a.wtol_A, a.eptol_eV),
                       match_lab(canon, lab2, 'Fe II', a.wtol_A, a.eptol_eV)],
                      ignore_index=True)
    lab_win, hfs_refused = resolve_lab_precedence(lab_m)

    print("\nLAB CROSSMATCH (lab line -> canonical row, wavelength AND EP)")
    for species, lab in (('Fe I', lab1), ('Fe II', lab2)):
        sub = canon[canon._species == species]
        nul = null_rate(lab.wavelength_air_A.to_numpy(float), lab.elo_eV.to_numpy(float),
                        sub.wavelength_air_A.to_numpy(float),
                        sub.excitation_potential_eV.to_numpy(float),
                        a.wtol_A, a.eptol_eV, a.null_trials, 945)
        g = lab_m[lab_m.species == species]
        rate = float((g.status == 'matched').mean())
        print(f"  {species:6s} {int((g.status == 'matched').sum()):4d}/{len(g)} matched "
              f"({100 * rate:.1f}%)   ambiguous {int((g.status == 'ambiguous').sum())}   "
              f"absent {int((g.status == 'absent').sum())}   "
              f"RANDOMISED NULL {100 * nul:.1f}%  excess {100 * (rate - nul):+.1f} pts")
        if rate <= nul:
            raise SystemExit(f"STOP: {species} match rate is at or below chance — this "
                             f"matcher is finding coincidences, not identifying lines")
        report['sources'][species] = {'matched': int((g.status == 'matched').sum()),
                                      'ambiguous': int((g.status == 'ambiguous').sum()),
                                      'absent': int((g.status == 'absent').sum()),
                                      'lab_lines': int(len(g)),
                                      'randomised_null_rate': round(nul, 4)}
    if len(hfs_refused):
        print(f"  HFS refusal      : {len(hfs_refused)} lab line(s) fell on a "
              f"MULTI-COMPONENT canonical row — not written (a per-transition lab gf "
              f"cannot replace a summed one)")

    nist_all = []
    for species in ('Fe I', 'Fe II'):
        n = load_nist(species)
        report.setdefault('nist_files', {})[species] = n.attrs['files']
        nist_all.append(match_nist(canon, n, species, a.wtol_A, a.eptol_eV))
    nist_m = pd.concat(nist_all, ignore_index=True)
    print(f"\nNIST ASD        : {len(nist_m)} rows; "
          f"{int(nist_m.log_gf.notna().sum())} with a log gf; "
          f"{int(nist_m.nist_grade.isin(NIST_ADOPT_GRADES).sum())} graded C or better; "
          f"{int((nist_m.canon_index >= 0).sum())} identify a unique canonical row")
    shared = int(((nist_m.canon_index >= 0) & (nist_m.canon_shared_with > 1)).sum())
    print(f"  {shared} ASD row(s) share a canonical row with another ASD row — "
          f"neither identifies it, both refused")

    plan = plan_updates(canon, lab_win, nist_m)
    print(f"  RYA-834 band     : {plan_updates.n_nist_recorded_not_adopted} NIST row(s) in "
          f"9199.9-12935 A RECORDED, NOT ADOPTED — there NIST agrees with the line list to "
          f"~0.0003 dex, so it is the same number twice and cannot referee it")
    dropped = plan[~plan.kept]
    plan_ok = plan[plan.kept]
    print(f"\nPLAN            : {len(plan)} candidate rewrites, {len(plan_ok)} kept, "
          f"{len(dropped)} refused as a tier DOWNGRADE or lateral move")
    if len(dropped):
        print("  refused:", dict(dropped.groupby(['from_tier', 'to_tier']).size()))
    for tier in (TIER_LAB, TIER_NIST):
        p = plan_ok[plan_ok.to_tier == tier]
        if not len(p):
            continue
        d = p.delta.dropna()
        print(f"  -> {tier:8s} {len(p):5d} rows   delta(new-old) median {d.median():+.4f} "
              f" MAD {1.4826 * (d - d.median()).abs().median():.4f}  "
              f"|max| {d.abs().max():.3f} dex   |delta|>0.02: {int((d.abs() > 0.02).sum())}")
        print(f"     from: {dict(p.from_tier.value_counts())}")

    before_diag = int(((canon._tier.isin([TIER_LAB, TIER_NIST]))
                       & (pd.to_numeric(canon.excitation_potential_eV,
                                        errors='coerce') >= DIAG_EP_MIN_EV)
                       & canon._species.isin(['Fe I', 'Fe II'])).sum())
    out = apply_plan(canon, plan)
    out, demoted = demote_unverified_nist(out)
    if len(demoted):
        print(f"\nDEMOTED         : {len(demoted)} pre-existing row(s) claimed the NIST "
              f"tier without a sigma this run could verify -> demoted_nist_rows.csv")
        for k, v in demoted.demoted_because.value_counts().items():
            print(f"    {v:2d}  {k}")
        demoted.to_csv(OUT_AUDIT / 'demoted_nist_rows.csv', index=False)
    out = flag_diagnostic(out)
    after_diag = int(out.is_diagnostic.sum())
    print(f"\nDIAGNOSTIC FLAG : {before_diag} -> {after_diag} Fe I+II lines "
          f"(lab or NIST C+, EP >= {DIAG_EP_MIN_EV} eV, not registry-excluded)")
    print(f"  registry exclusions applied: "
          f"{int(out.diagnostic_excluded_by_registry.sum())}")

    cov = coverage(out)
    print("\nCOVERAGE (rows by gf tier)\n" + cov.to_string())

    # TICKET ITEM 3, ANSWERED FROM THE DATA: do the rows canonical_gf ALREADY attributes
    # to a lab paper actually carry that paper's number? A reference code is a claim; the
    # delta is the check. RYA-799's SCALE-MISMATCH is exactly this claim failing.
    m = lab_m[lab_m.canon_index >= 0].copy()
    m['canon_reference'] = m.canon_reference.astype(str)
    ver = (m.groupby('canon_reference')
             .agg(n=('delta_lab_minus_canon', 'size'),
                  median_delta=('delta_lab_minus_canon', 'median'),
                  max_abs_delta=('delta_lab_minus_canon', lambda d: d.abs().max()),
                  n_agree_002=('delta_lab_minus_canon', lambda d: int((d.abs() <= 0.02).sum())))
             .sort_values('n', ascending=False))
    print("\nEXISTING ATTRIBUTION vs THE LAB PAPER (ticket item 3)")
    print(ver.head(12).to_string(float_format=lambda x: f"{x:+.3f}"))
    report['attribution_check'] = json.loads(ver.reset_index().to_json(orient='records'))

    lab_m.to_csv(OUT_AUDIT / 'lab_crossmatch.csv', index=False)
    nist_m.drop(columns=['ei_ek_raw'], errors='ignore').to_csv(
        OUT_AUDIT / 'nist_crossmatch.csv', index=False)
    # 🔴 NOTHING MOVES BY A DEX WITHOUT BEING ON A LIST. The overwhelming result of this
    # ingest is a NULL on the VALUES -- per TP reference the median delta is 0.000 with a
    # MAD of 0.000-0.002, which is RYA-834's "the same number arriving twice". What is NOT
    # null is the handful of lines where the graded source and the incumbent genuinely
    # disagree. Those are written out by themselves so a reader can audit every one.
    big = plan_ok[plan_ok.delta.abs() > OUTLIER_DEX].copy()
    big['wavelength_air_A'] = canon.loc[big.canon_index, 'wavelength_air_A'].to_numpy()
    big['species'] = canon.loc[big.canon_index, '_species'].to_numpy()
    big['excitation_potential_eV'] = canon.loc[
        big.canon_index, 'excitation_potential_eV'].to_numpy()
    big['old_log_gf'] = canon.loc[big.canon_index, 'log_gf'].to_numpy()
    big['old_reference'] = canon.loc[big.canon_index, 'loggf_reference'].to_numpy()
    big.sort_values('delta', key=lambda d: d.abs(), ascending=False).to_csv(
        OUT_AUDIT / 'outlier_ledger.csv', index=False)
    print(f"\nOUTLIER LEDGER  : {len(big)} row(s) move by more than {OUTLIER_DEX} dex "
          f"(of {len(plan_ok)} rewrites) -> outlier_ledger.csv")
    if len(big):
        show = big.sort_values('delta', key=lambda d: d.abs(), ascending=False).head(8)
        for r in show.itertuples():
            print(f"    {r.species:6s} {r.wavelength_air_A:9.3f} A  EP {r.excitation_potential_eV:5.2f}  "
                  f"{r.old_log_gf:+7.3f} ({str(r.old_reference)[:12]:12s}) -> "
                  f"{r.log_gf:+7.3f} ({r.lab_source_tag})  {r.delta:+.3f} dex")
    report['outliers_over_%.2f_dex' % OUTLIER_DEX] = int(len(big))

    plan.to_csv(OUT_AUDIT / 'update_plan.csv', index=False)
    spread = resolve_lab_precedence.last_spread
    if len(spread):
        print(f"\nINTER-SOURCE SPREAD: {len(spread)} canonical row(s) measured by more "
              f"than one lab; |spread| median "
              f"{spread.inter_source_spread_dex.abs().median():.3f} dex, max "
              f"{spread.inter_source_spread_dex.abs().max():.3f} dex")
        print("  " + ', '.join(f"{k}:{v}" for k, v in
                               spread.sources.value_counts().items()))
        spread.to_csv(OUT_AUDIT / 'lab_inter_source_spread.csv', index=False)
    report.update({
        'multi_lab_rows': int(len(spread)),
        'inter_source_spread_dex': {
            'median': float(spread.inter_source_spread_dex.abs().median()) if len(spread) else None,
            'max': float(spread.inter_source_spread_dex.abs().max()) if len(spread) else None},
        'canonical_rows': int(len(out)),
        'rewrites_kept': int(len(plan_ok)),
        'rewrites_refused_as_downgrade': int(len(dropped)),
        'hfs_refused_lab_lines': int(len(hfs_refused)),
        'nist_rows_sharing_a_canonical_row': shared,
        'diagnostic_before': before_diag, 'diagnostic_after': after_diag,
        'registry_excluded': int(out.diagnostic_excluded_by_registry.sum()),
        'demoted_unverified_nist': int(len(demoted)),
        'demoted_note': (
            'Pre-existing rows that claimed a NIST grade this run could not stand behind: '
            '5 sit on 2-component HFS clusters (RYA-822) and 4 carry a stored v5.11 letter '
            'with no verified sigma (RYA-354) — two of which RYA-852 already showed are '
            'stale by 4-7x. Demoted to OTHER: value and note kept, claim withdrawn.'),
        'lab_source_dois': LAB_DOI, 'lab_source_access': LAB_VIZIER,
        'nist_adopt_grades': list(NIST_ADOPT_GRADES),
        'nist_max_sigma_dex': NIST_MAX_SIGMA_DEX,
        'nist_recorded_not_adopted_rya834_band': int(
            plan_updates.n_nist_recorded_not_adopted),
        'ru_code_correction': (
            "The ticket calls the 5,113 `RU` rows a partial Ruffoni-2014 wiring. They are "
            "all Fe II, 4200-9199 A; Ruffoni 2014 is an Fe I paper. VALD's own reference "
            "block reads `Fe 2: Raassen & Uylings ... gf:RU`, so `RU` is Raassen & "
            "Uylings, a theoretical calculation. Ruffoni's Fe I values reach canonical_gf "
            "under the GES codes `GESHRL14` / `2014MNRAS.` instead."),
        'blend_flag_note': (
            "`blend_flag` is a column of the MEASURED EW table, set per observation and "
            "applied by abundances_derive._prefilter. A line list cannot carry it. The "
            "line-level criterion used here is the problem-children registry's standing "
            "`exclude` rows; the measurement-time filter still runs downstream."),
        'coverage': json.loads(cov.reset_index().to_json(orient='records')),
        'generated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    })
    (OUT_AUDIT / 'rya945_summary.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(f"\n  wrote {OUT_AUDIT}/lab_crossmatch.csv, nist_crossmatch.csv, "
          f"update_plan.csv, rya945_summary.json")

    if not a.apply:
        print("\n  [dry-run] re-run with --apply to write canonical_gf.csv")
        return
    out.drop(columns=['_species', '_tier']).to_csv(CANON, index=False)
    print(f"  wrote {CANON}  ({len(out)} rows, "
          f"{len(out.columns) - 2} columns)")


if __name__ == '__main__':
    main()
