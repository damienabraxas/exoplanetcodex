#!/usr/bin/env python3
"""
RYA-822 — pull NIST ASD Fe I 3000-3780 A: the primary gf source for the near-UV band.

⚠️ THE RAW-CURL RECIPE IN THE PROJECT'S NOTES IS DEAD, AND THE REASON IS NARROWER THAN
"the endpoint is down". Bisected 2026-08-15 against `lines1.pl`:

    spectra + format=3 + low_w/upp_w      -> ok
    ... + unit=0   (Angstrom, our recipe) -> HTTP 200 "Software error ... line 701"
    ... + unit=1   (nm)                   -> same crash
    ... + unit=2   (um)                   -> same crash
    ... + unit=3   (invalid)              -> input-error page
    ... omitting `unit` entirely          -> input-error page with an EMPTY message

So the fault is server-side in the `unit` parameter path, not in our query, and there is no
value of `unit` that works. Re-permuting the other parameters is wasted effort.

⚠️ BUT THE ASD LEG IS *NOT* BLOCKED, WHICH IS THE OPPOSITE OF WHAT THE NOTES CONCLUDED.
`astroquery.nist` reaches the same database and returns MORE than the curl recipe ever did:
`fik`, `Acc.` (the accuracy grade), `gi gk` (so log gf = log10(gi * fik)), and — the part
the notes recorded as impossible via this endpoint — the bibliographic reference codes
`TP` (transition probability) and `Line`.

🔴 `wavelength_type` MUST BE 'vac+air'. astroquery DEFAULTS TO 'vacuum', and our line
lists are AIR above 2000 A. The difference is ~0.95 A at 3400 A -- about five Fe I line
spacings in this band -- so the default silently compares two different wavelength scales.
It does not look like an error; it looks like "NIST and VALD disagree about which lines
exist". Caught by an offset scan whose match rate peaked at -1.00 A instead of 0.00, and
NOT adopted in the meantime only because the matcher also requires the excitation
potential to agree.

⚠️ AND A NIST CHECK IS NOT INDEPENDENT OF EVERY SOURCE. NIST's Fe I values in this band are
largely compilations; agreement proves only "no transcription error", and DISAGREEMENT is
the informative outcome. The per-line reference code is recorded so that a later reader can
tell which underlying measurement is actually being cited, instead of treating "NIST" as a
single authority.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'linelists' / 'primary_gf'

#: NIST accuracy ladder, worst-to-best percentage. Enumerated IN FULL because omitting the
#: '+' tiers once put a B+ line (<=7%) BELOW a B line (<=10%) — an inverted ladder that
#: silently demoted the better measurement (found on RYA-592).
NIST_ACC_PCT = {
    'AAA': 0.3, 'AA': 1.0, 'A+': 2.0, 'A': 3.0, 'B+': 7.0, 'B': 10.0,
    'C+': 18.0, 'C': 25.0, 'D+': 40.0, 'D': 50.0, 'E': 100.0,
}


def pull(lo_A: float, hi_A: float, step_A: float, species: str,
         pause_s: float) -> pd.DataFrame:
    from astroquery.nist import Nist
    import astropy.units as u

    frames, edges = [], np.arange(lo_A, hi_A, step_A)
    for i, a in enumerate(edges):
        b = min(a + step_A, hi_A)
        for attempt in range(3):
            try:
                t = Nist.query(a * u.AA, b * u.AA, linename=species,
                               wavelength_type='vac+air')
                break
            except Exception as e:                     # transient network, not data
                if attempt == 2:
                    print(f"  {a:.0f}-{b:.0f} FAILED after 3 tries: "
                          f"{type(e).__name__}: {str(e)[:120]}")
                    t = None
                    break
                time.sleep(2 + 3 * attempt)
        if t is None or len(t) == 0:
            print(f"  {a:7.1f}-{b:7.1f} A   0 rows")
            continue
        df = t.to_pandas()
        df['chunk_lo_A'], df['chunk_hi_A'] = a, b
        frames.append(df)
        print(f"  {a:7.1f}-{b:7.1f} A   {len(df):4d} rows "
              f"({i + 1}/{len(edges)})", flush=True)
        time.sleep(pause_s)
    if not frames:
        raise SystemExit('NIST returned nothing across the whole band — refusing to '
                         'write an empty pull that would read as "no data exists"')
    return pd.concat(frames, ignore_index=True)


def tidy(df: pd.DataFrame) -> pd.DataFrame:
    """Reduce to the columns an adjudication needs, and compute log gf HONESTLY.

    log gf = log10(g_i * f_ik). A row without BOTH f_ik and g_i has NO log gf — it is
    left as NaN and kept, because "NIST lists this line but publishes no oscillator
    strength for it" is a real and useful state, distinct from "NIST does not have it".
    """
    def _gi(s):
        try:
            return float(str(s).split('-')[0].strip())
        except (ValueError, AttributeError, IndexError):
            return np.nan

    out = pd.DataFrame({
        'wavelength_obs_A': pd.to_numeric(df.get('Observed'), errors='coerce'),
        'wavelength_ritz_A': pd.to_numeric(df.get('Ritz'), errors='coerce'),
        'fik': pd.to_numeric(df.get('fik'), errors='coerce'),
        'gi': df.get('gi   gk', pd.Series(dtype=object)).map(_gi),
        'aki_s-1': pd.to_numeric(df.get('Aki'), errors='coerce'),
        'nist_grade': df.get('Acc.'),
        'ref_transition_probability': df.get('TP'),
        'ref_line': df.get('Line'),
        # RYA-780: match on wavelength AND excitation potential. At this band's Fe I
        # density a 0.02 A wavelength coincidence happens ~20% of the time by chance,
        # so wavelength alone cannot identify a line well enough to ADOPT its value.
        'ei_ek_raw': df.get('Ei           Ek'),
    })
    def _ei(s):
        try:
            return float(str(s).split('-')[0].strip().replace('[','').replace(']','')
                         .replace('(','').replace(')','').replace('+x','').replace('?',''))
        except (ValueError, AttributeError, IndexError):
            return float('nan')
    # astroquery returns Ei/Ek ALREADY IN eV ('0.85899575  -   4.99127081'), so there
    # is no cm-1 conversion to do. Converting anyway produced 1e-4 eV excitation
    # potentials for Fe I -- physically impossible, which is how it was caught.
    out['ei_eV'] = out['ei_ek_raw'].map(_ei)
    out['log_gf'] = np.where(out.fik.notna() & out.gi.notna(),
                             np.log10(out.gi.astype(float) * out.fik.astype(float)),
                             np.nan)
    out['nist_grade'] = out['nist_grade'].astype(str).str.strip().replace(
        {'nan': None, '--': None, '': None})
    out['nist_acc_pct'] = out['nist_grade'].map(NIST_ACC_PCT)
    # the wavelength to match on: observed where it exists, else Ritz
    out['wavelength_A'] = out['wavelength_obs_A'].fillna(out['wavelength_ritz_A'])
    return out.dropna(subset=['wavelength_A']).sort_values('wavelength_A')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo-A', type=float, default=3000.0)
    ap.add_argument('--hi-A', type=float, default=3780.0)
    ap.add_argument('--step-A', type=float, default=20.0)
    ap.add_argument('--species', default='Fe I')
    ap.add_argument('--pause-s', type=float, default=0.5)
    a = ap.parse_args()

    import astroquery
    print(f"astroquery {astroquery.__version__}  ->  NIST ASD, {a.species} "
          f"{a.lo_A:.0f}-{a.hi_A:.0f} A in {a.step_A:.0f} A chunks")
    raw = pull(a.lo_A, a.hi_A, a.step_A, a.species, a.pause_s)
    tid = tidy(raw)

    with_gf = tid.log_gf.notna()
    graded = tid.nist_grade.notna()
    OUT.mkdir(parents=True, exist_ok=True)
    slug = a.species.replace(' ', '')
    dest = OUT / f'nist_asd_{slug}_{int(a.lo_A)}_{int(a.hi_A)}.tsv'
    tid.to_csv(dest, sep='\t', index=False)

    prov = {
        'ticket': 'RYA-822',
        'source': 'NIST Atomic Spectra Database (ASD), lines',
        'access': f'astroquery.nist {astroquery.__version__} (Nist.query)',
        'why_not_curl': (
            'the project recipe uses lines1.pl with unit=0; bisected 2026-08-15, the '
            '`unit` parameter crashes the CGI server-side for EVERY valid value '
            '(0=A, 1=nm, 2=um -> "Software error ... line 701") and omitting it returns '
            'an input-error page with an empty message. The fault is in the unit path, '
            'not the query. astroquery reaches the same database and additionally '
            'returns the TP/Line bibliographic codes the curl route never exposed.'),
        'species': a.species, 'band_A': [a.lo_A, a.hi_A], 'chunk_A': a.step_A,
        'n_rows_raw': int(len(raw)), 'n_rows_tidy': int(len(tid)),
        'n_with_log_gf': int(with_gf.sum()), 'n_graded': int(graded.sum()),
        'grade_counts': {k: int(v) for k, v in
                         tid.nist_grade.value_counts().items()},
        'caveat': ('NIST Fe I in this band is largely a COMPILATION. Agreement with a '
                   'compilation proves only "no transcription error"; DISagreement is '
                   'the informative outcome. Per-line TP/Line reference codes are kept '
                   'so a reader can see which measurement is actually cited.'),
        'pulled_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    (OUT / f'nist_asd_{slug}_{int(a.lo_A)}_{int(a.hi_A)}.prov.json').write_text(
        json.dumps(prov, indent=2))

    print(f"\n  rows (raw / tidy)      {len(raw)} / {len(tid)}")
    print(f"  with a log gf          {int(with_gf.sum())}  "
          f"({100 * with_gf.mean():.1f}%)")
    print(f"  with a NIST grade      {int(graded.sum())}")
    print(f"  grades: {dict(tid.nist_grade.value_counts())}")
    print(f"\n[out] {dest}")


if __name__ == '__main__':
    main()
