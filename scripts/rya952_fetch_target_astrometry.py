#!/usr/bin/env python3
"""
RYA-952 — vendor SIMBAD astrometry for the CRIRES target set.

WHY A CACHED FILE RATHER THAN A LIVE QUERY. Target identity is the thing this ticket must
not get wrong, so the reference it is judged against has to be (a) independent of the FITS
headers being judged, and (b) reproducible offline, including in CI on a runner with no
outbound network. A live SIMBAD call satisfies (a) and fails (b); this script pulls once and
commits the answer, and `audit_crires` reads the committed file.

🔴 THE HEADER IS NOT AN INDEPENDENT REFEREE FOR THE HEADER. `ESO TEL TARG ALPHA/DELTA/PMA/PMD`
carry the catalogue position the telescope was *pointed at*, which is excellent corroboration
but is written by the same observing block that wrote `OBJECT`. When `OBJECT` is wrong, those
fields are wrong together with it or right together with it — either way they cannot referee
it. SIMBAD can.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'reference' / 'crires_target_astrometry.csv'

#: The stars any CRIRES frame in this project could plausibly be, plus the ones a mislabel
#: would plausibly land on. `alf Cen B` and the RYA-423 quarantine cases are here ON PURPOSE:
#: a catalogue that only contains the answers you expect cannot produce a NEGATIVE result.
TARGETS = {
    'tau_cet': 'tau Cet', 'eps_eri': 'eps Eri', 'tau_boo': 'tau Boo',
    '55cnc_a': '55 Cnc', 'alpha_cen_a': 'alf Cen A', 'alpha_cen_b': 'alf Cen B',
    'procyon': 'Procyon',
}
#: Solar-system bodies have no fixed position and are identified by a different route
#: entirely (RYA-372's ephemeris). Recorded so the absence is deliberate, not an oversight.
MOVING = {
    'vesta': 'Vesta (minor planet 4) — ephemeris target, no catalogue position',
    'sun': ('the Sun — SIMBAD resolves the name but publishes NO coordinates for it, which '
            'is correct and is why it cannot live in a fixed-position catalogue. Solar '
            'CRIRES frames are Vesta reflected-solar and identify through Vesta.'),
}


def _aliases(simbad, name: str) -> set[str]:
    """SIMBAD's identifier list for one star, normalised for comparison."""
    try:
        t = simbad.query_objectids(name)
    except Exception:
        return set()
    if t is None:
        return set()
    col = t.colnames[0]
    return {str(v).strip().lower().replace(' ', '') for v in t[col]}


def main() -> None:
    import warnings
    warnings.filterwarnings('ignore')
    from astroquery.simbad import Simbad

    s = Simbad()
    s.add_votable_fields('pmra', 'pmdec', 'V')
    names = list(TARGETS.values())
    t = s.query_objects(names)
    df = t.to_pandas()
    # ⚠️ astroquery returns `user_specified_id` with LITERAL SINGLE QUOTES around the value
    # ("'tau Cet'", not "tau Cet"), so matching on it finds nothing and every target looks
    # missing. Caught only because this script refuses to write a catalogue with a hole in
    # it rather than dropping the unmatched rows.
    df['user_specified_id'] = df.user_specified_id.astype(str).str.strip().str.strip("'\"")

    rows = []
    for key, name in TARGETS.items():
        m = df[df.user_specified_id == name]
        if not len(m):
            raise SystemExit(f"SIMBAD returned no row for {name!r} — refusing to write a "
                             f"catalogue with a hole in it")
        r = m.iloc[0]
        rows.append({
            'star_id': key, 'query_name': name, 'simbad_main_id': str(r.main_id),
            'ra_deg_j2000': float(r.ra), 'dec_deg_j2000': float(r.dec),
            'pm_ra_cosdec_mas_yr': float(r.pmra), 'pm_dec_mas_yr': float(r.pmdec),
            'v_mag': float(r.V) if pd.notna(r.V) else None,
            # Every catalogue designation SIMBAD knows for this star. Without it, a frame
            # honestly labelled `HD 22049` reads as "OBJECT does not name this star", and
            # the genuinely interesting mislabels -- a ROLE like `STD`, an observing-run
            # placeholder like `Star S5` -- are buried in a list of false positives.
            'aliases': '|'.join(sorted(_aliases(s, name))),
        })
    out = pd.DataFrame(rows).sort_values('star_id')
    if out[['ra_deg_j2000', 'dec_deg_j2000']].isna().any().any():
        raise SystemExit('a target came back without coordinates — refusing to write it')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    prov = {
        'ticket': 'RYA-952', 'source': 'SIMBAD (CDS) via astroquery.simbad',
        'epoch': 'ICRS J2000; proper motions in mas/yr, pm_ra is mu_alpha* (cos-dec applied)',
        'why_cached': ('target identity is what this ticket must not get wrong, so its '
                       'referee must be independent of the headers being judged AND '
                       'reproducible offline, including in CI with no network'),
        'moving_targets_excluded': MOVING,
        'retrieved': time.strftime('%Y-%m-%d'),
        'regenerate': 'python3 scripts/rya952_fetch_target_astrometry.py',
    }
    (OUT.with_suffix('.prov.json')).write_text(json.dumps(prov, indent=2) + '\n')
    print(out.to_string(index=False))
    print(f"\n[out] {OUT}")


if __name__ == '__main__':
    main()
