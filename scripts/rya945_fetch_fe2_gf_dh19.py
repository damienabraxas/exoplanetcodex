#!/usr/bin/env python3
"""
RYA-945 — vendor Den Hartog et al. 2019 (ApJS 243, 33) Table 6: the PRIMARY laboratory
Fe II log gf set, 131 UV/blue lines.

WHY THIS IS THE Fe II REFEREE, AND WHY MELENDEZ & BARBUY 2009 IS NOT.
RYA-852 read MB09's own definition: lines flagged `S` there take their absolute gf from an
"inverse analysis based on the NSO solar flux spectrum", and even its `L` lines admit
"slight corrections (usually no larger than 0.1 dex) to reproduce the solar spectrum
better". That is lab-NORMALISED, not lab-MEASURED, and the `S` half is exactly the
reverse-solar-analysis the RYA-161 firewall exists to catch. DH19 is branching fractions
from FTS spectra times radiative lifetimes — no solar input anywhere in the chain.

🔴 THIS TABLE IS NOT ON VizieR AND THAT IS VERIFIED, NOT ASSUMED.
`Vizier.find_catalogs('J/ApJS/243/33')` returns a bucket whose `.description` is None and
`Vizier.get_catalogs('J/ApJS/243/33')` returns ZERO tables — CDS never ingested this paper.
`https://cdsarc.cds.unistra.fr/ftp/J/ApJS/243/33/` is a hard 404. The arXiv e-print
(1907.11760) is a PDF-only submission, so the Belmonte-2017 trick of parsing the LaTeX
source does not apply either. The published PDF IS the machine-readable source of last
resort, so it is parsed here with the paper's own stated line count as the control.

⚠️ pypdf renders this paper's '±' as the PRIVATE-USE-AREA pair U+F0A0 '±' U+F0A0, and its
minus sign as U+2212. A regex written against the characters you SEE in the extracted text
matches zero rows. `_normalise` maps the private-use block to a space before parsing, and
the row count control is what turned a silent zero-row parse into a loud one.

⚠️ THE PDF COLUMN ORDER IS λair | EU | JU | EL | JL | A ± σA | log gf ± σ. `EL` is the
LOWER level and it is the fourth column, not the second — reading EU as the excitation
potential would put every line 5-6 eV too high and quietly destroy the crossmatch. The
parser therefore checks EU > EL on every row and refuses the file if that ever inverts.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / 'data' / 'reference' / 'fe_gf_lab'
OUT_CSV = OUT_DIR / 'fe2_lab_loggf_dh19.csv'

#: The paper's own abstract/table caption. This is the CONTROL on the extraction: a PDF
#: parse that returns a different number of rows has silently dropped or duplicated lines.
STATED_N_LINES = 131

#: CODATA 2018. Same constant the RYA-799 Fe I pull used, so the two lab tables share a
#: single energy scale.
CM1_PER_EV = 8065.543937

PDF_DEFAULT = (Path.home() / 'Documents' / 'Exoplanet Codex' / 'Reference documents'
               / 'Den_Hartog_2019_ApJS_243_33.pdf')

#: λair  EU  JU  EL  JL  A ±σA  loggf ±σ  — commas inside the cm-1 values, unicode minus
#: signs, and a '±' that may or may not be preceded by a space.
_ROW = re.compile(
    r'^\s*(\d{4}\.\d{3})\s+'                       # lambda_air, always 3 decimals
    r'([\d,]+\.\d+)\s+(\d+\.?\d*)\s+'              # E_upper cm-1, J_upper
    r'([\d,]+\.\d+)\s+(\d+\.?\d*)\s+'              # E_lower cm-1, J_lower
    r'([\d.]+)\s*±\s*([\d.]+)\s+'                  # A-value (1e6 s-1) +- sigma
    r'([−\-]?[\d.]+)\s*±\s*([\d.]+)\s*$'           # log gf +- sigma (dex)
)


def parse_pdf(pdf: Path) -> pd.DataFrame:
    import pypdf
    reader = pypdf.PdfReader(str(pdf))
    rows, seen_caption = [], False
    for page in reader.pages:
        text = page.extract_text() or ''
        if 'Experimental Atomic Transition Probabilities' in text:
            seen_caption = True
        if not seen_caption:
            continue
        for line in text.split('\n'):
            m = _ROW.match(_normalise(line))
            if not m:
                continue
            lam, eu, ju, el, jl, a, sa, gf, sgf = m.groups()
            rows.append({
                'wavelength_air_A': float(lam),
                'eup_cm1': float(eu.replace(',', '')),
                'j_up': float(ju),
                'elo_cm1': float(el.replace(',', '')),
                'j_lo': float(jl),
                'aki_1e6_s-1': float(a),
                'e_aki_1e6_s-1': float(sa),
                'loggf': float(gf.replace('−', '-')),
                'e_loggf_dex': float(sgf),
            })
    if not seen_caption:
        raise SystemExit(f"{pdf.name}: Table 6's caption never appeared — this is not the "
                         f"Den Hartog 2019 paper, refusing to parse")
    d = pd.DataFrame(rows)

    # ---- controls -------------------------------------------------------------------
    if len(d) != STATED_N_LINES:
        raise SystemExit(f"parsed {len(d)} rows but the paper states {STATED_N_LINES} — "
                         f"the extraction is lossy; refusing to vendor it")
    if d.wavelength_air_A.duplicated().any():
        dup = d.wavelength_air_A[d.wavelength_air_A.duplicated()].tolist()
        raise SystemExit(f"duplicate wavelengths from the PDF parse: {dup}")
    bad = d[d.eup_cm1 <= d.elo_cm1]
    if len(bad):
        raise SystemExit(f"{len(bad)} row(s) with E_upper <= E_lower — the two energy "
                         f"columns have been read in the wrong order")

    d['elo_eV'] = d.elo_cm1 / CM1_PER_EV
    d['eup_eV'] = d.eup_cm1 / CM1_PER_EV
    d['source'] = 'DenHartog2019'
    d['species'] = 'Fe II'

    # The level splitting must reproduce the quoted wavelength. Fe II above 2000 A is
    # AIR in this table, so the comparison is made in VACUUM via Edlen (1966) -- the same
    # check the RYA-799 Fe I pull ran, and the reason its worst residual is 1e-6 eV.
    lam_vac = d.wavelength_air_A * _n_air(d.wavelength_air_A)
    hc_eV_A = 12398.419843320026
    resid = ((d.eup_eV - d.elo_eV) - hc_eV_A / lam_vac).abs()
    if resid.max() > 5e-4:
        raise SystemExit(f"level splitting disagrees with the quoted wavelength by up to "
                         f"{resid.max():.2e} eV — energies or wavelengths misparsed")
    d.attrs['level_resid_eV'] = float(resid.max())
    return d.sort_values('wavelength_air_A').reset_index(drop=True)


def _normalise(line: str) -> str:
    """Strip the PDF font's private-use glyphs so the row regex can see the columns."""
    return ''.join(' ' if '\uE000' <= ch <= '\uF8FF' else ch for ch in line)


def _n_air(lam_air_A):
    """Edlen (1966) refractive index of standard air, as a function of AIR wavelength."""
    s2 = (1e4 / lam_air_A) ** 2
    return 1 + 1e-8 * (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2))


def main() -> None:
    pdf = Path(sys.argv[1]) if len(sys.argv) > 1 else PDF_DEFAULT
    if not pdf.exists():
        raise SystemExit(f"source PDF not found: {pdf}")
    d = parse_pdf(pdf)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cols = ['source', 'species', 'wavelength_air_A', 'elo_cm1', 'eup_cm1', 'elo_eV',
            'eup_eV', 'j_lo', 'j_up', 'aki_1e6_s-1', 'e_aki_1e6_s-1', 'loggf',
            'e_loggf_dex']
    d[cols].to_csv(OUT_CSV, index=False)

    prov = {
        'ticket': 'RYA-945',
        'artifact': OUT_CSV.name,
        'source': 'Den Hartog, Lawler, Sneden, Cowan & Brukhovesky 2019, ApJS 243, 33, Table 6',
        'doi': '10.3847/1538-4365/ab322e',
        'bibcode': '2019ApJS..243...33D',
        'bibliography_key': 'denhartog2019',
        'method': 'FTS branching fractions x radiative lifetimes — NO solar input',
        'why_primary': (
            'The Fe II referee that clears RYA-161. Melendez & Barbuy 2009 does not: its '
            'own paper labels lines S when the absolute gf came from an inverse analysis '
            'of the NSO solar flux spectrum, and allows ~0.1 dex solar tweaks even on its '
            'L lines (RYA-852).'),
        'not_on_vizier': (
            'VERIFIED, not assumed: Vizier.get_catalogs("J/ApJS/243/33") returns 0 tables '
            'and cdsarc.cds.unistra.fr/ftp/J/ApJS/243/33/ is a 404. arXiv:1907.11760 is a '
            'PDF-only submission, so there is no LaTeX source to parse either. Extracted '
            'from the published PDF with pypdf.'),
        'local_pdf': str(pdf),
        'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
        'controls': [
            f'COUNT: {len(d)} rows parsed vs {STATED_N_LINES} stated by the table caption PASS',
            f'UNIQUE: 0 duplicate wavelengths PASS',
            f'ORDERING: 0 rows with E_upper <= E_lower PASS',
            f'LEVELS: max |(Eup-Elo) - hc/lambda_vac| = {d.attrs["level_resid_eV"]:.2e} eV PASS',
            f'SIGMA: {int(d.e_loggf_dex.notna().sum())}/{len(d)} rows carry a finite dex '
            f'uncertainty, range {d.e_loggf_dex.min():+.3f}..{d.e_loggf_dex.max():+.3f} PASS',
        ],
        'span_A': [float(d.wavelength_air_A.min()), float(d.wavelength_air_A.max())],
        'n_below_3300A': int((d.wavelength_air_A < 3300).sum()),
        'n_blue': int((d.wavelength_air_A >= 4000).sum()),
        'retrieved': time.strftime('%Y-%m-%d'),
        'regenerate': 'python3 scripts/rya945_fetch_fe2_gf_dh19.py',
    }
    (OUT_DIR / 'fe2_lab_loggf_dh19.prov.json').write_text(
        json.dumps(prov, indent=2) + '\n', encoding='utf-8')

    print(f"Den Hartog 2019 Table 6 — {len(d)} Fe II lines "
          f"({d.wavelength_air_A.min():.3f}-{d.wavelength_air_A.max():.3f} A)")
    for c in prov['controls']:
        print('  ' + c)
    print(f"  sigma median {d.e_loggf_dex.median():.3f} dex")
    print(f"  {prov['n_below_3300A']} below 3300 A, {prov['n_blue']} blueward-optical "
          f"(>= 4000 A)")
    print(f"[out] {OUT_CSV}")


if __name__ == '__main__':
    main()
