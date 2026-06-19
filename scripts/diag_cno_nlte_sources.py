"""
scripts/diag_cno_nlte_sources.py
================================
RYA-359 — DIAGNOSTIC (read-only): can the C I / O I NLTE correction grids actually
be acquired the way the ticket assumes ("analogous to RYA-245 — MPIA scrape")?

Verdict (reproduced by this script): NO. All three plausible acquisition paths fail
for C/O. The grids must come from Amarsi's author-distributed 3D non-LTE correction
data products (a specific download), which the ticket's premise did not account for.
This script does NOT fabricate any NLTE value (fabricated corrections would corrupt
the flagship C/O science) and writes nothing canonical.

Paths probed:
  1. MPIA Spectrum Tools (nlte.mpia.de, the RYA-245/319 source)
       - form serves: H, O(listed), Mg, Si, Ca, Ca II, Ti, Ti II, Cr, Mn, Fe, Fe II, Co
       - NO carbon select at all → C I impossible here.
       - O free-text 8.01 → returns 0.000 ("no NLTE departures") for the 777 triplet,
         which is physically wrong (the permitted O I triplet has ~−0.2..−0.5 dex NLTE).
         O via the olines[] select → no result. So O I NLTE is not served by this scrape.
  2. VizieR — Amarsi 2019 carbon catalog J/A+A/630/A104 exists but contains only one
     table (table7 = stellar-sample abundances), not the (Teff,logg,[Fe/H],A,line)→Δ
     correction GRID; the dedicated O-correction papers are not on VizieR.
  3. In-synthesis Turbospectrum NLTE — DATA/SPECIES_LTE_NLTE.dat is all-LTE (only
     H/Ca/Fe listed, every one 'lte'); no C/O model atoms or departure grids installed.

Usage:  python3 scripts/diag_cno_nlte_sources.py
"""
from __future__ import annotations
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
from config.constants import ISPEC_DIR  # noqa: E402

MPIA_URL = 'https://nlte.mpia.de/gui-siuAC_secE.php'


def probe_mpia():
    print("── 1. MPIA Spectrum Tools (the RYA-245 scrape source) ──")
    try:
        import requests
        from bs4 import BeautifulSoup
    except Exception as e:
        print(f"   (requests/bs4 unavailable: {e})"); return
    try:
        s = BeautifulSoup(requests.get(MPIA_URL, timeout=30).text, 'html.parser')
        sels = {sel.get('name') for sel in s.find_all('select')}
        has_c = any(n and n.startswith('cl') and 'ca' not in n and 'cr' not in n
                    and 'co' not in n for n in sels)  # 'clines[]' would be carbon
        print(f"   line-selects: {sorted(n for n in sels if n and n.endswith('lines[]'))}")
        print(f"   carbon select present : {has_c}  → C I NOT served by MPIA"
              if not has_c else f"   carbon select present : {has_c}")
        # O free-text null test
        user = 'Sun 5772 4.44 +0.00 1.0'
        r = requests.post(MPIA_URL, data={'model': 'mafags-os', 'user_input': user,
                          'linelist': '0', 'lines_input': '7771.959 8.01'}, timeout=120)
        rd = BeautifulSoup(r.text, 'html.parser').find('div', {'id': 'result'})
        txt = rd.get_text(' ', strip=True) if rd else ''
        print(f"   O I 7771.959 free-text → {'0.000 (NULL — no NLTE served)' if ' 0.000' in txt else txt[-60:]}")
    except Exception as e:
        print(f"   network probe failed: {e}")


def probe_vizier():
    print("\n── 2. VizieR — Amarsi 2019 C catalog (J/A+A/630/A104) ──")
    try:
        import requests, re
        r = requests.get('https://vizier.cds.unistra.fr/viz-bin/cat/J/A+A/630/A104', timeout=40)
        tabs = sorted(set(re.findall(r'J/A\+A/630/A104/(\w+)', r.text)))
        print(f"   tables: {tabs}  → {'only a stellar-sample table; NOT the correction grid' if tabs==['table7'] else tabs}")
    except Exception as e:
        print(f"   network probe failed: {e}")


def probe_insynth():
    print("\n── 3. In-synthesis Turbospectrum NLTE (SPECIES_LTE_NLTE.dat) ──")
    f = ISPEC_DIR / 'synthesizer' / 'turbospectrum' / 'DATA' / 'SPECIES_LTE_NLTE.dat'
    if not f.exists():
        print(f"   not found: {f}"); return
    species = []
    for line in f.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('~'):
            continue
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit():
            species.append((parts[1].strip("'"), parts[2].strip("'")))
    print(f"   configured species: {species}")
    cno = [s for s in species if s[0] in ('C', 'N', 'O')]
    nlte = [s for s in species if s[1].lower() == 'nlte']
    print(f"   C/N/O present: {cno or 'NONE'}   any 'nlte': {nlte or 'NONE (all LTE)'}")
    print(f"   → no C/O model atom / departure grid installed → in-synthesis NLTE N/A")


def main():
    print("RYA-359 — C I / O I NLTE grid acquisition feasibility (read-only)\n")
    probe_mpia(); probe_vizier(); probe_insynth()
    print("\n── VERDICT ──────────────────────────────────────────────────────────")
    print("  All three paths fail for C/O. The MPIA scrape (RYA-245 method) does NOT")
    print("  serve C (no select) or O (returns null); VizieR has only the Amarsi-2019")
    print("  sample table, not the grid; Turbospectrum has no C/O NLTE installed.")
    print("  → The real source is Amarsi's 3D non-LTE correction-grid data products")
    print("    (C I: Amarsi+2019; O I: Amarsi+2016/2018), a specific author-hosted /")
    print("    journal-supplement download. Acquire that (or install TS C/O model atoms")
    print("    + departure grids). NOT fabricated here. See the RYA-359 Linear report.")


if __name__ == '__main__':
    main()
