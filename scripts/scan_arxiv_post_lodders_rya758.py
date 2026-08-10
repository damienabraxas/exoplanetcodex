#!/usr/bin/env python3
"""
scripts/scan_arxiv_post_lodders_rya758.py
========================================
RYA-758 — post-Lodders-2025 arXiv scan.

Lodders/Bergemann/Palme 2025 was accepted 2025-02-07, so anything submitted on or
after 2025-03-01 is outside its evidence base.  Rather than one query per element
(52 requests -> arXiv 429s), this issues a handful of BROAD NLTE queries, then
filters the returned titles/abstracts locally for the outside-28 element names.
Every query outcome is recorded, successes and failures alike -- a query that
fails is written into the artifact as an ERROR row, never dropped.
"""
import json, re, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / 'data' / 'audit' / 'nlte_grid_inventory_beyond28_arxiv_scan.json'

NS = {'a': 'http://www.w3.org/2005/Atom'}
CUTOFF = '2025-03-01'

QUERIES = [
    'abs:"non-LTE" AND cat:astro-ph.SR',
    'abs:"NLTE" AND cat:astro-ph.SR',
    'abs:"model atom"',
    'abs:"departure coefficients"',
    'abs:"solar abundance"',
]

# Outside-28 universe (ticket section 2) + the element word forms to match on.
ELEMENTS = {
    'Be': ['beryllium'], 'B': ['boron'], 'F': ['fluorine'], 'Cl': ['chlorine'],
    'Ne': ['neon'], 'Ar': ['argon'], 'Ga': ['gallium'], 'Ge': ['germanium'],
    'Se': ['selenium'], 'Br': ['bromine'], 'Kr': ['krypton'], 'Rb': ['rubidium'],
    'Nb': ['niobium'], 'Mo': ['molybdenum'], 'Ru': ['ruthenium'],
    'Rh': ['rhodium'], 'Pd': ['palladium'], 'Ag': ['silver'], 'Cd': ['cadmium'],
    'In': ['indium'], 'Sn': ['tin'], 'Sb': ['antimony'], 'Te': ['tellurium'],
    'I': ['iodine'], 'Xe': ['xenon'], 'La': ['lanthanum'], 'Ce': ['cerium'],
    'Pr': ['praseodymium'], 'Nd': ['neodymium'], 'Sm': ['samarium'],
    'Gd': ['gadolinium'], 'Tb': ['terbium'], 'Dy': ['dysprosium'],
    'Ho': ['holmium'], 'Er': ['erbium'], 'Tm': ['thulium'], 'Yb': ['ytterbium'],
    'Lu': ['lutetium'], 'Hf': ['hafnium'], 'Ta': ['tantalum'],
    'W': ['tungsten'], 'Re': ['rhenium'], 'Os': ['osmium'], 'Ir': ['iridium'],
    'Pt': ['platinum'], 'Au': ['gold'], 'Hg': ['mercury'], 'Tl': ['thallium'],
    'Pb': ['lead'], 'Bi': ['bismuth'], 'Th': ['thorium'], 'U': ['uranium'],
}


def fetch(query, start, page=100, tries=6):
    url = 'http://export.arxiv.org/api/query?' + urllib.parse.urlencode({
        'search_query': query, 'start': start, 'max_results': page,
        'sortBy': 'submittedDate', 'sortOrder': 'descending'})
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return ET.fromstring(fh.read())
        except Exception as exc:
            last = exc
            wait = 30 * (attempt + 1)
            print(f'  retry {attempt + 1}/{tries} in {wait}s after {exc}', flush=True)
            time.sleep(wait)
    raise last


papers, queries_log = {}, []
for q in QUERIES:
    got, err, start = 0, None, 0
    try:
        while start < 300:                       # 3 pages per query is plenty
            root = fetch(q, start)
            entries = root.findall('a:entry', NS)
            if not entries:
                break
            stop = False
            for e in entries:
                pub = e.find('a:published', NS).text[:10]
                if pub < CUTOFF:
                    stop = True                  # descending order -> done
                    continue
                aid = e.find('a:id', NS).text.rsplit('/', 1)[-1]
                papers[aid] = {
                    'arxiv_id': aid, 'published': pub,
                    'title': ' '.join(e.find('a:title', NS).text.split()),
                    'summary': ' '.join(e.find('a:summary', NS).text.split()),
                }
                got += 1
            if stop:
                break
            start += 100
            time.sleep(8)
    except Exception as exc:
        err = f'{type(exc).__name__}: {exc}'
        print(f'[QUERY FAILED] {q}: {err}', flush=True)
    queries_log.append({'query': q, 'n_kept': got, 'error': err})
    print(f'{q!r}: kept {got} since {CUTOFF}  err={err}', flush=True)
    time.sleep(8)

# Local element filter over the union of returned papers.
by_element = {}
for sym, words in ELEMENTS.items():
    hits = []
    for p in papers.values():
        blob = (p['title'] + ' ' + p['summary']).lower()
        # whole-word match only: 'tin' must not fire on 'resulting', 'lead' not on 'leading'
        if any(re.search(rf'\b{w}\b', blob) for w in words):
            hits.append({k: p[k] for k in ('arxiv_id', 'published', 'title')})
    by_element[sym] = sorted(hits, key=lambda h: h['published'], reverse=True)

OUT.parent.mkdir(parents=True, exist_ok=True)
json.dump({'cutoff': CUTOFF, 'queries': queries_log,
           'n_papers_scanned': len(papers), 'by_element': by_element,
           'papers': sorted(papers.values(), key=lambda p: p['published'], reverse=True)},
          OUT.open('w'), indent=1)
print(f'scanned {len(papers)} papers; wrote {OUT}')
for sym, hits in by_element.items():
    if hits:
        print(f'  {sym}: {len(hits)}')
