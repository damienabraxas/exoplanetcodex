#!/usr/bin/env python3
"""Inventory existing Fe inputs and stamped terms; never infer a need to rerun."""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def build():
    products = json.loads((ROOT / 'data/products/solar/Fe.json').read_text())['products']
    rows = []
    for p in products:
        src = ROOT / p['provenance']['copied_to']
        bud = src.with_name(src.name.replace('_products.csv', '_budgets.txt'))
        if not src.name.endswith('_products.csv') or not bud.is_file():
            bud = None
        text = bud.read_text() if bud else ''
        rows.append({
            'key': '|'.join(str(p.get(k, '')) for k in
                ('element', 'ion', 'band', 'instrument', 'holding', 'tier', 'selector', 'route', 'treatment')),
            'measurement_artifact': str(src.relative_to(ROOT)),
            'measurement_artifact_present': src.is_file(),
            'measurement_sha256': hashlib.sha256(src.read_bytes()).hexdigest() if src.is_file() else None,
            'budget_artifact': str(bud.relative_to(ROOT)) if bud else None,
            'sigma_stat': p['sigma_stat'], 'sigma_syst': p['sigma_syst'],
            'sigma_xi': p.get('sigma_xi'), 'xi_state': p.get('xi_state'),
            'xi_source': p.get('xi_source'), 'xi_note': p.get('xi_note'),
            'recorded_budget_terms': [line.strip() for line in text.splitlines()
                                     if '[SYSTEMATIC]' in line or '[random]' in line],
            'evidence_status': 'existing_source_artifact_verified; not a new measurement request'
                               if src.is_file() else 'source_artifact_not_available_in_this_checkout',
            'mean3d_experimental': 'mean3D' in p['treatment'],
        })
    return {
        'ticket': 'RYA-587',
        'scope': 'Existing Fe evidence reconciliation after recovered Reference xi ingestion',
        'products': rows,
        'measurement_artifacts_present': sum(r['measurement_artifact_present'] for r in rows),
        'direct_budget_artifacts_present': sum(r['budget_artifact'] is not None for r in rows),
        'product_count': len(rows),
        'xi_states': dict(collections.Counter(p.get('xi_state') for p in products)),
        'interpretation': (
            'Missing new RYA-587 schema fields do not imply missing abundance runs. '
            'Source budget texts may predate separately stamped xi campaigns; do not treat '
            'those old UNMEASURED stellar labels as absence of completed xi work. '
            'Non-Reference xi labels must be reconciled against exact pools and historical '
            'artifacts; no blanket rerun is authorized or scheduled.'),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = build()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'products'}, indent=2))
