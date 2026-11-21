"""Canonical target/census screen; absence of active feed is never proof of blindness."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def screen(destination=None):
    targets_path=ROOT/'data/config/elements_master.json'
    census_path=ROOT/'data/audit/rya1353_line_census/v4/coverage_counts.csv'
    tracker_path=ROOT/'data/audit/element_status_tracker.csv'
    targets=json.loads(targets_path.read_text())['elements']
    census=list(csv.DictReader(census_path.open()))
    tracker=list(csv.DictReader(l for l in tracker_path.open() if not l.startswith('#')))
    inputs=[targets_path,census_path,tracker_path]
    rows=[]
    for target in targets:
        symbol=target['symbol']
        element,ion=('Fe','II') if symbol=='Fe II' else (symbol,'I')
        species=f'{element} {ion}'
        feed=ROOT/f'data/products/solar/{element}.json'
        products=[]
        if feed.exists():
            inputs.append(feed)
            products=[p for p in json.loads(feed.read_text())['products'] if p.get('ion')==ion]
        held=[r for r in census if r['species']==species]
        rows.append(dict(target=symbol,atomic_number=target['atomic_number'],priority=target['priority'],species=species,
                         active_feed_products=len(products),feed_state='ACTIVE_FEED_PRESENT' if products else 'NO_ACTIVE_FEED_PRODUCT_IN_THIS_SNAPSHOT',
                         processing_caveat='No active feed does not mean never processed; historical tracker is stale and cannot certify current survey completion',
                         historical_tracker=[{k:r[k] for k in ('ion','verdict','tier','refinement_debt')} for r in tracker if r['element']==element],
                         census_rows=sum(int(r['canonical_rows']) for r in held),
                         census_lab_label_rows=sum(int(r['canonical_rows']) for r in held if r['gf_tier']=='LAB'),
                         exposure='CANONICAL_AND_CENSUS_HOLDINGS_EXPOSED; LAB_LABEL_NOT_PRIMARY_VERIFICATION',
                         partition='Keep common upper levels, multiplets, shared lifetime/branching normalizations and experimental calibration scales together'))
    selected=[]
    for symbol,rationale in [('Li','Simplest canonical few-electron atom with no active feed; isotope/HFS and existing source exposure require explicit treatment'),
                             ('P','Priority-1 CHNOPS target with no active feed; primary experimental availability and many-electron feasibility require curator assessment'),
                             ('S','Priority-1 CHNOPS target with no active feed and oxygen-related valence structure; primary source/calibration independence remains unverified')]:
        row,=[r for r in rows if r['target']==symbol]
        if row['active_feed_products']:
            raise ValueError('Prospective no-feed shortlist premise changed')
        selected.append(dict(target=symbol,rationale=rationale,status='CANDIDATE_ONLY_NO_FROZEN_BENCHMARK',
                             primary_experiments_verified=False,blind=False))
    result=dict(schema='rya1354.prospective-canonical-screen.v2',scientific_admission='HOLD',
                inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(inputs))},
                rows=rows,selected=selected,
                previous_shortlist_note='Na and Mg have active products and do not serve the user-requested unfinished-feed branch; retain them only as exposed methodological controls.',
                next_gate='Distinct curator/human review must verify primary experiments, exposure and dependence groups, then freeze experiment/analysis and predictions before reveal. No genuine blindness is established.')
    (destination or ROOT/'data/audit/rya1354/source_curation_v2/prospective_screen.json').write_text(json.dumps(result,indent=2)+'\n')
    print([(r['target'],r['feed_state']) for r in rows if r['target'] in ('Li','P','S','Na','Mg')])


if __name__=='__main__':
    screen()
