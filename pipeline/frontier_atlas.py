"""Read-only, provenance-preserving inventory builder for the Frontier atlas.

No spectra, solver or abundance products are evaluated. Axis envelopes summarize
stored nodes and never authorize interpolation or a model/line combination.
"""
from __future__ import annotations
import copy,csv,gzip,hashlib,json,math,subprocess
from datetime import date
from pathlib import Path
from urllib.parse import urlsplit

ROOT=Path(__file__).resolve().parents[1]
SEED=Path('config/frontier/atlas.v1.json')
SCHEMA=Path('data/schemas/spectroscopic_frontier.schema.json')

def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))

def safe_path(root: Path, relative: str) -> Path:
    path=(root/relative).resolve()
    if path!=root.resolve() and root.resolve() not in path.parents:
        raise ValueError(f'Path outside repository: {relative}')
    if not path.is_file():raise ValueError(f'Missing source file: {relative}')
    return path

def axis(name,values,unit,source):
    if not values or not all(math.isfinite(v) for v in values):raise ValueError(f'Invalid axis {name}')
    return dict(name=name,min=min(values),max=max(values),unit=unit,source_id=source)

def derive(entry,sources,root):
    derivation=entry.pop('derivation')
    entry['coverage']={'transitions':[],'wavelength_air_nm':None,'bands':[],'atom_reach':None,'molecular_manifest':None,'source_labels':[],'labels_without_component_wavelength':[],'wavelength_native':None}
    if derivation is None:return
    source=derivation['source'];path=safe_path(root,sources[source]['path']);kind=derivation['type']
    if kind=='cno_grid':
        axes=[set() for _ in range(5)];lines=set()
        with gzip.open(path,'rt',encoding='utf-8') as stream:
            for line in stream:
                fields=line.split()
                if len(fields)!=7:raise ValueError(f'Invalid CDS row in {source}')
                for i in range(5):axes[i].add(float(fields[i]))
                lines.add(fields[5])
        entry['domain']['axes']=[axis(n,vs,u,source) for n,vs,u in zip(['teff','logg','feh','microturbulence','absolute_abundance'],axes,['K','log10(cm/s2)','dex','km/s','dex'])]
        selected=[]
        for line in safe_path(root,sources['cno-lines']['path']).read_text().splitlines():
            # Fixed-width schema is defined by the retained primary CDS ReadMe.
            if line[0:4].strip()!=derivation['species']:continue
            label=line[6:14].strip()
            if label in lines:selected.append({'source_label':label,'air_nm':float(line[17:25]),'vacuum_nm':float(line[28:36]),'source_id':'cno-lines'})
        if not selected:raise ValueError(f'No identified transitions for {entry["id"]}')
        entry['coverage']['source_labels']=sorted(lines)
        entry['coverage']['labels_without_component_wavelength']=sorted(lines-{r['source_label'] for r in selected})
        entry['coverage']['transitions']=selected
        wavelengths=[r['air_nm'] for r in selected]
        entry['coverage']['wavelength_air_nm']={'min':min(wavelengths),'max':max(wavelengths),'meaning':'sampled transition envelope, not continuous coverage'}
    elif kind=='fe_training':
        with path.open(encoding='utf-8') as stream:
            rows=[{**r,'source_row':str(i+2)} for i,r in enumerate(csv.DictReader(stream)) if r['ion']==derivation['ion']]
        if not rows:raise ValueError('Missing Fe training ion')
        entry['coverage']['transitions']=[{'source_label':r['network']+':CSV-row-'+r['source_row'],'air_nm':float(r['wavelength_air_A'])/10,'vacuum_nm':None,'source_id':source} for r in rows]
        vals=[r['air_nm'] for r in entry['coverage']['transitions']]
        entry['coverage']['wavelength_air_nm']={'min':min(vals),'max':max(vals),'meaning':'training-line envelope only; line and stellar-domain guards remain mandatory'}
        entry['domain']['axes']=[axis(name,{float(r[name]) for r in rows},unit,source) for name,unit in [('elo_eV','eV'),('eup_eV','eV'),('delta_E_eV','eV'),('loggf','dex')]]
    elif kind=='atom_reach':
        atom=read_json(path);stage=next(r for r in atom['stages'] if r['ion_stage']==derivation['ion_stage'])
        entry['coverage']['atom_reach']={'atom':atom['atom'],'sha256':atom['sha256'],'ion_stage':stage['ion_stage'],'n_levels':stage['n_levels'],'n_bound_bound':stage['n_bb_both'],'nlte_capable':stage['nlte_capable']}
    elif kind=='molecular_manifest':
        row=read_json(path)['molecules'][derivation['key']]
        entry['species']['isotopologues']=row.get('isotopologues',[])
        entry['coverage']['molecular_manifest']={k:row[k] for k in ['species','source','doi','version','origin','wavelength_coverage'] if k in row}
        bounds=row.get('wavelength_coverage')
        if bounds:
            entry['coverage']['wavelength_native']={'min':bounds['min_A'],'max':bounds['max_A'],'unit':'Angstrom','medium':'not established by this manifest','source_id':source}
            for name,lower,upper in [('UV',0,4000),('VIS',4000,7800),('IR',7800,1e13)]:
                if bounds['min_A'] < upper and bounds['max_A'] >= lower:entry['coverage']['bands'].append(name)
        # Preserve legacy source convention; do not relabel undocumented media as air.
        entry['domain']['scope']+=' Legacy manifest wavelength medium/frame are not inferred.'
    else:raise ValueError(f'Unknown derivation {kind}')
    envelope=entry['coverage']['wavelength_air_nm']
    if envelope:
        # Display bands are broad UI categories, not scientific validity boundaries.
        lo,hi=envelope['min'],envelope['max']
        entry['coverage']['bands']=[]
        for name,lower,upper in [('UV',0,400),('VIS',400,780),('IR',780,1e12)]:
            if lo < upper and hi >= lower:entry['coverage']['bands'].append(name)

def build(root=ROOT):
    seed=read_json(root/SEED);data=copy.deepcopy(seed);data['schema']='codex.spectroscopic_frontier/1';data['release_status']='DRAFT_SCIENTIFIC_REVIEW_REQUIRED';data['public_html_status']='DEFERRED_UNTIL_RYA_1370_APPROVED'
    with (root/'data/refs/bibliography.csv').open(encoding='utf-8') as stream:
        bib=list(csv.DictReader(stream))
    data['citations']={}
    for key,source in data['sources'].items():
        if 'path' in source:
            path=safe_path(root,source['path']);source['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            source['evidence_url']=f'https://github.com/damienabraxas/exoplanetcodex/blob/{data["science_base"]}/{source["path"]}'
            pinned=subprocess.check_output(['git','show',f'{data["science_base"]}:{source["path"]}'],cwd=root)
            if hashlib.sha256(pinned).hexdigest()!=source['sha256']:raise ValueError(f'Input differs from pinned Git evidence: {key}')
        else:source['evidence_url']=source['url']
        if source.get('bibliography_key'):
            bkey=source['bibliography_key'];row=next((r for r in bib if r['key']==bkey),None)
            if row is None:raise ValueError(f'Missing bibliography key {bkey}')
            data['citations'][bkey]={k:row[k] for k in ['authors','year','title','doi','url']}
    data['sources']['bibliography']={'path':'data/refs/bibliography.csv','sha256':hashlib.sha256((root/'data/refs/bibliography.csv').read_bytes()).hexdigest(),'evidence_url':f'https://github.com/damienabraxas/exoplanetcodex/blob/{data["science_base"]}/data/refs/bibliography.csv'}
    for entry in data['entries']:
        derive(entry,data['sources'],root)
        entry['upstream_citation_refs']=sorted({data['sources'][s]['bibliography_key'] for s in entry['source_ids'] if data['sources'][s].get('bibliography_key')})
    validate(data,root)
    return data

def validate(data,root=ROOT,check_hashes=True):
    import jsonschema
    jsonschema.Draft202012Validator(read_json(root/SCHEMA),format_checker=jsonschema.FormatChecker()).validate(data)
    sources=data['sources'];ids=[e['id'] for e in data['entries']]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate entry IDs')
    gap_ids={g['id'] for g in data['gaps']}
    entries_by_id={e['id']:e for e in data['entries']}
    canonical={r['symbol'].split()[0] for r in read_json(root/'data/config/elements_master.json')['elements']}
    for source in sources.values():
        if 'path' in source and check_hashes:
            if hashlib.sha256(safe_path(root,source['path']).read_bytes()).hexdigest()!=source['sha256']:raise ValueError('Source hash drift')
        parsed=urlsplit(source['evidence_url'])
        if parsed.scheme!='https' or parsed.username or parsed.password:raise ValueError('Unsafe evidence URL')
    for entry in data['entries']:
        if not set(entry['source_ids'])<=sources.keys():raise ValueError('Dangling source')
        if not set(entry['gap_ids'])<=gap_ids:raise ValueError('Dangling gap')
        if not set(entry['upstream_citation_refs'])<=data['citations'].keys():raise ValueError('Dangling citation')
        if entry['execution_approved']:raise ValueError('Atlas inventory cannot authorize execution')
        if entry['coverage_status']=='validated':raise ValueError('Draft atlas has no independent scientific sign-off')
        if entry['coverage_status']=='not_published' and not entry['literature_audit']['complete']:raise ValueError('Absence claim requires scoped completed search')
        if entry['literature_audit']['date'] and date.fromisoformat(entry['literature_audit']['date'])>date.fromisoformat(data['updated']):raise ValueError('Future audit date')
        if entry['component_type']=='line_list' and entry['treatment']!='NOT_APPLICABLE':raise ValueError('Line list is not an LTE/NLTE solver')
        species=entry['species']
        if species['kind']=='atom' and (species['element'] not in canonical or not species['ion'] or species['molecule'] is not None):raise ValueError('Invalid canonical atomic species')
        if species['kind']=='molecule' and (not species['molecule'] or species['element'] is not None or species['ion'] is not None):raise ValueError('Invalid molecular species')
        for key,kind in [('solver_ids','solver'),('atmosphere_ids','atmosphere_grid'),('model_atom_ids','model_atom')]:
            for ref in entry['components'][key]:
                if ref not in entries_by_id or entries_by_id[ref]['component_type']!=kind:raise ValueError('Invalid component reference/type')
        reach=entry['coverage']['atom_reach']
        if reach and not reach['nlte_capable'] and entry['coverage_status']!='incompatible':raise ValueError('Unsupported ion promoted')
        for ax in entry['domain']['axes']:
            if ax['min']>ax['max'] or not math.isfinite(ax['min']) or not math.isfinite(ax['max']):raise ValueError('Invalid domain bounds')
            if ax['source_id'] not in entry['source_ids']:raise ValueError('Untraced axis')
        for tr in entry['coverage']['transitions']:
            if tr['source_id'] not in entry['source_ids']:raise ValueError('Untraced transition')
    for gap in data['gaps']:
        if not set(gap['entry_ids'])<=set(ids) or not set(gap['source_ids'])<=sources.keys():raise ValueError('Dangling gap reference')
    table4=next(e for e in data['entries'] if e['id']=='fe-ii-amarsi2019-table4')
    if table4['treatment']!='3D_LTE':raise ValueError('Primary table4 is 3D LTE, not NLTE')
    return True

def matrix(data):
    lines=['# Spectroscopic Frontier coverage snapshot','',f'Draft, inspected {data["updated"]}; no scientific execution/adoption authorized.','', '| Entry | Species | Component | Treatment | Availability | Coverage / validation |','|---|---|---|---|---|---|']
    for e in data['entries']:
        s=e['species'];species=s['molecule'] or ' '.join(x for x in [s['element'],s['ion']] if x) or 'not species specific'
        lines.append('| '+' | '.join([e['id'],species,e['component_type'],e['treatment'],e['availability'],e['coverage_status']])+' |')
    lines+=['','Numeric envelopes, exact transitions, source hashes, citations, limits and gaps are','in the generated JSON. An envelope is not continuous coverage or an interpolation hull.','']
    return '\n'.join(lines)
