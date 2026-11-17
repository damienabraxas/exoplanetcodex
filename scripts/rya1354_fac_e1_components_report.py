"""Portable v44 report: re-verify every component context and decompose it.

Decodes the committed v44 archive (its v41 references restored and hash-checked)
and the v43 archive; writes a machine report. Diagnostics only: HOLD /
EXPOSED_NONBLIND, nothing adopted.
"""
import argparse
import json
import math
from pathlib import Path

from rya1354_fac_e1_components_audit import parse,verify,decompose,grid_context,configuration_shift
from rya1354_native_parts import load_parts

RUN='runs/fac_e1_components_v44/'

CASES=('ground_coherent10','ground_coherent100','equal_coherent10','equal_coherent100')


def native_amplitudes(raw_text):
    out={}
    for line in raw_text.splitlines():
        if line.startswith('#') or not line.strip():continue
        f=line.split()
        if int(f[5])==-1:out[int(f[1])]=dict(amplitude=float(f[12]),energy_au=float(f[11]),A_s1=float(f[14]))
    return out


def total_rate(raw_text):
    """Sum of every native channel rate from the upper level (all multipoles), s^-1."""
    return math.fsum(float(l.split()[14]) for l in raw_text.splitlines() if l.strip() and not l.startswith('#'))


def channels_for(files,prefix,name,top):
    raw=native_amplitudes(files[prefix+name+'/raw.tsv'].decode())
    header,contexts=parse(files[prefix+name+'/components.tsv'].decode());cut=float(header['ANGZCUT'])
    channels={}
    for c in contexts:
        n=raw[c['lower']]
        if c['energy_au']!=n['energy_au']:raise ValueError('Component energy differs from native raw energy')
        v=verify(c,n['amplitude'],cut);d=decompose(c,top)
        channels[str(c['lower'])]=dict(verification=v,grid=grid_context(c),decomposition=d,A_s1=n['A_s1'],energy_au=n['energy_au'])
    return channels,contexts


def case_report(files,name,v43,top):
    channels,contexts=channels_for(files,RUN,name,top)
    ref=contexts[0]
    if ref['lower']!=0:raise ValueError('Reference 3P2 context is not first')
    rows={((r['n1'],r['kappa1']),(r['n2'],r['kappa2'])):r for r in v43['cases'][name+'_own']['rows']}
    radial=[]
    for a in ref['A']:
        r=rows[a['pair']]
        if r['used_au']!=ref['energy_au']:raise ValueError('v43 own energy differs from reference transition energy')
        radial.append(dict(pair=[list(a['pair'][0]),list(a['pair'][1])],native_grid=a['radial'],v43_exact=r['radial_FR'],
            relative_difference=(a['radial']-r['radial_FR'])/r['radial_FR'] if r['radial_FR'] else None,
            term_shift=a['coeff']*(r['radial_FR']-a['radial'])))
    return dict(channels=channels,reference_radial_vs_v43=dict(rows=radial,
        amplitude_shift_if_exact=math.fsum(x['term_shift'] for x in radial),
        max_abs_relative=max(abs(x['relative_difference']) for x in radial if x['relative_difference'] is not None))),\
        {c['lower']:channels[str(c['lower'])]['decomposition'] for c in contexts}


def summarize(base,repo_root,top=25):
    files,manifest=load_parts(base/'fac_e1_components_complete_v44',repo_root=repo_root)
    v43_files,_=load_parts(base/'fac_radial_operator_complete_v43')
    receipt=json.loads(files[RUN+'receipt.json'])
    if receipt['status']!='FAC_E1_COMPONENTS_COMPLETED_REVIEW_REQUIRED':raise ValueError('Run not completed')
    v43=json.loads(v43_files['runs/fac_radial_operator_v43/receipt.json'])
    report=dict(schema='rya1354.e1-components-report.v44',run_receipt_status=receipt['status'],index_sha256=receipt['index_sha256'],
        archive=dict(stored_files=len(manifest['files']),referenced_files=len(manifest['referenced_files']),
                     native_paths=manifest['native_paths_preserved_byte_for_byte'],parts=len(manifest['parts']),
                     compressed_bytes=sum(v['bytes'] for v in manifest['parts'].values())),
        canary=receipt['canary']['byte_identical'] and len(receipt['canary']['byte_identical']),
        cases={},shifts={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    dec={}
    for case in CASES:
        for s in ('B','C'):
            name=case+'_'+s;report['cases'][name],dec[name]=case_report(files,name,v43,top)
            for low,ch in report['cases'][name]['channels'].items():
                want=dict(receipt['cases'][name]['components'][low]);grid=want.pop('grid')
                if ch['verification']!=want or ch['grid']!=grid:
                    raise ValueError('Portable verification differs from native receipt: '+name+'/'+low)
            report['cases'][name]['byte_identical_v41_files']=len(receipt['cases'][name]['compared_byte_identical'])
    for scale in ('coherent10','coherent100'):
        for s in ('B','C'):
            report['shifts'][f'ground_to_equal_{scale}_{s}']=configuration_shift(dec[f'ground_{scale}_{s}'][0],dec[f'equal_{scale}_{s}'][0])
    for case in CASES:
        report['shifts'][f'{case}_B_to_C']=configuration_shift(dec[case+'_B'][0],dec[case+'_C'][0])
    return report


V46='runs/fac_mixcut_control_v46/'


def summarize_mixcut(base,repo_root,top=25):
    """Full-eigenvector control vs the v44 (default mix-cut) decomposition of the same models."""
    files,manifest=load_parts(base/'fac_mixcut_control_complete_v46',repo_root=repo_root)
    failed,_=load_parts(base/'fac_mixcut_control_failed_v45',repo_root=repo_root)
    retained=json.loads(failed['runs/fac_mixcut_control_v45/receipt.json'])
    v44,_=load_parts(base/'fac_e1_components_complete_v44',repo_root=repo_root)
    receipt=json.loads(files[V46+'receipt.json'])
    if receipt['status']!='FAC_MIXCUT_CONTROL_COMPLETED_REVIEW_REQUIRED':raise ValueError('v46 run not completed')
    if retained['status']!='FAILED_RETAINED':raise ValueError('Retained v45 failure differs')
    report=dict(schema='rya1354.mixcut-control-report.v46',run_receipt_status=receipt['status'],index_sha256=receipt['index_sha256'],
        retained_failure=dict(status=retained['status'],error=retained['error'],canary_files=len(retained['canary']['byte_identical']),
            completed_native_case_reproduced=sorted(k for k,v in receipt['cases'].items() if 'retained_failure_outputs_identical' in v['compared_byte_identical'])),
        archive=dict(stored_files=len(manifest['files']),referenced_files=len(manifest['referenced_files']),
                     native_paths=manifest['native_paths_preserved_byte_for_byte'],parts=len(manifest['parts']),
                     compressed_bytes=sum(v['bytes'] for v in manifest['parts'].values())),
        canary_files=len(receipt['canary']['byte_identical']),cases={},
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    for name,entry in receipt['cases'].items():
        new,_=channels_for(files,V46,name,top);old,_=channels_for(v44,RUN,name,top)
        for low,ch in new.items():
            want=dict(entry['components']['contexts'][low]);grid=want.pop('grid')
            if ch['verification']!=want or ch['grid']!=grid:raise ValueError('Portable verification differs from v46 receipt: '+name+'/'+low)
        rows={}
        for low in new:
            a,b=old[low],new[low]
            rows[low]=dict(v44_amplitude=a['verification']['native_amplitude'],full_mixing_amplitude=b['verification']['native_amplitude'],
                relative_amplitude_change=b['verification']['native_amplitude']/a['verification']['native_amplitude']-1,
                v44_A_s1=a['A_s1'],full_mixing_A_s1=b['A_s1'],rate_ratio=b['A_s1']/a['A_s1'],
                stored_basis=dict(v44=a['decomposition']['stored_basis'],full_mixing=b['decomposition']['stored_basis']),
                cut=dict(v44=a['decomposition']['cut'],full_mixing=b['decomposition']['cut']),
                cancellation_index=dict(v44=a['decomposition']['cancellation_index'],full_mixing=b['decomposition']['cancellation_index']),
                configuration_shift=configuration_shift(a['decomposition'],b['decomposition']))
            rows[low]['configuration_shift']['rows']=rows[low]['configuration_shift']['rows'][:top]
        report['cases'][name]=dict(channels=rows,byte_identical_v41_files=len(entry['compared_byte_identical']['byte_identical']),
            transport=entry['compared_byte_identical']['transport'],
            total_rate_s1=dict(v44=total_rate(v44[RUN+name+'/raw.tsv'].decode()),full_mixing=total_rate(files[V46+name+'/raw.tsv'].decode())))
        tr=report['cases'][name]['total_rate_s1'];tr['ratio']=tr['full_mixing']/tr['v44']
    return report


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--base',type=Path,required=True);ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--version',choices=('v44','v46'),default='v44')
    a=ap.parse_args()
    report=summarize(a.base,a.repo_root) if a.version=='v44' else summarize_mixcut(a.base,a.repo_root)
    a.out.write_text(json.dumps(report,indent=1)+'\n')


if __name__=='__main__':main()
