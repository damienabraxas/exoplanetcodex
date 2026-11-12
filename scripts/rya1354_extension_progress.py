"""Replay an asynchronous n5 log checkpoint; never qualify the running state."""
import argparse
import hashlib
import json
from pathlib import Path
import re

from rya1354_native_parts import load_parts
from rya1354_frozen_extend import controls
from rya1354_radial_controls import initialization_input
from rya1354_scf_history import history
from rya1354_wave_probe import decode

ROOT=Path(__file__).resolve().parents[1]
RUN='runs/lifetime_frozen_prior_n5_v1/'
LOG='evidence/n5_ground_logs_v14/'


def summarize(base):
    files,manifest=load_parts(base/'n5_ground_snapshot_v14')
    initial,_=load_parts(base/'frozen_initial_v13')
    target,_=load_parts(base/'lifetime_core_n5_v4')
    source,_=load_parts(base/'frozen_complete_v12')
    capture=json.loads(files[LOG+'capture.json'])
    for name,record in capture['files'].items():
        data=files[LOG+name]
        if len(data)!=record['bytes'] or hashlib.sha256(data).hexdigest()!=record['sha256']:
            raise ValueError('Captured native log hash differs')
    rows=history(files[LOG+'rmcdhf.stdout'].decode())
    if rows!=capture['history']:
        raise ValueError('Captured cycle inventory differs')
    stage=RUN+'outer_n5_ground/'
    seed=RUN+'initial_ground/'
    for name in ('rcsf.inp','isodata','rwfn.inp'):
        if files[stage+name]!=files[stage+'fixed_load/'+name]:
            raise ValueError('Native loading inputs differ')
    if files[stage+'rwfn.inp']!=initial[seed+'state.w']:
        raise ValueError('N5 ground initializer differs')
    if (initial[seed+'previous.w']!=source['runs/lifetime_frozen_prior_n4_v1/ground/state.w']
            or initial[seed+'rwfnestimate.stdin'].decode()!=initialization_input(3,439,.03125)):
        raise ValueError('N5 source wave/initialization differs')
    if files[LOG+'rmcdhf.stdin'].decode()!=controls('ground'):
        raise ValueError('N5 optimization controls differ')
    baseline=stage+'fixed_load/'
    if files[baseline+'rmcdhf.stdin'].decode()!=controls('ground',True,True):
        raise ValueError('N5 loading controls differ')
    if 'All subshell radial wavefunctions are fixed' not in (files[baseline+'rmcdhf.stdout']+files[baseline+'rmcdhf.stderr']).decode():
        raise ValueError('Native fixed-load completion marker absent')
    summary=files[baseline+'rmcdhf.sum'].decode()
    if (float(re.search(r'^ H\s*=\s*([\d.D+\-]+)',summary,re.M)[1].replace('D','E'))!=.03125
            or int(re.search(r'^ N\s*=\s*(\d+)',summary,re.M)[1])!=439):
        raise ValueError('Native loading grid differs')
    for name,selected in (('rcsf.inp','state.c'),('isodata','isodata')):
        if files[stage+name]!=target['runs/lifetime_core_n5_v1/CV_2s_DC_ground/'+selected]:
            raise ValueError('N5 target CSF/nucleus differs')
    inventory=decode(files[baseline+'rwfn.out'])
    if len(inventory)!=25 or sum(k[0]<=4 for k in inventory)!=16:
        raise ValueError('Native loading inventory differs')
    complete=[r for r in rows if r['complete_printed_cycle']]
    return dict(schema='rya1354.n5-running-checkpoint.v14',captured_utc=capture['captured_utc'],
        last_complete_cycle=complete[-1],completed_cycle_count=len(complete),
        original_ACC=.0625**6,baseline_orbitals=25,declared_frozen_records=16,
        loading_control_verified=True,source_initializer_and_target_verified=True,
        scientific_admission='HOLD',adoption_allowed=False,accepted_lifetime=None,
        optimization_qualified=False,exposure='EXPOSED_NONBLIND',uncertainty=None,covariance=None,
        calibrated_sigma=None,parts=manifest['parts'],
        scope='ASYNCHRONOUS_LOG_PREFIX_AND_COMPLETED_NATIVE_LOADING_CONTROL_ONLY',
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
