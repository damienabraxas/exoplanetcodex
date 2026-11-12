"""Reproduce terminal n5 evidence and native n4-to-n5 model diagnostics."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_extension_audit import summarize
from rya1354_frozen_summary import summarize as previous
from rya1354_layer_change import decompose

ROOT=Path(__file__).resolve().parents[1]


def report(base):
    current=summarize(base,base/'n5_complete_v15')
    old=previous(base)
    return dict(schema='rya1354.n5-terminal-and-layer-diagnostic.v15',n5=current,
        dominant_E1_change=decompose(old['decay'],current['decay']),
        lifetime_fractional_change={g:current['decay']['model'][g]['model_one_photon_lifetime_ns']/
            old['decay']['model'][g]['model_one_photon_lifetime_ns']-1 for g in ('B','C')},
        scientific_admission='HOLD',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        interpretation='Exploratory layerwise frozen-prior-layer model sensitivity; gauge/model differences are not calibrated uncertainty. Full release failure and historical n3 component limitation retained.',
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(report(args.base),indent=2)+'\n')
