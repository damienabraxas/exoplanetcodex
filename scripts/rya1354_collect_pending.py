"""Collect terminal task-owned campaigns; never launch physics or adopt results.

Runs on Sirius. A heartbeat and terminal receipt/bundle paths survive a chat
checkpoint. Completed and failed native campaigns are kept distinct. Missing or
incomplete receipts never become accepted predictions.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import time

from rya1354_pack_native import pack
from rya1354_solver_smoke import sha

CAMPAIGNS=('lifetime_relax_n4_v1','lifetime_mr_n4_v1')
TERMINAL={'CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED',
          'BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED','FAILED_RETAINED'}


def terminal_receipt(path):
    if not path.exists():
        return None
    try:
        result=json.loads(path.read_text())
    except json.JSONDecodeError:
        # Native finally is still writing. Preserve it and inspect next heartbeat.
        return None
    if result.get('status') not in TERMINAL:
        raise ValueError('Unexpected campaign receipt status: '+str(result.get('status')))
    return result


def collect(root,out):
    """One observation/collection pass; return statuses without inferred success."""
    result={}
    for name in CAMPAIGNS:
        source=root/'runs'/name
        receipt_path=source/'receipt.json'
        receipt=terminal_receipt(receipt_path)
        target=out/name
        if receipt is None:
            result[name]=dict(collection='WAITING_FOR_TERMINAL_RECEIPT',scientific_admission='HOLD')
            continue
        if not target.exists():
            pack(root,[Path('runs')/name],target,omit_rci_restart_cache=True)
        elif not (target/'manifest.json').is_file():
            raise ValueError('Incomplete earlier collection retained; use a new explicitly reviewed output directory')
        manifest=json.loads((target/'manifest.json').read_text())
        if manifest['files'][f'runs/{name}/receipt.json']['sha256']!=sha(receipt_path):
            raise ValueError('Terminal native receipt changed after collection')
        result[name]=dict(collection='ARCHIVED_TERMINAL_OUTPUTS',native_status=receipt['status'],
                          native_error=receipt.get('error'),native_receipt_sha256=sha(receipt_path),
                          manifest_sha256=sha(target/'manifest.json'),bundle=str(target.relative_to(root)),
                          scientific_admission='HOLD',adoption_allowed=False,
                          verification='Packer verified byte-identical archive replay; physical-state/lifetime review still required')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--watch',action='store_true')
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Collect native task outputs on Sirius')
    root=args.root.resolve();out=args.out.resolve()
    if not out.is_relative_to(root/'evidence'):
        raise ValueError('Collection output must belong to the task evidence directory')
    out.mkdir(parents=True,exist_ok=False)
    while True:
        states=collect(root,out)
        record=dict(schema='rya1354.pending-native-collection.v1',observed_utc=datetime.now(timezone.utc).isoformat(),
                    scientific_admission='HOLD',adoption_allowed=False,campaigns=states,
                    scope='NATIVE_RECEIPT_AND_ARCHIVE_COLLECTION_ONLY_NO_SCIENTIFIC_ACCEPTANCE')
        temporary=out/'status.next.json';temporary.write_text(json.dumps(record,indent=2)+'\n')
        temporary.replace(out/'status.json')
        if not args.watch or all(s['collection']=='ARCHIVED_TERMINAL_OUTPUTS' for s in states.values()):
            print(json.dumps(record));return
        time.sleep(30)


if __name__=='__main__':
    main()
