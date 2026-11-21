"""Replay all previous full-space transport results before larger native CI."""
import argparse
import json
from pathlib import Path
from rya1354_fac_stream_transport import transport


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    receipt=json.loads((a.baseline/'receipt.json').read_text());results={}
    for name,item in receipt['cases'].items():
        p=a.baseline/name;result=transport(p/'o.basis',p/'o_export_000.c',p/'o_export_000.cm')
        if result!=item['transport']:raise ValueError('Streaming replay differs: '+name)
        results[name]=result
    a.out.write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))

if __name__=='__main__':main()
