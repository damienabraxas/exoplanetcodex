#!/usr/bin/env python3
"""Build/check review-only atlas JSON and Markdown without scientific execution."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pipeline.frontier_atlas import build,matrix
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
data=build();outputs={ROOT/'data/products/frontier/coverage.v1.json':json.dumps(data,indent=2,ensure_ascii=False)+'\n',ROOT/'docs/science/spectroscopic_frontier_matrix.md':matrix(data)}
for path,text in outputs.items():
 if args.check:
  if not path.is_file() or path.read_text()!=text:raise SystemExit(f'Stale atlas output: {path.relative_to(ROOT)}')
 else:path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
print('Frontier atlas schema, references, source hashes and generated snapshots: PASS')
