#!/usr/bin/env python3
"""RYA-1233 -- superseded by the generic builder: `python -m pipeline.line_sets [--element X]`.

The published sets are now DECLARED in data/reference/line_sets/SOURCES.yaml and built by
pipeline/line_sets.py (byte-identical output for every Si set). Kept as a thin entry point so
older notes that name this script still work.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.line_sets import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
