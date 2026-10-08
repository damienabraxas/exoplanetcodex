#!/usr/bin/env python3
"""RYA-1232 -- give the tracked VALD-built synthesis lists the damping their builder now writes.

    python3 scripts/rya1232_fill_unsold_damping.py [--check]

Our Turbospectrum DROPS an atomic line whose fdamp is 0.0 (it is not a "use the default"
value): C I 16890.38 never formed at any A(C), Mg I 15740.71 (45% deep) was absent from every
H synthesis. `pipeline.nearuv_linelist.to_ispec_array` now writes the Unsöld enhancement
GESv6 uses (`unsold_fdamp`) wherever VALD has no vdW. This applies the SAME mapping to the
lists already built through it -- only `turbospectrum_fdamp` on rows with waals == 0 and
fdamp == 0 changes; every other byte of each row is left as written.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.nearuv_linelist import unsold_fdamp  # noqa: E402

LISTS = ("ispec_ir_9200_13000", "ispec_h_15007_17494", "ispec_k_19452_24846")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="report only; exit 1 if any row is unfilled")
    a = ap.parse_args()
    bad = 0
    for name in LISTS:
        p = ROOT / "data/linelists" / name / "atomic_lines.tsv"
        lines = p.read_text().split("\n")
        hdr = lines[0].split("\t")
        ie, iw, ifd = hdr.index("element"), hdr.index("waals"), hdr.index("turbospectrum_fdamp")
        n = 0
        for k in range(1, len(lines)):
            if not lines[k]:
                continue
            f = lines[k].split("\t")
            if float(f[iw]) == 0.0 and float(f[ifd]) == 0.0:
                n += 1
                f[ifd] = f"{unsold_fdamp(f[ie]):g}"
                lines[k] = "\t".join(f)
        bad += n
        print(f"{name}: {n} row(s) with fdamp 0.0 -> Unsöld")
        if n and not a.check:
            p.write_text("\n".join(lines))
    return 1 if (a.check and bad) else 0


if __name__ == "__main__":
    raise SystemExit(main())
