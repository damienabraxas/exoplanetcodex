#!/usr/bin/env python3
"""RYA-1268 -- the per-line gf reference code of every line in our solar VALD extracts.

VALD long format follows each data line with two term lines and one reference line
("_  ... wl:XX  gf:CODE  ..."). The CODE names the gf's source; we hold no VALD BibTeX, so a
code is resolved to a paper only where the Asplund-cited tables name the gf source of a line
that VALD carries with that code (scripts/rya1268_coverage_audit.py). Output:
data/audit/rya1268_line_coverage/vald_solar_gf_refs.csv.
Run with CODEX_VALD_RAW set (Mac: config/hosts/mac.env).
"""
from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "data" / "audit" / "rya1268_line_coverage" / "vald_solar_gf_refs.csv"
FILES = ["vald_solar_nearuv_2000_3780_hfson_raw.txt", "vald_solar_raw.txt",
         "vald_solar_redopt_6910_9500_hfson_raw.txt", "vald_solar_ir_9500_17000_hfson_raw.txt",
         "vald_solar_ir_17000_25000_hfson_raw.txt"]
DATA = re.compile(r"^'([A-Z][a-z]?) (\d+)',\s+([\d.]+),\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+),"
                  r"\s*([-\d.]+),\s*([-\d.]+),")
GF = re.compile(r"\bgf:\s*([^\s]+)")
ROMAN = {1: "I", 2: "II", 3: "III"}


def main() -> None:
    from config.constants import codex_path
    raw = Path(os.environ.get("CODEX_VALD_RAW") or codex_path("data.vald_raw"))
    rows = []
    for name in FILES:
        lines = (raw / name).read_text(errors="replace").splitlines()
        for i, l in enumerate(lines):
            m = DATA.match(l)
            if not m or int(m.group(2)) not in ROMAN:
                continue
            depth = l.rstrip().rstrip(",").split(",")[-1].strip()
            ref = lines[i + 3] if i + 3 < len(lines) else ""
            g = GF.search(ref)
            rows.append({"species": f"{m.group(1)} {ROMAN[int(m.group(2))]}",
                         "wavelength_air_A": float(m.group(3)), "loggf": float(m.group(4)),
                         "ep_eV": float(m.group(5)), "j_lower": float(m.group(6)),
                         "ep_upper_eV": float(m.group(7)), "j_upper": float(m.group(8)),
                         "central_depth": depth,
                         "gf_ref_code": g.group(1) if g else "", "file": name})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(len(rows), "lines ->", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
