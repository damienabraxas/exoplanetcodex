#!/usr/bin/env python3
"""Engine-A IR capability, proved per element — RYA-713.

Ryan: "i want a double check though on every element for Engine A, and make sure
there is no IR capabilities."

WHY THIS EXISTS
---------------
"No IR NLTE" has been asserted four times and was wrong every time, always the same
way: a derived per-line CSV was read as the model's limit. It is not. The Amarsi/
Mallinson departure grids are indexed by

    (MARCS atmosphere index, depth index, [X/Fe], LEVEL index) -> b = n_NLTE / n_LTE

There is NO WAVELENGTH IN THE GRID. Any transition between two levels of the model
atom has departure coefficients, IR included. `Ti_Bergemann2011_MPIA.csv` stopping at
6861 A is a statement about that CSV and nothing else.

So this tool never infers capability from an extract. It reports the three things an
Engine-A correction actually requires and names which one is missing:

    1. the .grd            -- departure coefficients per level (wavelength-free)
    2. label_{El}.txt      -- level index -> species/config/term/J/energy, so a line's
                              upper and lower levels can be FOUND in the grid
    3. a level-identified linelist for the line (VALD wavelength + log gf is NOT
                              enough; the vendor warns VALD labels differ from the
                              NIST labels the grids are mapped to)

Miss any one and the line silently falls back to LTE.

THE DISTINCTION THIS TOOL REFUSES TO COLLAPSE
---------------------------------------------
"has no IR capability" and "is blocked on an unextracted label file" are not the same
statement and must never share a word. One is physics. The other is a `tar -x` that
was skipped. Reporting them identically is what produced the four wrong claims.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXTRACTS = ROOT / "data" / "nlte_grids"
OUT = ROOT / "data" / "audit" / "engine_a_ir"

SIRIUS = "sirius"
GRID_DIR = "/mnt/codex-data/grids/nlte/amarsi_galah"
GES = ("/mnt/codex-data/engines/Turbospectrum_NLTE/COM/linelists/"
       "nlte_ges_linelist_jmg17feb2022_I_II")

# Past this, HARPS cannot see the line at all -- the boundary that makes a span "IR".
HARPS_RED_EDGE_A = 6910.0

Z = {"H": 1, "Li": 3, "C": 6, "N": 7, "O": 8, "Na": 11, "Mg": 12, "Al": 13, "Si": 14,
     "P": 15, "S": 16, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24,
     "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Sr": 38, "Y": 39,
     "Zr": 40, "Ba": 56, "Eu": 63}

VERDICTS = ("IR-CAPABLE-NOW", "IR-CAPABLE-BLOCKED-ON-LABEL",
            "IR-CAPABLE-BLOCKED-ON-LINELIST", "PER-LINE-TABLE-ONLY", "NO-ENGINE-A")

# There are TWO Engine-A families and they are NOT interchangeable in what they let us do:
#
#   Amarsi / Mallinson (PySME .grd)  -- ships the LEVEL grid. Wavelength-free, so a new
#       line at any wavelength is a local computation once label_{El}.txt is present.
#
#   Bergemann MPIA / INSPECT         -- serves per-LINE corrections from a web tool
#       (nlte.mpia.de). What we hold is the delivered product, not an extract we chose to
#       truncate. Reaching a new wavelength needs a NEW UPSTREAM QUERY, not local compute.
#
# Calling the second one "no Engine A" would be the same collapse this tool exists to
# prevent -- Fe has 26,560 MPIA lines and is emphatically not un-modelled. But calling it
# level-extensible would be the mirror error. It gets its own verdict.
PERLINE_FAMILY_HINTS = ("MPIA", "INSPECT", "Korotin", "Mashonkina")


def _ssh(script: str) -> str:
    """Run a snippet on Sirius. Grids are Sirius-only and are never copied to the Mac."""
    r = subprocess.run(["ssh", "-o", "ConnectTimeout=30", SIRIUS, script],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"sirius unreachable ({r.returncode}); this audit cannot be "
                         f"faked from the Mac -- the grids only exist there.\n{r.stderr[:400]}")
    return r.stdout


def sirius_assets() -> tuple[dict[str, str], set[str]]:
    """Which .grd and which label_*.txt are actually staged."""
    out = _ssh(f"ls {GRID_DIR}/*.grd 2>/dev/null | xargs -r -n1 basename; "
               f"echo '###'; ls {GRID_DIR}/label_*.txt 2>/dev/null | xargs -r -n1 basename")
    grd_part, _, lab_part = out.partition("###")
    grds: dict[str, str] = {}
    for name in grd_part.split():
        # nlte_{El}_scatt_pysme.grd / nlte_{El}_pysme.grd / nlte_{El}_caliskan_..._pysme.grd
        m = re.match(r"nlte_([A-Z][a-z]?)_", name)
        if m:
            grds[m.group(1)] = name
    labels = {m.group(1) for n in lab_part.split()
              if (m := re.match(r"label_([A-Z][a-z]?)\.txt$", n))}
    return grds, labels


def ges_spans() -> dict[tuple[int, int], dict]:
    """Level-identified linelist reach, per (Z, stage). This is the linelist that
    carries the level identification NLTE needs -- not a plain VALD list."""
    script = f"""python3 - <<'PY'
import re, json
hdr = re.compile(r"^'\\s*(\\d+)\\.\\d\\d\\d\\s*'\\s+(\\d+)\\s+(\\d+)")
cur, span = None, {{}}
for ln in open({GES!r}, errors='replace'):
    m = hdr.match(ln)
    if m:
        cur = (int(m.group(1)), int(m.group(2)))
        span.setdefault(cur, [1e9, -1e9, 0])
        continue
    if cur is None:
        continue
    try:
        w = float(ln.split()[0])
    except Exception:
        continue
    s = span[cur]
    s[0] = min(s[0], w); s[1] = max(s[1], w); s[2] += 1
print(json.dumps({{f'{{z}}_{{st}}': dict(n=n, lo=round(a, 1), hi=round(b, 1))
                  for (z, st), (a, b, n) in span.items()}}))
PY"""
    raw = _ssh(script).strip().splitlines()[-1]
    return {tuple(int(x) for x in k.split("_")): v for k, v in json.loads(raw).items()}


def extract_spans() -> dict[str, list[dict]]:
    """Per-line CSV spans. These are EXTRACTS. Recorded so the audit can show what
    was previously mistaken for a ceiling -- never used to decide capability."""
    out: dict[str, list[dict]] = {}
    for f in sorted(EXTRACTS.glob("*.csv")):
        m = re.match(r"([A-Z][a-z]?)_", f.name)
        if not m:
            continue
        try:
            d = pd.read_csv(f)
        except Exception:
            continue
        wcol = next((c for c in d.columns if "wave" in c.lower()), None)
        if wcol is None or not len(d):
            continue
        out.setdefault(m.group(1), []).append(
            dict(file=f.name, n=len(d), lo=float(d[wcol].min()), hi=float(d[wcol].max())))
    return out


def audit() -> pd.DataFrame:
    grds, labels = sirius_assets()
    ges = ges_spans()
    ext = extract_spans()

    rows = []
    for el in sorted(set(Z) | set(grds) | set(ext)):
        z = Z.get(el)
        g = grds.get(el)
        has_label = el in labels
        # Neutral + singly-ionised reach in the level-identified linelist.
        gs = [ges[(z, st)] for st in (1, 2) if z is not None and (z, st) in ges]
        ges_hi = max((s["hi"] for s in gs), default=None)
        ges_n = sum(s["n"] for s in gs)

        exts = ext.get(el, [])
        ex_hi = max((e["hi"] for e in exts), default=None)
        ex_str = "; ".join(f"{e['file']}:{e['lo']:.0f}-{e['hi']:.0f}({e['n']})" for e in exts)

        perline = [e for e in exts
                   if any(h.lower() in e["file"].lower() for h in PERLINE_FAMILY_HINTS)]

        if not g and perline:
            verdict = "PER-LINE-TABLE-ONLY"
            src = "; ".join(f"{e['file']} ({e['n']} lines, {e['lo']:.0f}-{e['hi']:.0f} A)"
                            for e in perline)
            blocker = (f"Engine A EXISTS as a per-line table [{src}] but no level .grd is "
                       f"staged. This family serves corrections per line from upstream, so "
                       f"reaching a new wavelength is an UPSTREAM QUERY, not local compute. "
                       f"Not 'no Engine A' and not level-extensible either.")
        elif not g:
            verdict, blocker = "NO-ENGINE-A", "no departure-coefficient .grd staged on Sirius"
        elif not has_label:
            verdict = "IR-CAPABLE-BLOCKED-ON-LABEL"
            blocker = (f"{g} is staged but label_{el}.txt is NOT -- level index cannot be "
                       f"resolved, so no line can be located in the grid. This is an "
                       f"extraction from the source tarball, NOT a re-download.")
        elif ges_hi is None or ges_hi <= HARPS_RED_EDGE_A:
            verdict = "IR-CAPABLE-BLOCKED-ON-LINELIST"
            blocker = (f"grid+label staged, but no level-identified linelist reaches past "
                       f"{HARPS_RED_EDGE_A:.0f} A (GES hi={ges_hi})")
        else:
            verdict, blocker = "IR-CAPABLE-NOW", ""

        rows.append(dict(
            element=el, grd=g or "", label_staged=has_label,
            ges_lines=ges_n or 0, ges_hi_A=ges_hi,
            ges_reaches_ir=bool(ges_hi and ges_hi > HARPS_RED_EDGE_A),
            extract_hi_A=ex_hi,
            extract_already_ir=bool(ex_hi and ex_hi > HARPS_RED_EDGE_A),
            ir_verdict=verdict, blocker=blocker, extracts=ex_str))

    df = pd.DataFrame(rows)
    bad = set(df.ir_verdict) - set(VERDICTS)
    if bad:
        raise SystemExit(f"unknown verdict(s) {bad}")
    # A verdict that is not IR-CAPABLE-NOW must name its blocker. A bare "no" is the
    # failure mode this whole audit exists to prevent.
    silent = df[(df.ir_verdict != "IR-CAPABLE-NOW") & (df.blocker.str.strip() == "")]
    if len(silent):
        raise SystemExit(f"blocker unnamed for: {list(silent.element)}")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()

    df = audit()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "engine_a_ir_capability.csv", index=False)

    print("ENGINE-A IR CAPABILITY — RYA-713")
    print("The .grd is level-indexed. A per-line CSV's max wavelength is NOT a ceiling.\n")
    for v in VERDICTS:
        sub = df[df.ir_verdict == v]
        if not len(sub):
            continue
        print(f"  {v}  ({len(sub)})")
        for _, r in sub.iterrows():
            hi = f"{r.ges_hi_A:.0f}" if pd.notna(r.ges_hi_A) else "--"
            print(f"    {r.element:3s} linelist->{hi:>6s} A  {r.blocker[:96]}")
        print()

    ir_ext = df[df.extract_already_ir]
    if len(ir_ext):
        print("  Extracts that ALREADY reach past HARPS (the claim was falsified in-repo):")
        for _, r in ir_ext.iterrows():
            print(f"    {r.element:3s} {r.extracts}")
        print()
    print(f"  wrote {out / 'engine_a_ir_capability.csv'}")


if __name__ == "__main__":
    main()
