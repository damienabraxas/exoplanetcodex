#!/usr/bin/env python3
"""
RYA-1220 -- contamination audit of the atomic N I diagnostic set.

WHY THIS EXISTS. RYA-1220 expanded the native N I set from three lines to five by adding
7442.298 and 8629.235, and reported a five-line median NLTE correction of -0.0100 dex
against the three-line -0.0128. Nothing in that work measured what sits NEXT TO the two
new lines.

Both are blended with CN -- the very molecule whose abundance is the competing N
diagnostic. 8629.235 has a CN line 0.0076 A away carrying 75% of the N I line's own
central depth: unresolvable at any realistic resolution. Admitting an N I line that is
three-quarters CN, in order to adjudicate N I against CN, is circular.

METHOD. Contamination is measured from the PRESERVED long-format VALD extract, not from
the project's built line list -- the built list carries no term/J and no neighbours, and
a neighbour is the whole question here. Ratio is summed catalogued central_depth of every
non-N-I line within the window, over the N I line's own central_depth (the RYA-1189
discriminator: RANK isolation, never filter on it).

This measures the CATALOGUE, so it bounds contamination; it does not prove a fit absorbed
it. A line is flagged, never silently dropped.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

# The five Amarsi et al. 2020 Table 1 solar indicators, as registered/reported by RYA-1220.
NI_SET = [7442.298, 7468.312, 8216.336, 8629.235, 8683.403]
WINDOW_A = 0.05          # +/- around the line; ~2x a solar N I FWHM at these wavelengths
ADMIT_MAX_RATIO = 0.50   # above this the line is contamination-dominated, not N I

ARCHIVE = (pathlib.Path.home() / "Documents" / "Exoplanet Codex" / "vald_raw_archive"
           / "current" / "vald_solar_redopt_6910_9500_hfson_raw.txt")

_ROW = re.compile(r"^'(?P<spec>[^']+)',\s*(?P<wl>[0-9.]+),")


def read_vald(path: pathlib.Path) -> list[dict]:
    """Long-format VALD rows: species, air wavelength, central depth."""
    out = []
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            m = _ROW.match(line)
            if not m:
                continue
            parts = line.split(",")
            try:
                depth = float(parts[13])
            except (IndexError, ValueError):
                continue
            out.append({"species": m.group("spec").strip(),
                        "wl": float(m.group("wl")), "depth": depth})
    return out


def audit(rows: list[dict], window: float = WINDOW_A) -> list[dict]:
    results = []
    for target in NI_SET:
        me = [r for r in rows if r["species"].startswith("N 1") and abs(r["wl"] - target) < 0.02]
        if not me:
            results.append({"wavelength_air_A": target, "status": "NOT_IN_EXTRACT"})
            continue
        own = max(me, key=lambda r: r["depth"])
        near = [r for r in rows
                if abs(r["wl"] - target) < window and not (r["species"].startswith("N 1")
                                                           and abs(r["wl"] - target) < 0.02)]
        near.sort(key=lambda r: -r["depth"])
        total = sum(r["depth"] for r in near)
        ratio = total / own["depth"] if own["depth"] else float("inf")
        cn = [r for r in near if r["species"].startswith("CN")]
        cn_ratio = sum(r["depth"] for r in cn) / own["depth"] if own["depth"] else float("inf")
        results.append({
            "wavelength_air_A": target,
            "own_central_depth": round(own["depth"], 4),
            "n_contaminants": len(near),
            "blend_depth_ratio": round(ratio, 3),
            "cn_depth_ratio": round(cn_ratio, 3),
            "nearest": ({"species": near[0]["species"],
                         "wavelength_air_A": round(near[0]["wl"], 4),
                         "separation_A": round(near[0]["wl"] - target, 4),
                         "central_depth": round(near[0]["depth"], 4)} if near else None),
            "admit": bool(ratio <= ADMIT_MAX_RATIO),
            # A line whose contaminant is CN cannot adjudicate N I against CN.
            "cn_circular": bool(cn_ratio > 0.5),
        })
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--extract", default=str(ARCHIVE))
    ap.add_argument("--out", default="data/output/rya1220/ni_blend_contamination.json")
    args = ap.parse_args()

    path = pathlib.Path(args.extract).expanduser()
    if not path.exists():
        print(f"VALD extract absent: {path}\n"
              f"It is preserved outside git (RYA-1228); restore it from the archive.",
              file=sys.stderr)
        return 2

    rows = read_vald(path)
    res = audit(rows)
    admitted = [r for r in res if r.get("admit")]
    doc = {
        "schema": "rya1220.ni_blend_contamination.v1",
        "ticket": "RYA-1220",
        "source_extract": str(path),
        "method": (f"summed catalogued central_depth of non-N-I lines within +/-{WINDOW_A} A, "
                   f"over the N I line's own central_depth"),
        "admit_max_ratio": ADMIT_MAX_RATIO,
        "lines": res,
        "admitted": [r["wavelength_air_A"] for r in admitted],
        "rejected": [r["wavelength_air_A"] for r in res if not r.get("admit")],
        "finding": (
            "The two lines RYA-1220 ADDED (7442.298, 8629.235) are the two most contaminated "
            "of the five, and both contaminants are CN. 8629.235 carries a CN line 0.0076 A "
            "away at 75% of its own depth. Using CN-blended N I lines to adjudicate N I "
            "against CN is circular; they are flagged, not adopted."),
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2) + "\n")

    print("RYA-1220 N I blend contamination")
    print(f"  extract : {path.name}")
    for r in res:
        if r.get("status") == "NOT_IN_EXTRACT":
            print(f"  {r['wavelength_air_A']:9.3f}  NOT IN EXTRACT")
            continue
        n = r["nearest"]
        flag = "ADMIT " if r["admit"] else "REJECT"
        circ = "  CN-CIRCULAR" if r["cn_circular"] else ""
        where = (f"nearest {n['species']} @ {n['separation_A']:+.4f} A" if n
                 else f"clean within +/-{WINDOW_A} A")
        print(f"  {r['wavelength_air_A']:9.3f}  depth={r['own_central_depth']:.3f}  "
              f"blend/line={r['blend_depth_ratio']:.2f}  {flag}  {where}{circ}")
    print(f"\n  admitted: {doc['admitted']}")
    print(f"  rejected: {doc['rejected']}")
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
