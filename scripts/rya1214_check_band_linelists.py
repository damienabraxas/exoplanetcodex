#!/usr/bin/env python3
"""RYA-1214 Step 3 — does each band's SYNTHESIS line list already carry the CNO indicators?

The ticket's Step 3 says: "If a band-specific line list already exists and covers CNO,
confirm it; if not, build it." This is the confirmation, and it is done by looking in the
lists rather than by reasoning about what they ought to contain — `feedback: read the
files, not the filenames` (RYA-1190 read nm as Angstrom and scoped a fetch for data it
already had).

WHAT IT CHECKS. For every band in `config/synth_bands.yaml`, it opens the list that band
actually resolves to (including `ispec_ges_v6`, which is not a repo path but the
iSpec-vendored `GESv6_atom_hfs_iso.420_920nm` under `$ISPEC_DIR`) and asks two things:

  1. how many C I / C II / N I / N II / O I / O II lines it holds inside the band, and
  2. which of the AGSS21 Table 3 ATOMIC indicators are present, by name.

The second is the one that decides whether a product can be run, because a band list with
thousands of CNO lines and no [O I] 6300 cannot measure the oxygen the campaign is about.

⚠️ MUST RUN WHERE `$ISPEC_DIR` IS. The GES v6 list is not vendored in this repo — it
ships with iSpec, which lives on Sirius. Run there; the JSON it writes is the artifact.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.synth_bands import SYNTH_BANDS  # noqa: E402

OUT = ROOT / "data" / "audit" / "rya1214_cno_products"
CENSUS = ROOT / "data" / "audit" / "rya1136_cno_intake" / "atomic_source_census.csv"

#: The line list's own species spelling.
SPECIES = {"C I": "C 1", "C II": "C 2", "N I": "N 1", "N II": "N 2",
           "O I": "O 1", "O II": "O 2"}
#: Derived from the source's printed precision: the census quotes 0.01 nm = 0.1 A.
#: Same window Step 2 uses, for the same reason (RYA-1109).
TOL_A = 0.1


def list_path(cfg) -> Path:
    """Resolve a band's list, GES v6 included. `SynthBand.linelist` already does this,
    but it imports `pipeline.abundances_derive`, which imports iSpec — so the failure
    mode off-Sirius is an ImportError rather than a clear message about $ISPEC_DIR."""
    try:
        return Path(cfg.linelist)
    except Exception:
        ispec = os.environ.get("ISPEC_DIR")
        if not ispec:
            raise SystemExit(
                "this band resolves to the iSpec-vendored GES v6 list and $ISPEC_DIR is "
                "unset. That list is NOT vendored in this repo — run this on Sirius with "
                "ISPEC_DIR=/mnt/codex-data/engines/ispec_src.")
        return (Path(ispec) / "input" / "linelists" / "transitions"
                / "GESv6_atom_hfs_iso.420_920nm" / "atomic_lines.tsv")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cen = pd.read_csv(CENSUS)
    cen = cen[cen.element.isin(["C", "N", "O"])]

    bands, indicators = [], []
    for name, cfg in SYNTH_BANDS.items():
        p = list_path(cfg)
        if not p.exists():
            bands.append({"band": name, "linelist": str(p), "exists": False,
                          "verdict": "LIST ABSENT — regenerate before running this band",
                          "build_hint": cfg.build_hint})
            continue
        t = pd.read_csv(p, sep="\t", low_memory=False,
                        usecols=["element", "wave_A", "loggf", "lower_state_eV",
                                 "theoretical_depth"])
        el = t.element.astype(str).str.strip()
        inband = t[(t.wave_A >= cfg.lo_A) & (t.wave_A <= cfg.hi_A)]
        el_in = inband.element.astype(str).str.strip()
        counts = {sp: int((el_in == tok).sum()) for sp, tok in SPECIES.items()}
        bands.append({
            "band": name, "linelist": p.name, "exists": True,
            "lo_A": cfg.lo_A, "hi_A": cfg.hi_A,
            "list_span_A": [round(float(t.wave_A.min()), 2), round(float(t.wave_A.max()), 2)],
            "total_lines_in_band": int(len(inband)),
            "cno_lines_in_band": {k: v for k, v in counts.items() if v},
            "cno_total_in_band": int(sum(counts.values())),
            "use_molecules": bool(getattr(cfg, "use_molecules", False)),
        })
        # The indicators, by name.
        for _, r in cen.iterrows():
            w = float(r.wavelength_air_A)
            if not (cfg.lo_A <= w <= cfg.hi_A):
                continue
            tok = SPECIES.get(str(r.species))
            if tok is None:
                continue
            hit = t[(el == tok) & ((t.wave_A - w).abs() <= TOL_A)]
            indicators.append({
                "band": name, "element": r.element, "species": r.species,
                "line_label": r.line_label, "census_wavelength_A": w,
                "n_in_list": int(len(hit)),
                "list_wavelength_A": round(float(hit.wave_A.iloc[hit.loggf.values.argmax()]), 3)
                                     if len(hit) else None,
                "list_log_gf": round(float(hit.loggf.max()), 4) if len(hit) else None,
                "list_theoretical_depth": round(float(hit.theoretical_depth.max()), 3)
                                          if len(hit) else None,
                "verdict": "PRESENT" if len(hit) else "ABSENT FROM THE BAND LIST",
            })

    ind = pd.DataFrame(indicators)
    ind.to_csv(OUT / "band_linelist_indicator_coverage.csv", index=False)
    absent = ind[ind.verdict != "PRESENT"] if len(ind) else ind

    doc = {
        "ticket": "RYA-1214", "step": "3 — band synthesis line lists vs the CNO indicators",
        "match_tolerance_A": TOL_A,
        "match_tolerance_basis": "the census prints 0.01 nm = 0.1 A (RYA-1109)",
        "bands": bands,
        "indicators_checked": int(len(ind)),
        "indicators_absent": int(len(absent)),
        "verdict": ("NO NEW LINE LIST IS NEEDED for the bands checked"
                    if len(absent) == 0 else
                    f"{len(absent)} indicator(s) absent from their band's list"),
    }
    (OUT / "band_linelist_check.json").write_text(json.dumps(doc, indent=2) + "\n")

    print("=== RYA-1214 Step 3 — band synthesis line lists ===")
    for b in bands:
        if not b["exists"]:
            print(f"  {b['band']:12s} {b['verdict']}")
            continue
        print(f"  {b['band']:12s} {b['linelist']:34s} {b['total_lines_in_band']:7d} lines "
              f"in band, {b['cno_total_in_band']:5d} CNO   molecules={b['use_molecules']}")
        print(f"               {b['cno_lines_in_band']}")
    if len(ind):
        print(f"\n  AGSS21 atomic indicators found in their band's list: "
              f"{int((ind.verdict == 'PRESENT').sum())}/{len(ind)}")
        if len(absent):
            print(absent[["band", "element", "line_label",
                          "census_wavelength_A"]].to_string(index=False))
    print(f"\n  {doc['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
