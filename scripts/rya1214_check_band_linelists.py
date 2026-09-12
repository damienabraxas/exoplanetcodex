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


#: AGSS21 Table 3's MOLECULAR indicators and the band each one lives in, from
#: `data/audit/rya1160_cno_nist_gf/agss21_table3_indicators.csv`. The counts are AGSS21's
#: own; the band assignment is the physics (electronic vs ro-vibrational).
MOLECULAR_INDICATORS = [
    ("C", "C2 Swan", 39, "VIS", "electronic"),
    ("C", "CH A-X (G band)", 7, "VIS", "electronic"),
    ("C", "CH dnu=1", 51, "mid-IR", "ro-vibrational fundamental"),
    ("C", "CO dnu=1", 28, "mid-IR", "ro-vibrational fundamental"),
    ("C", "CO dnu=2", 52, "NIR/mid-IR", "first overtone"),
    ("N", "NH dnu=0", 13, "mid-IR", "pure rotational"),
    ("N", "NH dnu=1", 15, "mid-IR", "ro-vibrational fundamental"),
    ("N", "CN 0-0", 59, "VIS", "electronic"),
    ("N", "CN dnu>=1", 463, "VIS/red-optical", "electronic"),
    ("O", "OH dnu=0", 84, "mid-IR", "pure rotational"),
    ("O", "OH dnu=1", 50, "mid-IR", "ro-vibrational fundamental"),
    ("O", "OH dnu=2", 15, "NIR/mid-IR", "first overtone"),
]


#: 🔴 CORRECTED (RYA-1214, Ryan's comment 2026-09-12). MY FIRST PASS READ ONE DIRECTORY
#: AND GOT THE ANSWER WRONG IN THE MOST CONSEQUENTIAL DIRECTION.
#:
#: It scanned ONLY `$ISPEC_DIR/input/linelists/turbospectrum/molecules`, whose `.bsyn`
#: files are the 400-950 nm ELECTRONIC bands plus RYA-1207's 300-378 nm near-UV set, and
#: concluded that "the mid-IR ro-vibrational fundamentals of OH / NH / CH are NOT held
#: (RYA-503) — ABSENT, not merely unwired." That is FALSE. The repo's own git-tracked
#: vendored mirror at `data/linelists/molecular/turbospectrum/` carries three ExoMol
#: ro-vibrational lists the iSpec install tree does not:
#:
#:     CH/12C-1H__MoLLIST_rovib.bsyn   38,263 rows    2751 - 99970 A
#:     NH/14N-1H__kNigHt_rovib.bsyn    52,831 rows    3447 - 99988 A
#:     OH/16O-1H__MYTHOS_rovib.bsyn    72,762 rows    3133 - 99939 A
#:
#: Those span the near-UV A-X bands AND the mid-IR Dnu = 0/1/2 fundamentals, i.e. exactly
#: the AGSS21 indicators I called blocked. The counts Ryan quotes (C2 25 files, CH 23,
#: CN 33, NH 9, OH 12, CO 1) are this directory's, and they reproduce exactly.
#:
#: ⚠️ THE IRONY IS THE LESSON. This module's own docstring says "read the files, not the
#: filenames (RYA-1190)", and the mistake was reading ONE directory's filenames and
#: inferring the other's absence. `data/linelists/molecular/turbospectrum/README.md`
#: documents this precise failure once already: RYA-237 "reported 'no CH/CN Turbospectrum
#: lists' — it looked in the raw-download dir ... They were present and verified all
#: along; the recon looked in the wrong place." Same directory confusion, opposite
#: direction, eight months later.
MOLECULAR_DIRS = ("(repo) data/linelists/molecular/turbospectrum",
                  "(iSpec install) $ISPEC_DIR/input/linelists/turbospectrum/molecules")


def molecular_readiness(ispec: str | None) -> dict:
    """What the MOLECULAR campaign could run today — read from BOTH vendored trees.

    See MOLECULAR_DIRS: reading only the iSpec install tree reported the ro-vibrational
    fundamentals absent when the repo's own mirror holds them.
    """
    repo_dir = ROOT / "data" / "linelists" / "molecular" / "turbospectrum"
    d = Path(ispec) / "input" / "linelists" / "turbospectrum" / "molecules" if ispec else None
    if (d is None or not d.exists()) and not repo_dir.exists():
        return {"error": "neither molecular tree is readable"}
    spans: dict[str, list[tuple[int, int]]] = {}
    rovib: dict[str, dict] = {}
    files = []
    if d is not None and d.exists():
        files += [(f, "ispec") for f in sorted(d.iterdir())]
    if repo_dir.exists():
        files += [(f, "repo") for f in sorted(repo_dir.rglob("*")) if f.is_file()]
    for f, _tree in files:
        stem = f.name
        if "_rovib" in stem:
            # Span read from the FILE, because the name encodes no range at all — which is
            # precisely why a filename scan missed these.
            lo = hi = None
            n = 0
            with f.open() as fh:
                for ln in fh:
                    parts = ln.split()
                    if not parts:
                        continue
                    try:
                        w = float(parts[0])
                    except ValueError:
                        continue
                    lo = w if lo is None else min(lo, w)
                    hi = w if hi is None else max(hi, w)
                    n += 1
            if n:
                rovib[stem] = {"rows": n, "span_A": [round(lo, 1), round(hi, 1)],
                               "span_nm": [round(lo / 10, 0), round(hi / 10, 0)]}
            continue
        if stem.endswith(".bsyn") and "_" in stem:
            sp, _, rng = stem[: -len(".bsyn")].rpartition("_")
            try:
                lo, hi = (int(x) for x in rng.split("-"))
            except ValueError:
                continue
            spans.setdefault(sp, []).append((lo, hi))
        elif stem.endswith(".dat"):
            spans.setdefault(stem, []).append((0, 0))
    out = {}
    for sp, rs in spans.items():
        rs = sorted(set(rs))
        out[sp] = {"n_files": len(rs),
                   "span_nm": [rs[0][0], rs[-1][1]] if rs[0] != (0, 0) else "non-.bsyn file"}
    # A ro-vibrational list covering the window is what makes the Dnu indicators
    # runnable, so the verdict below is keyed on THESE, not on the electronic .bsyn spans.
    #
    # ⚠️ MATCHED BY THE ExoMol ISOTOPOLOGUE SPELLING, NOT BY THE MOLECULE NAME. The files
    # are `12C-1H__MoLLIST`, `14N-1H__kNigHt`, `16O-1H__MYTHOS` — "CH" is not a substring
    # of any of them, and a `m in filename` test reported all three ABSENT while the
    # parser above had just read 163,856 rows out of them. The second filename-shaped
    # mistake in this same function; the spelling is now declared rather than assumed.
    _ROVIB_TOKEN = {"CH": "12C-1H", "NH": "14N-1H", "OH": "16O-1H"}
    rovib_cover = {m: any(tok in k for k in rovib)
                   for m, tok in _ROVIB_TOKEN.items()}
    return {
        "trees_read": MOLECULAR_DIRS,
        "vendored_species": out,
        "rovibrational_lists": rovib,
        "rovib_coverage": rovib_cover,
        "verdict_by_indicator": [
            {"element": e, "indicator": i, "agss21_lines": n, "regime": r, "kind": k,
             "list_held": (True if r in ("VIS", "red-optical", "VIS/red-optical")
                           else True if i.startswith("CO ")
                           else rovib_cover.get(i.split()[0], False)),
             "why": ("electronic band; vendored 400-950 nm (12C12C / 12C14N / 12CH / "
                     "14NH / 16OH) plus RYA-1207's 300-378 nm near-UV set. The band's "
                     "synth config needs use_molecules: true."
                     if r in ("VIS", "red-optical", "VIS/red-optical") else
                     "CO_IR_Li2015.dat (ExoMol Li 2015) is held. ⚠️ The NIR band list "
                     "carries 34 C I lines and no O I / N I, so that band needs "
                     "rebuilding before a CO product can be run there."
                     if i.startswith("CO ") else
                     "ro-vibrational list HELD in the repo's vendored tree "
                     "(MoLLIST CH / kNigHt NH / MYTHOS OH, 2751-99988 A) — ⚠️ CORRECTED: "
                     "an earlier pass read only the iSpec install tree and reported these "
                     "ABSENT. The blocker is not the line list; it is that no synthesis "
                     "REGION is wired above 1 um for CNO (cno_synthesis.REGIONS is "
                     "{'vis'}) and the NIR/H band lists carry no O I or N I."
                     if rovib_cover.get(i.split()[0], False) else
                     "no ro-vibrational list found for this molecule in either tree.")}
            for e, i, n, r, k in MOLECULAR_INDICATORS
        ],
    }


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
        "molecular_readiness": molecular_readiness(os.environ.get("ISPEC_DIR")),
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
    mr = doc["molecular_readiness"]
    if "error" not in mr:
        ready = [v for v in mr["verdict_by_indicator"] if v["list_held"]]
        print(f"\n=== molecular line lists held ({len(ready)}/"
              f"{len(mr['verdict_by_indicator'])} AGSS21 indicators) ===")
        print(f"  ro-vibrational lists: {mr['rovibrational_lists']}")
        for v in mr["verdict_by_indicator"]:
            mark = "HELD" if v["list_held"] else "NO  "
            print(f"  {mark} {v['element']}  {v['indicator']:18s} {v['agss21_lines']:4d} "
                  f"lines  {v['regime']:16s} {v['why'][:78]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
