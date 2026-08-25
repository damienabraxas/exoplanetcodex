#!/usr/bin/env python3
"""RYA-1042: extract the solar ⟨3D⟩ Fe anchor slice from Amarsi+2016.

WHAT THIS ANCHOR IS, AND WHY IT IS THE RIGHT ONE
------------------------------------------------
`iron_abundancecorr.tar.gz` (Amarsi, Lind, Asplund, Barklem & Collet 2016, MNRAS 463,
1518; DOI 10.1093/mnras/stw2077 — verified at Crossref on acquisition, not cited from
memory) ships THREE grids, and its own readme says which is which:

    nmarcs_lmarcs.txt   1D non-LTE  −  1D  LTE     on 1D MARCS
    nmtd_lmarcs.txt     ⟨3D⟩ non-LTE −  1D  LTE     mixed atmospheres
    nmtd_lmtd.txt       ⟨3D⟩ non-LTE − ⟨3D⟩ LTE     on averaged-3D STAGGER   <-- THIS ONE

🔴 **`nmtd_lmtd` IS THE ONLY SCALE-MATCHED ONE, AND THE READ ME SAYS SO ITSELF**: *"Only
use this if you are correcting LTE results based on ⟨3D⟩ or 3D model atmospheres."* Our
RYA-1040 product is exactly that difference — ⟨3D⟩-NLTE minus ⟨3D⟩-LTE on the STAGGER
averaged atmosphere — so the two quantities are the same physical object computed by two
independent groups.

⚠️ **`nmtd_lmarcs` WOULD HAVE BEEN THE TRAP.** It is also a ⟨3D⟩ non-LTE number and it is
also Amarsi+2016, but its comparand is 1D LTE — so it carries the 1D→mean-3D ATMOSPHERE
shift inside it. Gating our differential against that one would fold an atmosphere
difference into an NLTE test, which is the RYA-542 confound the whole paired-product
design exists to avoid.

⚠️ **AND IT IS NOT Amarsi+2022.** RYA-817's MLP is FULL 3D. Gating a mean-⟨3D⟩ deck
against a full-3D reference measures mean-vs-full, not deck-vs-reference.

WHAT THE AUDIT FOUND (the file is self-documenting; nothing here is inferred)
----------------------------------------------------------------------------
Header, verbatim:

    T_eff/K, log10(g/cm s^-2), log10(eps_Fe), v_turb/km s^-1, Species,
    lambda_{Air,centre}/nm, E_low, log10(gf), lambda_{Air,min}/nm,
    lambda_{Air,max}/nm, Clean, Abundance correction

  * WAVELENGTH MEDIUM: **AIR**, stated in the column name. Our line lists are air too, so
    no conversion — but it is checked rather than assumed, because a medium mix-up is a
    ~1.4 Å error at 5000 Å and would silently mismatch every line.
  * UNITS: wavelength **nm** (our band arguments are Å); correction in **dex**.
  * ABUNDANCE AXIS is ABSOLUTE `log10(eps_Fe)`, not [Fe/H] — 2.50 to 8.00 in 0.25 steps.
    7.50 is a node, which is the Gerber deck's own A(Fe) (RYA-1035), so the anchor is read
    at the same abundance rather than interpolated to it.
  * GRID: Teff 4000–7000 / 250 K; logg 1.50–5.00 / 0.5; vturb 0.75, 1.50, 3.00.
  * ⚠️ **THE TITLE SAYS "metal-poor stars" AND THE GRID STILL COVERS SOLAR.** Paper III is
    a metal-poor analysis; the released grid is not restricted to it. Checked rather than
    assumed from the title — the solar node exists at (5750, 4.50, 7.50).
  * ⚠️ **THERE IS NO vturb = 1.0 NODE.** Solar vturb is 1.0 and the grid brackets it at
    0.75 and 1.50. The anchor is therefore a BRACKET, not a point, and this script emits
    both — an interpolated single number would hide that the grid was never asked.
  * `Clean` marks lines the authors consider unblended. Everything downstream filters on
    it; the counts are kept so the filter is visible rather than implied.

⚠️ **A SENTINEL LIVES IN THIS FILE.** At vturb 3.00 the Fe I minimum is exactly −4.0000,
which is a floor value and not a correction. Rows at exactly −4.0 are flagged, never
averaged in.

🔴 **THE ATMOSPHERE LEG COMES OUT OF THE SAME ARCHIVE, BY SUBTRACTION** (RYA-1042 scope
add). `nmtd_lmarcs` and `nmtd_lmtd` share the `<3D>`-non-LTE term exactly, so

    nmtd_lmarcs - nmtd_lmtd = (<3D>NLTE - 1D LTE) - (<3D>NLTE - <3D>LTE)
                            =  <3D>LTE - 1D LTE          <-- the ATMOSPHERE effect

That is Amarsi's own 1D->mean-3D shift for the same lines at the same node, computed by an
independent group, and it is what gates our `<3D>`-LTE leg. It is a SUBTRACTION of two
released columns, not a model of ours -- nothing is fitted and nothing is assumed.

⚠️ **THE JOIN MUST BE EXACT, NOT NEAREST.** Both files are the same grid over the same line
list, so a solar-node slice of each joins 540/540 on wavelength. A nearest-match join would
paper over a line-list difference if one ever appeared; this refuses instead. And a
sentinel in EITHER column poisons the difference, so a row is dropped if either end is -4.0.

    python3 scripts/rya1042_extract_amarsi2016_anchor.py \
        --source <nmtd_lmtd.txt> --source-lmarcs <nmtd_lmarcs.txt>
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "nlte_grids" / "amarsi2016_fe"

#: The solar node, in the grid's own coordinates. Teff/logg are the NEAREST NODES to the
#: pipeline's solar star (5772 / 4.438); they are named explicitly rather than searched, so
#: a change of node is a visible edit.
SOLAR_TEFF, SOLAR_LOGG, SOLAR_EPS = "5750", "4.50", "7.50"

#: A floor value in the released file, not a physical correction.
SENTINEL = -4.0

COLUMNS = ("teff_K", "logg", "log_eps_Fe", "vturb_kms", "species", "lambda_air_nm",
           "e_low", "log_gf", "lambda_air_min_nm", "lambda_air_max_nm", "clean",
           "delta_mean3d_nlte_minus_mean3d_lte",
           #: <3D>LTE - 1D LTE, by subtraction of the two released grids (see module
           #: docstring). Empty when the row is absent from `nmtd_lmarcs` or when either
           #: end is the -4.0 sentinel -- never 0.0, which would read as "no shift".
           "atmosphere_mean3d_lte_minus_1d_lte")


#: The join key for pairing the two grids. Wavelength alone would be ambiguous if two
#: transitions ever printed the same rounded lambda, so species and vturb ride along --
#: the same fields that make a row unique within one node.
def _key(r: dict) -> tuple:
    return (r["species"], r["vturb_kms"], r["lambda_air_nm"])


def _solar_rows(source: Path, n_fields: int = 12) -> list[list[str]]:
    rows = []
    with source.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) != n_fields:
                continue
            if f[0] == SOLAR_TEFF and f[1] == SOLAR_LOGG and f[2] == SOLAR_EPS:
                rows.append(f)
    return rows


def extract(source: Path, source_lmarcs: Path | None = None) -> list[dict]:
    base = _solar_rows(source)
    if not base:
        raise SystemExit(
            f"no rows at the solar node ({SOLAR_TEFF}, {SOLAR_LOGG}, {SOLAR_EPS}) in "
            f"{source} -- refusing to emit an empty anchor")
    rows = [dict(zip(COLUMNS, f)) for f in base]

    if source_lmarcs is None:
        for r in rows:
            r["atmosphere_mean3d_lte_minus_1d_lte"] = ""
        return rows

    # 🔴 The atmosphere leg. EXACT join -- see the module docstring.
    lm = {}
    for f in _solar_rows(source_lmarcs):
        d = dict(zip(COLUMNS[:12], f))
        lm[_key(d)] = float(d["delta_mean3d_nlte_minus_mean3d_lte"])
    if not lm:
        raise SystemExit(
            f"no rows at the solar node in {source_lmarcs} -- refusing to emit an "
            f"atmosphere column from an empty join")

    matched = 0
    for r in rows:
        a = lm.get(_key(r))
        b = float(r["delta_mean3d_nlte_minus_mean3d_lte"])
        # A sentinel at EITHER end poisons the difference. Empty, never 0.0 -- a zero
        # here would read as "the atmosphere does nothing", which is a claim.
        if a is None or a == SENTINEL or b == SENTINEL:
            r["atmosphere_mean3d_lte_minus_1d_lte"] = ""
            continue
        r["atmosphere_mean3d_lte_minus_1d_lte"] = f"{a - b:.4f}"
        matched += 1

    print(f"  atmosphere join: {matched} of {len(rows)} rows carry <3D>LTE - 1D LTE "
          f"({len(lm)} solar-node rows in {source_lmarcs.name})")
    return rows


def summarise(rows: list[dict]) -> dict:
    import statistics
    out = {}
    for species in ("Fe1", "Fe2"):
        for vturb in ("0.75", "1.50", "3.00"):
            sel = [float(r["delta_mean3d_nlte_minus_mean3d_lte"]) for r in rows
                   if r["species"] == species and r["vturb_kms"] == vturb
                   and r["clean"] == "yes"]
            sentinels = [d for d in sel if d == SENTINEL]
            sel = [d for d in sel if d != SENTINEL]
            if not sel:
                continue
            atm = [float(r["atmosphere_mean3d_lte_minus_1d_lte"]) for r in rows
                   if r["species"] == species and r["vturb_kms"] == vturb
                   and r["clean"] == "yes"
                   and r.get("atmosphere_mean3d_lte_minus_1d_lte") not in (None, "")]
            out[f"{species}_vturb{vturb}"] = {
                "n_clean": len(sel), "n_sentinel_excluded": len(sentinels),
                "median": round(statistics.median(sel), 4),
                "mean": round(statistics.fmean(sel), 4),
                "min": round(min(sel), 4), "max": round(max(sel), 4),
                "atmosphere_n": len(atm),
                "atmosphere_median": round(statistics.median(atm), 4) if atm else None,
                "atmosphere_mean": round(statistics.fmean(atm), 4) if atm else None}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True,
                    help="path to nmtd_lmtd.txt (Sirius-only; the release is 748 MB)")
    ap.add_argument("--source-lmarcs", default=None,
                    help="path to nmtd_lmarcs.txt -- enables the ATMOSPHERE column "
                         "(<3D>LTE - 1D LTE) by subtraction. Omit and that column is "
                         "emitted EMPTY rather than absent, so a consumer can tell "
                         "'not extracted' from 'no shift'.")
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    a = ap.parse_args()

    rows = extract(Path(a.source),
                   Path(a.source_lmarcs) if a.source_lmarcs else None)
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "amarsi2016_mean3d_solar_anchor.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    doc = {
        "ticket": "RYA-1042",
        "extracted_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "quantity": "<3D> non-LTE minus <3D> LTE abundance correction, on averaged-3D "
                    "STAGGER model atmospheres (the file's own readme)",
        "node": {"teff_K": SOLAR_TEFF, "logg": SOLAR_LOGG, "log_eps_Fe": SOLAR_EPS},
        "n_rows": len(rows),
        "summary_clean_only": summarise(rows),
    }
    (out_dir / "amarsi2016_mean3d_solar_anchor.json").write_text(
        json.dumps(doc, indent=2) + "\n")

    print(f"wrote {csv_path} ({len(rows)} rows at the solar node)")
    for k, v in doc["summary_clean_only"].items():
        print(f"  {k:<16} n={v['n_clean']:>4} median={v['median']:+.4f} "
              f"[{v['min']:+.4f}, {v['max']:+.4f}]"
              + (f"  atm={v['atmosphere_median']:+.4f} (n={v['atmosphere_n']})"
                 if v.get("atmosphere_median") is not None else "")
              + (f"  ⚠️ {v['n_sentinel_excluded']} sentinel rows excluded"
                 if v["n_sentinel_excluded"] else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
