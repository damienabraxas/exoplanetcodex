#!/usr/bin/env python3
"""Which published ⟨3D⟩ aux tables carry the RYA-1035 metallicity-zeroing defect?

RYA-1035 found it in Fe's aux and had to decide, for Fe alone, which of the two shipped
files to register against. That decision is needed **once per element**, and discovering it
one element at a time is exactly the pattern this project keeps paying for — so this sweeps
all of them at once, from the source.

## The defect

The vendor's aux tables key each departure record to an atmosphere by (Teff, logg, [Fe/H]).
On some decks the `[Fe/H]` COLUMN reads +0.00 on the seven Teff=5777 rows — the STAGGER
solar member, whose model NAMES (`…m00 / m05 / m10 / m20 / m30 / m40 / p05`) span the full
metallicity axis. The names are right and the column is wrong: on every clean deck, and on
every non-solar row of a defective one, the two agree exactly.

Left unrefereed it is **not a crash but a wrong star**: all seven tie at the solar node and
file order wins, which is the `[Fe/H] = −1.0` model. `gerber_nlte._parse_aux_text` referees
the column with the name and marks the row SUSPECT; `read_deck_node` then refuses it.

## And the `_marcs_names` sibling can be UNRECOVERABLE

`convert_3d_grid_to_marcs_names.py` builds the new name **from the `[Fe/H]` column**. On a
defective deck it therefore propagates the zeroing into the NAME, collapsing seven distinct
atmospheres into one byte-identical string — after which nothing can referee it, because
name and column agree and both are wrong. **A defective deck must be registered against the
PLAIN aux; a clean one may use either.** Al uses `_marcs_names`; Fe must not.

    python3 scripts/rya1035_mean3d_aux_defect_sweep.py            # live, all elements
    python3 scripts/rya1035_mean3d_aux_defect_sweep.py --offline  # re-render from the artifact
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT_DIR = ROOT / "data" / "results" / "rya1035"
KEEPER = "6eaecbf95b88448f98a4"
API = f"https://keeper.mpdl.mpg.de/api/v2.1/share-links/{KEEPER}/dirents/"
DL = f"https://keeper.mpdl.mpg.de/d/{KEEPER}/files/?p="

#: The pipeline's solar star (IAU) against the deck's STAGGER solar member (5777/4.44).
SOLAR = dict(teff=5772.0, logg=4.438, feh=0.0)


def _get(url: str, timeout: int = 120) -> bytes:
    with urllib.request.urlopen(
            urllib.request.Request(url, headers={"Accept": "*/*"}), timeout=timeout) as fh:
        return fh.read()


def _fetch_text(element: str, name: str) -> str:
    """⚠️ MUST follow redirects: the share 302s to its file server, and a bare request
    returns 0 bytes. An empty aux parses to zero rows and reads as 'no defect found' --
    an absence manufactured by a fetch bug."""
    data = _get(DL + urllib.parse.quote(f"/dep-grids/{element}/{name}") + "&dl=1")
    text = data.decode("utf-8", "replace")
    if not text.strip():
        raise RuntimeError(f"{element}/{name}: empty response -- refusing to score it")
    return text


def list_aux(element: str) -> dict:
    """The ⟨3D⟩ aux pair for one element: the plain table and its converted sibling."""
    listing = json.loads(_get(API + "?path=" + urllib.parse.quote(f"/dep-grids/{element}/")))
    names = [e.get("file_name") for e in listing["dirent_list"] if not e.get("is_dir")]
    aux = [n for n in names if n and n.startswith("auxData_")
           and ("STAGGER" in n or "mean3D" in n or "mean3d" in n)]
    plain = [n for n in aux if "_marcs_names" not in n]
    conv = [n for n in aux if "_marcs_names" in n]
    return {"plain": plain[0] if plain else None, "marcs_names": conv[0] if conv else None}


def score(text: str) -> dict:
    """Referee one aux table with `gerber_nlte`'s own parser -- the same code the pipeline
    runs, so this measures the shipped behaviour rather than a re-implementation of it."""
    from pipeline import gerber_nlte as G
    rows = G._parse_aux_text(text)
    if not rows:
        raise RuntimeError("no parseable rows -- refusing to score")
    over = [r for r in rows if r["feh_from_name"]]
    solar = [r for r in rows
             if abs(r["teff"] - SOLAR["teff"]) <= G._NODE_TOL_TEFF
             and abs(r["logg"] - SOLAR["logg"]) <= G._NODE_TOL_LOGG
             and abs(r["feh"] - SOLAR["feh"]) <= G._NODE_TOL_FEH]
    at_5777 = [r["id"] for r in rows if r["teff"] == 5777.0]
    return {
        "n_rows": len(rows),
        "n_feh_overridden": len(over),
        "overridden_teff_logg": sorted({(r["teff"], r["logg"]) for r in over}),
        "n_rows_at_5777": len(at_5777),
        "n_distinct_names_at_5777": len(set(at_5777)),
        "solar_node_distinct_ids": len({r["id"] for r in solar}),
        "solar_node_distinct_abundances": len({round(r["abundance"], 3) for r in solar}),
        "unparseable_names": sum(1 for r in rows
                                 if G._node_from_model_name(r["id"]) is None),
    }


def verdict(plain: dict, conv: dict | None) -> tuple[str, str]:
    """(aux to register against, why). Derived from the two scores, never chosen."""
    if plain["n_feh_overridden"] == 0:
        return "either", ("clean: the [Fe/H] column agrees with every model name, so the "
                          "converted sibling is faithful and either file works")
    # defective -- is the conversion recoverable?
    if conv is None:
        return "plain", "defective, and no converted sibling is published"
    if conv["n_distinct_names_at_5777"] == 1 and conv["n_rows_at_5777"] > 1:
        return "plain", (
            f"DEFECTIVE and the conversion is UNRECOVERABLE: {plain['n_feh_overridden']} "
            f"rows carry a zeroed [Fe/H], and the converter propagated it into the NAME -- "
            f"{conv['n_rows_at_5777']} rows at Teff=5777 collapse to ONE name, after which "
            f"nothing can referee it. Register against the plain aux (Fe/Mn pattern).")
    return "plain", "defective; the plain aux keeps the names that referee the column"


ELEMENTS = ("Al", "Ba", "Ca", "Co", "Cr", "Eu", "Fe", "H", "Mg", "Mn", "Na", "Ni", "O",
            "Si", "Sr", "Ti", "Y")


def build() -> dict:
    out = {"ticket": "RYA-1035",
           "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "share": KEEPER, "elements": {}}
    for el in ELEMENTS:
        aux = list_aux(el)
        if not aux["plain"]:
            out["elements"][el] = {"error": "no plain <3D> aux published"}
            continue
        p = score(_fetch_text(el, aux["plain"]))
        c = score(_fetch_text(el, aux["marcs_names"])) if aux["marcs_names"] else None
        use, why = verdict(p, c)
        out["elements"][el] = {"plain_file": aux["plain"], "marcs_names_file": aux["marcs_names"],
                               "plain": p, "marcs_names": c, "register_against": use,
                               "why": why}
        print(f"  {el:<3} rows={p['n_rows']:>5} overridden={p['n_feh_overridden']:>4}  "
              f"-> register against: {use}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true",
                    help="re-render the CSV from the committed JSON instead of fetching")
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    a = ap.parse_args()

    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    jf = out_dir / "mean3d_aux_defect_sweep.json"

    doc = json.loads(jf.read_text()) if a.offline else build()
    if not a.offline:
        jf.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")

    with open(out_dir / "mean3d_aux_defect_sweep.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "plain_file", "n_rows", "n_feh_overridden",
                    "marcs_names_collapsed", "register_against", "why"])
        for el, d in doc["elements"].items():
            if "error" in d:
                w.writerow([el, "", "", "", "", "", d["error"]])
                continue
            c = d.get("marcs_names")
            collapsed = bool(c and c["n_distinct_names_at_5777"] == 1
                             and c["n_rows_at_5777"] > 1)
            w.writerow([el, d["plain_file"], d["plain"]["n_rows"],
                        d["plain"]["n_feh_overridden"], collapsed,
                        d["register_against"], d["why"]])

    bad = [el for el, d in doc["elements"].items()
           if "error" not in d and d["plain"]["n_feh_overridden"]]
    print(f"\nDEFECTIVE: {bad or 'none'}  ({len(bad)} of {len(doc['elements'])})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
