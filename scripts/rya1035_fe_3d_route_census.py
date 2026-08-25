#!/usr/bin/env python3
"""RYA-1035 Step 0 — what Fe 3D-NLTE already exists as PUBLIC-FETCHABLE.

The ticket's Step 0: *"do not build what we can consume"*. Verify against a LIVE fetch,
on the RYA-1015 availability-matrix pattern, before any build is scoped.

🔴 IT CAME BACK **HAVE**, FOUR TIMES OVER — and the register said the opposite.
`CODEX_STATE_REGISTER.md` v117 states *"Fe <3D> is unaffected — Fe has NO deck
(T2_BUILD_OWED, build-our-own only)"*. That was read off OUR DISK and never checked against
the SOURCE: `NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin` has been sitting on the same MPG
Keeper share we already fetch every Gerber deck from since 2021, and TSFitPy's own canonical
downloader lists it as `[Fe] 3d_bin_link`. We fetched Al, Cr, Eu and Y from that share and
never looked in the Fe folder.

WHAT THIS SCRIPT CHECKS
  1. the Keeper share, live, for every element's <3D> deck (is Fe really there?);
  2. Amarsi's homepage data index, live, for the Fe products it publishes;
  3. the vendored Amarsi+2022 Fe MLP — does it actually run on main, or is the matrix's
     `BROKEN` cell stale? (this leg needs no network);
  4. the two Fe <3D> aux tables for the metallicity defect, against the committed copies.

Legs 1-2 need the network; `--offline` skips them and runs 3-4 only.

    python3 scripts/rya1035_fe_3d_route_census.py
    python3 scripts/rya1035_fe_3d_route_census.py --offline
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
AUX_DIR = ROOT / "data" / "nlte_grids" / "gerber_ts" / "fe_mean3d_aux"

KEEPER_TOKEN = "6eaecbf95b88448f98a4"
KEEPER_API = f"https://keeper.mpdl.mpg.de/api/v2.1/share-links/{KEEPER_TOKEN}/dirents/"
AMARSI_HOME = "https://www.astro.uu.se/~amarsi/"

#: The solar star as the pipeline knows it (IAU), against the deck's STAGGER solar member.
SOLAR = dict(teff=5772.0, logg=4.438, feh=0.0)


def _get(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return fh.read()


def _head_ok(url: str, timeout: int = 60) -> tuple[bool, int | None]:
    """(reachable, content-length). Follows redirects; Keeper 302s to its file server."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as fh:
            cl = fh.headers.get("Content-Length")
            return True, int(cl) if cl else None
    except Exception:
        return False, None


# ── leg 1: the Keeper share, live ────────────────────────────────────────────

def probe_keeper() -> dict:
    """Which elements ship a STAGGERmean3D deck, and is Fe among them?"""
    out: dict[str, dict] = {}
    listing = json.loads(_get(KEEPER_API + "?path=/dep-grids/"))["dirent_list"]
    elements = sorted(e["folder_name"] for e in listing if e.get("is_dir"))
    for el in elements:
        try:
            files = json.loads(
                _get(KEEPER_API + "?path=" + urllib.parse.quote(f"/dep-grids/{el}/"))
            )["dirent_list"]
        except Exception as exc:                                  # pragma: no cover
            out[el] = {"error": str(exc)}
            continue
        names = [f.get("file_name") for f in files if not f.get("is_dir")]
        mean3d = [n for n in names if n and "STAGGERmean3D" in n and ".bin" in n]
        out[el] = {
            "mean3d_bin": mean3d[0] if mean3d else None,
            "mean3d_bytes": next((f["size"] for f in files
                                  if f.get("file_name") in mean3d), None),
            "has_marcs_names_aux": any(n and "_marcs_names" in n for n in names),
            "n_files": len(names),
        }
    return out


# ── leg 2: Amarsi's homepage index ───────────────────────────────────────────

def probe_amarsi() -> dict:
    """💡 THE HOMEPAGE IS THE INDEX, not the paper and not CDS (RYA-1013 finding, reused).

    It lists three separate Fe products, two of which are 3D and one of which is <3D>.
    """
    import re
    html = _get(AMARSI_HOME).decode("utf-8", "replace")
    links = re.findall(r'href="([^"]+)"[^>]*>(.*?)</a>', html, flags=re.S | re.I)
    fe = []
    for href, label in links:
        text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", label)).strip()
        if re.search(r"\bFe\s?[12]\b", text):
            url = href if href.startswith("http") else urllib.parse.urljoin(AMARSI_HOME, href)
            ok, size = _head_ok(url)
            fe.append({"label": text, "url": url, "reachable": ok, "bytes": size})
    return {"page": AMARSI_HOME, "fe_entries": fe}


# ── leg 3: is the vendored Amarsi+2022 Fe MLP actually broken? ────────────────

def probe_mlp() -> dict:
    """🔴 RE-CHECK THE `BROKEN` CELL. `model_availability_matrix` carries Fe FULL_3D as
    BROKEN (RYA-923): *"the MLP returns NaN for EVERY in-domain line"*. Measured here on
    the committed RYA-817 in-domain pool, on main, pinned and unpinned."""
    import warnings

    import numpy as np
    import pandas as pd
    warnings.filterwarnings("ignore")
    from pipeline import amarsi3d

    per_line = ROOT / "data" / "results" / "rya817" / "rya817_3dnlte_per_line.csv"
    df = pd.read_csv(per_line)
    kw = dict(teff=SOLAR["teff"], logg=SOLAR["logg"], vmic=1.0)
    out = {}
    for (band, ion), grp in df.groupby(["band", "ion"]):
        dom = grp[grp.in_domain.astype(str).str.lower().isin(["true", "1"])]
        dom = dom.dropna(subset=["elo_eV", "eup_eV", "loggf", "a_1dlte"])
        if not len(dom):
            out[f"{band} Fe {ion}"] = {"in_domain": 0}
            continue
        leg = {"in_domain": int(len(dom))}
        for label, pin in (("unpinned_per_line", None), ("pinned_7.46", 7.46)):
            vals = []
            for r in dom.itertuples():
                try:
                    ab, _ = amarsi3d.aberr_for_line(
                        str(r.ion), float(r.elo_eV), float(r.eup_eV), float(r.loggf),
                        float(r.a_1dlte), afe3n_axis=pin, **kw)
                except Exception:
                    ab = np.nan
                vals.append(ab)
            v = np.array(vals, dtype=float)
            leg[label] = {
                "n_finite": int(np.isfinite(v).sum()),
                "median_aberr": (round(float(np.nanmedian(v)), 4)
                                 if np.isfinite(v).any() else None),
            }
        out[f"{band} Fe {ion}"] = leg
    return out


# ── leg 4: the Fe <3D> aux metallicity defect ────────────────────────────────

def probe_aux() -> dict:
    """The wiring blocker Step 0 turned up: which of the two shipped aux tables can
    address the solar node at all."""
    from pipeline import gerber_nlte as G

    out = {}
    for tag, path in (("plain", AUX_DIR / "auxData_Fe_STAGGERmean3D_May-21-2021.txt"),
                      ("marcs_names",
                       AUX_DIR / "auxData_Fe_STAGGERmean3D_May-21-2021_marcs_names.txt")):
        rows = G._parse_aux_text(path.read_text())
        over = [r for r in rows if r["feh_from_name"]]
        solar = [r for r in rows
                 if abs(r["teff"] - SOLAR["teff"]) <= G._NODE_TOL_TEFF
                 and abs(r["logg"] - SOLAR["logg"]) <= G._NODE_TOL_LOGG
                 and abs(r["feh"] - SOLAR["feh"]) <= G._NODE_TOL_FEH]
        out[tag] = {
            "file": path.name,
            "n_rows": len(rows),
            "n_feh_overridden": len(over),
            "overridden_ids": sorted(r["id"] for r in over),
            "solar_node_candidates": sorted({r["id"] for r in solar}),
            "solar_node_addressable": len(solar) == 1,
        }
    return out


# ── the census ───────────────────────────────────────────────────────────────

#: One row per Fe 3D/<3D> route Step 0 found. `verdict` is what the route IS today, never
#: what it might become. RYA-161 firewall: routes are described on physics and access, and
#: never ranked by how close their answer sits to a reference value.
ROUTES = [
    dict(
        route="gerber-mean3d-deck",
        scale="<3D>-NLTE",
        product="NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin (+ aux)",
        source="MPG Keeper share 6eaecbf95b88448f98a4, /dep-grids/Fe/",
        citation="Gerber, Bergemann et al. 2023, A&A 669, A43",
        access="PUBLIC_FETCHABLE",
        coverage="Teff 4000-7000, logg 1.5-5.0, [Fe/H] -4..+0.5; 189 nodes; dedicated "
                 "solar member 5777/4.44; departure deck, so NOT line-list limited",
        verdict="HAVE_NOT_STAGED",
        blocker="72 MB zip not yet on Sirius; vendor aux [Fe/H] defect (see aux_check)",
        note="The Al route transfers UNCHANGED: same layout, same direct reader, same "
             "<3D> atmosphere. Simpler than Al -- one A(X) per node, no abundance axis.",
    ),
    dict(
        route="amarsi2022-mlp",
        scale="3D-NLTE (full)",
        product="1L-3NErrors MLP (fe1_model_gt02.p / fe1_model_lt02.p / fe2_model.p)",
        source="https://github.com/sliljegren/1L-3NErrors (MIT)",
        citation="Amarsi, Liljegren & Nissen 2022, A&A 668, A68",
        access="VENDORED",
        coverage="Teff 5000-6500, logg 4.0-4.5, [Fe/H] -3..0; OPTICAL ONLY "
                 "4787.83-6810.26 A, Jofre+2014 golden lines; A(Fe) axis ceiling 7.5",
        verdict="HAVE_RUNS",
        blocker="none for the pinned-axis path; the A(Fe)=7.5 ceiling rails our VIS "
                "1D-LTE zero point (7.586) -- recorded, sensitivity 0.0066 dex",
        note="The matrix's BROKEN cell is STALE -- see mlp_check.",
    ),
    dict(
        route="amarsi2016-corrections",
        scale="1D & <3D>-NLTE",
        product="iron_abundancecorr.tar.gz",
        source="https://www.astro.uu.se/~amarsi/iron_abundancecorr.tar.gz",
        citation="Amarsi et al. 2016 (1D/mean-3D non-LTE Fe I and Fe II corrections)",
        access="PUBLIC_FETCHABLE",
        coverage="not yet characterised -- 182 MB, unfetched",
        verdict="AVAILABLE_UNFETCHED",
        blocker="not fetched; an INDEPENDENT <3D> cross-check on the Gerber deck",
        note="A second <3D> Fe product from a different group. Its value is as a REFEREE "
             "for the Gerber route, not as a substitute.",
    ),
    dict(
        route="amarsi2019-cofe",
        scale="3D-NLTE",
        product="cofe_tools.tar.gz (C1, O1, Fe2)",
        source="https://www.astro.uu.se/~amarsi/cofe_tools.tar.gz",
        citation="Amarsi, Nissen & Skuladottir 2019, A&A 630, A104",
        access="PUBLIC_FETCHABLE",
        coverage="Fe II leg; 3D below Teff 6500 K",
        verdict="PARTIALLY_HELD",
        blocker="we hold data/nlte_grids/amarsi2019_cno/ for C/N/O; the Fe II leg of the "
                "SAME release is not registered against Fe",
        note="Fe II only. THREED_HOLDINGS maps this release to C/N/O and not to Fe.",
    ),
    dict(
        route="mpia-siu-mean3d",
        scale="<3D>-NLTE",
        product="MPIA Spectrum Tools online abundance-correction service",
        source="http://nlte.mpia.de/gui-siuAC_secE.php",
        citation="Bergemann et al. / MPIA SIU",
        access="UNAVAILABLE",
        coverage="n/a",
        verdict="ADVERTISED_DISABLED",
        blocker="the `3dmean` radio is COMMENTED OUT in the live page source; only 1D "
                "models are selectable (re-verified live, RYA-1035)",
        note="The page still LISTS 'mean 3D' in its model text, so a reader sees an "
             "option that cannot be selected.",
    ),
]


def build(offline: bool) -> dict:
    doc: dict = {
        "ticket": "RYA-1035",
        "step": "Step 0 -- fetch-check before build",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "offline": offline,
        "verdict": "HAVE -- the build in RYA-1035's 'only if Step 0 comes back empty' "
                   "branch is NOT scoped. The Fe <3D> row is a fetch-and-wire job.",
        "routes": ROUTES,
    }
    doc["aux_check"] = probe_aux()
    doc["mlp_check"] = probe_mlp()
    if not offline:
        doc["keeper_check"] = probe_keeper()
        doc["amarsi_check"] = probe_amarsi()
    return doc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true",
                    help="skip the live Keeper and Amarsi probes; run the in-repo legs only")
    ap.add_argument("--out-dir", default=str(OUT_DIR),
                    help="where to write the census artifacts")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    doc = build(args.offline)

    (out_dir / "fe_3d_route_census.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False) + "\n")

    cols = ["route", "scale", "product", "source", "citation", "access", "coverage",
            "verdict", "blocker", "note"]
    with open(out_dir / "fe_3d_route_census.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in ROUTES:
            w.writerow(r)

    print(f"wrote {out_dir}/fe_3d_route_census.{{csv,json}}")
    for r in ROUTES:
        print(f"  {r['verdict']:<20} {r['route']:<24} {r['scale']}")
    print(f"\nVERDICT: {doc['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
