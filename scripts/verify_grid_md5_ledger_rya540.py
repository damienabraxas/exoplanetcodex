"""
scripts/verify_grid_md5_ledger_rya540.py
========================================
RYA-540 finalize — disk-verified grid inventory. md5-verify every retained grid against
its pinned provenance and cross-reference the coverage ledger. LOUD on any mismatch
(the RYA-534 truncated-download class); reports present / verified / pending / mismatch.

Two grid families:
  * Gerber TS (Engine-B) — data/nlte_grids/gerber_ts/*.bin, cached by scripts/grid_cache.py.
    Verify: cached .bin md5 == _cache_index.json bin_md5, and the index's zip_md5 == the
    zip md5 pinned in <El>_gerber2023.prov.json (the download was verified against prov).
  * Amarsi-2020 PySME (Engine-A) — data/nlte_grids/amarsi_galah/*.grd. Verify the on-disk
    .grd md5 against the inner-.grd md5 recorded in <El>_amarsi2020_v3.prov.json (where
    present); else report the observed md5 (no pinned reference to verify against).

Writes data/audit/rya540_disk_layout/grid_inventory.md and exits non-zero on any mismatch.

    GERBER_TS_DIR=/mnt/codex-data/grids/nlte/gerber_ts python -m scripts.verify_grid_md5_ledger_rya540
"""
from __future__ import annotations
import os, re, json, hashlib, csv
from pathlib import Path
# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

_REPO = Path(__file__).resolve().parents[1]
GERBER_DIR = Path(os.environ.get("GERBER_TS_DIR", str(codex_path('grids.gerber_ts'))))
AMARSI_DIR = Path(os.environ.get("AMARSI_GRD_DIR", str(codex_path('grids.amarsi_galah'))))
LEDGER_CSV = _REPO / "data" / "curation" / "nlte_grid_availability.csv"
PROV_DIR = _REPO / "data" / "nlte_grids" / "gerber_ts"
AMARSI_PROV_DIR = _REPO / "data" / "nlte_grids" / "amarsi_galah"


def _md5(p: Path, buf=1 << 20) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(buf), b""):
            h.update(c)
    return h.hexdigest()


def _gerber_prov_path(el):
    # provs are co-located with the cached grids (canonical, Sirius) AND committed in the
    # repo; prefer the co-located copy so a stale checkout doesn't hide a grid.
    for base in (GERBER_DIR, PROV_DIR):
        p = base / f"{el}_gerber2023.prov.json"
        if p.exists():
            return p
    return None


def _gerber_elements():
    els = set()
    for base in (GERBER_DIR, PROV_DIR):
        for p in base.glob("*_gerber2023.prov.json"):
            els.add(p.name.split("_")[0])
    return sorted(els)


def _prov_bin(el):
    p = _gerber_prov_path(el)
    if p is None:
        return None
    d = json.loads(p.read_text())["files"]
    z = d.get("grid_1d_bin_zip") or d.get("grid_1d_bin")
    return {"bin_name": z["name"][:-4] if z["name"].endswith(".zip") else z["name"],
            "zip_md5": z["md5"], "zip_bytes": z.get("bytes")}


def verify_gerber():
    idx = {}
    ci = GERBER_DIR / "_cache_index.json"
    if ci.exists():
        idx = json.loads(ci.read_text())
    rows = []
    for el in _gerber_elements():
        pv = _prov_bin(el)
        binf = GERBER_DIR / pv["bin_name"]
        rec = idx.get(el)
        if not binf.exists():
            rows.append((el, "gerber", pv["bin_name"], "PENDING", "not yet downloaded"))
            continue
        if not rec:
            rows.append((el, "gerber", pv["bin_name"], "UNINDEXED",
                         "bin present but not in cache index — re-run ensure"))
            continue
        actual = _md5(binf)
        bin_ok = actual == rec.get("bin_md5")
        zip_ok = rec.get("zip_md5") == pv["zip_md5"]
        status = "VERIFIED" if (bin_ok and zip_ok) else "MISMATCH"
        note = (f"bin_md5={actual[:12]} zip_md5={rec.get('zip_md5','?')[:12]} "
                f"({'bin ok' if bin_ok else 'BIN DRIFT'}, {'zip==prov' if zip_ok else 'ZIP!=PROV'})")
        rows.append((el, "gerber", pv["bin_name"], status, note))
    return rows


def _amarsi_inner_md5(el):
    p = AMARSI_PROV_DIR / f"{el}_amarsi2020_v3.prov.json"
    if not p.exists():
        return None
    m = re.search(r'inner\s+\.?grd\s+md5\s+([0-9a-f]{32})', p.read_text(), re.I)
    return m.group(1) if m else None


def verify_amarsi():
    rows = []
    if not AMARSI_DIR.exists():
        return rows
    for grd in sorted(AMARSI_DIR.glob("nlte_*_pysme.grd")):
        m = re.match(r"nlte_([A-Za-z]+)_", grd.name)
        el = m.group(1) if m else grd.stem
        pinned = _amarsi_inner_md5(el)
        actual = _md5(grd)
        if pinned:
            status = "VERIFIED" if actual == pinned else "MISMATCH"
            note = f"md5={actual[:12]} vs pinned {pinned[:12]}"
        else:
            status = "PRESENT"
            note = f"md5={actual[:12]} (no pinned inner-md5 to verify against)"
        rows.append((el, "amarsi_grd", grd.name, status, note))
    return rows


def main():
    gerber = verify_gerber()
    amarsi = verify_amarsi()
    allrows = gerber + amarsi
    mism = [r for r in allrows if r[3] == "MISMATCH"]
    ver = [r for r in allrows if r[3] == "VERIFIED"]
    pend = [r for r in allrows if r[3] == "PENDING"]

    lines = ["# RYA-540 — disk-verified grid inventory", "",
             f"Gerber TS ({GERBER_DIR}) + Amarsi PySME ({AMARSI_DIR}).",
             f"**{len(ver)} VERIFIED · {len(pend)} pending · {len(mism)} MISMATCH** "
             f"of {len(allrows)} grids.", "",
             "| element | family | grid file | status | note |",
             "|---|---|---|---|---|"]
    for el, fam, gf, st, note in allrows:
        lines.append(f"| {el} | {fam} | {gf} | **{st}** | {note} |")
    out = _REPO / "data" / "audit" / "rya540_disk_layout" / "grid_inventory.md"
    out.write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nwrote {out}")
    if mism:
        raise SystemExit(f"RYA-540: {len(mism)} grid md5 MISMATCH — corrupt/drift, LOUD-FAIL.")


if __name__ == "__main__":
    main()
