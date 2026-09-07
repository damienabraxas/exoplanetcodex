#!/usr/bin/env python3
"""RYA-1181 (A9) — the NUV transitions are ON DISK; the intake recorded them "not published".

`rejected_indicator_ledger.csv` carried NH A-X (~340 nm), OH A-X (~320 nm) and CN B-X
(~390 nm) with `count = NOT_PUBLISHED` and reason "individual list not published". That
conflates two different facts:

    what is unpublished   WHICH SUBSET Amarsi used for the abundance
    what is held          THE TRANSITIONS THEMSELVES, in this repo, acquired and unread

"not published" reads as "not available", and it is not. This script reads what is held and
states the negative precisely, so the ledger stops understating the repo's own holdings.

🔴 AND THE UV LISTS CARRY THE INTAKE'S OWN BLOCKER. RYA-1136 names missing rotational
identity as the reason the crossmatch cannot close. `NH-A-X-linelist.csv` publishes J', J'',
symmetry, branch, v', v'', N', N'', E_upper, E_lower, f AND A -- richer rotational identity
than any list the intake currently parses.

⚠️ THIS STAGES, IT DOES NOT RE-MATCH. Adding these transitions to `inventory()` would make
them candidates for the A-X targets and would move `join_status` on real rows -- a value
change, which this ticket forbids. What wiring them WOULD do is measured and reported at
the bottom instead, so the decision is Ryan's and is made on a number.
"""
from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

PRIMARY = ROOT / "data" / "reference" / "cno_molecular_primary"
AUDIT = ROOT / "data" / "audit" / "rya1136_cno_intake"
OUT = AUDIT / "held_uv_transitions_rya1181.csv"

#: The two loose UV line lists that sit beside the archives the ingest opens, and that the
#: ingest never touches: `parse_brooke_xx` reads only the X-X member INSIDE each zip.
NH_AX = PRIMARY / "nh_brooke2014" / "NH-A-X-linelist.csv"
OH_AX = PRIMARY / "oh_brooke2016" / "OH-A-X-linelist-final.csv"

NUV_A = (2000.0, 4000.0)          # the near-UV bin, in Angstrom


def read_nh_ax() -> list[dict]:
    """NH A-X. Comma-separated with a BOM; `Position(angair)` is air Angstrom."""
    rows = []
    with NH_AX.open(encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                lam = float(r["Position(angair)"])
            except (TypeError, ValueError, KeyError):
                continue
            rows.append({
                "species": "NH", "system": "A-X",
                "wavelength_air_A": f"{lam:.4f}",
                "j_lower": (r.get('J"') or "").strip(),
                "j_upper": (r.get("J'") or "").strip(),
                "branch": (r.get("Branch") or "").strip(),
                "sym_lower": (r.get('Sym"') or "").strip(),
                "e_upper": (r.get("Eupper") or "").strip(),
                "e_lower": (r.get("Elower") or "").strip(),
                "f_value": (r.get("f-value") or "").strip(),
                "einstein_A": (r.get("A") or "").strip(),
                "source": str(NH_AX.relative_to(ROOT)),
            })
    return rows


def read_oh_ax() -> list[dict]:
    """OH A-X. A prose preamble precedes the table; the header row is the one naming
    `wl_air`, and that column is NANOMETRES (309.37 nm), not Angstrom."""
    lines = OH_AX.read_text(errors="replace").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Molecule,Manifold"))
    rows = []
    for r in csv.DictReader(lines[start:]):
        try:
            lam_nm = float(r["wl_air"])
        except (TypeError, ValueError, KeyError):
            continue
        rows.append({
            "species": "OH", "system": "A-X",
            "wavelength_air_A": f"{lam_nm * 10.0:.4f}",
            "j_lower": (r.get('J"') or "").strip(),
            "j_upper": (r.get("J'") or "").strip(),
            "branch": "",
            "sym_lower": (r.get('sym"') or "").strip(),
            "e_upper": (r.get("Eupper") or "").strip(),
            "e_lower": (r.get("Elower") or "").strip(),
            "f_value": (r.get("f_ocs") or "").strip(),
            "einstein_A": (r.get("A") or "").strip(),
            "source": str(OH_AX.relative_to(ROOT)),
        })
    return rows


def parsed_but_unqueried() -> dict:
    """Systems the ingest PARSES and then never looks up.

    The match index is keyed on (species, system) and Amarsi's Table 2 carries no B-X or
    C-X row, so those transitions are read off disk into memory and silently go nowhere.
    They are not lost data -- but nothing counts them, which is how ~46k CN B-X lines can
    be in hand while a ledger records the band as "not published".
    """
    import ingest_cno_molecular_primary_rya1136 as M
    have = collections.Counter((t.species, t.system) for t in M.inventory())
    wanted = {(r["species"], r["system"])
              for r in csv.DictReader((AUDIT / "primary_molecular_crossmatch.csv").open())}
    return {f"{sp} {sysn}": n for (sp, sysn), n in sorted(have.items())
            if (sp, sysn) not in wanted}


def main() -> int:
    nh, oh = read_nh_ax(), read_oh_ax()
    rows = nh + oh
    AUDIT.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"=== RYA-1181 A9 — UV transitions HELD ON DISK ===\n")
    for name, rs, path in (("NH A-X", nh, NH_AX), ("OH A-X", oh, OH_AX)):
        lam = [float(r["wavelength_air_A"]) for r in rs]
        nuv = [x for x in lam if NUV_A[0] <= x <= NUV_A[1]]
        withj = sum(1 for r in rs if r["j_lower"])
        print(f"  {name:7s} {len(rs):6d} transitions   {min(lam):8.1f}-{max(lam):8.1f} A"
              f"   NUV {len(nuv):6d}   with J'': {withj}")
        print(f"          read from {path.relative_to(ROOT)} — never opened by the ingest")

    unq = parsed_but_unqueried()
    print(f"\n  PARSED BY THE INGEST AND NEVER QUERIED (index keyed on species+system,")
    print(f"  and Amarsi's Table 2 carries no such row):")
    for k, n in unq.items():
        print(f"    {k:10s} {n:7d} transitions")
    print(f"    {'TOTAL':10s} {sum(unq.values()):7d}")

    print(f"\n  staged -> {OUT.relative_to(ROOT)}  ({len(rows)} rows)")
    print(f"\n  ⚠️ NOT wired into inventory(): these would become candidates for the A-X")
    print(f"     targets and move join_status on live rows. That is a value change and")
    print(f"     belongs to a measurement ticket, not this one.")
    (AUDIT / "held_uv_summary_rya1181.json").write_text(json.dumps({
        "ticket": "RYA-1181",
        "nh_ax_transitions": len(nh), "oh_ax_transitions": len(oh),
        "parsed_but_unqueried": unq,
        "total_held_uv_not_previously_counted": len(rows) + sum(unq.values()),
        "negative_restated": ("Amarsi's UV SELECTION is unpublished; the TRANSITIONS are "
                              "held in this repo and are now counted and staged."),
    }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
