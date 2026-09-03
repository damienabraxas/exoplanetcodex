#!/usr/bin/env python3
"""RYA-1055 — measure which IONISATION STAGES a Gerber model atom can serve in NLTE.

    python3 scripts/rya1055_atom_ion_reach.py --atom atom.fe607a \
        --out data/results/rya1055/atom_ion_reach.json

WHAT THIS ANSWERS, AND WHY IT IS NOT THE QUESTION ANYONE WAS ASKING
-------------------------------------------------------------------
bsyn applies departures PER LINE and falls back to departure = 1 for any line whose two
levels are not identified in the atom. So the question "can this deck do NLTE for
{element} {ion}" is decided by ONE property of the atom file: does it declare any
bound-bound transition BOTH of whose levels belong to that ionisation stage?

RYA-1055's first version asked a LINE-LIST question instead — "how many Fe II lines in
our production list carry NLTE level labels" — measured 0 of 11 on the pool, and proposed
labelling a VALD Fe II list against the atom. That project could not have worked, and
this script is the cheap bounding computation that says why in one number.

`atom.fe607a`:  607 levels (548 Fe I, 58 Fe II, 1 Fe III), 12,635 bound-bound
transitions, and every single one of them is Fe I -> Fe I. The 58 Fe II levels are a
pure IONISATION RESERVOIR — the targets of Fe I photoionisation, which in turn ionise to
Fe III. There is nothing in the atom for an Fe II line to point at, so no line list,
VALD or otherwise, can enable Fe II NLTE on this deck.

WHY BOTH HALVES ARE REPORTED
-----------------------------
`n_bb_both` (both levels in the stage) is what decides a product. `n_bb_any` (at least
one level in the stage) is reported beside it because the two numbers can differ — a
transition spanning stages would be a parse fault, not physics — and a bare zero on the
first is not distinguishable from "the parser found nothing at all". `n_parsed` vs the
header's declared count is the third guard: a SHORT READ would report zero for every
stage and look exactly like this result. It refuses instead.

⚠️ THE ATOM IS A VENDOR DECK, NOT A REPO FILE. It is resolved through the RYA-810 path
register (`grids.gerber_ts`), so this runs wherever the deck is staged and fails loudly
where it is not. The committed artifact is the record for machines without it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import codex_path, require_codex_path  # noqa: E402

#: A level line: `  energy_cm   g   'label'   ion_stage`.
_LEVEL = re.compile(r"^\s*(-?[\d.]+)\s+([\d.eE+-]+)\s+'([^']*)'\s+(\d+)\s*$")
#: A bound-bound line: `  UP  LOW  f  nq  qmax  q0  iw  ga  gv  gq  PROFILE`.
_TRANS = re.compile(r"^\s*(\d+)\s+(\d+)\s+([-\d.eE+]+)\s+")


class AtomParseError(RuntimeError):
    """Raised rather than returning a partial atom. A short read reports zero for every
    stage, which is indistinguishable from the real finding this script exists to make."""


def measure(path: Path) -> dict:
    """Ion-stage reach of one Gerber `atom.*` file. Every number is read, none assumed."""
    text = path.read_text(errors="replace").splitlines()
    element = text[0].strip().split()[0] if text and text[0].strip() else ""
    if not element:
        raise AtomParseError(f"{path}: first line carries no element name")

    n_lev = n_trans = n_cont = None
    hdr = -1
    for i, line in enumerate(text):
        parts = line.split()
        if len(parts) >= 3 and all(p.lstrip("-").isdigit() for p in parts[:3]):
            n_lev, n_trans, n_cont = (int(parts[0]), int(parts[1]), int(parts[2]))
            hdr = i
            break
    if n_lev is None:
        raise AtomParseError(f"{path}: no counts header (n_levels n_transitions n_continua)")

    # ── levels ────────────────────────────────────────────────────────────────────
    stage: dict[int, int] = {}
    energy_eV: dict[int, float] = {}
    j = hdr + 1
    while len(stage) < n_lev and j < len(text):
        m = _LEVEL.match(text[j].rstrip())
        if m:
            idx = len(stage) + 1
            stage[idx] = int(m.group(4))
            energy_eV[idx] = float(m.group(1)) / 8065.543937
        j += 1
    if len(stage) != n_lev:
        raise AtomParseError(
            f"{path}: header declares {n_lev} levels, parsed {len(stage)}. Refusing: a "
            f"short read reports zero transitions for every stage.")

    # ── bound-bound transitions ───────────────────────────────────────────────────
    both: Counter = Counter()
    any_: Counter = Counter()
    cross = 0
    n_parsed = 0
    max_index = 0
    for line in text[j:]:
        if not line.strip() or line.lstrip().startswith(("*", "#")):
            continue
        m = _TRANS.match(line)
        if not m:
            continue
        up, lo = int(m.group(1)), int(m.group(2))
        if not (1 <= up <= n_lev and 1 <= lo <= n_lev):
            continue
        n_parsed += 1
        max_index = max(max_index, up, lo)
        su, sl = stage[up], stage[lo]
        if su == sl:
            both[su] += 1
        else:
            cross += 1
        any_[su] += 1
        if sl != su:
            any_[sl] += 1
        if n_parsed >= n_trans:
            break
    if n_parsed != n_trans:
        raise AtomParseError(
            f"{path}: header declares {n_trans} bound-bound transitions, parsed "
            f"{n_parsed}. Refusing rather than reporting a partial reach.")

    stages = sorted(set(stage.values()))
    per_stage = []
    for s in stages:
        idx = [i for i, v in stage.items() if v == s]
        per_stage.append({
            "ion_stage": s,
            "n_levels": len(idx),
            "level_index_lo": min(idx),
            "level_index_hi": max(idx),
            "energy_eV_lo": round(min(energy_eV[i] for i in idx), 4),
            "energy_eV_hi": round(max(energy_eV[i] for i in idx), 4),
            "n_bb_both": int(both.get(s, 0)),
            "n_bb_any": int(any_.get(s, 0)),
            # THE VERDICT, and it is `n_bb_both` alone that decides it: a line needs two
            # levels of ONE stage. A stage with levels but no transitions among them is a
            # reservoir closing the ionisation balance, not a species that can be
            # synthesised.
            "nlte_capable": bool(both.get(s, 0)),
        })
    return {
        "atom": path.name,
        "element": element,
        # md5 BESIDE sha256 deliberately: md5 is the identifier the deck's own
        # `Fe_gerber2023.prov.json` and RYA-1035 already pin, so this artifact can be tied
        # to the STAGED bytes without re-hashing a 6.7 MB vendor file by hand.
        "md5": hashlib.md5(path.read_bytes()).hexdigest(),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n_levels_declared": n_lev,
        "n_bb_declared": n_trans,
        "n_bb_parsed": n_parsed,
        "n_continua_declared": n_cont,
        "n_cross_stage_bb": cross,
        "highest_level_index_in_any_bb": max_index,
        "stages": per_stage,
        "measured_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def _resolve_atom(name: str) -> Path:
    """The deck dir through the RYA-810 register, tolerating a root OVERRIDE.

    ⚠️ `require_codex_path` checks the `grids` root's REMOVABLE mount + sentinel before
    anything else, and those are Sirius literals (`/mnt/codex-ext`) read from the register
    rather than from the resolved root. So a machine that legitimately overrides
    `CODEX_GRID_ROOT` — the Mac, where the deck is staged locally — is told the volume is
    unmounted while the file is sitting right there. Resolve first, and fall back to
    `require_codex_path` ONLY when the resolved directory is genuinely absent, so its
    stale-mount diagnosis is still what a Sirius reader gets.
    """
    d = codex_path("grids.gerber_ts")
    if not d.is_dir():
        require_codex_path("grids.gerber_ts")       # raises with the right explanation
    p = Path(name)
    if not p.is_absolute():
        p = d / name
    if not p.exists():
        raise SystemExit(f"model atom not found: {p}\n"
                         f"  deck dir: {d}  (key 'grids.gerber_ts')\n"
                         f"  There is NO local fallback (RYA-567): stage it, or run on "
                         f"Sirius.")
    return p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atom", default="atom.fe607a",
                    help="atom file NAME inside the gerber_ts deck (default atom.fe607a)")
    ap.add_argument("--out", default="data/results/rya1055/atom_ion_reach.json")
    ap.add_argument("--check", action="store_true",
                    help="re-measure and compare against the committed artifact; "
                         "exit 1 on any difference other than the timestamp")
    a = ap.parse_args(argv)

    got = measure(_resolve_atom(a.atom))
    out = ROOT / a.out
    if a.check:
        have = json.loads(out.read_text())
        g = {k: v for k, v in got.items() if k != "measured_at"}
        h = {k: v for k, v in have.items() if k != "measured_at"}
        if g != h:
            print("DIFFERS from the committed artifact:", file=sys.stderr)
            for k in sorted(set(g) | set(h)):
                if g.get(k) != h.get(k):
                    print(f"  {k}: committed={h.get(k)!r} measured={g.get(k)!r}",
                          file=sys.stderr)
            return 1
        print(f"{a.atom}: matches {a.out}")
        return 0

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(got, indent=2) + "\n")
    print(f"{got['atom']}  {got['element']}  "
          f"{got['n_levels_declared']} levels / {got['n_bb_parsed']} bound-bound")
    for s in got["stages"]:
        print(f"  stage {s['ion_stage']}: {s['n_levels']:>4} levels "
              f"(idx {s['level_index_lo']}-{s['level_index_hi']}, "
              f"{s['energy_eV_lo']:.3f}-{s['energy_eV_hi']:.3f} eV)  "
              f"bound-bound both={s['n_bb_both']:>6}  any={s['n_bb_any']:>6}  "
              f"NLTE-capable={s['nlte_capable']}")
    print(f"  cross-stage bound-bound: {got['n_cross_stage_bb']}")
    print(f"  highest level index in any bound-bound transition: "
          f"{got['highest_level_index_in_any_bb']}")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
