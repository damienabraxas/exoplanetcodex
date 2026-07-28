#!/usr/bin/env python3
"""
scripts/gen_state_register_targets.py  (RYA-516)
================================================
Regenerate the "Targets & stellar parameters" MIRROR table in
CODEX_STATE_REGISTER.md **from** config/constants.py (STAR_PARAMS <- stars.yaml).

The Codex State Register distinguishes NATIVE rows (the register IS the source of
truth) from MIRROR rows (a fast view over another single source). Stellar
parameters are a MIRROR of constants.py: this script is the mechanism that keeps
them a *generated view*, never a hand-typed second copy that can silently drift.

Usage:
    python scripts/gen_state_register_targets.py            # print table to stdout
    python scripts/gen_state_register_targets.py --write    # splice into the register
    python scripts/gen_state_register_targets.py --check    # exit 1 if the register is stale

The generated block is delimited in the register by:
    <!-- BEGIN GENERATED: targets (scripts/gen_state_register_targets.py) -->
    <!-- END GENERATED: targets -->
Everything outside those markers (NATIVE caveats, ratification flags) is hand-kept.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import STAR_PARAMS  # noqa: E402  (needs sys.path first)

REGISTER = ROOT / "CODEX_STATE_REGISTER.md"
BEGIN = "<!-- BEGIN GENERATED: targets (scripts/gen_state_register_targets.py) -->"
END = "<!-- END GENERATED: targets -->"

# Presentation only (display names + intended science order). NOT parameter data:
# every numeric value below is read live from STAR_PARAMS, nothing is transcribed.
DISPLAY = {
    "solar": "Sun (Star Zero)",
    "procyon": "Procyon",
    "alpha_cen_a": "α Cen A",
    "alpha_cen_b": "α Cen B",
    "55cnc_a": "55 Cnc A",
    "synthetic_no_logg": "synthetic_no_logg (test fixture)",
}
TARGET_ORDER = ["solar", "procyon", "alpha_cen_a", "alpha_cen_b", "55cnc_a"]


def _fmt_feh(rec: dict) -> str:
    return f"{rec['feh_ref']:+.2f}"


def _fmt_xi(rec: dict) -> str:
    if "xi" in rec:
        return f"{rec['xi']:.2f} (pinned)"
    if "xi_init" in rec:
        xc = rec.get("xi_xcheck")
        xc_s = f"; x-check {xc[0]}–{xc[1]}" if xc else ""
        return f"{rec['xi_init']:.2f} (solved{xc_s})"
    return "—"


def _fmt_policy(rec: dict) -> str:
    pin = ",".join(rec.get("pin") or []) or "—"
    solve = ",".join(rec.get("solve") or []) or "—"
    return f"pin: {pin} · solve: {solve}"


def _ordered_keys() -> list[str]:
    keys = [k for k in TARGET_ORDER if k in STAR_PARAMS]
    keys += [k for k in STAR_PARAMS if k not in keys]  # any extras (fixtures) last
    return keys


def build_table() -> str:
    lines = [
        "| Star | Teff (K) | log g | [Fe/H] ref | ξ (km/s) | pin/solve policy | Source (verbatim from stars.yaml) |",
        "|---|---|---|---|---|---|---|",
    ]
    for key in _ordered_keys():
        rec = STAR_PARAMS[key]
        name = DISPLAY.get(key, key)
        teff = f"{rec['teff']:.0f} ± {rec.get('e_teff', 0):.0f}"
        logg = f"{rec['logg']}"  # natural float repr (4.0, 4.438) — no zero-strip to "4"
        src = str(rec.get("source", "—")).replace("\n", " ").replace("|", "/")
        lines.append(
            f"| {name} | {teff} | {logg} | {_fmt_feh(rec)} | {_fmt_xi(rec)} "
            f"| {_fmt_policy(rec)} | {src} |"
        )
    return "\n".join(lines)


def render_block() -> str:
    from config import constants  # for the source stamp

    stamp = getattr(constants, "__version__", "unknown")
    return (
        f"{BEGIN}\n"
        f"<!-- generated from config/constants.py STAR_PARAMS (stars.yaml, __version__={stamp}); "
        f"regenerate with: python scripts/gen_state_register_targets.py --write -->\n"
        f"{build_table()}\n"
        f"{END}"
    )


def splice(text: str, block: str) -> str:
    if BEGIN not in text or END not in text:
        raise SystemExit(
            f"markers not found in {REGISTER.name}; add the BEGIN/END pair first."
        )
    pre = text[: text.index(BEGIN)]
    post = text[text.index(END) + len(END) :]
    return pre + block + post


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="splice the block into the register")
    ap.add_argument("--check", action="store_true", help="exit 1 if the register block is stale")
    args = ap.parse_args()

    block = render_block()

    if not args.write and not args.check:
        print(block)
        return 0

    text = REGISTER.read_text()
    new = splice(text, block)

    if args.check:
        if new != text:
            print("STALE: CODEX_STATE_REGISTER.md Targets block differs from constants.py.")
            print("Run: python scripts/gen_state_register_targets.py --write")
            return 1
        print("OK: Targets block is in sync with constants.py.")
        return 0

    REGISTER.write_text(new)
    print(f"wrote generated Targets block into {REGISTER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
