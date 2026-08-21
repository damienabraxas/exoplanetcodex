#!/usr/bin/env python3
"""Emit an explicit, resumable execution PLAN from a run request (RYA-767).

This is the seam between judgment and execution. Everything that requires
knowing something -- which holdings cover a band, which method a band permits,
which interpreter the engine needs, what order the steps go in, what proves a
step worked -- happens HERE, deterministically, and is written down. What comes
out is a plan an executor can follow while deciding nothing at all.

That split is the point. The eventual executor is a local model on Sirius, and a
local model asked to decide "does Kurucz 2005 cover the NIR?" will answer
confidently and wrongly. Asked to run step 3 of 24 and check that the named file
exists with at least one row, it will be right every time. So the plan carries
every decision already made, and every step carries the postcondition that
proves it -- never "did the command exit 0", which RYA-682 showed is not the
same question (iSpec writes a zero-row artifact and exits 0).

Usage:
    python scripts/rya767_plan.py --element Fe --ion I --ion II \
        --interpreter /mnt/codex-data/venv312/bin/python \
        --ispec-dir /mnt/codex-data/engines/ispec_src --out plan.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import band_policy                      # noqa: E402
from pipeline.run_descriptor import RunDescriptor, resolve  # noqa: E402


def candidate_runs(elements, ions, instruments):
    """Every (holding x band) an instrument actually offers -- enumerated from the
    harness's own table and the band policy, never from a list typed here."""
    from measure_band_ew import _INSTRUMENT_HOLDINGS
    for element in elements:
        for ion in ions:
            for instrument in instruments:
                for spec in _INSTRUMENT_HOLDINGS.get(instrument, ()):
                    for policy in band_policy.POLICIES:
                        yield RunDescriptor(element, ion, instrument, spec.holding_id,
                                            policy.lo_A, policy.hi_A)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element", action="append", required=True)
    ap.add_argument("--ion", action="append", required=True)
    ap.add_argument("--instrument", action="append", default=None)
    ap.add_argument("--interpreter", required=True,
                    help="MUST have numpy < 2.3 (RYA-682); the plan records it and the "
                         "executor is required to verify it before the first step")
    ap.add_argument("--ispec-dir", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--include-blocked", action="store_true",
                    help="keep runs that cannot happen, WITH their reason. On by default "
                         "in the report: a run we cannot do is a finding, and dropping it "
                         "silently is how 'we do not hold this' becomes indistinguishable "
                         "from 'nobody looked' (RYA-833).")
    args = ap.parse_args()

    instruments = args.instrument or ["harps", "kpno_solar_atlas"]
    runnable, blocked = [], []
    for descriptor in candidate_runs(args.element, args.ion, instruments):
        resolved = resolve(descriptor, interpreter=args.interpreter,
                           ispec_dir=args.ispec_dir)
        (runnable if resolved.runnable else blocked).append(resolved.as_dict())

    plan = {
        "ticket": "RYA-767",
        "elements": args.element, "ions": args.ion, "instruments": instruments,
        "interpreter": args.interpreter, "ispec_dir": args.ispec_dir,
        "executor_contract": {
            "decide_nothing": "Every decision in this plan is already made. An executor "
                              "that re-derives one has reintroduced the judgment this "
                              "layer exists to remove.",
            "verify_interpreter_first": "Check numpy < 2.3 before step 1. Above the "
                                        "ceiling iSpec writes a zero-row artifact and "
                                        "EXITS 0, so a clean exit proves nothing.",
            "postconditions_not_exit_codes": "A step succeeded when its postcondition "
                                             "holds. Exit code 0 is necessary and not "
                                             "sufficient.",
            "sequential": "Sirius is SHARED with concurrent sessions. Run steps in "
                          "order, one at a time, and never pattern-kill a process.",
            "resumable": "Skip any step whose postcondition already holds. Re-running a "
                         "completed step wastes an hour of synthesis.",
            "blocked_is_a_result": "The blocked list is part of the output, not an "
                                   "error. Report it; do not retry it.",
        },
        "n_runnable": len(runnable), "n_blocked": len(blocked),
        "runs": runnable,
        "blocked": blocked if args.include_blocked or True else [],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(plan, indent=2) + "\n")
    print(f"{len(runnable)} runnable, {len(blocked)} blocked -> {args.out}")
    for r in blocked:
        print(f"  BLOCKED {r['descriptor']['key'][:62]:<64} {r['blocked_reason'][:70]}")


if __name__ == "__main__":
    main()
