#!/usr/bin/env python3
"""Run a resolved plan, verifying POSTCONDITIONS rather than exit codes (RYA-767/941).

Two drivers, one contract:

  --driver local   plain Python. The reference implementation and the control.
  --driver qwen    the same steps, with a local model on Sirius asked to read
                   each step and return a verdict. The model does not choose
                   what runs, in what order, or what "done" means -- all of that
                   is already in the plan.

That split is not stylistic. Measured on Sirius (RYA-941): the model runs at
**3.0 tokens/s** warm, with a ~110 s cold load. A driver that asked it to reason
about coverage or compose commands would cost minutes per step and would be
wrong on exactly the questions that matter -- a 7B model answers "does Kurucz
2005 cover the near-infrared?" confidently and incorrectly. Asked instead to
confirm a named file exists with rows, it emits ~20 tokens in ~7 s and is right.

A step is done when its POSTCONDITION holds. Never when it exited 0: RYA-682 is
the standing counter-example, where iSpec writes a zero-row artifact and exits 0,
so a clean exit is evidence of nothing.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Long enough that the ~110 s cold load is paid once per session, not per step.
OLLAMA_KEEP_ALIVE = "60m"
OLLAMA_URL = "http://localhost:11434/api/generate"


def postcondition_holds(step: dict, *, search_roots: list[Path]) -> tuple[bool, str]:
    """Did this step actually produce something? Checked on the FILESYSTEM.

    Deliberately not a model's opinion and not a return code. The plan names the
    artifact; either it is there with content, or the step is not done.
    """
    produced = step.get("produces")
    if not produced:
        return False, "step declares no artifact, so nothing can prove it ran"
    name = Path(produced).name
    for root in search_roots:
        for hit in root.rglob(name):
            try:
                rows = sum(1 for _ in hit.open()) - 1        # minus the header
            except OSError as exc:
                return False, f"{hit} unreadable: {exc}"
            if rows < 1:
                return False, (f"{hit} exists but holds {rows} data rows. A zero-row "
                               f"artifact is a FAILURE, not an empty result (RYA-682).")
            return True, f"{hit} ({rows} rows)"
    return False, f"no file named {name} under {[str(r) for r in search_roots]}"


def ask_qwen(model: str, step: dict, evidence: str, ok: bool) -> dict:
    """Give the model the verdict to confirm, and a strict shape to answer in.

    It is shown what was checked and what was found. It is NOT asked to decide
    whether the science is right, which holding to use, or what to run next --
    those were resolved before this plan was written.
    """
    prompt = (
        "You are verifying one step of a pre-computed pipeline plan. "
        "You do not choose what runs. Answer ONLY with compact JSON.\n\n"
        f"step: {step['name']}\n"
        f"expected artifact: {step.get('produces')}\n"
        f"postcondition: {step.get('postcondition')}\n"
        f"filesystem check result: {'PASS' if ok else 'FAIL'}\n"
        f"evidence: {evidence}\n\n"
        'Reply exactly: {"verdict":"pass"} or {"verdict":"fail","why":"<8 words>"}'
    )
    body = json.dumps({
        "model": model, "prompt": prompt, "stream": False,
        "keep_alive": OLLAMA_KEEP_ALIVE,
        "options": {"num_predict": 40, "temperature": 0},
    }).encode()
    started = time.time()
    request = urllib.request.Request(OLLAMA_URL, data=body,
                                     headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=600) as response:
        payload = json.load(response)
    text = (payload.get("response") or "").strip()
    parsed: dict | None = None
    start, end = text.find("{"), text.rfind("}")
    if start >= 0 and end > start:
        try:
            parsed = json.loads(text[start:end + 1])
        except json.JSONDecodeError:
            parsed = None
    return {"raw": text[:200], "parsed": parsed,
            "elapsed_s": round(time.time() - started, 1),
            "eval_count": payload.get("eval_count"),
            "tok_per_s": round((payload.get("eval_count") or 0)
                               / max(payload.get("eval_duration", 1) / 1e9, 1e-9), 2)}


def run_step(step: dict, *, repo: Path, dry: bool) -> int:
    cmd = [step.get("interpreter") or sys.executable, step["script"], *step["args"]]
    if dry:
        print(f"      DRY: {' '.join(cmd)}")
        return 0
    import os
    env = dict(os.environ, **{k: v for k, v in (step.get("env") or {}).items() if v})
    return subprocess.run(cmd, cwd=repo, env=env).returncode


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("plan", type=Path)
    ap.add_argument("--driver", choices=["local", "qwen"], default="local")
    ap.add_argument("--model", default="qwen2.5-coder:latest")
    ap.add_argument("--repo", type=Path, default=ROOT)
    ap.add_argument("--search-root", type=Path, action="append", default=None)
    ap.add_argument("--only", default=None,
                    help="substring filter on the run key, e.g. FeII")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report", type=Path, default=None)
    args = ap.parse_args()

    plan = json.loads(args.plan.read_text())
    roots = args.search_root or [args.repo / "data" / "measured" / "band_ew",
                                 args.repo / "data" / "results"]
    runs = [r for r in plan["runs"]
            if not args.only or args.only in r["descriptor"]["key"]]
    print(f"plan: {len(plan['runs'])} runnable, {plan['n_blocked']} blocked; "
          f"executing {len(runs)} (driver={args.driver})")
    for b in plan.get("blocked", []):
        print(f"  BLOCKED (reported, not retried) {b['descriptor']['key'][:58]}")

    report = {"driver": args.driver, "model": args.model if args.driver == "qwen" else None,
              "runs": []}
    for i, run in enumerate(runs, 1):
        key = run["descriptor"]["key"]
        print(f"\n[{i}/{len(runs)}] {key}  ({run['method']})")
        entry = {"key": key, "steps": []}
        for step in run["steps"]:
            ok, evidence = postcondition_holds(step, search_roots=roots)
            if ok:
                # Resumability. Re-running a satisfied step wastes an hour of synthesis.
                print(f"    SKIP {step['name']} -- already satisfied: {evidence}")
                entry["steps"].append({"name": step["name"], "action": "skipped",
                                       "evidence": evidence})
                continue
            print(f"    RUN  {step['name']}")
            rc = run_step(step, repo=args.repo, dry=args.dry_run)
            ok, evidence = postcondition_holds(step, search_roots=roots)
            record = {"name": step["name"], "action": "ran", "returncode": rc,
                      "postcondition_holds": ok, "evidence": evidence}
            if args.driver == "qwen":
                verdict = ask_qwen(args.model, step, evidence, ok)
                record["qwen"] = verdict
                said = (verdict["parsed"] or {}).get("verdict")
                print(f"    QWEN {said!r}  ({verdict['elapsed_s']}s, "
                      f"{verdict['tok_per_s']} tok/s)")
                if said == "pass" and not ok:
                    # The filesystem wins. The model is a second opinion that can be
                    # wrong; it is never allowed to certify an artifact into existence.
                    record["disagreement"] = ("model said pass, postcondition FAILED -- "
                                              "filesystem is authoritative")
                    print("    !! model disagreed with the filesystem; filesystem wins")
            entry["steps"].append(record)
            if not ok:
                entry["stopped"] = f"postcondition failed: {evidence}"
                print(f"    STOP {evidence}")
                break
        report["runs"].append(entry)

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(f"\nreport -> {args.report}")


if __name__ == "__main__":
    main()
