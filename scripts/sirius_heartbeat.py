#!/usr/bin/env python3
"""Heartbeat for long-running Sirius jobs — is it PROGRESSING, STUCK, or GONE?

WHY THIS EXISTS (RYA-821, 2026-08-23). A Fortran interpolator ran for minutes with no
output and no log (Fortran buffers stdout), and there was no way to tell a working job
from a hung one. Worse, the obvious watcher — `until ! pgrep -f <binary>` — MATCHED ITS
OWN ssh command line, so it never exited and I read the self-match as "still running".
Silence looked identical to progress, and a false signal looked identical to a real one.

TWO RULES THIS ENCODES:

1. **Match on `comm`, never `pgrep -f`.** `ps -eo comm` holds only the executable name,
   so a shell or ssh wrapper carrying the pattern can never appear in it. `pgrep -f` is
   the trap; `comm` is the cure. See [[feedback_waiter_self_match]].

2. **One sample cannot tell progress from a hang.** A process pinned at 99.9% CPU and a
   process spinning on a bad loop look identical in a single `ps`. The verdict comes from
   DELTAS between samples: cpu-seconds burned, bytes read, output bytes written. So this
   keeps a rolling history and reports change, not state.
"""
from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "heartbeat" / "sirius_heartbeat.json"
HOST = "sirius"

#: Executable-name prefixes worth watching, matched against `comm` (15-char truncated).
WATCH = ("interpol", "babsma", "bsyn", "python3", "turbospec", "moog", "molecfit")

#: Paths whose growth proves a job is producing something, not just burning CPU.
DEFAULT_ARTIFACTS = ("/tmp/multi_t/Testout", "/tmp/rya821/Testout")

_PS = (
    "ps -eo pid,comm,etimes,times,pcpu,rss --no-headers"
)


def _ssh(cmd: str, timeout: int = 60) -> str:
    r = subprocess.run(["ssh", "-o", "ConnectTimeout=15", HOST, cmd],
                       capture_output=True, text=True, timeout=timeout)
    return r.stdout


def sample(artifacts: tuple[str, ...]) -> dict:
    """One observation of every watched process plus the artifact sizes."""
    procs = []
    for line in _ssh(_PS).splitlines():
        f = line.split(None, 5)
        if len(f) < 6:
            continue
        pid, comm, etimes, times, pcpu, rss = f
        if not comm.startswith(WATCH):
            continue
        # io is best-effort: /proc/<pid>/io is unreadable for other users' processes.
        io = _ssh(f"grep -E '^(rchar|wchar)' /proc/{pid}/io 2>/dev/null || true", 30)
        rchar = wchar = None
        for l in io.splitlines():
            k, _, v = l.partition(":")
            if k == "rchar":
                rchar = int(v)
            elif k == "wchar":
                wchar = int(v)
        procs.append(dict(pid=int(pid), comm=comm, etimes=int(etimes),
                          cpu_s=int(times), pcpu=float(pcpu), rss_kb=int(rss),
                          rchar=rchar, wchar=wchar))

    art = {}
    for a in artifacts:
        q = shlex.quote(a)
        n = _ssh(f"du -sb {q} 2>/dev/null | cut -f1 || true", 30).strip()
        c = _ssh(f"find {q} -type f 2>/dev/null | wc -l || true", 30).strip()
        art[a] = dict(bytes=int(n) if n.isdigit() else None,
                      files=int(c) if c.isdigit() else None)

    return dict(t=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                procs=procs, artifacts=art)


def verdict(prev: dict | None, cur: dict) -> list[dict]:
    """PROGRESSING / STUCK / IDLE per process — from DELTAS, never one sample.

    A single observation cannot distinguish real work from a spin loop, so anything
    reported here without a previous sample is explicitly UNKNOWN rather than guessed.
    """
    out = []
    pmap = {p["pid"]: p for p in (prev or {}).get("procs", [])}
    for p in cur["procs"]:
        q = pmap.get(p["pid"])
        if q is None:
            state, why = "UNKNOWN", "first sample — a verdict needs two"
        else:
            dcpu = p["cpu_s"] - q["cpu_s"]
            drch = (p["rchar"] or 0) - (q["rchar"] or 0)
            dwch = (p["wchar"] or 0) - (q["wchar"] or 0)
            if dcpu > 0 or drch > 0 or dwch > 0:
                state = "PROGRESSING"
                why = f"+{dcpu}s cpu, +{drch}B read, +{dwch}B written since last sample"
            else:
                state = "STUCK"
                why = ("no cpu burned, no bytes moved since last sample — alive but "
                       "doing nothing")
        out.append(dict(pid=p["pid"], comm=p["comm"], elapsed_s=p["etimes"],
                        cpu_s=p["cpu_s"], pcpu=p["pcpu"], state=state, why=why))
    if not cur["procs"]:
        out.append(dict(pid=None, comm=None, state="GONE",
                        why="no watched process on the host — finished, or never started"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifact", action="append", default=None,
                    help="path whose growth proves output (repeatable)")
    ap.add_argument("--keep", type=int, default=200, help="samples to retain")
    ap.add_argument("--watch", type=int, default=0,
                    help="seconds between samples; 0 = take one and exit")
    args = ap.parse_args()
    artifacts = tuple(args.artifact or DEFAULT_ARTIFACTS)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    hist = []
    if OUT.exists():
        try:
            hist = json.loads(OUT.read_text()).get("samples", [])
        except (json.JSONDecodeError, OSError):
            hist = []

    while True:
        cur = sample(artifacts)
        v = verdict(hist[-1] if hist else None, cur)
        hist = (hist + [cur])[-args.keep:]
        OUT.write_text(json.dumps(
            dict(generated=cur["t"], host=HOST, generator="scripts/sirius_heartbeat.py",
                 verdict=v, samples=hist), indent=2) + "\n")
        for x in v:
            print(f"  {x['state']:12s} {x.get('comm') or '-':16s} "
                  f"pid={x.get('pid')} {x['why']}")
        for a, d in cur["artifacts"].items():
            print(f"  artifact {a}: {d['files']} files, {d['bytes']} bytes")
        if not args.watch:
            return
        time.sleep(args.watch)


if __name__ == "__main__":
    main()
