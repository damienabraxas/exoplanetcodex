"""
pipeline/intake_debug.py  (RYA-318 MVP -- the diagnostics-layer trace)
=====================================================================
A LOCAL, git-ignored, per-run trace of what happens when an element is taken
into the pipeline for a given product. It answers one question after the fact:
*where did this run silently default?*

WHY THIS AND NOT THE EXISTING LEDGERS
-------------------------------------
rejection_ledger.py (RYA-429), model_attempt_ledger.py (RYA-695) and
provenance_honesty.py (RYA-653) are COMMITTED science accounting -- git-tracked
on purpose. This is NOT one of those. It is an EPHEMERAL DEBUG TRACE written to
debug/intake/ (git-ignored, per Ryan): a chronological event stream for a single
run -- stage entered, asset expected/found, decision taken, and every FALLBACK
as it fires -- so a wrong number becomes a timeline you can read.

INVARIANTS
----------
1. NUMERICALLY INERT. The tracer only reads and records. A run with tracing OFF
   and the same run with tracing ON must produce byte-identical science outputs.
2. OPT-IN. Off unless explicitly started (start_trace / CODEX_INTAKE_TRACE=1).
   The deep-site helpers no-op when no trace is active, so instrumentation left
   in place costs nothing in production.
3. LOUD-BUT-SAFE ON ITS OWN FAILURE. Observability must not break science: if the
   tracer cannot write, it warns ONCE and disables itself -- it never raises into
   the pipeline. (The one correct swallow: of the tracer's OWN error, announced.)
4. DE-SILENCING FALLBACKS IS ITS WHOLE JOB. Every .fallback() is the pipeline
   admitting it defaulted; the severity says how much that should worry you.
"""
from __future__ import annotations

import contextvars
import json
import math
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TRACE_DIR = ROOT / 'debug' / 'intake'   # git-ignored; see .gitignore

INFO, WARN, ERROR = 'INFO', 'WARN', 'ERROR'
_RANK = {INFO: 0, WARN: 1, ERROR: 2}
_ICON = {INFO: 'ok', WARN: 'WARN', ERROR: 'ERR'}

_CURRENT: 'contextvars.ContextVar[IntakeTrace | None]' = \
    contextvars.ContextVar('codex_intake_trace', default=None)


def _jsonable(v):
    """JSON-safe. NaN/inf -> string tag (bare NaN is invalid JSON, breaks readers)."""
    if isinstance(v, float):
        if math.isnan(v):
            return 'NaN'
        if math.isinf(v):
            return 'Infinity' if v > 0 else '-Infinity'
        return v
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (str, int, bool)) or v is None:
        return v
    return str(v)   # paths, numpy scalars -> repr, never crash


class IntakeTrace:
    """One run's event stream at debug/intake/{run_id}.jsonl. Build via start_trace()."""

    def __init__(self, element, product, star, *, run_id=None, extra=None):
        self.element = str(element)
        self.product = str(product)
        self.star = str(star)
        ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        safe = f"{self.element}_{self.product}_{self.star}".replace('/', '-').replace(' ', '')
        self.run_id = run_id or f"{safe}_{ts}"
        self.started = time.time()
        self.path = TRACE_DIR / f"{self.run_id}.jsonl"
        self._fh = None
        self._disabled = False
        self._counts = {INFO: 0, WARN: 0, ERROR: 0}
        self._fallbacks = []
        self._missing = []
        self._extra = dict(extra or {})
        try:
            TRACE_DIR.mkdir(parents=True, exist_ok=True)
            self._fh = self.path.open('w', encoding='utf-8')
        except OSError as exc:
            self._degrade(f"cannot open trace file {self.path}: {exc}")
        self.event('run', 'start', INFO, phase='preflight', name_over=self.run_id,
                   detail=f"{self.element} / {self.product} / {self.star}",
                   data=self._extra)

    def _degrade(self, why):
        if not self._disabled:
            msg = f"[intake_debug] tracing disabled: {why}"
            # BOTH channels, deliberately. `warnings.warn` alone is not loud enough to
            # satisfy invariant #3: the first driver wired to this tracer
            # (scripts/rya527_two_engine_run.py) calls warnings.filterwarnings('ignore')
            # at import, and several harnesses do the same, so a warn-only degrade is
            # SWALLOWED -- the tracer would fail silently, which is the exact failure
            # class this module exists to de-silence. Found by RYA-765 proof (c): the
            # run completed and produced correct output, and nothing said tracing had
            # stopped. stderr is not filterable by a warnings filter; the warn is kept
            # so pytest.warns / -W still see it.
            print(msg, file=sys.stderr, flush=True)
            warnings.warn(msg, stacklevel=2)
        self._disabled = True
        self._fh = None

    def event(self, event, name, severity=INFO, *, phase='inline', stage=None,
              ion=None, detail='', data=None, name_over=None):
        if severity not in _RANK:
            severity = INFO
        self._counts[severity] = self._counts.get(severity, 0) + 1
        rec = dict(
            t=round(time.time() - self.started, 4),
            run_id=self.run_id, phase=phase, event=event, severity=severity,
            stage=stage, element=self.element, ion=ion, product=self.product,
            star=self.star, name=str(name_over if name_over is not None else name),
            detail=str(detail), data=_jsonable(data) if data else None,
        )
        if self._disabled or self._fh is None:
            return rec
        try:
            self._fh.write(json.dumps(rec, ensure_ascii=True) + '\n')
            self._fh.flush()   # a mid-run crash must leave the trace up to the failure
        except (OSError, ValueError) as exc:
            self._degrade(f"write failed: {exc}")
        return rec

    def asset(self, name, present, *, path=None, detail='', phase='preflight'):
        """Expected per-element asset (grd, label, linelist, nlteinfo, ...).
        Missing = ERROR and loud -- a skipped asset is the #1 intake break."""
        sev = INFO if present else ERROR
        if not present:
            self._missing.append(name)
        return self.event('asset', name, sev, phase=phase,
                          detail=detail or ('present' if present else
                                            'MISSING -- expected asset not found'),
                          data={'present': bool(present), 'path': path})

    def decision(self, name, choice, *, detail='', phase='inline', **data):
        return self.event('decision', name, INFO, phase=phase,
                          detail=detail or f"chose {choice}",
                          data={'choice': choice, **data})

    def fallback(self, name, detail, *, severity=WARN, phase='inline', ion=None, **data):
        """THE point: the pipeline admitting it defaulted (silent LTE, NaN gf,
        all-zero iSpec, inherited default broadening)."""
        self._fallbacks.append((name, severity))
        return self.event('fallback', name, severity, phase=phase, ion=ion,
                          detail=detail, data=data or None)

    def check(self, name, ok, *, severity_if_fail=ERROR, phase='postflight',
              detail='', **data):
        sev = INFO if ok else severity_if_fail
        return self.event('check', name, sev, phase=phase,
                          detail=detail or ('pass' if ok else 'FAIL'),
                          data={'ok': bool(ok), **data})

    def stage(self, name):
        return _Stage(self, name)

    def summary(self, *, print_report=True):
        health = 'red' if self._counts[ERROR] else ('yellow' if self._counts[WARN] else 'green')
        s = dict(
            run_id=self.run_id, element=self.element, product=self.product,
            star=self.star, health=health, n_error=self._counts[ERROR],
            n_warn=self._counts[WARN], n_info=self._counts[INFO],
            missing_assets=list(self._missing),
            fallbacks=[n for n, _ in self._fallbacks],
            elapsed_s=round(time.time() - self.started, 3), path=str(self.path),
        )
        self.event('run', 'summary', INFO, phase='postflight', data=s)
        if self._fh is not None:
            try:
                self._fh.close()
            except OSError:
                pass
            self._fh = None
        if print_report:
            _print_health(s)
        return s


class _Stage:
    def __init__(self, trace, name):
        self.trace, self.name, self.t0 = trace, name, None
    def __enter__(self):
        self.t0 = time.time()
        self.trace.event('stage_enter', self.name, INFO, stage=self.name)
        return self.trace
    def __exit__(self, exc_type, exc, tb):
        dt = round(time.time() - (self.t0 or time.time()), 3)
        if exc_type is not None:
            self.trace.event('stage_exit', self.name, ERROR, stage=self.name,
                             detail=f"raised {exc_type.__name__}: {exc}",
                             data={'elapsed_s': dt})
        else:
            self.trace.event('stage_exit', self.name, INFO, stage=self.name,
                             data={'elapsed_s': dt})
        return False   # never swallow a pipeline exception


def start_trace(element, product, star, **kw) -> IntakeTrace:
    """Begin a trace and make it current for this context. Call once at the top
    of an element run. Retrievable via current()."""
    tr = IntakeTrace(element, product, star, **kw)
    _CURRENT.set(tr)
    return tr


def current():
    return _CURRENT.get()


def active() -> bool:
    return _CURRENT.get() is not None


def end_trace(*, print_report=True):
    tr = _CURRENT.get()
    if tr is None:
        return None
    s = tr.summary(print_report=print_report)
    _CURRENT.set(None)
    return s


# Deep-site helpers: safe anywhere, no-op when no trace is active, so they can
# live PERMANENTLY in nlte_corrections / gf_resolver / etc. at zero prod cost.
def trace_fallback(name, detail, *, severity=WARN, **data):
    tr = _CURRENT.get()
    if tr is not None:
        tr.fallback(name, detail, severity=severity, **data)

def trace_decision(name, choice, *, detail='', **data):
    tr = _CURRENT.get()
    if tr is not None:
        tr.decision(name, choice, detail=detail, **data)

def trace_asset(name, present, *, path=None, detail=''):
    tr = _CURRENT.get()
    if tr is not None:
        tr.asset(name, present, path=path, detail=detail)

def trace_check(name, ok, *, severity_if_fail=ERROR, detail='', **data):
    tr = _CURRENT.get()
    if tr is not None:
        tr.check(name, ok, severity_if_fail=severity_if_fail, detail=detail, **data)


def _print_health(s):
    tag = {'green': 'GREEN', 'yellow': 'YELLOW', 'red': 'RED'}[s['health']]
    print(f"\nRUN HEALTH [{tag}] -- {s['element']} / {s['product']} / {s['star']}")
    print(f"  errors {s['n_error']}  warn {s['n_warn']}  info {s['n_info']}"
          f"   ({s['elapsed_s']}s)")
    if s['missing_assets']:
        print(f"  missing assets: {', '.join(s['missing_assets'])}")
    if s['fallbacks']:
        print(f"  fallbacks fired: {', '.join(s['fallbacks'])}")
    print(f"  trace: {s['path']}")


def show(path):
    """Print a run's trace as a readable timeline + its health summary."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"no such trace: {p}")
    summary = None
    print(f"# {p.name}")
    for line in p.read_text(encoding='utf-8').splitlines():
        if not line.strip():
            continue
        try:
            r = json.loads(line)
        except ValueError:
            print(f"  ??  unparseable: {line[:80]}")
            continue
        if r.get('event') == 'run' and r.get('name') == 'summary':
            summary = r.get('data')
            continue
        sev = _ICON.get(r.get('severity'), '?')
        stg = f"[{r['stage']}] " if r.get('stage') else ''
        det = f" -- {r['detail']}" if r.get('detail') else ''
        print(f"  {r.get('t'):>8}  {sev:<4} {r.get('event'):<12} {stg}{r.get('name')}{det}")
    if summary:
        _print_health(summary)


def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description='Codex intake debug trace reader (RYA-318).')
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--show', metavar='PATH', help='print a specific trace file')
    g.add_argument('--last', action='store_true', help='print the most recent trace')
    g.add_argument('--list', action='store_true', help='list traces, newest last')
    a = ap.parse_args(argv)
    if a.list:
        if not TRACE_DIR.exists():
            raise SystemExit(f"no trace dir yet: {TRACE_DIR}")
        for q in sorted(TRACE_DIR.glob('*.jsonl'), key=lambda q: q.stat().st_mtime):
            print(q)
        return
    if a.last:
        if not TRACE_DIR.exists():
            raise SystemExit(f"no trace dir yet: {TRACE_DIR}")
        files = sorted(TRACE_DIR.glob('*.jsonl'), key=lambda q: q.stat().st_mtime)
        if not files:
            raise SystemExit(f"no traces in {TRACE_DIR}")
        show(files[-1])
        return
    show(a.show)


if __name__ == '__main__':
    _main()
