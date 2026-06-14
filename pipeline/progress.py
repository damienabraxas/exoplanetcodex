"""
pipeline/progress.py
Lightweight terminal progress + ETA reporter for long pipeline runs
(synthesis abundance runs are hours per star on Sirius/Orion).

Headless/SSH-friendly: discrete plain-text lines (no redrawing bar that
spams logs), throttled to <= one update per `min_interval` s, plus a
provisional ETA after a short warmup so the run can be left unattended.
Provisional ETA assumes roughly uniform per-item cost; it self-corrects
on every refresh as slower (e.g. escalated synthesis) lines come through.
All output goes to stderr so stdout/results stay clean.
"""
import sys
import time
from datetime import datetime, timedelta


def _fmt_dur(seconds: float) -> str:
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


class ProgressReporter:
    """
    prog = ProgressReporter(total=len(linemasks), label="Fe synthesis")
    for lm in linemasks:
        ... slow per-line work ...
        prog.update(item_label=f"{lm['element']} {lm['wave_A']:.2f}")
    prog.finish()
    """

    def __init__(self, total, label="run", stream=None,
                 min_interval=5.0, warmup=3):
        self.total = int(total) if total else None
        self.label = label
        self.stream = stream if stream is not None else sys.stderr
        self.min_interval = float(min_interval)
        self.warmup = int(warmup)
        self.done = 0
        self.t0 = time.perf_counter()
        self._last_print = 0.0
        self._warmed = False
        if self.total:
            self._emit(f"[{self.label}] starting — {self.total} items")
        else:
            self._emit(f"[{self.label}] starting — total unknown")

    def _emit(self, msg):
        self.stream.write(msg + "\n")
        self.stream.flush()

    def update(self, item_label="", n=1):
        self.done += n
        now = time.perf_counter()
        elapsed = now - self.t0
        first_estimate = (not self._warmed and self.total
                          and self.done >= self.warmup)
        due = (now - self._last_print) >= self.min_interval
        last = self.total and self.done >= self.total
        if not (first_estimate or due or last):
            return
        self._last_print = now
        per_item = elapsed / self.done if self.done else 0.0
        if self.total:
            remaining = max(self.total - self.done, 0)
            eta_s = remaining * per_item
            finish_at = datetime.now() + timedelta(seconds=eta_s)
            pct = 100.0 * self.done / self.total
            tag = "  (provisional)" if first_estimate else ""
            self._emit(
                f"[{self.label}] {self.done}/{self.total} ({pct:4.1f}%)  "
                f"elapsed {_fmt_dur(elapsed)}  avg {per_item:.1f}s/item  "
                f"ETA {_fmt_dur(eta_s)} (~{finish_at:%H:%M})"
                f"{('  ' + item_label) if item_label else ''}{tag}"
            )
            self._warmed = True
        else:
            self._emit(
                f"[{self.label}] {self.done} done  elapsed {_fmt_dur(elapsed)}  "
                f"avg {per_item:.1f}s/item"
                f"{('  ' + item_label) if item_label else ''}"
            )

    def finish(self):
        elapsed = time.perf_counter() - self.t0
        self._emit(f"[{self.label}] complete — {self.done} items in {_fmt_dur(elapsed)}")
