"""Read native statistical-equilibrium update flags without replacing the solver."""
import json
import time
import numpy as np


class StatEquilTrace:
    """Forward a Context and record returned native flags, including invalid retries.

    The original solver flushes one history row per completed outer iteration.
    A streaming reader therefore joins each native call to its actual outer row;
    an invalid-update retry shares the row and receives a new attempt index.
    No native update, population, rate, convergence control or return value changes.
    """
    def __init__(self, context, path, history):
        self._context = context
        self._stream = path.open('x')
        self._history_path = history
        self._reader = None
        self._completed = 0
        self._calls = 0
        self._outer = 0
        self._attempt = 0
        self.overhead_seconds = 0.

    def __getattr__(self, name):
        return getattr(self._context, name)

    def stat_equil(self, *args, **kwargs):
        before = time.perf_counter()
        if self._reader is None: self._reader = self._history_path.open('r')
        self._completed += sum(1 for _ in self._reader)
        outer = self._completed + 1
        self._attempt = self._attempt + 1 if outer == self._outer else 1
        self._outer = outer
        self._calls += 1
        self.overhead_seconds += time.perf_counter()-before
        update = self._context.stat_equil(*args, **kwargs)
        before = time.perf_counter()
        flags = update.ngAccelerated
        active = self._context.activeAtoms
        if (not update.updatedPops or len(flags) != len(active)
                or not all(isinstance(v, (bool, np.bool_)) for v in flags)):
            raise ValueError('Native update lacks exact per-active-species Ng flags')
        row = dict(schema='codex.bride_stat_equil_update/1', call_index=self._calls,
            outer_iteration=outer, attempt_in_outer=self._attempt,
            active_species=[str(a.atomicModel.element) for a in active],
            ngAccelerated=[bool(v) for v in flags], updatedPops=bool(update.updatedPops),
            populations_finite=[bool(np.isfinite(a.n).all()) for a in active],
            populations_positive=[bool(np.all(a.n > 0)) for a in active])
        self._stream.write(json.dumps(row, allow_nan=False)+'\n'); self._stream.flush()
        self.overhead_seconds += time.perf_counter()-before
        return update

    def __enter__(self): return self

    def __exit__(self, *args):
        self._stream.close()
        if self._reader is not None: self._reader.close()
