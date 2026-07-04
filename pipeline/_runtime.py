"""
pipeline/_runtime.py — centralized runtime init (RYA-514). Import FIRST, before numpy.

One platform-agnostic fix for the RYA-506 class of failures. Import this at the top of
every pipeline entry point (and any module that touches iSpec or multiprocessing) so the
process is configured identically no matter which entry point launched it.

Two coupled settings, both required:

1. **Force the multiprocessing start method to `fork`.** iSpec runs its SPECTRUM synthesis
   in a child process that imports the C `synthesizer` extension and relies on FORK-inherited
   C-extension state (the pipeline was developed on Linux, where fork was the default). Under
   a re-importing start method the child dies WITHOUT enqueueing its result and iSpec
   silently returns a zero-initialised array (RYA-506). BOTH platform defaults break it:
   macOS Python 3.8+ = `spawn`; Linux Python 3.14 = `forkserver`. Forcing `fork` overrides
   both — no per-OS branching. (The Codex cannot make iSpec run in-process without patching
   iSpec, so forcing fork is the fix.)

2. **Pin BLAS/LAPACK to a single thread** (OMP/OPENBLAS/MKL/VECLIB = 1). Forced fork + a
   THREADED BLAS is a classic silent-corruption hazard (a thread pool forked mid-operation
   leaves the child in an undefined state) — and the RYA-506 all-zero guard does NOT catch
   that. Single-thread BLAS removes the hazard and also makes results bit-deterministic
   across machines (the RYA-511 cross-system floor). Set as early as possible so the BLAS
   library reads them before it spins up its pool.

Idempotent + guarded: safe to import many times; never raises.
"""
import os

# BLAS/LAPACK single-thread pins — set before numpy/scipy import their backends.
_BLAS_ENV = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _v in _BLAS_ENV:
    os.environ[_v] = '1'

# Force fork — overrides macOS spawn and Linux-3.14 forkserver alike.
import multiprocessing as _mp  # noqa: E402
try:
    if _mp.get_start_method(allow_none=True) != 'fork':
        _mp.set_start_method('fork', force=True)
except (RuntimeError, ValueError):
    # A start method was already fixed by an earlier import, or the platform has no fork.
    # force=True normally prevents the RuntimeError; if we still land here the guard in
    # _fe2_theoretical_ew (RYA-506) is the loud backstop.
    pass

START_METHOD = _mp.get_start_method(allow_none=True)


def assert_fork_and_pins() -> dict:
    """Smoke helper: confirm the start method is fork and the BLAS pins are in-process."""
    sm = _mp.get_start_method(allow_none=True)
    env = {v: os.environ.get(v) for v in _BLAS_ENV}
    assert sm == 'fork', f"start method is {sm!r}, expected 'fork' (RYA-514 not applied)"
    assert all(val == '1' for val in env.values()), f"BLAS pins not all 1: {env}"
    return {'start_method': sm, 'blas_env': env}
