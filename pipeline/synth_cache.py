"""Per-line synthesis cache under `abundances_derive._synth_flux_at_abund` (RYA-1232).

WHY. A synthesis-route fit calls the synthesis ~20 times per line, and the orchestrator
repeats that fit for every holding and every RYA-587 budget leg -- 9 passes x 13 holdings of
the same physics. What a line looks like at a trial abundance depends on the line, the model
and the broadening, NOT on which spectrum it is compared with. Measured on Asplund's 34 VIS
Fe I lines (HARPS, ts-lte): a warm pass (a continuum leg, or a second holding) needs 16 new
syntheses instead of ~970 and runs in ~55 s instead of ~600 s.

HOW, without changing the fit.
  * Each call is answered from abundance-grid NODES (STEP_A dex apart) synthesised on a fixed
    fine wavelength grid (FINE_NM) covering the caller's window, by 4-node cubic Lagrange
    interpolation in A, then linear interpolation onto the caller's pixels. The fitter, its
    bracket and its RYA-847 constraint metrics are untouched.
  * Nodes are synthesised with R = 0 -- rotation and macroturbulence applied, the
    instrumental Gaussian NOT -- and the caller's R is applied afterwards with iSpec's own
    `convolve_spectrum`, exactly the step `apply_post_fundamental_effects` runs LAST. So one
    node serves every instrument resolution: HARPS, Kitt Peak and IAG share it.
  * The key is a digest of EVERY other argument (atmosphere, model parameters, line list,
    abundances, NLTE departures, molecules, ...) plus the window: a changed input is a new
    key, never a stale answer. tmp_dir is excluded (a scratch path).

Synthesising on the fine grid instead of on the observed pixels is MORE accurate where the
pixels under-sample the line: on HARPS chi2r fell 3-8x on several lines and the per-line
abundances agree to a median 0.003 dex (max 0.05 on an under-resolved line).

Cache dir: $CODEX_SYNTH_CACHE_DIR, default data/cache/synth (gitignored). Disable with
CODEX_SYNTH_CACHE=0. Writes are atomic (temp + rename): worker processes share it safely.
"""
from __future__ import annotations

import hashlib
import inspect
import os
import pickle
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STEP_A = 0.05            # dex between abundance nodes
FINE_NM = 0.0002         # 0.002 A synthesis grid
PAD_NM = 0.02            # grid margin beyond the caller's window
STATS = {"calls": 0, "synth": 0}
#: per-process synthesis scratch dirs seen by the cache; removed at exit (RYA-1232: 2,039 of
#: them, 50 GB, filled the Mac's disk mid-run)
_SCRATCH: set = set()
_INSTALLED = False


def cache_dir() -> Path:
    d = Path(os.environ.get("CODEX_SYNTH_CACHE_DIR") or ROOT / "data" / "cache" / "synth")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _digest(x, h=None):
    top = h is None
    h = h or hashlib.sha1()
    if isinstance(x, np.ndarray):
        if x.dtype.names:                     # structured (line list, atmosphere records)
            h.update(repr(x.dtype.descr).encode())
        else:
            h.update(str(x.dtype).encode())
        h.update(str(x.shape).encode())
        h.update(np.ascontiguousarray(x).tobytes())
    elif isinstance(x, (list, tuple)):
        h.update(b"[")
        for e in x:
            _digest(e, h)
        h.update(b"]")
    elif isinstance(x, dict):
        h.update(b"{")
        for k in sorted(x, key=str):
            h.update(str(k).encode())
            _digest(x[k], h)
        h.update(b"}")
    else:
        h.update(repr(x).encode())
    return h.hexdigest() if top else h


def install() -> bool:
    """Wrap `abundances_derive._synth_flux_at_abund` once per process. Returns True when on."""
    global _INSTALLED
    if _INSTALLED:
        return True
    # the live MPIA query gets its own disk cache + retry (pipeline/mpia_cache.py) -- installed
    # here because this is the one hook the fitter calls at start-up
    try:
        from pipeline import mpia_cache as _mc
        _mc.install()
    except Exception as exc:                                     # noqa: BLE001
        print(f"WARNING: MPIA cache not installed: {type(exc).__name__}: {exc}")
    if os.environ.get("CODEX_SYNTH_CACHE", "1") == "0":
        return False
    import pipeline.abundances_derive as AD
    orig = AD._synth_flux_at_abund
    sig = inspect.signature(orig)
    from ispec.spectrum import convolve_spectrum, create_spectrum_structure

    def cached(*a, **k):
        STATS["calls"] += 1
        b = sig.bind(*a, **k)
        b.apply_defaults()
        kw = dict(b.arguments)
        w = np.asarray(kw.pop("waveobs_nm"), float)
        trial = float(kw.pop("trial_A"))
        R = kw.pop("R")
        tmp = kw.pop("tmp_dir", None)
        if tmp:
            _SCRATCH.add(str(tmp))
        # NLTE departures (Gerber decks): the tuple carries the trial abundance as a STAMP
        # that bsyn checks against the synthesis abundance (gerber_nlte.coefficients). For a
        # deck with ONE A(X) node (Fe) the coefficients are node-fixed, so the stamp is
        # restamped per grid node and the key ignores it. A deck WITH an abundance axis has
        # coefficients that depend on A: the cache steps aside and synthesises directly.
        nd = kw.get("nlte_departures")
        if nd:
            from pipeline import gerber_nlte as _gn
            if any(_gn.has_abundance_axis(el) for el in nd):
                return orig(*a, **k)
            kw["nlte_departures"] = {el: tuple(t[:5]) + (None,) + tuple(t[6:]) for el, t in nd.items()}
        lo = np.floor((w.min() - PAD_NM) / 0.01) * 0.01
        hi = np.ceil((w.max() + PAD_NM) / 0.01) * 0.01
        fine = np.arange(lo, hi + FINE_NM / 2, FINE_NM)
        key = _digest([round(float(lo), 4), round(float(hi), 4), FINE_NM, kw])
        d = cache_dir()

        def node(a_):
            f = d / f"{key}_{a_:+.3f}.pkl"
            if f.exists():
                with f.open("rb") as fh:
                    return pickle.load(fh)
            extra = {"tmp_dir": tmp} if tmp is not None else {}
            kw_ = dict(kw)
            if kw_.get("nlte_departures"):
                kw_["nlte_departures"] = {el: tuple(t[:5]) + (float(a_),) + tuple(t[6:])
                                          for el, t in kw_["nlte_departures"].items()}
            flux = orig(fine, trial_A=float(a_), R=0, **kw_, **extra)
            STATS["synth"] += 1
            fd, tp = tempfile.mkstemp(dir=d, suffix=".tmp")
            with os.fdopen(fd, "wb") as fh:
                pickle.dump(flux, fh)
            os.replace(tp, f)
            return flux

        i = int(np.floor(trial / STEP_A))
        nodes = [round((i + j) * STEP_A, 3) for j in (-1, 0, 1, 2)]
        fl = [node(a_) for a_ in nodes]
        if any(f is None for f in fl):
            return None
        out = np.zeros(fine.size)
        for j, (aj, fj) in enumerate(zip(nodes, fl)):
            lj = 1.0
            for m, am in enumerate(nodes):
                if m != j:
                    lj *= (trial - am) / (aj - am)
            out += lj * np.asarray(fj, float)
        if R is not None and R > 0:
            # iSpec's own last broadening step (synth/effects.apply_post_fundamental_effects)
            out = convolve_spectrum(create_spectrum_structure(fine, out), R,
                                    from_resolution=None, frame=None)["flux"]
        return np.interp(w, fine, out)

    AD._synth_flux_at_abund = cached
    import atexit
    import shutil

    def _bye():
        for t in _SCRATCH:
            shutil.rmtree(t, ignore_errors=True)
        print(f"synth_cache: calls={STATS['calls']} new_syntheses={STATS['synth']} "
              f"dir={cache_dir()}", flush=True)
    atexit.register(_bye)
    _INSTALLED = True
    return True
