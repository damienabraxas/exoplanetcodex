"""
pipeline/telluric/esorex_runtime.py
===================================
RYA-963 — the ONE place that knows where the ESO pipeline lives and how to invoke it.

RYA-939 found the esorex path baked into `scripts/rya931_molecfit_model.py` and fixed it
there. It was baked into `pipeline/crires_telluric.py` too (`_ESOREX =
"/opt/homebrew/bin/esorex"`), so the CRIRES+ telluric leg stayed Mac-only while the HARPS
leg had already been freed — a refuted rule that died at one site only. This module is the
shared implementation both call; there is no second copy to drift.

What it owns:

* `resolve_esorex()` — $ESOREX → PATH → the `eso_pipelines` path-register root → loud
  failure naming everything tried. A missing engine must never read as an uncorrected
  product.
* `esorex_env()` — the ESO source-kit prefix needs its own `lib/` on `LD_LIBRARY_PATH`
  and its `bin/` on `PATH`; the Homebrew tap needs neither, and adding them is harmless.
* `SUPPRESS_PREFIX` — **always pass `--suppress-prefix=TRUE` explicitly** (RYA-939).
  esorex otherwise takes it from `~/.esorex/esorex.rc`, and the two machines disagree:
  the Mac has no rc (compiled default TRUE), a fresh source-kit install writes FALSE.
  Under FALSE the recipe SUCCEEDS and writes every product as `out_0000.fits …`, so
  `BEST_FIT_MODEL.fits` is simply absent and the run reads as the self-contradicting
  `failed (rc=0)`.
* `gdas_dirs()` — the installed telluriccorr GDAS profile directories, including the
  registered ESO prefix. `pipeline.telluric.gdas_fetch` globbed only Homebrew and
  `/usr/share/esopipes`, so on Sirius the per-night profile was unreachable and the
  loud-fail would have fired on a profile that is sitting on the disk.
"""
from __future__ import annotations

import os
from pathlib import Path
from shutil import which

#: Always pass this to esorex. See the module docstring (RYA-939).
SUPPRESS_PREFIX = "--suppress-prefix=TRUE"


class EsorexNotAvailable(RuntimeError):
    """esorex/molecfit is not installed. Raising this is the correct outcome — a
    telluric correction that did not run must never be reported as one that did."""


def eso_root() -> "Path | None":
    """The registered ESO pipeline install prefix, or None if the register is
    unavailable (the resolver still has $ESOREX and PATH to try)."""
    try:
        from config.constants import codex_root
        return Path(codex_root("eso_pipelines"))
    except Exception:
        return None


def _candidates() -> tuple:
    """Known install prefixes, most specific first: the registered root (Sirius's ESO
    source kit) then the Mac's Homebrew tap."""
    out = []
    root = eso_root()
    if root is not None:
        out.append(root / "bin" / "esorex")
    out.append(Path("/opt/homebrew/bin/esorex"))
    out.append(Path("/usr/local/bin/esorex"))
    return tuple(out)


def resolve_esorex(required: bool = True) -> "str | None":
    """Return the esorex executable path. Order: `$ESOREX`, `PATH`, the registered
    `eso_pipelines` root, the Homebrew tap. With `required=False` return None instead of
    raising (the availability probe)."""
    explicit = os.environ.get("ESOREX")
    if explicit:
        if not Path(explicit).is_file():
            raise EsorexNotAvailable(f"$ESOREX={explicit!r} is not a file")
        return explicit
    found = which("esorex")
    if found:
        return found
    for cand in _candidates():
        if cand.is_file():
            return str(cand)
    if not required:
        return None
    raise EsorexNotAvailable(
        "esorex not found. Tried $ESOREX, PATH, and "
        + ", ".join(str(c) for c in _candidates())
        + ". Install the ESO molecfit kit or set $ESOREX. Failing here rather than "
          "reporting an uncorrected product — a missing engine is not a correction.")


def esorex_available() -> bool:
    return resolve_esorex(required=False) is not None


def esorex_env(esorex: str, base: "dict | None" = None) -> dict:
    """Environment for an esorex subprocess: the install's `bin/` on PATH and, for an
    ESO source-kit prefix, its `lib/` on LD_LIBRARY_PATH."""
    env = dict(base if base is not None else os.environ)
    bindir = Path(esorex).parent
    env["PATH"] = f"{bindir}:{env.get('PATH', '')}"
    libdir = bindir.parent / "lib"
    if libdir.is_dir():
        env["LD_LIBRARY_PATH"] = f"{libdir}:{env.get('LD_LIBRARY_PATH', '')}"
    return env


def gdas_dirs() -> tuple:
    """Candidate telluriccorr GDAS profile directories (glob patterns), registered root
    first. Consumed by `pipeline.telluric.gdas_fetch`."""
    pats = []
    root = eso_root()
    if root is not None:
        pats.append(str(root / "share" / "molecfit" / "data" / "profiles" / "gdas"))
        pats.append(str(root / "share" / "telluriccorr" / "profiles" / "gdas"))
    pats += [
        "/opt/homebrew/Cellar/telluriccorr/*/share/molecfit/data/profiles/gdas",
        "/usr/local/Cellar/telluriccorr/*/share/molecfit/data/profiles/gdas",
        "/usr/share/esopipes/datastatic/telluriccorr*/profiles/gdas",
    ]
    return tuple(pats)
