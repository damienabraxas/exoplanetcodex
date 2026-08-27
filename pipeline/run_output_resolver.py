"""RYA-1074 — resolve what a run ACTUALLY wrote. Never rebuild a path from intent.

🔴 THE DEFECT THIS EXISTS TO KILL. `derive_band_products` narrows a requested band to the
arm's real extent (RYA-1046 — correct science: a window half outside the data is not a
measurement) and then builds every output stem from the NARROWED values. So a run asked
for 4200-6910 A on HARPS writes `FeI_4200_6908_harps_...`.

On 2026-08-27 a verification step globbed the REQUESTED stem, `FeI_4200_6910_harps_*`,
matched files left by an earlier run, and reported "EXACT reproduction of published
1D-LTE and ENGINE-A". The check could not have failed: it compared the store against
files the run never wrote. A self-confirming check, reported as verification.

⚠️ AND THE OBVIOUS FIX IS WRONG. "Just glob more carefully" or "derive the trimmed stem
yourself" still reconstructs a path from intent — and the reconstruction is not even
well-defined: the trim logs `%.1f` (6909.0) while the stem uses `int()` (6908), so the
log and the filename disagree by rounding and NEITHER is the number. The only durable
answer is that the run records what it produced and verification reads that record.

THE RULE: given a requested band, this module resolves the effective stem via the run's
own `_runinfo.json`, or RAISES. It never falls back to the requested stem, because the
fallback is exactly the bug — a stale file sitting at the requested path makes the
fallback SUCCEED, silently, with the wrong answer.
"""
from __future__ import annotations

import json
from pathlib import Path


class StemResolutionError(AssertionError):
    """Asked to verify against a stem no run in scope recorded writing."""


def requested_stem(element: str, ion: str, lo: float, hi: float, instrument: str) -> str:
    """The stem the CALLER asked for. Not necessarily one that was ever written."""
    return f"{element}{ion}_{int(lo)}_{int(hi)}_{instrument}"


def load_runinfos(out_dir: Path) -> list[dict]:
    infos = []
    for p in sorted(Path(out_dir).glob("*_runinfo.json")):
        try:
            infos.append(json.loads(p.read_text()))
        except Exception as exc:                       # a corrupt record is not "absent"
            raise StemResolutionError(
                f"{p.name} is unreadable ({type(exc).__name__}). A run manifest that "
                f"cannot be parsed must not be skipped — skipping it would let "
                f"verification fall back to globbing, which is the RYA-1074 defect.")
    return infos


def resolve_stem(out_dir: Path, *, element: str, ion: str, lo: float, hi: float,
                 instrument: str, holding: str | None = None) -> str:
    """The stem a run ACTUALLY wrote for this request. Raises if no run recorded one.

    🔴 THERE IS NO FALLBACK TO THE REQUESTED STEM, BY DESIGN. If nothing recorded writing
    this request, the honest answer is "no run in scope produced it" — not "here is the
    path it would have had", which is how a stale file gets mistaken for a fresh result.
    """
    want = requested_stem(element, ion, lo, hi, instrument)
    infos = load_runinfos(out_dir)
    if not infos:
        raise StemResolutionError(
            f"no *_runinfo.json in {out_dir}: nothing recorded writing anything, so the "
            f"stem for {want!r} cannot be resolved. Re-run the deriver (it writes one), "
            f"and do NOT fall back to globbing {want}_* — a file sitting there may be "
            f"from an earlier run with a different effective band (RYA-1074).")
    cands = [i for i in infos
             if i.get("requested_stem") == want or i.get("stem") == want]
    if holding is not None:
        cands = [i for i in cands if i.get("holding") in (None, holding)] or cands
    if not cands:
        seen = sorted({i.get("requested_stem", "?") for i in infos})
        raise StemResolutionError(
            f"no run recorded writing {want!r} in {out_dir}. Runs present asked for: "
            f"{seen}. Refusing to glob {want}_* — matching a pre-existing file there is "
            f"exactly the vacuous check RYA-1074 exists to prevent.")
    newest = max(cands, key=lambda i: i.get("written_utc", ""))
    return newest["stem"]


def resolve_products(out_dir: Path, **kw) -> dict[str, Path]:
    """`treatment -> products.csv path`, from the run's own file list.

    Reads `files_written`, so a file that appeared in the directory AFTER the run — by any
    other route — cannot enter the comparison.
    """
    out_dir = Path(out_dir)
    stem = resolve_stem(out_dir, **kw)
    info = next(i for i in load_runinfos(out_dir) if i.get("stem") == stem)
    got = {}
    for name in info.get("files_written", []):
        if not name.endswith("_products.csv"):
            continue
        mid = name[len(stem):-len("_products.csv")].strip("_")
        got[mid or "__base__"] = out_dir / name
    if not got:
        raise StemResolutionError(
            f"run {stem!r} recorded no *_products.csv in files_written. It produced no "
            f"comparable product; say so rather than searching the directory for one.")
    return got
