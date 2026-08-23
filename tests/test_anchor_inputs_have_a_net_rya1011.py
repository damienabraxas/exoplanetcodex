"""A load-bearing product must be recoverable by something other than luck (RYA-1011).

On 2026-08-23 the 163-line graded anchor's per-line inputs were silently wiped by a
stem-sharing flag (`--local-renorm` on pre-RYA-1006 code). They were recovered ONLY
because Sirius happened to hold copies whose mtimes predated the overwrite. RYA-1006
closed that particular collision; this closes the shape of it.

The exposure was structural, not incidental. Every other per-line artifact in the
project lives in a tracked ticket subdirectory (`data/results/rya847/`, `rya880/`,
`rya925/` ...). The anchor's inputs lived only in `data/results/band_products/` -- the
one results directory that is gitignored -- so they had no git safety net at all, while
13 less important per-line files did.

So the invariant is not "band_products must be tracked" (it is ignored on purpose, to
stop bulk accidental tracking). It is narrower and it is about consequence: **anything
`anchor_pools` reads by name must be recoverable** -- committed, and able to say how it
was made. 96 KB buys that for all three.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def _anchor_parts() -> list[str]:
    """Every filename anchor_pools load-bears on, read from the module itself.

    Read from ANCHORS rather than hardcoded here: a new anchor must inherit the
    protection automatically, or this test is just a second place to forget.
    """
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline import anchor_pools
    parts: list[str] = []
    for pool in anchor_pools.ANCHORS.values():
        parts.extend(pool.parts)
    return parts


def _tracked() -> set[str]:
    out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True)
    return {Path(p).name for p in out.stdout.split("\n") if p}


def _registered() -> str:
    return (ROOT / "data" / "results" / "GENERATORS.yaml").read_text(encoding="utf-8")


def test_there_is_at_least_one_anchor_to_protect():
    """Guard the guard: an empty ANCHORS makes every assertion below vacuously true."""
    assert _anchor_parts(), "anchor_pools.ANCHORS is empty — this test would pass on nothing"


@pytest.mark.parametrize("name", _anchor_parts())
def test_every_anchor_input_is_committed(name):
    """🔴 A gitignored load-bearing file has NO safety net — losing it is unrecoverable
    on the machine that lost it."""
    assert name in _tracked(), (
        f"{name} is load-bearing (anchor_pools reads it by name) but is NOT git-tracked. "
        f"band_products/ is gitignored, so this file can be silently overwritten with no "
        f"way back. Commit it with `git add -f` — the anchor inputs are ~40 KB each."
    )


@pytest.mark.parametrize("name", _anchor_parts())
def test_every_anchor_input_can_say_how_it_was_made(name):
    """Committed is not enough: the RYA-1011 answer is regenerability, so the file must
    also carry the invocation that rebuilds it."""
    assert name in _registered(), (
        f"{name} is load-bearing but has no data/results/GENERATORS.yaml entry, so it "
        f"cannot state how to rebuild itself. Being committed protects the bytes; the "
        f"generator record is what makes a LOSS recoverable rather than terminal."
    )
