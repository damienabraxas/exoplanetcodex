"""
RYA-1096 — `copied_to` is a REPO-RELATIVE path, and an absolute one is a machine fact.

🔴 THE DEFECT, MEASURED. `publish_product` builds its sources with `Path(x).resolve()` and
stored that resolved path in `provenance.copied_to`. Publishing from inside the repo
therefore recorded e.g. `/Users/<me>/codex/rya1089/data/results/band_products/...`, which is
inside the repository ON THE MACHINE THAT PUBLISHED and nowhere else.

`feed_repo_reconciliation` reads `copied_to` to check the artifact is COMMITTED. On any
other checkout the absolute path is not under the repo root, so it reports
COPIED_TO_OUTSIDE_REPO. Eight freshly published products passed every check on the Mac and
turned that guard RED on Sirius: three tests that pass locally and fail on the runner, which
is the same shape as a test that silently skips.

⚠️ THIS TEST EXISTS BECAUSE THE MAC COULD NOT SEE IT. Nothing local was wrong -- the path
really was inside the local repo. Only the set-diff against Sirius exposed it, and a guard
that only one machine can evaluate is not a guard.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / "data" / "products"


def _feeds():
    return sorted(STORE.glob("*/[A-Z]*.json"))


@pytest.mark.parametrize("feed", _feeds(), ids=lambda p: p.name)
def test_every_copied_to_is_repo_relative_and_present(feed):
    doc = json.loads(feed.read_text())
    absolute, missing = [], []
    for pool in ("products", "quarantine", "superseded", "archive"):
        for p in doc.get(pool) or []:
            ct = (p.get("provenance") or {}).get("copied_to")
            if not ct:
                continue
            if Path(ct).is_absolute():
                absolute.append((pool, ct))
            elif not (ROOT / ct).exists():
                missing.append((pool, ct))
    assert not absolute, (
        f"{len(absolute)} record(s) store an ABSOLUTE copied_to. It resolves only on the "
        f"machine that published and makes feed_repo_reconciliation report "
        f"COPIED_TO_OUTSIDE_REPO everywhere else:\n  "
        + "\n  ".join(f"{pool}: {ct}" for pool, ct in absolute[:6]))
    assert not missing, (
        f"{len(missing)} record(s) name a copied_to that is not in this checkout:\n  "
        + "\n  ".join(f"{pool}: {ct}" for pool, ct in missing[:6]))


def test_publish_product_REFUSES_a_source_outside_the_repo():
    """POSITIVE CONTROL for the other half. `copied_to` asserts the artifact is committed
    HERE, so publishing from a scratch directory must refuse rather than record a path no
    other checkout can resolve -- which is exactly what I did, from a scratchpad, before
    the guard caught it."""
    import ast
    src = (ROOT / "scripts" / "publish_product.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "_repo_relative")
    body = ast.get_source_segment(src, fn)
    assert "relative_to(ROOT)" in body
    assert "SystemExit" in body, "an outside-repo source must REFUSE, not fall back"

    prov = next(n for n in ast.walk(ast.parse(src))
                if isinstance(n, ast.Call)
                and getattr(n.func, "id", "") == "_repo_relative")
    assert prov is not None, "the provenance builder must call _repo_relative"
