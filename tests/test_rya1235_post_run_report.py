"""RYA-1235: zero-token report poster. Routing only through linear_dossiers.csv; ASCII
bodies; an unchanged report is never re-posted; no key, no post."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import post_run_report as prr  # noqa: E402


def test_every_canonical_target_routes_to_exactly_one_dossier():
    master = [e["symbol"] for e in
              json.loads((ROOT / "data/config/elements_master.json").read_text())["elements"]]
    r = prr.routes("solar")
    assert sorted(r) == sorted(master)
    assert len(set(r.values())) == len(r)


def test_no_issue_identifier_is_hardcoded_in_the_script():
    src = (ROOT / "scripts" / "post_run_report.py").read_text()
    import re
    assert not re.search(r"RYA-\d+", src.split('"""', 2)[2]), "an issue id is in the code"


def test_bodies_are_ascii_and_the_hash_ignores_commit_and_clock():
    doc = {"code_commit": "abc", "generated_at": "2026-10-03T00:00:00Z",
           "counts": {"DONE": 1}, "cells": [],
           "verdict": {"element_verdict": "PASS", "cells": []}, "problem_lines": []}
    t1, h1 = prr.body("solar", "Fe", doc)
    t2, h2 = prr.body("solar", "Fe", {**doc, "code_commit": "def",
                                      "generated_at": "2026-10-04T00:00:00Z"})
    assert t1.isascii() and h1 == h2
    _, h3 = prr.body("solar", "Fe", {**doc, "verdict": {"element_verdict": "REVIEW",
                                                       "cells": []}})
    assert h3 != h1


def test_ascii_only_transliterates():
    assert prr.ascii_only("Å ± σ") .isascii()


def test_missing_key_fails_loud_and_posts_nothing(monkeypatch, capsys):
    monkeypatch.delenv("LINEAR_API_KEY", raising=False)
    monkeypatch.setattr(prr, "post", lambda *a, **k: pytest.fail("posted without a key"))
    assert prr.main(["--star", "solar", "--element", "Si"]) == 2
    assert "LINEAR_API_KEY" in capsys.readouterr().err


def test_an_unchanged_report_is_skipped_not_reposted(monkeypatch, tmp_path):
    rep = tmp_path / "solar_Si_latest.json"
    rep.write_text(json.dumps({"counts": {"HELD": 3}, "cells": [], "problem_lines": [],
                               "verdict": {"element_verdict": "INCOMPLETE", "cells": []}}))
    monkeypatch.setattr(prr, "REPORT_DIR", tmp_path)
    monkeypatch.setattr(prr, "POSTED", tmp_path / "posted_ledger.json")
    monkeypatch.setenv("LINEAR_API_KEY", "x")
    posts = []
    monkeypatch.setattr(prr, "post", lambda issue, text, key: (posts.append(issue), (True, "ok"))[1])
    assert prr.main(["--star", "solar", "--element", "Si"]) == 0
    assert prr.main(["--star", "solar", "--element", "Si"]) == 0
    assert posts == ["RYA-725"], posts


def test_an_unrouted_element_is_posted_nowhere(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(prr, "DOSSIERS", tmp_path / "d.csv")
    (tmp_path / "d.csv").write_text("star,element,ion,issue_identifier\n")
    (tmp_path / "solar_Si_latest.json").write_text(json.dumps({"counts": {}, "cells": []}))
    monkeypatch.setattr(prr, "REPORT_DIR", tmp_path)
    monkeypatch.setattr(prr, "post", lambda *a, **k: pytest.fail("posted an unrouted element"))
    assert prr.main(["--star", "solar", "--element", "Si", "--dry-run"]) == 0
    assert "UNROUTED 1" in capsys.readouterr().out
