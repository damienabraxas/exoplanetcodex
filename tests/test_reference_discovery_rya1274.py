"""Reference discovery must be portable and fail without fabricating holdings."""
from pathlib import Path

import pytest

from config.constants import codex_path, require_codex_path
from scripts.generate_sources_page import audit_library


def test_existing_library_override_is_resolved_without_copying(tmp_path, monkeypatch):
    library = tmp_path / "held literature"
    library.mkdir()
    paper = library / "paper.txt"
    paper.write_text("Synthetic held evidence")
    monkeypatch.setenv("CODEX_REFERENCE_LIBRARY", str(library))
    assert codex_path("refs.library") == library
    assert require_codex_path("refs.library") == library
    assert list(library.iterdir()) == [paper]
    assert paper.read_text() == "Synthetic held evidence"


def test_missing_library_override_fails_without_creating_a_fallback(tmp_path, monkeypatch):
    missing = tmp_path / "missing holding"
    monkeypatch.setenv("CODEX_REFERENCE_LIBRARY", str(missing))
    with pytest.raises(FileNotFoundError, match="refs.library"):
        require_codex_path("refs.library", context="RYA-1274 synthetic probe")
    assert not missing.exists()


def test_library_audit_rejects_missing_resource_and_unindexed_evidence(tmp_path):
    with pytest.raises(SystemExit, match="reference library not found"):
        audit_library([], tmp_path / "missing")
    (tmp_path / "unindexed.txt").write_text("Synthetic evidence without a citation")
    assert audit_library([], tmp_path) == 1
