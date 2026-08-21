"""RYA-939: the telluric leg must not depend on WHICH MACHINE it runs on.

Three defects, all of the same shape — a property of the environment being
treated as a property of the code:

1. the esorex path was a Mac literal, so the correction leg was Mac-only while
   the project rule is that compute happens on Sirius;
2. product FILENAMES came from ``~/.esorex/esorex.rc``, which the two machines
   disagree about, so a wholly successful run reported ``failed (rc=0)``;
3. ``scripts/rya929_full_sweep.py`` still read the Kurucz 2005 product as air
   after RYA-938 measured it to be vacuum — queued reruns would have inherited
   a line table displaced by ~200 sampled pixels.

The tests pin the invariants, not the paths: any machine, resolved not assumed,
and a failure that names what is actually wrong.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = _load("rya931_molecfit_model")


def test_no_machine_specific_esorex_path_is_baked_in():
    """The one literal that made this leg Mac-only must not come back."""
    source = (ROOT / "scripts" / "rya931_molecfit_model.py").read_text()
    assert 'ESOREX = "/opt/homebrew' not in source
    assert "resolve_esorex" in source


def test_esorex_env_override_wins(tmp_path, monkeypatch):
    fake = tmp_path / "esorex"
    fake.write_text("#!/bin/sh\n")
    monkeypatch.setenv("ESOREX", str(fake))
    assert runner.resolve_esorex() == str(fake)


def test_esorex_env_override_must_exist(tmp_path, monkeypatch):
    """A wrong $ESOREX must fail loudly, not fall through to another machine's."""
    monkeypatch.setenv("ESOREX", str(tmp_path / "nope"))
    with pytest.raises(SystemExit) as excinfo:
        runner.resolve_esorex()
    assert "ESOREX" in str(excinfo.value)


def test_absent_esorex_names_everything_it_tried(tmp_path, monkeypatch):
    """A missing engine must never read as an uncorrected product."""
    monkeypatch.delenv("ESOREX", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    # Neutralise BOTH install prefixes, not just the registered one: the Mac's
    # Homebrew path is a real second candidate, and on a developer Mac it exists.
    monkeypatch.setattr(runner, "_esorex_candidates", tuple)
    with pytest.raises(SystemExit) as excinfo:
        runner.resolve_esorex()
    message = str(excinfo.value)
    assert "$ESOREX" in message and "PATH" in message
    assert "not a correction result" in message


def test_the_eso_prefix_is_registered_not_hard_coded():
    """Sirius's install prefix is declared in ONE place (RYA-810 discipline)."""
    import sys
    sys.path.insert(0, str(ROOT))
    from config.constants import codex_root
    assert codex_root("eso_pipelines").name == "molecfit"
    os.environ["CODEX_ESO_ROOT"] = "/tmp/elsewhere"
    try:
        assert str(codex_root("eso_pipelines")) == "/tmp/elsewhere"
    finally:
        del os.environ["CODEX_ESO_ROOT"]


def test_product_naming_is_pinned_by_the_invocation():
    """Not by whoever's home directory the recipe happened to run in."""
    source = (ROOT / "scripts" / "rya931_molecfit_model.py").read_text()
    assert "--suppress-prefix=TRUE" in source


def test_a_successful_recipe_with_no_product_is_not_reported_as_a_failure():
    """`failed (rc=0)` is a self-contradiction that cost real debugging time."""
    source = (ROOT / "scripts" / "rya931_molecfit_model.py").read_text()
    assert "SUCCEEDED (rc=0) but" in source, (
        "the missing-product branch must name what is actually wrong")
    assert "Produced:" in source, "it must list what the recipe did write"


def test_kurucz_2005_is_read_as_vacuum_everywhere():
    """RYA-938 measured the medium; every reader must honour it.

    RYA-933/934 are queued against this script. A reader that silently treats
    vacuum as air produces a line table that looks entirely plausible and is
    displaced by about 200 sampled pixels.
    """
    source = (ROOT / "scripts" / "rya929_full_sweep.py").read_text()
    assert "vacuum=True" in source, "the Kurucz 2005 call site must declare its medium"
    kurucz_call = source[source.index("kurucz_2005"):]
    assert "vacuum=True" in kurucz_call[:400]


def test_the_normalised_solar_products_are_tracked_with_generators():
    """A landed artifact without a committed generator is the RYA-559 hole."""
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import yaml
    manifest = yaml.safe_load((ROOT / "data" / "results" / "GENERATORS.yaml").read_text())
    by_artifact = {e["artifact"]: e for e in manifest["artifacts"]}
    for name in ("data/processed/solar_normalized.csv",
                 "data/processed/solar_normalized_tellcorr.csv"):
        entry = by_artifact[name]
        assert entry.get("status", "COMMITTED") == "COMMITTED"
        assert (ROOT / entry["generator"]).is_file(), entry["generator"]
