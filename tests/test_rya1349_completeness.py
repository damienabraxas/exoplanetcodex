"""RYA-1349: the orchestrator answers "did we do ENGINE-A? CRIRES? every band?" itself.

Ryan 2026-11-07, after Ca shipped without ENGINE-A on some cells: the checklist must prove an
element complete. Three holes it closes in the RYA-1233 coverage table:
  - a missing ENGINE-A product was HELD whenever a 1D-LTE product existed, without asking
    whether the grid holds the pool's lines (it did, on Ca I near-UV Reference Grade);
  - a published ENGINE-A product that corrected no line counted as done;
  - STALE moved on any edit in a band-wide window -- another species' gf source tag, or the
    code -- so every live Ti product read STALE and nobody could use the table.
"""
import json
import subprocess
import sys

from pipeline import coverage_report as cr
from pipeline import run_matrix as rm


def _p(holding, A, sigma, **kw):
    return dict(ion="II", band="NIR", selector="GRADED", treatment="1D-LTE",
                holding=holding, A=A, sigma_reported=sigma, **kw)


def test_disagreement_flags_the_outlier_not_the_holdings_that_agree():
    """Ti II NIR before RYA-1315's quarantine: KP-Kurucz 5.30 +/- 0.38 against IAG 4.81 and
    KP-molecfit 4.84 on the same two lines. Judged by its OWN sigma it passes (and the two
    agreeing holdings get blamed); judged by the others' precision it is the one flagged."""
    flags = cr.disagreements([_p("solar_iag", 4.808, 0.092),
                              _p("solar_kpno_molecfit_corrected", 4.841, 0.056),
                              _p("solar_kpno_kurucz2005_corrected", 5.300, 0.379)])
    assert [f["holding"] for f in flags] == ["solar_kpno_kurucz2005_corrected"]
    assert flags[0]["others_median"] == 4.825


def test_two_holdings_are_never_flagged():
    """With two holdings there is no majority to say which one is wrong."""
    assert cr.disagreements([_p("solar_iag", 4.8, 0.05), _p("solar_harps", 5.8, 0.05)]) == []


def test_ti_has_no_cited_3d_nlte_table_and_ca_does():
    """ENGINE-A-3DNLTE is HELD for Ti because no table exists in the literature (A21: owed),
    not because a 1D-LTE product does; Ca I 4512.27 is in A21 Table A.1."""
    n, why = cr.correction_reach("Ti", "I", "ENGINE-A-3DNLTE", [5866.45])
    assert n == 0 and "no cited 3D NLTE table for Ti I" in why
    n, _ = cr.correction_reach("Ca", "I", "ENGINE-A-3DNLTE", [4512.268])
    assert n == 1


def test_engine_a_reach_reads_the_registered_grid_and_its_ions():
    n, src = cr.correction_reach("Ti", "II", "ENGINE-A", [4501.27, 9999.0])
    assert src == "Ti_Mallinson2024_PySME.csv" and n == 1


def test_every_published_ti_engine_a_product_corrected_its_lines():
    feed = json.loads((rm.ROOT / "data/products/solar/Ti.json").read_text())
    ea = [p for p in feed["products"] if p["treatment"] == "ENGINE-A"]
    assert len(ea) == 55
    for p in ea:
        moved, n = cr.correction_applied(p)
        assert moved > 0 and moved == n, (p["band"], p["holding"], moved, n)


def test_every_live_product_of_the_orchestrator_era_has_its_lines_file():
    """lines_artifact follows the feed's own `copied_to`; the 10 Fe II products RYA-1232
    restored to superseded/ lost theirs and RYA-1349 put them back, sha-matched."""
    for el in ("Ca", "Mg", "Na", "Ti"):
        feed = json.loads((rm.ROOT / f"data/products/solar/{el}.json").read_text())
        missing = [p["provenance"]["copied_to"] for p in feed["products"]
                   if cr.lines_artifact(p) is None]
        assert missing == [], (el, missing[:3])
    feed = json.loads((rm.ROOT / "data/products/solar/Fe.json").read_text())
    fe2 = [p for p in feed["products"] if str(p.get("ion")) == "II"]
    assert fe2 and all(cr.lines_artifact(p) for p in fe2)


def test_code_drift_is_reported_not_owed():
    rec = {"inputs": [{"kind": "stage_script", "name": "scripts/derive_band_products.py",
                       "digest": "a"},
                      {"kind": "spectrum", "name": "x.fits", "digest": "s"}],
           "code_commit": "0" * 40}
    now = [{"kind": "stage_script", "name": "scripts/derive_band_products.py", "digest": "b"},
           {"kind": "spectrum", "name": "x.fits", "digest": "s"}]
    out, why = cr._staleness(rec, now, "Ti", (5000.0, 5100.0))
    assert out == "CODE_DRIFT" and "derive_band_products.py" in why
    now[1]["digest"] = "t"
    assert cr._staleness(rec, now, "Ti", (5000.0, 5100.0))[0] == "STALE"


def test_a_dirty_build_commit_is_used_only_when_its_file_proves_out():
    """A build from a dirty tree records `<sha>-dirty`. The committed file is read in its
    place only when the build's own scoped fingerprint of it equals the recorded digest."""
    rel = "data/registry/problem_children.csv"
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=rm.ROOT, capture_output=True,
                          text=True).stdout.strip()
    good = rm.registry_fingerprint(3770, 6920)
    assert cr._built_table("registry_scoped", rel, head + "-dirty", good, 3770.0, 6920.0) is not None
    assert cr._built_table("registry_scoped", rel, head + "-dirty", "f" * 64, 3770.0, 6920.0) is None


def test_a_dry_run_never_becomes_the_elements_report(tmp_path):
    doc = rm.run("solar", "Si", dry_run=True, interpreter=sys.executable,
                 ispec_dir="/nonexistent/ispec", echo=False, report_dir=tmp_path)
    assert doc["dry_run"] and not doc["scope"]
    assert (tmp_path / "solar_Si_latest_dryrun.json").exists()
    assert not (tmp_path / "solar_Si_latest.json").exists()


def test_ti_bookkeeping_is_committed():
    """Ti's runs recorded 355 builds in a Sirius checkout and none reached git, and its
    `_latest` was the 10-29 dry run: coverage read every live Ti product as unbuilt."""
    led = json.loads(rm.LEDGER.read_text())
    assert sum("|Ti|" in k for k in led) >= 330
    rep = json.loads((rm.REPORT_DIR / "solar_Ti_latest.json").read_text())
    assert not rep["dry_run"]


def test_deep_cells_without_a_product_are_held_by_policy():
    cov = cr.coverage("solar", "Ti")
    deep = [r for r in cov["rows"] if r["pool"] == "deep" and r["A"] is None]
    assert deep and all(r["outcome"] == "HELD" for r in deep)
    assert any(r["reason"] == cr.DEEP_POLICY for r in deep)


def test_ti_owes_no_cell_and_ca_owes_exactly_its_engine_a_cells():
    """Host-independent: NOT_RUN is decided from the feed, the grid and the run report, never
    from spectrum fingerprints (STALE needs the data host)."""
    ti = cr.coverage("solar", "Ti")
    assert [r for r in ti["rows"] if r["outcome"] == "NOT_RUN"] == []
    ca = cr.coverage("solar", "Ca")
    owed = [r for r in ca["rows"] if r["outcome"] == "NOT_RUN"]
    assert owed and all(r["treatment"] == "ENGINE-A" and "Ca_Amarsi2020_PySME.csv covers" in r["reason"]
                        for r in owed)
    assert {r["band"] for r in owed} == {"red-optical"}
    # near-UV Reference Grade: the pool is quarantined (RYA-1311), so its ENGINE-A is HELD
    nuv = [r for r in ca["rows"] if r["band"] == "near-UV" and r["pool"] == "reference"
           and r["treatment"] == "ENGINE-A" and r["ion"] == "I"]
    assert nuv and all(r["outcome"] == "HELD" and "quarantined" in r["reason"] for r in nuv)
