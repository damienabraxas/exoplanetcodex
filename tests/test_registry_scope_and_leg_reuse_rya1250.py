"""RYA-1250: culling one line re-runs only the cells that can see it.

The cull registry was hashed whole into every cell's identity, and budget legs were kept or
wiped for a whole element at once, so culling Mg I 24566.30 (K) meant re-running every Mg cell
and all 1548 legs. The registry is now scoped to the cell's window, and legs are kept per unit.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

from pipeline import element_publish as ep
from pipeline import run_matrix as rm

ROOT = Path(__file__).resolve().parents[1]
HEAD = ("species,lambda_or_scope,problem_class,required_treatment,observed_in,amplifies_with,"
        "severity,governing_tickets,status,population_source,notes\n")


def _reg(tmp_path, *rows):
    p = tmp_path / "problem_children.csv"
    p.write_text(HEAD + "".join(r + "\n" for r in rows))
    return p


def test_a_k_cull_moves_k_cells_only(tmp_path):
    base = ["Mg I,15748.988,ATOMIC_BLEND,exclude,Sun,all,medium,1250,active,curated,x",
            "Si I,gf pool,GF,flag,Sun,all,low,1,owed,curated,y"]
    a = _reg(tmp_path, *base)
    k0, v0 = rm.registry_fingerprint(19442, 24856, a), rm.registry_fingerprint(3770, 6920, a)
    a.write_text(a.read_text() + "Mg I,24566.296,ATOMIC_BLEND,exclude,Sun,all,medium,1250,active,curated,z\n")
    assert rm.registry_fingerprint(19442, 24856, a) != k0
    assert rm.registry_fingerprint(3770, 6920, a) == v0


def test_a_scope_row_without_a_wavelength_reaches_every_window(tmp_path):
    a = _reg(tmp_path, "Mg I,15748.988,ATOMIC_BLEND,exclude,Sun,all,medium,1250,active,curated,x")
    v0 = rm.registry_fingerprint(3770, 6920, a)
    a.write_text(a.read_text() + "Ti I,gf pool,GF,flag,Sun,all,low,1,owed,curated,y\n")
    assert rm.registry_fingerprint(3770, 6920, a) != v0


def _legs(root, n):
    for i in range(n):
        for leg in ("nominal", "xi_minus"):
            d = root / f"unit{i:02d}" / leg
            d.mkdir(parents=True)
            (d / "DONE").write_text(f"u{i}")


def test_legs_are_kept_per_unit_and_follow_their_unit(tmp_path):
    legs, prev = tmp_path / "legs", tmp_path / "units.prev.json"
    _legs(legs, 3)
    prev.write_text(json.dumps(["a", "b", "c"]))
    kept = ep.reuse_legs(legs, prev, ["a", "c2", "b"])
    assert kept == [0, 2]
    assert (legs / "unit00/nominal/DONE").read_text() == "u0"
    assert (legs / "unit02/nominal/DONE").read_text() == "u1"      # b moved from index 1 to 2
    assert not (legs / "unit01").exists()                            # c changed: re-priced
    assert json.loads(prev.read_text()) == ["a", "c2", "b"]


def test_no_per_unit_record_keeps_nothing(tmp_path):
    legs = tmp_path / "legs"
    _legs(legs, 2)
    assert ep.reuse_legs(legs, tmp_path / "units.prev.json", ["a", "b"]) == []
    assert not legs.exists()


def test_a_duplicated_signature_is_never_reused(tmp_path):
    legs, prev = tmp_path / "legs", tmp_path / "units.prev.json"
    _legs(legs, 2)
    prev.write_text(json.dumps(["a", "a"]))
    assert ep.reuse_legs(legs, prev, ["a", "a"]) == []


def test_migration_reproduces_the_scoped_identity():
    spec = importlib.util.spec_from_file_location(
        "rescope", ROOT / "scripts/rya1250_rescope_registry_fingerprint.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    whole = rm._file_fingerprint(rm.PROBLEM_CHILDREN)
    rows = [{"kind": "ledger_scoped", "name": "data/linelists/canonical_gf.csv[19442-24856]",
             "digest": "x"},
            {"kind": "ledger", "name": m.REG, "digest": whole}]
    entry = {"inputs": rows,
             "inputs_hash": hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest()}
    new, why = m.rescope(entry, whole)
    assert why == "" and new is not None
    reg = [r for r in entry["inputs"] if r["kind"] == "registry_scoped"]
    assert reg == [{"kind": "registry_scoped", "name": f"{m.REG}[19442-24856]",
                    "digest": rm.registry_fingerprint(19442, 24856)}]
    assert new == hashlib.sha256(json.dumps(entry["inputs"], sort_keys=True).encode()).hexdigest()
    # a cell built against an older registry stays stale
    stale = {"inputs": [dict(r, digest="old") if r["kind"] == "ledger" else r for r in rows]}
    stale["inputs_hash"] = hashlib.sha256(json.dumps(stale["inputs"], sort_keys=True).encode()).hexdigest()
    assert m.rescope(stale, whole)[0] is None


def test_mg_k_24566_is_culled_for_the_sun():
    from pipeline import problem_children as pc
    d = pc.disposition_for_line("Mg", "I", 24566.296, table=pc.line_dispositions())
    assert pc.aggregate_action(d) == "exclude"


def test_an_unpairable_product_is_held_not_fatal(tmp_path, monkeypatch):
    """Mg run 3 died in assembly on Mg II K Codex (no accepted line shared by core and
    nominal) after all 1548 legs had run; that one product is held, the rest assemble."""
    import sys
    from pipeline.paired_differential import UnpairableProducts
    spec = importlib.util.spec_from_file_location("budget", ROOT / "scripts/rya1230_cno_budget.py")
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    legs = tmp_path / "legs"
    (legs / "unit00" / "nominal").mkdir(parents=True)
    (legs / "unit00" / "nominal" / "MgII_1_2_x_lines.csv").write_text("a\n")
    (legs / "units.txt").write_text("python3 scripts/derive_band_products.py --element Mg\n")

    def boom(*a, **k):
        raise UnpairableProducts("no line is in-aggregate in BOTH products")
    monkeypatch.setattr(b, "build", boom)
    monkeypatch.setattr(b, "GF", ROOT / "data/registry/problem_children.csv")   # any small csv
    rep = tmp_path / "report.json"
    monkeypatch.setattr(sys, "argv", ["x", "--legs", str(legs), "--stage", str(tmp_path / "stage"),
                                      "--report", str(rep)])
    b.main()
    (o,) = json.loads(rep.read_text())["products"]
    assert o["skip"].startswith("UNPAIRABLE")


def _budget_module():
    spec = importlib.util.spec_from_file_location("budget", ROOT / "scripts/rya1230_cno_budget.py")
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    return b


def test_the_budget_gf_join_is_the_synthesis_join():
    """15765.839 (list) sits 0.014 A from its canonical row 15765.825; the synthesis resolves
    it (gf_resolver, 0.02 A) so the budget must too -- it HELD 21 Mg products."""
    import pandas as pd
    from pipeline import gf_resolver as gfr
    g = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    g = g[g.species == "Mg I"]
    for w, ep in ((15765.839, 5.9328), (15886.188, 5.9459), (9983.188, 5.9315), (10811.122, 5.9459)):
        m = g[((g.wavelength_air_A - w).abs() <= gfr.WTOL) & ((g.excitation_potential_eV - ep).abs() <= gfr.EPTOL)]
        assert len(m) >= 1, w
    src = (ROOT / "scripts/rya1230_cno_budget.py").read_text()
    assert "<= WTOL" in src and "< 0.01)" not in src


def test_a_strong_split_line_takes_its_measured_collapse_bound(tmp_path, monkeypatch):
    import pandas as pd
    b = _budget_module()
    f = tmp_path / "c.csv"
    pd.DataFrame({"physical_id": ["p1", "p2"], "dA_collapse_dex": [0.004, -0.012]}).to_csv(f, index=False)
    monkeypatch.setattr(b, "HFS_COLLAPSE", f)
    monkeypatch.setattr(b, "_hfs_thin_evidence", lambda rows: None)
    rows = pd.DataFrame({"physical_id": ["p1", "p2"], "wavelength_air_A": [1.0, 2.0]})
    assert b._hfs_collapse_bound(rows)["sigma_dex"] == 0.012
    assert b._hfs_collapse_bound(rows.assign(physical_id=["p1", "p3"])) is None   # unmeasured: HOLD
