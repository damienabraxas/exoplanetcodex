"""RYA-470 (Ryan 2026-11-18): the IR hyperfine free-continuum special case is wired into the
ORCH as a second cell BESIDE the standard one -- never instead of it. iSpec-free (Mac-runnable)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]


# ── the ORCH wiring (RYA-470: trigger = IR band AND a hyperfine line; BESIDE the standard) ──

def _canon(tmp_path):
    f = tmp_path / "canonical_gf.csv"
    f.write_text("species,wavelength_air_A,hfs_n_components\n"
                 "V I,15924.821,21\nV I,5737.031,9\nFe I,15964.865,1\n")
    return f


def test_the_trigger_needs_an_ir_band_and_a_hyperfine_line(tmp_path):
    from pipeline import free_continuum as fc
    c = _canon(tmp_path)
    assert fc.applies("V", "I", "H", 15007, 17493, c)
    assert not fc.applies("V", "I", "VIS", 4800, 6900, c)          # hyperfine, but not IR
    assert not fc.applies("Fe", "I", "H", 15007, 17493, c)        # IR, but no hyperfine line
    assert not fc.applies("V", "I", "K", 19455, 24842, c)         # no V hyperfine line in K here


def test_the_special_case_never_shares_the_standard_identity():
    from pipeline import run_matrix as rm
    from pipeline.run_descriptor import RunDescriptor
    base = dict(element="V", ion="I", instrument="crires_plus",
                holding="solar_crires_plus_h_rya1232", lo_A=15007.11, hi_A=17493.69,
                method="synthesis", pool="codex")
    std, free = RunDescriptor(**base), RunDescriptor(**base, continuum="free")
    assert std.key == free.key                                   # same window, same lines
    assert rm.pool_tier(std) == "GRADED"                         # standard cell unchanged
    assert rm.pool_tier(free) == "GRADED-FREECONT"
    assert rm.ledger_key("solar", std, "H") != rm.ledger_key("solar", free, "H")
    assert rm.cell_prefix_key(std) != rm.cell_prefix_key(free)


def test_the_flag_reaches_only_the_special_case_cell():
    from pipeline.run_descriptor import RunDescriptor, resolve
    base = dict(element="V", ion="I", instrument="crires_plus",
                holding="solar_crires_plus_h_rya1232", lo_A=15007.11, hi_A=17493.69,
                method="synthesis", pool="codex")

    def args(d):
        return next(s for s in resolve(d).steps if s["name"] == "derive_products")["args"]
    assert "--free-continuum" not in args(RunDescriptor(**base))
    assert "--free-continuum" in args(RunDescriptor(**base, continuum="free"))


def test_each_cell_takes_only_its_own_files(tmp_path, monkeypatch):
    import time
    from pipeline import run_matrix as rm
    from pipeline.run_descriptor import RunDescriptor
    monkeypatch.setattr(rm, "ARTIFACT_ROOT", tmp_path)
    out = tmp_path / "data/results/band_products"
    out.mkdir(parents=True)
    stem = "VI_15007_17493_crires_plus_solar_crires_plus_h_rya1232_SYNTH_GRADED"
    (out / f"{stem}_products.csv").write_text("x")
    (out / f"{stem}_FREECONT_products.csv").write_text("x")
    base = dict(element="V", ion="I", instrument="crires_plus",
                holding="solar_crires_plus_h_rya1232", lo_A=15007.11, hi_A=17493.69,
                method="synthesis", pool="codex")
    t0 = time.time() - 60
    std = [a["path"] for a in rm.artifacts_written(RunDescriptor(**base), t0)]
    free = [a["path"] for a in rm.artifacts_written(RunDescriptor(**base, continuum="free"), t0)]
    assert std == [f"data/results/band_products/{stem}_products.csv"]
    assert free == [f"data/results/band_products/{stem}_FREECONT_products.csv"]


def test_the_budget_publishes_beside_the_standard_selector_on_the_same_tier():
    import rya1230_cno_budget as b
    std = ["--lines-tier", "graded"]
    assert b._selector(std) == "GRADED" and b._tier(std) == "GRADED"
    assert b._selector(std + ["--free-continuum"]) == "GRADED-FREECONT"
    assert b._tier(std + ["--free-continuum"]) == "GRADED"
    s = ["--lines-from-set", "V_IR=x.csv", "--free-continuum"]
    assert b._selector(s) == "SET-V_IR-FREECONT"


def test_the_special_case_grades_as_the_pool_it_measured():
    from pipeline.cno_grade import grade_for
    for sel in ("GRADED", "REFERENCE", "DEEPGRADED"):
        assert grade_for({"element": "V", "selector": sel}) == \
            grade_for({"element": "V", "selector": sel + "-FREECONT"})
