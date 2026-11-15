"""RYA-1388: a patch refits the listed lines and rebuilds every other line's fit result from
its stored row, so the run's own gates and aggregate reproduce a full rerun."""
import csv

import numpy as np
import pytest

from pipeline import line_patch as lp


def _write(path, rows):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return path


def _row(w, ep, a="7.40", reason="", level="0.99", method="RYA-1230 model-guided; ok"):
    return {"wavelength_air_A": w, "ep_eV": ep, "abundance": a, "excluded_reason": reason,
            "sigma_A": "0.05", "frac_rise_weaker": "4.2", "edge_distance_dex": "2.5",
            "red_chi2": "", "continuum_level": level, "continuum_method": method}


def test_the_stored_row_gives_back_what_the_fit_loop_reads():
    res = lp.res_from_row(_row("5000.0", "2.0"), "holding=x")
    assert res["a_synth"] == 7.40 and res["status"] == "ok" and res["reason"] == ""
    assert res["frac_rise_weaker"] == 4.2 and res["obs_source"] == "holding=x"
    assert np.isnan(res["red_chi2"])          # '' was None; as_float_or_none maps NaN back


def test_a_failed_fit_keeps_its_status_behind_a_registry_prefix():
    r = _row("5000.0", "2.0", a="", reason="REGISTRY-BLEND: exclude/active per x | "
                                           "SYNTHESIS: failed too few observed pixels")
    res = lp.res_from_row(r, "")
    assert (res["status"], res["reason"]) == ("failed", "too few observed pixels")
    assert np.isnan(res["a_synth"])
    # the loop rebuilds exactly the text it wrote
    assert f"{lp.SYNTH_PREFIX}{res['status']} {res['reason']}".strip() == r[
        "excluded_reason"].split(" | ")[1]


def test_a_moved_continuum_is_drift_and_an_unmoved_one_is_not():
    r = _row("5000.0", "2.0")
    assert lp.continuum_drift(r, 0.99, "RYA-1230 model-guided; ok") == ""
    assert "level" in lp.continuum_drift(r, 0.97, "RYA-1230 model-guided; ok")
    assert "method" in lp.continuum_drift(r, 0.99, "RYA-1230 model-guided; refused")
    blank = _row("5000.0", "2.0", level="", method="")
    assert lp.continuum_drift(blank, None, "") == ""


def test_refit_lines_and_stored_fits_come_together(tmp_path):
    assert lp.LinePatch.from_args(None, None, "s", "Fe I") is None
    with pytest.raises(lp.LinePatchError):
        lp.LinePatch.from_args(str(_write(tmp_path / "r.csv", [{"wavelength_air_A": 1}])),
                               None, "s", "Fe I")


def test_a_wavelength_not_in_the_pool_is_refused(tmp_path):
    refit = _write(tmp_path / "r.csv", [{"species": "Fe I", "wavelength_air_A": "5000.003",
                                         "loggf": "-9"},
                                        {"species": "Fe I", "wavelength_air_A": "5010.5"}])
    p = lp.LinePatch.from_args(str(refit), str(tmp_path), "s", "Fe I")
    with pytest.raises(lp.LinePatchError, match="5010.500"):
        p.bind(np.array([5000.0, 5005.0]), np.array([2.0, 3.0]))


def test_other_species_rows_are_left_to_their_own_run(tmp_path):
    refit = _write(tmp_path / "r.csv", [{"species": "Fe II", "wavelength_air_A": "6000.0"},
                                        {"species": "Fe 1", "wavelength_air_A": "5005.001"}])
    p = lp.LinePatch.from_args(str(refit), str(tmp_path), "s", "Fe I", ("Fe 1",))
    p.bind(np.array([5000.0, 5005.0]), np.array([2.0, 3.0]))
    assert p.refit_wl == {5005.0}


def test_listed_lines_refit_others_reuse_and_a_missing_store_refits(tmp_path):
    _write(tmp_path / "s_1D-LTE_lines.csv", [_row("5000.0", "2.0"), _row("5005.0", "3.0")])
    refit = _write(tmp_path / "r.csv", [{"wavelength_air_A": "5005.0"}])
    p = lp.LinePatch.from_args(str(refit), str(tmp_path), "s", "Fe I")
    p.bind(np.array([5000.0, 5005.0, 5009.0]), np.array([2.0, 3.0, 1.0]))
    assert p.prior("1D-LTE", 5000.0, 2.0)["abundance"] == "7.40"
    assert p.prior("1D-LTE", 5005.0, 3.0) is None                      # listed
    assert p.prior("1D-LTE", 5009.0, 1.0) is None                      # never stored
    assert p.prior("1D-LTE", 5000.0, 4.5) is None                      # EP says another line
    assert p.prior("ENGINE-B-NLTE", 5000.0, 2.0) is None               # no store: refit
    led = p.ledger["1D-LTE"]
    assert led["reused"] == 1 and led["refit"] == [5005.0]
    assert led["refit_no_prior"] == [5009.0, 5000.0]
    assert p.write(tmp_path).name == "s_patch.json"


def test_drift_on_a_reused_line_refuses_the_patch(tmp_path):
    p = lp.LinePatch.from_args(str(_write(tmp_path / "r.csv", [{"wavelength_air_A": 1}])),
                               str(tmp_path), "s", "Fe I")
    p.assert_no_drift()
    p.note_drift("1D-LTE", 5000.0, "level 0.99 -> 0.97")
    with pytest.raises(lp.LinePatchError, match="5000.000"):
        p.assert_no_drift()


def _leg(store, sig, leg, args, extra=(), done=True):
    d = store / sig / leg
    d.mkdir(parents=True)
    import shlex
    d.joinpath("argv.txt").write_text(shlex.join(
        ["py", "scripts/derive_band_products.py", *args, *extra, "--out", str(d)]) + "\n{}\n")
    if done:
        d.joinpath("DONE").write_text("0\n")


def test_a_unit_finds_its_own_earlier_legs_by_command(tmp_path):
    from pipeline.element_publish import leg_priors
    a = ["--element", "Fe", "--ion", "I", "--lines-tier", "reference"]
    b = ["--element", "Fe", "--ion", "II", "--lines-tier", "reference"]
    _leg(tmp_path, "s1", "cscale", a)
    _leg(tmp_path, "s1", "core", a, extra=("--half-width-A", "0.25"))
    _leg(tmp_path, "s2", "cscale", b, done=False)                 # killed leg: no prior
    units = [("u0", "python3 scripts/derive_band_products.py " + " ".join(a), "ts", "k0"),
             ("u1", "python3 scripts/derive_band_products.py " + " ".join(b), "ts", "k1"),
             ("u2", "python3 scripts/derive_band_products.py " + " ".join(a), "ts", "k2")]
    assert leg_priors(tmp_path, units, kept=[2]) == {0: str(tmp_path / "s1")}


def test_two_filed_builds_of_one_command_are_ambiguous(tmp_path):
    from pipeline.element_publish import leg_priors
    a = ["--element", "Fe", "--ion", "I"]
    _leg(tmp_path, "s1", "cscale", a)
    _leg(tmp_path, "s2", "cscale", a)
    units = [("u0", "python3 scripts/derive_band_products.py " + " ".join(a), "ts", "k0")]
    assert leg_priors(tmp_path, units, kept=[]) == {}


def test_a_line_list_runs_only_as_a_synthesis_patch():
    from pipeline.run_request import RunRequest
    ok = dict(star="solar", element="Fe", mode="select", ions=["I"], routes=["synthesis"],
              lines="x.csv", patch=True)
    assert RunRequest.from_dict(ok).problems() == []
    assert RunRequest.from_dict({**ok, "patch": False}).problems()
    assert RunRequest.from_dict({**ok, "routes": ["profile-fit"]}).problems()
    assert RunRequest.from_dict({**ok, "lines": ""}).problems()
