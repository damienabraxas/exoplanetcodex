import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_run_is_exact_source_set_and_never_substitutes():
    tree = ast.parse((ROOT / "scripts/run_si_asplund_ew_rya1169.py").read_text())
    text = ast.unparse(tree)
    assert "src.asplund_grade == 'asplund'" in text
    assert "published_loggf" in text
    assert "SCOTT_EW_MA[centre]" in text
    assert "no missing line substituted" in text


def test_run_checks_tellurics_before_kp_measurement():
    text = (ROOT / "scripts/run_si_asplund_ew_rya1169.py").read_text()
    assert "exclusion(centre, \"kpno_solar_atlas\")" in text
