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


def test_run_uses_only_corrected_1984_and_2005():
    text = (ROOT / "scripts/run_si_asplund_ew_rya1169.py").read_text()
    assert "load_kp1984_composite_window" in text
    assert "_read_kurucz2005_residual" in text
    assert "load_kp_window(" not in text
    assert '"raw_1984_policy":"QUARANTINED; never loaded by this runner"' in text
