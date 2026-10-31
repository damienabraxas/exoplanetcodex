"""RYA-1315 -- ENGINE-A corrects every stage its table carries (`nlte_ions`), not only the
registry's reported-ion lock: Ti II was left 1D-LTE although the Mallinson table holds it."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import nlte_corrections as N                # noqa: E402
from config.constants import NLTE_CORRECTION_ELEMENTS     # noqa: E402


def test_ti_corrects_both_stages_and_keeps_its_reported_ion_lock():
    assert N.nlte_ions("Ti") == (1, 2)
    assert NLTE_CORRECTION_ELEMENTS["Ti"]["ion"] == 1
    assert N.nlte_ions("Ca") == (int(NLTE_CORRECTION_ELEMENTS["Ca"]["ion"]),)


def test_each_stage_reads_only_its_own_rows():
    one = N._load_mpia_element_grid("Ti", "I")
    two = N._load_mpia_element_grid("Ti", "II")
    assert one["ion"] == 1 and two["ion"] == 2
    near = lambda waves: any(abs(w - 5336.77) < 0.06 for w in waves)
    assert near(two["waves"]) and not near(one["waves"])
    # the Asplund-set Ti II line gets its own (small, negative) solar delta
    d = N._mpia_element_delta("Ti", 5336.77, 5772, 4.44, 0.0, tol=0.06, ion="II")
    assert -0.05 < d < 0.0


def test_default_call_is_the_registry_ion():
    assert N._load_mpia_element_grid("Ti")["ion"] == 1


def test_ti_model_axis_is_stated_not_inferred():
    import derive_band_products as d
    assert d.engine_a_model("Ti") == "amarsi"
