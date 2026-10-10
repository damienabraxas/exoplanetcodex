"""RYA-1233: the synthesis route honours problem_children exclusions, as the EW route does.

`_stamp` (RYA-807) applied the registry on the EW route only; synthesis_route returned
before it, so a registered cull was ignored by every synthesis product. Si's Codex-graded
KP NIR pool aggregated 7 lines culled for the Sun.
"""
import inspect
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))


def test_the_synthesis_fit_loop_consults_the_registry():
    import derive_band_products as d
    src = inspect.getsource(d.synthesis_route)
    assert "_pc.aggregate_action" in src and "line_dispositions" in src


def test_a_culled_si_line_resolves_to_exclude():
    from pipeline import problem_children as pc
    t = pc.line_dispositions()
    for w in (10660.973, 12031.504, 10827.088):
        assert pc.aggregate_action(pc.disposition_for_line("Si", "I", w, table=t)) == "exclude"
