"""RYA-1315 -- Ti's Gerber decks: 1D MARCS and <3D> STAGGER on atom.ti503b, both ions
NLTE-capable (measured), and Ti labels beside every VALD-built synthesis list."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import gerber_nlte as G                     # noqa: E402


def test_both_ti_decks_load_the_same_atom_on_their_own_atmospheres():
    one, three = G.DECKS["Ti"], G.DECKS["Ti@mean3D"]
    assert one["atom"] == three["atom"] == "atom.ti503b"
    assert one["grid"] != three["grid"] and "STAGGERmean3D" in three["grid"]
    assert three.get("read_via") == "direct" and one.get("read_via") is None
    assert three["atmos"].endswith("sun_avg3d_stagger.mod")


def test_atom_ti503b_serves_ti_i_and_ti_ii():
    assert G.nlte_ion_capability("Ti", "I")[0] is True
    assert G.nlte_ion_capability("Ti", "II")[0] is True
    assert G.nlte_ion_capability("Ti@mean3D", "II")[0] is True


def test_ti_sidecar_rows_label_both_levels():
    n = 0
    for p in (ROOT / "data" / "linelists").glob("ispec_*/nlte_labels.tsv"):
        s = pd.read_csv(p, sep="\t")
        ti = s[s.element.str.startswith("Ti ")]
        n += len(ti)
        assert ti.nlte_level_low.notna().all() and ti.nlte_level_up.notna().all(), p
    assert n > 1000
