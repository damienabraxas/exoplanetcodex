"""RYA-731 -- a lab line is matched to its own TRANSITION, not to whichever synthesis row
sits nearest in wavelength (Ni I 3566.372 / 3413.476 took a coincident VALD3 line)."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import derive_band_products as D  # noqa: E402

# Ni I 3566.372: an unrelated EP 3.679 line first, then the lab line's hfs rows (EP 0.423).
W = np.array([3566.371, 3566.372, 3566.372, 3566.372, 3566.372, 3566.372])
EP = np.array([0.4228, 3.6785, 0.4228, 0.4228, 0.4228, 0.4228])
IDX = np.arange(len(W))


def test_without_ep_the_nearest_row_is_kept_as_before():
    keep, missing = D._match_into_list(np.array([3566.372]), W, IDX)
    assert keep == [1] and not missing          # the EP 3.679 row: the old behaviour


def test_with_ep_the_row_of_the_wanted_transition_wins():
    keep, _ = D._match_into_list(np.array([3566.372]), W, IDX,
                                 want_ep=np.array([0.4228]), ep_all=EP)
    assert EP[keep[0]] == 0.4228
    assert W[keep[0]] == 3566.372               # nearest of the agreeing rows


def test_a_nearer_row_of_another_transition_does_not_win():
    # Ni I 3413.476 (EP 0.165): the EP 3.539 row at .479 vs own components at .472/.479
    w = np.array([3413.472, 3413.479, 3413.479])
    ep = np.array([0.1652, 3.5388, 0.1652])
    keep, _ = D._match_into_list(np.array([3413.476]), w, np.arange(3),
                                 want_ep=np.array([0.1652]), ep_all=ep)
    assert ep[keep[0]] == 0.1652


def test_no_agreeing_ep_falls_back_to_the_nearest_row():
    keep, missing = D._match_into_list(np.array([3566.372]), W, IDX,
                                       want_ep=np.array([2.0]), ep_all=EP)
    assert keep == [1] and not missing
