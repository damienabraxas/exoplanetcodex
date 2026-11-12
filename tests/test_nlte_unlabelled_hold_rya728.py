"""RYA-728: derive's "NOT ONE pooled line carries an NLTE label" refusal is a HOLD, not a FAILED.

8 Cr cells (HARPS Cr II VIS, Cr I red-optical; gerber-nlte / gerber-mean3d, graded pool of
1-2 lines) were refused correctly by derive and then filed as FAILED. The orchestrator
recognises the refusal by a marker string, so the marker must stay in derive's message.
"""
from pathlib import Path

from pipeline import run_matrix as rm

ROOT = Path(__file__).resolve().parents[1]
DERIVE = (ROOT / "scripts" / "derive_band_products.py").read_text(encoding="utf-8")


def test_every_marker_is_in_derives_refusal_message():
    for m in rm.NLTE_UNLABELLED_POOL_MARKERS:
        # derive builds the message across f-string pieces; the marker sits in one piece
        assert m in DERIVE, m


def test_the_marker_does_not_also_classify_as_an_empty_pool():
    msg = ("derive_products: exit 1: NLTE was requested for Cr I but NOT ONE of the 2 pooled "
           "lines carries an NLTE label, so bsyn would apply departure = 1")
    assert any(m in msg for m in rm.NLTE_UNLABELLED_POOL_MARKERS)
    assert not any(m in msg for m in rm.EMPTY_POOL_MARKERS)
