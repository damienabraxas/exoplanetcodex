"""
tests/test_ace_co_reopen_derisk_rya443.py
=========================================
RYA-443 — guard the de-risk decision arithmetic (docs/decisions/ace_co_reopen_derisk.md).
Pinned PRIMARY-SOURCE numbers; the verdict logic must stay internally consistent.
No synthesis here (the 8.646 reproducibility is the 441 CLI smoke test).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import SOLAR_ASPLUND2021  # noqa: E402

A_REF = float(SOLAR_ASPLUND2021['C'])            # 8.46 (single source)
OUR_1D_MU1 = 8.646                               # RYA-441 ATLAS9 mu=1 (bit-reproducible)

# Amarsi 2021 (A&A 656, A113), 12C16O Delta-nu=2 overtone, Li2015 gf (= our gf):
AMARSI_CO_1D_MARCS = 8.608
AMARSI_CO_1D_ATMO = 8.621
AMARSI_CO_3D = 8.479
AMARSI_CO_DIFF_3D_1D = AMARSI_CO_3D - AMARSI_CO_1D_MARCS      # -0.129
AMARSI_CH_DIFF_3D_1D = +0.045                                  # CH A-X, opposite sign
FAITHFUL_BAR = 0.07
NEAR_REF = 0.10


def test_step1_pipeline_faithful_same_gf():
    # Our 1D CO matches a modern same-gf (Li2015) 1D CO analysis within the bar.
    assert abs(OUR_1D_MU1 - AMARSI_CO_1D_MARCS) < FAITHFUL_BAR
    assert abs(OUR_1D_MU1 - AMARSI_CO_1D_ATMO) < FAITHFUL_BAR


def test_step2_published_co_differential_lands_in_band():
    # The PUBLISHED CO-specific 1D->3D differential, applied to OUR faithful 1D value.
    corrected = OUR_1D_MU1 + AMARSI_CO_DIFF_3D_1D
    assert round(corrected, 3) == 8.517
    assert abs(corrected - A_REF) <= NEAR_REF


def test_step3_co_negative_ch_positive_opposite_signs():
    assert AMARSI_CO_DIFF_3D_1D < 0          # CO toward reference (canonical)
    assert AMARSI_CH_DIFF_3D_1D > 0          # CH away (modern STAGGER)
    assert (AMARSI_CO_DIFF_3D_1D < 0) != (AMARSI_CH_DIFF_3D_1D < 0)   # opposite


def test_verdict_corroborated():
    faithful = abs(OUR_1D_MU1 - AMARSI_CO_1D_MARCS) < FAITHFUL_BAR
    diff_in_band = abs((OUR_1D_MU1 + AMARSI_CO_DIFF_3D_1D) - A_REF) <= NEAR_REF
    co_negative = AMARSI_CO_DIFF_3D_1D < 0
    assert faithful and diff_in_band and co_negative   # all three -> CORROBORATED


def test_decision_record_present_with_verdict_and_gf_check():
    doc = ROOT / 'docs' / 'decisions' / 'ace_co_reopen_derisk.md'
    assert doc.exists()
    txt = doc.read_text()
    assert 'CORROBORATED' in txt
    norm = ' '.join(txt.split())                      # collapse newlines/spaces
    assert 'gf cross-check PASSES' in norm or 'gf cross-check passes' in norm
    assert 'Li et al. 2015' in txt or 'Li2015' in txt
    assert 'OPPOSITE' in txt and 'CH' in txt          # sign resolution stated
    assert '8.646' in txt and '8.517' in txt
