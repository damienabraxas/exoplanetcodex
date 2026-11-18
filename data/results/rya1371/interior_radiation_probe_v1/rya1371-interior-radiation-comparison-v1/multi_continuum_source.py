"""Explicit MULTI wavelength-table semantics for the optional LW backend.

MULTI2.3 ATOM/FREQC consume absolute tabulated wavelengths when QMAX<0;
they do not discard points beyond the level-energy threshold. Its Milne factor
uses the actual photon frequency and LTE population ratio (mul23_sub.f, 1010ff).
Source audited at MingjieJian/pymulti commit 4d6aa3c39d634b38da772e469a5692435eba1354,
pymulti/files/multi/mul23/source/mul23_sub.f:540-572,2325-2376.
This preserves an existing model; it does not validate its atomic physics.
"""
from dataclasses import dataclass
import numpy as np
from lightweaver.atomic_model import ExplicitContinuum


@dataclass(eq=False,repr=False)
class MultiExplicitContinuum(ExplicitContinuum):
    """Use the supplied table extent, while retaining energy-derived lambda0.

    LW's ordinary ExplicitContinuum clips to an energy-derived lambdaEdge and
    can append a fabricated table endpoint. MULTI's explicit-table branch does
    neither. Retain both definitions for audit, with unchanged source values.
    """
    @property
    def lambdaEdge(self):
        return float(np.max(self.wavelengthGrid))

    @property
    def energy_edge_nm(self):
        return super().lambdaEdge
