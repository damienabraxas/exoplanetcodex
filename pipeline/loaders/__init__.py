"""
pipeline/loaders/
=================
Instrument-specific spectrum loaders.

Each loader ingests a raw FITS file and returns a standardised
SpectrumData object: wave_A (Å, barycentric), flux, err, meta dict.

Supported instruments:
  SPIRou   — CFHT/APERO _t.fits (telluric-corrected)  ← this session
  HARPS    — ESO S1D_A                                  ← future
  HARPS-N  — TNG/IA2 S1D_A                              ← future
  HIRES    — Keck/KOA                                   ← future
  CRIRES+  — ESO/cr2res + molecfit                      ← future
"""

from .spirou_loader import SPIRouLoader
from .base_loader import SpectrumData

__all__ = ['SPIRouLoader', 'SpectrumData']
