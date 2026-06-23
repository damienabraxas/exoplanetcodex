"""
pipeline/acen_orbit.py
======================
RYA-423 — the alpha Cen AB visual+spectroscopic orbit, for predicting the heliocentric RV of
A and B at any epoch (the model-independent PRIMARY star-ID discriminator for IR frames the
optical classifier cannot reach).

Elements: Kervella, Thevenin & Lovis 2016 (A&A 594, A107, "Close stellar conjunctions of
alpha Cen A and B until 2050"), Table 1, which refines Pourbaix & Boffin 2016 (A&A 586, A90)
with 11 yr of HARPS RVs. Pulled from the source, NOT memory (RYA-423 acceptance).

  P     = 79.91 yr           e   = 0.5179
  T     = 1875.88 (Besselian periastron)
  a     = 17.571 arcsec      omega = 231.81 deg (component B relative to A)
  i     = 79.20 deg          Omega = 139.67 deg     parallax = 747.17 mas
  gamma = -22.3930 km/s      M_A = 1.1055 Msun       M_B = 0.9343 Msun

K_A / K_B are DERIVED from (a, i, P, e, masses), giving K_A ~ 4.6, K_B ~ 5.5 km/s. The A/B
assignment of the two RV branches (which star is more negative now) is pinned by the
INDEPENDENT RYA-384 spectral-type classification (J-band line depth + calibrated flux) — the
confirmed-A NIRPS frames sit at the more-negative branch (~-26.7), so alpha Cen A = the
omega=231.81 branch. This is NOT circular: it uses the abundance-independent spectral type,
not the CO/abundance the IR ID is meant to protect.
"""
from __future__ import annotations
import numpy as np
from astropy.time import Time

# --- Kervella 2016 Table 1 ---
P_YR = 79.91
ECC = 0.5179
T_PERI_BYEAR = 1875.88
A_ARCSEC = 17.571
OMEGA_DEG = 231.81           # component B relative to A
INC_DEG = 79.20
PARALLAX_AS = 0.74717
GAMMA = -22.3930             # km/s, systemic
M_A, M_B = 1.1055, 0.9343    # Msun
SOURCE = "Kervella, Thevenin & Lovis 2016 (A&A 594, A107) Table 1; Pourbaix & Boffin 2016 (A&A 586, A90)"

_AUYR_KMS = 4.740470         # 1 AU/yr in km/s
_MTOT = M_A + M_B
_A_AU = A_ARCSEC / PARALLAX_AS
# semi-major axes of each star about the barycentre
_aA_AU = _A_AU * M_B / _MTOT
_aB_AU = _A_AU * M_A / _MTOT
_sini = np.sin(np.radians(INC_DEG))
K_A = (2 * np.pi / P_YR) * (_aA_AU * _sini) / np.sqrt(1 - ECC ** 2) * _AUYR_KMS
K_B = (2 * np.pi / P_YR) * (_aB_AU * _sini) / np.sqrt(1 - ECC ** 2) * _AUYR_KMS

# A occupies the omega=231.81 branch (pinned by the independent spectral-type-confirmed-A
# NIRPS frames landing at the more-negative predicted RV ~-26.7); B is the opposite branch.
_OMEGA_A = np.radians(OMEGA_DEG)
_OMEGA_B = np.radians(OMEGA_DEG) + np.pi


def _true_anomaly(byear: float) -> float:
    M = 2 * np.pi * (((byear - T_PERI_BYEAR) / P_YR) % 1.0)
    E = M
    for _ in range(60):
        E = E - (E - ECC * np.sin(E) - M) / (1 - ECC * np.cos(E))
    return 2 * np.arctan2(np.sqrt(1 + ECC) * np.sin(E / 2), np.sqrt(1 - ECC) * np.cos(E / 2))


def _rv(byear: float, K: float, omega: float) -> float:
    nu = _true_anomaly(byear)
    return GAMMA + K * (np.cos(nu + omega) + ECC * np.cos(omega))


def predicted_rv(mjd: float) -> dict:
    """Predicted heliocentric RV (km/s) of alpha Cen A and B at a given MJD, plus their
    separation. delta_AB = RV_B - RV_A (the binary RV split the frames must reflect)."""
    by = Time(mjd, format='mjd').byear
    rva = _rv(by, K_A, _OMEGA_A)
    rvb = _rv(by, K_B, _OMEGA_B)
    return {'mjd': float(mjd), 'byear': float(by), 'rv_A': float(rva), 'rv_B': float(rvb),
            'delta_AB': float(rvb - rva), 'sep_kms': float(abs(rvb - rva))}


if __name__ == '__main__':
    print(f"alpha Cen AB orbit — {SOURCE}")
    print(f"derived K_A={K_A:.3f} K_B={K_B:.3f} km/s; gamma={GAMMA}")
    for yr in (2022.0, 2023.3, 2024.2, 2025.2):
        mjd = Time(yr, format='byear').mjd
        p = predicted_rv(mjd)
        print(f"  {yr}: RV_A={p['rv_A']:+.2f}  RV_B={p['rv_B']:+.2f}  |sep|={p['sep_kms']:.2f} km/s")
