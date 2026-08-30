"""
scripts/rya1118_multi3d_read.py
==============================
RYA-1118 (the Bride, M2/M3) — read a Multi3D-format STAGGER snapshot into SI.

This is OURS. Lightweaver and the STAGGER cube are external and stay out of this
repo (G3, as with VoronoiRT in RYA-1012); the evaluation tree lives at
`~/bride_rya1118/` on the Mac. What is committed here is the ingest — the one piece
the Bride actually needs from us, and the piece a later rung (M4/RYA-1119) reuses.

WHY THIS IS A READER AND NOT AN EOS CONVERSION. The obvious fear about feeding a
STAGGER cube to an RT code is the equation of state: raw STAGGER `atm3d` carries
density and internal energy, and temperature / electron density have to be looked up
through an EOS table whose python is broken-as-shipped. **That work is already done
in these files.** The snapshot's own `info.` header says so:

    staggtom3d
    Reference depth scale: Rosseland
    Input EOS table: EOSrhoe.tab
    Output Multi3D atmosphere: atm3d.t5777g44m0005_00091_80x80x240

so what we hold is the POST-EOS, Multi3D-format product on a Rosseland depth scale —
the "synthesis-ready processed snapshot" the ticket asks for. This module reads it and
nothing more.

LAYOUT, taken from Muspel.jl's `read_atmos_multi3d` (Pereira; the reference
implementation for this format) rather than reverse-engineered:

  mesh file : whitespace-numeric — nx, x[nx], ny, y[ny], nz, z[nz]   (cm)
  atm3d file: SIX consecutive Float32 arrays of shape (nx, ny, nz), in this order
              electron_density [cm^-3] · temperature [K] · vx · vy · vz [km/s] ·
              rho [g cm^-3]

Byte arithmetic checks the layout independently: 6 x 80 x 80 x 240 x 4 = 36,864,000,
which is exactly the file size. A wrong variable count would not land on it.

Hydrogen is split by the Saha ionisation fraction, as Muspel does, using the same
grams-per-hydrogen constant (2.380491e-24 g). SI is returned throughout, because that
is what Lightweaver's `Atmosphere` wants and a unit seam is where this would break.
"""
from __future__ import annotations

import numpy as np

#: grams of atmosphere per hydrogen atom — Muspel's `grph` default, same value.
GRPH_G = 2.380491e-24
CM_TO_M = 1e-2
KM_TO_M = 1e3

#: Hydrogen ionisation energy [eV] and the Saha prefactor pieces, in SI.
_CHI_H_EV = 13.5984
_EV_TO_J = 1.602176634e-19
_K_B = 1.380649e-23
_H_PLANCK = 6.62607015e-34
_M_E = 9.1093837015e-31


def read_mesh(path):
    """(nx, ny, nz, x, y, z) in METRES. The file is whitespace-numeric, counts inline."""
    vals = np.fromstring(open(path).read().replace("\n", " "), sep=" ", dtype=np.float64)
    i = 0
    nx = int(vals[i]); i += 1
    x = vals[i:i + nx]; i += nx
    ny = int(vals[i]); i += 1
    y = vals[i:i + ny]; i += ny
    nz = int(vals[i]); i += 1
    z = vals[i:i + nz]; i += nz
    # A SHORT read and a LONG read are different faults and must say so separately.
    # numpy slicing silently truncates, so `z` can come back shorter than nz while the
    # leftover arithmetic goes NEGATIVE -- which reads as gibberish rather than as the
    # real problem (the mesh claims more depth points than it carries).
    if z.size != nz:
        raise ValueError(
            f"{path}: header claims nz={nz} but only {z.size} depth values follow. "
            "The mesh is truncated or its counts are wrong — refusing to guess.")
    if i != vals.size:
        raise ValueError(
            f"{path}: {vals.size - i} numbers left over after nx={nx} ny={ny} nz={nz}. "
            "The mesh does not describe this file — refusing to guess its shape.")
    return nx, ny, nz, x * CM_TO_M, y * CM_TO_M, z * CM_TO_M


def h_ionfrac_saha(temperature, electron_density):
    """Fraction of hydrogen that is ionised, from Saha. Mirrors Muspel's helper.

    Kept explicit rather than pulled from a library so the one place hydrogen is split
    is visible: an error here moves every line's opacity and would look like a bad atom.
    """
    t = np.asarray(temperature, dtype=np.float64)
    ne = np.asarray(electron_density, dtype=np.float64)
    # 2 (2 pi m_e k T / h^2)^{3/2} / n_e * exp(-chi/kT), with g_II/g_I = 1/2
    fac = 2.0 * (2.0 * np.pi * _M_E * _K_B * t / _H_PLANCK ** 2) ** 1.5
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        ratio = fac / ne * np.exp(-_CHI_H_EV * _EV_TO_J / (_K_B * t))
        frac = ratio / (1.0 + ratio)
    return np.clip(np.nan_to_num(frac, nan=1.0, posinf=1.0), 0.0, 1.0)


def read_atmos_multi3d(mesh_file, atmos_file, *, grph_g: float = GRPH_G):
    """Read a Multi3D snapshot into SI arrays shaped (nz, ny, nx).

    Returns a dict with z/y/x, temperature, electron_density, hydrogen_density
    (NEUTRAL H), proton_density, and vx/vy/vz.
    """
    nx, ny, nz, x, y, z = read_mesh(mesh_file)

    n = nx * ny * nz
    expect = 6 * n * 4
    import os
    got = os.path.getsize(atmos_file)
    if got != expect:
        raise ValueError(
            f"{atmos_file}: {got} bytes, expected {expect} for six float32 "
            f"({nx}x{ny}x{nz}) arrays. The mesh and the cube disagree; reading it "
            "anyway would silently reinterpret one variable as another.")

    raw = np.fromfile(atmos_file, dtype=np.float32, count=6 * n)
    # (nx, ny, nz) on disk -> (nz, ny, nx), matching Muspel's permutedims(3, 2, 1)
    def _f(k):
        return np.transpose(raw[k * n:(k + 1) * n].reshape(nx, ny, nz), (2, 1, 0)).astype(np.float64)

    electron_density = _f(0) / CM_TO_M ** 3          # cm^-3 -> m^-3
    temperature = _f(1)
    vx, vy, vz = (_f(k) * KM_TO_M for k in (2, 3, 4))
    rho = _f(5)                                       # g cm^-3

    n_h_total = rho / (grph_g * CM_TO_M ** 3)         # -> m^-3
    ionfrac = h_ionfrac_saha(temperature, electron_density)
    proton_density = n_h_total * ionfrac
    hydrogen_density = n_h_total * (1.0 - ionfrac)

    return dict(nx=nx, ny=ny, nz=nz, x=x, y=y, z=z,
                temperature=temperature, electron_density=electron_density,
                hydrogen_density=hydrogen_density, proton_density=proton_density,
                vx=vx, vy=vy, vz=vz)
