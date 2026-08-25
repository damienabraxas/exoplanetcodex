"""The ⟨3D⟩ STAGGER atmosphere, in the shape iSpec's Turbospectrum path needs.

RYA-821/1035 made the ⟨3D⟩ DECKS readable. This module makes the ⟨3D⟩ ATMOSPHERE those
departures were solved on usable by the synthesis — the step the register has carried as
"still owed, and a PHYSICS one rather than plumbing" for both Al and Fe.

It turns out to be neither exactly: it is an INTERFACE CONTRACT, and every clause of it is
readable in iSpec's own source. The reason it looked like physics is that the obvious route
— hand iSpec the model and let it write a MARCS file — is genuinely impossible, and the
reason is structural rather than a bug:

    ispec/atmospheres.py:274   writes lgTauR, layer[7], layer[8], layer[1], layer[9],
                               layer[2], layer[5]   =   lgTau5, Depth, T, Pe, Pg, Prad

A ⟨3D⟩ mul23 model carries **five** columns — log τ₅₀₀, T, nₑ, V, v_mic — and has no
Depth, no Pₑ, no P_g, no P_rad. `write_atmosphere` cannot produce a MARCS model from it,
and inventing P_g to satisfy the format would be fabricating the structure the whole
exercise exists to measure. (⚠️ iSpec already fabricates one field there itself: `lgTauR`
is a counter stepped by 0.1/0.2 with the comment *"only needed by Turbospectrum to read
the model atmosphere but it has no effect"*.)

## The contract, derived from `ispec/synth/turbospectrum.py`

Pass the ⟨3D⟩ model as a FILE and give the array only what iSpec actually reads off it:

| line | code | what it forces |
| --- | --- | --- |
| 110 | `if atmosphere_layers_file is None:` | supplying the file SKIPS `write_atmosphere` |
| 114 | `is_marcs_model = len(atmosphere_layers[0]) == 11` | the array must **not** have 11 columns, so babsma is called `MARCS-FILE=.false.` |
| 159 | `interpolated_tau = atmosphere_layers[:, 7]` | the array must have **≥ 8** columns, and column **7** must be lgTau5 |
| 85–89 | `radius = …[0][-1]`; spherical needs `nvalues == 11` | ncols ≠ 11 ⇒ plane-parallel, which is correct for this model |
| 331 | `if remove_tmp_atm_file: os.remove(…)` | that flag is only set when iSpec MADE the file, so a caller-supplied path is **not** deleted |

⇒ **8 ≤ ncols ≤ 10, and column 7 is log τ₅₀₀.** This module emits 10.

`MARCS-FILE=.false.` is the mode RYA-442 measured as correct for a `TAU5000 SCALE` model
(101 depths, log τ −5→+5 and T round-tripping exactly); `.true.` fails loudly with *"This
model is probably not a MARCS model!"*, so the wrong mode is audible rather than silent.

## Why the τ overwrite at line 159 is a no-op here, and not a fudge

iSpec overwrites the departure file's τ with the atmosphere's to dodge Turbospectrum's
"tau scales differ" error. That is normally a papering-over. For a ⟨3D⟩ deck read against
its OWN atmosphere it is an identity: measured on the real Fe deck against this model,
**max |Δ log τ| = 5.0e-5 over all 101 depths**, and the residual is the model file's
four-decimal text precision (deck −3.09995 printed as −3.0999). `assert_tau_consistent`
turns that into a gate rather than a remark.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

#: iSpec reads `atmosphere_layers[:, 7]` as lgTau5 (turbospectrum.py:159), matching the
#: MARCS writer's own column map (atmospheres.py:274).
TAU_COLUMN = 7

#: 11 columns is iSpec's "this is a MARCS model" signal, which would send babsma to
#: `MARCS-FILE=.true.` and make it reject a TAU5000 model. Anything from 8 to 10 satisfies
#: the tau column without tripping it; 10 leaves room and keeps the last column (which
#: `radius` reads) unambiguous.
N_COLUMNS = 10

#: The ⟨3D⟩ mul23 column order, as written by the STAGGER→TSv20 conversion.
MUL23_COLUMNS = ("log_tau500", "T_K", "n_e_cm3", "V", "vmic_kms")


class Mean3DAtmosphereError(RuntimeError):
    """Raised loudly. A ⟨3D⟩ atmosphere that is silently wrong yields a well-formed
    spectrum for a structure nobody chose, which is the whole failure class here."""


def read_mul23(path: str | Path) -> dict:
    """Parse a mul23 `TAU5000 SCALE` model into its header and its five columns.

    Refuses anything it cannot identify rather than guessing: the file is small, the
    format is fixed, and a mis-parse here is a wrong atmosphere that still synthesises.
    """
    p = Path(path)
    text = p.read_text().splitlines()

    scale = next((ln.strip() for ln in text if ln.strip().endswith("SCALE")), None)
    if scale != "TAU5000 SCALE":
        raise Mean3DAtmosphereError(
            f"{p}: depth scale is {scale!r}, not 'TAU5000 SCALE'. This module only handles "
            f"the τ₅₀₀ form; a Rosseland-scaled model is a DIFFERENT depth variable and "
            f"pairing it with τ₅₀₀ departures would misplace every layer (RYA-1013).")

    # NDEP is declared in the header, and it is checked against what we actually read --
    # a truncated file otherwise yields a shorter, perfectly plausible atmosphere.
    ndep_declared = None
    for i, ln in enumerate(text):
        if ln.strip().lstrip("* ").upper().startswith("NDEP"):
            for nxt in text[i + 1:]:
                if nxt.strip() and not nxt.lstrip().startswith("*"):
                    try:
                        ndep_declared = int(nxt.split()[0])
                    except ValueError:
                        pass
                    break
            break

    rows, started = [], False
    for ln in text:
        parts = ln.split()
        if len(parts) != 5:
            if started:
                break
            continue
        try:
            rows.append([float(x) for x in parts])
            started = True
        except ValueError:
            if started:
                break
    if not rows:
        raise Mean3DAtmosphereError(f"{p}: no 5-column depth rows found")

    data = np.asarray(rows, dtype=float)
    if ndep_declared is not None and len(data) != ndep_declared:
        raise Mean3DAtmosphereError(
            f"{p}: header declares NDEP={ndep_declared} but {len(data)} depth rows were "
            f"read. A truncated model still synthesises -- refusing.")

    tau = data[:, 0]
    if not np.all(np.diff(tau) > 0):
        raise Mean3DAtmosphereError(
            f"{p}: log τ₅₀₀ is not strictly increasing. The depth ordering is the one thing "
            f"every downstream index depends on.")

    head = [ln.strip() for ln in text[:3]]
    try:
        teff, logg, feh = float(head[0]), float(head[1]), float(head[2])
    except (IndexError, ValueError) as e:
        raise Mean3DAtmosphereError(f"{p}: unreadable Teff/logg/[Fe/H] header: {e}") from e

    model_id = next((ln.strip() for ln in text[:12]
                     if ln.strip() and not ln.startswith("*")
                     and not ln.strip().endswith("SCALE")
                     and not _is_number(ln.strip())), "")

    return dict(path=str(p), teff=teff, logg=logg, feh=feh, model_id=model_id,
                ndep=len(data), columns=MUL23_COLUMNS, data=data)


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def as_ispec_layers(model: dict) -> np.ndarray:
    """The array iSpec needs beside the FILE — not a substitute for it.

    🔴 THIS ARRAY DOES NOT DESCRIBE THE ATMOSPHERE, AND MUST NOT BE READ AS THOUGH IT DID.
    babsma reads the real structure out of the mul23 file; iSpec touches this array in
    exactly three places (the column count, the last column, and column 7), so only those
    carry meaning. Every other entry is NaN **on purpose**: a zero would look like a
    physical value, and the first caller to read one would get a plausible number for a
    quantity this model does not carry.
    """
    data = model["data"]
    layers = np.full((len(data), N_COLUMNS), np.nan, dtype=float)
    layers[:, TAU_COLUMN] = data[:, 0]          # lgTau5 -- iSpec writes this as departure τ
    layers[:, N_COLUMNS - 1] = 0.0              # `radius`: plane-parallel (see module docs)
    return layers


def assert_tau_consistent(parsed_deck: dict, model: dict,
                          tol: float = 2e-4) -> float:
    """The ⟨3D⟩ deck and the ⟨3D⟩ atmosphere must describe the SAME depths.

    🔴 THIS IS THE GATE THAT MAKES iSpec's τ OVERWRITE HONEST. `turbospectrum.py:159`
    replaces the departure τ with the atmosphere's, which normally hides a disagreement.
    Here the two come from the same STAGGER model, so they should agree — and if they ever
    do not, the overwrite would silently apply departures at the wrong depths.

    Returns the measured max |Δ log τ|.

    ⚠️ THE TOLERANCE IS THE FILE'S PRINT PRECISION, DERIVED NOT CHOSEN. The mul23 model
    stores log τ to four decimals while the deck stores float64, so a half-digit
    disagreement (5e-5) is expected and means nothing physical. 2e-4 is 4× that — tight
    enough that a real regrid (0.1 dex spacing) cannot hide, loose enough that text
    rounding cannot fire it. Measured on the real Fe deck: **5.0e-5**.
    """
    deck_tau = np.asarray(parsed_deck["tau"], dtype=float)
    mod_tau = model["data"][:, 0]
    if len(deck_tau) != len(mod_tau):
        raise Mean3DAtmosphereError(
            f"depth count differs: deck has {len(deck_tau)}, atmosphere has {len(mod_tau)}. "
            f"iSpec pairs them by INDEX, so unequal lengths mean depth i of the departures "
            f"lands on a different physical layer.")
    dmax = float(np.abs(deck_tau - mod_tau).max())
    if dmax > tol:
        worst = int(np.argmax(np.abs(deck_tau - mod_tau)))
        raise Mean3DAtmosphereError(
            f"the deck and the atmosphere are on DIFFERENT depth scales: max "
            f"|Δ log τ₅₀₀| = {dmax:.3e} > {tol:.0e}, worst at depth {worst} "
            f"(deck {deck_tau[worst]:.6f} vs atmosphere {mod_tau[worst]:.6f}). iSpec "
            f"overwrites the departure τ with the atmosphere's, so this would be applied "
            f"silently and the departures would sit at the wrong depths.")
    return dmax


def assert_not_marcs_shaped(layers: np.ndarray) -> None:
    """🔴 ELEVEN COLUMNS IS A SWITCH, NOT A SIZE. `turbospectrum.py:114` sets
    `is_marcs_model = len(atmosphere_layers[0]) == 11` and hands it to `calculate_opacities`
    as babsma's `MARCS-FILE` flag. An 11-column array here would tell babsma to read a
    TAU5000 model as native MARCS, which RYA-442 measured as a LOUD failure -- but loud in
    babsma's log, several layers below anything that reports a result."""
    ncols = int(layers.shape[1])
    if ncols == 11:
        raise Mean3DAtmosphereError(
            "the layer array has 11 columns, which iSpec reads as 'this is a MARCS model' "
            "and turns into MARCS-FILE=.true. -- the mode that rejects a TAU5000 model.")
    if ncols <= TAU_COLUMN:
        raise Mean3DAtmosphereError(
            f"the layer array has {ncols} columns; iSpec reads column {TAU_COLUMN} as the "
            f"departure τ (turbospectrum.py:159), so it needs at least {TAU_COLUMN + 1}.")


def load(path: str | Path) -> tuple[np.ndarray, dict]:
    """`(layers, model)` ready to pass beside `atmosphere_layers_file=model['path']`."""
    model = read_mul23(path)
    layers = as_ispec_layers(model)
    assert_not_marcs_shaped(layers)
    return layers, model
