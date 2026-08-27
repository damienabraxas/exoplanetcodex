"""RYA-1077 — a STABLE identity for a canonical_gf row.

🔴 THE DEFECT THIS REPLACES. `canonical_gf.line_id` is `gf_000000`, `gf_000001`, … —
assigned by ROW POSITION. It is a row number wearing an identifier's clothes, and it moves
whenever a block of rows is removed or replaced (appending alone is safe; RYA-1052 finding
RYA-1047's lab lines unselectable and RYA-1053 re-extending is not).

Measured on the committed artifacts at the time of writing:

    artifact rows carrying a gf_ id and a wavelength : 6,979
    id NOW POINTS AT A DIFFERENT LINE                : 1,739   (25%)

    gf_177842  recorded against 21553.299 A  ->  now resolves to 21816.566 A
    gf_177933  recorded against 21735.457 A (Fe I) -> now Ni I at 21993.985 A

Not near-misses: different lines, and in places a different SPECIES. The damage is
concentrated exactly where you would expect from a block replacement — rya1058 28%,
rya1059 66%, and rya1060 0.1% because it was generated afterwards.

🔴 WHY THERE IS NO ROUNDING HERE, WHICH IS THE WHOLE DESIGN.
The obvious fix is to hash a rounded physical key. RYA-1033 already paid for that lesson:
a rounded number is not an identity, and `round(6136.615, 2)` is 6136.61 in Python and
6136.62 in numpy/pandas. And the data does not permit it — canonical_gf's wavelengths carry
a MEDIAN of 3 decimals but 13.3% carry four or more (26.0% for the excitation potential),
so any fixed precision would fuse lines that genuinely differ.

So the key is built from `repr(float(x))` — the shortest string that round-trips a float64
EXACTLY. The precision is therefore DERIVED from the stored value rather than chosen, and
two rows collide only if they are the same transition to the last bit.

Measured on the current table: **178,680 rows, 178,680 distinct keys, 0 collisions.**

⚠️ A collision is REFUSED, never broken by a tiebreak. Two rows sharing a physical key are
either a duplicate that needs removing or a key that is too coarse — both want a human, and
silently keeping the first is how a join starts lying (RYA-1033 again).
"""
from __future__ import annotations

import hashlib

#: Length of the hex digest kept. 12 hex chars = 48 bits; at ~1.8e5 rows the birthday
#: probability of any collision is ~3e-10, and a collision RAISES rather than being
#: resolved, so the cost of the unlikely case is a loud failure and not a wrong join.
_DIGEST_CHARS = 12
PREFIX = "pk_"


class PhysicalKeyError(ValueError):
    """Two rows share a physical key, or a row has no usable one."""


def physical_key(species, wavelength_air_A, excitation_potential_eV) -> str:
    """The exact, unrounded identity of a transition.

    ⚠️ `repr` and not `format`: `repr(float)` is the shortest representation that round
    trips, so it neither invents precision the value does not have nor discards precision
    it does. `f"{x:.3f}"` would do both, in different rows of the same file.
    """
    try:
        w = float(wavelength_air_A)
        e = float(excitation_potential_eV)
    except (TypeError, ValueError) as exc:
        raise PhysicalKeyError(
            f"cannot build a physical key from species={species!r} "
            f"wavelength={wavelength_air_A!r} ep={excitation_potential_eV!r}") from exc
    if w != w or e != e:                      # NaN
        raise PhysicalKeyError(
            f"species={species!r} has a NaN wavelength or excitation potential; a row "
            f"with no physical key cannot be given a stable identity")
    return "|".join((str(species).strip(), repr(w), repr(e)))


def physical_id(species, wavelength_air_A, excitation_potential_eV) -> str:
    """A stable id for one transition. Same physics in, same id out, forever."""
    k = physical_key(species, wavelength_air_A, excitation_potential_eV)
    return PREFIX + hashlib.sha256(k.encode("utf-8")).hexdigest()[:_DIGEST_CHARS]


def add_physical_ids(df, *, species_col="species", wave_col="wavelength_air_A",
                     ep_col="excitation_potential_eV", out_col="physical_id"):
    """Stamp `out_col` onto a canonical_gf-shaped frame. Refuses on collision."""
    ids = [physical_id(r[species_col], r[wave_col], r[ep_col])
           for _, r in df.iterrows()]
    n_dup = len(ids) - len(set(ids))
    if n_dup:
        from collections import Counter
        worst = [i for i, c in Counter(ids).most_common(3) if c > 1]
        raise PhysicalKeyError(
            f"{n_dup} row(s) share a physical id — e.g. {worst}. That is either a "
            f"duplicate transition or a key too coarse to separate two real ones. "
            f"REFUSING rather than keeping the first: a silent tiebreak here is how a "
            f"join starts lying (RYA-1033).")
    out = df.copy()
    out[out_col] = ids
    return out
