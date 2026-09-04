"""
pipeline/gerber_nlte.py
=======================
RYA-798 — the adapter that lets the PRODUCTION flux-fit path use the TS-native Gerber
departure decks. Engine-B NLTE had a validated deck (RYA-785) and no route to it.

WHAT WAS MISSING, precisely
---------------------------
iSpec's Turbospectrum wrapper has always accepted NLTE:
`ispec.synth.turbospectrum.generate_spectrum(..., nlte_departure_coefficients=...)` writes
the departure + `nlteinfo` files and drives `bsyn_lu`. Two things kept us out of it:

  * `_synth_flux_at_abund` — the ONE generator both the v1 EW path and the v2 flux fit call
    — had no NLTE parameter, so there was nowhere to hand them in;
  * iSpec's own interpolator reads `input/dep-grid/{El}_nlte_grid_data.h5`, and that store
    is EMPTY on Sirius. Our decks are TS-native (`atom.*` + `auxData_*.dat` + the `.bin`),
    consumed only by `scripts/ts_gerber_gate.py`'s direct `interpol_modeles_nlte` + `bsyn`.

So this module produces, from the TS-native deck, exactly the tuple iSpec builds from its
own HDF5 grids:

    {element_name: (departures, tau, ndep, nk, Z, abundance, atom_model_path)}

⚠️ THE DECK HAS NO ABUNDANCE AXIS — MEASURED, NOT ASSUMED
----------------------------------------------------------
`interpol_modeles_nlte` takes an abundance argument, so it looks as though departures can
follow a trial abundance. They cannot. Running the interpolator at A = 7.36 / 7.46 / 7.56
gives departure files that are **byte-identical except for line 8, the abundance stamp
itself** (444836 bytes each; 1 differing line of 132). The Gerber aux table says why:
`A(X)` is not an axis at all, it is exactly `7.50 + [Fe/H]` across all 15229 nodes. The
grid's axes are Teff / logg / [Fe/H] / [alpha/Fe] / vturb.

Two consequences, both load-bearing:

  1. **The departures are a property of the ATMOSPHERE NODE, not of the trial abundance.**
     A flux fit varies A(X) to minimise chi2; the departure coefficients cannot track it.
     We therefore hold them FIXED across the fit and let only the stamp follow the trial
     value — which is what the deck does anyway, not an approximation we introduced. It is
     still an approximation *of the physics* (more Fe means more line opacity means a
     different radiation field means different departures), and it is second-order over the
     ~0.3 dex a fit explores, so it is STATED on the product rather than buried here.
  2. **The interpolation can be cached.** One interpolator run per (element, node) serves
     every trial abundance, instead of one per chi2 evaluation.

⚠️ SOLAR NODE ONLY, TODAY
--------------------------
`ts_gerber_gate` passes ONE concrete MARCS model eight times as the eight interpolation
corners — the TS idiom for a degenerate box, i.e. "use this model, do not interpolate".
That is correct for a solar gate and is why the departure file is unchanged by the
[Fe/H] argument. Serving another star needs real corner selection, so `for_node()` refuses
any node it cannot prove it has the model for, rather than silently returning solar
departures for Procyon. The eight corner paths are recorded in the departure file's own
footer, so the check is against the file rather than against our intent.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np

from config.constants import codex_path  # RYA-810 path register

ROOT = Path(__file__).resolve().parents[1]

# Resolved through the RYA-810 register, never a literal. These stay STRINGS: they are
# interpolated into f-strings, passed to subprocess and fed to os.path.basename, so a Path
# here would be a behaviour change rather than a refactor.
INTERP = str(codex_path('engines.ts_interpolator'))
GT = str(codex_path('grids.gerber_ts'))
PROV_DIR = ROOT / "data" / "nlte_grids" / "gerber_ts"
MARCS_SOLAR = str(codex_path('grids.marcs_standard')
                  / "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_"
                    "o+0.00_r+0.00_s+0.00.mod")

# element -> (Z, atom, aux, grid). Mirrors ts_gerber_gate.ELEMENTS but carries only what
# the adapter needs; the GATE remains the authority on whether a deck may be used at all.
DECKS = {
    "Fe": dict(Z=26, atom="atom.fe607a",
               aux="auxData_Fe_MARCS_May-07-2021.dat",
               grid="NLTEgrid4TS_Fe_MARCS_May-07-2021.bin"),
    # RYA-1005. ⚠️ Al's grid is NOT shaped like Fe's — see `abundance_axis()`. Fe resolves
    # ONE A(X) at fixed [Fe/H]; Al resolves 31 (4.43-7.43 at the solar node). Adding this
    # entry alone, without the axis handling below, would have served every trial
    # abundance the same cached departures and thrown away a real dimension of the grid
    # while still returning a perfectly good-looking file.
    "Al": dict(Z=13, atom="atom.al_qmh",
               aux="auxData_Al_MARCS_Jul-25-2023.dat",
               grid="NLTEgrid_Al_MARCS_Jul-25-2023.bin"),
    # RYA-821/1019 — Al's <3D> deck. A SECOND deck for the SAME element, differing only
    # in the ATMOSPHERE its departures were solved on. That is the atmosphere axis, and
    # it is why a registry keyed by element alone cannot hold both.
    #
    # 🔴 THE NODE COORDINATES ARE STAGGER'S, NOT MARCS'S. This deck is keyed at
    # Teff 5777 / logg 4.44 (the STAGGER solar member). The MARCS decks are keyed at
    # 5750 / 4.5. Interpolating this deck against MARCS_SOLAR asks for a node that DOES
    # NOT EXIST in it — measured: 0 rows at (5750, 4.5, 0.0), 31 rows at (5777, 4.44, 0.0).
    # So the atmosphere is carried per-deck, never assumed.
    #
    # The A(Al) axis is the same 31 values (4.43–7.43) as the MARCS deck, so the RYA-1005
    # abundance-axis handling applies unchanged — do not re-derive it.
    "Al@mean3D": dict(Z=13, atom="atom.al_qmh",
                      aux="auxData_Al_STAGGERmean3D_Aug-05-2023_marcs_names.txt",
                      grid="NLTEgrid_Al_STAGGERmean3D_Aug-05-2023.bin",
                      atmos=str(ROOT / "data" / "atmospheres" / "stagger_avg3d_rya442"
                                / "sun_avg3d_stagger.mod"),
                      atmosphere_family="<3D> STAGGER (Magic et al. 2013)",
                      node_coords="Teff 5777 / logg 4.44 — STAGGER, NOT MARCS 5750/4.5",
                      # RYA-821: no vendor binary can read a <3D> deck -- see
                      # `read_deck_node` for the traced reason and the verified layout.
                      read_via="direct"),
    # RYA-710 (route found + staged under RYA-1035). Fe's <3D> deck -- the second element
    # to take the RYA-821 direct-read route, and the SIMPLER of the two.
    #
    # 🔴 THE PLAIN AUX, NOT THE `_marcs_names` SIBLING -- THE OPPOSITE OF Al ABOVE, AND
    # DELIBERATELY SO. The vendor's aux zeroes [Fe/H] on the seven Teff=5777 rows, and its
    # `convert_3d_grid_to_marcs_names.py` builds the MARCS-style name FROM that column --
    # so the converted file collapses all seven distinct atmospheres to ONE byte-identical
    # string and the solar node becomes unaddressable. The plain file keeps STAGGER's own
    # names, which still encode the metallicity, so `_parse_aux_text` can referee the
    # column and resolve the Sun to `p5777g44m00`. Measured both ways: `_marcs_names`
    # REFUSES the solar node, plain returns it. ⚠️ Only the plain aux is staged on Sirius,
    # so the wrong one cannot be picked up by accident. See RYA-1035.
    #
    # 🔴 NO ABUNDANCE AXIS, so the v118 machinery Al needed does NOT apply here: A(X) is
    # exactly 7.50 + [Fe/H] across all 189 nodes, one value per atmosphere. The pre-v118
    # hoisting of departures out of the chi2 loop was always correct for Fe and stays
    # correct -- `has_abundance_axis` decides that per deck and answers False for this one.
    #
    # Node coords are STAGGER's (5777 / 4.44), like Al's, so the deck carries its own
    # <3D> atmosphere rather than borrowing MARCS_SOLAR.
    "Fe@mean3D": dict(Z=26, atom="atom.fe607a",
                      aux="auxData_Fe_STAGGERmean3D_May-21-2021.txt",
                      grid="NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin",
                      atmos=str(ROOT / "data" / "atmospheres" / "stagger_avg3d_rya442"
                                / "sun_avg3d_stagger.mod"),
                      atmosphere_family="<3D> STAGGER (Magic et al. 2013)",
                      node_coords="Teff 5777 / logg 4.44 — STAGGER, NOT MARCS 5750/4.5",
                      read_via="direct"),
}


# ─────────────────────────────────────────────────────────────────────────────────────
# RYA-1055 — WHICH IONISATION STAGES A DECK CAN ACTUALLY SERVE, AND THE STATED LIMIT
# ─────────────────────────────────────────────────────────────────────────────────────
# 🔴 THE MEASUREMENT. bsyn applies departures PER LINE and falls back to departure = 1
# for any line whose two levels are not identified in the atom. So whether a deck can do
# NLTE for {element} {ion} is decided by ONE property of the atom file: does it declare a
# bound-bound transition BOTH of whose levels belong to that stage?
#
# `atom.fe607a`, measured by `scripts/rya1055_atom_ion_reach.py` on the staged bytes
# (md5 d08dc8232ed68eec65f9bb6631e82ea8, the file `Fe_gerber2023.prov.json` pins):
#
#     levels        548 Fe I (idx 1-548, 0.000-7.505 eV)
#                    58 Fe II (idx 549-606, 7.950-23.753 eV)
#                     1 Fe III (idx 607, 24.102 eV)
#     bound-bound  12,635 declared, 12,635 parsed, 12,635 Fe I -> Fe I
#                       0 with BOTH levels Fe II
#                       0 with EVEN ONE level Fe II
#     highest level index appearing in any radiative transition : 546  (< 549)
#
# ⇒ **The 58 Fe II levels carry no bound-bound transitions at all.** They are a pure
# ionisation reservoir: the targets of Fe I photoionisation, which in turn ionise to
# Fe III. `atom.fe607a` is an Fe I model atom with an Fe II reservoir attached, not a
# two-stage model.
#
# 🔴 THIS IS A DECK LIMIT, NOT A LINE-LIST GAP, AND THE DIFFERENCE DECIDES WHAT TO DO
# ABOUT IT. RYA-1055's first version measured the LINE LIST (0 of 11 lines labelled in
# our Fe II pool, 854 of 8,870 in the window) and proposed labelling a VALD Fe II list
# against the atom. There is nothing in the atom to label against, so that project could
# not have worked. The labelling gap is a SYMPTOM. Enabling Fe II NLTE needs a DIFFERENT
# MODEL ATOM, one carrying Fe II bound-bound transitions — a fetch-or-build of the same
# shape as RYA-1035's <3D> deck hunt, DEFERRED by Ryan (2026-09-03) to the off-solar /
# Beta program where the Fe I/Fe II balance does the log g work.
#
# ⚠️ AND THE FIELD THAT WOULD HAVE REVEALED IT IS THE ONE NOBODY FILTERS ON.
# `turbospectrum_species` COLLAPSES the ionisation stages — every Fe row carries 26.0 —
# while `element`, `ion`, `spectrum_moog_species` and `width_species` all keep them. A
# Z-only filter reports "49% of Fe labelled": the average of a 99% population and a 10%
# population that have nothing to do with each other.
#
# ⚠️ WHAT THIS LIMIT IS **NOT**. Three neighbouring things stay available and must not be
# swept up by a reader who remembers only the headline:
#   * **Fe II ENGINE-A NLTE is a DIFFERENT MECHANISM and is unaffected.** That leg reads
#     the MPIA/Bergemann per-line delta grid (`data/nlte_grids/Fe_Bergemann_MPIA.csv`,
#     6,400 Fe II rows over 80 lines, 3805.5-6586.7 A), never this atom. Its solar Fe II
#     corrections are real and SMALL BUT NOT NEGLIGIBLE-BY-INSPECTION: measured on the
#     grid at the solar node (Teff 5800, logg 4.3/4.5, [Fe/H] 0) the Fe II deltas run
#     -0.002 to +0.016 dex with a median of +0.000, and 146 of 160 entries sit inside
#     +-0.005 — the tail is the three strong multiplet-42 lines 4923.932 / 4924.921 /
#     5169.033, all at +0.015. (Fe I control at the same node: median +0.011, up to
#     +0.040.) What our live Fe II products ACTUALLY applied is -0.001 to -0.002 dex,
#     recorded per line in each product's own `*_ENGINE-A_lines.csv` (`nlte_delta_dex` /
#     `nlte_source`, RYA-880) -- the run speaking, not a label being read.
#     ⚠️ TWO DRAFTS OF THIS COMMENT WERE WRONG, both by generalising from a neighbouring
#     number. The first said "+-0.002 across the grid", read off the first five rows of a
#     truncated summary -- wrong by 8x at the tail. The second said "+0.001 dex on
#     6147.7341 / 6238.3859 / 6247.5570", which is a real measurement of the RYA-489
#     replication product's 11-line 5256-6456 A Fe II set -- a DIFFERENT POOL from the
#     live band product's 4233.162 / 4303.170 / 4583.829, and of the OPPOSITE SIGN. Quote
#     the artifact the product was actually built from.
#   * **Fe II <3D>-LTE is available** via the Frankenstein LTE leg (model 5,
#     `synth-mean3D-LTE-gerber-stagger`): the <3D>-mean atmosphere is ION-AGNOSTIC, so
#     Fe II synthesised on it in LTE gets a real <3D>-LTE treatment. RYA-1135 owns that
#     run. It is NOT the 3D axis that is missing for Fe II — it is the NLTE one.
#   * **Fe I NLTE is unaffected** — 12,635 transitions is the whole Fe I term system.
#
# 🔴 AND THE CONSEQUENCE FOR THE MANDATORY <3D> PAIR: model 5 (LTE) and model 6 (NLTE)
# MUST come out byte-identical ON Fe II, because every Fe II departure is 1 by
# construction. If they ever DIFFER, an Fe II line has inherited a departure from a
# mis-identified Fe I level — a PHANTOM departure. That is a live detector, and it is
# RYA-1135's standing test, not this module's.
#
# Keyed by ATOM, not by deck key, because reach is a property of the atom file: `Fe` and
# `Fe@mean3D` share `atom.fe607a`, so one entry states the limit for both decks (and so
# for BOTH model-registry NLTE members, 4 and 6) and the two can never drift apart.
#
# Al's atom is NOT listed. It has not been measured, and an unmeasured atom must read as
# UNMEASURED rather than inherit Fe's answer — a missing key returns "unknown", never
# "capable" and never "refused".
ATOM_ION_REACH = {
    "atom.fe607a": {
        # ion stage -> can this atom form a line of that stage in NLTE?
        1: True,      # 548 levels, 12,635 bound-bound
        2: False,     # 58 levels, ZERO bound-bound — an ionisation reservoir
        3: False,     # 1 level, the continuum reservoir
    },
}

#: The measurement behind `ATOM_ION_REACH`, regenerated by
#: `scripts/rya1055_atom_ion_reach.py`. Named so a test can hold the table and the
#: artifact against each other -- a constant and a committed artifact that describe the
#: same fact drift silently otherwise, and the stale side is the one that passes.
ATOM_ION_REACH_ARTIFACT = "data/results/rya1055/atom_ion_reach.json"

#: The one-sentence limit, in the words Ryan ratified (2026-09-03). Stamped beside the
#: deck registration AND into every Fe II product's `science_provenance`, because it is a
#: property of the DECK WE SHIP, not of any one run.
FE_II_NLTE_LIMIT = (
    "Fe II NLTE unavailable on the Gerber deck: atom.fe607a carries 58 Fe II levels but "
    "ZERO Fe II bound-bound transitions (12,635 of 12,635 are Fe I -> Fe I), so bsyn "
    "applies departure = 1 to every Fe II line and no line list can enable it. Reported "
    "LTE. The effect is literature-negligible at solar -- Amarsi's grid gives -0.001 dex "
    "(vturb 0.75) and -0.0008 (1.50), and Lind et al. 2017 section 2.2.4 states 'For the "
    "Sun, NLTE effects on Fe II lines are insignificant' -- and our own MPIA/Bergemann "
    "Fe II deltas agree IN THE MEDIAN: +0.000 dex at the solar node, 146 of 160 grid "
    "entries inside +-0.005 -- though the tail is the three strong multiplet-42 lines "
    "4923.932 / 4924.921 / 5169.033, all at +0.015 (max single grid entry +0.016), so "
    "'negligible' is a median statement and not a bound (Fe I control at the same "
    "node: median +0.011, up to +0.040). A two-stage Fe atom is DEFERRED to the off-solar "
    "programme, where the Fe I/Fe II balance does the log g work (Ryan, 2026-09-03). "
    "Fe II <3D>-LTE IS available (model 5, the <3D>-mean atmosphere is ion-agnostic); "
    "model 5 == model 6 on Fe II by construction, which IS the phantom-departure check. "
    "This does NOT cover ENGINE-A: that leg reads the MPIA/Bergemann per-line grid, not "
    "this atom, and its Fe II corrections are real. RYA-1055."
)


_ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV"}


def atom_ion_reach(atom: str) -> "dict | None":
    """{ion_stage: bool} for an atom file, or None where it has never been measured.

    None is not False. An unmeasured atom is a gap in OUR knowledge; returning False
    would state a physical fact we have not established, and callers branch on exactly
    that difference (see `nlte_ion_capability`).
    """
    return ATOM_ION_REACH.get(atom)


def nlte_ion_capability(element: str, ion) -> tuple:
    """(capable, note) — can the registered deck for {element} form {ion} lines in NLTE?

    `capable` is True / False / None, where None means UNMEASURED (see `atom_ion_reach`).
    `note` is the stated limit when the answer is False, "" otherwise.

    An ACCESSOR rather than a lookup table exposed at the call sites, so the ion-stage
    normalisation happens once and a caller cannot spell the stage a second way. Reach
    hangs off the ATOM, so this is the same answer for `Fe` and `Fe@mean3D` — as it must
    be, since they load the same file.
    """
    # ⚠️ `parse_ion` returns the STAGE AS AN INTEGER (1, 2), not the roman spelling --
    # verified, not assumed; reading it as roman is what made the first draft of this
    # accessor answer None for every input while looking perfectly correct.
    from pipeline.species import parse_ion

    cfg = DECKS.get(element)
    if cfg is None:
        return (None, "")
    reach = atom_ion_reach(cfg["atom"])
    if reach is None:
        return (None, "")
    try:
        stage = int(parse_ion(ion))
    except (KeyError, ValueError, TypeError):
        return (None, "")
    ok = reach.get(stage)
    if ok is None:
        return (None, "")
    if ok:
        return (True, "")
    note = FE_II_NLTE_LIMIT if (element == "Fe" and stage == 2) else (
        f"{element} {_ROMAN.get(stage, stage)} NLTE unavailable on {cfg['atom']}: the "
        f"atom carries no bound-bound transition with both levels in that ionisation "
        f"stage, so bsyn applies departure = 1 to every line of it. RYA-1055.")
    return (False, note)

def deck_atmosphere(element: str) -> str:
    """The atmosphere this deck's departures were solved on.

    Defaults to MARCS_SOLAR so every pre-existing deck behaves byte-identically; a <3D>
    deck overrides it. Applying <3D> departures on a 1D MARCS structure (or vice versa)
    is a physics mismatch that produces a perfectly well-formed file, which is exactly
    the class of error worth a guard rather than a comment.
    """
    if element not in DECKS:
        raise GerberDeckError(f"no deck registered for {element!r}")
    return DECKS[element].get("atmos", MARCS_SOLAR)


class GerberDeckError(RuntimeError):
    """Raised loudly. There is no degraded mode: silent LTE is the failure to avoid."""


_AXIS_CACHE: dict[str, tuple] = {}
_AUX_ROW_CACHE: dict[str, list] = {}


def abundance_axis(element: str) -> tuple:
    """The distinct A(X) values this deck resolves at solar [Fe/H], sorted.

    🔴 RYA-1005 — DECKS ARE NOT ALL SHAPED LIKE Fe's, AND THE DIFFERENCE IS INVISIBLE.
    This module's original design rests on a property of the Fe grid that it measured and
    stated: the departures do not depend on A(X), the aux table has ONE abundance at fixed
    [Fe/H], so one interpolation serves every trial and only the stamp follows the fitted
    value. That is true for Fe and FALSE for Al:

        Fe   15,229 aux rows,   1 distinct A(X) at [Fe/H]=0
        Al  454,466 aux rows,  31 distinct A(X) at [Fe/H]=0  (4.43 - 7.43)

    Measured, not inferred: interpolating Al at A = 6.20 / 6.43 / 6.70 returns THREE
    DIFFERENT departure blocks (sha256 d54fb26c / 58e03eb4 / ef995777), where the same
    test on Fe returns one. So for Al the grid carries real abundance information, and
    caching on (element, teff, logg, feh) alone would serve every trial the first
    interpolation and silently discard a dimension — while still returning a perfectly
    well-formed file, which is what makes it worth a guard rather than a comment.

    Returns a 1-tuple for a single-abundance deck, so callers need no special case.

    🔴 RYA-1035 -- IT READS `_aux_rows`, SO THE [Fe/H] COLUMN IS ALREADY REFEREED BY THE
    MODEL NAME. It selects on "[Fe/H] == 0", which is exactly the column the Fe <3D> aux
    zeroes on seven rows spanning -4.0 to +0.5. Parsing the file a second time here would
    have swept all seven metallicities into the solar abundance axis; had their A(X)
    differed it would have MANUFACTURED an abundance axis for a deck that has none, and
    `has_abundance_axis` gates whether departures may be hoisted out of the chi2 loop.
    """
    if element in _AXIS_CACHE:
        return _AXIS_CACHE[element]
    vals = {round(r["abundance"], 4) for r in _aux_rows(element)
            if abs(r["feh"]) < 1e-9}
    if not vals:
        raise GerberDeckError(
            f"{GT}/{DECKS[element]['aux']}: no rows at [Fe/H]=0 — cannot establish this "
            f"deck's abundance axis, and guessing it would be the silent-LTE failure one "
            f"level up")
    out = tuple(sorted(vals))
    _AXIS_CACHE[element] = out
    return out


def has_abundance_axis(element: str) -> bool:
    """True when the deck resolves more than one A(X) — i.e. the abundance must be part
    of the interpolation and of the cache key."""
    return len(abundance_axis(element)) > 1


def deck_abundance(element: str) -> float:
    """The abundance the deck's departures were computed at, from its provenance record.

    NEVER a fitted or reference value. Fe's grid was computed at **A(Fe) = 7.50**: that is
    what `atom.fe607a` declares on its own second line (`7.50  55.85`, md5-matched to our
    staged copy), what BOTH Fe aux tables encode as A(X) = 7.50 + [Fe/H] exactly (15,229
    MARCS rows and 183 clean ⟨3D⟩ rows), and what Turbospectrum's own MARCS reader
    hardcodes as solar iron (`metal = abund(15) - 7.50`, interpol_modeles_nlte.f:1177).

    🔴 RYA-1035 — THIS RECORD SAID 7.46, AND THAT NUMBER WAS OUR OWN INPUT COMING BACK.
    The evidence for 7.46 was a bsyn message, *"NLTE departure coeff calculated for
    abundance = 7.46 while it is 7.50 here"*, read as a statement about the grid. Traced
    through the vendor source, it is a statement about us:

        interpol_modeles_nlte.f:206   read(*,*) abu_ref        <- OUR stdin
        interpol_modeles_nlte.f:761   write(27,1971) abu_ref   <- verbatim into the file
        read_departure.f              -> abundance_nlte
        bsyn.f:988                    -> prints it back at us

    `abu_ref` is a LABEL the caller supplies; the deck never asserts it. This module's own
    docstring already carried the measurement that proves it — running the interpolator at
    A = 7.36 / 7.46 / 7.56 gives departure files **byte-identical except the stamp**. So
    the record was sourced from an echo of `deck_abundance()`'s own previous value: a loop
    with no external referee in it. ⚠️ 7.46 is also the Asplund solar A(Fe), i.e. exactly
    the number that looks right on arrival — which is why the loop went unchallenged.

    The ⟨3D⟩ direct-read path never had the problem: `read_deck_node` reports the AUX's own
    A(X), so it self-declares 7.50 whatever it is asked for. Only the interpolator path
    echoes the caller, and the disagreement between the two paths *was* the ambiguity.
    """
    # A deck key may carry an atmosphere suffix (`Fe@mean3D`); the provenance record is per
    # ELEMENT, and both Fe decks were solved with the same model atom at the same A(Fe).
    base = element.split("@", 1)[0]
    p = PROV_DIR / f"{base}_gerber2023.prov.json"
    if not p.exists():
        raise GerberDeckError(
            f"no provenance record at {p} — an unregistered deck has not passed the "
            f"RYA-534/785 gate and must not be used as an Engine-B leg")
    d = json.loads(p.read_text())
    try:
        a_prov = float(d["deck_abundance"]["a_sun"])
    except (KeyError, TypeError, ValueError) as e:
        raise GerberDeckError(f"{p} carries no usable deck_abundance.a_sun: {e}") from e

    # 🔴 THE RECORD IS CROSS-EXAMINED BY THE GRID, whenever the grid is reachable. A
    # provenance value that contradicts the deck's own aux is how 7.46 survived: nothing
    # ever compared the two, so a number that came from us was passed BACK to the vendor
    # binary as though it had come from the vendor. `abundance_axis` raises when the aux is
    # Sirius-only and absent, and there the record is genuinely all we have.
    try:
        axis = abundance_axis(element)
    except GerberDeckError:
        return a_prov
    if len(axis) == 1 and abs(axis[0] - a_prov) > 1e-6:
        raise GerberDeckError(
            f"{element}: provenance records A(X) = {a_prov:.2f} but the deck's own aux "
            f"declares {axis[0]:.2f} at [Fe/H] = 0. The AUX is the deck speaking; the "
            f"record is us. Do not paper over this by editing whichever is convenient — "
            f"establish what the model atom header says and fix the record to match "
            f"(RYA-1035).")
    return a_prov


def read_departure_file(path: str | os.PathLike) -> dict:
    """Parse `interpol_modeles_nlte`'s output.

    Layout, established by reading a real file rather than from documentation:

        0-7                 eight `#` comment lines (the interpolation weights)
        8                   abundance A(X) the file is STAMPED with
        9                   ndep   (depth points)
        10                  nk     (model-atom levels; 607 for atom.fe607a)
        11 .. 10+ndep       ndep log-tau values
        next ndep lines     nk departure coefficients each  -> (ndep, nk)
        final 8 lines       the eight MARCS corner model paths (provenance footer)
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    if len(lines) < 12:
        raise GerberDeckError(f"{path}: too short to be a departure file")
    try:
        abundance = float(lines[8])
        ndep = int(lines[9])
        nk = int(lines[10])
    except ValueError as e:
        raise GerberDeckError(f"{path}: unreadable header ({e})") from e

    tau = np.array([float(x) for x in lines[11:11 + ndep]], dtype=float)
    if tau.size != ndep:
        raise GerberDeckError(f"{path}: expected {ndep} tau values, got {tau.size}")

    body = lines[11 + ndep:]
    rows = [r.split() for r in body if len(r.split()) == nk]
    if len(rows) != ndep:
        raise GerberDeckError(
            f"{path}: expected {ndep} departure rows of {nk} values, found {len(rows)}")
    dep = np.array(rows, dtype=float)

    corners = [r.strip() for r in body[len(rows):] if r.strip()]
    return dict(abundance=abundance, ndep=ndep, nk=nk, tau=tau,
                departures=dep, corners=corners)


# ── reading a deck DIRECTLY, without the vendor interpolator (RYA-821) ───────
#
# 🔴 WHY THIS EXISTS. `interpol_modeles_nlte` cannot consume a <3D> deck at all, and the
# reason is structural rather than a bug to report upstream: it reads ONLY native MARCS,
# which requires tau_Rosseland and P_g per depth, and BOTH public <3D> STAGGER archives
# ship the 5-column mul23 TAU5000 form (log tau500, T, n_e, V, v_mic) carrying neither.
# The deck's aux names MARCS files the distribution does not contain -- Gerber's group
# evidently converted <3D> models to MARCS internally and shipped the deck, not the
# models. The sibling `interpol_multi_nlte` returns ALL-ZERO b-values and then corrupts
# the heap (glibc malloc assertion) on the very record this reader parses correctly.
#
# So the vendor binary is not a dependency worth having here. The deck's own layout is
# fully determined, and reading it directly ALSO dissolves the corner-model problem: the
# node comes from the aux table, so no eight MARCS corners are needed to name it.
#
# ⚠️ THIS IS A LOOKUP, NOT AN INTERPOLATION, AND THE DIFFERENCE IS THE POINT. It returns
# the departures stored AT a grid node. `interpol_modeles_nlte` interpolates BETWEEN
# eight corners. For a node the deck actually contains -- which is the solar case, and the
# only case we consume today -- those agree by construction and the lookup is the exact
# answer rather than an approximation of it. For an OFF-NODE star they do not, and this
# function refuses rather than pretending; see `read_deck_node`.

#: Record layout, VERIFIED against the file rather than taken from documentation:
#:   500 bytes  atmosphere id (space-padded)
#:     4 bytes  n_dep   int32
#:     4 bytes  n_lev   int32
#:  n_dep*8     log tau           float64
#:  n_dep*nlev*8  departure coefficients, (n_dep, n_lev)  float64
#: For the Al <3D> deck that is 500 + 4 + 4 + 101*8 + 101*354*8 = 287,348 bytes, and the
#: aux table's own consecutive pointers differ by exactly that (1001 -> 288349), so the
#: layout is confirmed twice by independent evidence.
_DECK_ID_BYTES = 500


def _parse_aux_text(text: str) -> list[dict]:
    """The aux table's rows, with the [Fe/H] column REFEREED BY THE MODEL NAME.

    🔴 RYA-1035 -- THE VENDOR'S Fe <3D> AUX HAS A ZEROED METALLICITY COLUMN, AND IT IS
    ZEROED ON EXACTLY THE ROWS A SOLAR RUN SELECTS. Its seven Teff=5777 members are named
    `p5777g44m00` / `m05` / `m10` / `m20` / `m30` / `m40` / `p05` -- the full metallicity
    axis at solar Teff -- but the file's [Fe/H] column reads +0.00 for ALL SEVEN. The
    other 182 rows are correct (measured: name and column agree exactly, 182/182), so the
    column is wrong and the name is right, not the other way round. The Al <3D> aux is the
    positive control: 0 disagreements in 6345 rows, so this correction cannot move Al.

    Left alone the consequence is not a crash, it is a WRONG STAR: every one of the seven
    ties at the solar node, `read_deck_node` breaks the tie on A(X) -- identical at 7.50
    across all seven -- and takes the first, which is `p5777g44m10`, [Fe/H] = -1.0. The
    RYA-821 record-vs-aux check then fires and says "the pointer is wrong", which is a
    correct refusal with a misleading diagnosis: the pointer is right and the aux is
    wrong, and the true solar record (`p5777g44m00`, 6th of the seven) is unreachable.

    So the name referees the column. The override is recorded per row (`feh_aux`,
    `feh_from_name`) and surfaced in provenance -- a silently corrected input is the same
    class of defect as the one being corrected.

    🔴 AND THE OVERRIDE CONDEMNS THE ROW RATHER THAN REPAIRING IT. [Fe/H] is not the only
    metallicity-dependent field in an aux row: A(X) is one too, and on those same six rows
    it is wrong the same way. The deck's own relation A(X) = 7.50 + [Fe/H] holds EXACTLY
    on all 183 other rows and is violated on exactly the six, which shipped A(X) = 7.50 at
    [Fe/H] = -4.0 .. +0.5. The name can referee [Fe/H] because it encodes it; nothing
    referees A(X). So an overridden row is marked SUSPECT and `read_deck_node` REFUSES it
    -- the override exists to get the row OUT of the candidate set, not to repair it.
    Guessing A(X) from the relation would be inventing vendor data.

    The Sun is unaffected: `p5777g44m00` is the one row of the seven the column happens to
    get right, so it is never overridden and never suspect.
    """
    out = []
    for ln in text.splitlines():
        if ln.lstrip().startswith("#"):
            continue
        f = ln.split()
        if len(f) < 9:
            continue
        try:
            row = dict(id=f[0].strip("'"), teff=float(f[1]), logg=float(f[2]),
                       feh=float(f[3]), abundance=float(f[7]), pointer=int(f[8]),
                       feh_aux=float(f[3]), feh_from_name=False)
        except ValueError:
            continue
        named = _node_from_model_name(row["id"])
        # _NODE_TOL_FEH is defined below, next to the node lookup it also governs: ONE
        # tolerance decides both "is this row the node you asked for" and "does this row
        # disagree with its own name", because they are the same question.
        if named is not None and abs(named[2] - row["feh"]) > _NODE_TOL_FEH:
            row["feh"] = named[2]
            row["feh_from_name"] = True
        out.append(row)
    return out


def aux_metallicity_overrides(element: str) -> list[dict]:
    """Rows whose [Fe/H] column disagreed with their own model name (RYA-1035).

    Empty for a clean deck. Non-empty is a VENDOR defect, not ours -- report it, do not
    hide it. Exposed so a status surface can read the defect out of the deck rather than
    carry a hand-written note about it.
    """
    return [r for r in _aux_rows(element) if r["feh_from_name"]]


def _aux_rows(element: str) -> list[dict]:
    """Every aux row as (id, teff, logg, feh, abundance, pointer). Cached per element."""
    if element in _AUX_ROW_CACHE:
        return _AUX_ROW_CACHE[element]
    if element not in DECKS:
        raise GerberDeckError(f"no deck registered for {element!r}")
    aux = f"{GT}/{DECKS[element]['aux']}"
    if not os.path.exists(aux):
        raise GerberDeckError(f"missing Gerber aux table (Sirius-only): {aux}")
    with open(aux) as fh:
        out = _parse_aux_text(fh.read())
    if not out:
        raise GerberDeckError(f"{aux}: no parseable rows")
    _AUX_ROW_CACHE[element] = out
    return out


def _node_from_model_name(name: str) -> tuple[float, float, float] | None:
    """(Teff, logg, [Fe/H]) parsed from ANY of the THREE naming conventions in use.

        't5777g44m00'                             -> (5777, 4.4, 0.0)   STAGGER, long Teff
        'p50g25m40'                               -> (5000, 2.5, -4.0)  STAGGER, SHORT Teff
        's5000_g+2.5_m1.0_t02_st_z-4.00_a+0.00_…' -> (5000, 2.5, -4.0)  MARCS-style

    ⚠️ In the STAGGER forms `g44` is logg*10 and `m00` is [Fe/H]*10 SIGNED BY A LETTER
    ('m' minus / 'p' plus), so `m00` is -0.0 and `m10` is -1.0 -- NOT 0.0 and 1.0. Reading
    those as unsigned is how a [Fe/H] sign gets lost, which is a defect this project has
    already met once in the STAGGER GRID_STATUS table.

    🔴 RYA-1035 -- IN THE MARCS-STYLE FORM THE METALLICITY IS `z`, AND `m` IS THE MASS.
    The two forms put an `m` field in the same place and mean different things by it:
    STAGGER's `m00` is [Fe/H], MARCS's `m1.0` is the stellar mass and its [Fe/H] lives in
    `z-4.00` further along. Reading the MARCS `m` as metallicity returns 0.0 or 1.0 for
    EVERY MARCS-style name -- correct at the solar node by coincidence (mass 0.0, z+0.00)
    and wrong everywhere else, which is exactly why the previous version passed its test:
    that test pinned the solar EXAMPLE, where the two readings agree. The Al <3D> aux is
    entirely MARCS-style, so this is measured on real data, not hypothesised: 123 of its
    node rows disagree between the two readings.

    🔴 RYA-1035 -- THE SHORT-Teff FORM IS NOT COSMETIC, IT IS THE MAJORITY OF A DECK.
    The Fe <3D> deck writes Teff/100 for every model EXCEPT its seven Teff=5777 members:
    182 of its 189 rows are `p50g25m40`-shaped. A parser that only knew the 4-digit form
    returned None for all 182, and `read_deck_node` refuses a record it cannot identify --
    so Fe <3D> would have been unusable at every node except the solar one. Two digits
    are unambiguous against four because STAGGER Teff is always >= 1000 K.

    `s` is the third model-type letter (MARCS spherical, used for logg < 3), alongside
    `p` (plane-parallel) and STAGGER's `t`.
    """
    import re
    n = name.strip()

    # MARCS-style alias -- checked FIRST, because its `_m1.0_` field would otherwise be
    # captured by the STAGGER pattern below and read as a metallicity.
    m = re.match(r"^[tps](\d{4})_g([+-]?\d+(?:\.\d+)?)_.*_z([+-]?\d+(?:\.\d+)?)", n)
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(3))

    # STAGGER's own form: Teff as 4 digits or as Teff/100, logg*10, then [Fe/H]*10 signed
    # by the letter.
    m = re.match(r"^[tp](\d{4}|\d{2})g(\d{2})([mp])(\d{2})$", n)
    if m:
        traw = m.group(1)
        teff = float(traw) if len(traw) == 4 else float(traw) * 100.0
        sign = -1.0 if m.group(3) == "m" else 1.0
        return teff, float(m.group(2)) / 10.0, sign * float(m.group(4)) / 10.0
    return None


#: "Is the requested star AT this grid node?" -- DERIVED from the measured node spacing,
#: not chosen. The Al <3D> deck's Teff nodes are 4000/4500/5000/5500/**5777**/6000/6500/
#: 7000: an otherwise 500 K grid with STAGGER's solar member dropped in, so 5777's nearest
#: neighbour is 223 K away. Our solar star is 5772/4.438 (IAU) against the deck's
#: 5777/4.44 (STAGGER) -- a 5 K CONVENTION difference for the same physical star, not a
#: different one. 25 K gives 5x margin on that and ~9x clearance from the next node, so it
#: cannot reach a neighbour. At Teff 5777 the deck holds exactly ONE logg (4.44).
_NODE_TOL_TEFF = 25.0
_NODE_TOL_LOGG = 0.05
_NODE_TOL_FEH = 0.01


def read_deck_node(element: str, teff: float, logg: float, feh: float,
                   abundance: float, *, tol_teff: float = _NODE_TOL_TEFF,
                   tol_logg: float = _NODE_TOL_LOGG,
                   tol_feh: float = _NODE_TOL_FEH) -> dict:
    """Departures at ONE grid node, read straight out of the deck binary.

    Returns the same shape as `read_departure_file`, so nothing downstream changes.

    🔴 IT REFUSES AN OFF-NODE REQUEST rather than silently returning the nearest node.
    That refusal is the whole safety property: a wrong-node lookup produces a perfectly
    well-formed departure block for a DIFFERENT star, which is exactly the failure this
    module was built to make impossible (and is what feeding a 5750/4.50 MARCS corner to
    the vendor binary would have done -- rc=0 and a plausible file).
    """
    rows = _aux_rows(element)
    node_rows = [r for r in rows
                 if abs(r["teff"] - teff) <= tol_teff
                 and abs(r["logg"] - logg) <= tol_logg
                 and abs(r["feh"] - feh) <= tol_feh]
    if not node_rows:
        near = sorted({(r["teff"], r["logg"], r["feh"]) for r in rows},
                      key=lambda t: (abs(t[0] - teff), abs(t[1] - logg)))[:4]
        raise GerberDeckError(
            f"{element}: deck holds NO node at Teff={teff} logg={logg} [Fe/H]={feh}. "
            f"Nearest nodes: {near}. This deck is keyed at STAGGER coordinates, not "
            f"MARCS -- the solar member is 5777/4.44, NOT 5750/4.50. Refusing to "
            f"substitute a neighbouring node: that returns a well-formed departure block "
            f"for a different star.")

    hit = min(node_rows, key=lambda r: abs(r["abundance"] - abundance))
    if hit["feh_from_name"]:
        # RYA-1035. Reached only if a caller asks for one of the rows whose metallicity
        # column disagreed with its own model name -- which means the row was written
        # wrong, and A(X) (unrefereeable) is wrong on it the same way.
        a_sun = hit["abundance"] - hit["feh_aux"]     # A(X) the row implies at [Fe/H]=0
        raise GerberDeckError(
            f"{element}: aux row {hit['id']!r} is SUSPECT and will not be served. Its "
            f"[Fe/H] column read {hit['feh_aux']:+.2f} while its own model name says "
            f"{hit['feh']:+.2f}. The name referees the metallicity; nothing referees "
            f"A(X), and it is wrong on this row the same way -- it ships "
            f"{hit['abundance']:.2f} where the deck's A(X) = A_sun + [Fe/H] relation "
            f"(exact on every clean row) gives {a_sun + hit['feh']:.2f}. Repairing that "
            f"would be inventing vendor data. RYA-1035.")
    if abs(hit["abundance"] - abundance) > 1e-6 and has_abundance_axis(element):
        axis = abundance_axis(element)
        raise GerberDeckError(
            f"{element}: deck holds no departures at A(X)={abundance:.4f} for this node. "
            f"This deck RESOLVES ABUNDANCE ({len(axis)} values, {axis[0]}-{axis[-1]}), so "
            f"serving the nearest one would discard a real dimension of the grid while "
            f"still returning a well-formed block (RYA-1005). Nearest stored value is "
            f"{hit['abundance']:.4f}.")

    binf = f"{GT}/{DECKS[element]['grid']}"
    if not os.path.exists(binf):
        raise GerberDeckError(f"missing Gerber deck binary (Sirius-only): {binf}")

    # The aux pointer is 1-BASED -- the first record starts at 1001, i.e. after a
    # 1000-byte file header. Off by one here would shift every float by a byte and
    # produce numbers that are finite, plausible and wrong.
    offset = hit["pointer"] - 1
    with open(binf, "rb") as fh:
        fh.seek(offset)
        head = fh.read(_DECK_ID_BYTES + 8)
        if len(head) < _DECK_ID_BYTES + 8:
            raise GerberDeckError(
                f"{binf}: record at pointer {hit['pointer']} is truncated")
        model_id = head[:_DECK_ID_BYTES].decode("ascii", "replace").strip()
        ndep, nk = np.frombuffer(head[_DECK_ID_BYTES:], dtype="<i4", count=2)
        ndep, nk = int(ndep), int(nk)
        if not (0 < ndep < 10_000 and 0 < nk < 10_000):
            raise GerberDeckError(
                f"{binf}: record at pointer {hit['pointer']} reports ndep={ndep} nk={nk}, "
                f"which is not a plausible model -- the pointer or the layout is wrong, "
                f"and a plausible-looking block from a wrong offset is the danger here")
        tau = np.frombuffer(fh.read(ndep * 8), dtype="<f8", count=ndep)
        # 🔴 RYA-1040 — THE BLOCK IS STORED (nlev, ndep), NOT (ndep, nlev).
        # This read `.reshape(ndep, nk)` and every check RYA-821 applied passed anyway,
        # because all of them are orientation-blind: the record SIZE and the aux pointer
        # stride are identical either way (it is the same 101x607 doubles), "100% finite"
        # is true either way, "zero all-zero ROWS" is true either way, and even "median
        # b -> 1.0 at the deepest layer" comes out 1.0000 under BOTH. The layout really
        # was verified twice; the ORIENTATION was never tested, because nothing that had
        # been tested could tell the two apart.
        #
        # What tells them apart is PHYSICS, and three independent lines agree:
        #
        #   1. THE LTE LIMIT. Deep layers must thermalise, so nearly every level should
        #      sit at b~1 there. As-read: 50% of levels within 5% of LTE at depth.
        #      Transposed: 100% -- exactly like the validated 1D deck. And at the SURFACE,
        #      where NLTE bites, as-read gives median b = 0.9905 (i.e. almost no NLTE at
        #      all) against the 1D deck's 0.2026; transposed gives 0.2764.
        #   2. THE ZEROS BECOME AN OBJECT. The record's last 101 doubles are zero, and
        #      101 == ndep. As-read they scatter across one depth over levels 506-606 --
        #      not a thing. Transposed they are LEVEL 606 at every depth: the top level,
        #      unpopulated, which is a thing.
        #   3. THE SYNTHESIS BECOMES PHYSICAL. As-read, the solar Fe I 5240-5260 A window
        #      comes back with 86 pixels ABOVE the continuum and max flux 1.1213 -- a zero
        #      b divides into the source function. Transposed: max 0.9998, zero pixels
        #      above continuum, and an NLTE effect small enough to be Fe I's.
        #
        # ⚠️ THIS AFFECTS `Al@mean3D` TOO, merged since v117. Nothing published moves:
        # no band product has ever been derived through EITHER <3D> deck, which is exactly
        # what v117 and v120 recorded as still owed -- and is why a reader this wrong could
        # sit in main for days looking healthy.
        dep = np.frombuffer(fh.read(ndep * nk * 8), dtype="<f8",
                            count=ndep * nk).reshape(nk, ndep).T

    # The aux row NAMES the model and the record CARRIES its own name -- and they use
    # DIFFERENT CONVENTIONS, which is the whole reason the vendor binary cannot match them
    # and the reason this aux file is called `_marcs_names`:
    #
    #     record : 't5777g44m00'                       <- STAGGER's own naming
    #     aux    : 'p5777_g+4.4_m0.0_t02_st_z+0.00...' <- MARCS-style alias
    #
    # So they are compared on the PHYSICS they encode, never as strings. A string test
    # here is wrong in both directions: it rejects the correct record (which is what it
    # did on the first run of this reader), and it would accept a wrong one that happened
    # to share a prefix.
    node = _node_from_model_name(model_id)
    if node is None:
        raise GerberDeckError(
            f"{binf}: record at pointer {hit['pointer']} carries the unparseable name "
            f"{model_id!r}. Refusing to trust an offset whose record cannot be identified.")
    r_teff, r_logg, r_feh = node
    if (abs(r_teff - hit["teff"]) > tol_teff or abs(r_logg - hit["logg"]) > 0.06
            or abs(r_feh - hit["feh"]) > tol_feh):
        raise GerberDeckError(
            f"{binf}: pointer {hit['pointer']} lands on record {model_id!r} "
            f"(Teff {r_teff} logg {r_logg} [Fe/H] {r_feh}) but the aux row names "
            f"{hit['id']!r} (Teff {hit['teff']} logg {hit['logg']} [Fe/H] {hit['feh']}"
            + (f", column said {hit['feh_aux']} — overridden by the name, RYA-1035"
               if hit["feh_from_name"] else "")
            + "). The pointer is wrong; refusing to return another model's departures.")

    # 🔴 THE DECK STORES LINEAR TAU; `read_departure_file` RETURNS LOG TAU. Converted here
    # so both paths honour ONE contract. Measured on the Al <3D> deck: the stored array is
    # 1e-5 .. 1e5 ascending, i.e. log tau -5.000 .. +5.000 -- exactly this atmosphere's
    # range. Today nothing downstream reads these values (iSpec overwrites the departure
    # tau with the atmosphere's own and `assert_depth_match` checks only the COUNT), which
    # is precisely why the mismatch would have sat here undetected until the first caller
    # that did read them got a number 5 orders of magnitude off.
    tau = np.asarray(tau, float)
    if np.any(tau <= 0):
        raise GerberDeckError(
            f"{binf}: record at pointer {hit['pointer']} has non-positive tau values, so "
            f"it cannot be the linear tau scale this layout expects")
    return dict(abundance=hit["abundance"], ndep=ndep, nk=nk, tau=np.log10(tau),
                departures=np.asarray(dep, float),
                corners=[f"{model_id} (DIRECT deck read, RYA-821 -- node lookup, no "
                         f"interpolation, no vendor binary)"
                         + (f" [RYA-1035: aux [Fe/H] column read {hit['feh_aux']:+.2f} "
                            f"for this row and was overridden to {hit['feh']:+.2f} by the "
                            f"model name -- vendor aux defect, corrected not hidden]"
                            if hit["feh_from_name"] else "")])


def _interpolate(element: str, node: str, teff: float, logg: float, feh: float,
                 abundance: float, workdir: str) -> str:
    cfg = DECKS[element]
    binf, aux = f"{GT}/{cfg['grid']}", f"{GT}/{cfg['aux']}"
    for p in (INTERP, binf, aux):
        if not os.path.exists(p):
            raise GerberDeckError(f"missing Gerber asset (Sirius-only): {p}")
    with open(aux) as fh:
        nrows = sum(1 for ln in fh if not ln.lstrip().startswith("#"))

    os.makedirs(f"{workdir}/work", exist_ok=True)
    os.makedirs(f"{workdir}/Testout", exist_ok=True)
    dep = f"{workdir}/Testout/{node}_{element}_coef.dat"
    if os.path.exists(dep):
        os.remove(dep)          # never inherit a previous run's file (RYA-785 stale guard)

    atmos = deck_atmosphere(element)          # RYA-821: per-deck, not global
    stdin = "\n".join([f"'{atmos}'"] * 8 + [
        f"'{workdir}/Testout/{node}_{element}.interpol'",
        f"'{workdir}/Testout/{node}_{element}.alt'", f"'{dep}'",
        f"'{binf}'", f"'{aux}'", str(nrows),
        f"{teff:.0f}", f"{logg:+.2f}", f"{feh:+.2f}", f"{abundance:.2f}",
        ".false.", ".false.", "'none'", ""])
    r = subprocess.run([INTERP], input=stdin, capture_output=True, text=True,
                       cwd=f"{workdir}/work")
    if not os.path.exists(dep):
        raise GerberDeckError(
            f"interpol_modeles_nlte produced no departure file for {element}:\n"
            f"{r.stdout[-1200:]}\n{r.stderr[-400:]}")
    return dep


_CACHE: dict[tuple, dict] = {}


def departures_at_abundance(element: str, teff: float, logg: float, feh: float,
                            abundance: float) -> dict:
    """Departures at an ARBITRARY A(X), interpolated along the deck's abundance axis.

    🔴 WHY THIS IS NEEDED AT ALL. `read_deck_node` deliberately refuses a non-node
    abundance -- serving the nearest stored value would discard a real dimension of the
    grid (RYA-1005 measured Al's departures genuinely differing across its 31 values).
    But the chi2 loop evaluates CONTINUOUS trial abundances, so a deck with an axis needs
    a defined answer between nodes. Refusing there would make the deck unusable; snapping
    to the nearest node would quantise the fit and bias it toward whichever node the
    optimiser happened to sit near.

    So: LINEAR INTERPOLATION IN A(X) between the two bracketing nodes. That is what the
    vendor `interpol_modeles_nlte` does for the MARCS decks, so the two routes answer the
    same question the same way rather than diverging on method. The axis is UNIFORM at
    0.1000 dex (measured, 31 nodes over 4.43-7.43), so no node is ever far away.

    Degenerates to an EXACT node read when `abundance` lands on one, which the tests pin.
    Extrapolation past either end is refused -- that is `for_node`'s existing rule and it
    is not relaxed here.
    """
    axis = abundance_axis(element)
    a = float(abundance)
    if not (axis[0] - 1e-9 <= a <= axis[-1] + 1e-9):
        raise GerberDeckError(
            f"A({element}) = {a:.3f} is outside this deck's abundance axis "
            f"[{axis[0]:.2f}, {axis[-1]:.2f}] -- extrapolating departures off the grid is "
            f"not a correction, it is an invention.")

    exact = [v for v in axis if abs(v - a) <= 1e-9]
    if exact:
        return read_deck_node(element, teff, logg, feh, exact[0])

    lo = max(v for v in axis if v <= a)
    hi = min(v for v in axis if v >= a)
    d_lo = read_deck_node(element, teff, logg, feh, lo)
    d_hi = read_deck_node(element, teff, logg, feh, hi)
    if d_lo["ndep"] != d_hi["ndep"] or d_lo["nk"] != d_hi["nk"]:
        raise GerberDeckError(
            f"{element}: bracketing nodes A={lo} and A={hi} have different shapes "
            f"({d_lo['ndep']}x{d_lo['nk']} vs {d_hi['ndep']}x{d_hi['nk']}). Interpolating "
            f"across them would pair unrelated depths and levels.")

    w = (a - lo) / (hi - lo)
    dep = (1.0 - w) * d_lo["departures"] + w * d_hi["departures"]
    return dict(abundance=a, ndep=d_lo["ndep"], nk=d_lo["nk"], tau=d_lo["tau"],
                departures=dep,
                corners=[f"{d_lo['corners'][0]} | LINEAR in A(X) between {lo:.2f} and "
                         f"{hi:.2f}, w={w:.4f} (RYA-821)"])


def for_node(element: str, teff: float, logg: float, feh: float,
             node: str = "solar", workdir: str = "/tmp/rya798_gerber",
             abundance: float | None = None) -> dict:
    """Interpolated departures for one atmosphere node.

    Cached. ⚠️ RYA-1005 — the cache key includes the ABUNDANCE for a deck that resolves
    one (`has_abundance_axis`), and does not for a deck that does not. Fe has a single
    A(X), so its key and its behaviour are unchanged; Al has 31 and its departures really
    do differ between them, so keying without the abundance would hand every trial in the
    chi2 loop the first interpolation.
    """
    if element not in DECKS:
        raise GerberDeckError(
            f"no TS-native Gerber deck registered for {element!r} in this adapter "
            f"(registered: {sorted(DECKS)}). Staging a deck is RYA-710; validating it is "
            f"the RYA-534/785 gate. Both must happen before it can be used here.")
    # 🔴 RYA-1049 — `node` IS NOT A DECK SELECTOR AND NEVER WAS.
    # It is a scratch-FILENAME tag: `_interpolate` spends it on
    # `{workdir}/Testout/{node}_{element}_coef.dat` and nothing else. The physics is
    # identical whatever is passed, and on a direct-read deck it is not referenced at all.
    #
    # But the NAME reads like "which node of the grid", so the natural way to reach for a
    # ⟨3D⟩ deck is `for_node("Fe", ..., node="Fe@mean3D")` — which is ACCEPTED, runs
    # clean, and returns the 1D MARCS deck, byte-identical to plain `for_node("Fe", ...)`.
    # The deck key belongs in the FIRST positional argument.
    #
    # That silently served the wrong deck to a ⟨3D⟩ diagnostic, whose numbers were then
    # published as having cleared the ⟨3D⟩ deck. The conclusion happened to survive
    # re-measurement, but the evidence was from a different object and nothing said so.
    # ⚠️ The two decks agree in every orientation-blind respect — both (ndep, 607), both
    # A(X) 7.50, both atom.fe607a, both b→1 at depth, both b<1 at the surface. Only ndep
    # (56 vs 101) and the tau range discriminate. A check that passes on both candidates
    # cannot tell them apart, which is exactly how RYA-821's "verified twice" missed the
    # transpose.
    if node in DECKS:
        raise GerberDeckError(
            f"`node={node!r}` names a registered DECK, but `node` is a scratch-filename "
            f"tag and selects nothing — you would silently get {element!r}'s deck. Pass "
            f"the deck as the FIRST argument instead: for_node({node!r}, teff, logg, feh).")
    axis = abundance_axis(element)
    if has_abundance_axis(element):
        if abundance is None:
            raise GerberDeckError(
                f"{element}'s deck resolves {len(axis)} abundances "
                f"({axis[0]:.2f}-{axis[-1]:.2f}) and its departures DIFFER between them "
                f"(RYA-1005), so an abundance is required. Passing none would silently "
                f"pin the whole chi2 loop to one arbitrary value.")
        if not (axis[0] - 1e-9 <= abundance <= axis[-1] + 1e-9):
            raise GerberDeckError(
                f"A({element}) = {abundance:.3f} is outside this deck's abundance axis "
                f"[{axis[0]:.2f}, {axis[-1]:.2f}] — extrapolating departures off the grid "
                f"is not a correction, it is an invention.")
    key = (element, round(teff, 1), round(logg, 3), round(feh, 3),
           round(float(abundance), 3) if has_abundance_axis(element) else None)
    if key in _CACHE:
        return _CACHE[key]

    # SOLAR NODE ONLY — see the docstring. Refuse rather than quietly hand back solar
    # departures for another star.
    if not (abs(teff - 5750.0) <= 60.0 and abs(logg - 4.50) <= 0.10
            and abs(feh - 0.0) <= 0.10):
        raise GerberDeckError(
            f"the Gerber departure path is solar-node-only today: it passes ONE MARCS "
            f"model ({os.path.basename(deck_atmosphere(element))}) eight times as the eight "
            f"interpolation corners, which is a degenerate box. Asked for "
            f"teff={teff:.0f} logg={logg:.2f} feh={feh:+.2f}. Serving another star needs "
            f"real corner selection — do that rather than accepting solar departures.")

    # A single-abundance deck interpolates at ITS OWN value (Fe: 7.50 — RYA-1035 corrected
    # this from 7.46, which was our own input echoing back through bsyn's message; see
    # `deck_abundance`). An axis deck interpolates at the value actually being synthesised.
    # ⚠️ This feeds the interpolator's `abu_ref`, which is a LABEL: the departures do not
    # move with it (measured), and `as_ispec_tuple` stamps the TRIAL abundance into what
    # bsyn finally reads — which is why bsyn's abundance STOP does not fire during a fit,
    # and why it MUST be left alone rather than downgraded to a warning.
    a_deck = float(abundance) if has_abundance_axis(element) else deck_abundance(element)

    # RYA-821 -- a <3D> deck is READ DIRECTLY; there is no vendor binary that can consume
    # one. See `read_deck_node`. The MARCS decks keep the interpolator path byte-for-byte.
    if DECKS[element].get("read_via") == "direct":
        parsed = (departures_at_abundance(element, teff, logg, feh, a_deck)
                  if has_abundance_axis(element)
                  else read_deck_node(element, teff, logg, feh, a_deck))
        # 🔴 RYA-1040 — THIS BRANCH RETURNED WITHOUT THE THREE FIELDS `as_ispec_tuple`
        # READS, SO NO <3D> DECK COULD EVER REACH A SYNTHESIS. The interpolator branch
        # below sets `atom_path`, `Z` and `deck_abundance` after parsing; the direct
        # reader returned early and set none of them, so the very next call in the chain
        # died with `KeyError: 'Z'`.
        #
        # ⚠️ THAT MEANS Al@mean3D HAS BEEN `T2_CONSUME_WIRED` SINCE v117 AND COULD NOT
        # HAVE SYNTHESISED EITHER. The deck read is genuinely fine -- which is what was
        # tested, and it passed -- but reading a deck is not using one. Nothing caught it
        # because no band product had ever been derived through the route, which is
        # precisely the step v117 and v120 both recorded as STILL OWED. The gap between
        # "the route runs" and "a measurement came out of it" was exactly this wide.
        parsed["atom_path"] = f"{GT}/{DECKS[element]['atom']}"
        parsed["Z"] = DECKS[element]["Z"]
        parsed["deck_abundance"] = a_deck
        # RYA-1049: say WHICH deck this is. The discriminating fields (ndep, tau range)
        # differ between a 1D and a ⟨3D⟩ deck while everything else matches, so a caller
        # that wants to be sure should be able to ask rather than infer.
        parsed["deck_key"] = element
        _CACHE[key] = parsed
        return parsed

    dep_path = _interpolate(element, node, teff, logg, feh, a_deck, workdir)
    parsed = read_departure_file(dep_path)

    # Verify against the file's OWN footer, not against our intent.
    if parsed["corners"] and not all(os.path.basename(deck_atmosphere(element)) in c
                                     for c in parsed["corners"]):
        raise GerberDeckError(
            f"departure file corners are not the requested solar model: "
            f"{parsed['corners'][:2]}")
    parsed["atom_path"] = f"{GT}/{DECKS[element]['atom']}"
    parsed["Z"] = DECKS[element]["Z"]
    parsed["deck_abundance"] = a_deck
    parsed["deck_key"] = element          # RYA-1049 — see the direct-read branch above
    _CACHE[key] = parsed
    return parsed


def as_ispec_tuple(parsed: dict, abundance: float) -> tuple:
    """Pack a parsed departure file into iSpec's tuple.

    ⚠️ ORIENTATION. `interpol_modeles_nlte` writes the block as (ndep, nk) — one line per
    depth, nk values along it — and iSpec's writer wants **(nk, ndep)**. It checks, and
    raises `Shape of depart_coeffs is (56, 607), but expected (607, 56)`, so this is a
    caught error rather than a silent one; the transpose lives here so exactly one place
    knows about the difference.
    """
    return (parsed["departures"].T, parsed["tau"], parsed["ndep"], parsed["nk"],
            parsed["Z"], float(abundance), parsed["atom_path"])


def coefficients(element: str, teff: float, logg: float, feh: float,
                 abundance: float, **kw) -> dict:
    """The dict `ispec.generate_spectrum(nlte_departure_coefficients=...)` expects.

    `abundance` is the value the SYNTHESIS will run at. bsyn checks it against the stamp in
    the departure file and STOPs on a mismatch, so the stamp must follow the trial value.
    The coefficients themselves are node-fixed (module docstring) — this is the stamp
    following the synthesis, not the physics following the abundance.
    """
    kw.setdefault("abundance", float(abundance) if has_abundance_axis(element) else None)
    return {element: as_ispec_tuple(for_node(element, teff, logg, feh, **kw),
                                    abundance)}


def assert_depth_match(parsed: dict, atmosphere) -> None:
    """iSpec OVERWRITES the departure tau with the atmosphere's own (`atmosphere_layers[:, 7]`)
    to dodge Turbospectrum's "tau scales differ" error. That silently assumes the two have
    the same number of layers — if they do not, depth i of the departures is paired with a
    different physical depth of the atmosphere and the correction is quietly wrong.
    """
    n = int(len(atmosphere))
    if n != parsed["ndep"]:
        raise GerberDeckError(
            f"depth mismatch: the departure grid has {parsed['ndep']} layers, the model "
            f"atmosphere has {n}. iSpec overwrites the departure tau with the "
            f"atmosphere's, so these MUST match or the departures are applied at the "
            f"wrong depths.")


def _select_species_rows(linelist, Z: int, ion=None,
                         wave_lo_A: float | None = None,
                         wave_hi_A: float | None = None):
    """Rows for ONE species in ONE window — ionisation stage included when it is known.

    🔴 RYA-1050. `turbospectrum_species` is the ONE field in this list that COLLAPSES the
    ionisation stages: every Fe row carries 26.0, so a filter built on it cannot tell Fe I
    from Fe II. Both callers here used exactly that field, and the consequence is not
    cosmetic. Measured on GESv6 420-920 nm:

        4200-6910 A   Fe 1   8194 of 8243 labelled   99.4%
        4200-6910 A   Fe 2    854 of 8870 labelled    9.6%
                      ALL Fe 9048 of 18366           49.3%   <- what a Z-only filter sees

    Fe I is essentially COMPLETE and Fe II is essentially ABSENT, and the "49%" is an
    average of two populations that have nothing to do with each other. Worse, an Fe II
    NLTE product would pass a Z-only check ON Fe I's LABELS and then synthesise ~90% of its
    own lines in LTE under an NLTE label. That is the case this guard's docstring already
    recorded ("every line in the Fe II 6910-9199 window is unlabelled" -- measured here at
    12.0%, so the substance holds).

    The stage IS available: `element` carries 'Fe 1'/'Fe 2', and so do `ion`,
    `spectrum_moog_species` and `width_species`. Only the field in use lost it.

    ⚠️ `ion=None` keeps the old element-global behaviour deliberately, for callers that
    genuinely mean "this element, any stage".
    """
    species = np.floor(np.asarray(linelist["turbospectrum_species"], dtype=float))
    sel = species == Z
    names = linelist.dtype.names or ()
    if ion is not None and "element" in names:
        from pipeline.species import parse_ion, z_symbol
        want = f"{z_symbol(Z)} {parse_ion(ion)}"
        el = np.asarray([str(x).strip() for x in linelist["element"]])
        sel = sel & (el == want)
    if wave_lo_A is not None and wave_hi_A is not None:
        if "wave_A" in names:
            w = np.asarray(linelist["wave_A"], dtype=float)
        elif "wave_nm" in names:
            w = np.asarray(linelist["wave_nm"], dtype=float) * 10.0
        else:
            w = None
        if w is not None:
            sel = sel & (w >= wave_lo_A) & (w <= wave_hi_A)
    return linelist[sel]


def label_coverage(linelist, Z: int, element: str,
                   wave_lo_A: float, wave_hi_A: float, ion=None) -> tuple[int, int]:
    """(labelled, total) {element} lines in the window — coverage as a FRACTION.

    🔴 RYA-1050. `assert_linelist_supports_nlte` answers "is there at least one", which is
    the only question an assert can usefully ask: bsyn's fallback is departure = 1 PER
    LINE, so partial coverage never raises — it DILUTES. A band at half coverage emits a
    well-formed NLTE product carrying roughly half the effect, and nothing records it.

    Measured on GESv6 420–920 nm, which is our own production list: **9,048 of 18,366 Fe
    lines are labelled in 4200–6910 Å — 49%**. The 67 lines of the graded pool are
    themselves 67/67 labelled, so the TARGETS are fine; it is the surrounding BLEND lines
    that are half in LTE.

    A bare count cannot show that. 9,048 looks reassuring until it is set beside 18,366.

    ⚠️ NO FLOOR IS IMPOSED HERE and none should be without deriving one (RYA-161). This
    reports; the reader judges.
    """
    n_lab = 0
    try:
        rows = _select_species_rows(linelist, Z, ion, wave_lo_A, wave_hi_A)
        if not len(rows):
            return (0, 0)
        flagged = rows[rows["nlte"] == "T"]
        n_lab = int(np.sum((flagged["nlte_label_low"] != "none")
                           | (flagged["nlte_label_up"] != "none")))
        return (n_lab, int(len(rows)))
    except Exception:
        # A coverage REPORT must never be the thing that breaks a run; the assert above is
        # what refuses. Returning (n_lab, 0) makes "not measurable" distinguishable from
        # "zero coverage", which 0/0 would not be.
        return (n_lab, 0)


def pool_label_coverage(linelist, Z: int, waves_A, ion=None,
                        tol_A: float = 0.05) -> tuple[int, int]:
    """(labelled, total) over THE LINES ABOUT TO BE FITTED — the only question that decides
    whether a fit is actually NLTE.

    🔴 RYA-1050, THE THIRD AND SHARPEST NARROWING. Element-global asks the wrong question;
    window-local asks a better one; but a run does not synthesise "the window", it
    synthesises a POOL. Measured, and this is why the distinction is not academic:

        Fe II 4200-6910 A window : 854 of 8870 lines labelled  (9.6%)  -> a window check PASSES
        RYA-877's measured Fe II pool :   0 of 11 lines labelled       -> every line is LTE

    So an Fe II NLTE run on that pool passes a window check on 854 lines it will never
    touch, and then synthesises ALL ELEVEN of the lines it does touch in LTE, under an
    NLTE label. Nothing raises, nothing looks wrong, and the product is a lie.

    ⚠️ This does NOT replace the window check. The window number is what tells you about
    BLENDS -- unlabelled neighbours inside the fitting window contribute to the profile in
    LTE whether or not they are targets. Both numbers are real and they answer different
    questions, so both are reported.
    """
    try:
        rows = _select_species_rows(linelist, Z, ion)
        if not len(rows):
            return (0, 0)
        names = rows.dtype.names or ()
        w = (np.asarray(rows["wave_A"], dtype=float) if "wave_A" in names
             else np.asarray(rows["wave_nm"], dtype=float) * 10.0)
        flag = np.asarray([str(x).strip() for x in rows["nlte"]]) == "T"
        lo = np.asarray([str(x).strip() for x in rows["nlte_label_low"]])
        up = np.asarray([str(x).strip() for x in rows["nlte_label_up"]])
        labelled = flag & ((lo != "none") | (up != "none"))
        n_lab = 0
        n_tot = 0
        for x in np.asarray(waves_A, dtype=float):
            m = np.abs(w - x) <= tol_A
            if not m.any():
                continue
            n_tot += 1
            n_lab += bool(labelled[m].any())
        return (n_lab, n_tot)
    except Exception:
        return (0, 0)


def assert_linelist_supports_nlte(linelist, Z: int, element: str,
                                  wave_lo_A: float | None = None,
                                  wave_hi_A: float | None = None,
                                  ion=None) -> int:
    """iSpec fails SOFT — replicate its skip test and RAISE instead.

    `generate_spectrum` appends an element to `nlte_ignored` and continues **in LTE without
    raising** when no line of that element carries NLTE labels. That is the RYA-764 defect
    (2,644 in-band Fe I lines at `nlte_label_up='none'` running silently in LTE), and it is
    the single worst outcome available here: an LTE synthesis shipped under an NLTE label.
    iSpec does not report `nlte_available` back to the caller, so the check has to happen
    on this side of the call.

    🔴 RYA-1135 — THE ATOM'S REACH IS CHECKED FIRST, AND IT IS A REFUSAL, NOT A FOOTNOTE.
    Everything below keys on "how many of these lines carry NLTE labels", which is the right
    question only when the ATOM has transitions for the ion. For Fe II it does not:
    `atom.fe607a` declares 12,635 bound-bound transitions and not one involves an Fe II
    level (RYA-1055), so bsyn applies departure = 1 to EVERY Fe II line whatever the list
    says about it. Labels are irrelevant to the outcome.

    That made the old placement a live hole rather than a style point. `nlte_ion_capability`
    already returns (False, limit) for Fe II — but it was consulted ONLY inside the
    `n == 0` branch, to phrase the message. A pool with even ONE labelled Fe II line never
    reaches that branch, so the verdict was never applied, and the run proceeded to emit an
    LTE product under an NLTE label. That pool is not hypothetical: 854 of the 8,870 Fe II
    lines in the 4200-6910 A window ARE labelled, and only the accident that RYA-1135's
    9-line DEEPGRADED pool contains none of them made `derive_band_products` refuse.

    So the capability verdict now gates the function, before any counting. A species whose
    atom cannot serve the ion is refused no matter how well its line list is labelled —
    which is the standing phantom-departure detector RYA-1135 asks for, placed where a
    phantom is actually born rather than downstream of it.
    """
    capable, limit = nlte_ion_capability(element, ion) if ion is not None else (None, "")
    if capable is False:
        raise GerberDeckError(
            f"{element} {ion}: NLTE was requested, but this deck's model atom carries NO "
            f"bound-bound transitions for this ionisation stage, so bsyn would apply "
            f"departure = 1 to every line and the product would be LTE under an NLTE label "
            f"(RYA-764/RYA-1050). NO LINE LIST CAN FIX THIS — label coverage is irrelevant "
            f"when there is nothing in the atom to label against. {limit}")

    try:
        # RYA-1050: one selection, shared with `label_coverage`, so the number the guard
        # asserts on and the number the provenance reports cannot disagree.
        # ⚠️ WINDOW-LOCAL AND STAGE-LOCAL, NOT ELEMENT-GLOBAL. bsyn applies departures
        # per LINE and falls back to departure = 1 for any line whose levels are
        # unidentified, so "Fe is labelled somewhere in the list" says nothing about the
        # lines actually being synthesised here. Fe answers the global question yes 15706
        # times while every line in the Fe II 6910-9199 window is unlabelled -- measured
        # at 12.0% there, so the substance of that recorded case holds.
        # `_select_species_rows` applies BOTH narrowings, and is shared with
        # `label_coverage` so the asserted number and the reported number cannot disagree.
        rows = _select_species_rows(linelist, Z, ion, wave_lo_A, wave_hi_A)
        if len(rows) == 0:
            raise GerberDeckError(
                f"no {element} lines in this window's linelist — NLTE cannot engage")
        flagged = rows[rows["nlte"] == "T"]
        n = int(np.sum((flagged["nlte_label_low"] != "none")
                       | (flagged["nlte_label_up"] != "none")))
    except GerberDeckError:
        raise
    except Exception as e:
        raise GerberDeckError(f"cannot read NLTE labels from the linelist: {e}") from e
    if n == 0:
        where = (f" in {wave_lo_A:.1f}-{wave_hi_A:.1f} A"
                 if wave_lo_A is not None and wave_hi_A is not None else "")
        # 🔴 RYA-1055 — THE DIAGNOSIS WAS HARDCODED, AND FOR Fe II IT WAS THE WRONG ONE.
        # This message used to close, unconditionally, with "This is a LINE-LIST coverage
        # gap, not a deck failure — the deck is fine". That sentence is true when the atom
        # HAS the transitions and the list does not point at them, and it is false when
        # the atom has none: `atom.fe607a` declares 12,635 bound-bound transitions and not
        # one of them is Fe II, so an Fe II refusal here is the DECK's limit and no
        # labelling work could ever clear it. A reader who acted on the old sentence would
        # have gone and labelled a VALD Fe II list against an atom with nothing to label
        # against — which is precisely the project RYA-1055 stopped. So the diagnosis is
        # now LOOKED UP rather than asserted.
        capable, limit = nlte_ion_capability(element, ion) if ion is not None \
            else (None, "")
        if capable is False:
            why = (f" 🔴 AND NO LINE LIST CAN FIX THIS. {limit}")
        elif capable is True:
            why = (" This is a LINE-LIST coverage gap, not a deck failure — the atom "
                   "carries transitions for this ionisation stage, these lines have no "
                   "level identification.")
        else:
            # UNMEASURED. Naming which it is would be inventing the answer; say so.
            why = (f" Which side this is — a line-list gap or a deck limit — is NOT "
                   f"established for this species: no ion-stage reach has been measured "
                   f"for its atom (pipeline.gerber_nlte.ATOM_ION_REACH, "
                   f"scripts/rya1055_atom_ion_reach.py). Measure it before labelling "
                   f"anything: RYA-1055 measured Fe II at ZERO bound-bound transitions "
                   f"after a labelling project had already been scoped against it.")
        raise GerberDeckError(
            f"{element}: {len(rows)} lines{where} but NONE carry NLTE level labels. bsyn "
            f"sets departure = 1 for an unidentified level and iSpec drops an unlabelled "
            f"element into `nlte_ignored` — either way the synthesis runs in LTE WITHOUT "
            f"RAISING (RYA-534 Co/Ni, RYA-764). Refusing: an LTE spectrum under an NLTE "
            f"label is worse than no product." + why)
    return n
