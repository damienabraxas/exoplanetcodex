"""Telluric-applied determination at INTAKE, from the file headers (RYA-806).

THE AXIS THIS OWNS
------------------
`telluric_basis` (RYA-786, `instrument_catalog.csv`) answers **does this band NEED
correction** — a property of the instrument. This module answers the orthogonal
question: **has correction ALREADY been APPLIED to THIS product** — a property of the
downloaded product level, which varies dataset to dataset. The same exposures can exist
as a raw cr2res IDP and as a molecfit-corrected ADP; the instrument axis cannot tell
them apart, and collapsing the two is what RYA-806 exists to prevent.

WHY THIS IS A HEADER QUESTION AND NOT A PREFERENCE
--------------------------------------------------
RYA-805: RYA-370 asserted "no telluric correction" for the CRIRES+ Vesta set, RYA-373's
spec repeated it, and **neither ever showed the keyword** — four tickets inherited the
claim. It happened to be true. That was not knowable without looking, so intake looks.

WHAT COUNTS AS EVIDENCE
-----------------------
A pipeline records a telluric step in one of three places, and this module checks all
three before concluding anything:

  1. **The recipe chain** — `ESO PRO REC*n* ID` walked to whatever depth it goes. A
     molecfit / telluric / `corr_tell` recipe anywhere in the chain is direct evidence
     of APPLIED. A chain that ends at extraction is evidence of NOT-APPLIED.
  2. **A transmission extension** — a `TRANS` / `RECON` / `TELL` HDU. ⚠️ Its *presence*
     is not proof: a transmission array that is **all 1.0 means the correction was NOT
     applied** (the model exists, nothing was divided by it). The values decide.
  3. **A second flux column** — `FLUX_TELL` beside `FLUX`, the corrected-vs-raw pair.

⚠️ TWO TRAPS, both of which fake an answer rather than erroring:

  * **`HIERARCH` (RYA-791).** astropy strips the prefix, so looking up
    `HIERARCH ESO PRO CATG` returns empty and **manufactures an absence**. When the
    verdict *is* an absence, that bug is indistinguishable from the finding. Every
    lookup here uses the bare `ESO ...` form.
  * **Vocabulary false positives (RYA-805).** A naive telluric regex over the CRIRES+
    Vesta headers returns **162 hits, all 162 spurious** — `ESO DET DEV1 BOARD*n* TRANS`
    detector shift registers, the `ESO OBS AMBI TRANS` sky-transparency *constraint*,
    and FITS boilerplate. Counting them would have inverted the verdict, so they are
    filtered explicitly and the filter is tested.

THE THIRD VALUE IS LOAD-BEARING
-------------------------------
`unknown` is a real answer, not a failure to try. A product whose headers do not speak
to telluric state must be refused, never defaulted — defaulting to `applied` fabricates a
correction (the RYA-786 forbidden move) and defaulting to `not-applied` silently sends a
corrected product through a second correction. Consumers refuse on `unknown`; see
`telluric_policy.gate_holding`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

APPLIED = "applied"
NOT_APPLIED = "not-applied"
UNKNOWN = "unknown"
VALUES = (APPLIED, NOT_APPLIED, UNKNOWN)

# A telluric-correction step, wherever a pipeline records one.
_TELLURIC_RECIPE = re.compile(
    r"molecfit|telluric|corr_tell|calctrans|skycorr|wapiti|atmospheric_correction",
    re.IGNORECASE,
)
# Cards that MATCH the vocabulary but are not telluric evidence (RYA-805).
_FALSE_POSITIVE_KEY = re.compile(
    r"^(COMMENT|HISTORY)$"
    r"|ESO DET .*TRANS"          # detector board shift registers
    r"|ESO OBS AMBI TRANS",      # the OB sky-transparency CONSTRAINT, not a correction
    re.IGNORECASE,
)
_TRANSMISSION_HDU = re.compile(r"TRANS|RECON|TELL|MTRANS", re.IGNORECASE)
_TELLURIC_COLUMN = re.compile(r"FLUX_TELL|TELL|MTRANS|TRANSMISSION|RECON", re.IGNORECASE)
# A FLUX column carrying the correction (NIRPS/ESPRESSO DRS: FLUX_TELL_EL, ...) vs one
# that does not (FLUX, FLUX_EL, FLUX_CAL). Anchored on FLUX so ERR_*/QUAL_* siblings do
# not count as flux, and `_TELL` is what separates the two.
_CORRECTED_FLUX_COLUMN = re.compile(r"^FLUX_TELL(_|$)", re.IGNORECASE)
_RAW_FLUX_COLUMN = re.compile(r"^FLUX(_(?!TELL)[A-Z0-9_]+)?$", re.IGNORECASE)


@dataclass
class TelluricEvidence:
    """What the headers actually said. `value` is the verdict; the rest is the citation."""

    value: str
    path: str = ""
    recipe_chain: tuple[str, ...] = ()
    telluric_recipes: tuple[str, ...] = ()
    transmission_hdus: tuple[str, ...] = ()
    telluric_columns: tuple[str, ...] = ()
    transmission_all_unity: bool | None = None
    pro_catg: str = ""
    # ⚠️ Set when the correction lives in a SEPARATE column while an UNCORRECTED flux
    # column survives in the same file (NIRPS `S1D_FINAL_A`). `applied` is then true of
    # the product but NOT of its default column, so a consumer that reads `FLUX` gets
    # uncorrected flux out of a holding the registry calls corrected. Callers must be
    # told which column to read; see `required_column`.
    corrected_columns: tuple[str, ...] = ()
    uncorrected_columns: tuple[str, ...] = ()
    reasons: list[str] = field(default_factory=list)

    @property
    def applied_in_column_only(self) -> bool:
        return bool(self.corrected_columns and self.uncorrected_columns)

    @property
    def required_column(self) -> str:
        """The column a consumer MUST read to get corrected flux, or '' if in place."""
        return self.corrected_columns[0] if self.applied_in_column_only else ""

    def citation(self) -> str:
        """One line fit to paste into the holdings registry `notes`."""
        bits = [f"telluric_applied={self.value}"]
        if self.pro_catg:
            bits.append(f"PRO CATG={self.pro_catg}")
        if self.recipe_chain:
            bits.append(f"PRO REC chain={'|'.join(self.recipe_chain)} "
                        f"({len(self.recipe_chain)} step"
                        f"{'s' if len(self.recipe_chain) != 1 else ''})")
        if self.telluric_recipes:
            bits.append(f"telluric recipe(s)={'|'.join(self.telluric_recipes)}")
        if self.transmission_hdus:
            bits.append(f"transmission HDU={'|'.join(self.transmission_hdus)}"
                        + (" ALL 1.0 => NOT applied"
                           if self.transmission_all_unity else ""))
        if self.applied_in_column_only:
            bits.append(
                f"CORRECTION IS IN A COLUMN, NOT IN PLACE: read {self.required_column}; "
                f"{'/'.join(self.uncorrected_columns)} in the SAME file "
                f"{'are' if len(self.uncorrected_columns) > 1 else 'is'} UNCORRECTED")
        elif self.telluric_columns:
            bits.append(f"telluric column={'|'.join(self.telluric_columns)}")
        return "; ".join(bits + self.reasons)


def _bare(key: str) -> str:
    """astropy strips `HIERARCH `; normalise so a caller's prefixed key still matches."""
    return key[9:] if key.upper().startswith("HIERARCH ") else key


def from_headers(path: str | Path) -> TelluricEvidence:
    """Determine `telluric_applied` for ONE product file, from its headers alone.

    Returns the verdict WITH its evidence. Never guesses: a file that says nothing about
    telluric state comes back `unknown`, which consumers refuse.
    """
    from astropy.io import fits

    path = str(path)
    ev = TelluricEvidence(value=UNKNOWN, path=path)
    with fits.open(path) as hdul:
        h0 = hdul[0].header
        ev.pro_catg = str(h0.get("ESO PRO CATG", h0.get("PRODCATG", ""))).strip()

        # ── 1. the recipe chain, walked to whatever depth it goes ────────────
        chain, n = [], 1
        while f"ESO PRO REC{n} ID" in h0:
            chain.append(str(h0[f"ESO PRO REC{n} ID"]).strip())
            n += 1
        ev.recipe_chain = tuple(chain)
        ev.telluric_recipes = tuple(r for r in chain if _TELLURIC_RECIPE.search(r))

        # any other card naming a telluric step (some pipelines log it outside REC*)
        vocab_hits = []
        for hdu in hdul:
            for k, v in hdu.header.items():
                if _FALSE_POSITIVE_KEY.search(_bare(k)):
                    continue
                if _TELLURIC_RECIPE.search(f"{_bare(k)} = {v}"):
                    vocab_hits.append(f"{_bare(k)}={str(v).strip()}")

        # ── 2. a transmission extension (presence is NOT proof; values decide) ─
        trans_hdus, all_unity = [], None
        for hdu in hdul:
            if hdu.name and _TRANSMISSION_HDU.search(hdu.name):
                trans_hdus.append(hdu.name)
                data = getattr(hdu, "data", None)
                if data is not None and getattr(data, "size", 0):
                    try:
                        import numpy as np
                        arr = np.asarray(data, dtype=float).ravel()
                        arr = arr[np.isfinite(arr)]
                        if arr.size:
                            unity = bool(np.allclose(arr, 1.0))
                            all_unity = unity if all_unity is None else (all_unity and unity)
                    except (TypeError, ValueError):
                        pass          # structured/non-numeric ext: presence only
        ev.transmission_hdus = tuple(trans_hdus)
        ev.transmission_all_unity = all_unity

        # ── 3. a corrected-vs-raw column pair ───────────────────────────────
        cols, all_cols = [], []
        for hdu in hdul[1:]:
            if getattr(hdu, "columns", None) is not None:
                names = [c.name for c in hdu.columns]
                all_cols += names
                cols += [c for c in names if _TELLURIC_COLUMN.search(c)]
        ev.telluric_columns = tuple(cols)

        # Does an UNCORRECTED flux column survive alongside the corrected one? If so the
        # product is corrected but its DEFAULT column is not, and saying only "applied"
        # would hand a consumer uncorrected flux under a corrected label.
        corrected_flux = [c for c in all_cols
                          if _CORRECTED_FLUX_COLUMN.match(c)]
        if corrected_flux:
            raw_flux = [c for c in all_cols
                        if _RAW_FLUX_COLUMN.match(c) and c not in corrected_flux]
            ev.corrected_columns = tuple(corrected_flux)
            ev.uncorrected_columns = tuple(raw_flux)

    # ── the verdict ─────────────────────────────────────────────────────────
    if ev.transmission_all_unity is True:
        ev.value = NOT_APPLIED
        ev.reasons.append(
            "a transmission extension is present but is ALL 1.0 — the model was computed "
            "and nothing was divided by it, so the correction was NOT applied")
        return ev

    if ev.telluric_recipes or vocab_hits or ev.telluric_columns or (
            ev.transmission_hdus and ev.transmission_all_unity is False):
        ev.value = APPLIED
        if ev.telluric_recipes:
            ev.reasons.append(
                f"a telluric recipe is IN the PRO REC chain: {'|'.join(ev.telluric_recipes)}")
        elif vocab_hits:
            ev.reasons.append(f"telluric keyword(s) present: {'; '.join(vocab_hits[:3])}")
        elif ev.telluric_columns:
            ev.reasons.append(
                f"a corrected-flux/transmission column is present: "
                f"{'|'.join(ev.telluric_columns)}")
        else:
            ev.reasons.append(
                f"a transmission extension is present and is NOT all-unity: "
                f"{'|'.join(ev.transmission_hdus)}")
        return ev

    if ev.recipe_chain:
        # The chain is legible and contains no telluric step. That is a POSITIVE
        # finding, not a failed search — we can see every recipe that ran.
        ev.value = NOT_APPLIED
        ev.reasons.append(
            f"the PRO REC chain is legible and complete ({len(ev.recipe_chain)} step"
            f"{'s' if len(ev.recipe_chain) != 1 else ''}: {'|'.join(ev.recipe_chain)}) "
            f"and contains NO telluric step; no transmission HDU and no FLUX_TELL column")
        return ev

    ev.value = UNKNOWN
    ev.reasons.append(
        "the headers carry no PRO REC chain, no transmission extension and no telluric "
        "column, so they do not speak to telluric state. UNKNOWN is refused by "
        "telluric_policy.gate_holding — it is never defaulted either way (RYA-806)")
    return ev


def from_many(paths) -> tuple[str, list[TelluricEvidence]]:
    """Fold a holding's files into ONE value. A holding is only as certain as its worst
    file: any `unknown` makes the holding `unknown`, and a mixed applied/not-applied
    holding is `unknown` too, because a single flag cannot describe it honestly."""
    evs = [from_headers(p) for p in paths]
    if not evs:
        return UNKNOWN, evs
    vals = {e.value for e in evs}
    if vals == {APPLIED}:
        return APPLIED, evs
    if vals == {NOT_APPLIED}:
        return NOT_APPLIED, evs
    return UNKNOWN, evs
