"""RYA-1388 -- refit a SUBSET of a pool's lines and reuse every other line's stored fit.

A product is the median of independently fitted lines, so a defect that touches some lines
(RYA-1350: the continuum gate moved 108 Fe HARPS VIS lines) needs only those refitted. This
module supplies the per-line fit result for every OTHER line from the product's own stored
per-line artifact, and `derive_band_products` runs everything after the fit -- the per-line
gates, the aggregate, the budget, every treatment -- through its normal code. A patch
therefore reproduces a full rerun by construction rather than by a second implementation.

WHAT IS REUSED IS THE RAW FIT RESULT, NOT THE ROW. The row carries decisions taken after the
fit (constraint verdict, physical bound, registry, no-data); those are taken again, so a
registry edit since the stored run still applies. `res_from_row` inverts exactly what the
synthesis fit loop reads from `fit_one`'s result: `a_synth`, `status`/`reason` (recovered
from the loop's own `SYNTHESIS: <status> <reason>` text), and the four RYA-847 metrics.

THE PATCH CHECKS ITS OWN PREMISE. Reusing a fit asserts that nothing the fit depends on
changed for that line. The defect class this exists for is the continuum, so the caller
re-places every reused line's continuum (the loader runs for it anyway) and `continuum_drift`
refuses the patch if any reused line's level or method moved: the refit list was incomplete.

Atomic data is never read from the patch list (Ryan, 2026-11-15: "we will want our values
for things like gf"); it selects lines, and each must be a line of the pool being patched.
"""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from pipeline.line_match import MATCH_TOL_A, match

SYNTH_PREFIX = "SYNTHESIS: "
#: the continuum is placed by the same code on the same flux; any difference is a change
CONTINUUM_REL_TOL = 1e-9


class LinePatchError(SystemExit):
    """The patch cannot be applied as asked. SystemExit, like derive's other refusals."""


def read_refit_list(path: str | Path) -> list[dict]:
    """The lines to refit: `wavelength_air_A` required, `species` / `ep_eV` optional.
    Any other column (gf included) is ignored by design."""
    with Path(path).open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "wavelength_air_A" not in rows[0]:
        raise LinePatchError(f"refit list {path}: needs a `wavelength_air_A` column and at "
                             f"least one row")
    out = []
    for r in rows:
        out.append({"species": (r.get("species") or "").strip(),
                    "wavelength_air_A": float(r["wavelength_air_A"]),
                    "ep_eV": float(r["ep_eV"]) if (r.get("ep_eV") or "").strip() else None})
    return out


def _num(v):
    """A stored float back to what `fit_one` returned: '' (written from None) -> NaN, which
    `as_float_or_none` maps back to None exactly as it did the first time."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def res_from_row(row: dict, obs_source: str) -> dict:
    """The `fit_one` result the synthesis loop consumed for this stored row."""
    seg = next((s for s in (row.get("excluded_reason") or "").split(" | ")
                if s.startswith(SYNTH_PREFIX)), None)
    status, reason = "ok", ""
    if seg is not None:
        status, _, reason = seg[len(SYNTH_PREFIX):].partition(" ")
    return {"a_synth": _num(row.get("abundance")), "status": status, "reason": reason,
            "obs_source": obs_source,
            "sigma_A": _num(row.get("sigma_A")),
            "frac_rise_weaker": _num(row.get("frac_rise_weaker")),
            "edge_distance_dex": _num(row.get("edge_distance_dex")),
            "red_chi2": _num(row.get("red_chi2"))}


def continuum_applied(method: str) -> bool:
    """Whether the recorded continuum was divided out: the loop writes `<method>; <reason>`
    and the reason of an applied level is `applied`."""
    return bool(method) and method.rsplit("; ", 1)[-1].startswith("applied")


def continuum_drift(row: dict, level, method: str) -> str:
    """'' when the re-placed continuum divides the flux exactly as the stored one did, else
    what moved. The LEVEL and whether it was APPLIED decide the fitted flux; the method's
    prose does not (RYA-1350 appended "sag corroborated: ..." to unchanged levels)."""
    old_l, new_l = _num(row.get("continuum_level")), _num(level)
    same_l = ((np.isnan(old_l) and np.isnan(new_l))
              or (np.isfinite(old_l) and np.isfinite(new_l)
                  and abs(new_l - old_l) <= CONTINUUM_REL_TOL * max(1.0, abs(old_l))))
    old_a = continuum_applied(row.get("continuum_method") or "")
    new_a = continuum_applied(method or "")
    if same_l and old_a == new_a:
        return ""
    return (f"level {row.get('continuum_level') or '-'} -> "
            f"{level if level is not None else '-'}"
            + ("" if old_a == new_a else f"; applied {old_a} -> {new_a}"))


@dataclass
class LinePatch:
    """One derive run's patch: which candidate lines are refitted, where the stored fits
    of the others live, and what happened per treatment (written as `<stem>_patch.json`)."""
    refit: list[dict]
    prior_dir: Path
    stem: str
    species: str
    aliases: tuple = ()
    refit_wl: set = field(default_factory=set)
    ledger: dict = field(default_factory=dict)
    drift: list = field(default_factory=list)

    @classmethod
    def from_args(cls, refit_csv, prior_dir, stem: str, species: str,
                  aliases: tuple = ()) -> "LinePatch | None":
        if not refit_csv and not prior_dir:
            return None
        if not (refit_csv and prior_dir):
            raise LinePatchError("--refit-lines and --reuse-fits-from go together: a patch "
                                 "needs both the lines to refit and the stored fits of the rest")
        return cls(read_refit_list(refit_csv), Path(prior_dir), stem, species,
                   (species, *aliases))

    def bind(self, cand_wl, cand_ep) -> None:
        """Tie every refit row to exactly one candidate line of THIS pool (EP when the list
        gives it), or refuse. A row of another species is not this run's business."""
        mine = [r for r in self.refit if not r["species"] or r["species"] in self.aliases]
        if not mine:
            return
        ep = [r["ep_eV"] for r in mine]
        use_ep = all(e is not None for e in ep)
        m = match([r["wavelength_air_A"] for r in mine], list(cand_wl),
                  want_ep=ep if use_ep else None, src_ep=list(cand_ep) if use_ep else None,
                  require_ep=use_ep)
        bad = [mine[i]["wavelength_air_A"] for i in np.flatnonzero(~m.resolved)]
        bad += list(m.ambiguous)
        if bad:
            raise LinePatchError(
                f"refit list: {len(bad)} {self.species} line(s) are not lines of this pool "
                f"(no candidate within {MATCH_TOL_A} A, or ambiguous): "
                + ", ".join(f"{float(w):.3f}" for w in bad[:20])
                + ". Pick lines from the pool's own list; a mistyped wavelength is refused.")
        self.refit_wl = {round(float(cand_wl[i]), 3) for i in m.index[m.resolved]}

    def _prior_rows(self, treatment: str) -> list[dict] | None:
        p = self.prior_dir / f"{self.stem}_{treatment}_lines.csv"
        if not p.exists():
            return None
        if treatment not in self.ledger:
            self.ledger[treatment] = {"prior": str(p),
                                      "prior_sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                                      "reused": 0, "refit": [], "refit_no_prior": []}
        with p.open(newline="") as fh:
            return list(csv.DictReader(fh))

    def prior(self, treatment: str, wave_A: float, ep_eV: float) -> dict | None:
        """The stored row to reuse for this line, or None when it must be fitted."""
        rows = self._prior_rows(treatment)
        led = self.ledger.setdefault(treatment, {"prior": None, "reused": 0, "refit": [],
                                                 "refit_no_prior": []})
        if round(float(wave_A), 3) in self.refit_wl:
            led["refit"].append(round(float(wave_A), 3))
            return None
        if rows:
            m = match([wave_A], [float(r["wavelength_air_A"]) for r in rows],
                      want_ep=[ep_eV], src_ep=[_num(r.get("ep_eV")) for r in rows],
                      require_ep=True)    # a lone row of another level is not this line
            if m.n_resolved:
                led["reused"] += 1
                return rows[int(m.index[0])]
        led["refit_no_prior"].append(round(float(wave_A), 3))
        return None

    def note_drift(self, treatment: str, wave_A: float, what: str) -> None:
        self.drift.append((treatment, round(float(wave_A), 3), what))

    def assert_no_drift(self) -> None:
        if self.drift:
            raise LinePatchError(
                f"patch refused: {len(self.drift)} reused line(s) no longer place the "
                f"continuum they were fitted on, so their stored fits are not current. Add "
                f"them to the refit list: "
                + "; ".join(f"{t} {w:.3f} ({d})" for t, w, d in self.drift[:20]))

    def write(self, out_dir: Path) -> Path:
        p = Path(out_dir) / f"{self.stem}_patch.json"
        p.write_text(json.dumps({"stem": self.stem, "species": self.species,
                                 "prior_dir": str(self.prior_dir),
                                 "refit_requested": len(self.refit_wl),
                                 "treatments": self.ledger}, indent=1) + "\n")
        return p
