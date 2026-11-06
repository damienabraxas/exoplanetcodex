"""RYA-1315 -- one run request, two modes. The CLI and any later front end (Ryan, 2026-10-31:
"this software will eventually live in an online GUI, where a person could select a band,
an instrument, an engine, an element, and do a run") build the SAME object and call the
SAME `execute`, so a selection can never mean one thing on the command line and another
behind a button.

    full    every cell of the element (the first run; publish prices and files the element)
    select  only the cells matching the chosen ions / bands / instruments / engines / routes /
            pools; `force` re-runs them even where the ledger calls them current, and
            `publish` prices, files and withdraws inside that selection only

`options(star, element)` lists every value a selection may take, read from the matrix
itself (`run_matrix.expand`), so a front end offers only what the orchestrator can run.

Not here yet: a run on a spectrum the USER supplies. Every holding is a registered,
prepared (telluric-corrected, frozen) spectrum with a wired reader; a user upload needs
that registration step first, and is the next piece of this interface.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

MODES = ("full", "select")
#: the selection axes, as `run_matrix.run` names its keyword arguments
AXES = ("ions", "bands", "instruments", "engines", "routes", "pools")


@dataclass(frozen=True)
class RunRequest:
    star: str
    element: str
    mode: str = "full"
    ions: tuple[str, ...] = ()
    bands: tuple[str, ...] = ()
    instruments: tuple[str, ...] = ()       # instrument_id or holding_id
    engines: tuple[str, ...] = ()
    routes: tuple[str, ...] = ()            # 'profile-fit' / 'synthesis'
    pools: tuple[str, ...] = ()             # 'set' / 'reference' / 'codex' / 'deep'
    force: bool = False
    publish: bool = False
    dry_run: bool = False
    jobs: int = 1
    extra: dict = field(default_factory=dict, compare=False)

    @classmethod
    def from_dict(cls, d: dict) -> "RunRequest":
        unknown = set(d) - {f for f in cls.__dataclass_fields__}
        if unknown:
            raise ValueError(f"run request: unknown field(s) {sorted(unknown)}")
        d = dict(d)
        for ax in AXES:
            v = d.get(ax) or ()
            d[ax] = (v,) if isinstance(v, str) else tuple(v)
        return cls(**d)

    @classmethod
    def from_json(cls, path: str | Path) -> "RunRequest":
        return cls.from_dict(json.loads(Path(path).read_text()))

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        return {k: (list(v) if isinstance(v, tuple) else v) for k, v in d.items()}

    def selection(self) -> dict[str, list[str]]:
        return {ax: list(getattr(self, ax)) for ax in AXES if getattr(self, ax)}

    def problems(self, opts: dict | None = None) -> list[str]:
        """Why this request cannot run, or [] -- checked BEFORE anything is dispatched."""
        out = []
        if self.mode not in MODES:
            out.append(f"mode {self.mode!r} is not one of {MODES}")
        if self.mode == "full" and self.selection():
            out.append(f"mode 'full' runs every cell; a selection {self.selection()} needs "
                       f"mode 'select'")
        if self.mode == "full" and self.force:
            out.append("force re-runs work the ledger calls current; it needs mode 'select', "
                       "never the whole element")
        if self.mode == "select" and not self.selection():
            out.append("mode 'select' needs at least one of " + ", ".join(AXES))
        if self.publish and self.dry_run:
            out.append("a dry run builds nothing, so it cannot publish")
        if opts:
            for ax, vals in self.selection().items():
                allowed = set(opts.get(ax, ()))
                if ax == "instruments":
                    allowed |= set(opts.get("holdings", ()))
                bad = [v for v in vals if v not in allowed]
                if bad:
                    out.append(f"{ax}: {bad} not available for {self.star} {self.element} "
                               f"(choose from {sorted(allowed)})")
        return out


def options(star: str, element: str) -> dict:
    """Every value each selection axis can take for (star, element), with the cell count."""
    from pipeline import run_matrix as rm
    ds = rm.expand(star, element)
    pools = {("set" if (d.pool or "").startswith(rm.SET_POOL_PREFIX) else d.pool)
             for d in ds if d.pool}
    return {"star": star, "element": element, "cells": len(ds),
            "ions": sorted({d.ion for d in ds}),
            "bands": sorted({rm._band_of(d, star).name_declared for d in ds}),
            "instruments": sorted({d.instrument for d in ds}),
            "holdings": sorted({d.holding for d in ds}),
            "engines": sorted({d.engine_deck for d in ds}),
            "routes": sorted({d.method for d in ds if d.method}),
            "pools": sorted(pools)}


def execute(req: RunRequest, *, interpreter: str | None = None,
            ispec_dir: str | None = None, check_options: bool = True) -> dict:
    """Run the request: the matrix for its cells, then (if asked) publish those cells.
    Returns {"report": run report, "publish": publish summary or None}."""
    from pipeline import element_publish, run_matrix as rm
    why = req.problems(options(req.star, req.element) if check_options else None)
    if why:
        raise ValueError("run request refused: " + "; ".join(why))
    sel = req.selection()
    report = rm.run(req.star, req.element, ions=sel.get("ions"), bands=sel.get("bands"),
                    instruments=sel.get("instruments"), engines=sel.get("engines"),
                    methods=sel.get("routes"), pools=sel.get("pools"),
                    dry_run=req.dry_run, interpreter=interpreter, ispec_dir=ispec_dir,
                    jobs=req.jobs, prepare_first=True, force=req.force)
    report["request"] = req.to_dict()
    pub = None
    if req.publish and not req.dry_run:
        pub = element_publish.publish(req.star, req.element, report=report,
                                      jobs=max(req.jobs, 1), interpreter=interpreter)
    return {"report": report, "publish": pub}
