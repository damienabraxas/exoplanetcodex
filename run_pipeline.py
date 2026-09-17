"""
run_pipeline.py
===============
Single, honest entry point for the Exoplanet Codex abundance pipeline.

This is a THIN DRIVER (RYA-541): it only orchestrates the REAL stage run()s in
order and surfaces their own loud guards. It does NOT generate data, does NOT
hold star parameters, and does NOT paper over unimplemented stages.

The real, validated sequence (the path that produces the banked solar numbers):

    spectra_normalize.run(star)   → continuum-normalize the FITS spectrum
    lines_fit.run(star)           → measure equivalent widths
   [params_stellar.run(star)]     → CONDITIONAL: only for stars that must SOLVE a
                                     spectroscopic-equilibrium parameter (Teff/logg/ξ);
                                     pinned calibrators (solar/Procyon/αCen) skip it
    abundances_derive.run(star)   → EW → A(X) abundances (params pinned from stars.yaml)
    uncertainty_stack.run(star)   → Type A + Type B uncertainty budget
    ratios_interpret.run(star)    → [X/Fe], Mg/Si, C/O, science outputs

There is NO synthetic-acquire step: the real path reads FITS directly inside
spectra_normalize. Linelist loading is internal to abundances_derive (via
data/linelists/loader.py), so there is no separate lines_load stage in the driver.

Star parameters come ONLY from config/stars.yaml via get_star_params() — the
driver never hardcodes them (RYA-292/298).

RYA-1222 adds a SECOND, independent mode on the same entry point. `--element X`
does not run the whole-star chain above; it delegates to `pipeline.run_matrix`,
which expands (star, element) into the full (band x holding x ion x engine) matrix
and drives the RYA-767 resolver over every cell. The two modes share nothing but
this file: the whole-star path is unchanged and is what runs when `--element` is
absent, so nothing that worked before this ticket behaves differently.

Usage:
    python run_pipeline.py --list-stars
    python run_pipeline.py --star solar
    python run_pipeline.py --star solar --validate-only   # stop after the params stage

    # RYA-1222 -- the (band x holding x ion x engine) matrix for one element
    python run_pipeline.py --star solar --element Si --dry-run
    python run_pipeline.py --star solar --element Si
    python run_pipeline.py --star solar --element Fe --band VIS --engine ts-lte

Naming convention — pipeline scripts use the subject_action pattern:
    spectra_*  lines_*  params_*  abundances_*  uncertainty_*  ratios_*
"""

from __future__ import annotations

import argparse
import sys

from config.constants import get_star_params

# The REAL stages the driver orchestrates. (spectra_acquire / lines_load are NOT
# imported here: acquire's run() is a demo inspector only, and linelist loading is
# internal to abundances_derive — see RYA-541.)
#
# 🔴 RYA-1222 — imported INSIDE the whole-star path, not at module load. They used
# to be top-level, which meant every mode of this driver paid for every stage: on
# Python 3.9 `pipeline/abundances_derive.py:835` evaluates `str | None` at class-body
# time and raises TypeError, so `run_pipeline.py --help` and the new --element matrix
# (which loads NONE of these six) both died on an unrelated stage. One broken stage
# taking down a mode that never calls it is not thin-driver behaviour. The whole-star
# sequence itself is untouched: same six stages, same order, same guards — they are
# just resolved when that path runs instead of when the file is read.
def _whole_star_stages():
    from pipeline.spectra_normalize import run as normalize
    from pipeline.lines_fit         import run as fit_lines
    from pipeline.params_stellar    import run as solve_params
    from pipeline.abundances_derive import run as derive
    from pipeline.uncertainty_stack import run as uncertainty
    from pipeline.ratios_interpret  import run as interpret
    return normalize, fit_lines, solve_params, derive, uncertainty, interpret

# ── The star ladder ───────────────────────────────────────────────────────────
# The observing/analysis ladder, in order. Every id MUST resolve via
# get_star_params (config/stars.yaml is the single source of truth). We validate
# that HERE so a stale or unbacked id is reported at startup rather than crashing
# mid-run — this is what killed the old RYA-42 map (it listed hd89307 / gl581 /
# a bare HD id that no longer resolve). Ids without a stars.yaml record are shown
# as "pending" and are NOT runnable (no fabricated parameters).
_PROGRAM_LADDER = [
    'solar', 'procyon', 'alpha_cen_a', 'alpha_cen_b', 'tau_boo', '55cnc_a',
]
STARS: list[str] = []
UNAVAILABLE: list[str] = []
for _sid in _PROGRAM_LADDER:
    try:
        get_star_params(_sid)
        STARS.append(_sid)
    except KeyError:
        UNAVAILABLE.append(_sid)

# Spectroscopic-equilibrium parameters. A star needs the params_stellar SOLVER
# only when its stars.yaml `solve` set contains one of these; [Fe/H] is always
# an output of abundances_derive itself, not something params_stellar solves.
_EQUILIBRIUM_PARAMS = {'teff', 'logg', 'xi'}


def _needs_param_solve(star_id: str) -> bool:
    """True iff the star must SOLVE a spectroscopic-equilibrium parameter (so the
    params_stellar stage applies). Pinned calibrators (solve = [feh] only) → False."""
    solve = set(get_star_params(star_id).get('solve', []))
    return bool(solve & _EQUILIBRIUM_PARAMS)


def _stage(name: str, fn, star_id: str, ticket: str | None = None):
    """Run a stage, converting an unimplemented-stage NotImplementedError into a
    clear, deliberate STOP (never a bare stub traceback). The stage's own loud
    prerequisite guards (e.g. 'Run spectra_normalize.py first') propagate as-is."""
    print(f"\n▶ {name}")
    try:
        return fn(star_id)
    except NotImplementedError as exc:
        ref = f" — see {ticket}" if ticket else ""
        raise SystemExit(
            f"\nSTOP: pipeline stage '{name}' is not yet implemented{ref}.\n"
            f"run_pipeline.py is a thin driver; it will not fabricate this stage's "
            f"output. Implement the stage, then re-run.\n(underlying: {exc})"
        ) from exc


#: RYA-1222. The matrix-only flags, and the whole-star flags they must not be
#: silently mixed with. `--dry-run` without `--element` would otherwise run the
#: REAL whole-star chain while the operator believed nothing was executing, which
#: is a silent-wrong of the worst kind: it writes.
_MATRIX_ONLY = ('ion', 'band', 'instrument', 'engine', 'dry_run')


def _refuse_matrix_flags_without_element(parser, args) -> None:
    given = [f'--{f.replace("_", "-")}' for f in _MATRIX_ONLY if getattr(args, f, None)]
    if given and not args.element:
        parser.error(
            f"{', '.join(given)} {'is' if len(given) == 1 else 'are'} matrix-mode "
            f"flag{'' if len(given) == 1 else 's'} and only mean something with "
            f"--element. Without it this driver runs the whole-star chain, which would "
            f"have ignored them and executed for real.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Exoplanet Codex — stellar abundance pipeline (thin driver)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--star', choices=STARS,
                        help='Target star id (resolves via config/stars.yaml)')
    parser.add_argument('--validate-only', action='store_true',
                        help='Run only through the stellar-parameter stage, then stop')
    parser.add_argument('--list-stars', action='store_true',
                        help='List runnable star ids and exit')

    # ── RYA-1222: the element matrix ─────────────────────────────────────────
    # Validated against data/config/elements_master.json inside run_matrix, not by
    # an argparse `choices` list here: a second copy of the 27-element vocabulary
    # in this file is exactly the drift RYA-914 sweeps for.
    m = parser.add_argument_group(
        'element matrix (RYA-1222)',
        'Given --element, drive the full (band x holding x ion x engine) matrix for '
        'that element instead of the whole-star chain. Every cell gets a terminal '
        'status and a reason; one cell failing never stops the rest.')
    m.add_argument('--element',
                   help='Element symbol as data/config/elements_master.json spells it '
                        '(e.g. Si, Fe, "Fe II"). Drives the matrix mode.')
    m.add_argument('--ion', action='append', metavar='ION',
                   help='Limit to one ion (repeatable). Default: every ion the graded '
                        'pool holds a line for.')
    m.add_argument('--band', action='append', metavar='BAND',
                   help='Limit to one band (repeatable). Default: all.')
    m.add_argument('--instrument', action='append', metavar='ID',
                   help='Limit to one instrument_id or holding_id (repeatable). '
                        'Default: all.')
    m.add_argument('--engine', action='append', metavar='DECK',
                   help='Limit to one engine deck (repeatable). Default: all live decks '
                        'bound to data/catalog/model_registry.csv.')
    m.add_argument('--dry-run', action='store_true',
                   help='Expand, resolve and report intended statuses; execute nothing.')
    args = parser.parse_args()

    if args.list_stars:
        print("\nRunnable stars (resolve via config/stars.yaml):")
        for sid in STARS:
            p = get_star_params(sid)
            kind = 'SOLVE ' + ','.join(sorted(set(p.get('solve', [])) & _EQUILIBRIUM_PARAMS)) \
                if _needs_param_solve(sid) else 'pinned'
            print(f"  {sid:<14} Teff={p['teff']:.0f}  ({kind})")
        if UNAVAILABLE:
            print("\nPending (in the program ladder but no stars.yaml record yet — not runnable):")
            for sid in UNAVAILABLE:
                print(f"  {sid:<14} (add a config/stars.yaml entry before running)")
        print()
        sys.exit(0)

    _refuse_matrix_flags_without_element(parser, args)

    if not args.star:
        parser.print_help()
        sys.exit(1)

    star_id = args.star

    # ── RYA-1222: matrix mode. Delegates entirely; decides nothing itself. ────
    if args.element:
        from pipeline import run_matrix
        try:
            report = run_matrix.run(
                star_id, args.element, ions=args.ion, bands=args.band,
                instruments=args.instrument, engines=args.engine,
                dry_run=args.dry_run)
        except run_matrix.MatrixError as exc:
            raise SystemExit(f"\nSTOP: the matrix could not be built.\n{exc}") from exc
        # A cell that FAILED is a real failure and the exit code must say so, or a
        # CI job reads a green run off a report full of broken cells. BLOCKED and
        # NOT_READY are honest empties and do NOT fail the run -- they are the
        # documented output this layer exists to produce (RYA-1187).
        sys.exit(1 if report["counts"][run_matrix.FAILED] else 0)

    (normalize, fit_lines, solve_params,
     derive, uncertainty, interpret) = _whole_star_stages()

    print(f"\n{'='*62}")
    print(f"  Exoplanet Codex Pipeline  |  {star_id}")
    print(f"{'='*62}")

    # 1–2. Normalize the spectrum, then measure EWs. (Both carry their own loud
    #      'run the previous stage first' guards if an upstream output is stale.)
    _stage('spectra_normalize', normalize, star_id)
    _stage('lines_fit', fit_lines, star_id)

    # 3. Stellar parameters — CONDITIONAL. Only stars that must solve an equilibrium
    #    parameter (55 Cnc ξ; synthetic_no_logg) invoke params_stellar; pinned
    #    calibrators skip it entirely. While params_stellar is the RYA-537 stub, a
    #    solve-star STOPS here with a clear message (never a silent proceed).
    #    (When RYA-537 lands, feed its solved params into abundances_derive via
    #    stellar_params_override at the call below.)
    if _needs_param_solve(star_id):
        _stage('params_stellar', solve_params, star_id, ticket='RYA-537')
    else:
        print("\n▶ params_stellar — SKIPPED (parameters pinned in stars.yaml)")

    if args.validate_only:
        print(f"\n✓  --validate-only: stopped after the parameter stage for {star_id}.")
        return

    # 4. Abundances (the validated EW → A(X) stage; params pinned from stars.yaml).
    _stage('abundances_derive', derive, star_id)

    # 5–6. Post-processing stages. Both are currently stubs (flagged in RYA-536);
    #      the driver STOPS with a clear message rather than skipping silently.
    _stage('uncertainty_stack', uncertainty, star_id, ticket='RYA-536')
    _stage('ratios_interpret', interpret, star_id, ticket='RYA-536')

    print(f"\n✓  Pipeline complete for {star_id}.")


if __name__ == '__main__':
    main()
