#!/usr/bin/env python3
"""RYA-1075 — correct the canonical_gf rows inflated by log10(n_isotopes).

    python3 scripts/rya1075_isotope_inflation_fix.py            # classify + report
    python3 scripts/rya1075_isotope_inflation_fix.py --apply     # write the correction
    python3 scripts/rya1075_isotope_inflation_fix.py --verify    # re-check what was written

THE DEFECT. The GES v6 HFS/ISO list is *form (A)* in RYA-684's terms: isotope-coded
blocks whose components carry the FULL oscillator strength, with the engine applying
``isotopfrac`` afterwards. So each isotope's component set sums to the whole transition
gf, and Sum(f_i) = 1 makes the physics come out right at synthesis time.

``cluster_physical_lines`` is isotope-blind. It groups on species + EP + wavelength gap
— and the components of two isotopes of the SAME transition share the lower level and sit
milli-Angstroms apart, so they land in one cluster. Every consumer that then writes
``log10(sum 10**gf)`` over that cluster gets ``n_isotopes x`` the physical gf.

Demonstrated on Eu II 6645, end to end, not inferred from the offset:

    11 GES components, isotopes 151 (7) and 153 (4)
    isotope 151 sums to  +0.1199        <- the transition total
    isotope 153 sums to  +0.1198        <- the same total, as form (A) requires
    naive sum over all   +0.4208        == canonical_gf's published log_gf, exactly
    gf-weighted centroid  6645.0905 A   == canonical_gf's wavelength_air_A, exactly
    offset               +0.3010        == log10(2)

🔴 THIS IS NOT RYA-684's DEFECT, AND ITS CORRECTION TERM IS DIFFERENT.
RYA-684 measured an ENGINE-side double-application — ``isotopfrac`` applied on top of
already-folded gf — whose offset is ``-log10(sum f_i^2)`` and which lives in the
TSFitPy ``linelist_vald`` files. That module's docstring warns, correctly, not to reuse
log10(2) as a correction. This is a CONSUMER-side aggregation error over a
correctly-formed source, and its offset is ``log10(n_isotopes)`` — a count, with no
dependence on the abundances at all.

**La II settles which is which.** La is 99.911% La-139, so RYA-684's term is
``-log10(sum f_i^2) = +0.0008`` while a count gives ``log10(2) = +0.3010``. Nd separates
them too: +0.7258 against +0.8451. The measured offsets are +0.3010 and +0.8451. It is
the count. Applying RYA-684's term here would have left the defect essentially untouched.

🔴 CORRECTING THE STORED VALUE ALONE WOULD INJECT A LIVE ERROR.
``gf_resolver.apply_to_synth_array`` shifts each cluster by ``canon_total - cur_total``,
and it computes ``cur_total`` with the same isotope-blind sum. Today both numbers are
equally inflated, so the shift is 0.0000 and the delivered components pass through
untouched — which is why RYA-684 correctly found this reached no live value. De-inflate
canonical_gf on its own and that shift becomes ``-log10(n)``: every component would be
scaled DOWN by a factor n and the synthesis would be wrong where it is currently right.
So this ticket also makes that aggregation isotope-aware (``gf_resolver.physical_total``).
The two changes are one change; neither is correct alone.

SELECTION. Per row, from the source — never a blanket species rule:

  1. the row's GES cluster spans more than one isotope;
  2. the per-isotope totals agree (form A — the physics test);
  3. the published ``log_gf`` equals the naive all-component sum (the row was in fact
     built by that sum);
  4. the resulting offset equals ``log10(n_isotopes)``, where n is counted from the
     component set actually present, NOT from the element's catalogued isotope count.
     Those differ: Ba II 4934 codes 5 isotopes in GES while makeabund.f lists 7, and
     the measured offset is log10(5).

Anything failing (2) or (4) is a DIFFERENT problem and is reported, never touched.
Li I 6707 is the built-in positive control: two isotopes whose sets do NOT agree
(spread 0.301), offset +0.3266 != log10(2). It must survive untouched.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import codex_path                                   # noqa: E402
from pipeline.gf_resolver import cluster_physical_lines, species_key      # noqa: E402

CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
OUT = ROOT / "data" / "audit" / "rya1075_isotope_inflation"
SIDECAR = ROOT / "data" / "linelists" / "canonical_gf_isotope_corrections_rya1075.csv"

#: How close a canonical row must sit to its reconstructed cluster to BE that row. These
#: are identity tolerances, not physics: the migrator wrote the gf-weighted centroid and
#: the mean EP straight out of the same clustering, so a true match is exact to rounding.
#: An ambiguous match is refused rather than resolved by proximity (RYA-1034).
MATCH_WL_A = 0.02
MATCH_EP_EV = 0.02

#: The published value must BE the naive sum for the row to have been built by it. The
#: canonical table stores 4 decimals, so this is that quantum with a little room.
PUBLISHED_EQ_NAIVE_DEX = 0.001

#: |offset - log10(n_isotopes)|. Measured separation on the real data is a factor ~1400:
#: every form-A cluster lands within 1.8e-05, and the one non-form-A cluster (Li I 6707)
#: sits at 0.0266. Anything in between does not occur; this threshold is nowhere near a
#: boundary and `report` prints both sides of the gap so that stays checkable.
OFFSET_TOL_DEX = 0.002

#: Per-isotope totals must agree for the source to be form (A). Components print 3-4
#: decimals and a cluster carries up to ~25 of them, so summation rounding alone can reach
#: ~0.01 dex; this is that, doubled. Observed maximum among form-A clusters: 0.0097.
FORM_A_SPREAD_DEX = 0.02

CORRECTION_NOTE = (
    "RYA-1075: de-inflated by log10(n_isotopes). The GES v6 source is RYA-684 form (A) "
    "(isotope-coded, fraction-free), so each isotope's component set already sums to the "
    "full transition gf; the isotope-blind cluster sum that built this row counted it "
    "n_isotopes times. Found by RYA-1070's linemake cross-reference."
)


class FixError(RuntimeError):
    """Loud stop. This script never proceeds on a source it could not read."""


# ── the source of truth for isotope identity ─────────────────────────────────
def parse_isotopfrac(path: Path) -> dict[int, dict[int, float]]:
    """``isotopfrac(Z,A)=f`` straight out of the engine source.

    Used ONLY to cross-check the count taken from the component set and to report where
    the two differ. It is NOT the correction's input: the defect counts the isotopes that
    were summed, and the list codes fewer than nature has for some species.
    """
    if not path.exists():
        raise FixError(f"makeabund.f not found at {path} — isotope catalogue unavailable")
    pat = re.compile(r"^\s*isotopfrac\((\d+),(\d+)\)\s*=\s*([0-9.eEdD+-]+)")
    out: dict[int, dict[int, float]] = {}
    for line in path.read_text(errors="replace").splitlines():
        m = pat.match(line)
        if not m:
            continue
        z, a = int(m.group(1)), int(m.group(2))
        f = float(m.group(3).replace("d", "e").replace("D", "e"))
        if a and f > 0:
            out.setdefault(z, {})[a] = f
    if not out:
        raise FixError(f"parsed zero isotopfrac assignments from {path}")
    return out


def reconstruct_clusters(arr) -> pd.DataFrame:
    """Every multi-isotope physical-line cluster in the GES synth source."""
    iso = arr["spectrum_synthe_isotope"]
    keys = [species_key(arr["element"][i], arr["ion"][i], arr["molecule"][i])
            for i in range(len(arr))]
    wls = arr["wave_A"].astype(float)
    eps = arr["lower_state_eV"].astype(float)
    gf = arr["loggf"].astype(float)

    rows = []
    for cl in cluster_physical_lines(keys, wls, eps):
        isotopes = sorted({int(iso[i]) for i in cl if int(iso[i]) != 0})
        if len(isotopes) < 2:
            continue
        w = 10.0 ** gf[cl]
        per = {a: float(np.log10(np.sum(10.0 ** gf[[i for i in cl if int(iso[i]) == a]])))
               for a in isotopes}
        vals = list(per.values())
        rows.append(dict(
            key=keys[cl[0]], Z=keys[cl[0]][0], ion=keys[cl[0]][1],
            centroid_A=float((wls[cl] * w).sum() / w.sum()),
            ep_mean_eV=float(np.mean(eps[cl])),
            n_components=len(cl), n_isotopes=len(isotopes), isotopes=isotopes,
            naive_total=float(np.log10(w.sum())),
            per_isotope_total=float(np.mean(vals)),
            per_isotope_spread=float(max(vals) - min(vals)),
        ))
    if not rows:
        raise FixError("no multi-isotope clusters found — the source or the clustering "
                       "changed shape; refusing to conclude 'nothing to fix'")
    return pd.DataFrame(rows)


def classify(clusters: pd.DataFrame, canon: pd.DataFrame,
             catalogue: dict[int, dict[int, float]]) -> pd.DataFrame:
    """One row per reconstructed cluster: which canonical row it is, and the verdict."""
    canon = canon[canon.key_z.astype(str).str.isdigit()].copy()
    canon["Z"] = canon.key_z.astype(int)
    out = []
    for c in clusters.itertuples():
        m = canon[(canon.Z == c.Z) & (canon.ion == c.ion)
                  & ((canon.wavelength_air_A - c.centroid_A).abs() <= MATCH_WL_A)
                  & ((canon.excitation_potential_eV - c.ep_mean_eV).abs() <= MATCH_EP_EV)]
        base = dict(centroid_A=c.centroid_A, ep_mean_eV=c.ep_mean_eV,
                    n_components=c.n_components, n_isotopes=c.n_isotopes,
                    isotopes=";".join(str(a) for a in c.isotopes),
                    n_isotopes_catalogued=len(catalogue.get(c.Z, {})),
                    naive_total=round(c.naive_total, 6),
                    per_isotope_total=round(c.per_isotope_total, 6),
                    per_isotope_spread=round(c.per_isotope_spread, 6))
        if len(m) != 1:
            out.append(dict(base, line_id="", species="", verdict="NO_UNIQUE_CANONICAL_ROW",
                            reason=f"{len(m)} canonical rows within "
                                   f"{MATCH_WL_A} A / {MATCH_EP_EV} eV — refusing to guess"))
            continue
        r = m.iloc[0]
        published = float(r.log_gf)
        offset = published - c.per_isotope_total
        log10n = math.log10(c.n_isotopes)
        form_a = c.per_isotope_spread <= FORM_A_SPREAD_DEX
        built_by_naive = abs(published - c.naive_total) <= PUBLISHED_EQ_NAIVE_DEX
        offset_ok = abs(offset - log10n) <= OFFSET_TOL_DEX

        rec = dict(base, line_id=r.line_id,
                   #  RYA-1077: `line_id` is a ROW INDEX and rots when a block of rows is
                   #  replaced -- it moved 1,739 committed references, 25% of them. The
                   #  sidecar is exactly such a reference, so it is keyed on the STABLE
                   #  `physical_id` and carries line_id only as a human convenience.
                   physical_id=(r.physical_id if "physical_id" in canon.columns else ""),
                   species=r.species,
                   wavelength_air_A=float(r.wavelength_air_A),
                   excitation_potential_eV=float(r.excitation_potential_eV),
                   published_log_gf=published,
                   sibling_gf_linelist_vald=(float(r.gf_linelist_vald)
                                             if pd.notna(r.gf_linelist_vald) else np.nan),
                   loggf_reference=r.loggf_reference, gf_tier=r.gf_tier,
                   offset_dex=round(offset, 6), log10_n_isotopes=round(log10n, 6),
                   offset_minus_log10n=round(offset - log10n, 8),
                   form_A=form_a, published_equals_naive=built_by_naive)

        if not form_a:
            rec.update(verdict="NOT_FORM_A",
                       reason=f"per-isotope totals disagree by {c.per_isotope_spread:.4f} dex "
                              f"(> {FORM_A_SPREAD_DEX}); the isotopes do not each carry the "
                              f"full gf, so a count correction is not defined here")
        elif not built_by_naive:
            rec.update(verdict="NOT_BUILT_BY_NAIVE_SUM",
                       reason=f"published {published:+.4f} is not the naive cluster sum "
                              f"{c.naive_total:+.4f}; this row took another column and is "
                              f"already correct or differs for another reason")
        elif not offset_ok:
            rec.update(verdict="OFFSET_NOT_LOG10N",
                       reason=f"offset {offset:+.4f} != log10({c.n_isotopes}) = {log10n:+.4f}; "
                              f"a DIFFERENT problem — not corrected")
        else:
            new = c.per_isotope_total
            sib = rec["sibling_gf_linelist_vald"]
            rec.update(verdict="CORRECT", corrected_log_gf=round(new, 4),
                       correction_term=round(-log10n, 6),
                       sibling_residual=(round(new - sib, 6) if np.isfinite(sib) else np.nan))
        out.append(rec)
    return pd.DataFrame(out)


def load_synth_array():
    """Read the GES synth line list through iSpec, as every consumer of it does."""
    try:
        import ispec                                              # noqa: PLC0415
    except ImportError as exc:                                    # pragma: no cover
        raise FixError(
            "iSpec is not importable — this script reads the GES source the same way the "
            "pipeline does, so it cannot run without it. Set PYTHONPATH to the ispec_src "
            f"checkout (see docs). Original error: {exc}") from exc
    import pipeline.abundances_derive as ad                       # noqa: PLC0415
    path = ad._SYNTH_LINELIST_FILE
    if not Path(path).exists():
        raise FixError(f"GES synth line list not found at {path}")
    return ispec.read_atomic_linelist(path), str(path)


def apply_corrections(rep: pd.DataFrame, canon_path: Path) -> tuple[int, str]:
    """Write the corrected log_gf, per row, and stamp the adjudication status.

    The file is rewritten with pandas but the ONLY cells that change are `log_gf` and
    `adjudication_status` on the corrected rows — `verify()` re-reads the result and
    proves that, column by column, against the pre-image.
    """
    todo = rep[rep.verdict == "CORRECT"]
    if todo.empty:
        raise FixError("nothing qualifies — refusing to rewrite the table for no change")
    df = pd.read_csv(canon_path, dtype=str, keep_default_na=False)
    by_id = {lid: i for i, lid in enumerate(df.line_id)}
    prior: dict[str, str] = {}
    for r in todo.itertuples():
        i = by_id.get(r.line_id)
        if i is None:
            raise FixError(f"{r.line_id} vanished from the table between classify and apply")
        if abs(float(df.at[i, "log_gf"]) - r.published_log_gf) > PUBLISHED_EQ_NAIVE_DEX:
            raise FixError(f"{r.line_id} log_gf moved under us "
                           f"({df.at[i, 'log_gf']} vs {r.published_log_gf}) — STOP")
        #  Two of these rows carry a real adjudication (RYA-354 stamped Eu II 6645 to
        #  LWHS, RYA-466 stamped Cu I 5782 to KR1968). Those adjudications are still
        #  CORRECT — they picked the right paper; the number attached to it was the
        #  inflated one. So the status is restamped to this ticket, and the one it
        #  replaces is preserved in the sidecar rather than lost.
        prior[r.line_id] = str(df.at[i, "adjudication_status"])
        df.at[i, "log_gf"] = f"{r.corrected_log_gf:.4f}"
        df.at[i, "adjudication_status"] = "isotope_rya1075"
    df.to_csv(canon_path, index=False)
    import hashlib
    return len(todo), hashlib.sha256(canon_path.read_bytes()).hexdigest(), prior


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="write the corrections")
    ap.add_argument("--verify", action="store_true", help="re-check an applied correction")
    args = ap.parse_args()

    print(f"RYA-1075 — HFS-per-isotope inflation in canonical_gf\n{'=' * 78}")
    ts_src = Path(str(codex_path("engines.ts_source")))
    catalogue = parse_isotopfrac(ts_src / "makeabund.f")
    print(f"isotope catalogue   makeabund.f — {len(catalogue)} elements")

    arr, synth_path = load_synth_array()
    print(f"GES synth source    {synth_path}\n                    {len(arr)} rows")

    clusters = reconstruct_clusters(arr)
    print(f"multi-isotope       {len(clusters)} physical-line clusters span >1 isotope")

    canon = pd.read_csv(CANONICAL, low_memory=False)
    rep = classify(clusters, canon, catalogue)

    print(f"\nVERDICTS\n{'-' * 78}")
    for v, n in rep.verdict.value_counts().items():
        print(f"  {v:26s} {n:4d}")

    ok = rep[rep.verdict == "CORRECT"]
    #  Re-running after a successful apply is the NORMAL case, not an error: every row
    #  now fails `published == naive` because it carries the corrected value. Say so and
    #  stop, rather than crashing on an empty frame and looking like a real failure.
    if ok.empty:
        already = int((rep.verdict == "NOT_BUILT_BY_NAIVE_SUM").sum())
        print(f"\nNOTHING TO CORRECT — {already} clusters no longer match the naive sum, "
              f"which is what an APPLIED correction looks like. Re-run with --verify to "
              f"check the table against the committed sidecar.")
        if args.apply:
            raise FixError("--apply with nothing to correct; the table is already fixed")
    print(f"\nSELECTOR SEPARATION (this is what makes it specific, not a blanket rule)")
    print(f"  |offset - log10(n)| among CORRECT rows : max {rep.loc[rep.verdict == 'CORRECT', 'offset_minus_log10n'].abs().max():.2e}")
    others = rep[rep.verdict.isin({"NOT_FORM_A", "OFFSET_NOT_LOG10N"})]
    if len(others):
        print(f"  |offset - log10(n)| among rejected     : min {others.offset_minus_log10n.abs().min():.2e}")
        print(f"  -> the two populations are separated by a factor "
              f"{others.offset_minus_log10n.abs().min() / max(rep.loc[rep.verdict == 'CORRECT', 'offset_minus_log10n'].abs().max(), 1e-12):.0f}")

    print(f"\nPER SPECIES\n{'-' * 78}")
    for sp, g in ok.groupby("species"):
        n_iso = int(g.n_isotopes.iloc[0])
        cat = int(g.n_isotopes_catalogued.iloc[0])
        note = "" if n_iso == cat else f"  (source codes {n_iso}; makeabund catalogues {cat})"
        print(f"  {sp:7s} {len(g):3d} rows   n_isotopes={n_iso}  "
              f"correction {-math.log10(n_iso):+.4f} dex{note}")

    sib = (ok.dropna(subset=["sibling_residual"]) if "sibling_residual" in ok
           else ok.iloc[:0])
    print(f"\nCORROBORATION vs the gf_linelist_vald sibling ({len(sib)} of {len(ok)} have one)")
    if len(sib):
        print(f"  |corrected - sibling|  median {sib.sibling_residual.abs().median():.2e}  "
              f"95th {sib.sibling_residual.abs().quantile(0.95):.2e}  "
              f"max {sib.sibling_residual.abs().max():.4f}")
        far = sib[sib.sibling_residual.abs() > 0.01]
        for r in far.itertuples():
            print(f"  ⚠️  {r.line_id} {r.species} {r.wavelength_air_A:.4f}: corrected "
                  f"{r.corrected_log_gf:+.4f} vs sibling "
                  f"{r.sibling_gf_linelist_vald:+.4f} ({r.sibling_residual:+.4f} dex)")
    nosib = ok[ok.sibling_residual.isna()] if "sibling_residual" in ok else ok.iloc[:0]
    if len(nosib):
        print(f"  {len(nosib)} corrected rows have NO sibling column — corrected on the "
              f"source reconstruction alone:")
        for r in nosib.itertuples():
            print(f"      {r.line_id} {r.species} {r.wavelength_air_A:.4f}")

    print(f"\nREJECTED — reported, never touched\n{'-' * 78}")
    for r in rep[rep.verdict != "CORRECT"].itertuples():
        print(f"  {r.verdict:24s} {r.line_id or '(unmatched)':10s} {r.species or '':7s} "
              f"{r.centroid_A:10.4f}\n      {r.reason}")

    OUT.mkdir(parents=True, exist_ok=True)
    #  Never overwrite the record with a post-apply re-run: once corrected, every row
    #  classifies as NOT_BUILT_BY_NAIVE_SUM and writing that would ERASE the evidence of
    #  what was changed and why.
    if not ok.empty:
        rep.to_csv(OUT / "classification.csv", index=False)
    else:
        print("  (classification.csv left as committed — a post-apply re-run must not "
              "overwrite the record of what was corrected)")

    if args.apply:
        n, sha, prior = apply_corrections(rep, CANONICAL)
        cols = ["physical_id", "line_id", "species", "wavelength_air_A",
                "excitation_potential_eV",
                "n_components", "n_isotopes", "isotopes", "published_log_gf",
                "correction_term", "corrected_log_gf", "sibling_gf_linelist_vald",
                "sibling_residual", "naive_total", "per_isotope_total",
                "per_isotope_spread", "loggf_reference", "gf_tier"]
        side = ok[cols].copy()
        side["prior_adjudication_status"] = (
            side.physical_id.map(prior) if "physical_id" in side else side.line_id.map(prior))
        side["correction_note"] = CORRECTION_NOTE
        side["source"] = synth_path
        side.to_csv(SIDECAR, index=False)
        print(f"\nAPPLIED {n} corrections")
        print(f"  canonical_gf.csv sha256 {sha}")
        print(f"  provenance sidecar {SIDECAR.relative_to(ROOT)}")

    if args.verify:
        #  Verify the TABLE against the committed SIDECAR, not against a fresh
        #  classification — after a successful apply the classifier finds nothing, and
        #  "0/0 verified" is not a check, it is a vacuous pass.
        if not SIDECAR.exists():
            raise FixError(f"{SIDECAR} missing — nothing to verify against")
        side = pd.read_csv(SIDECAR)
        cur = pd.read_csv(CANONICAL, low_memory=False)
        kc = "physical_id" if ("physical_id" in cur.columns
                               and "physical_id" in side.columns) else "line_id"
        cur = cur.set_index(kc)
        bad = []
        for r in side.itertuples():
            key = getattr(r, kc)
            if key not in cur.index:
                bad.append((key, "ABSENT", r.corrected_log_gf))
                continue
            got = float(cur.at[key, "log_gf"])
            if abs(got - r.corrected_log_gf) > 1e-9:
                bad.append((key, got, r.corrected_log_gf))
        print(f"\nVERIFY  {len(side) - len(bad)}/{len(side)} sidecar rows carry the "
              f"corrected value in canonical_gf (keyed on {kc})")
        for b in bad:
            print(f"  MISMATCH {b}")
        if bad:
            return 1

    (OUT / "provenance.json").write_text(json.dumps(dict(
        ticket="RYA-1075", synth_source=synth_path,
        isotopfrac_source=str(ts_src / "makeabund.f"),
        multi_isotope_clusters=int(len(clusters)),
        verdicts={k: int(v) for k, v in rep.verdict.value_counts().items()},
        thresholds=dict(match_wl_A=MATCH_WL_A, match_ep_eV=MATCH_EP_EV,
                        published_eq_naive_dex=PUBLISHED_EQ_NAIVE_DEX,
                        offset_tol_dex=OFFSET_TOL_DEX,
                        form_a_spread_dex=FORM_A_SPREAD_DEX),
        correction="log_gf := per-isotope total (== naive - log10(n_isotopes_in_cluster))",
        n_isotopes_source="counted from the isotopes present in the GES component set, "
                          "NOT from the element's catalogued isotope count",
        not_rya684="RYA-684's offset is -log10(sum f_i^2) on an ENGINE-side double "
                   "application; this is log10(n) on a CONSUMER-side aggregation. La II "
                   "separates them: +0.0008 vs +0.3010",
        paired_change="pipeline/gf_resolver.physical_total — correcting the stored value "
                      "without it would turn a cancelling error into a live one",
    ), indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)}/classification.csv + provenance.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
