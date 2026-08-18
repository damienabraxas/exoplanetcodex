#!/usr/bin/env python3
"""
RYA-870 — the replication-grade per-line data product (RYA-489 Section 6).

One generated CSV per (element x star): every constant needed to reproduce a per-line
abundance in the same stack, the per-line A(X), and the problem children in the SAME file
with a coded reason. It is the download the appendix links to.

🔴 THIS IS A JOIN, NOT A MEASUREMENT. Every number is projected from an artifact that is
already committed. Nothing here fits, infers, or fills. That is what makes it cheap to
regenerate and what makes it safe to publish: if a number here is wrong, its source is
wrong, and the source is named in the header.

THE FIVE SOURCES, and the rule each one carries:

  linelist_<star>.csv        atomic constants — EP, damping rad/stark/vdW
  canonical_gf.csv           🔴 gf OVERRIDES the linelist wherever canonical has the line
                             (RYA-834 single source). Never read a duplicated gf literal.
  band products *_lines.csv  measurement + per-line A(X), one row per (line x instrument
                             x engine). ⚠️ ENGINE IS A PER-ROW AXIS (RYA-489 change 1) and
                             engines are NEVER combined (RYA-712).
  problem_children.csv       status + reason_code + reason_note. The vocabulary is
                             IMPORTED from pipeline.problem_children, never retyped
                             (RYA-463); an excluded line is a ROW, never a drop (RYA-844).
  STAR_PARAMS + litscan      header. Cited, never remembered.

⚠️ NO SILENT FALLBACK. A missing required column or artifact RAISES and names the file.
The ONLY legal blank is a measured line absent from problem_children, which means
`status=in_aggregate` and an empty reason — an explicit state, not a gap.

⚠️ WHAT THE ROW COUNT MUST SATISFY. Emitted rows == measured rows, exactly. The product is
an ACCOUNTING of the measured set, so a line that was excluded still appears. `build` fails
if the two disagree, because a per-line product that quietly drops lines is the defect
RYA-844 exists to prevent, wearing a data-product costume.
"""
from __future__ import annotations

import datetime as _dt
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from config import constants
from pipeline import litscan
from pipeline.problem_children import CURATED_CLASSES, AUTO_CLASSES, REGISTRY_CSV

ROOT = Path(__file__).resolve().parents[1]
LINELIST_DIR = ROOT / "data" / "linelists"
CANONICAL_GF = LINELIST_DIR / "canonical_gf.csv"
GOLD_CURRENT = ROOT / "data" / "reference" / "solar" / "elements" / "_current.json"
LAB_GF = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"
PRODUCTS_DIR = ROOT / "data" / "products"

#: Matching tolerances. EP is required alongside wavelength — RYA-780/852 both found that
#: a wavelength-only window returns a high-excitation neighbour as if it were the line.
WAVE_TOL_A = 0.05
EP_TOL_EV = 0.05

#: The status vocabulary of RYA-489 Section 6.4. `in_aggregate` is the absence of a
#: problem-children row; the other three come from the registry's `required_treatment`.
STATUS_IN_AGGREGATE = "in_aggregate"
_TREATMENT_TO_STATUS = {
    "exclude": "excluded",
    "investigate": "flagged_kept",
}

ROW_COLUMNS = [
    # identity
    "element", "ion", "wavelength_air_A", "excitation_potential_eV",
    # atomic constants — the replication payload
    "log_gf", "gf_source", "gf_grade", "gf_sigma_dex",
    "damping_rad", "damping_stark", "damping_vdW", "damping_form", "hfs_isotope_note",
    # measurement
    "instrument", "arm", "method", "ew_mA", "reduced_ew", "red_chi2",
    # the per-line STATISTICAL uncertainty (see the note on _sigma_A below)
    "sigma_A_dex", "sigma_A_basis",
    # engine / model context
    "engine", "scale",
    # result
    "A_X_line",
    # status / provenance
    "status", "reason_code", "reason_note",
]


class PerLineProductError(RuntimeError):
    """A required artifact or column is missing. Never downgraded to a warning."""


@dataclass
class PerLineProduct:
    header: dict
    rows: pd.DataFrame
    accounting: dict = field(default_factory=dict)

    def to_csv(self, path: Path) -> Path:
        """Write the commented header block, then the rows."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="") as fh:
            fh.write("# RYA-489 Section 6 replication-grade per-line data product\n")
            fh.write("# Every value below is PROJECTED from a committed artifact; this\n")
            fh.write("# file measures nothing. Regenerate with the generator named below.\n")
            for k, v in self.header.items():
                fh.write(f"# {k}: {v}\n")
            fh.write("#\n")
            self.rows.to_csv(fh, index=False)
        return path


def _require(df: pd.DataFrame, cols, artifact: str) -> None:
    """Loud-fail on a missing column, naming the artifact (RYA-419/632)."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise PerLineProductError(
            f"{artifact} is missing required column(s) {missing}. This product refuses to "
            f"fill a replication constant with a default — fix the source or widen the "
            f"binding deliberately.")


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True, check=True).stdout.strip()


def _assert_inputs_committed(paths) -> str:
    """RYA-654 pattern: refuse to stamp a product with a SHA that does not describe it.

    A file generated off an uncommitted artifact carries a commit_sha that points at
    something else — the stamp becomes a lie that looks like provenance.
    """
    dirty = _git("status", "--porcelain", "--", *[str(p) for p in paths]).splitlines()
    if dirty:
        raise PerLineProductError(
            "refusing to generate: these inputs are uncommitted, so `commit_sha` would "
            "not describe the file it stamps —\n  " + "\n  ".join(dirty))
    return _git("rev-parse", "HEAD")


def _build_header(star: str, element: str, sources: dict, input_commit: str) -> dict:
    if star not in constants.STAR_PARAMS:
        raise PerLineProductError(
            f"{star!r} is not in STAR_PARAMS; known: {sorted(constants.STAR_PARAMS)}")
    p = constants.STAR_PARAMS[star]
    for k in ("teff", "logg", "xi", "source"):
        if k not in p:
            raise PerLineProductError(f"STAR_PARAMS[{star!r}] has no {k!r}")

    gold = json.loads(GOLD_CURRENT.read_text()) if GOLD_CURRENT.exists() else {}
    gold_keys = sorted(k for k in gold if k.split("_")[0] == element)
    rng = litscan.literature_range(element)

    return {
        "star": star,
        "element": element,
        "Teff": p["teff"], "log_g": p["logg"], "Fe_H": p.get("feh_ref"),
        "xi_microturb": p["xi"],
        "params_source": p["source"],
        "reference_stack": "py3.12 + numpy 2.x (RYA-517 reference stack)",
        "damping_source": sources.get("damping_source", "synthesis"),
        "replication_grade": sources.get("replication_grade", "yes"),
        "gold_version": ", ".join(f"{k}={gold[k]}" for k in gold_keys) or "none frozen",
        "commit_sha": input_commit,
        "input_commit": input_commit,
        "litscan_anchor": (
            f"best_external {rng.best_external} sigma_ext {rng.sigma_external}"
            if rng is not None else "no litscan for this element"),
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "generator_script": "scripts/generate_perline_product.py",
        "sources": "; ".join(f"{k}={v}" for k, v in sources.items()),
    }


def _load_linelist(star: str) -> pd.DataFrame:
    f = LINELIST_DIR / f"linelist_{star}.csv"
    if not f.exists():
        raise PerLineProductError(f"no linelist for {star!r}: {f}")
    df = pd.read_csv(f)
    _require(df, ["element", "ion", "wavelength_air_A", "excitation_potential_eV",
                  "damping_rad", "damping_stark", "damping_vdW"], str(f))
    return df


def _load_canonical_gf() -> pd.DataFrame:
    df = pd.read_csv(CANONICAL_GF)
    _require(df, ["species", "wavelength_air_A", "excitation_potential_eV", "log_gf",
                  "loggf_reference", "nist_grade"], str(CANONICAL_GF))
    return df


def _grade_label(nist_letter, gf_source: str) -> str:
    """RYA-711: a grade must NAME ITS SUBJECT.

    🔴 The schema line says "RYA-711 MQ-A/B/C/D, never NIST's bare letters". Those are two
    different subjects and only one of them lives here: `mq_grade` is a MEASUREMENT-quality
    score computed from line_score in abundances_derive, while what canonical_gf carries is
    the ATOMIC-DATA grade of the gf itself. A projection cannot manufacture the former, and
    republishing the latter bare is exactly the collision RYA-711 was opened about.

    So the letter is NAMESPACED rather than dropped: `NIST:B` cannot be mistaken for `MQ-B`
    wherever the column travels, and the information survives. Flagged for RYA-489 to settle
    the wording; this is the reading that loses nothing and confuses nothing.
    """
    letter = "" if nist_letter is None else str(nist_letter).strip()
    if letter and letter.lower() not in ("nan", "none"):
        return f"NIST:{letter}"
    src = (gf_source or "").lower()
    if "kurucz" in src or src.startswith("k") and src[1:3].isdigit():
        return "ungraded_kurucz"
    return "ungraded"


def load_synthesis_damping():
    """The damping constants the SYNTHESIS actually used, from the GES list itself.

    🔴 WHY THIS EXISTS. The first cut of this product published damping from
    `linelist_<star>.csv` because that is what the ticket's source binding named. The
    RYA-870 guard then caught the consequence: for 6094.372 the product published a
    classical log gamma of -7.179 while the synthesis that produced its A(X) used the GES
    list's ABO packed sigma.alpha of 914.276 — a DIFFERENT PHYSICAL QUANTITY in a
    different form — and re-inverting on the published value missed by 0.011 dex.

    A replication-grade row must carry the constants that produced ITS OWN number. So the
    binding moves to the synthesis linelist, and the ticket's own instruction covers it:
    "bind each source-column name to the actual artifact at implementation time".

    ⚠️ THIS MAKES THE GENERATOR SIRIUS-ONLY. The GES list is an iSpec resource, and iSpec
    lives on Sirius. That is a real cost and it is taken deliberately rather than papered
    over with a fallback: a `damping_vdW` whose meaning depended on which machine ran the
    generator would be the silent-fallback defect this module refuses everywhere else.
    """
    from pipeline.abundances_derive import _load_synth_resources   # iSpec import
    ll, _, _ = _load_synth_resources()
    return ll


def _damping_form(value) -> str:
    """ABO packed sigma.alpha, or a classical log gamma — stated, per RYA-489 §6.4.

    The convention is the sign: a packed ABO value is a positive sigma.alpha composite
    (914.276 = sigma 914, alpha 0.276); a classical broadening constant is a negative log
    gamma. Publishing the number without the form makes it unusable, which is why the
    schema asks for the form explicitly.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "unknown"
    return "ABO packed sigma.alpha" if v > 0 else "classical log gamma"


def _load_lab_sigma() -> pd.DataFrame:
    """RYA-850: the pool's OWN cited per-line laboratory sigma, where one exists."""
    if not LAB_GF.exists():
        return pd.DataFrame(columns=["wavelength_air_A", "elo_eV", "e_loggf_dex", "source"])
    df = pd.read_csv(LAB_GF)
    _require(df, ["wavelength_air_A", "elo_eV", "e_loggf_dex", "source"], str(LAB_GF))
    return df


def _load_problem_children() -> pd.DataFrame:
    df = pd.read_csv(REGISTRY_CSV)
    _require(df, ["species", "lambda_or_scope", "problem_class", "required_treatment",
                  "status", "notes"], str(REGISTRY_CSV))
    known = CURATED_CLASSES | AUTO_CLASSES
    unknown = sorted({c for cls in df["problem_class"].astype(str)
                      for c in cls.split("+")} - known)
    if unknown:
        raise PerLineProductError(
            f"{REGISTRY_CSV} carries problem_class values outside the RYA-463 vocabulary: "
            f"{unknown}. The reason code is imported, never retyped — declare it in "
            f"pipeline/problem_children.py or fix the row.")
    return df


def _match(df: pd.DataFrame, wl: float, ep: float | None, species: str | None = None):
    """Wavelength AND excitation potential. Never wavelength alone (RYA-780/852)."""
    m = (df["wavelength_air_A"] - wl).abs() <= WAVE_TOL_A
    if species is not None and "species" in df.columns:
        m &= df["species"].astype(str).str.strip() == species
    if ep is not None and "excitation_potential_eV" in df.columns:
        m &= (df["excitation_potential_eV"] - ep).abs() <= EP_TOL_EV
    hit = df[m]
    return None if hit.empty else hit.iloc[0]


def _disposition(pc: pd.DataFrame, species: str, wl: float, in_aggregate: bool,
                 excluded_reason: str):
    """status / reason_code / reason_note for one line (RYA-463 + 844).

    Two independent sources have to agree here, and where they do not, the DERIVER wins on
    membership and the registry wins on the reason: the registry says what is known about a
    line, the deriver says what actually happened to it in this product.
    """
    row = None
    for _, r in pc.iterrows():
        if str(r["species"]).strip() != species:
            continue
        scope = str(r["lambda_or_scope"])
        head = scope.split()[0].split("(")[0]
        try:
            if abs(float(head) - wl) <= WAVE_TOL_A:
                row = r
                break
        except ValueError:
            continue

    if row is None:
        if in_aggregate:
            return STATUS_IN_AGGREGATE, "", ""
        # excluded by the deriver with no registry row: the reason still has to be stated,
        # so it is carried verbatim from the deriver rather than invented here.
        return "excluded", "DERIVER_EXCLUDED", (excluded_reason or "").strip()

    code = str(row["problem_class"])
    note = str(row["notes"])
    if in_aggregate:
        return "flagged_kept", code, note
    return _TREATMENT_TO_STATUS.get(str(row["required_treatment"]), "excluded"), code, note


def _lab_sigma(lab: pd.DataFrame, wl: float, ep: float | None):
    if lab.empty:
        return ""
    m = (lab["wavelength_air_A"] - wl).abs() <= WAVE_TOL_A
    if ep is not None:
        m &= (lab["elo_eV"] - ep).abs() <= EP_TOL_EV
    hit = lab[m]
    return "" if hit.empty else float(hit.iloc[0]["e_loggf_dex"])


#: RYA-870 — THE TWO PER-LINE UNCERTAINTIES ARE DIFFERENT AXES AND ARE NEVER SUMMED HERE.
#:
#:   sigma_A_dex    STATISTICAL. 1 sigma on THIS line's abundance from the local curvature
#:                  of its own chi2 (`fit_constraint.curvature_sigma`) — the published
#:                  sigma_stat for CNO's N and O (RYA-848), hoisted by RYA-847 so one
#:                  definition serves every synthesis path. It already carries two
#:                  corrections worth not re-learning: chi2 is rescaled to red_chi2 == 1
#:                  first, because Delta-chi2 = 1 is one sigma only for a CALIBRATED chi2
#:                  and ours runs 2.6-1226 (unrescaled, solar O was understated 7.7x); and
#:                  the probe is two-sided and unclamped, because the old one-sided clamped
#:                  probe reported sigma = 0.000 for a fit sitting ON the bracket — the
#:                  tightest possible bar on the worst possible fit.
#:   gf_sigma_dex   SYSTEMATIC. How well this line's oscillator strength is known, from the
#:                  pool's own cited laboratory sigma where one exists (RYA-850).
#:
#: A reader may combine them; this file will not. They answer different questions, the
#: budget combines them at PRODUCT level, and RYA-850 established that a cited sigma
#: REPLACES the generic bound rather than joining it.
#:
#: ⚠️ BLANK ON THE EW ROUTE, AND THAT IS A STATEMENT. Curvature sigma needs a chi2 surface;
#: an EW inversion does not produce one. `sigma_A_basis` says so per row rather than
#: leaving a bare blank nobody can distinguish from "nobody computed it" (RYA-833).
def _sigma_A(row):
    v = row.get("sigma_A")
    if v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() in ("", "nan"):
        return "", ("no chi2 surface on this route — EW inversion, "
                    "see ew_mA for the measured quantity")
    return float(v), "chi2 curvature at the minimum, rescaled to red_chi2==1 (RYA-847/848)"


def _damp(synth_ll, ll, wl: float, ep, synth_col: str, list_col: str):
    """The damping the SYNTHESIS used, falling back to the linelist ONLY in the
    explicitly-labelled non-replication-grade mode."""
    if synth_ll is not None:
        import numpy as _np
        near = _np.abs(synth_ll["wave_A"] - wl) < WAVE_TOL_A
        if ep is not None:
            near = near & (_np.abs(synth_ll["lower_state_eV"] - ep) < EP_TOL_EV)
        idx = _np.where(near)[0]
        if idx.size and synth_col in synth_ll.dtype.names:
            return float(synth_ll[synth_col][idx[0]])
        return None
    return (ll[list_col] if ll is not None else None)


def _arm_of(wl_A: float) -> str:
    """RYA-489 Section 2 boundaries, as used by the band products."""
    if wl_A < 3780:
        return "near-UV"
    if wl_A < 6910:
        return "VIS"
    if wl_A < 9199:
        return "red-optical"
    return "NIR"


def build_perline_product(star: str, element: str, band_products,
                          damping_source: str = "synthesis") -> PerLineProduct:
    band_products = [Path(p) for p in band_products]
    # 🔴 A CONTROL EXPERIMENT IS NOT A PRODUCT. RYA-877 committed before/after per-line
    # files under data/results/rya877/control/ to evidence a same-inputs diff; a plain
    # rglob swept them in and the Fe II lines were counted THREE times (1039 -> 1104 rows).
    # The row-accounting check could not see it — it compares emitted against what was
    # READ, and both were inflated together. Directories that hold a deliberate re-run of
    # the same measurement are excluded by name, and the duplicate guard below is the
    # backstop for the ones I have not thought of.
    EXPERIMENT_DIRS = {"control", "before", "after", "sweep", "_regen"}
    line_files = sorted(
        f for root in band_products for f in root.rglob("*_lines.csv")
        if not (EXPERIMENT_DIRS & set(f.relative_to(root).parts[:-1])))
    if not line_files:
        raise PerLineProductError(
            f"no *_lines.csv under {[str(p) for p in band_products]} — this product is a "
            f"projection of the band products and cannot be built without them")

    linelist = _load_linelist(star)
    canonical = _load_canonical_gf()
    pc = _load_problem_children()
    lab = _load_lab_sigma()

    synth_ll = None
    if damping_source == "synthesis":
        try:
            synth_ll = load_synthesis_damping()
        except Exception as e:                                   # pragma: no cover
            raise PerLineProductError(
                f"cannot read the synthesis linelist ({type(e).__name__}: {e}). The "
                f"damping constants a replication-grade row must publish are the ones the "
                f"SYNTHESIS used, and those live in the GES list, which needs iSpec — run "
                f"this on Sirius under venv312. To emit anyway, pass "
                f"--damping-source linelist, which produces a file that is explicitly NOT "
                f"replication-grade and says so in its header.")

    inputs = [LINELIST_DIR / f"linelist_{star}.csv", CANONICAL_GF, REGISTRY_CSV, *line_files]
    input_commit = _assert_inputs_committed(inputs)

    rows, n_measured = [], 0
    for f in line_files:
        d = pd.read_csv(f)
        _require(d, ["element", "ion", "wavelength_air_A", "instrument", "treatment",
                     "in_aggregate", "abundance"], str(f))
        d = d[d["element"].astype(str).str.strip() == element]
        if d.empty:
            continue
        deck = f.parent.name if f.parent.name in ("ts-lte", "gerber-nlte") else ""
        for _, r in d.iterrows():
            n_measured += 1
            wl = float(r["wavelength_air_A"])
            species = f"{r['element']} {r['ion']}".strip()
            ep = (float(r["ep_eV"]) if "ep_eV" in d.columns and pd.notna(r.get("ep_eV"))
                  else None)
            ll = _match(linelist.assign(species=linelist["element"].astype(str) + " "
                                        + linelist["ion"].astype(str)), wl, ep, species)
            if ep is None and ll is not None:
                ep = float(ll["excitation_potential_eV"])
            gf = _match(canonical, wl, ep, species)

            in_agg = str(r["in_aggregate"]) == "True"
            status, code, note = _disposition(pc, species, wl, in_agg,
                                              str(r.get("excluded_reason", "")))
            engine = str(r["treatment"])
            rows.append({
                "element": r["element"], "ion": r["ion"],
                "wavelength_air_A": wl,
                "excitation_potential_eV": ep,
                # 🔴 canonical_gf OVERRIDES the linelist (RYA-834). Where canonical does
                # not carry the line the linelist value stands, and gf_source says so.
                "log_gf": (float(gf["log_gf"]) if gf is not None
                           else (float(ll["log_gf"]) if ll is not None
                                 and "log_gf" in linelist.columns else None)),
                "gf_source": (str(gf["loggf_reference"]) if gf is not None
                              else ("linelist (not in canonical_gf)" if ll is not None
                                    else "UNRESOLVED")),
                "gf_grade": _grade_label(
                    gf["nist_grade"] if gf is not None else None,
                    str(gf["loggf_reference"]) if gf is not None else ""),
                # RYA-850: the cited lab sigma where the line is in a primary-lab pool.
                # Blank is the honest value for a line nobody measured — never a default.
                "gf_sigma_dex": _lab_sigma(lab, wl, ep),
                "damping_rad": _damp(synth_ll, ll, wl, ep, "rad", "damping_rad"),
                "damping_stark": _damp(synth_ll, ll, wl, ep, "stark", "damping_stark"),
                "damping_vdW": _damp(synth_ll, ll, wl, ep, "waals", "damping_vdW"),
                # STATE THE FORM (RYA-489 §6.4) — per line, because the GES list mixes
                # them: some vdW entries are ABO packed sigma.alpha, others classical.
                "damping_form": _damping_form(
                    _damp(synth_ll, ll, wl, ep, "waals", "damping_vdW")),
                "hfs_isotope_note": ("HFS components in canonical_gf: "
                                     f"{gf.get('hfs_n_components')}" if gf is not None
                                     and "hfs_n_components" in canonical.columns else ""),
                "instrument": r["instrument"],
                "arm": _arm_of(wl),
                "method": ("ew_integration" if str(r.get("ew_inversion")) == "True"
                           else "synthesis_fit"),
                "ew_mA": r.get("ew_mA"),
                "reduced_ew": r.get("rew"),
                "red_chi2": r.get("red_chi2"),
                "sigma_A_dex": _sigma_A(r)[0],
                "sigma_A_basis": _sigma_A(r)[1],
                "engine": engine + (f" ({deck})" if deck else ""),
                "scale": ("1D-NLTE" if "NLTE" in engine.upper() else "1D-LTE"),
                "A_X_line": r.get("abundance"),
                "status": status, "reason_code": code, "reason_note": note,
            })

    out = pd.DataFrame(rows, columns=ROW_COLUMNS)
    # 🔴 ONE ROW PER (line x instrument x engine) — RYA-489 §6.4 and RYA-712. A repeated
    # key means the same physical measurement entered twice, which silently doubles its
    # weight for any consumer that groups by engine. Loud, with the offenders named.
    key = ["element", "ion", "wavelength_air_A", "instrument", "engine"]
    dup = out[out.duplicated(subset=key, keep=False)]
    if not dup.empty:
        show = (dup.groupby(key).size().sort_values(ascending=False).head(5).to_dict())
        raise PerLineProductError(
            f"{len(dup)} rows share a (line x instrument x engine) key — the same "
            f"measurement was read more than once, probably from a control or re-run "
            f"directory. Worst offenders: {show}")
    # 🔴 THE ACCOUNTING. Emitted must equal measured — a projection that filters is not a
    # projection (RYA-844).
    if len(out) != n_measured:
        raise PerLineProductError(
            f"row accounting failed: {n_measured} measured rows read, {len(out)} emitted. "
            f"A per-line product must account for every measured line, never filter one.")

    header = _build_header(star, element,
                           {"damping_source": damping_source,
                            "replication_grade": ("yes" if damping_source == "synthesis"
                                                  else "NO — damping is NOT the constant "
                                                       "the synthesis used"),
                            "band_products": ", ".join(str(p) for p in band_products),
                            "linelist": f"linelist_{star}.csv",
                            "canonical_gf": CANONICAL_GF.name,
                            "problem_children": REGISTRY_CSV.name},
                           input_commit)
    accounting = {
        "n_measured_rows_read": n_measured,
        "n_emitted_rows": len(out),
        "n_line_files": len(line_files),
        "by_status": out["status"].value_counts().to_dict(),
        "by_engine": out["engine"].value_counts().to_dict(),
        "n_gf_from_canonical": int((~out["gf_source"].isin(
            ["UNRESOLVED", "linelist (not in canonical_gf)"])).sum()),
        "by_gf_grade": out["gf_grade"].value_counts().to_dict(),
        "n_with_cited_lab_sigma": int((out["gf_sigma_dex"].astype(str) != "").sum()),
        "n_with_statistical_sigma": int((out["sigma_A_dex"].astype(str) != "").sum()),
    }
    return PerLineProduct(header=header, rows=out, accounting=accounting)
