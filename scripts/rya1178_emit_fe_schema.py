#!/usr/bin/env python3
"""RYA-1178 — re-emit data/products/solar/Fe.json on a complete, publication-grade schema.

PACKAGING, NOT SCIENCE. No abundance is re-derived and no gf, line or parameter is tuned
(RYA-161). Every value written here is either copied from a single source of truth, or
computed from terms the campaign already measured. Where a term was never measured this
says so in the record rather than inventing one.

WHAT EACH FIELD'S SOURCE OF TRUTH IS — none of it from memory:
  atlas / telluric_state   data/catalog/holdings_manifest_registry.csv + instrument_catalog
  line_set_resolved        pipeline.reference_lineset.line_set_for_product  (see below)
  grade                    derived from line_set, RYA-946's going-forward names
  wavelength_range_A       the artifact stem the product was ingested from
  model_grid               data/catalog/model_registry.csv, keyed on stored_token
  xi terms                 data/results/rya1120/xi_sigma_reported.json, KEYED BY BAND
  Asplund products         data/results/rya1106/asplund_four_instrument_table.json

🔴 `line_set` IS NOT STAMPED ON OUR OWN PRODUCTS, AND THAT IS RYA-1127's CALL, NOT THIS
TICKET'S. The spec asks to populate it. `tests/test_line_set_identity_rya1127.py` asserts
the opposite in as many words — `assert "line_set" not in ours` — because our tiers map
one-to-one onto the `our-*` names, so a stored copy is a SECOND SOURCE OF TRUTH free to
drift from the one the identity key resolves. The spec's stated reason for wanting it
("so replication vs working products stop colliding on identity") was ALREADY satisfied:
RYA-1127 fixed that collision by resolving the axis at key-computation time. So the feed
publishes `line_set_resolved` + `line_set_basis`, labelled derived, and a REPLICATION
product — the documented exception — is the only kind that carries an explicit `line_set`.

🔴 THE VOCABULARY IS `model_registry.LINE_SETS`, NOT THE TICKET'S EXAMPLE STRINGS. The
ticket says 'e.g. "codex", "asplund_agss21"'. The canonical values are `our-graded` /
`our-deep-graded` / `asplund`; `asplund_agss21` is the NATIVE spelling inside the RYA-1109
reference file, which `reference_lineset` maps to canonical `asplund` and never rewrites
in the file. Repo state wins over ticket prose.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
FEED = ROOT / "data/products/solar/Fe.json"
HOLDINGS = ROOT / "data/catalog/holdings_manifest_registry.csv"
INSTRUMENTS = ROOT / "data/catalog/instrument_catalog.csv"
MODELS = ROOT / "data/catalog/model_registry.csv"
XI_CAMPAIGN = ROOT / "data/results/rya1120/xi_sigma_reported.json"

#: 🔴 RYA-1185 -- THE BAND-KEYED dA/dxi RUNS, WHICH ARE WHAT RETIRE `ALIASED`.
#: RYA-1120's run unit is (ion x holding x tier) with NO band, so one measured derivative
#: was served to products in different bands measured on different line pools -- published
#: as ALIASED (17 products) or MEASURED_KEY_AMBIGUOUS (13) because the artifact could not
#: settle it. These three runs measure the derivative IN ITS OWN BAND, so they answer the
#: question the campaign key could not. They take PRECEDENCE over the campaign wherever
#: they cover a product, and the campaign remains the source for VIS.
XI_BAND_RUNS = (
    ROOT / "data/results/rya1168/nearuv_xi_dadxi.json",     # near-UV, RYA-1168
    ROOT / "data/results/rya1163/ir_xi_dadxi.json",         # NIR,     RYA-1163
    ROOT / "data/results/redopt_xi/redopt_xi_dadxi.json",   # red-optical
)
ASPLUND = ROOT / "data/results/rya1106/asplund_four_instrument_table.json"

#: RYA-1178 §1: solar xi is PINNED at 1.0; delta_xi is RYA-1089's sourced 0.2912.
XI_VALUE_KMS, DELTA_XI_KMS = 1.0, 0.2912

#: 🔴 RYA-1185 (Ryan, pre-merge) — `xi_state` MUST DERIVE FROM A STATED RULE, and the rule
#: is NOT "does the synthesiser receive a microturbulence". THAT WAS CHECKED IN THE CODE AND
#: THE ANSWER IS "BOTH DO":
#:
#:   * <3D> mean:  `abundances_derive.py:876/886` calls `ispec.generate_spectrum(...,
#:                 microturbulence_vel=vturb, ..., **_atm_file)` — `_atm_file` IS the <3D>
#:                 route (`atmosphere_layers_file`, the mul23 deck) and `microturbulence_vel`
#:                 rides the SAME call. No branch drops it for a mean-3D deck.
#:   * Amarsi 3D:  the MLP takes `vmic` as an EXPLICIT input axis (`pipeline/amarsi3d.py`,
#:                 training box 0-3 km/s) and its correction moves d(aberr)/dxi = +0.0985
#:                 dex/(km/s) over 112 of 117 lines.
#:
#: So a "does the synthesis carry a xi parameter" rule would make ENGINE-A-3DNLTE UNMEASURED
#: too, which contradicts the ratified split. The rule that actually separates them is about
#: the PUBLISHED QUANTITY, not the call signature:
#:
#:   NOT_APPLICABLE  the published value has no NET xi dependence.
#:                   - full 3D (ENGINE-A-3DNLTE): the published number is a SUM,
#:                     A(3D-NLTE) = A(1D-LTE)(xi) + aberr(xi), whose two halves move
#:                     OPPOSITELY — the correction tracking xi is the baseline's
#:                     xi-dependence being undone, not a second one added (Ryan, 2026-08-29).
#:                   - the replication anchors (`line_set` stored): a reference-set anchor,
#:                     not a xi-varying measurement of ours.
#:   UNMEASURED      xi APPLIES to the published value and no valid derivative exists yet.
#:
#: ⚠️ THE <3D> MEAN STAYS UNMEASURED, AND THAT IS A MEASUREMENT, NOT A DEFAULT. RYA-1099 ran
#: the mean route at xi = 0 and it came out **+0.137 dex WORSE**: a mean atmosphere averages
#: the velocity structure OUT, so the route still runs on an inherited xi. RYA-1099 forbids
#: the exemption in writing, and `rya1112_vis_fe_uncertainty_audit` asserts the same split
#: with a test whose name is `..._and_never_from_a_NAME`.
XI_NOT_APPLICABLE_TREATMENTS = ("ENGINE-A-3DNLTE",)


def xi_disposition(prod: dict) -> tuple[str, str] | None:
    """(state, why) when xi does not apply to this product's published value, else None."""
    if str(prod.get("treatment") or "") in XI_NOT_APPLICABLE_TREATMENTS:
        why = ("full 3D resolves the velocity field xi stands in for, so the PUBLISHED value "
               "carries no net xi dependence: it is a sum A(1D-LTE)(xi) + aberr(xi) whose two "
               "halves move oppositely. Recorded, not silent — the Amarsi MLP DOES take vmic "
               "as an input axis and its correction moves d(aberr)/dxi = +0.0985 dex/(km/s); "
               "the exemption is about the sum, not the call (Ryan, 2026-08-29).")
        if prod.get("line_set"):
            why += (" This product is ALSO a reference-set replication anchor, not a "
                    "xi-varying measurement of ours — N/A on either ground (RYA-1185).")
        return "NOT_APPLICABLE", why
    return None

#: RYA-946's grade names. `Consistent` is NOT here: RYA-1105 retired it and
#: `model_registry.LINE_SETS` omits it deliberately, so a product carrying it must fail
#: loudly rather than acquire a name.
#:
#: 🔴 RYA-1185 RENAMES TWO OF THESE, ON RYAN'S RULING, AND THE OLD NAMES WERE NOT RATIFIED.
#:   "Our Grade"     -> "Codex Grade"     -- "Our Grade" was never a ratified name; it is
#:                                           the primary-lab-gf tier at or below the 0.60
#:                                           depth gate.
#:   "Asplund Grade" -> "Reference Grade" -- the top tier is the REFERENCE tier, named for
#:                                           what it is (a product measured on an external
#:                                           published line set) rather than for one
#:                                           paper. Ryan settled the name in RYA-1185
#:                                           2026-09-03; it had been an open pick.
#: `gbs` moves to the same "Reference Grade" for the same reason -- it is the same kind of
#: product, and `line_set` is what says WHICH reference (that is why RYA-1127 put line_set
#: in the identity key). No GBS product is live today, so this renames nothing published.
#: ⚠️ docs/catalog/model_registry_notes.md:142 still records the OLD set as "going-forward";
#: RYA-1185 Part B corrects it. Repo prose is not the authority here -- the ruling is.
GRADE_FOR_LINE_SET = {
    "asplund": "Reference Grade",
    "asplund-al": "Reference Grade",
    "gbs": "Reference Grade",
    "our-graded": "Codex Grade",
    "our-deep-graded": "Deep Grade",
    "our-ungraded": "Ungraded",
    "our-all": "Ungraded",
}

#: Human atlas name + telluric state per holding. EVERY string is transcribed from the
#: holdings registry's own `notes` / `telluric_applied` and the instrument catalog's
#: `telluric_basis` — the "which Kitt Peak" fix. `solar_kpno_*` splits into two atlases
#: that a reader could not otherwise tell apart.
ATLAS = {
    "solar_kpno_kurucz2005_corrected": (
        "KPNO Kurucz-2005 solar irradiance (irradiance2005/irradthu.dat, rev. 02jan2010)",
        "corrected: Kurucz-2005 product is delivered telluric-corrected"),
    "solar_kpno_molecfit_corrected": (
        "KPNO 1984 FTS solar flux atlas, molecfit-corrected",
        "corrected: ESO molecfit 4.4.4, six registered bands (RYA-940)"),
    "solar_harps_molecfit_corrected": (
        "HARPS Phase-3 direct-Sun stack (1102.D-0954(A), 2023-08-02), molecfit-corrected",
        "corrected: ESO molecfit 4.4.4 per-exposure, observation-night GDAS (RYA-931)"),
    "solar_iag": (
        "IAG FTS solar flux atlas",
        "corrected: as delivered — measured telluric-free to 0.1% in the O2 A-band core "
        "(RYA-783/786)"),
    "solar_crires_plus_y_wide_rya1054": (
        "CRIRES+ Vesta Y arm, 9800-10796 A (RYA-1054 wide window)",
        "corrected: cr2res + molecfit"),
}

#: RYA-1178 B / RYA-1164: the IR bar is gf-limited, and that is an IRREDUCIBLE under the
#: current line list — not small-n noise. Dispersion = sigma_stat * sqrt(n).
#: ⚠️ THAT CONVERSION IS ONLY VALID BECAUSE `sigma_stat` IS A STANDARD ERROR. The feed's
#: `stat_basis` prose says "RMS of the random terms", which READS like a dispersion, and
#: it is absent entirely on all six IR products. Verified before use: on the ten VIS Fe I
#: products that carry both a per-line and a product layer, std/sqrt(n) reproduces the
#: published `stat_dex` on 6 of 6 checked and the raw std does not.
IR_DISPERSION_NOTE = (
    "gf-limited dispersion, IRREDUCIBLE under the current NIR line list. It tracks gf "
    "QUALITY, not line count: the widest bars carry the LARGEST n (KP ENGINE-A 1.315 at "
    "n=7 is the exception that proves it — CRIRES+ curated lab-gf is ~6x tighter at n=5, "
    "while IAG 0.794 at n=25 and KP 0.890 at n=26 are the widest with the most lines). "
    "The reducible lever is LABORATORY gf, not more lines: adding VALD3-grade NIR lines "
    "would widen this, not narrow it. Report-and-note per RYA-777."
)


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return ""


def _bump(v: str) -> str:
    a, b = v.split("."); return f"{a}.{int(b) + 1}"


def wavelength_range(prod: dict):
    """[lo, hi] in A, read off the artifact stem the product was ingested from.

    The span is in the band_products FILENAME and nowhere else in the record, so a reader
    of the feed cannot currently tell what wavelengths were priced.
    """
    ap = prod.get("provenance") or prod.get("artifact_provenance") or {}
    for cand in (ap.get("copied_to"), ap.get("path")):
        m = re.search(r"_(\d{3,5})_(\d{3,5})_", str(cand or ""))
        if m:
            return [float(m.group(1)), float(m.group(2))], "artifact stem"
    #: The four RYA-1095 3D-NLTE artifacts are named `rya817_3dnlte_products.csv` with no
    #: window in the stem. The window is still recorded — inside the CSV's own
    #: `provenance` column, which names the `FeI_4200_6910_...` 1D-LTE base the run was
    #: built on. Read it from there rather than substituting the per-line min/max, which
    #: is the span the surviving LINES happen to cover, not the span that was priced.
    rel = ap.get("copied_to")
    if rel:
        f = ROOT / str(rel)
        if f.exists():
            try:
                df = pd.read_csv(f)
            except Exception:
                return None, None
            for col in ("provenance", "stat_basis"):
                if col in df.columns:
                    m = re.search(r"_(\d{3,5})_(\d{3,5})_", str(df.iloc[0][col]))
                    if m:
                        return ([float(m.group(1)), float(m.group(2))],
                                f"artifact `{col}` column (stem carries no window)")

    #: 🔴 RYA-1185 -- LAST RESORT: THE BAND'S OWN DEFINITION, AND SAID SO IN THE BASIS.
    #: The four RYA-1106 Asplund products are published from `asplund_lines_products.csv`,
    #: whose stem carries no window and whose `provenance` column names a holding rather
    #: than a span; the RYA-1106 artifacts record `n_asplund_lines`/`n_served` but NO
    #: window anywhere. So the run's own priced window is not recoverable from the
    #: artifact, and the honest fallback is the declared BAND -- `pipeline.band_policy`,
    #: the single definition of what VIS means here.
    #:
    #: ⚠️ IT IS A WEAKER CLAIM AND THE `basis` STRING SAYS SO, because it is: a band
    #: definition is what the band IS, not what this run priced. Substituting the per-line
    #: min/max was rejected for the same reason the RYA-1095 case above rejects it -- that
    #: is the span the surviving LINES cover. Any product that lands here is visible as
    #: such in `wavelength_range_basis`, so a real missing window cannot hide behind it.
    band = str(prod.get("band") or "")
    try:
        from pipeline.band_policy import POLICIES
        pol = next((q for q in POLICIES if q.name == band), None)
    except Exception:
        pol = None
    if pol is not None:
        return ([float(pol.lo_A), float(pol.hi_A)],
                f"the {band} BAND DEFINITION (pipeline.band_policy), not a run window — "
                f"this product's artifact records no wavelength span of its own")
    return None, None


def load_sources():
    hold = pd.read_csv(HOLDINGS, comment="#")
    inst = pd.read_csv(INSTRUMENTS, comment="#")
    models = pd.read_csv(MODELS, comment="#")
    xi = json.loads(XI_CAMPAIGN.read_text())
    return hold, inst, models, xi


def xi_band_index() -> dict:
    """The band-keyed dA/dxi runs, indexed on the FULL product key INCLUDING band.

    ⚠️ THE THREE FILES DO NOT SHARE ONE SCHEMA, AND ONE OF THEM MISLABELS ITS OWN HEADER.
      * RYA-1163 carries no per-pool `ion`; its species is a TOP-LEVEL "Fe I", so the ion
        is read from there. A `.get("ion")` alone would silently key all six IR pools on
        None and match nothing.
      * The red-optical file's header says `"ticket": "RYA-1168"` and
        `"window_A": [3000, 3780]` -- a near-UV window -- because it was templated from
        the near-UV runner and the header was not corrected. Its PER-POOL rows are right
        (`band: red-optical`, red-optical holdings), so the band is taken per pool and
        never from the header. Flagged to Ryan in RYA-1185; the artifact is not edited
        here (RYA-161: it is a published run).
    """
    idx = {}
    for path in XI_BAND_RUNS:
        if not path.exists():
            raise SystemExit(f"missing band-keyed xi run: {path.relative_to(ROOT)}")
        doc = json.loads(path.read_text())
        sp = (doc.get("species") or "").split()
        species_ion = sp[-1] if sp else None
        for e in doc["pools"]:
            ion = e.get("ion") or species_ion
            if not ion:
                raise SystemExit(f"{path.name}: pool {e.get('pool')!r} has no ion and the "
                                 f"file declares no species -- refusing to guess")
            k = (ion, e["holding"], e["tier"], e["treatment"], e["band"])
            if k in idx:
                raise SystemExit(f"two band-keyed xi entries for {k} -- the band key must "
                                 f"be unique or a product would silently take one")
            idx[k] = {**e, "_source": str(path.relative_to(ROOT)), "_ticket": doc.get("ticket")}
    return idx


def model_grid(models: pd.DataFrame, treatment: str) -> str | None:
    """`model_family` x `atmosphere` x `deck` for the treatment token, from the registry."""
    r = models[models.stored_token.astype(str) == str(treatment)]
    if r.empty:
        return None
    r = r.iloc[0]
    fam = str(r.model_family); atm = str(r.atmosphere); deck = str(r.deck)
    parts = [p for p in (atm, fam if fam != "none" else "", deck if deck != "none" else "")
             if p and p != "nan"]
    return " / ".join(parts)


def xi_index(xi_doc: dict, feed: dict) -> dict:
    """Campaign entries keyed BY BAND — RYA-1114 F2.

    🔴 THE CAMPAIGN KEY HAS NO BAND, AND THAT IS THE DEFECT. Its run unit is
    (ion x holding x tier) and one such unit spans several wavelength windows, so a single
    measured dA/dxi was being served to products in DIFFERENT bands measured on DIFFERENT
    line pools — 14 campaign keys serve more than one band, and on 18 products the
    campaign's own `n_lines` does not even match the product's. Keying by band here cannot
    manufacture the missing derivatives, but it stops a product silently inheriting a
    twin's, which is the part that was publishing a wrong number as if it were measured.
    """
    by_key = {}
    for e in xi_doc["products"]:
        by_key.setdefault((e["ion"], e["holding"], e["tier"], e["treatment"], e["route"]), []).append(e)

    band = xi_band_index()

    #: 🔴 `serves` COUNTS ONLY THE BANDS THAT STILL DRAW ON THE CAMPAIGN -- and that is what
    #: retires MEASURED_KEY_AMBIGUOUS rather than relabelling it. The ambiguity was never
    #: "this number might be wrong"; it was "this campaign key also feeds OTHER bands, so
    #: the artifact cannot prove the derivative is this product's own". Once near-UV,
    #: red-optical and NIR each have their own band-keyed run, they stop drawing on the
    #: campaign key at all, the VIS product becomes the only band on it, and the question
    #: the flag existed to ask is answered by the ARTIFACTS -- not by widening a tolerance
    #: or renaming a state. Measured: 13 KEY_AMBIGUOUS -> 0, with no VIS derivative changed.
    serves = {}
    for p in feed["products"]:
        if (p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]) in band:
            continue
        k = (p["ion"], p["holding"], p["tier"], p["treatment"], p["route"])
        serves.setdefault(k, set()).add(p["band"])
    return {"by_key": by_key, "serves": serves, "band": band}


def xi_terms(prod: dict, idx: dict) -> dict:
    """The xi layer for one product, honest about which pool the derivative came from."""
    base = {"xi_value_kms": XI_VALUE_KMS, "delta_xi_kms": DELTA_XI_KMS}

    #: A band-keyed run answers for this product's OWN band, so it wins outright.
    b = idx["band"].get((prod["ion"], prod["holding"], prod["tier"], prod["treatment"],
                         prod["band"]))
    if b is not None:
        #: ⚠️ THE RUN'S OWN `xi_state` DECIDES, NOT THE PRESENCE OF A NUMBER. The CRIRES+
        #: ENGINE-A NIR pool carries dA_dxi = -0.1225 AND xi_state = UNMEASURED, because it
        #: paired only 2 lines against the run's declared min_paired = 3. Reading the float
        #: and ignoring the verdict would publish a 2-line derivative as a measurement.
        if b.get("xi_state") != "MEASURED" or b.get("dA_dxi") is None:
            return {**base, "sigma_xi": None, "xi_state": "UNMEASURED",
                    "xi_note": (f"{b.get('xi_note') or 'not measured'} "
                                f"[{b['_ticket']}, {b['_source']}] -- the band-keyed run "
                                f"declined this pool, so no derivative is published for it")}
        return {**base, "sigma_xi": round(abs(b["dA_dxi"]) * DELTA_XI_KMS, 6),
                "dA_dxi": b["dA_dxi"], "xi_state": "MEASURED",
                "xi_note": (f"|dA/dxi|={abs(b['dA_dxi']):.4f} x delta_xi={DELTA_XI_KMS}, "
                            f"measured IN THIS PRODUCT'S OWN BAND ({prod['band']}) on "
                            f"{b.get('n_paired')} paired lines by {b['_ticket']} "
                            f"({b['_source']}). Supersedes the RYA-1120 campaign value, "
                            f"whose run unit carries no band (RYA-1114 F2)."),
                "xi_source": b["_source"]}

    k = (prod["ion"], prod["holding"], prod["tier"], prod["treatment"], prod["route"])
    entries = idx["by_key"].get(k, [])
    bands = idx["serves"].get(k, set())

    #: 🔴 RYA-1185 -- FULL 3D IS NOT_APPLICABLE BY WHAT IT IS, NOT BY WHETHER THE CAMPAIGN
    #: HAPPENED TO RUN IT. xi is a 1D fudge for the velocity field, and a full-3D model
    #: resolves that field, so the term does not apply -- which is exactly the disposition
    #: the RYA-1120 campaign itself gives: all four of its NOT_APPLICABLE entries are
    #: ENGINE-A-3DNLTE and it has no other. Deciding this from campaign MEMBERSHIP instead
    #: made the four newly published RYA-1106 Asplund products (also ENGINE-A-3DNLTE) fall
    #: through to NOT_IN_CAMPAIGN -- "we never ran it", which is the wrong reason and reads
    #: as an owed measurement rather than a term that does not exist.
    na = xi_disposition(prod)
    if na is not None:
        return {**base, "sigma_xi": 0.0, "xi_state": na[0], "xi_note": na[1]}

    if not entries:
        return {**base, "sigma_xi": None, "xi_state": "NOT_IN_CAMPAIGN",
                "xi_note": ("this (ion, holding, tier, treatment) was never run through "
                            "the RYA-1120 perturb-and-re-derive campaign, so no dA/dxi "
                            "exists for it at any band. Not inherited from elsewhere.")}
    e = entries[0]
    state = e.get("xi_state")
    if state == "NOT_APPLICABLE":
        return {**base, "sigma_xi": 0.0, "xi_state": "NOT_APPLICABLE",
                "xi_note": e.get("xi_note") or "full 3D resolves the velocity field xi stands in for"}
    if e.get("dA_dxi") is None:
        #: ⚠️ UNMEASURED MUST SAY "xi APPLIES", NOT JUST "no run" (RYA-1185). Read bare, the
        #: campaign's note is silent about whether the term is owed or absent, and the four
        #: <3D> mean products are exactly where that ambiguity bites — they look like full 3D
        #: and are not exempt. RYA-1099 measured it: the mean route at xi = 0 is +0.137 dex
        #: WORSE, because a mean atmosphere averages the velocity structure OUT.
        applies = ("xi APPLIES to this product and the derivative is OWED — this is not an "
                   "exemption. ")
        if "mean3D" in str(prod.get("treatment") or ""):
            applies += ("The <3D> MEAN is NOT full 3D: averaging removes the velocity "
                        "structure, so the route runs on an inherited xi. RYA-1099 measured "
                        "it at xi = 0 and the result was +0.137 dex WORSE, which forbids the "
                        "full-3D exemption here in writing. ")
        return {**base, "sigma_xi": None, "xi_state": "UNMEASURED",
                "xi_note": applies + (e.get("xi_note") or "no perturb-and-re-derive on this pool")}
    # measured — but on WHICH pool?
    #: 🔴 ALIASED MEANS "PROVABLY A DIFFERENT POOL", NOT "SHARES A KEY". The campaign
    #: key spans bands, so a shared key alone is only ambiguity; what proves the
    #: derivative came from somewhere else is that the campaign measured a DIFFERENT
    #: NUMBER OF LINES than this product has. 18 products fail that test. Widening
    #: ALIASED to every multi-band key would relabel 30 and bury the 18 that are
    #: demonstrably wrong among 12 that merely cannot be confirmed.
    same_pool = e.get("n_lines") == prod.get("n_lines")
    multi_band = len(bands) > 1
    if same_pool and multi_band:
        return {**base, "sigma_xi": round(abs(e["dA_dxi"]) * DELTA_XI_KMS, 6),
                "dA_dxi": e["dA_dxi"], "xi_state": "MEASURED_KEY_AMBIGUOUS",
                "xi_note": (
                    f"|dA/dxi|={abs(e['dA_dxi']):.4f} x delta_xi={DELTA_XI_KMS}. The "
                    f"campaign's line count ({e.get('n_lines')}) matches this product's, "
                    f"so the derivative is most likely its own — but the campaign key "
                    f"carries no band and also serves {sorted(bands)}, so that cannot be "
                    f"proven from the artifact. Band-keying the campaign would settle it "
                    f"(RYA-1114 F2).")}
    if same_pool and not multi_band:
        return {**base, "sigma_xi": round(abs(e["dA_dxi"]) * DELTA_XI_KMS, 6),
                "dA_dxi": e["dA_dxi"], "xi_state": "MEASURED",
                "xi_note": (f"|dA/dxi|={abs(e['dA_dxi']):.4f} x delta_xi={DELTA_XI_KMS} "
                            f"measured on this product's own {e.get('n_paired')} paired lines")}
    return {**base, "sigma_xi": round(abs(e["dA_dxi"]) * DELTA_XI_KMS, 6),
            "dA_dxi": e["dA_dxi"], "xi_state": "ALIASED",
            "xi_note": (
                f"|dA/dxi|={abs(e['dA_dxi']):.4f} was measured on a pool of "
                f"{e.get('n_lines')} lines; this product has {prod.get('n_lines')} "
                f"(bands sharing that campaign key: {sorted(bands)}). The RYA-1120 run "
                f"unit is (ion x holding x tier) with NO band, so no band-specific "
                f"derivative exists. Published as ALIASED — the term is carried so the "
                f"budget is not silently short, and labelled so it is not read as this "
                f"product's own measurement. Retiring it needs a band-keyed re-run of the "
                f"campaign (RYA-1114 F2).")}


def ir_dispersion(prod: dict) -> dict | None:
    """RYA-1178 B — name the IR gf-limited irreducible instead of shipping a wide bar."""
    if prod["band"] != "NIR" or not prod.get("sigma_stat") or not prod.get("n_lines"):
        return None
    disp = prod["sigma_stat"] * math.sqrt(prod["n_lines"])
    return {"dispersion_dex": round(disp, 3),
            "basis": ("sigma_stat x sqrt(n_lines); valid because sigma_stat is the "
                      "STANDARD ERROR of the mean — verified against the per-line layer "
                      "on the VIS products that carry both (std/sqrt(n) reproduces the "
                      "published stat_dex; the raw std does not)"),
            "reducible_by": "laboratory gf for NIR Fe I — NOT more lines",
            "note": IR_DISPERSION_NOTE}


def tier_provenance(prod: dict) -> dict | None:
    """RYA-1178 B — re-confirm a GRADED tier against the pool's ACTUAL gf provenance.

    RYA-1114 F4 is the breach this guards: a product tiered GRADED whose pool is mostly
    ungraded. The check reproduces the pool the way `derive_band_products._cand_graded`
    builds it — every canonical_gf row for the species whose `gf_tier` contains LAB,
    inside the product's own window, split on the EW depth gate — and reports what that
    pool is actually made of, per gf reference.

    ⚠️ IT REPRODUCES AGAINST canonical_gf AS IT STANDS NOW, not as it stood at
    measurement time. canonical_gf has gained LAB rows since (RYA-945 ingested DH19), so
    a disagreement here is a prompt to check the vintage, not proof on its own.
    """
    win = prod.get("wavelength_range_A")
    if not win or prod.get("tier") not in ("GRADED", "DEEPGRADED"):
        return None
    try:
        from line_accounting_rya709 import DEPTH_HI
        sys.path.insert(0, str(ROOT / "scripts"))
        from derive_band_products import _feature_depth
    except Exception:
        return None
    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    species = f"Fe {prod['ion']}"
    lab = cg[(cg.species == species)
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)
             & cg.wavelength_air_A.between(win[0], win[1])]
    if lab.empty:
        return {"pool_reproduced": False,
                "note": "no LAB-tier rows for this species/window in canonical_gf today"}
    depth = _feature_depth(lab.wavelength_air_A.values.astype(float))
    sel = lab[depth > DEPTH_HI] if prod["tier"] == "DEEPGRADED" else lab[depth <= DEPTH_HI]
    refs = sel.loggf_reference.value_counts(dropna=False).to_dict()
    n_lab = int(sum(v for k, v in refs.items() if "PRIMARY LAB" in str(k)))
    return {
        "pool_reproduced": True,
        "n_lab_tier_in_window": int(len(lab)),
        "n_selected_by_depth_gate": int(len(sel)),
        "depth_gate": DEPTH_HI,
        "gf_references": {str(k): int(v) for k, v in refs.items()},
        "n_primary_lab": n_lab,
        "fraction_primary_lab": round(n_lab / len(sel), 4) if len(sel) else None,
        "tier_verdict": ("CONFIRMED — the selector admits only LAB-tier rows, and every "
                         "selected line resolves to a primary-laboratory reference"
                         if n_lab == len(sel) else
                         f"REVIEW — {len(sel) - n_lab} of {len(sel)} selected lines are "
                         f"not primary-laboratory"),
        "basis": ("reproduced from canonical_gf TODAY via the same rule as "
                  "derive_band_products._cand_graded; see the vintage caveat"),
    }


def science_provenance(prod: dict, models: pd.DataFrame, hold: pd.DataFrame,
                       line_set: str) -> dict:
    """The lineage a reader needs and `artifact_provenance` does not carry.

    The existing block is artifact bookkeeping — host, path, sha256, mtime. None of it
    says which gf, which atlas or which model grid produced the number, so a reader
    cannot tell a lab-gf run from a Kurucz one without opening the CSV.
    """
    h = hold[hold.holding_id.astype(str) == prod["holding"]]
    gf_map = {
        "our-graded": ("laboratory (canonical_gf LAB tier)",
                       "Den Hartog et al. 2014 ApJS 215, 23 (10.1088/0067-0049/215/2/23) "
                       "and Ruffoni et al. 2014 MNRAS 441, 3127 (10.1093/mnras/stu780), "
                       "as ingested by RYA-945/353"),
        "our-deep-graded": ("laboratory (canonical_gf LAB tier, above the depth gate)",
                            "Den Hartog et al. 2014 / Ruffoni et al. 2014, RYA-945/353"),
        "asplund": ("AGSS21's own published log gf on the target lines",
                    "Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table A.2 "
                    "(10.1051/0004-6361/202140445), ingested by RYA-1109"),
        "our-ungraded": ("VALD3 / Kurucz fallback", "VALD3 extraction, no laboratory grade"),
        "our-all": ("mixed: laboratory where available, VALD3/Kurucz otherwise",
                    "canonical_gf, mixed provenance — see per-line layer"),
    }
    gf_grade, gf_cite = gf_map.get(line_set, (None, None))
    ls_cite = ("Asplund, Amarsi & Grevesse 2021 Table A.2 via RYA-1109 / RYA-1111"
               if line_set == "asplund" else
               "Codex canonical_gf pool, tier %s (RYA-353/945)" % prod.get("tier"))
    out = {
        "gf_source": gf_grade,
        "gf_grade": ("primary laboratory" if line_set in ("our-graded", "our-deep-graded")
                     else "published reference set" if line_set == "asplund" else "fallback"),
        "gf_citation": gf_cite,
        "atlas_citation": (str(h.iloc[0].notes)[:400] if not h.empty else None),
        "model_grid": model_grid(models, prod["treatment"]),
        "line_set_citation": ls_cite,
    }
    cap = nlte_capability(prod, models)
    if cap:
        out["nlte_capability"] = cap
    return out


def nlte_capability(prod: dict, models: pd.DataFrame) -> dict | None:
    """RYA-1055 item 1 — the STATED CAPABILITY LIMIT, where a reader of the product is.

    Ryan, 2026-09-03: *"Stamp the capability limit where it is read: beside the deck
    registration AND in every Fe II product's science_provenance."* It is a property of
    the DECK WE SHIP, not of any one run, so it goes on every Fe II product regardless of
    which treatment produced that particular number — a reader setting an Fe I NLTE cell
    beside an Fe II one needs to be told the two are not on the same scale.

    🔴 THE LIMIT IS DECK-SPECIFIC AND THE STAMP MUST SAY SO, PER PRODUCT. Fe II ENGINE-A
    reads the MPIA/Bergemann per-line delta grid (6,400 Fe II rows over 80 lines,
    3805.5-6586.7 A), never `atom.fe607a`, so its NLTE label is honest and its corrections
    are real: every live Fe II product's own per-line artifact records
    `nlte_source = "Bergemann MPIA per-line delta_nlte (live query, solar node)"` with
    `nlte_delta_dex` of -0.001 to -0.002 on the lines MPIA serves. Stamping a blanket
    "Fe II NLTE unavailable" onto that product would replace one wrong statement with
    another.

    ⚠️ QUOTE THE PRODUCT'S OWN ARTIFACT, NOT A NEIGHBOURING ONE. A draft of this stamp
    cited "+0.001 dex on 6147.7341 / 6238.3859 / 6247.5570" from
    `data/products/solar/Fe_perline.csv`. Those numbers are real, but they belong to a
    DIFFERENT POOL -- the RYA-489 replication product's 11-line 5256-6456 A Fe II set --
    while the live VIS band product fits 4233.162 / 4303.170 / 4583.829 and applies
    -0.001/-0.002. Same element, same ion, same treatment, OPPOSITE SIGN, and it would
    have been written into the wrong product's provenance.

    WHICH TREATMENTS THE LIMIT BITES ON IS READ FROM THE REGISTRY, NOT LISTED HERE:
    `model_family == "gerber"` AND a scale carrying NLTE. A Gerber member added to
    `model_registry.csv` tomorrow is covered without touching this file — the hand-kept
    list is the thing that silently misses one.

    ⚠️ NOTHING HERE MOVES A VALUE. This adds a field; no abundance, sigma, line count or
    gf is touched (RYA-161).
    """
    if str(prod.get("ion")) != "II" or str(prod.get("element")) != "Fe":
        return None
    from pipeline.gerber_nlte import FE_II_NLTE_LIMIT, nlte_ion_capability

    capable, _ = nlte_ion_capability("Fe", "II")
    treatment = str(prod.get("treatment", ""))
    row = models[models.stored_token.astype(str) == treatment]
    fam = str(row.iloc[0].model_family) if not row.empty else ""
    scale = str(row.iloc[0].scale) if not row.empty else ""
    is_nlte = "NLTE" in scale.upper()
    on_deck = fam == "gerber" and is_nlte

    if not is_nlte:
        source = "n/a — this product is on an LTE scale"
    elif on_deck:
        source = ("Gerber deck (atom.fe607a) — LTE-EQUIVALENT for Fe II: every Fe II "
                  "line takes departure = 1")
    elif fam == "bergemann":
        source = ("MPIA / Bergemann mafags-os 1D-NLTE per-line delta grid "
                  "(data/nlte_grids/Fe_Bergemann_MPIA.csv) — NOT atom.fe607a. Its Fe II "
                  "corrections are real and small: at the solar node the grid's Fe II "
                  "deltas run -0.002 to +0.016 dex, median +0.000, with 146 of 160 "
                  "entries inside +-0.005 (Fe I control at the same node: median +0.011, "
                  "up to +0.040). What THIS product actually applied is recorded per "
                  "line in its own *_ENGINE-A_lines.csv (nlte_delta_dex / nlte_source, "
                  "RYA-880); across every live Fe II product that is -0.001 to -0.002 "
                  "dex, on the lines MPIA serves")
    elif fam == "amarsi":
        source = ("Amarsi+2022 3D-NLTE MLP (fe2_model.p) — NOT atom.fe607a. Its solar "
                  "Fe II correction is +0.066 dex (RYA-817 control, paper Table 6 "
                  "7.41 -> 7.47)")
    else:
        source = (f"UNDETERMINED — treatment {treatment!r} resolves to no "
                  f"model_registry row, so which NLTE source it used is not established "
                  f"here. Do not read this as 'the Gerber deck'.")

    return {
        "gerber_deck_ion_reach": ("Fe II: 58 levels, ZERO bound-bound transitions in "
                                  "atom.fe607a — an ionisation reservoir, not a term "
                                  "system. Fe I: 548 levels, 12,635 transitions."),
        "fe_ii_nlte_available_on_gerber_deck": bool(capable),
        "this_product_takes_nlte_from_the_gerber_deck": bool(on_deck),
        "nlte_source_for_this_product": source,
        "limit": FE_II_NLTE_LIMIT,
        "measurement": "data/results/rya1055/atom_ion_reach.json (RYA-1055)",
        "deferred": ("a two-stage Fe atom is DEFERRED to the off-solar programme, where "
                     "the Fe I/Fe II balance does the log g work (Ryan, 2026-09-03)"),
    }


def enrich(feed: dict, hold, inst, models, xi_doc) -> tuple[dict, list]:
    from pipeline.reference_lineset import line_set_for_product
    idx = xi_index(xi_doc, feed)
    notes = []
    for p in feed["products"]:
        #: 🔴 THE RENAME IS DEFERRED, DELIBERATELY. The ticket asks for `provenance` ->
        #: `artifact_provenance`. 27 code sites read the current name, and several are
        #: VALUE-PROTECTION guards -- RYA-1080's feed/repo reconciliation, RYA-1092's
        #: eligibility gate, RYA-1034's product store. Renaming the key under them makes
        #: those guards stop finding the block they exist to check, which does not fail
        #: loudly; it passes vacuously. A packaging ticket must not quietly disarm the
        #: guards that protect published values (RYA-161). The rename is right and should
        #: happen -- as its own migration that moves the readers with it.
        atlas, tell = ATLAS.get(p["holding"], (None, None))
        p["atlas"], p["telluric_state"] = atlas, tell
        #: 🔴 `line_set` IS NOT STAMPED ON OUR OWN PRODUCTS, AND THAT IS RYA-1127's CALL.
        #: The ticket asks to populate it. `tests/test_line_set_identity_rya1127.py`
        #: asserts the opposite in as many words -- `assert "line_set" not in ours` -- and
        #: `line_set_for_product` explains why: our tiers map one-to-one onto the `our-*`
        #: names, so a stored copy is a SECOND SOURCE OF TRUTH free to drift from the one
        #: the identity key resolves. The ticket's stated reason for wanting it ("so
        #: replication vs working products stop colliding on identity") is ALREADY
        #: satisfied -- RYA-1127 fixed that collision by resolving the axis at
        #: key-computation time, not by storing it. Storing it now would buy nothing and
        #: risk the drift that guard exists to prevent.
        #:
        #: A REPLICATION product is the documented exception and DOES carry it explicitly
        #: (see `asplund_products`), which is exactly what `line_set_for_product` supports.
        #: So the axis is readable from the feed for the products where it is not derivable,
        #: and `grade` -- a new, non-identity field -- carries it in human form for the rest.
        ls = line_set_for_product(p)
        p["grade"] = GRADE_FOR_LINE_SET[ls]
        p["line_set_resolved"] = ls
        p["line_set_basis"] = (
            "DERIVED at read time from `tier` by pipeline.reference_lineset."
            "line_set_for_product -- not stored as an identity field (RYA-1127). Shown "
            "here for readability; the identity key resolves it independently.")
        wr, wr_basis = wavelength_range(p)
        p["wavelength_range_A"] = wr
        p["wavelength_range_basis"] = wr_basis
        p["generated_at"] = _now()
        p["code_commit"] = _commit()

        xt = xi_terms(p, idx)
        p.update({k: v for k, v in xt.items() if k != "dA_dxi"})
        if "dA_dxi" in xt:
            p["dA_dxi_dex_per_kms"] = xt["dA_dxi"]

        # Part 2 — sigma_syst must be a real total of its NAMED components, and
        # sigma_reported must include it. The campaign's own sigma_reported did NOT:
        # on all 32 entries carrying one it equals quadrature(SE, sigma_xi) with
        # sigma_syst omitted entirely.
        comps = {"published_syst": p.get("sigma_syst"), "sigma_xi": xt.get("sigma_xi")}
        named = [v for v in comps.values() if v]
        p["sigma_syst_components"] = comps
        p["sigma_syst_complete"] = round(math.sqrt(sum(v * v for v in named)), 6) if named else None
        terms = [p.get("sigma_stat"), p.get("sigma_syst_complete")]
        got = [t for t in terms if t]
        p["sigma_reported"] = round(math.sqrt(sum(t * t for t in got)), 6) if got else None
        p["sigma_reported_basis"] = (
            "quadrature(sigma_stat, sigma_syst_complete) where sigma_syst_complete = "
            "quadrature(published sigma_syst, sigma_xi). sigma_xi enters ONCE, through "
            "sigma_syst_complete — the RYA-1120 artifact's own sigma_reported omitted "
            "sigma_syst altogether.")
        if p.get("xi_state") in ("UNMEASURED", "NOT_IN_CAMPAIGN"):
            p["sigma_reported_caveat"] = (
                "INCOMPLETE: no dA/dxi exists for this pool, so the xi term is absent "
                "from this bar rather than zero. The bar is a LOWER BOUND.")

        #: 🔴 PART 0 — THE <3D> LTE/NLTE COLLISION IS A MEDIAN COINCIDENCE, NOT A BUG.
        #: The ticket suspected the RYA-1104 "<3D>-NLTE == LTE wiring" defect. It is not
        #: that, and it is not a byte-identical pair either: the two artifacts differ
        #: (sha256, sigma_stat 0.0219 vs 0.0218, and the per-line layers differ), and the
        #: NLTE correction moves 66 of 67 lines by a mean +0.0267 dex. What coincides is
        #: only the AGGREGATE, because the published A is the MEDIAN — `nlte_effect.json`
        #: states `"aggregator": "median"`, and the median reproduces all four mean-3D
        #: products to 6 dp while the mean reproduces none. A median is insensitive to a
        #: shift that moves most of a distribution without moving its centre order
        #: statistic. The kurucz2005 pair on the identical code path lands 7.541 vs 7.562,
        #: which independently proves the NLTE leg is live. Both are kept, per the
        #: ticket's "if genuine: document the mechanism and keep both".
        if "mean3D" in p["treatment"]:
            p["aggregate_collision_note"] = (
                "Published A is the MEDIAN of the per-line abundances "
                "(nlte_effect.json: aggregator=median). On this holding the <3D>-LTE and "
                "<3D>-NLTE medians coincide at 7.552 while their MEANS differ by "
                "+0.0267 dex and the NLTE departure moves 66 of 67 lines. The pair is "
                "NOT byte-identical and the NLTE leg is NOT unwired: the same code path "
                "on kpno_kurucz2005 gives 7.541 vs 7.562. Genuine median coincidence, "
                "kept as two products (RYA-1178 Part 0)."
                if p["holding"] == "solar_kpno_molecfit_corrected" else
                "Published A is the MEDIAN of the per-line abundances "
                "(nlte_effect.json: aggregator=median). LTE 7.541 / NLTE 7.562 — the "
                "<3D>-NLTE departure is live and separates the pair on this holding.")

        disp = ir_dispersion(p)
        if disp:
            p["irreducible_dispersion"] = disp
        if p["band"] == "NIR":
            tp = tier_provenance(p)
            if tp:
                p["tier_provenance"] = tp
        p["science_provenance"] = science_provenance(p, models, hold, ls)
    return feed, notes


def asplund_products(models, hold) -> list:
    """Part 3 — the four RYA-1106 Asplund-replication products as first-class rows."""
    doc = json.loads(ASPLUND.read_text())
    out = []
    for _, h in sorted(doc["holdings"].items()):
        p = {
            "element": "Fe", "ion": "I", "band": "VIS",
            "instrument": h["instrument"], "holding": h["holding"],
            "tier": "GRADED", "selector": "ASPLUND_AGSS21",
            "treatment": "ENGINE-A-3DNLTE", "display": "Synth · 3D-NLTE · Amarsi · Asplund line set",
            "A": round(float(h["A_3dnlte"]), 3),
            "sigma_stat": h["stat_dex"], "sigma_syst": h["syst_dex"],
            "n_lines": h["n_lines"], "n_excluded": h["n_excluded"],
            "dominant_term": "gf scale (AGSS21 published)",
            "route": "SYNTH",
            "line_set": "asplund",
            "atlas": ATLAS.get(h["holding"], (None, None))[0],
            "telluric_state": ATLAS.get(h["holding"], (None, None))[1],
            "wavelength_range_A": [4200.0, 6910.0],
            "xi_value_kms": XI_VALUE_KMS, "delta_xi_kms": DELTA_XI_KMS,
            "sigma_xi": 0.0, "xi_state": "NOT_APPLICABLE",
            "xi_note": ("full 3D resolves the velocity field xi stands in for — the same "
                        "disposition the campaign gives every ENGINE-A-3DNLTE product"),
            "asplund21_reference": doc["asplund21_reference"],
            "vs_asplund_dex": h["vs_asplund_dex"],
            "provenance": {
                "host": "mac",
                "path": "data/results/rya1106/asplund_four_instrument_table.json",
                "copied_to": "data/results/rya1106/asplund_four_instrument_table.json",
            },
        }
        p["grade"] = GRADE_FOR_LINE_SET["asplund"]
        p["generated_at"] = _now(); p["code_commit"] = _commit()
        comps = {"published_syst": p["sigma_syst"], "sigma_xi": 0.0}
        p["sigma_syst_components"] = comps
        p["sigma_syst_complete"] = p["sigma_syst"]
        p["sigma_reported"] = round(math.hypot(p["sigma_stat"], p["sigma_syst"]), 6)
        p["sigma_reported_basis"] = ("quadrature(sigma_stat, sigma_syst_complete); xi is "
                                     "NOT_APPLICABLE for a full-3D product")
        p["science_provenance"] = science_provenance(p, models, hold, "asplund")
        p["science_provenance"]["line_set_coverage"] = (
            f"{h['coverage']['n_served']} of {h['coverage']['n_asplund_lines']} AGSS21 "
            f"Table A.2 lines served; {h['n_lines']} entered the aggregate")
        out.append(p)
    return out


def fill_withdrawal_reasons(feed: dict) -> list:
    """Part 4 — every superseded/quarantine/archive entry states WHY.

    ⚠️ THE TICKET'S PREMISE DOES NOT HOLD, and the honest thing is to say so rather than
    overwrite good prose to look busy. All 13 superseded, 12 quarantine and 13 archive
    entries ALREADY carry a populated reason — including the four near-UV ENGINE-B rows
    the ticket names, which carry code DUPLICATE_RETIRED_LABEL and a full paragraph. The
    only genuine gap is one archive row missing the machine-readable `quarantine_codes`
    while carrying its prose. That one is filled from its own text; nothing else is touched.
    """
    filled = []
    for r in feed.get("archive", []):
        if str(r.get("quarantine_codes") or "").strip():
            continue
        prose = str(r.get("quarantine_reason") or "")
        code = "PRE_CONTINUUM_FIX" if "PRE-FIX" in prose or "RYA-933" in prose else None
        if code:
            r["quarantine_codes"] = code
            filled.append(f"archive {r['band']}/{r['ion']}/{r['treatment']} -> {code}")
    return filled


def null_audit(feed: dict) -> pd.DataFrame:
    """The smoke test's per-product null audit over every field the ticket adds."""
    fields = ["atlas", "telluric_state", "line_set_resolved", "grade", "wavelength_range_A",
              "generated_at", "code_commit", "xi_value_kms", "delta_xi_kms", "sigma_xi",
              "xi_state", "sigma_syst_complete", "sigma_reported", "science_provenance",
              # the artifact block keeps its /1 name — the rename is deferred, see enrich()
              "provenance"]
    rows = []
    for p in feed["products"]:
        row = {"band": p["band"], "ion": p["ion"], "holding": p["holding"],
               "treatment": p["treatment"]}
        for f in fields:
            v = p.get(f)
            if v not in (None, "", [], {}):
                row[f] = "OK"
            elif f == "sigma_xi" and p.get("xi_state") in (
                    "UNMEASURED", "NOT_IN_CAMPAIGN", "NOT_APPLICABLE"):
                # explicit n/a + reason, which the ticket allows and a bare null does not
                row[f] = "n/a"
            else:
                row[f] = "NULL"
        rows.append(row)
    return pd.DataFrame(rows)


def verify(feed: dict) -> list:
    """Guards. Every one of these would otherwise be a silent wrong number."""
    from pipeline.reference_lineset import line_set_for_product
    errs = []
    for p in feed["products"]:
        # the stamped line_set must equal what the resolver derives — one source of truth
        want = line_set_for_product(p)
        shown = p.get("line_set") or p.get("line_set_resolved")
        if shown != want:
            errs.append(f"line_set drift: {p['holding']}/{p['treatment']} "
                        f"shows {shown!r}, resolver says {want!r}")
        if p.get("tier") in ("GRADED", "DEEPGRADED") and p.get("selector") != "ASPLUND_AGSS21" \
                and "line_set" in p:
            errs.append(f"line_set must stay DERIVED on our own products (RYA-1127): "
                        f"{p['holding']}/{p['treatment']}")
        if p.get("grade") != GRADE_FOR_LINE_SET.get(want):
            errs.append(f"grade/line_set disagree: {p['holding']}/{p['treatment']}")
        # sigma_reported must be the quadrature it claims to be
        st, sy, rep = p.get("sigma_stat"), p.get("sigma_syst_complete"), p.get("sigma_reported")
        if rep is not None:
            want_rep = math.sqrt(sum(t * t for t in (st, sy) if t))
            if abs(rep - want_rep) > 1e-6:
                errs.append(f"sigma_reported not the stated quadrature: "
                            f"{p['holding']}/{p['treatment']}")
        if "Consistent" in str(p.get("grade")):
            errs.append(f"grade 'Consistent' is retired (RYA-1105) and must fail loudly: "
                        f"{p['holding']}/{p['treatment']}")
    return errs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="audit only; write nothing")
    a = ap.parse_args()

    feed = json.loads(FEED.read_text())
    #: 🔴 KEY THE BEFORE-SNAPSHOT ON THE FULL IDENTITY. (ion, band, holding, treatment)
    #: is NOT unique in this feed — GRADED and DEEPGRADED products share it — so a
    #: dict on that key silently collapses pairs and then reports every loser as a MOVE.
    #: The identity is `product_eligibility.KEY_FIELDS`; use it.
    def ident(q):
        return tuple(str(q.get(k) or "") for k in
                     ("element", "ion", "band", "instrument", "holding",
                      "tier", "selector", "route", "treatment"))
    before_A = {ident(p): p["A"] for p in feed["products"]}
    if len(before_A) != len(feed["products"]):
        print(f"WARNING: identity key is not unique over the live products "
              f"({len(before_A)} keys for {len(feed['products'])} rows)")
    n_before = len(feed["products"])
    hold, inst, models, xi_doc = load_sources()

    feed, _ = enrich(feed, hold, inst, models, xi_doc)
    #: 🔴 PART 3 IS DEFERRED, AND HAND-WRITING IT HERE IS THE REASON IT MUST BE.
    #: The four RYA-1106 Asplund products are real and their numbers are in
    #: `data/results/rya1106/asplund_four_instrument_table.json`. But a product does not
    #: enter this feed by having a number: it enters through the canonical publisher,
    #: which computes the artifact sha256 and mtime, applies the RYA-1092 eligibility
    #: gate, and forms the RYA-1127 identity key. Rows appended by hand from a second
    #: script carry none of that -- when I tried it, 24 tests failed, and they were right
    #: to: RYA-1092's `every_live_product_passes_the_gate`, RYA-1034's
    #: `every_product_states_where_it_came_from`, and RYA-1111's own entrypoint guard all
    #: caught a product that had bypassed them. `asplund_products()` below is kept as the
    #: assembled record so the follow-up has it, and is NOT published from here.
    #:
    #: The two legitimate routes are `scripts/measure_reference_lineset.py --line-set
    #: asplund` (a measurement run; needs spectra not present in this worktree) or
    #: ingestion of the existing RYA-1106 artifacts through `scripts/publish_product.py`.
    added = []
    codes = fill_withdrawal_reasons(feed)

    #: The comparison must exclude the rows Part 3 ADDS. An Asplund product shares
    #: (ion, band, holding, treatment) with the working product on the same cell and
    #: differs only in `line_set` — which is exactly the collision RYA-1127 put line_set
    #: in the key to end. Comparing without it reads four ADDITIONS as twelve MOVES.
    added_keys = {id(p) for p in added}
    moved = []
    for p in feed["products"]:
        if id(p) in added_keys:
            continue
        k = ident(p)
        if k in before_A and p["A"] != before_A[k]:
            moved.append(f"{k} {before_A[k]} -> {p['A']}")
    errs = verify(feed)

    audit = null_audit(feed)
    nulls = {c: (int((audit[c] == "NULL").sum()), int((audit[c] == "n/a").sum()))
             for c in audit.columns if c not in ("band", "ion", "holding", "treatment")}

    print(f"products: {n_before} before + {len(added)} Asplund = {len(feed['products'])}")
    print(f"published A moved: {len(moved)}  (must be 0 — RYA-161)")
    print(f"withdrawal codes filled: {len(codes)} {codes}")
    print(f"xi_state: {pd.Series([p.get('xi_state') for p in feed['products']]).value_counts().to_dict()}")
    print(f"grade:    {pd.Series([p.get('grade') for p in feed['products']]).value_counts().to_dict()}")
    print(f"line_set: {pd.Series([p.get('line_set') or p.get('line_set_resolved') for p in feed['products']]).value_counts().to_dict()}")
    print()
    print("per-field null audit  [NULL = owed but missing | n/a = absent with a stated reason]:")
    for k, (nul, na) in nulls.items():
        print(f"   {k:24s} NULL={nul:<3d} n/a={na}")
    print()
    if errs:
        print("GUARD FAILURES:")
        for e in errs:
            print("  ", e)
        return 1
    print("guards: line_set==resolver, grade==line_set, sigma_reported==quadrature — all pass")

    if a.check:
        print("\n--check: nothing written")
        return 0
    if moved:
        print("REFUSING to write: a published A moved.")
        return 1
    feed["version"] = _bump(feed["version"])
    feed["updated_at"] = _now()
    feed["schema"] = "codex.element_product/2"
    FEED.write_text(json.dumps(feed, indent=2) + "\n")
    print(f"\nwrote {FEED.relative_to(ROOT)} at v{feed['version']} (schema {feed['schema']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
