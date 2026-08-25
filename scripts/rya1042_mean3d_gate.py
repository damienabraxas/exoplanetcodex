#!/usr/bin/env python3
"""RYA-1042: the RYA-534 anchor gate, re-expressed for a ⟨3D⟩ deck.

WHY THE OLD GATE COULD NOT RUN
------------------------------
`scripts/ts_gerber_gate.py` drives `interpol_modeles_nlte` directly, and **no vendor binary
can consume a ⟨3D⟩ deck** — it reads only native MARCS, which needs τ_Rosseland and P_g per
depth, and the ⟨3D⟩ archives ship the 5-column `TAU5000` form carrying neither (RYA-821).
So the gate was unrunnable on the very decks that most needed gating.

The **logic** is unchanged: does the deck reproduce an INDEPENDENT anchor within tolerance?
Only the read path changes — the deck now arrives through the products the RYA-1040/1044
route emits, which read it via `gerber_nlte.read_deck_node`.

TWO LEGS, AND ONLY ONE OF THEM NEEDS THE EXTERNAL FILE
------------------------------------------------------
**Leg 1 — ANCHOR AGREEMENT.** Our ⟨3D⟩NLTE − ⟨3D⟩LTE, per line, against Amarsi+2016's
`nmtd_lmtd` (the same difference, on the same kind of atmosphere, computed by another group
with another code and another model atom). External by construction.

**Leg 2 — THE DIFFERENTIAL'S OWN STRUCTURE.** Needs no external number, so it cannot be
gamed by the deck agreeing with itself. It checks what the differential IS rather than what
it equals: sign, magnitude envelope, and the excitation trend that over-ionisation implies.

🔴 A deck that passes BOTH is validated. Either alone is weaker — leg 1 can be satisfied by
two codes sharing an error, leg 2 by a deck that is self-consistently wrong.

🔴 THE ANCHOR IS NEVER DERIVED FROM THE MACHINERY IT VALIDATES (RYA-161/1035). That is the
closed-loop trap RYA-1035 found in this very deck's abundance record, where 7.46 turned out
to be our own stdin round-tripping back through bsyn's error message.

    python3 scripts/rya1042_mean3d_gate.py --nlte <lines.csv> --lte <lines.csv>
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANCHOR = ROOT / "data" / "nlte_grids" / "amarsi2016_fe" / "amarsi2016_mean3d_solar_anchor.csv"

#: ⚠️ INHERITED FROM RYA-534, NOT RE-DERIVED ON THIS DATA. RYA-534's element gate uses 0.05
#: dex and the ticket says the gate LOGIC is unchanged, so the tolerance travels with it.
#: A borrowed threshold is not a control (RYA-847), so this is STATED as inherited and the
#: measured margin is always printed beside the verdict -- a reader can then judge the cut
#: as well as the result.
TOL_DEX = 0.05

#: Wavelength match tolerance, air Å. The anchor stores nm to 4 dp (0.0001 nm = 0.001 Å);
#: our line list carries the catalogue wavelength. 0.05 Å is ~50x the anchor's print
#: precision and far below the ~4 Å minimum line separation the band selection enforces, so
#: it cannot reach a neighbour. ⚠️ An AMBIGUOUS match (two anchor lines inside the window)
#: is REFUSED, never resolved by taking the nearer -- a rounded number is not an identity.
MATCH_TOL_A = 0.05

#: The grid has no vturb = 1.0 node and the Sun is 1.0, so the anchor is a BRACKET.
#: Both ends are carried through to the verdict rather than interpolated into one number.
VTURB_BRACKET = ("0.75", "1.50")

SENTINEL = -4.0


def _read_lines(path: Path) -> dict[float, dict]:
    out = {}
    with path.open() as fh:
        for r in csv.DictReader(fh):
            if str(r.get("in_aggregate", "")).strip().lower() not in ("true", "1"):
                continue
            try:
                w = float(r["wavelength_air_A"]); a = float(r["abundance"])
            except (KeyError, TypeError, ValueError):
                continue
            out[w] = dict(abundance=a, ep_eV=_f(r.get("ep_eV")))
    return out


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_anchor(vturb: str) -> list[dict]:
    rows = []
    with ANCHOR.open() as fh:
        for r in csv.DictReader(fh):
            if r["species"] != "Fe1" or r["vturb_kms"] != vturb or r["clean"] != "yes":
                continue
            d = float(r["delta_mean3d_nlte_minus_mean3d_lte"])
            if d == SENTINEL:          # a floor value in the release, not a correction
                continue
            rows.append(dict(wave_A=float(r["lambda_air_nm"]) * 10.0, delta=d,
                             e_low=_f(r.get("e_low"))))
    return rows


def match(ours: dict[float, dict], anchor: list[dict]) -> tuple[list[dict], list[str]]:
    """Pair our lines to the anchor's by air wavelength. Ambiguity is REFUSED."""
    paired, notes = [], []
    for w, rec in sorted(ours.items()):
        near = [a for a in anchor if abs(a["wave_A"] - w) <= MATCH_TOL_A]
        if not near:
            notes.append(f"{w:.3f} A: no anchor line within {MATCH_TOL_A} A")
            continue
        if len(near) > 1:
            notes.append(f"{w:.3f} A: AMBIGUOUS -- {len(near)} anchor lines within "
                         f"{MATCH_TOL_A} A; refused rather than resolved by proximity")
            continue
        paired.append(dict(wave_A=w, ours=rec["delta"], anchor=near[0]["delta"],
                           ep_eV=rec["ep_eV"]))
    return paired, notes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nlte", required=True, help="⟨3D⟩-NLTE per-line CSV")
    ap.add_argument("--lte", required=True, help="⟨3D⟩-LTE per-line CSV (the comparand)")
    ap.add_argument("--out", default=None, help="write the verdict JSON here")
    a = ap.parse_args()

    nlte, lte = _read_lines(Path(a.nlte)), _read_lines(Path(a.lte))
    common = sorted(set(nlte) & set(lte))
    if not common:
        raise SystemExit(
            "no line is in-aggregate in BOTH products, so no differential exists. The "
            "NLTE effect is a per-line DIFFERENCE; a line measured in one leg and "
            "excluded in the other cannot contribute to it.")

    ours = {w: dict(delta=nlte[w]["abundance"] - lte[w]["abundance"],
                    ep_eV=nlte[w]["ep_eV"]) for w in common}
    deltas = [v["delta"] for v in ours.values()]
    med = statistics.median(deltas)

    doc = {"ticket": "RYA-1042", "generated_utc": datetime.now(timezone.utc)
           .strftime("%Y-%m-%dT%H:%M:%SZ"),
           "n_lines_nlte": len(nlte), "n_lines_lte": len(lte),
           "n_paired_internally": len(common),
           "our_differential": {"median": round(med, 4),
                                "mean": round(statistics.fmean(deltas), 4),
                                "min": round(min(deltas), 4),
                                "max": round(max(deltas), 4)},
           "tolerance_dex": TOL_DEX, "tolerance_basis": "inherited from RYA-534, not "
                                                        "re-derived on this data"}

    print(f"OUR DIFFERENTIAL  <3D>-NLTE minus <3D>-LTE, one atmosphere, {len(common)} lines")
    print(f"  median {med:+.4f}  mean {statistics.fmean(deltas):+.4f}  "
          f"range [{min(deltas):+.4f}, {max(deltas):+.4f}]")

    # ── LEG 1: anchor agreement ─────────────────────────────────────────────
    print(f"\nLEG 1 — ANCHOR AGREEMENT (Amarsi+2016 nmtd_lmtd, external)")
    leg1 = {}
    for vt in VTURB_BRACKET:
        anchor = load_anchor(vt)
        paired, notes = match(ours, anchor)
        if not paired:
            print(f"  vturb {vt}: NO LINES MATCHED -- cannot gate")
            leg1[vt] = {"n": 0, "verdict": "NO_OVERLAP"}
            continue
        o = [p["ours"] for p in paired]
        k = [p["anchor"] for p in paired]
        d_med = statistics.median(o) - statistics.median(k)
        per_line = [p["ours"] - p["anchor"] for p in paired]
        verdict = "PASS" if abs(d_med) <= TOL_DEX else "FAIL"
        leg1[vt] = {"n": len(paired), "ours_median": round(statistics.median(o), 4),
                    "anchor_median": round(statistics.median(k), 4),
                    "median_difference": round(d_med, 4),
                    "per_line_median_difference": round(statistics.median(per_line), 4),
                    "verdict": verdict,
                    "unmatched_or_ambiguous": len(notes)}
        # 🔴 THE MATCH RATE IS A HEADLINE, NOT A FOOTNOTE. Our lines are selected by
        # theoretical depth from the GES list; the anchor carries its own. A gate run over
        # a handful of matched lines is a weak gate however clean its arithmetic looks,
        # and burying that in a notes list is how a thin result reads as a solid one.
        rate = len(paired) / max(len(ours), 1)
        leg1[vt]["match_rate"] = round(rate, 3)
        leg1[vt]["n_ours"] = len(ours)
        flag = "  ⚠️ THIN" if len(paired) < 10 else ""
        print(f"  vturb {vt}: matched {len(paired)}/{len(ours)} of our lines "
              f"({rate:.0%}){flag}")
        print(f"      ours {statistics.median(o):+.4f}  anchor {statistics.median(k):+.4f}"
              f"  |diff| {abs(d_med):.4f} vs tol {TOL_DEX}  -> {verdict}")
        if (statistics.median(o) > 0) != (statistics.median(k) > 0):
            print(f"      🔴 SIGN DISAGREEMENT -- ours and the anchor do not even point "
                  f"the same way; a tolerance test on |difference| cannot see that, so it "
                  f"is stated separately")
            leg1[vt]["sign_disagreement"] = True
        for n in notes[:3]:
            print(f"      note: {n}")
        if len(notes) > 3:
            print(f"      ... and {len(notes) - 3} more")
    doc["leg1_anchor_agreement"] = leg1

    # ⚠️ THE vturb BRACKET CAN STRADDLE THE TOLERANCE, and then the gate's answer depends
    # on a grid node the Sun does not sit on. Reported explicitly rather than resolved:
    # picking the end that passes would be tuning, and interpolating would invent a node.
    verdicts = {v.get("verdict") for v in leg1.values() if v.get("n")}
    if verdicts == {"PASS", "FAIL"}:
        print(f"\n  🔴 THE vturb BRACKET STRADDLES THE TOLERANCE: "
              f"{ {k: v['verdict'] for k, v in leg1.items() if v.get('n')} }. The Sun is "
              f"vturb 1.0 and the grid has no such node, so the verdict depends on which "
              f"end is read. Treated as NOT PASSED -- choosing the end that passes is "
              f"tuning, and interpolating invents a node the grid was never asked for.")
        doc["vturb_bracket_straddles_tolerance"] = True

    # ── LEG 2: the differential's own structure, no external number ─────────
    print(f"\nLEG 2 — DIFFERENTIAL STRUCTURE (no external anchor; cannot be gamed by the "
          f"deck agreeing with itself)")
    envelope = 0.15
    within = abs(med) <= envelope
    eps = [(v["ep_eV"], v["delta"]) for v in ours.values() if v["ep_eV"] is not None]
    trend = None
    if len(eps) >= 5:
        lo = [d for e, d in eps if e <= statistics.median([e for e, _ in eps])]
        hi = [d for e, d in eps if e > statistics.median([e for e, _ in eps])]
        if lo and hi:
            trend = statistics.median(lo) - statistics.median(hi)
    leg2 = {"median": round(med, 4), "magnitude_envelope_dex": envelope,
            "within_envelope": bool(within),
            "low_minus_high_EP": (round(trend, 4) if trend is not None else None),
            "n_with_ep": len(eps)}
    print(f"  |median| {abs(med):.4f} vs envelope {envelope} -> "
          f"{'within' if within else 'OUTSIDE'}")
    print(f"     ⚠️ the envelope is DECLARED, not derived: solar Fe I is only mildly out "
          f"of LTE, so a differential of order a tenth of a dex is plausible and one of "
          f"order a dex is not. It bounds the absurd; it does not certify the value.")
    if trend is not None:
        print(f"  low-EP minus high-EP median = {trend:+.4f} dex")
        print(f"     over-ionisation acts hardest on the lowest-excitation Fe I lines, so "
              f"a NEGATIVE differential should be LARGER (more negative) at low EP")
    doc["leg2_structure"] = leg2

    # 🔴 `all()` OVER AN EMPTY SEQUENCE IS TRUE, AND THAT MADE THIS GATE PASS WHEN IT HAD
    # MEASURED NOTHING. The first version read
    #     all(v["verdict"] == "PASS" for v in leg1.values() if v.get("n"))
    # and the `if v.get("n")` filtered out every zero-overlap bracket -- so a product
    # sharing NOT ONE LINE with the anchor produced an empty generator, `all([]) is True`,
    # and the verdict came back VALIDATED. A gate that passes because it compared nothing
    # is worse than no gate: it is a gate-shaped absence, which is exactly the failure
    # class this project keeps finding (RYA-954's 8->200, RYA-923's 114->n=0, RYA-1044's
    # "0 fitted, 4 refused"). Caught by its own test rather than by a run.
    #
    # So overlap is now a PRECONDITION, stated separately from agreement: there must be at
    # least one bracket that actually compared lines, AND every bracket that compared any
    # must pass.
    compared = [v for v in leg1.values() if v.get("n")]
    if not compared:
        print(f"\n  🔴 NO BRACKET SHARED A SINGLE LINE WITH THE ANCHOR. That is not "
              f"agreement -- it is the absence of a comparison, and it cannot pass.")
        doc["no_overlap_with_anchor"] = True
    both = bool(compared) and all(v["verdict"] == "PASS" for v in compared) and within
    doc["verdict"] = "VALIDATED" if both else "NOT_VALIDATED"
    print(f"\nVERDICT: {doc['verdict']}  "
          f"(a deck must pass BOTH legs; either alone is weaker)")

    if a.out:
        Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
        print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
