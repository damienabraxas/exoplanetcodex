#!/usr/bin/env python3
"""RYA-679 section 3E — re-run every affected adjudication under the ratified rule.

Re-fitting the syntheses is not required and would be misleading if done: the ratified
rule is a PURE FUNCTION of (dEW_dA, railed, red_chi2), all three of which are already
recorded per line in the committed result artefacts, and RYA-679 changed no fitting
maths (the only edit inside the chi2 was naming the literal 0.01 as SIGMA_FLUX_ASSUMED,
which is numerically identical — see `verify_fit_maths_unchanged`). Re-running
Turbospectrum would therefore reproduce the same inputs at some cost and some risk of
drift from an unrelated grid change, and would obscure rather than demonstrate the
effect of the rule change.

So this script does the honest thing: it reads the committed artefacts, applies BOTH
the old per-harness rule and the ratified rule to the same stored fit outputs, and
prints the before/after `reliable` table with any flip called out. That isolates the
decision, which is the entire question.

    python3 scripts/rya679_readjudicate.py [--json OUT]
"""
import argparse
import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
from solar_profile_fit import (RCHI2_REVIEW, RELIABLE_DEWDA,  # noqa: E402
                               SIGMA_FLUX_ASSUMED, assess_reliability)

RESULTS = os.path.join(_ROOT, "data", "results")

# The rule each harness applied BEFORE RYA-679, as a red_chi2 ceiling (None = no
# ceiling). Every harness also required (not railed) AND dEW_dA >= 40.0.
OLD_CEILING = {'Sr II': None, 'Zr II (plain)': None, 'Zr II (deblend)': 5.0,
               'Ba II': 60.0, 'Eu II': 15.0, 'Co I': 60.0}


def _old_rule(dewda, railed, rchi2, ceiling):
    if dewda is None or railed is None:
        return None
    ok = (not railed) and dewda >= 40.0
    if ceiling is not None and rchi2 is not None:
        ok = ok and rchi2 <= ceiling
    return bool(ok)


def _rows():
    """(species, line, arm, dEW_dA, railed, red_chi2) from the committed artefacts."""
    out = []

    def push(sp, line, arm, d):
        if not isinstance(d, dict) or 'red_chi2' not in d:
            return
        out.append((sp, line, arm, d.get('dEW_dA_mA_dex', d.get('dEW_dA')),
                    d.get('railed'), d.get('red_chi2')))

    p = os.path.join(RESULTS, "sr2_synthesis_rya551.json")
    if os.path.exists(p):
        for wl, v in json.load(open(p)).items():
            for arm in ('harps', 'iag'):
                push('Sr II', wl, arm, (v or {}).get(arm))

    p = os.path.join(RESULTS, "zr2_synthesis_rya560.json")
    if os.path.exists(p):
        for wl, v in json.load(open(p)).items():
            if wl == '_meta':
                continue
            for arm in ('harps', 'iag'):
                push('Zr II (plain)', wl, arm, (v or {}).get(arm))

    p = os.path.join(RESULTS, "zr2_deblend_rya585.json")
    if os.path.exists(p):
        for wl, v in json.load(open(p)).items():
            if wl == '_meta':
                continue
            for arm in ('harps', 'iag'):
                push('Zr II (deblend)', wl, arm, (v or {}).get(arm))

    p = os.path.join(RESULTS, "solar_ba_deblend_rya581.json")
    if os.path.exists(p):
        d = json.load(open(p))
        for arm, v in (d.get('per_arm') or {}).items():
            for hw, f in ((v or {}).get('deblended') or {}).items():
                if hw == '0.6':          # the headline half-window
                    push('Ba II', '5853.668', arm, f)

    p = os.path.join(RESULTS, "eu2_synthesis_rya565.json")
    if os.path.exists(p):
        d = json.load(open(p))
        for wl, v in (d.get('lines') or {}).items():
            for leg, legv in ((v or {}).get('legs') or {}).items():
                for arm in ('harps', 'iag', 'kp'):
                    push('Eu II', f"{wl}/{leg}", arm, (legv or {}).get(arm))

    p = os.path.join(RESULTS, "co_synthesis_rya564.json")
    if os.path.exists(p):
        for wl, v in json.load(open(p)).items():
            if wl.startswith('_'):
                continue
            for arm in ('harps', 'iag'):
                push('Co I', wl, arm, (v or {}).get(arm))
    return out


def verify_fit_maths_unchanged():
    """RYA-679 touched the chi2 only by naming its denominator. Assert the constant is
    what the literal was, so the before/after comparison is about the RULE alone."""
    assert SIGMA_FLUX_ASSUMED == 0.01, (
        "SIGMA_FLUX_ASSUMED moved away from the historical literal 0.01 — the committed "
        "red_chi2 values are no longer comparable and the syntheses MUST be re-run.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', help='write the table as JSON')
    args = ap.parse_args()

    verify_fit_maths_unchanged()
    print("RYA-679 section 3E — re-adjudication under the ratified reliability rule")
    print(f"  ratified: reliable = (not railed) AND dEW_dA >= {RELIABLE_DEWDA} mA/dex")
    print(f"  red_chi2: REPORTED, review-flagged above {RCHI2_REVIEW}, NEVER gating")
    print(f"  sigma_flux assumed by every chi2: {SIGMA_FLUX_ASSUMED}\n")

    hdr = (f"{'species':16s} {'line':>10s} {'arm':5s} {'dEW/dA':>8s} {'railed':>6s} "
           f"{'rchi2':>9s} {'old bar':>8s} {'before':>7s} {'after':>6s} {'review':>7s}")
    print(hdr)
    print('-' * len(hdr))

    rows, flips = [], []
    for sp, line, arm, dewda, railed, rchi2 in _rows():
        ceiling = OLD_CEILING[sp]
        before = _old_rule(dewda, railed, rchi2, ceiling)
        new = assess_reliability(dict(dEW_dA=dewda, railed=railed, red_chi2=rchi2))
        after, review = new['reliable'], new['rchi2_review']
        flag = '  FLIP' if before is not None and before != after else ''
        print(f"{sp:16s} {line:>10s} {arm:5s} {str(dewda):>8s} {str(railed):>6s} "
              f"{str(rchi2):>9s} {str(ceiling):>8s} {str(before):>7s} {str(after):>6s} "
              f"{str(review):>7s}{flag}")
        rec = dict(species=sp, line=line, arm=arm, dEW_dA=dewda, railed=railed,
                   red_chi2=rchi2, old_ceiling=ceiling, reliable_before=before,
                   reliable_after=after, rchi2_review=review)
        rows.append(rec)
        if before is not None and before != after:
            flips.append(rec)

    print(f"\nFLIPS: {len(flips)}")
    for f in flips:
        print(f"  {f['species']} {f['line']} {f['arm']}: "
              f"{f['reliable_before']} -> {f['reliable_after']}")
    if not flips:
        print("  none — no disposition moves under the ratified rule.")

    rev = [r for r in rows if r['rchi2_review']]
    live = [r for r in rev if r['reliable_after']]
    print(f"\nreviewed rows: {len(rev)} of {len(rows)} exceed the review trigger; "
          f"{len(rev) - len(live)} of those are already not reliable (they fail "
          f"sensitivity or are railed, so the flag adds nothing).")
    print(f"\nACTIONABLE — RELIABLE but the in-window model does not reproduce the "
          f"window ({len(live)}):")
    for r in sorted(live, key=lambda x: -x['red_chi2']):
        rms = (r['red_chi2'] ** 0.5) * SIGMA_FLUX_ASSUMED
        print(f"  {r['species']:16s} {r['line']:>10s} {r['arm']:5s} "
              f"red_chi2={r['red_chi2']:>7} (residual RMS {rms * 100:.1f}% of continuum)")
    if not live:
        print("  none")

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(ticket='RYA-679', rule=dict(
                reliable_dEW_dA_floor=RELIABLE_DEWDA,
                red_chi2_review_trigger=RCHI2_REVIEW,
                red_chi2_gates_reliable=False,
                sigma_flux_assumed=SIGMA_FLUX_ASSUMED),
                old_ceilings=OLD_CEILING, rows=rows, flips=flips), fh, indent=2)
        print(f"\nWrote {args.json}")


if __name__ == '__main__':
    main()
