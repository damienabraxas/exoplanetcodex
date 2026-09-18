"""
Vintage supersession for band-keyed dA/dxi entries -- RYA-1225.

🔴 THE DEFECT. A dA/dxi is a property of a SYNTHESIS, not just of a line pool. When the
synthesis for a band changes, every derivative measured before that change describes a
model the products no longer use -- while still looking like a live measurement, because
nothing in the entry says which vintage it belongs to.

That is not hypothetical. RYA-1207 (merged 2026-09-09) turned molecular opacity on for the
near-UV band -- `config/synth_bands.yaml` sets `use_molecules: true` on near-UV ALONE --
and near-UV Fe II's published slope disagreed 2.0-2.8x between the run before it
(RYA-1168, 2026-09-03, dA/dxi = -0.0475) and the run after it (RYA-1213, 2026-09-12,
-0.1100), on a byte-identical 12-line pool. Red-optical, whose synthesis did not move,
agreed to 1.03-1.08x across the same two runs.

🔴 STALENESS IS MEASURED, NEVER LISTED. A run that records the feed `A` it measured against
carries its own expiry date: if that `A` is not the product's `A` today, the derivative was
taken on a superseded synthesis. There is no hardcoded ticket list here, so a future
synthesis change is caught by the same rule rather than needing this module edited.

⚠️ AND A SUPERSEDED ENTRY IS NOT AUTOMATICALLY REPLACEABLE. Substituting another tier's
slope is only honest when the two tiers are the SAME POOL, and that is a measured fact
owned by RYA-1213's `same_pool_ledger.json` (delta_dex == 0 AND same_n), not an assumption
from equal `n_lines`. Where the ledger does not certify identity the entry is REPORTED and
left alone -- near-UV Fe I is exactly that case (its pools differ by RYA-1209's line drop,
delta 0.002-0.007, same_n False), and silently swapping there would move published bars on
a pool nobody proved was the same.
"""

from __future__ import annotations

from typing import Any

#: Tolerance on the abundance comparison. Feed A is published to 3-4 dp, so anything
#: above this is a real difference and not a rounding artifact.
A_TOL_DEX = 1e-9


def certified_same_pool(ledger: dict) -> set[tuple[str, str, str, str]]:
    """(ion, holding, treatment, band) whose two tiers RYA-1213 MEASURED as one pool.

    The ledger's own reading: "a row with delta 0.000 AND the same n is the same-pool
    property holding exactly." Both halves are required -- equal n alone does not make
    two pools the same pool (RYA-1204's 40-line and 58-line pools shared 2 lines).
    """
    out = set()
    for r in ledger.get("rows", []):
        if r.get("delta_dex") == 0.0 and r.get("same_n") is True:
            ion = str(r.get("ion", "")).split()[-1]
            out.add((ion, r["holding"], r["treatment"], r["band"]))
    return out


def is_stale(entry: dict, product: dict | None) -> bool:
    """Was this derivative measured against an abundance the product no longer has?

    Only entries that RECORD the `A` they ran against can be judged. An entry without one
    is not assumed fresh and is not assumed stale -- it is simply not judgeable here, and
    says so by returning False rather than guessing.
    """
    if product is None or entry.get("A") is None:
        return False
    return abs(float(product["A"]) - float(entry["A"])) > A_TOL_DEX


#: 🔴 DETECT BROADLY, ACT NARROWLY. The staleness rule is general on purpose -- it finds
#: every band-keyed entry measured against a superseded abundance. Acting on all of them
#: is a different matter: substituting a slope changes a PUBLISHED uncertainty, and that
#: is only defensible where the divergence has actually been adjudicated and its physical
#: cause written down. RYA-1225 adjudicated near-UV Fe II (cause: RYA-1207's near-UV
#: molecular opacity). The same rule also flags stale red-optical and NIR entries, whose
#: cause is NOT established here -- RYA-1191's telluric re-run moved those pools (RYA-515),
#: which is a different mechanism. Those are REPORTED as owed and left untouched, because
#: swapping them would be this ticket quietly re-pricing bars nobody diagnosed.
ADJUDICATED: frozenset[tuple[str, str]] = frozenset({("II", "near-UV")})


def apply_vintage_supersession(
    band_idx: dict[tuple, dict],
    products: list[dict],
    ledger: dict,
    adjudicated: frozenset[tuple[str, str]] = ADJUDICATED,
) -> tuple[dict[tuple, dict], list[dict]]:
    """Replace superseded band entries with the same-pool measurement taken after the change.

    Returns the (possibly rewritten) index and a report row per superseded-or-flagged
    entry. Nothing is mutated in place: the caller gets a new dict.
    """
    live = {(p["ion"], p["holding"], p["tier"], p["treatment"], p["band"]): p
            for p in products}
    certified = certified_same_pool(ledger)
    out: dict[tuple, dict] = dict(band_idx)
    report: list[dict] = []

    for key, entry in band_idx.items():
        ion, holding, tier, treatment, band = key
        if not is_stale(entry, live.get(key)):
            continue
        product = live[key]
        row = {
            "key": {"ion": ion, "holding": holding, "tier": tier,
                    "treatment": treatment, "band": band},
            "stale_entry": {"dA_dxi": entry.get("dA_dxi"), "A_it_ran_against": entry.get("A"),
                            "source": entry.get("_source"), "ticket": entry.get("_ticket")},
            "live_product_A": product["A"],
        }

        if (ion, band) not in adjudicated:
            row["action"] = "FLAGGED_NOT_ADJUDICATED"
            row["reason"] = ("measured against a superseded abundance, but no RCA has "
                             "established the cause for this band/ion, so no slope is "
                             "substituted here -- reported as owed, value left untouched")
            report.append(row)
            continue

        if (ion, holding, treatment, band) not in certified:
            row["action"] = "FLAGGED_NOT_SUPERSEDED"
            row["reason"] = ("the same-pool ledger does not certify this cell (delta_dex "
                             "!= 0 or same_n false), so no other tier's slope may stand "
                             "in for it -- reported as owed, value left untouched")
            report.append(row)
            continue

        #: The donor must be the SAME pool in another tier, and must not itself be stale.
        donor_key = next(
            (k for k in band_idx
             if k[0] == ion and k[1] == holding and k[3] == treatment and k[4] == band
             and k[2] != tier
             and not is_stale(band_idx[k], live.get(k))
             and band_idx[k].get("xi_state") == "MEASURED"
             and band_idx[k].get("dA_dxi") is not None),
            None)
        if donor_key is None:
            row["action"] = "FLAGGED_NO_DONOR"
            row["reason"] = ("certified same-pool, but no contemporaneous MEASURED entry "
                             "exists in another tier to supersede it")
            report.append(row)
            continue

        donor = band_idx[donor_key]
        out[key] = {**donor, "tier": tier,
                    "_superseded": {
                        "replaced_ticket": entry.get("_ticket"),
                        "replaced_source": entry.get("_source"),
                        "replaced_dA_dxi": entry.get("dA_dxi"),
                        "replaced_A_it_ran_against": entry.get("A"),
                        "donor_tier": donor_key[2],
                        "cause": ("measured before the band's synthesis changed; the "
                                  "recorded A is not the product's A today (RYA-1225)"),
                    }}
        row.update({"action": "SUPERSEDED",
                    "donor": {"tier": donor_key[2], "dA_dxi": donor.get("dA_dxi"),
                              "ticket": donor.get("_ticket"),
                              "source": donor.get("_source")},
                    "ratio_new_over_old": None if not entry.get("dA_dxi") else round(
                        abs(donor["dA_dxi"] / entry["dA_dxi"]), 3)})
        report.append(row)

    return out, report
