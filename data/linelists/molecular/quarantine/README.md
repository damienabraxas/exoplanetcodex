# VALD3 molecular redistribution — QUARANTINED, NOT ADOPTED (RYA-1182)

`vald3_redistribution_quarantine_rya1182.csv` holds **8,977 molecular rows that were
sitting inside the atomic canonical store**, `data/linelists/canonical_gf.csv`, until
RYA-1182 moved them here.

| species | rows |
| --- | --- |
| CN  | 2766 |
| C2  | 2654 |
| CH  | 2265 |
| MgH |  594 |
| SiH |  583 |
| NH  |  114 |
| OH  |    1 |
| **total** | **8977** |

All 8,977 carry `seed_source=linelist(VALD3)`, `loggf_reference=VALD3`, `gf_tier=VALD3`,
and all 8,977 lie in **VIS only, 3780.19–6909.76 Å** — one VALD3 extract, redistributed.

## Why they are here and not in the RYA-1130 molecular store

RYA-1130 defines a provenance-first molecular transition schema
(`data/schemas/canonical_molecular_transition.schema.json`) that **requires**
`isotopologue`, `transition_identity`, `thermochemistry`, `normalization`, `precedence`
and `band_reach`. These rows carry none of it — they are a bare species/λ/EP/log gf
redistribution. Promoting them into that schema would mean inventing the fields it exists
to demand, so they are quarantined with their provenance intact instead.

**They are kept, not deleted** (RYA-946 keep-not-delete): every row is byte-identical to
the line that was removed from `canonical_gf.csv`, in the original order, under the same
26-column header. A reader can still find them; nothing adopts them.

## Why they must not live in `canonical_gf.csv`

`canonical_gf.csv` is the **atomic** laboratory-gf store. RYA-1136's CNO intake spent five
primary archives specifically avoiding VALD3 redistribution for molecular data — and these
rows sat one join away from being picked up by the next CNO measurement that reached for
canonical_gf. RYA-1142 check B2 found them; the consequence was live regardless of origin.

## The invariant that now holds

`pipeline/canonical_gf_invariants.py` states the atomic contract positively: **every
`canonical_gf.csv` row has an integer `key_z` in 1–92 and a non-empty `ion`.** That is a
closed rule — it does not depend on a list of known molecule names, so the next TiO or FeH
cannot slip in the way these did. Enforced by
`tests/test_canonical_gf_molecular_separation_rya1182.py`.

If a molecular transition ever *legitimately* belongs in the atomic store, the guard
requires it to carry vetted non-VALD3 provenance and to say so — see the invariant module.
