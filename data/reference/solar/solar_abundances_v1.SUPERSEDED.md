# solar_abundances_v1.csv — SUPERSEDED by v2 (RYA-522)

`solar_abundances_v1.csv` is retained **immutable** (RYA-469 write-once) as a historical
record. It is **no longer the authoritative solar gold reference** — the `CURRENT`
pointer names **v2**.

**Reason superseded:** v1 carried `C, I = 10.26` (XH +1.8) — the saturated **C I 5380.337**
EW artifact (RYA-520): a saturated line on the flat curve-of-growth turned a tiny EW error
into a ~1.8-dex abundance error, and the RYA-469 freeze had sourced the gold from the raw
EW channel. v2 is regenerated **entirely from the phase_c verdict channel** (RYA-521,
synthesis + HFS + atlas + NLTE + curation), tiered by row-confidence (RYA-522): **C = 8.491**
(CH G-band synthesis).

See `docs/audit/solar_gold_v2_ratification_rya522.md` for the full v2-vs-v1-vs-Asplund
ratification table.

**Do not read v1 as authoritative.** Consumers read `read_solar_reference('CURRENT')` → v2.
