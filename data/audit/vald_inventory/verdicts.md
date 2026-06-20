# RYA-376 — VALD line-list inventory: coverage matrix + near-IR readiness

Read-only audit of all VALD assets under `data/linelists/`.
Regenerate: `python scripts/audit_vald_inventory_rya376.py` → `assets.csv`, `coverage_matrix.csv`.
Bands (air Å): optical <9500 · Y 9500–11000 · J 11000–14000 · H 14000–18000 · K 18000–25000.

## Headline verdict — NEAR-IR EXTRACTION NEEDED (solar) + K-band is quarantine-grade

- **Y/J/H (0.95–1.80 µm) atomic VALD data IS held, HFS-on**, for **55 Cnc A, α Cen A, α Cen B** (per-star lists `linelist_{55cnc,alpha_cen_a,alpha_cen_b}.csv`, 1150–17000 Å) and **Procyon to 1.10 µm (Y only)**. The named diagnostics — **P I 1.05 µm (Y), S I 1.045 µm (Y), K I 1.17 µm (J), near-IR C I + Fe-peak, OH 1.5–1.8 µm (H)** — are present for those stars.
- **K-band (1.80–2.50 µm, the CO-overtone region) atomic context exists ONLY in `linelist_full.csv`, sourced from the HFS-OFF quarantine** (see defect 1). No HFS-on K-band anywhere.
- **Solar has ZERO near-IR — optical only (3780–6910 Å).** Since the reflected-asteroid IR (Vesta, RYA-373) is **reflected sunlight**, it needs a **solar-params** NIR atomic list → **does not exist → extraction needed.** (The 2.3 µm CO overtone itself is *molecular* — covered by the RYA-236 CO_IR list, not VALD atomic — so VALD-NIR is needed only to validate **atomic** diagnostics beyond CO.)

## Per-delivery intake verdicts (raw VALD3 deliveries)

| delivery | star | req Å | actual hi Å | n lines | trunc | >100k cap | HFS | bands | verdict |
|---|---|---|---:|---:|---|---|---|---|---|
| vald_solar_raw | solar | 3780–6910 | 6910 | 108,969 | no | **yes (complete)** | ? | optical | ACCEPT |
| vald_55cnc_raw | 55cnc | 3780–6910 | 6910 | 125,615 | no | **yes (complete)** | ? | optical | ACCEPT |
| vald_55cnc_nir_raw | 55cnc | 6910–17000 | 16997 | 5,218 | no | no | on* | Y/J/H | ACCEPT |
| vald_55cnc_nir_5k30k_hfsoff_quarantine | 55cnc | 5000–30000 | 29995 | 21,312 | no | no | **off** | Y/J/H/**K** | **QUARANTINE** |
| vald_55cnc_uv_{a,b,raw} | 55cnc | 1150–3780 | 3780 | 24.5k/90.5k/77.5k | no | no | ? | UV | ACCEPT |
| vald_55cnc_uv_019509_hfsoff_quarantine | 55cnc | 1150–3780 | 3780 | 77,458 | no | no | **off** | UV | **QUARANTINE** |
| vald_alpha_cen_a_nir_raw | αCenA | 6910–17000 | 16997 | 3,707 | no | no | ? | Y/J/H | ACCEPT |
| vald_alpha_cen_b_nir_raw | αCenB | 6910–17000 | 16997 | 5,152 | no | no | ? | Y/J/H | ACCEPT |
| vald_alpha_cen_{a,b}_{optical,uv1,uv2}_raw | αCen | 1150–6910 | — | — | no | no | ? | UV/opt | ACCEPT |
| vald_procyon_nir_hfson_raw | procyon | 6910–11000 | 10984 | 1,147 | no | no | **on** | Y only | ACCEPT |
| vald_procyon_{optical,uv}_hfson_raw | procyon | 1150–6910 | — | — | no | no | on | UV/opt | ACCEPT |

17 ACCEPT, 2 QUARANTINE (HFS-off). **No delivered file is truncated.** Truncation-trap note: the historical full-range extraction **019387 (1150–30000 Å) HIT the 100k cap** (per README) and was re-extracted in sub-ranges. Two *optical* deliveries (solar 108,969; 55 Cnc 125,615) **exceed 100k yet delivered complete** (n_data == n_selected, no warning) — the web cap is therefore not a hard universal 100k here; **flag for the curator to confirm no silent server-side cap on future large extractions.**

## Assembled line lists

| list | star | range Å | n | NIR Y/J/H/K | note |
|---|---|---|---:|---|---|
| linelist_solar | solar | 3780–6910 | 108,971 | 0/0/0/0 | **optical-only** |
| linelist_master | 55cnc | 3780–6910 | 125,617 | 0/0/0/0 | optical-only |
| linelist_55cnc | 55cnc | 1150–16997 | 245,878 | 419/957/1790/0 | UV+opt+NIR→H, **HFS-on** |
| linelist_alpha_cen_a | αCenA | 1150–16997 | 141,940 | 323/656/1259/0 | UV+opt+NIR→H, HFS-on |
| linelist_alpha_cen_b | αCenB | 1150–16997 | 155,900 | 411/951/1753/0 | UV+opt+NIR→H, HFS-on |
| linelist_procyon | procyon | 1150–10984 | 108,563 | 209/0/0/0 | UV+opt+NIR→Y only |
| **linelist_full** | 55cnc | 3780–29995 | 140,483 | 1217/1442/3677/3032 | opt + **NIR→K from the HFS-OFF quarantine** (defect 1) |
| canonical_gf (RYA-353) | — | 3780–9199 | 145,886 | 0/0/0/0 | **single-source gf is optical-only** (defect 2) |

## Provenance defects (RYA-332 single-source)

1. **`linelist_full.csv` NIR = the HFS-OFF quarantined extraction.** Its NIR band counts (Y/J/H/K = 1217/1442/3677/3032) are **identical** to `vald_55cnc_nir_5k30k_hfsoff_quarantine.txt`, *not* the HFS-on `vald_55cnc_nir_raw` (419/957/1790/0) its README claims. So a curator-quarantined (HFS-off) source leaked into the assembled "full" list — every NIR line, incl. all K-band. **HFS-off mis-partitions log_gf for odd-Z hyperfine elements (Mn, Co, Cu, V, Sc, Na, Al, K)** → `linelist_full` NIR is not science-ready for those. README provenance is stale/wrong.
2. **`canonical_gf.csv` (the RYA-353 single-source gf table) stops at 9199 Å** → **all NIR gf is outside the single-source table** and outside the RYA-355 stewardship guard. Divergent-gf instance: the same NIR transition carries different `log_gf` HFS-on (per-star lists) vs HFS-off (`linelist_full`) — the HFS partition vs the merged value — an un-reconciled duplication.

## Scoped extraction (NEEDED), per the codex-vald-extraction recipe

- **Solar** (the gap that blocks the reflected-solar Vesta IR): Extract Stellar at Teff 5772 / logg 4.44 / [M/H] 0.0 / vmic 1.0, **HFS splitting ON**, 27 elements (`elements_master.json`), **0.95–2.5 µm**. Split to dodge the cap: **9500–14000 (Y+J), 14000–18000 (H), 18000–25000 (K)** — finer if any sub-range > ~100k.
- **K-band HFS-on** for 55 Cnc to replace the quarantined HFS-off K data now in `linelist_full` (18000–25000 Å, HFS-on).
- **Extend Procyon NIR** beyond 1.10 µm (currently Y-only) to J/H/K if Procyon IR is wanted.
- After extraction: route NIR gf into `canonical_gf.csv` so the single-source guard (RYA-353/355) covers the NIR.
