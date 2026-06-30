# RYA-376 — VALD line-list inventory: coverage matrix + near-IR readiness

Read-only audit of all VALD assets under `data/linelists/`.
Regenerate: `python scripts/audit_vald_inventory_rya376.py`
→ `assets.csv`, `coverage_matrix.csv`, `ir_triage.csv`, `groupA_depth_confirm.csv`.
Bands (air Å): optical <9500 · Y 9500–11000 · J 11000–14000 · H 14000–18000 · K 18000–25000.

> **Re-baselined 2026-06-30.** The original audit (cb43ef4, 2026-06-20) reported
> *NEAR-IR EXTRACTION NEEDED (solar)*. Since then **RYA-381** (assemble non-optical
> solar VALD → `linelist_solar` 1150–25000 Å) and **RYA-387** (re-extract the solar
> wings at central-depth 0.001, HFS-on) **closed that gap**. This re-run reflects the
> current holdings AND adds the ticket's two new deliverables: the **27-element IR
> triage** (IR-PRIMARY / IR-CROSSCHECK / IR-NONE) and the **Group-A depth confirmation**.

## Headline verdict — NEAR-IR READY (solar), HFS-on, full Y/J/H/K to 2.5 µm

- **Solar now holds the complete near-IR atomic list, HFS-on.** `linelist_solar.csv`
  spans **1150–24998 Å** (354,496 lines) with **NIR Y/J/H/K = 3600 / 4966 / 10207 / 7988**,
  assembled from five HFS-on solar wing deliveries (`vald_solar_{fuv,nearuv,redopt,
  ir_9500_17000,ir_17000_25000}_hfson`). **The reflected-asteroid IR (Vesta = reflected
  sunlight, RYA-373) now has a solar-params NIR atomic list — the gap the prior audit
  flagged is CLOSED.** (The 2.3 µm CO overtone itself is *molecular* — RYA-236 CO_IR list,
  not VALD atomic; VALD-NIR validates the *atomic* diagnostics beyond CO.)
- **Y/J/H (0.95–1.80 µm) HFS-on atomic data is also held** for **55 Cnc A, α Cen A,
  α Cen B** (per-star lists to 1.70 µm) and **Procyon to 1.10 µm (Y only)**.
- **Residual extraction items remain (see below):** 55 Cnc **K-band** is still
  quarantine-grade (HFS-off, defect 1); **Procyon** is Y-only (no J/H/K); and the
  **single-source gf table is optical-only** (defect 2) so NIR gf is unguarded.

## 27-element IR triage (IR-PRIMARY / IR-CROSSCHECK / IR-NONE)

Each target element tagged by the VALUE the IR adds, then CONFIRMED against actual
holdings (full table: `ir_triage.csv`). All Group-A/B elements have NIR coverage held;
Group-C heavies are barren or low-value in the IR — confirming they belong to the
UV/near-blue leg, not the IR extraction.

| group | tag | elements | confirmation |
|---|---|---|---|
| **A** | **IR-PRIMARY** | C, N, O, S, P, K, Na, Mg, Al | all CONFIRMED — NIR held + diagnostics present (below) |
| **B** | **IR-CROSSCHECK** | Si, Ca, Ti, Cr, Mn, V, Co, Ni, Sc, Cu, **Fe**, **Zn** | all CONFIRMED — NIR held; optical owns value |
| **C** | **IR-NONE** | Ba, Y, Zr, Sr, Eu, Li | CONFIRMED — IR barren/low-value → UV/near-blue leg |

- **IR-PRIMARY (Group A) — IR is a new or best channel.** S I 1.045 µm ≫ the lone
  weak optical S I 6757 (RYA-488); P I 1.05 µm is essentially the ONLY real P channel;
  K I 1.17 µm sidesteps the O₂ A-band telluric on optical 7699; Na/Mg/Al near-IR
  subordinate lines are UNSATURATED where the optical resonance lines saturate (the
  RYA-465 problem); C/N/O flagship (C I near-IR, N I red + NH, O I 777 + OH).
- **IR-CROSSCHECK (Group B) — optical owns the value, IR confirms.** Deliberate
  exception **Fe = the IR SIGMA RULER**: Fe carries **7581 NIR lines — by far the most
  of any atomic element** (H-band alone 2968, of which 845 with central depth >0.05),
  so it is the right element to characterize the IR scatter/systematic, then apply that
  to the sparse-line CNOPS measurements. **Zn is the 27th element, not named in the
  ticket triage** → classified IR-CROSSCHECK (optical Zn I 4722/4810/6362 owns the
  value; IR lines are subordinate).
- **IR-NONE (Group C) — IR does not help.** The neutron-capture heavies (Ba, Y, Zr,
  Sr, Eu) live in the blue/near-UV as ionized lines; the few NIR lines the solar list
  carries for them are sparse/weak (`nir_held` true but low-value). Li stays optical
  (6707). Group C confirms the UVES Ceres NUV leg is where the heavies belong.

## Group-A IR diagnostic depth confirmation (held solar list)

Reach + depth come from `linelist_solar.csv`, never from memory (full: `groupA_depth_confirm.csv`).

| diagnostic | species | target Å | band | matched Å | depth | note |
|---|---|---:|---|---:|---:|---|
| S I 1.045 µm triplet | S I | 10455.45 | Y | 10455.470 | **0.281** | far stronger than optical S I 6757 (~0.01–0.03) |
| S I 1.045 µm triplet | S I | 10456.76 | Y | 10456.790 | 0.173 | |
| S I 1.045 µm triplet | S I | 10459.41 | Y | 10459.460 | 0.249 | |
| P I 1.05 µm quartet | P I | 10511.59 | Y | 10511.588 | 0.023 | intrinsically weak — the ONLY real P channel |
| P I 1.05 µm quartet | P I | 10529.52 | Y | 10529.524 | 0.045 | |
| P I 1.05 µm quartet | P I | 10581.58 | Y | 10581.577 | **0.062** | deepest P line |
| P I 1.05 µm quartet | P I | 10596.90 | Y | 10596.903 | 0.020 | |
| K I 1.17 µm doublet | K I | 11690.22 | J | 11690.219 | **0.292** | HFS-split (6 comps); dodges 7699 O₂ telluric |
| K I 1.17 µm doublet | K I | 11769.64 | J | 11769.638 | 0.152 | HFS-split (9 comps) |
| C I near-IR | C I | 10691.25 | Y | 10691.250 | **0.325** | |
| O I 777 triplet | O I | 7771.94 | optical | 7771.944 | 0.166 | |
| O I near-IR | O I | 11302.38 | J | 11302.378 | 0.019 | weak |
| Mg I 1.57 µm | Mg I | 15740.71 | H | 15740.706 | **0.479** | strong unsaturated subordinate (RYA-465) |
| Al I 1.67 µm | Al I | 16718.96 | H | 16718.911 | **0.471** | strong unsaturated subordinate, HFS (6 comps) |

**14/14 named diagnostics covered.** Molecular bands cited by the ticket — **OH 1.5–1.8 µm,
NH 3360, CO 2.3 µm** — are NOT in the VALD atomic list and are out of this audit's scope
(handled by the molecular lists, e.g. RYA-236 CO_IR).

## Per-delivery intake verdicts — 24 raw deliveries: 22 ACCEPT, 2 QUARANTINE

| delivery | star | req Å | actual hi Å | n lines | trunc | >100k cap | HFS | bands | verdict |
|---|---|---|---:|---:|---|---|---|---|---|
| vald_solar_raw | solar | 3780–6910 | 6910 | 108,969 | no | **yes (complete)** | ? | optical | ACCEPT |
| vald_solar_fuv_1150_2000_hfson | solar | 1150–2000 | 2000 | 42,886 | no | no | on | UV | ACCEPT |
| vald_solar_nearuv_2000_3780_hfson | solar | 2000–3780 | 3780 | 161,526 | no | **yes (complete)** | on | UV | ACCEPT |
| vald_solar_redopt_6910_9500_hfson | solar | 6910–9500 | 9500 | 14,352 | no | no | on | optical | ACCEPT |
| vald_solar_ir_9500_17000_hfson | solar | 9500–17000 | 16998 | 16,793 | no | no | **on** | Y/J/H | ACCEPT |
| vald_solar_ir_17000_25000_hfson | solar | 17000–25000 | 24998 | 9,968 | no | no | **on** | H/**K** | ACCEPT |
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

**No delivered file is truncated.** Truncation-trap note: the historical full-range
extraction **019387 (1150–30000 Å) HIT the 100k cap** (per README) and was re-extracted
in sub-ranges. **Three deliveries exceed 100k yet delivered complete** (solar 108,969;
solar near-UV 161,526; 55 Cnc 125,615 — all `n_data == n_selected`, no warning): the web
cap is therefore **not a hard universal 100k** here — flag for the curator to confirm no
silent server-side cap on future large extractions.

## Assembled line lists

| list | star | range Å | n | NIR Y/J/H/K | note |
|---|---|---|---:|---|---|
| linelist_solar | solar | 1150–24998 | 354,496 | 3600/4966/10207/7988 | **full UV+opt+NIR→K, HFS-on** (RYA-381/387) |
| linelist_master | 55cnc | 3780–6910 | 125,617 | 0/0/0/0 | optical-only |
| linelist_55cnc | 55cnc | 1150–16997 | 245,878 | 419/957/1790/0 | UV+opt+NIR→H, HFS-on |
| linelist_alpha_cen_a | αCenA | 1150–16997 | 141,940 | 323/656/1259/0 | UV+opt+NIR→H, HFS-on |
| linelist_alpha_cen_b | αCenB | 1150–16997 | 155,900 | 411/951/1753/0 | UV+opt+NIR→H, HFS-on |
| linelist_procyon | procyon | 1150–10984 | 108,563 | 209/0/0/0 | UV+opt+NIR→Y only |
| **linelist_full** | 55cnc | 3780–29995 | 140,483 | 1217/1442/3677/3032 | opt + **NIR→K from the HFS-OFF quarantine** (defect 1) |
| canonical_gf (RYA-353) | — | 3780–9199 | 145,886 | 0/0/0/0 | **single-source gf is optical-only** (defect 2) |

## Provenance defects (RYA-332 single-source)

1. **`linelist_full.csv` NIR = the HFS-OFF quarantined extraction (STILL OPEN, 55 Cnc).**
   Its NIR band counts (Y/J/H/K = 1217/1442/3677/3032) are **identical** to
   `vald_55cnc_nir_5k30k_hfsoff_quarantine.txt`, *not* the HFS-on `vald_55cnc_nir_raw`
   (419/957/1790/0) its README claims. A curator-quarantined (HFS-off) source leaked into
   the assembled "full" list — every K-band line included. **HFS-off mis-partitions
   log_gf for odd-Z hyperfine elements (Mn, Co, Cu, V, Sc, Na, Al, K)** → `linelist_full`
   NIR is not science-ready for those. README provenance is stale/wrong. **Note:** the new
   *solar* K-band (`vald_solar_ir_17000_25000_hfson`) is HFS-**on** and is NOT affected —
   the defect is now scoped to 55 Cnc's `linelist_full` K-band only.
2. **`canonical_gf.csv` (the RYA-353 single-source gf table) stops at 9199 Å (STILL OPEN).**
   **All NIR gf is outside the single-source table** and outside the RYA-355 stewardship
   guard — now more pressing because the solar NIR list is in active use. Divergent-gf
   instance: the same NIR transition carries different `log_gf` HFS-on (per-star/solar
   lists) vs HFS-off (`linelist_full`) — an un-reconciled duplication.

## Residual extraction NEEDED (hand to RYA-427), per the codex-vald-extraction recipe

- **55 Cnc K-band HFS-on** (18000–25000 Å, HFS-on) to replace the quarantined HFS-off
  K data now in `linelist_full`.
- **Extend Procyon NIR** beyond 1.10 µm (currently Y-only) to J/H/K if Procyon IR is wanted.
- **Route NIR gf into `canonical_gf.csv`** so the single-source guard (RYA-353/355) covers
  the NIR (closes defect 2).
- α Cen A/B already reach H (1.70 µm); extend to K only if K-band α Cen IR is scoped.

## Element-list double-check (the ticket's second ask)

Confirmed against actual holdings + line depth: **Group A reach is real** — every named
IR-PRIMARY diagnostic (S I 1.045, P I 1.05, K I 1.17, C I near-IR, O I, Mg I 1.57,
Al I 1.67) is present in the held solar list with sane depth. **Group C belongs on the
UV/near-blue leg** — the heavies are barren/low-value in the IR. **Fe is the IR sigma
ruler** (7581 NIR lines). The 27th element **Zn** (unlisted in the ticket) is tagged
IR-CROSSCHECK.

*Audit only — no merge.*
