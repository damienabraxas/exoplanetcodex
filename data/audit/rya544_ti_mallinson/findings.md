# RYA-544 — Mallinson-2024 Ti grid: acquired + derived, gate STOP (Ti I gf-scale is the blocker)

**Date:** 2026-07-10 · **Branch:** `ryandamienschmitt/rya-544-...` (off origin/main 73088bc, register v18) · **NO MERGE.**
Script: `scripts/rya544_ti_mallinson.py`. Runs on Sirius (PySME/MARCS, venv_pysme).

## What this ticket did (per the RYA-544 spec-rewrite comment)

RYA-542 adjudication: the correct solar Ti I NLTE correction is ~+0.04 (Mallinson 2022/2024 + Sitnova
2020, ab-initio H collisions), NOT our Engine-A MAFAGS-OS +0.108 or Engine-B MARCS +0.20 (both on the
outdated Bergemann-2011 scaled-Drawin atom). RYA-544 wires the *grid* (option 1, infra-consistent) and
gates on OUR ionization balance — reference-blind, NOT "matches +0.04".

## 1. Grid acquired (RYA-540 persistent M.2 cache)

Mallinson-2024 Ti grid from **Zenodo 10753497** (DOI 10.5281/zenodo.10753497), PySME `.grd`:
- `nlte_Ti_pysme.grd` 11.38 GB, **md5 80ad1bb5ac676910de8bb905eee4bca1** (verified)
- `atmos_Ti.txt` (MARCS node table), `label_Ti.txt` (NIST level table)
in `/srv/codex/grids/nlte/amarsi_galah/` (M.2), md5-pinned + retained. Prov:
`data/nlte_grids/amarsi_galah/nlte_Ti_pysme.grd.prov.json`. 587-level atom (459 Ti I + 127 Ti II +
Ti III), Grumer & Barklem 2020 ab-initio H collisions; 3756 MARCS models, solar node well inside.

**Reader gotcha:** our Amarsi-format `.grd` reader (`read_amarsi_grid`) MIS-PARSES the Mallinson
grid's internal level table (same v1.10 container, different layout → 813 garbage levels). Level
labels are read from the ASCII `label_Ti.txt` (energy-matched, grid-native NIST); PySME's native
`set_nlte` reads the departures for synthesis (works).

## 2. OUR Ti I NLTE correction — DERIVED from the grid (+0.0506) ✓

`scripts/rya544_ti_mallinson.py --derive` (PySME NLTE-vs-LTE COG, MARCS solar node 5772/4.44/0.0/1.0,
A(Ti)=4.97), departures engaged:

| Ti I line | delta = A_NLTE − A_LTE |
|---|---|
| 5689.460 | +0.0503 |
| 5648.565 | +0.0542 |
| 5662.150 | +0.0506 |
| **median** | **+0.0506** |

OUR +0.0506 reproduces Mallinson-2024's reported solar Ti I correction **+0.052** (and Mallinson-2022
+0.03 / Sitnova-2020 +0.03) — **validate-don't-tune, firewall held**: derived from the grid, NOT
fitted to +0.04. This RESOLVES the RYA-542 NLTE-correction question: the correct ab-initio value is
~+0.05; our old engines' +0.108 (MAFAGS-OS) and +0.20 (MARCS) were both inflated by the outdated
Bergemann-2011 scaled-Drawin atom. **The NLTE correction is a differential (gf cancels) → this value
is robust regardless of §3.**

## 3. Ionization-balance ACCEPTANCE GATE → **STOP** (reference-blind)

`--gate`: LTE EW-invert OUR solar Ti I + Ti II pools (measured EWs `sol_ew_results_v1.csv`, canonical
gf) on the SAME MARCS atmosphere (weak/moderate lines, EW 5–60 mA = linear COG, well-conditioned;
railed inversions dropped, not clamped), apply the delta, require |A(Ti I)_NLTE − A(Ti II)| <
`FE_IONISATION_GATE` (0.05).

| ion | n | A_LTE | MAD | note |
|---|---|---|---|---|
| Ti I | 9/19 | **5.598** | **0.508** | half-dex line-to-line scatter |
| Ti II | 3/3 | **5.008** | 0.032 | clean, ≈ solar (Asplund 4.97) |

- A(Ti I)_NLTE = 5.598 + 0.0506 = **5.649**;  A(Ti II)_NLTE = **5.006**
- **NLTE ionization balance = +0.643  ≫  0.05  →  STOP.**

**Diagnosis — the blocker is OUR Ti I gf-scale, NOT the NLTE grid:**
- Ti I MAD = **0.508 dex**: real abundances don't scatter half a dex; a bad/ungraded gf-scale does.
  This independently reproduces the RYA-521/v1 finding (Ti I `curation_residual +0.502`, "thin graded
  pool / LOW_CONFIDENCE"). Many Ti I canonical gf are K10 synth-gf, single-source, ungraded.
- Ti II is clean (MAD 0.03) and sits at the correct solar value → the reliable ion.
- The NLTE correction (+0.0506) is ~13× smaller than the imbalance (+0.643) — it cannot be the cause,
  and (being differential) is unaffected by the gf-scale.

## 4. Verdict + actions (STOP branch)

- **Do NOT wire** the Mallinson grid into production. `constants.py` Ti stays `Ti_Bergemann2011_MPIA.csv`;
  register Ti stays **CHECK**; the RYA-534 Ti strict-xfail is **NOT flipped**. (All per the spec STOP branch.)
- **Blocker reassigned:** RYA-542 asked "what is the right Ti I NLTE correction?" → **answered, +0.05
  ab-initio (Mallinson-2024, derived here)**. The remaining Ti blocker is now the **Ti I gf-scale**
  (RYA-521), surfaced here via ionization balance before it propagated into a wired abundance.
- **Grid is banked + derivation validated** — ready to wire the moment the Ti I gf-scale is fixed and
  the gate re-runs to ACCEPT.

## 5. Follow-on (opened)

1. **Ti I gf-scale fix + re-gate** — NIST-grade the Ti I line pool (retire the ungraded K10 synth-gf
   lines driving the 0.5-dex scatter), re-run the RYA-544 ionization-balance gate; on ACCEPT, wire the
   Mallinson grid (`Ti_Mallinson2024_PySME.csv`), register Ti forward/v2, flip the RYA-534 test.
2. **Model-atom vintage audit** (from the RYA-544 spec flag) — `atom.ti503b` = Bergemann-2011
   scaled-Drawin was the first-caught instance; audit the H-collision vintage across the Engine-A/B
   atom set (Mn/Cr/Co/Ni and any pre-Grumer-Barklem-2020 scaled-Drawin atoms), which systematically
   inflate over-ionization corrections.

## Sources
- Mallinson et al. 2024, A&A 687 A77 (grid: Zenodo 10753497) — solar Ti I +0.052, ab-initio GB2020.
- Mallinson et al. 2022, A&A 668 A103; Sitnova et al. 2020 (DETAIL) — corroborating ~+0.03.
- RYA-521 / RYA-395 / RYA-398 — the Ti I graded-pool +0.502 residual (now the ionization-balance blocker).
