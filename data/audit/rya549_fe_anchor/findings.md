# RYA-549 — Fe anchor: reconcile code≠register + ab-initio vintage assessment

_The RYA-524 headline finding, run to ground. Firewall: choose the Fe grid on provenance, never to hit a number; Fe II = untouched validator. Read at `origin/main` b69f7a6._

## Gate 0 — forensic (which grid is the real intent: code or register?)

**Verdict: the CODE (`Fe_Bergemann_MPIA`) reflects the last real decision; the register's "Amarsi-2022-NN" row was STALE doc.**

Evidence (git + code + provenance, not memory):
- `pipeline/nlte_corrections.py` docstring (lines 6–32) states plainly: **"LIVE Fe NLTE source (RYA-319): MPIA Bergemann MAFAGS-OS 1D NLTE grid … covers the full benchmark box (Procyon Teff, α Cen A/B, 55 Cnc). The Amarsi 3D-NLTE MLP below (Teff≤6500 ceiling) is ARCHIVED — kept only for the solar 3D-vs-1D cross-check (RYA-283)."**
- Git: `5d42702` **"integrate(rya-319): switch Fe NLTE to MPIA grid-lookup (single source); solar revalidated"** (+ `8b355bc` "data(rya-319): MPIA Fe I+II NLTE grid", `ff27604` "MPIA Fe scrape (gate passed)"), merged to main `472a237`. This is the last real Fe-grid decision.
- The Amarsi-2022 MLP files (`vendor/1L-3NErrors/fe1_model_*.p`, `fe2_model.p`) are present but ARCHIVED; coverage (docstring 22): **Teff 5000–6500, logg 4.0–4.5, [Fe/H] −3.0 to 0.0** → excludes Procyon (6554 K) and 55 Cnc ([Fe/H] +0.32). RYA-319 switched to MPIA precisely to get single-methodology coverage across the benchmark box; the ≲0.05 dex 3D→1D trade was documented and accepted (Amarsi 2022 Fig 7).
- The stale row: `CODEX_STATE_REGISTER.md:104` `| Fe I | Amarsi 2022 neural-network NLTE … | confirm bounds cover 55 Cnc … | RYA-251 Phase-3 / RYA-247 | SETTLED (selection); coverage per-target [confirm] |` — it even carried an unresolved `[confirm]` coverage flag (which RYA-319 answered by switching away). It predates/missed RYA-319.
- The wiring matrix `nlte_grid_availability.csv` was ALREADY correct (`Fe_Bergemann_MPIA.csv`, production). So the only code≠doc gap was this one register row.

**Integrity fix:** reconcile the register selection row → MPIA (with the RYA-319 rationale). Now **code = register = matrix**. No code change — the code was right.

## Gate 1 — is a covering ab-initio Fe grid banked? → NO

- The best-available ab-initio Fe NLTE grid is **Amarsi-2022** (NN on STAGGER 3D-NLTE, ab-initio Barklem H-collisions) — banked in-repo but coverage-ceilinged at Teff≤6500 / [Fe/H]≤0.0 → **cannot cover the metal-rich / hot benchmark targets.**
- No Fe `.grd` in the Sirius `amarsi_galah` set. No other ab-initio Fe grid covering [Fe/H]+0.32 is banked.
- → **Ab-initio migration is a sourcing blocker AND low-stakes (see vintage).** Deferred to RYA-550.

## Vintage — is the MPIA Fe grid the Ti/Mn/Cr inflated scaled-Drawin class? → NO (confirmed benign)

The MPIA Fe NLTE correction (read directly from `Fe_Bergemann_MPIA.csv`, median over 252 Fe I / ~75 Fe II lines per node):

| node | Fe I median δ | Fe II median δ |
|---|---|---|
| Sun (5800/4.5/0.0) | **+0.010** | +0.000 |
| 55 Cnc (5400/4.5/+0.20) | **+0.004** | −0.001 |
| 55 Cnc (5400/4.5/+0.35) | +0.004 | −0.001 |

- The Fe I correction is **+0.01** — an **order of magnitude smaller** than the Ti/Mn/Cr inflated class (+0.108). Fe is the dominant ion, ionization-balance-CALIBRATED (Bergemann's SH is fit, not pure Drawin=0), and empirically **ionization-balance-gated in production** (RYA-406: Fe I−Fe II −0.007…−0.015 PASS; RYA-407 scatter 0.138 honest floor). An inflated over-ionization correction would BREAK Fe I/II balance — it doesn't.
- The **+0.056** solar A(Fe) vs Asplund 7.46 is the documented **1D-NLTE-vs-3D-true scale offset** (RYA-336), NOT the NLTE correction (which is +0.01).
- **55 Cnc residual-risk** (RYA-546 Part-B logic on the anchor): the param-dependent, non-cancelling part of [Fe/H] = δ(55 Cnc) − δ(Sun) = +0.004 − 0.010 = **−0.006 dex — negligible.** Even a maximal ab-initio-vs-SH disagreement is bounded by the ~+0.01 total correction, i.e. sub-0.01 dex on [Fe/H].

**Conclusion:** the RYA-524 "Fe anchor on inflated scaled-Drawin" headline is **confirmed-but-benign** — the correction is tiny and balance-gated. The real defect was the integrity gap (now closed). Fe stays on MPIA.

## Validation (per spec)
- **Solar A(Fe I) vs 7.46:** 7.516 (+0.056), within tol; the offset is the RYA-336 1D-vs-3D scale, not a re-wire (grid unchanged).
- **Fe I/II NLTE balance < gate:** PASS — the production arbiter (RYA-406, Fe I−Fe II −0.007…−0.015); Fe is line-rich so NOT precision-limited (contrast Ti/Cr).
- **55 Cnc metal-rich δ:** characterized above (+0.004; non-cancelling part −0.006 dex).

## Disposition
- **Reconciled** (register + matrix → MPIA; code=register=matrix). **Vintage benign.** **Migration deferred → RYA-550.**
- **Unblocks RYA-527** — the Fe anchor is confirmed correct as-is; safe to re-freeze. (The re-freeze must not fold in any Fe change — there is none.)
- STOP conditions: none tripped (solar A(Fe) unchanged; grid not re-wired; the "intended" ab-initio grid is banked-but-coverage-limited, tracked as RYA-550).
