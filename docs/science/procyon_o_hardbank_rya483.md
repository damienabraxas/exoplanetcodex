# Procyon oxygen — hard-bank attempt (RYA-483)

**Verdict: NOT hard-bankable.** Procyon O I 777 [O/H] +0.085 is **banked provisionally with a
named residual**, not hard-banked. Two findings block a tight bank, both larger than RYA-478
expected. Primary indicator stays **O I 777 (1D-NLTE)**; differential is vs **our measured Sun
(8.735)**, never Asplund. The method does not move the Sun (the RYA-478 guard holds).

## The number (provisional)

| quantity | value |
|---|---|
| Procyon A(O) 1D-NLTE (O I 777, UVES) | **8.82** |
| [O/H] vs our Sun (8.735) | **+0.085** |
| 1σ (continuum zero-point dominated) | **± 0.186** |
| O I 777 fit χ²ᵣ | 46.7 (still ≫ 1) |
| solar control (KPNO O I 777, 1D-NLTE) | 8.736 = our Sun + 0.001 ✓ unmoved |

## Part A — cross-instrument zero-point: **~0.18 dex, not 0.04–0.08**

The differential is A(O, Procyon, **UVES**, 777) − A(O, Sun, **KPNO**, 777) — two instruments.
Characterized via the continuum-localization lever measured **per arm** (the same RYA-478
local-continuum method, applied global-vs-local on each):

| arm | A(O) global | A(O) local | shift |
|---|---|---|---|
| Sun / KPNO | 8.956 | 8.907 | **−0.049** |
| Procyon / UVES | 9.591 | 9.360 | **−0.231** |

- **Differential continuum zero-point = −0.231 − (−0.049) = −0.182 dex.** The UVES O I 777
  triplet sits on a sloped pseudo-continuum that the KPNO solar flux atlas does not — so the
  local-vs-global continuum choice swings Procyon A(O) by **0.23 dex**, vs 0.05 on the Sun.
- **Crucially, χ²ᵣ barely moves (46.4 global vs 46.7 local)** — the data does *not* discriminate
  which continuum is right, yet A(O) swings 0.23 dex. This is a real, under-constrained systematic.
- Reference accuracy: KPNO O I 777 reproduces our Sun to **+0.001** → the reference is on-scale;
  the solar KPNO-vs-ESPRESSO spread (~0.04 dex, RYA-460) bounds the reference-instrument term.
- **Disposition: fold into σ (≈0.186, the continuum-differential in quadrature with the 0.04
  reference spread).** It is *not* a stable additive offset — there is no shared star+diagnostic
  (no Procyon-KPNO, no Sun-UVES at 777) to correct one. Reported, never silently absorbed.

> This contradicts RYA-478's 0.04–0.08 estimate: the cross-instrument zero-point is **~0.18 dex**,
> and it is the dominant uncertainty on Procyon O.

## Part B — independent [O I] 6300: **terminal, leg unavailable**

Fit Procyon HARPS [O I] 6300 independently (continuum-localized, Ni I 6300.34 pinned at the
Johansson-2003 gf via gf_resolver, RYA-365) to de-anchor the RYA-348 UVES-vs-HARPS leg (which
was Caffau-anchor-dominated).

- **Solar [O I] 6300 control railed to 7.256** (vs our Sun 8.735, |Δ| 1.48) — a failed/terminal
  fit, exactly the RYA-447/448/455 lesson: [O I] 6300 is continuum + Ni-blend terminal, not
  science-grade. Procyon came out 8.228, but with the solar control broken the line provides **no
  usable de-anchored number**.
- **The corroboration leg does not materialize.** [O I] 6300 cannot confirm or refute O I 777.
  This is a finding, not a bank — and it means O I 777 has **no independent in-house cross-check**
  on Procyon today.

## Bank decision

**Provisional, not hard-banked.** A(O) 8.82, [O/H] +0.085 ± 0.186 (1σ continuum-zero-point
dominated), primary = O I 777 1D-NLTE, vs our Sun 8.735. This is the number α Cen A / 55 Cnc O
get compared against — carried **with its σ and the two blockers**, not as a tight value.

**Blockers to a tight (hard) bank, named:**
1. **The UVES O I 777 continuum-localization lever (0.23 dex, no χ² discrimination)** — needs an
   independent continuum determination or a **second O I 777 epoch / instrument** for Procyon to
   break the degeneracy. χ²ᵣ 46.7 also says the 777 fit is not yet clean.
2. **[O I] 6300 is terminal** (RYA-447/448/455 confirmed on Procyon too) — it cannot serve as the
   corroboration leg; an independent O cross-check must come from elsewhere (UV C/O, a CO arm).

Solar control unmoved (8.736); validate-don't-tune held throughout — nothing massaged toward 8.69
or 8.735. STOP at verdict; not merged.

**Related:** RYA-478 (provisional measurement), RYA-348 (Phase 2 parent), RYA-447/448/455 ([O I]
6300 terminal-precision), RYA-365 (Ni gf), RYA-359 (Amarsi NLTE), RYA-460 (KPNO-ESPRESSO zero-point).
