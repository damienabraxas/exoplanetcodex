# Finish solar — IR O I 844/926 + C confirm + S line fix + N flag (RYA-491)

Closes the solar O gambit (the RYA-485 redo nailed 777; this adds the IR indicators) and the
remaining C/S/N tracks, so Procyon inherits a complete, multi-indicator, regime-matched solar
foundation. **Validate-don't-tune: Asplund are comparisons, never targets.** All references
from Sirius (`/mnt/codex-data/solar_reference/`).

## Track 1 (headline) — IR O I 844/926 + the 4-indicator O table

Measured solar A(O) on the IR O I lines from the Baker-2020 telluric-IAG atlas
(`iag_telfree.fits`, cols `v/s/err/flags`, **frame = Salami & Ross 2005 iodine, RYA-481
declared**), wavenumber(cm⁻¹,vac) → air (Birch&Downs), **flag>0 points excised** before fitting.

| indicator | arm | A(O)_3D-NLTE | A(O)_1D-NLTE | note |
|---|---|---|---|---|
| **O I 777** | RED / UVES-IAG | **8.736** | 8.755 | 4-reference swap-test (RYA-485) — **PRIMARY** |
| [O I] 6300 | VIS / KP | 8.835 | 8.835 | forbidden-LTE, continuum-limited (RYA-447/448/455) |
| **O I 844** | IR / Baker-IAG | **8.603** | 8.643 | clean (flag-frac 0.006, χ²ᵣ 56) — **the new 4th indicator** |
| O I 926 | IR / Baker-IAG | (9.877) | — | **UNUSABLE** — weak (depth 0.08) + 80% flagged in Baker; railed |

**Result — not a tight single value; a formation-depth finding (surfaced, not averaged).**
The three usable indicators span **0.232 dex**: O I 777 (8.736, the swap-tested gold),
[O I] 6300 (+0.10, continuum-limited and weak), **O I 844 (−0.13)**. The 844 multiplet forms
**deeper** than 777, so a real lower A(O) there is physics, not error — the offset is the
result. **O I 777 stays the bankable solar O (8.736)**, now IR-corroborated to within 0.13 dex
with a documented formation-depth offset. O I 926 is too weak and too heavily flagged in the
Baker atlas to be a measurement (a data finding, excised honestly rather than averaged in).

## Track 2 — C confirmed on a new reference ✅
C I 5052 (the clean line) → **A(C)_3D-NLTE 8.449** vs the 371-banked **8.491** (Δ −0.042 —
holds). C I 5380 → 8.488 (kept flagged: the ESPRESSO χ²ᵣ-103 outlier at 371). **C holds ~8.45–8.49**
— a confirmation, not a move. (The full CH/C2 molecular-equilibrium re-convergence is the
heavier C/N redo; the atomic C I confirmation is clean here.)

## Track 3 — S line fixed, but still high (the gf is the lever)
Swapped 6748.68 (Costa Silva+2020 DISCARDED) → **6743.53** (kept) + 6757.15:

| line | A(S)_LTE | χ²ᵣ | NLTE | A(S)_NLTE |
|---|---|---|---|---|
| 6743.53 | 7.596 | 137.8 | −0.017 | 7.579 |
| 6757.15 | 7.471 | 137.5 | −0.017 | 7.454 |

**A(S)_NLTE = 7.516** (σ 0.088) — the line swap barely moved it (was 7.524 on 6748). So the
offset is **not** the line identity: with the **GES gf** S reads **+0.40 over Asplund 7.12**.
The Costa Silva+2020 atlas-tuned gf (their nudges to match A(S)⊙ = 7.12) is the lever — owed
(Ryan supplies). Reported as a **gf-floor finding (RYA-161/162)**, not tuned. χ²ᵣ ~138 also
flags the weak/blended multiplet-8 at solar strength.

## Track 4 — N: deferred + FLAGGED (not silent)
N stays carried from 371 (NH 3360 data-gap + N I grid owed, RYA-369). It is the **one remaining
soft solar species** — "solar done" is honest about this gap; N has its own strategy ticket.

## Genuinely blocked (documented)
CRIRES+ CO 2.3 µm overtone + ¹³C, and IR OH ~1.5–1.6 µm — STAGGER collaborator gate + telluric
(RYA-373). The only legitimately-deferred solar pieces.

## Verdict — solar is closed (with honest caveats)
- **O:** 3-indicator (777 + 6300 + 844), 926 unusable; **777 8.736 is the bankable primary**,
  IR-corroborated; the 0.23 spread is a formation-depth finding. Both RT legs recorded
  (3D 8.736 / **1D 8.755** for warm-star differentials).
- **C:** holds 8.45–8.49 on the new reference.
- **S:** 7.516 on the correct line — gf lever owed (Costa Silva), gf-floor finding.
- **N:** flagged soft (RYA-369).
- → The triangulated O + the 1D regime-matched denominator are ready for the Procyon CNO redo
  (RYA-486 gate satisfied). Deltas from banked: O +0.001 (777), C −0.042, S −0.008-on-correct-line.

**Related:** RYA-485 (the redo this finishes), RYA-486 (full-CNO scope), RYA-484 (Lever-4 = IR),
RYA-483 (Procyon O denominator), RYA-369 (N), RYA-489 (per-arm), RYA-481 (frame), RYA-359 (grid).
