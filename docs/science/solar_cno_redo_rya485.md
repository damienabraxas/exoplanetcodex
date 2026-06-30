# Solar CNO redo — O I 777 swap-test + regime legs + S (RYA-485/486)

Re-derive the measured solar O denominator across the three new references + verified RT,
record both RT legs for regime-matched differentials, and measure solar S — through the
RYA-485 iterate-til-right loop. **Validate-don't-tune: the goal is OUR measured value with
honest spreads; the Asplund numbers are comparisons, never fit targets.** All reference
atlases sourced from Sirius (`/mnt/codex-data/solar_reference/`, RYA-481/477).

## Headline — O I 777 four-reference swap-test (RYA-484 Lever-3) ✅ **continuum is NOT reference-driven**

Solar A(O) on O I 777, same method per reference (RYA-478 continuum-localize + Amarsi-2019
3D-NLTE), each FTS atlas converted wavenumber(cm⁻¹,vac) → air (Birch&Downs):

| reference | A(O)_LTE | χ²ᵣ | **A(O)_3D-NLTE** | A(O)_1D-NLTE |
|---|---|---|---|---|
| Kurucz-1984 KPNO (current) | 8.907 | 48.7 | **8.736** | 8.755 |
| Reiners-2016 IAG | 8.921 | 58.3 | **8.750** | 8.768 |
| Baker-2020 IAG-telluric | 8.915 | 54.7 | **8.744** | 8.762 |
| Wallace-2011 KPNO | 8.886 | 51.1 | **8.715** | 8.735 |

**Per-reference spread = 0.035 dex** (Kurucz 8.736 · IAG mean 8.747 · Wallace 8.715).

**The diagnostic resolves clean:** all four independent references agree on solar A(O) to
within ~0.035 dex. The O I 777 continuum is **not reference-driven** — so the RYA-483 Procyon
continuum lever (0.18 dex, UVES) is **intrinsic to the UVES arm**, not a KPNO reference
artifact. This *closes* the RYA-483 open question (one reference couldn't tell): the σ-inflation
on Procyon O stands as a real, instrument-specific systematic, not something a better solar
reference fixes. (IAG sits ~0.01 above Kurucz, Wallace ~0.02 below — all within the noise; no
reference is an outlier worth correcting.)

## Regime-matched solar O denominators (RYA-484 / the Issue-2 fix at its source)

The Sun is the 3D-NLTE leg (5772 K < the 6500 K STAGGER ceiling), but warm stars (Procyon,
6554 K) only have the 1D-NLTE leg — so their differential must use the **1D solar denominator**.
Recorded both:

- **solar_A(O)_3D-NLTE = 8.736** — the Sun's own best value (= the 371-banked 8.735, Δ +0.001).
- **solar_A(O)_1D-NLTE = 8.755** — the regime-matched denominator for >6500 K stars.

Using 8.755 for Procyon shifts [O/H] +0.085 → +0.066 (the ~0.02 dex RYA-485 Issue-2 fix). The
solar O denominator is re-confirmed and reference-robust at **8.736**, unchanged from banked.

## Sulphur — S I multiplet-8 6748.68 / 6757.15 (synthesis, HARPS)

MEASURED (the EW-pool 7.753 was curation-owed), Amarsi-2025 NLTE (small):

| line | A(S)_LTE | χ²ᵣ | NLTE δ | A(S)_NLTE |
|---|---|---|---|---|
| 6748.68 | 7.608 | 135.6 | −0.014 | 7.594 |
| 6757.15 | 7.471 | 137.5 | −0.017 | 7.454 |

**A(S)_NLTE = 7.524** (σ 0.099, n=2). Supersedes the EW-pool 7.753. **A finding, not tuned:**
it sits **+0.40 over Asplund 7.12** with 0.14 dex line-to-line scatter and high χ²ᵣ (~136) —
the multiplet-8 lines are weak/blended at solar strength. Reported as our measured value with
its honest spread; not massaged toward 7.12. The high-χ²ᵣ + line scatter is the curation owed
(stronger S lines / better continuum) — S stays a finding, off the EW-pool floor.

## New-vs-371 deltas
- **O: 8.736 (4-ref combined 3D-NLTE) vs 8.735 banked — Δ +0.001.** Denominator confirmed + now
  reference-triangulated and regime-split.
- **S: 7.524 (synthesis) vs 7.753 (EW-pool banked) — Δ −0.229**, and now NLTE-corrected with a real method.

## Deferred (named, not silently dropped)
- **Full C/N 3D-1D legs across the molecular arms (CH/CN/NH).** The Amarsi grid has the 3D/1D
  legs for atomic **C I** and **O I**; molecular C/N and N I post-hoc are a separate regime
  split — do alongside the C/N redo (the atomic O leg, the load-bearing one, is done here).
- **IR O I 844/926 nm (Lever-4).** In the Amarsi grid + in the IAG spvis/spnir range, but needs
  telluric-gated IR solar data reaching the line (the permanent IR rule); a follow-on channel.

## Validate-don't-tune
Corrections are cited-grid (Amarsi 3D-NLTE for O, Amarsi-2025 for S), never fitted; Asplund
8.69/7.12 are comparisons. Solar O lands at its measured 8.736 across four independent
references — earned, not snapped. S sits high — a finding, surfaced.

**Related:** RYA-485 (the 3-reference acquisition + RT verification), RYA-483 (the continuum
lever this closes), RYA-484 (the levers), RYA-486 (CNO-redo scope), RYA-371 (371 baseline),
RYA-489 (per-arm/per-reference convention), RYA-481 (Sirius references).
