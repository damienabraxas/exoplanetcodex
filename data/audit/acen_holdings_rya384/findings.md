# RYA-384 — alpha Cen A & B holdings audit (VALD + spectra)

Read-only readiness audit for the first science targets (Solar -> Procyon -> **alpha Cen A/B** -> 55 Cnc); alpha Cen A carries the Aug-2026 JWST deadline (RYA-116). A and B audited separately. Do NOT merge.

## Part A — VALD line-list coverage (per star)

| | alpha Cen A | alpha Cen B |
|---|---|---|
| assembled span | 1150.3–16997.2 Å | 1150.3–16997.2 Å |
| total lines | 141,940 | 155,900 |
| bands present | UV + optical + Y/J/H | UV + optical + Y/J/H |
| **K-band (>1.70 µm)** | **ABSENT (0 lines)** | **ABSENT (0 lines)** |
| HFS state | **ON** | **ON** |

- **HFS-on confirmed by data** (not just filename): the Mn I 6013/6016/6021 triplet carries ~14 hyperfine components per feature in the optical raw, and Co I in the NIR raw shows ~6 components per feature — for **both** stars. No HFS-off leak (unlike the 55 Cnc NIR quarantine). The 376 audit reports HFS="unknown" only because the alpha Cen filenames lack the `hfson`/`hfsoff` token — the data are HFS-on.
- **Red-optical diagnostics present** (HFS-on): O I 7771-5 (A:3 / B:2 lines), N I 8216 (1/1), N I 8680 (2/3), K I 7699 (4/4, HFS-split).
- **Near-UV molecular bands present**: NH ~3360 (A:229 / B:272), CN ~3883 (A:130 / B:173).
- **Element x band**: the 10 registered-NLTE + Fe/C/N/O diagnostics are ✓ across optical/Y/J/H for both stars. Minor per-band sparsity (V in J/H for A; P/Zn/Cu/Sr thin in some bands) = VALD line density, not a defect.

**K-band owed: YES** — a K-band (>1.70 µm) HFS-on extraction is owed for the CO arm (same as solar/55 Cnc), for both stars.

## Part B — spectral data inventory + binary contamination

Header-truth target split (HD128620 = alpha Cen A; HD128621 = alpha Cen B). See spectral_inventory.csv.

**HEADLINE — the alpha Cen folders are NOT cleanly separated by star.** UVES and the entire HST/STIS UV set are DUPLICATED into both the "Alpha Centauri A" and "Alpha Centauri B" folders, and each is internally a MIX of A and B by header target. **Selection MUST be by header target (OBJECT/TARGNAME), never by folder name** — feeding "Alpha Centauri A/UVES" wholesale would put ~75% alpha Cen B into an "A" reduction.

| instrument | files | header truth | usable as |
|---|---|---|---|
| **HARPS** (A folder) | 88 | **all HD128620** | **clean alpha Cen A** — primary optical 378–691 nm, SNR~320, BERV in header |
| GIRAFFE (A) | 22 | 9 A + 13 coord-named | mostly A, narrow 533–562 nm, low SNR — minor |
| FEROS (A) | 2 | Gl-559 (system, ambiguous) | optical+red 352–922 nm, SNR~173 — A/B disambiguation owed |
| UVES (A folder) | 75 | **19 A / 56 B** | MIXED — duplicate of B/ESO; mislabeled folder |
| UVES "B ESO" (B folder) | 79 | 19 A / 56 B (+4) | MIXED — same bundle + 4 extra |
| STIS uv-e140m (both) | 42 ea | 25 A / 17 B | MIXED FUV 1140–1730 Å |
| STIS uv-e140h (both) | 22 ea | 3 A / 19 B | MIXED FUV (high-res) |
| STIS nuv-e230h (both) | 43 ea | 20 A / 23 B | MIXED NUV 2568–2846 Å |
| STIS nuv-e230m (both) | 4 ea | 4 B | clean B NUV |
| **CHIRON** (B folder) | 76 | CTIO format (target in ext) | alpha Cen B primary (2021 CTIO), optical — format-aware read owed to confirm |
| Phoenix (A) | 10 | raw Gemini, sample all "Dark" | **unusable** — raw unreduced IR frames, no wavelength solution |

- **HST STIS** adds the FUV/NUV arm (1140–1730 + 2568–2846 Å) for UV C/N/O — but mixed A/B and only in echelle windows (not continuous UV).
- **GDAS (RYA-380)**: FEROS red (to 922 nm) and CHIRON (red-optical) need per-night GDAS telluric handling.
- **No reduced NIR (Y/J/H) spectra exist** — Phoenix is raw darks. The VALD NIR line list is ahead of the spectra.

## Part C — lines-vs-spectra alignment (per star)

| band | VALD lines | spectra | aligned? |
|---|---|---|---|
| FUV/NUV (1150–3780 Å) | ✓ HFS-on | HST/STIS windows only (1140–1730 + 2568–2846) | PARTIAL (echelle windows, mixed A/B) |
| optical (3780–6910) | ✓ | HARPS (A), CHIRON (B), FEROS/UVES | **aligned** |
| red-optical (6910–9500) | ✓ (O/N/K present) | FEROS→922 nm (A), CHIRON (B); HARPS stops at 691 | aligned via FEROS/CHIRON (not HARPS); GDAS owed |
| Y/J/H (9500–17000) | ✓ (~1900 lines) | **none reduced** (Phoenix raw) | **MISMATCH — lines without spectra** |
| K (>17000) | absent | none | double gap |

## What's owed before the alpha Cen pipeline run

1. **De-contaminate by header target (CRITICAL).** Re-select UVES + all HST/STIS by OBJECT/TARGNAME (HD128620=A, HD128621=B), not folder. The folders are duplicated and internally mixed; HARPS is the only cleanly-separated optical set (A).
2. **Confirm CHIRON = alpha Cen B** via a CTIO-format-aware read (target is in the extension), then it is B's optical primary.
3. **K-band HFS-on VALD extraction** (>1.70 µm) for the CO arm — both stars (owed, same as solar/55 Cnc).
4. **Reduced NIR spectra** (Y/J/H) if NIR diagnostics are needed — Phoenix is raw/unusable; another IR source or a reduction is owed.
5. **Per-night GDAS** (RYA-380) for FEROS-red (A) and CHIRON (B).
6. **FEROS A/B disambiguation** (Gl-559 is the system name).

**Bottom line:** alpha Cen A is in good shape for an OPTICAL run on the clean HARPS set (+ FEROS red with GDAS); the blocker is the duplicated/mixed UVES+STIS folders (select by header) and the absent reduced NIR. alpha Cen B's optical primary (CHIRON) needs a format-aware confirm. VALD coverage is strong (1150–16997 Å, HFS-on) with the K-band extraction the one line-list gap.
