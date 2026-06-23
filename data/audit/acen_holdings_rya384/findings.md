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

---

## ADDENDUM (2026-06-23) — systemic binary mislabeling + spectral-type re-sort

After deeper vetting (and the user dropping more ESO downloads), the headline hardened: **for this close visual binary the ESO archive target labels are systematically unreliable across EVERY instrument.** The only trustworthy discriminator is the spectrum itself — G2V (alpha Cen A) shallow lines vs K1V (alpha Cen B) ~2x deeper lines, validated by RV-lag cross-correlation against templates built from the HARPS set.

### Confirmed by spectral type (high confidence)
- **HARPS "A folder" (88) = 75 A + 13 B.** The 13 highest-SNR (SN0 ~300, single night 2010-04-02) are spectroscopically K1V = alpha Cen B, mislabeled HD128620. (CC margins: known-A CC_A=1.0/CC_B=0.52; known-B CC_A=0.5/CC_B=0.96.)
- **HARPS archive "HD128621" query** returns mostly A — selecting by SNR biases to the brighter star (A). Harvested 12 genuine B from 150 downloaded.
- **ESPRESSO (15) labeled "alf Cen B" = alpha Cen A** (depth 0.18, firmly G-dwarf; resolution-matched). Excellent high-res A optical (3772-7900 A), just mislabeled.

### Could NOT be auto-classified (cross-instrument CC too weak; routed by label, FLAGGED UNVERIFIED)
- **UVES** (CC ~0.22/0.00 vs both templates — unreliable), **CHIRON** (CTIO IDP, no usable window -> _REVIEW), **FEROS**, **GIRAFFE**, **HST/STIS** (UV), **NIRPS** (YJH).
These need dedicated per-instrument vetting (resolution-matched templates or RV cross-match to simultaneous optical).

### IR data (the CO / NIR gaps — both stars, mixed targets)
- **CRIRES** in both folders: Y/J/H/K settings. **All K-band settings (K2148/66/92/2217, 1.95-2.49 um) cover the 12CO (2.2935 um) AND 13CO (2.3448 um) bandheads** = the CO arm. Targets mixed; vetted tentatively by CO-bandhead depth (strong->B). The "Star S5" K-band CO is at alpha Cen coords.
- **NIRPS**: 322 processed YJH (0.97-1.92 um) spectra exist; 8 best downloaded (SNR up to 403). Labels mixed; A/B in the NIR needs a dedicated discriminator.
- **GEMINI = TReCS mid-IR (10-20 um)** -> NOT-FOR-ABUNDANCES; **Phoenix** = raw Gemini darks -> NOT-FOR-ABUNDANCES.

### Vetted folder tree (NON-DESTRUCTIVE copy; originals preserved)
`Alpha Centauri (vetted)/Alpha Cen A` (302) + `Alpha Cen B` (147) + `_REVIEW` (102: CHIRON/CRIRES-YJH/FEROS/GIRAFFE/UVES) + `_NEEDS-REVIEW` (HARPS-borderline, NIRPS) + `_NOT-FOR-ABUNDANCES` (TReCS, Phoenix). Built by scripts/reorg_acen_by_spectral_type_rya384.py; per-file routing + method/CC in reorg_by_spectraltype_manifest.csv. **alpha Cen B real optical = 25 HARPS (confirmed) + UVES + CHIRON; high-res B optical still owed (the ESPRESSO turned out to be A).**

### Standing recommendation
Going forward, do NOT trust HD128620/HD128621/"alf Cen B" labels. Pull a generous batch per instrument and sort by spectral type. Confirmed-reliable method exists for HARPS-resolution optical; UVES/IR/UV need per-instrument templates.

---

## PROOF: alpha Cen A/B labels swapped (NIRPS) — two independent discriminators (2026-06-23)

Demanded a rigorous proof. Built one with TWO PHYSICALLY INDEPENDENT discriminators that must
agree, anchored to the validated HARPS ground truth. Figure: plots/acen_label_swap_proof.png
(4 panels: HARPS A/B overlay; HARPS bimodal depth histogram = the 13 mislabeled B; NIRPS
flux-vs-depth 2D; NIRPS spectrum overlay). Script: scripts/prove_acen_label_swap_rya384.py.

Quantitative (decisive):
- HARPS GROUND TRUTH: confirmed-A line-depth = 0.224 (n=75); confirmed-B = 0.506 (n=13).
- NIRPS 'AlphaCenB'-labeled: J-depth 0.132 (SHALLOW=G2) + flux 3.9e-10 erg/cm2/s/A (BRIGHT)
  -> alpha Cen A on BOTH axes.
- NIRPS 'alf Cen A'/'Star S5'-labeled: J-depth 0.797 (DEEP=K1) + flux 4.1e-12 (FAINT, ~95x
  dimmer) -> alpha Cen B on BOTH axes.
The spectroscopic axis (line depth, 6x ratio) and the photometric axis (flux, 95x ratio) are
INDEPENDENT and agree perfectly: the 'AlphaCenB' label = actual alpha Cen A. LABELS SWAPPED.
(NIRPS flux is calibrated: BUNIT = erg.cm-2.s-1.angstrom-1.)

Acted on it: re-queried ESO for the swapped-label B candidates (label 'alf Cen A' + 'Star S5')
and verified each by J-depth -> 22 confirmed-B NIRPS routed to Alpha Cen B/NIRPS (up from 2);
15 -> Alpha Cen A/NIRPS. alpha Cen B NIR (YJH) is now well-populated.

---

## IR target ID on all accounts + NIRPS-A sufficiency (2026-06-23)

NIRPS (proof-grade, see plot): labels SWAPPED, confirmed by 2 independent axes. Acted ->
20 confirmed-A + 20 confirmed-B NIRPS (YJH 0.97-1.92um) routed to Alpha Cen A|B/NIRPS by
J-band line-depth (clean split). alpha Cen A NIRPS = MORE THAN ENOUGH (293 available in the
ESO archive; the 20 grabbed include the top-SNR set 232-403; co-add -> SNR>1000). alpha Cen B
NIRPS now well-populated (20, archive SNR up to 179).

CRIRES (Y/J/H/K) target ID — CANNOT be proof-graded by the quick two-axis method, and here is
why (honest):
- FLUX axis dead: CRIRES BUNIT = ADU (counts, NOT flux-calibrated) -> exposure-dependent, so
  the A-label vs B-label "fluxes" are ~equal and carry no brightness information.
- DEPTH axis weak: the reduced CRIRES line depths are uniformly ~0.05 across Y/J/H/K for the
  'alf Cen A'/'alf Cen B' pairs (over-flattened continuum / telluric-dominated) -> no clean
  G2-vs-K1 separation.
- ONLY solid IR signal: 'Star S5' K-band has deep CO (~0.18 vs ~0.05) = a real K-dwarf =
  alpha Cen B's CO arm (confident).
So the CRIRES 'alf Cen A/B' Y/J/H/K target labels are UNVERIFIED (flagged, not pretended);
resolving them needs telluric correction + NIR (G2 vs K1) templates — a dedicated reduction
task, not quick vetting. The CO arm (Star S5 K-band) for alpha Cen B is confirmed.
