# CNO method step-back — 3 corrections before RYA-362 (RYA-485)

Diagnose-only pass (no production fixes this round). Three methodology issues, sequenced
Issue 2 → 1 → 3.

## Issue 2 — VERIFY solar O is 3D/synthesis, not silent-LTE ✅ **CORRECT (no silent-LTE)**

Traced the live `pipeline.nlte_cno` path (`scripts/verify_solar_o_modeling_rya485.py`),
not assumptions. The solar O is the differential denominator, so this was the
correctness-critical check.

| check | result |
|---|---|
| **O I 777 leg at solar Teff** | Sun (5772 K) → **3D-NLTE leg** (table3), δ = **−0.171** ✓. Procyon (6554 K) → 1D-NLTE leg (table6), δ −0.540; 3D returns NaN (above the 6500 K STAGGER ceiling). The Sun is **not** silently 1D. |
| **[O I] 6300 method** | `forbidden_blend` synthesis ([O I] 6300.30 + Ni I 6300.34 joint, Ni gf = Johansson 2003), flagged `lte_forbidden_insensitive`. The grid **has** a 630.0 nm node — its 3D differential at solar params is **+0.001 dex** (1D-leg exactly 0). So the forbidden-LTE treatment is the **correct** method, not a silent skip of a real 3D term (RYA-447 confirmed quantitatively). |
| **per-arm solar O spread** | O I 777 (3D-NLTE) ≈ 8.78 vs [O I] 6300 (forbidden-LTE) 8.835 → spread **≈ 0.05 dex**. |

**The spread does NOT mirror Procyon's 0.18.** Ryan's hypothesis was that the solar O
spread mirroring Procyon's would signal an upstream solar bug. It doesn't: solar O is
internally consistent to ~0.05, while Procyon's 0.18 was the **O I 777 continuum-localization
lever on the UVES arm** (solar 777 continuum shift was only −0.049; Procyon's −0.231). So
**Procyon's continuum sensitivity is UVES-instrument-specific, not inherited from the solar
reference.** Confirmed, not assumed.

### The real finding (named, NOT fixed this pass): differential regime-mismatch
The Procyon [O/H] differences **Procyon-1D-NLTE** (8.82, forced to the 1D leg because 6554 K
> the 6500 K 3D ceiling) against **Sun-3D-NLTE** (8.736). The two sides of the differential
are in **different RT regimes**. The solar 3D-minus-1D term at 777 is **+0.019 dex**, so a
regime-matched differential (Procyon-1D vs **Sun-1D-NLTE** = 8.755) gives **[O/H] +0.065**
instead of **+0.085** — a −0.019 dex shift.

**Fix to apply later (RYA-484/362):** form the Procyon (and any >6500 K star's) O I 777
differential against the **Sun on the same leg** (1D-NLTE), not the Sun's own best-regime
(3D-NLTE) value. Each star's *own* abundance still uses its best leg; only the *differential
denominator* must be regime-matched. Small (~0.02 dex) but it's a real, citable systematic
that otherwise rides silently in every >6500 K O differential.

## Issue 1 — per-arm ARM-MATCHED literature comparison (re-frame)

The error: comparing our **combined multi-arm carbon spread (0.253 dex)** against Bruntt+2010's
**optical-C I-only −0.04**. That violates RYA-307/237 ("one product per star, regime as
provenance layer, per-region combined by inverse-variance"). The rule, restated + enforced:

- Report abundances **per-arm with per-arm spread** (VIS / UV / IR), THEN combine by
  inverse-variance into one product.
- **Cross-compare to literature ARM-MATCHED only.** Bruntt measured Procyon C from **optical
  C I lines** — no CH, C2, UV, or IR. So the only valid Bruntt comparand is **our optical
  C I (5052 / 5380) vs their optical C I.** The molecular CH/C2 lines have **no Bruntt
  comparand** and must not be folded into a "vs Bruntt" spread.
- **Our combined multi-arm CNO spread has no external comparand** (no one else does multi-arm
  CNO) — it is our **novel product**, reported as such, never benchmarked against single-arm
  studies (which would make careful multi-arm work look worse than work that never touched the
  hard lines).
- **Action (for the RYA-484 / Procyon C/O writeup):** split the finding into (a) internal VIS
  C I spread — the Bruntt-comparable number — and (b) the molecular-vs-atomic spread — our
  finding, no comparand. Same for O (O I 777 is the IR-ish primary; [O I] 6300 the optical
  cross-check; no single study spans both).

## Issue 3 — non-Kitt-Peak solar 777 reference (acquire) ⏳ **fetching to Sirius**

**Inventory (on disk now):** the registered `SOLAR_REFERENCE_SPECTRA` are KPNO Kurucz-1984
flux atlas (296–1300 nm — **the only 777 reference**, covers 777 + 6300), CALSPEC UV composite
(low-res), ACE-FTS IR (2289–2349 nm), IRTF (deferred). **HARPS/Vesta cannot help — O I 777 is
outside HARPS range (≤6910 Å).** So today the 777 continuum lever (RYA-483) can't be isolated:
one reference, can't tell if the slope is KPNO-specific. **IAG was absent** (Sirius `sol/` empty).

**Acquired now (to Sirius — the RYA-477/481 pattern, `scripts/sirius_fetch_solar_atlas_rya485.py`,
no silent partials, `PROVENANCE.json` per record).** Three independent references for a
**2-vs-1 triangulation** of the 777 continuum lever, not just one:

| # | reference | source | integrity | dest (Sirius) |
|---|---|---|---|---|
| 1 | **Baker+2020 telluric-corrected IAG** (0.5–1.0 µm, the cleanest 777) | Zenodo 10.5281/zenodo.3598136 (13 files, 551 MB) | **upstream md5 ✓ (13/13)** | `solar_reference/iag_baker2020/` |
| 2 | **Reiners 2016 IAG base** (VIS 4050–10650 Å, R≈10⁶) | CDS J/A+A/587/A65 (spvis/spnir .gz, ~115 MB) | gunzip-integrity + computed md5 (CDS publishes none) | `solar_reference/iag_reiners2016/` |
| 3 | **Wallace 2011 telluric-removed KPNO** (sptr.reg1–6) | NSO `nispdata.nso.edu/.../Wallace_2011_solar_flux_atlas/` (~80 MB) | computed md5 (NSO publishes none) | `solar_reference/wallace2011_kpno/` |

**The triangulation logic (RYA-484 Lever-3):** our current 777 reference is Kurucz-1984 KPNO.
Swap it for each of the three and re-fit solar (and Procyon) 777:
- IAG (1+2) disagrees with our-KPNO → Wallace-2011 (3) breaks the tie: Wallace ≈ Kurucz ⇒
  **KPNO-family** systematic; Wallace ≈ IAG ⇒ specific to the **Kurucz reduction** we happen to use.
- All three agree ⇒ the 777 continuum is **not** reference-driven; the σ-inflation is intrinsic
  (the RYA-483 finding stands, σ stays).

The stitched **`iag_telfree_solaratlas.fits`** (Baker) is the primary 777+6300 swap-test file.
**RYA-481 discipline:** each new reference is *declared* explicitly before any swap-test, never
silently substituted for Kitt Peak. Note (3) was **NOT** what the A&A HTML's `ftp://vso.nso.edu`
implied for IAG — that path is the NSO KPNO atlases; IAG is CDS. NSO FTP (port 21) is firewalled,
so all three fetch over **https** (CDS + `nispdata.nso.edu`), reachable from Sirius.

## Gate to RYA-362 / RYA-484
- ✅ solar O modeling verified correct (no silent-LTE);
- ⚠️ one named fix owed (differential regime-match, ~0.02 dex) — apply in RYA-484/362;
- ✅ comparison re-framed (arm-matched);
- ⏳ 2nd solar 777 reference (IAG) landing on Sirius for the Lever-3 swap-test.

**Related:** RYA-307/237 (per-arm presentation), RYA-484 (O-resolution levers), RYA-483
(provisional O bank), RYA-447/448 ([O I] 6300 path), RYA-362 (NLTE cross-check), RYA-459
(solar atlases), RYA-481 (no silent substitution), RYA-477 (Sirius staging).
