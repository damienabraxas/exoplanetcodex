# RYA-684 — isotope-fraction double-application audit

**Verdict: the mechanism is real, confirmed at engine-source level, and present on five
species of the shipped VALD "for-grid" lists — and it reaches NO live abundance.**
Every reported value is measured on a surface that is fraction-free and therefore correct.
The +0.300 dex RYA-565 saw is a property of a control leg it deliberately ran and
deliberately did not adopt.

Regenerate with `scripts/rya684_isotope_gf_audit.py` on Sirius; the machine-readable
record is `data/audit/rya684_isotope_gf_audit.json`.

---

## 1. The mechanism, end to end

Four links, each read off the engine or the shipped file rather than inferred:

1. **The engine applies the fraction.** `bsyn.f:1350-1351` multiplies the absorber
   population of every isotope-coded species by `isotopfrac(atom, isotope)`:

   ```fortran
   ntot(k) = ntot(k) * isotopfrac(atom(i), isotope(i))
   ntt(k)  = ntt(k)  * isotopfrac(atom(i), isotope(i))
   ```

2. **The fractions are the engine's own defaults.** `makeabund.f` hardcodes
   `isotopfrac(63,151) = 0.478` and `isotopfrac(63,153) = 0.522`, and every bsyn harness
   in `scripts/` emits `'ISOTOPES : ' '0'` — zero overrides — so those defaults are what
   runs. (The iSpec path differs: `ispec/synth/turbospectrum.py:260` writes an explicit
   `ISOTOPES :` block from `input/isotopes/SPECTRUM.lst`.)

3. **`makeabund.f` states the convention that goes with it**, in its own words:

   > `isotopfrac(x,0)==1.0`; this is used in the case of no isotopes wanted in
   > calculation. It will ensure that the total abundance of the species will be
   > used.(if for example isotopic factors were included in gf-values)

   So a list either isotope-codes and ships fraction-free gf, or writes `Z.000` and folds
   the fraction into gf. Doing both applies it twice.

4. **The shipped VALD list does both.** Ba II 5853 in
   `vald-5800-6300-for-grid.list` is written as seven isotope-coded blocks whose log gf
   values reconstruct to one physical number once the fraction is divided back out:

   | isotope | shipped log gf | `isotopfrac` | shipped − log₁₀f |
   |---|---|---|---|
   | 138 | −1.144 | 0.7170 | −0.9995 |
   | 137 | −1.950 | 0.1123 | −1.0004 |
   | 135 | −2.181 | 0.06592 | −1.0003 |

   All seven land on **−1.000** — which is exactly the "total loggf −1.0" the RYA-559
   Ba verdict text quotes. The fraction is unmistakably already in the shipped gf, and
   `bsyn` then applies it again.

The audit does not take that on one species' word: it runs the reconstruction over every
isotope-coded feature on every surface and reports which hypothesis collapses.

### The size of the error

For **co-located** structure (atomic isotope/HFS components inside one blended feature)
the feature should carry `Σᵢ fᵢ·gf = gf` and instead carries `Σᵢ fᵢ²·gf`, so the fitted
abundance is high by

    ΔA = −log₁₀( Σᵢ fᵢ² )

For Eu II that is **+0.30020**, against RYA-565's measured **+0.300**. Note this is *not*
log₁₀2 = 0.30103 — Eu's isotopes are 0.478/0.522, not 0.500/0.500. The two agree to three
decimals by coincidence of Eu's near-even split; Ba's seven isotopes give +0.26938, which
log₁₀2 would have misstated by 0.03 dex.

For **separated** structure (molecular isotopologues — ¹²CH and ¹³CH are different lines at
different wavelengths) each line should carry `fᵢ·gfᵢ` and instead carries `fᵢ²·gfᵢ`, so
the offset is `−log₁₀ fᵢ` **per isotopologue**. Using the co-located formula there would
have called ¹³CH a 0.01 dex effect when it is a 1.96 dex one.

---

## 2. Per-species exposure table

`exposed = YES` means the shipped log gf carries the isotope fraction **and** the block is
isotope-coded, so `bsyn` applies it a second time.

| species | surface | isotopes coded | shipped log gf | exposed | dex if folded |
|---|---|---|---|---|---|
| Li I | ges(iSpec path) | 6,7 | fraction-free (full gf) | no | — |
| Li I | ts-nlte-ges(Engine-B) | 6,7 | fraction-free (full gf) | no | — |
| **Li I** | **vald(bsyn harnesses)** | 6,7 | **fraction folded in** | **YES** | +0.0339 (⁷Li) … +1.1249 (⁶Li) |
| **Ca II** | **vald(bsyn harnesses)** | 40,42,43,44,46,48 | **fraction folded in** | **YES** | +0.0135 (⁴⁰Ca) … +4.3979 (⁴⁶Ca) |
| Cu I | ges(iSpec path) | 63,65 | fraction-free (full gf) | no | — |
| Cu I | ts-nlte-ges(Engine-B) | 63,65 | fraction-free (full gf) | no | — |
| **Cu I** | **vald(bsyn harnesses)** | 63,65 | **fraction folded in** | **YES** | +0.2415 (blended feature) |
| Ba II | ges(iSpec path) | 134–138 | fraction-free (full gf) | no | — |
| Ba II | ts-nlte-ges(Engine-B) | 134–138 | fraction-free (full gf) | no | — |
| **Ba II** | **vald(bsyn harnesses)** | 130–138 | **fraction folded in** | **YES** | +0.2694 (blended feature) |
| La II | ges / Engine-B | 139 | undecidable (single isotope) | ≤0.0004 | bounded, nil |
| Pr II | ges / Engine-B | 141 | undecidable (single isotope) | 0.0000 | mono-isotopic |
| Nd II | ges / Engine-B | 142–150 | fraction-free (full gf) | no | — |
| Sm II | ges / Engine-B | 144–154 | fraction-free (full gf) | no | — |
| Eu II | ges(iSpec path) | 151,153 | fraction-free (full gf) | no | — |
| Eu II | ts-nlte-ges(Engine-B) | 151,153 | fraction-free (full gf) | no | — |
| **Eu II** | **vald(bsyn harnesses)** | 151,153 | **fraction folded in** | **YES** | **+0.3002 (blended feature)** |
| **CH** | **vald(bsyn harnesses)** | ¹²C,¹³C | **fraction folded in** | **YES** | +0.0048 (¹²CH) … +1.9586 (¹³CH) |
| OH / MgH / SiH / C₂ / CN / TiO | vald(bsyn harnesses) | single isotopologue | undecidable by this test | bounded | +0.0010 … +0.1319 |
| CO | vald(bsyn harnesses) | 12016,12017,12018,13016 | ambiguous | bounded | +0.0058 … +3.4250 (¹³CO) |

### Two premise corrections to the ticket

* **HFS is not isotope structure.** The ticket names Mn and Co as plausibly exposed. They
  cannot be: **Mn-55, Co-59 and Sc-45 are 100 % abundant and V-51 is 99.75 %**, so
  `Σf² = 1` and the offset is identically zero no matter how hyperfine-split the feature
  is. Neither element is isotope-coded anywhere in the pipeline. *Any* reasoning of the
  form "this element has HFS, therefore it is exposed" is invalid.
* **The two surfaces that feed live values are clean.** The GES v6 HFS/ISO list that
  `pipeline/abundances_derive` hands to iSpec, and the Gerber NLTE deck that Engine-B
  synthesises on, are both fraction-free on every species tested — i.e. correctly form (A).
  Only the TSFitPy `linelist_vald` "for-grid" files are defective.

### Single-isotope coverage gaps (separate, small, worth recording)

Two GES-surface species do not code a complete isotope set, so the engine applies less than
unity: **Ba II Σf = 0.99793** (¹³⁰Ba/¹³²Ba absent, −0.009 dex) and **Sm II Σf = 0.97237**
(−0.0122 dex). Both are far below anything this project acts on and neither is a
double-application; recorded so they are not rediscovered as one.

---

## 3. Live values: correction magnitude and disposition impact

**No reported abundance is affected. Nothing needs re-deriving and no disposition changes.**
Traced element by element:

| element | live value | route | exposed? | correction | disposition |
|---|---|---|---|---|---|
| **Eu II** | none (owed-no-value, PR#188) | RYA-565 ran BOTH legs | the VALD leg is, by +0.3002 | n/a — the GES leg is the adopted one | **unchanged.** Eu is owed because `cog_linearity` 0.953 puts it on the linear CoG, not because of gf |
| **Ba II** | 2.237, PASS (PR#190, merged) | `rya581_ba2_deblend_sirius.py` writes its OWN Ba HFS block as **`56.000`** and drops all seven VALD Ba II blocks from the blend | **no** | none | **unchanged** — the PASS flip stands |
| **Cu I** | RYA-466 HFS synthesis | GES/iSpec path (`GESv6_atom_hfs_iso`) | **no** | none | unchanged |
| **Li I** | UPPER LIMIT | EW pool → GES path; Li I is not isotope-coded on the surfaces used | **no** | none | unchanged |
| **Ca** | owed (EW pool, Ca **I**) | Ca I is not isotope-coded on any surface; only Ca **II** is coded, and no harness targets it | **no** | none | unchanged |
| **Mn, Co, Sc, V** | Mn 5.466 PASS, Co 4.960 | mono-isotopic — structurally immune | **no** | none | unchanged |

### Blend-model exposure (the second-order path)

Every bsyn harness copies the in-window VALD blocks verbatim as its blend model, so an
exposed block can bias a target it does not belong to. The audit scores each fit window by
the fraction of its line strength that is missing, weighting each line by `gf·10^(−θ·EP)`
with θ = 5040/5777:

| window | missing | note |
|---|---|---|
| RYA-551 Sr II 4305.443 | **12.376 %** | CH G-band head. This line is **EXCLUDED** from the Sr result (dEW/dA = 0.2 mÅ/dex) — touches nothing |
| RYA-560/585 Zr II 4258.041 | **3.798 %** | Zr is **owed** with `reliable_lines: []` — touches nothing |
| RYA-581 Ba II 5853.668 | 2.451 % | entirely the target's own VALD blocks, which the harness **drops** |
| RYA-565 Eu II 6645.064 | 1.965 % | the deliberate control leg |
| RYA-551 Sr II 4215.519 | **0.517 %** | a live CROSS-CHECK line — the largest live-adjacent number in the sweep |
| RYA-551 Sr II 4077.709 | 0.0003 % | the line the reported **A(Sr) 2.759** is measured on |
| all other windows | < 0.021 % | |

The reported Sr value comes from 4077.709, where 0.0003 % of window absorption is missing.
The one number worth naming is the **Sr II 4215.519 cross-check at 0.517 %**; it does not
set the value and is recorded rather than acted on.

The mechanism there is ¹³CH: folded twice, a ¹³CH line is synthesised 1.96 dex weak, i.e.
effectively deleted from the blend model. That matters only where CH dominates the window,
which is the 4300 Å G band.

---

## 4. The enforceable convention

Stated in `docs/CONVENTIONS.md` and enforced by
`pipeline/isotope_gf_convention.assert_target_convention(linelist, Z, ion)`, which a bsyn
harness calls immediately before handing its line list to the engine. It raises if the
species being FITTED is isotope-coded on a surface known to ship folded gf, and records —
without failing — exposed blocks that are only in the blend model, which is the distinction
the measurements above justify.

The guard does not hardcode which species are folded. It reads them from
`data/audit/rya684_isotope_gf_audit.json`, which `scripts/rya684_isotope_gf_audit.py`
regenerates by measuring the shipped files, so re-vendoring a line list and re-running the
audit keeps the guard true. `tests/test_isotope_gf_convention_rya684.py` pins the measured
verdicts, asserts the two live-value surfaces stay fraction-free, and fails if a third
target species ever starts being fitted against its own exposed block.

---

## 5. Open follow-ups (not done here, not blocking)

1. **Molecular isotopologue convention is undecided for OH, MgH, SiH, C₂, CN, TiO.** Each
   codes a single isotopologue, so the cross-isotope reconstruction cannot discriminate.
   The bound is small for most (≤0.035 dex) but **MgH +0.1024** and **TiO +0.1319** are not
   negligible if either is ever used quantitatively. CH is decided (folded) and CO is
   ambiguous with a ¹³CO bound of +3.425 dex, which matters to any future ¹²C/¹³C work
   (RYA-503).
2. **Fix the shipped lists, or keep routing around them.** Nothing today needs the VALD
   for-grid isotope blocks, so the cheapest correct answer is the guard. If a future
   measurement wants them, divide the fractions out at vendor time rather than at use time.
3. **Ca II's coded components sit 219 mÅ apart** in the VALD list, which is larger than
   published Ca isotope shifts. Unused today; flagged as a line-data question, not a
   convention one.
