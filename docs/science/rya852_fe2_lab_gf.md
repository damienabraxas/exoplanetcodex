# RYA-852 — can the Fe II arbiter lines be graded?

**Status:** audit. No product changed. **Not merged — Ryan reviews.**

```
python3 scripts/rya852_fe2_lab_gf_audit.py     # Sirius (needs astroquery in venv_ci)
```

---

## The answer is no, and the reason is worse than a coverage gap

**None of the three Fe II arbiter lines has a confirmed primary-laboratory gf** — and the
grade already stored for two of them is wrong.

### Which lines — confirmed from the products, not the guess

The ticket guessed *"6247.6 / 6432.7 / 6456.4 or similar; confirm from linelist"*. The
three-line set is the **Fe II VIS ENGINE-A aggregate**:

```
6147.734   6238.386   6247.557        (EP 3.889 / 3.889 / 3.892 eV)
```

Only 6247.6 was in the guess; **6432.7 and 6456.4 belong to the wider 11-line VIS 1D-LTE
pool**, not the trio.

---

## 🔴 The stored grade disagrees with NIST on both the grade and the value

| line | `canonical_gf` says | log gf stored | NIST log gf | Δ | **NIST accuracy** |
|---|---|---|---|---|---|
| 6147.734 | `RU` (Raassen & Uylings) | −2.827 | −2.796 | −0.031 | **E** (0.301 dex) |
| 6149.246 | **"NIST ASD v5.11 grade B"** | −2.724 | −2.854 | **+0.130** | **E** (0.301 dex) |
| 6238.386 | `2009A&A…` (Meléndez & Barbuy) | −2.600 | −2.754 | **+0.154** | **D** (0.176 dex) |
| 6247.557 | **"NIST ASD v5.11 grade B"** | −2.329 | −2.444 | **+0.115** | **D** (0.176 dex) |

Grade B is 10% → 0.041 dex. NIST actually grades these **D (50% → 0.176)** and
**E (100% → 0.301)**:

- **6149.246** — stored `B`, NIST says `E`: **understated 7.3×**
- **6247.557** — stored `B`, NIST says `D`: **understated 4.3×**

And in both cases the stored *value* isn't NIST's either, sitting +0.115/+0.130 dex above.
So the label "NIST ASD v5.11 grade B" is wrong on **both** counts.

**Why this matters now:** RYA-850 keys its graded gf term on exactly this metadata. Wiring
Fe II into the graded pool on the stored grade would publish **0.041 dex** on lines whose
own source says **0.176–0.301** — the precise failure this ticket forbids.

### Verdict

The trio stays **UNGRADED**, with the reason stated (*"do not fabricate coverage"*), and the
stored grade B is a **defect to correct**. Their cited accuracies span **0.176–0.301 dex**,
worse than the ~0.1 dex floor the ticket anticipated — ionized Fe really is hard to measure.

---

## The pool sits ~0.1 dex above NIST as a whole — and that lands on the ionization balance

This is not confined to the two mislabelled lines. Across all 9 matched pool lines:

**median offset +0.106 dex, 7 of 9 above NIST** (range −0.031 … +0.154) — *including its
plain-VALD3 members*. A coherent offset across an entire pool is a **scale** difference,
not nine independent errors.

A gf that is too high yields an abundance that is too low by about the same amount, so:

| | Fe I VIS | Fe II VIS | Fe I − Fe II |
|---|---|---|---|
| current scale | 7.586 | 7.568 | **+0.018** |
| on NIST's scale | 7.586 | ~7.674 | **−0.088** |

⚠️ **The current scale is what makes the solar ionization balance work.** Adopting NIST's
values would break it by ~0.1 dex. Given that one pool member is labelled MB09 outright and
MB09 is partly built by reverse solar analysis, *"a gf scale that makes the solar ionization
balance come out right"* is exactly what a solar-fitted scale would look like — the RYA-161
circularity, showing up as a property of the data rather than a choice anyone made.

**That is a hypothesis, not a finding.** It needs MB09's own S/L flags to test, which this
audit could not obtain.

---

## Two traps, both of which nearly produced a wrong answer

⚠️ **`astroquery.nist` defaults to `wavelength_type='vacuum'`.** Queried that way, **none**
of the three arbiter wavelengths appear in the Fe II list, and the obvious conclusion is
"NIST doesn't cover these lines". They are all there — air→vacuum is **+1.71 Å** at 6150 Å,
so 6147.734 is listed at 6149.435. Same class as RYA-846's Wallace atlas: a unit convention
manufacturing an absence.

⚠️ **Match on wavelength *and* EP.** A ±0.05 Å window alone returns EP **13.436 eV** for
6149.246 and **10.930 eV** for 6432.676, against true values of 3.889 and 2.891 — high-
excitation neighbours (RYA-780). I flagged this for the `canonical_gf` side and then **made
the same mistake on the NIST side**, which produced a nonsense −2.298 dex offset for
6432.676 until EP matching was added there too. The guard has to be on **both** sides of a
cross-match.

---

## What this does NOT establish

1. **Den Hartog 2019's optical subset could not be reached.** `Vizier.find_catalogs`
   returned nothing for three phrasings. That is an absence in **the search**, not in the
   source (RYA-833) — whether these three lines are among its ten optical lines is
   **unverified**, and it remains the most promising route to a genuine grade.
2. **MB09's S/L flags were not obtained**, so the RYA-161 firewall cannot be applied line by
   line. 6238.386 is MB09-labelled and is a firewall **candidate**, not a confirmed
   exclusion.
3. The consistent +0.115…+0.154 dex offset across the three non-`RU` lines *looks* like one
   common scale (plausibly MB09), but that is inference from the pattern.

---

## Recommended next steps

1. **Correct the stored grade** for 6149.246 and 6247.557 — `B` → the NIST-cited `E`/`D`.
   This is a data defect independent of everything else here, and RYA-850 will consume it.
2. **Obtain Den Hartog 2019** (ApJS) directly rather than through VizieR search; its ten
   optical Fe II lines are the one plausible path to a graded Fe II arbiter.
3. **Get MB09's S/L flags** to settle the firewall question and test the scale hypothesis.
4. Treat the **+0.106 dex pool-wide offset** as a finding in its own right — it is larger
   than the Fe I − Fe II balance it underwrites.
