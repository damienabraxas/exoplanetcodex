# RYA-853 — do the stored gf grades and values match the source they cite?

**Status:** audit. Nothing corrected yet. **Not merged — Ryan reviews.**

```
python3 scripts/rya853_gf_grade_integrity_audit.py     # Sirius (astroquery in venv_ci)
```

---

## The headline: it is not two lines, it is most of the file

The two hand-maintained NIST extracts (`nist_reference.csv`, `nist_crosscheck.csv`) carry
60 rows. Verified against **live NIST ASD**, in air, matched on wavelength **and** EP:

| | |
|---|---|
| rows | 60 |
| uniquely matched to NIST | **44** |
| ambiguous (>1 NIST row) — **not judged** | 8 |
| no NIST row / query failed | 8 |
| **stored grade ≠ NIST accuracy** | **31 of 44 (70%)** |
| **stored log gf ≠ NIST (>0.02 dex)** | **14 of 44** |
| **the two files disagree with each other** | **13** (offline, no NIST needed) |

`nist_reference.csv` describes itself as the *"Type B uncertainty anchor"*. Seven in ten of
its verifiable grades disagree with the source they cite.

### It is a known defect class, found once and fixed for two rows

The file's own header records it:

> *RYA-592 re-verification (2026-07-26, NIST ASD v5.12, **Mg I rows 5528/5711 ONLY**): the
> log_gf values RE-DERIVE EXACTLY from the live source, but **two columns were wrong** … the
> **accuracy CODE** (both rows carried A; ASD reports B+ and B) and aki_s-1 … **The OTHER
> rows in this file are NOT re-verified — the drift found here means they should not be
> assumed current.**"

This audit is what that warning asked for. The Fe II rows RYA-852 found are simply the next
instance.

---

## Worst offenders (understatement of the cited uncertainty)

| line | stored | NIST | stored dex | NIST dex | understated |
|---|---|---|---|---|---|
| Mg I 4571.096 | `A` | **D** | 0.013 | 0.176 | **13.7×** |
| Fe I 5576.089 | `A` | **D** | 0.013 | 0.176 | **13.7×** |
| Ca I 5857.451 | `A` | **D** | 0.013 | 0.176 | **13.7×** |
| Ca I 6439.075 | `A` | **D** | 0.013 | 0.176 | **13.7×** |
| O I 6300.304 | `A+` | **C+** | 0.009 | 0.072 | **8.4×** |
| Si I 5793.07 | `B` | **E** | 0.041 | 0.301 | **7.3×** |
| Fe II 6149.25 | `B` | **E** | 0.041 | 0.301 | **7.3×** |
| Ca I 6102/6122/6162 | `A` | **C** | 0.013 | 0.097 | **7.5×** |
| Fe II 6247.56 | `B` | **D** | 0.041 | 0.176 | **4.3×** |

**Why this is not cosmetic:** RYA-850 keys `graded_gf_term` on exactly this metadata. A
stored `A` on a line NIST grades `D` publishes **0.013 dex** where the source says **0.176**.
A grade is a claim about how well a number is known; a wrong one is an understated error bar
with a citation attached.

⚠️ Three lines run the *other* way — Li I 6707.76/6707.91 stored `A`/`A+` where NIST says
`AAA`, i.e. our bar is **too wide**. The defect is not one-directional, which rules out a
simple "someone rounded optimistically" story.

---

## 🔴 The two files disagree with each other — and that needs no network

13 disagreements on lines both files carry, detectable offline on any laptop:

| line | `nist_reference` | `nist_crosscheck` |
|---|---|---|
| **Mg I 5711.090** | `B` | `A` |
| Fe I 5576.090 / 6065.49 / 6136.62 | `A` | `B` |
| O I 6300.304 | `A` | `A+` |
| Li I 6707.760 | log gf **+0.177** | **−0.002** |
| S I 6052.670 | log gf **+0.370** | **−0.550** |
| Ni I 6300.336 | log gf **−2.110** | **−2.310** |

**Mg I 5711.090 is the smoking gun**: that is precisely the row RYA-592 corrected — in
`nist_reference` only. The twin file was left stale, so the same line now carries two
different grades depending on which file you read.

Two consequences worth flagging on their own:

- **S I 6052.670 differs by 0.92 dex between our own two files.**
- **Ni I 6300.336 differs by 0.20 dex** — that is the Ni blend under [O I] 6300, the line
  RYA-365 rerouted to gate the solar oxygen abundance.

---

## Scope 4 — the cross-match guards

| site | air/vac | EP |
|---|---|---|
| `rya822_pull_nist_nearuv.py` | **OK** (`vac+air`) | **OK** (RYA-780, both sides) |
| `rya853_gf_grade_integrity_audit.py` | OK | OK |
| `rya347_fe2_atomic_data_audit.py` | n/a (local file) | 🔴 **MISSING** — wavelength-only ±0.1 Å |

⚠️ **My first diagnosis was wrong and is recorded because it is the obvious one.** The
missing EP guard in `rya347` looks like the cause of the Fe II defect. **It is not** — the
EP in these files is *correct* (3.889 / 3.892), the right physical line was matched, and the
stored grade and value simply are not NIST's. The loose ±0.1 Å window is a latent risk worth
closing separately; it did not produce this.

⚠️ **And I had to guard my own audit.** A wavelength+EP window is *not* a unique line
identifier — several transitions share a wavelength and a lower level while differing in the
upper level. Taking the first match manufactured 12-dex "defects" (Mg I 5183.604 stored
+0.180 against a NIST row at −11.908; O I 6300.304 against −12.201). Those 8 rows are now
reported **AMBIGUOUS and not judged**, because the stored files carry no upper level. Every
number above is from a *uniquely* matched row.

---

## Scope 2 — the Fe I lab pool: materially CLEAN, and that is the important result

The NIST pass covered ~60 rows of two hand-maintained extracts. It did **not** touch the
pool RYA-850 actually promotes: the 465 Fe I lines graded on primary lab sources. Those
grades come from the lab papers, not NIST, so a NIST cross-check could never have cleared
them — and after a 70% failure next door there was no basis for assuming they were exempt.

The papers are on the Mac (`Reference documents/`), read off disk. Refereed against the
tables **present in the PDFs**:

| source | refereed | of pool | coverage |
|---|---|---|---|
| Belmonte 2017 | 118 | 120 | **98.3%** |
| Den Hartog 2014 | 64 | 203 | 31.5% |
| Ruffoni 2014 | 35 | 142 | 24.6% |
| **total** | **217** | **465** | **46.7%** |

### The verdict

**207 of 217 refereed lines (95.4%) reproduce their source paper's log gf within 0.02 dex.**

That is the opposite of the NIST extracts. The pool RYA-850 promotes is sound where it can
be checked.

**10 value mismatches**, nine of them ≤0.10 dex. One real outlier:

> 🔴 **Belmonte 3935.307 — ours −2.199, paper −1.820 (Δ +0.379)**, and its stored σ is wrong
> too (0.070 against the paper's 0.180). Both axes disagree on the same line, which is the
> signature of a bad row rather than rounding.

**44 cited-σ mismatches**, almost all at the ±0.01 rounding level. Den Hartog's run
consistently **ours = paper + 0.01** (0.030 vs 0.020, 0.040 vs 0.030) — our bar is
*conservative* there, not understated. Belmonte's go both ways.

⚠️ **53% of the pool is unverified, not clean.** Ruffoni's full line list is Table 3, whose
caption reads *"Only the first upper level is shown here. The full [table is online]"*, and
Den Hartog 2014 is the same. The unrefereed remainder is absent from this check, **not shown
to be correct** (RYA-833).

### 🔴 Three encoding traps, each of which produced a confident wrong answer

1. **U+2212 MINUS SIGN, not ASCII `-`.** An ASCII `-?` silently fails, the capture starts at
   the digit, and every negative log gf comes back **positive**. First run: **95 of 99 lines
   "mismatched" by exactly twice their value.** A defect rate that high with a regular
   signature is a parser bug, not a finding.
2. **U+F0A0 (Private Use Area) padding around Belmonte's `±`** — `0.43\uf0a0±\uf0a00.02`. A
   `\s*±\s*` pattern matches nothing, and Belmonte looked like it covered **none** of our
   120 lines. A manufactured absence caused by a font quirk.
3. **`pdftotext -layout` is required.** pypdf's `extract_text` collapses the columns and the
   tables become unparseable — which is why the Mac run found 0 Belmonte rows.

Belmonte's wavelengths are also in **nanometres** while the Wisconsin tables are in air
Ångströms; the conversion is explicit and coverage is reported so a wrong convention shows
up as zero matches rather than a clean pass.

---

## Scope 3 — the DH19 referee: INCONCLUSIVE, and the question was malformed

Ryan supplied Den Hartog 2019's ten optical Fe II lines (BF × LIF lifetimes, pure lab, no
solar normalisation) — the referee RYA-852 needed.

| comparison | median | 95% CI | sd | n |
|---|---|---|---|---|
| **ours − DH** | **+0.020** | **[−0.070, +0.160]** | 0.179 | 10 |
| NIST − DH | +0.124 | [−0.019, +0.213] | 0.145 | 10 |
| ours − NIST *(same 10 lines)* | **−0.066** | [−0.190, +0.080] | 0.220 | 10 |
| ours − NIST *(RYA-852, red pool)* | +0.106 | — | — | 9 |

**The counter-evidence is confirmed in direction.** NIST − DH = **+0.124**: NIST's Fe II gf
sit *above* pure lab, which gives the *lower* abundance Ryan quoted (NIST 7.31 vs DH 7.46).
Higher gf → lower abundance, so the signs are consistent.

### But the test cannot choose, and saying so is the result

`ours − DH = +0.020` looks like a clean "our scale is the lab scale". It is not: the 95% CI
**[−0.070, +0.160] covers both the pure-lab prediction (~0) and the solar-fitted one
(~+0.13)**. Ten lines spanning −0.360…+0.260 dex cannot separate them.

⚠️ My first pass thresholded on `|median| ≤ 0.05` and returned **"LEGITIMATE — hypothesis
REFUTED"**. That was overreading a median with an uncertainty three times its size. The CI
is what caught it, and it is now the verdict rule.

### 🔴 And the premise itself does not hold: the offset is band-dependent

| band | ours − NIST |
|---|---|
| blue overlap 4173–4584 Å | **−0.066** |
| red pool 5256–6456 Å (RYA-852) | **+0.106** |

**Swing +0.172 dex, and the sign flips.** So *"the Fe II pool sits +0.106 above NIST"* is
not a scale property — RYA-852 measured it on the **red pool alone**, and the blue lines
carry different provenance (`T83av`/`PGHcor`/`BSScor`/`KK` against `VALD3`/`RU`/MB09).
There is no single "our Fe II scale" to referee, which is why a ten-line overlap in one band
was never going to settle a claim made in another.

### What would settle it

The **full DH19 Table 6 — 131 lines**, VizieR `J/ApJS/243/33`. ⚠️ That ID returns **0 tables**
from Sirius via astroquery (tried both `find_catalogs` and `Vizier(catalog=…).get_catalogs`),
so the ten values here are **transcribed from the PDF** and this run is not network-
reproducible. A direct CDS fetch would give 131 lines and a per-band offset with real
statistics.

⚠️ **The arbiter lines are not refereed by this.** DH19's optical set stops at 4584 Å; the
three Fe II arbiter lines (6147.734 / 6238.386 / 6247.557) are redward of it. This tests the
*scale* on the blue overlap, not those lines.

**Status of the two hypotheses: both still open.** RYA-852's solar-fitting reading is neither
supported nor refuted, and the NIST-is-the-outlier reading is confirmed in *direction*
(NIST − DH > 0) but not in magnitude against our own pool.

## Recommended order

1. **Correct the 31 grade mismatches and 14 value mismatches**, and reconcile the two files
   against each other. Regenerate anything keyed on them.
2. **Add the offline extract-vs-extract check to CI** — it caught RYA-592's half-applied fix
   with no network and would have caught it a month ago.
3. **Fetch DH19 directly** and run the referee test; until then RYA-850's graded bars stay
   provisional, and the RYA-852 scale hypothesis stays open.
4. Close the `rya347` ±0.1 Å wavelength-only match on general principle.
