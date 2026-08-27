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

## Scope 2 — the Fe I lab pool: CLEAR. This is what unblocks RYA-850.

The NIST pass never touched the pool RYA-850 promotes. Refereed against the **machine-
readable CDS tables** (the full lists the PDFs only excerpt), vendored into
`data/reference/fe_gf_lab/cds/` so the audit reproduces offline:

| source | refereed | of pool | coverage | value mismatches | σ mismatches |
|---|---|---|---|---|---|
| **Ruffoni 2014** *(CDS `J/MNRAS/441/3127` table3)* | 142 | 142 | **100%** | **0** | **0** |
| **Den Hartog 2014** *(CDS `J/ApJS/215/23` table4)* | 203 | 203 | **100%** | **0** | **0** |
| Belmonte 2017 *(PDF excerpt — no CDS table)* | 118 | 120 | 98.3% | 8 | 9 |
| **total** | **463** | **465** | **99.6%** | **8** | **9** |

### 455 of 463 lines (98.3%) reproduce their source paper exactly

**Ruffoni and Den Hartog are perfect — 345 lines, zero defects on either value or cited σ.**
That is the opposite of the NIST extracts (70% wrong) and it clears the pool RYA-850
promotes.

**Every mismatch is Belmonte**, and Belmonte is the only source refereed from a *typeset
PDF table* rather than a machine-readable one (`J/ApJ/848/125` is a 404 on CDS). Its table
prints 2 decimals while our stored values carry 3, so the small deltas may be extraction or
rounding artifacts rather than defects. One is not:

> 🔴 **Belmonte 3935.307 — ours −2.199, paper −1.820 (Δ +0.379)**, with its σ wrong too
> (0.070 vs 0.180). Too large for rounding, and wrong on both axes — a genuine bad row.

⚠️ The other seven (≤0.10 dex) and the nine σ differences (±0.01–0.02) **cannot be
adjudicated without Belmonte's machine-readable table.** They are flagged, not condemned.

### 🔴 Four parser bugs, every one of which produced a confident wrong answer

Each was caught only by checking a value that was already known:

1. **U+2212 MINUS SIGN, not ASCII `-`.** An ASCII `-?` fails silently, the capture starts at
   the digit, and every negative log gf returns **positive** — first run: **95 of 99 lines
   "mismatched" by exactly twice their value.**
2. **U+F0A0 (Private Use Area) padding around Belmonte's `±`** — `0.43␣±␣0.02`. A `\s*±\s*`
   pattern matched nothing and Belmonte looked like it covered **none** of our 120 lines.
3. **`pdftotext -layout` is mandatory** — pypdf collapses the columns, which is why the Mac
   run found zero Belmonte rows.
4. **Den Hartog column off-by-one.** Field 7 is `A_ul`, not `log gf`. Reading it as log gf
   reported "paper +35.503" against our −0.310 and turned **211 of 463 lines** into
   mismatches. Caught by spot-checking one line whose value was already known from the PDF.

A defect rate that is suddenly enormous, with a regular signature, is a parser bug — not a
finding. Three of these four would have condemned a clean pool.

⚠️ Belmonte's wavelengths are in **nanometres**; the Wisconsin tables are in air Ångströms.

## Scope 3 — the DH19 referee: **LEGITIMATE**, at n=22, after the test was repaired

> **Re-run 2026-08-27.** The August run of this scope reported INCONCLUSIVE on ten lines
> transcribed from the DH19 PDF. Both the input and the verdict have changed, and the reason
> the verdict changed is not that the evidence moved — it is that the test had broken and had
> to be fixed first. That story is the important part of this section.

### 🔴 The referee had come to sit on both sides of its own comparison

[RYA-945] (`b545dc6`, 2026-08-21) ingested Den Hartog 2019 Table 6 into `canonical_gf`.
Correct for the line list. Fatal for this experiment: the referee read *ours* from
`canonical_gf`, so by 2026-08-27 all ten overlap lines cited `PRIMARY LAB DenHartog2019` and
the script was comparing DH19 against DH19.

| | Aug 17 (committed) | Aug 27 (re-run, unrepaired) |
|---|---|---|
| ours − DH, median | +0.020 | **+0.000** |
| 95% CI | [−0.070, +0.160] | **[+0.000, +0.000]** |
| sd | 0.179 | **0.000** |
| our source, all 10 lines | VALD3 / T83av / PGHcor / BSScor / KK | **`PRIMARY LAB DenHartog2019`** |
| verdict | INCONCLUSIVE | **LEGITIMATE … REFUTED** |

The verdict flipped from INCONCLUSIVE to REFUTED with nobody touching RYA-853. **The
bootstrap CI added in August specifically to stop the median being overread did not catch
it** — a zero-width CI passes every "is this well determined?" test there is. A guard against
overreading a *noisy* statistic is silent about a *collapsed* one.

### The repair: freeze the scale being refereed, and refuse a comparison that is not a test

1. **`ours` is now a frozen artifact**, `data/reference/fe_gf_lab/fe2_pre945_scale_snapshot.csv`
   — Fe II as it stood at `b545dc6^`, before the ingest. That is the scale that underwrote
   the ionization balance, and the only one the referee can test. Post-945 the overlap lines
   *are* Den Hartog's by adoption, and asking whether they agree with Den Hartog is not a
   question.
2. **Guard 1 — self-reference.** Any line whose stored reference cites the referee is excluded
   and counted. On the frozen snapshot: 0 excluded. On live `canonical_gf`: **22 of 22**.
3. **Guard 2 — degenerate estimator.** A verdict is refused when the CI width or the sd
   collapses, or fewer than 5 lines survive. `--source live` reproduces the broken run on
   purpose and now returns **`UNUSABLE — the comparison is not a test`**.
4. **The referee is read from the repo, not transcribed.** RYA-945 vendored the
   machine-readable Table 6 with a recorded PDF sha256, so this is **131 lines, reproducible
   offline** — the limitation the August run declared is closed.

### The result, on 22 EP-matched lines

| comparison | median | 95% CI | sd | n |
|---|---|---|---|---|
| **ours − DH** | **−0.035** | **[−0.050, +0.030]** | 0.135 | 22 |
| NIST − DH | +0.021 | [−0.038, +0.060] | 0.136 | 22 |
| ours − NIST *(same lines)* | −0.036 | [−0.092, +0.054] | 0.166 | 22 |
| ours − DH, **near-UV** 3003–3277 Å | −0.035 | [−0.055, −0.003] | 0.085 | 12 |
| ours − DH, **blue** 4173–4584 Å | +0.020 | [−0.070, +0.160] | 0.179 | 10 |

**VERDICT: LEGITIMATE — our scale IS the pure-lab scale.** The CI **excludes the
solar-fitted prediction (~+0.13)**, which is what a verdict requires; a median inside a
threshold is not enough, and that was the August mistake.

⚠️ **The blue arm did not change its mind — it is still INCONCLUSIVE on its own** (CI covers
both readings, exactly as in August). The twelve near-UV lines that RYA-945's full table
added are what carry the verdict. Pinned in the tests so this cannot be misread later.

### 🔴 The premise still does not hold: the offset is band-dependent

| band | ours − NIST | n |
|---|---|---|
| blue overlap 4173–4584 Å | **−0.066** | 10 |
| red pool 5257–6456 Å | **+0.107** | 8 |

**Swing +0.173 dex, sign flips.** Both arms are now *measured on this run against live NIST
from the same source* — the previous revision carried RYA-852's `+0.106` as a hardcoded
literal, and quoting a number taken before 852/877/945 touched the pool beside a freshly
measured one compares two pools on two dates and calls the difference a finding. (The
re-measured red arm, +0.107, happens to land on RYA-852's +0.106 — but it is now a
measurement, and three lines are dropped as ambiguous rather than argmin-matched.)

### What this does and does not clear

✅ The Fe II **gf scale** is on the laboratory scale, so *sitting above NIST is NIST being
low*, as Den Hartog 2019 reports, and RYA-852's solar-fitting hypothesis is **refuted for the
scale**.

⚠️ **It does not clear the arbiter lines.** DH19 stops at 4584 Å; 6147.734 / 6238.386 /
6247.557 are redward of it. And `canonical_gf` still carries the fabricated `NIST ASD v5.11
grade B` on **6149.246** and **6247.557** — the defect that opened this ticket is unfixed.

⚠️ **Independence is established only as deep as the reference string.** The twelve near-UV
lines carrying the verdict are plain `VALD3`, `single_source`. VALD3 is an aggregator, so
that label does not prove independence from Den Hartog — only that we did not adopt him
directly. The values differ by up to 0.250 dex so this is not a laundered copy, but a shared
upstream ancestor cannot be excluded from what we record. Closing it needs VALD3's per-line
source for those twelve.

🔴 **It does not clear the ionization balance, because the measurement leg is pre-continuum-fix.**
The balance the verdict speaks to is Fe I 7.586 / Fe II 7.568 (RYA-783), measured on
`kpno_solar_atlas` PROFILEFIT products dated **2026-08-16…18**. Every continuum fix postdates
them:

| fix | landed |
|---|---|
| RYA-911 — HARPS Fe II EW legs 0.34 dex low, continuum placement a partial cause | 2026-08-19 |
| RYA-913 — ENGINE-B hardcoded `load_kp_window`, produced the retracted 7.486 | 2026-08-19 |
| RYA-1000 / RYA-1006 — `--local-renorm` and the conditioning axis in the artifact stem | 2026-08-23 |
| RYA-1026 — KP class ratified `pre_normalised`, `prenormalised_guard` wired | 2026-08-24 |
| RYA-1030 — normalisation detected off the flux, not the flag | 2026-08-24 |

**The Fe II VIS cell has never been re-run.** There is no Fe II product in
`data/results/band_products` at all — every current Fe product is Fe I. The most recent Fe II
numbers are the 2026-08-18 RYA-877/880 artifacts, which already moved the 1D-LTE cell
7.568 → **7.542** when the circular MB09-S line 5991.371 was dispositioned; the RYA-852 test
still hardcodes 7.568. A gf scale cleared here is one leg of the balance, and the other leg
is stale.

## Recommended order

1. ~~Fetch DH19 directly and run the referee test~~ — **done 2026-08-27**, and the test had
   to be repaired before it could be run at all. See scope 3.
2. **Re-run the Fe II VIS cell.** The gf scale is cleared; the measurement leg is not. Nothing
   downstream of the ionization balance should be treated as settled until an Fe II product
   exists that postdates RYA-1026/1030, and none does.
3. **Correct the 31 grade mismatches and 14 value mismatches**, and reconcile the two extract
   files against each other. Regenerate anything keyed on them. These rows feed published
   anchors (O I 6300, Li I 6707, Ba II 5853, Ni I 6300.336), so this is a reviewable change
   with a real blast radius — hence tabulated in `rya853_corrections.csv` and not applied.
4. **Fix `canonical_gf` Fe II 6149.246 and 6247.557**, which still read `NIST ASD v5.11
   grade B` where NIST grades them **E** and **D** — the defect that opened this ticket.
5. **Add the offline extract-vs-extract check to CI** — it caught RYA-592's half-applied fix
   with no network and would have caught it a month earlier.
6. Close the `rya347` ±0.1 Å wavelength-only match on general principle.
7. Fix the Belmonte Fe I **3935.307** row (stored −2.199 / σ 0.070; paper −1.820 / 0.180).

## Provenance of this document

| artifact | what |
|---|---|
| `scripts/rya853_freeze_pre945_fe2_scale.py` | freezes the pre-RYA-945 Fe II scale from `b545dc6^` |
| `data/reference/fe_gf_lab/fe2_pre945_scale_snapshot.csv` | that snapshot, 15,280 rows, 0 citing the referee |
| `scripts/rya853_dh19_scale_referee.py` | scope 3, with both guards; `--source live` shows them fire |
| `data/results/rya853/rya853_dh19_referee.json` | the repaired run (verdict LEGITIMATE, n=22) |
| `data/results/rya853/rya853_dh19_referee_live.json` | the degenerate control (verdict UNUSABLE) |
| `data/results/rya853/rya853_red_arm_vs_nist.csv` | the red arm, measured rather than quoted |
