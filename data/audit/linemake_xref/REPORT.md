# RYA-1070 — is `linemake` a good gf cross-reference for the lines we already hold?

READ-ONLY audit. linemake `61af45bd8293c91059360cbb0cfa2d16d34fc21e`, `canonical_gf.csv` sha256 `3f87cfe8c64fa489…` (178,680 rows). No line list was edited, no gf changed, no linemake line imported, no synthesis run.

## The answer

**Yes for corroboration, and more usefully as a POINTER to lab papers we are not citing.** linemake's curated database holds 13,902 comparable transitions across 77 species; 74 of those species appear in our pool, covering 137,130 of our 178,680 canonical_gf lines. Of those, **10,555 matched** on the full physical key — species AND wavelength AND excitation potential — for a coverage of **7.7%**. Of the 10,491 matched lines that carry a directly comparable value, **7,767 = 74.0% agree** — the two lists carry the same number to the precision both of them printed (median threshold 0.0100 dex, derived per line; 6,716 of the matched deltas are *exactly* zero). That is the corroboration answer.

The verdict buckets split those comparable lines differently, and deliberately: `AGREES` (719) is reserved for lines where we already sit at a comparable tier, while **9,539 lines land in `LINEMAKE_STRONGER_SOURCE`** — matched lines where we are still on a Kurucz/VALD/other fallback and linemake carries a value tagged to a primary laboratory measurement. Those are the point of this audit, and they are counted separately from agreement precisely so that a line agreeing with linemake is not mistaken for a line that is already well sourced.

The strongest internal check passes. Our Fe I LAB-tier lines and linemake's Fe I list ultimately trace the *same* Wisconsin papers (Ruffoni 2014, Den Hartog 2014, Belmonte 2017), so they have to agree; 396 of them matched with a median |Δlog gf| of **0.0000 dex**. A large systematic offset there would have meant a match-key or air/vacuum bug rather than a real disagreement, and there is none. The wavelength scale was measured, not assumed: under the air hypothesis linemake produces 10,031 EP-consistent coincidences within 5 mÅ of our air wavelengths, against 3 under the vacuum hypothesis — linemake is on the air scale above 2000 Å, and its 214 lines below that boundary were excluded from the numeric match rather than matched raw.

## What we could promote — `LINEMAKE_STRONGER_SOURCE` (9539 lines)

These are lines where our `gf_tier` is a compilation fallback (KURUCZ / VALD3 / OTHER) and linemake carries a value tagged to a primary laboratory measurement. **This ticket does not promote them** — RYA-161 validate-don't-tune. It records where a promotion is available, via the paper linemake points at, as future work.

`linemake primary source` lists every paper the README's cell for that species cites, verbatim and in order — including papers cited for context rather than as the adopted source (the `V I` cell, for instance, names Holmes et al. 2016 while explaining that Wood et al. showed the Lawler values are correct). The full unedited cell for every species is preserved in `linemake_readme_sources.csv`; the adopted value for any one line is better identified by its per-line `linemake_source_tag` in `per_line_xref.csv`.

| species | lines | linemake primary source |
|---|---:|---|
| Ce II | 1320 | Lawler et al. (2009, ApJS, 182, 51) |
| Sm II | 931 | Lawler et al. (2006, ApJS, 162, 227); Roederer & Lawler (2012, ApJ, 750, 76) |
| Ti I | 845 | Lawler et al. (2013, ApJS, 205, 11) |
| V I | 754 | Lawler et al. (2014, ApJS, 215, 20); Wood et al. (2018, ApJS, 234, 25); Holmes et al. (2016, ApJS, 224, 35) |
| Nd II | 706 | Den Hartog et al. (2003, ApJS, 148, 543); Roederer et al. (2008, ApJ, 675, 723) |
| Co I | 609 | Lawler et al. (2015, ApJS, 220, 13) |
| Gd II | 497 | Den Hartog et al. (2006, ApJS, 167, 292) |
| Dy II | 329 | Sneden et al. (2009, ApJS, 182, 80) |
| Er II | 307 | Lawler et al. (2008, ApJS, 178, 71) |
| Ti II | 305 | Wood et al. (2013, ApJS, 208, 27) |
| Ni I | 247 | Wood et al. (2014, ApJS, 211, 20) |
| Cr II | 230 | Lawler et al. (2017, ApJS, 228, 10); Nilsson et al. (2006, A&A, 445, 1165); Gurell et al. (2010, A&A, 511, A68); Ward et al. (2023, ApJ, 959, 8) |
| Zr II | 215 | Ljung et al. (2006, A&A, 456, 1181); Malcheva et al. (2006, MNRAS, 367, 754) |
| Cr I | 213 | Sobeck et al. (2007, ApJ, 667, 1267) |
| Pr II | 200 | Sneden et al. (2009, ApJS, 182, 80); Ivarsson et al. (2001, Phys. Scr., 64, 455); Li et al. (2007, Phys. Scr., 76, 577); Biemont et al. (2003, Eur. Phys. J. D, 27, 33) |
| Tb II | 181 | Lawler et al. (2001, ApJS, 137, 341) |
| Sc I | 179 | Lawler et al. (2019, ApJS, 241, 21) |
| Th II | 154 | Nilsson et al. (2002, A&A, 382, 368) |
| Nb II | 150 | Nilsson et al. (2010, A&A, 511, A16) |
| Ca I | 147 | Den Hartog et al. (2021, ApJS, 255, 227) |
| Tm II | 138 | Sneden et al. (2009, ApJS, 182, 80); Wickliffe & Lawler (1997, JOSA B, 14, 737); Den Hartog et al. (2024, ApJS, 274, 9); Kebapci et al. (2024, ApJ, 970, 23) |
| V II | 114 | Wood et al. (2014, ApJS, 214, 18) |
| Fe II | 111 | Den Hartog et al. (2019, ApJS, 243, 33); Meléndez & Barbuy (2009, A&A, 497, 611) |
| La II | 83 | Lawler et al. (2001, ApJ, 556, 452) |
| Os I | 78 | Quinet et al. (2006, A&A, 448, 1207) |
| Y II | 75 | Hannaford et al. (1982, ApJ, 261, 736); Biémont et al. (2011, MNRAS, 414, 3350) |
| Si I | 62 | Den Hartog et al. (2023, ApJS, 265, 42); Pehlivan Rhodin et al. (2024, A&A, 682, A184) |
| Sc II | 57 | Lawler et al. (2019, ApJS, 241, 21) |
| Pt I | 48 | Den Hartog et al. (2005, ApJ, 619, 639) |
| Mn I | 38 | Den Hartog et al. (2011, ApJS, 194, 35) |
| U II | 37 | Nilsson et al. (2002, A&A, 381, 1090) |
| Fe I | 35 | Ruffoni et al. (2014, MNRAS, 441, 3127); Den Hartog (2014, ApJS, 215, 23); Belmonte et al. (2017, ApJ, 848, 126); O'Brian et al. (1991, JOSAB, 8, 1185) |
| Eu II | 24 | Lawler et al. (2001, ApJ, 563, 1075) |
| Ho II | 22 | Lawler et al. (2004, ApJ, 604, 850) |
| Lu II | 19 | Sneden et al. (2009, ApJS, 182, 80); Roederer et al. (2010, ApJL, 714, L123); Roederer et al. (2012, ApJS, 203, 27); Den Hartog et al. (2020, ApJS, 248, 10); Quinet et al. (1999, MNRAS, 307, 934) |
| Ca II | 18 | Den Hartog et al. (2021, ApJS, 255, 227) |
| Mg I | 14 | Pehlivan Rhodin et al. (2017, A&A, 598, A102) |
| Zn I | 10 | Roederer & Lawler (2012, ApJ, 750, 76) |
| Mn II | 9 | Den Hartog et al. (2011, ApJS, 194, 35) |
| Co II | 6 | Lawler et al. (2018, ApJS, 238, 7); Ding & Pickering (2020, ApJS, 251, 24); Roederer et al. 2022, ApJS, 260, 27 |
| Ge I | 6 | Li et al. (1999, PRA, 60, 198); Biémont et al (1999, MNRAS, 303, 721) |
| Ir I | 3 | Xu et al. (2007, JQSRT, 104, 52); Roederer et al. 2022, ApJS, 260, 27; Cowan et al. (2005, ApJ, 627, 238) |
| Os II | 3 | Quinet et al. (2006, A&A, 448, 1207) |
| Pb I | 3 | Roederer & Lawler (2012, ApJ, 750, 76) |
| Si II | 3 | Den Hartog et al. (2023, ApJS, 265, 42); Pehlivan Rhodin et al. (2024, A&A, 682, A184) |
| Ag I | 2 | Hansen et al. (2012, A&A, 545, 31) |
| Ta II | 1 | Quinet et al. (2009, A&A, 493, 711); Morton (2000, ApJS, 130, 403) |
| Yb II | 1 | Sneden et al. (2009, ApJS, 182, 80); Kedzierski et al. (2010, Spectrochimica Acta B, 65, 248) |

Twelve examples, largest |Δ| first — the size of Δ is how much a promotion would actually move the line:

| line_id | species | λ_air (Å) | EP (eV) | our log gf | our source | linemake log gf | Δ | linemake tag |
|---|---|---:|---:|---:|---|---:|---:|---|
| `gf_043533` | Cr II | 6421.699 | 4.476 | -4.987 | RU | -1.720 | -3.267 | Lawler17 |
| `gf_147940` | Ge I | 3067.021 | 2.029 | +0.490 | VALD3 | -2.410 | +2.900 | Li99exp |
| `gf_133676` | V I | 4270.315 | 1.804 | -3.158 | K09 | -0.500 | -2.658 | LAWNOHFS |
| `gf_117406` | Si I | 5597.940 | 4.920 | -2.278 | K07 | -4.861 | +2.583 | Peh24_49% |
| `gf_161716` | Co I | 3562.097 | 2.280 | -3.436 | VALD3 | -0.900 | -2.536 | LAWNHFS |
| `gf_162870` | Ca I | 3607.908 | 1.886 | -2.860 | VALD3 | -5.312 | +2.452 | LAWLER |
| `gf_152297` | Os II | 3213.315 | 1.637 | +1.280 | VALD3 | -1.030 | +2.310 | Quinet06 |
| `gf_148870` | Co I | 3096.404 | 1.785 | -3.631 | VALD3 | -1.500 | -2.131 | LAWNHFS |
| `gf_117514` | Si I | 5872.707 | 4.954 | -2.733 | K07 | -4.841 | +2.108 | Peh24_36% |
| `gf_147202` | Os II | 3042.745 | 1.421 | +1.080 | VALD3 | -0.980 | +2.060 | Quinet06 |
| `gf_117507` | Si I | 5859.201 | 4.954 | -2.407 | K07 | -4.435 | +2.028 | Peh24_99% |
| `gf_016702` | Ca I | 8588.417 | 4.131 | -5.463 | K07 | -7.421 | +1.958 | LAWLER |

## Disagreements for later adjudication — |Δ| > 0.10 dex (56 lines)

By species: Si I (22), Al I (14), Fe II (7), Fe I (6), Mg I (4), Ba II (1), Hf II (1), Sr II (1).

Full list in `disagreements_over_0p10dex.csv`. The 30 largest:

| line_id | species | λ_air (Å) | EP (eV) | our log gf | our source | our tier | linemake log gf | Δ | linemake tag |
|---|---|---:|---:|---:|---|---|---:|---:|---|
| `gf_000028` | Al I | 6906.287 | 4.022 | -1.111 | K75 | KURUCZ | -2.480 | +1.369 | NIST/E |
| `gf_117782` | Si I | 6553.883 | 5.964 | -3.161 | K07 | KURUCZ | -2.450 | -0.711 | NIST/E |
| `gf_117962` | Si I | 7097.473 | 5.984 | -3.918 | K07 | KURUCZ | -3.212 | -0.706 | NIST/E |
| `gf_001465` | Ba II | 4934.074 | 0.000 | +0.542 | 1992A&A... | OTHER | -0.157 | +0.699 | GALNOHFS |
| `gf_118233` | Si I | 8029.166 | 6.083 | -2.286 | K07 | KURUCZ | -1.677 | -0.609 | NIST/E |
| `gf_117937` | Si I | 7016.787 | 5.964 | -2.243 | K07 | KURUCZ | -1.658 | -0.585 | NIST/E |
| `gf_118294` | Si I | 8353.719 | 6.223 | -2.279 | K07 | KURUCZ | -1.750 | -0.529 | NIST/D |
| `gf_117992` | Si I | 7193.553 | 6.079 | -2.284 | K07 | KURUCZ | -1.777 | -0.507 | NIST/E |
| `gf_118038` | Si I | 7289.603 | 6.083 | -1.863 | K07 | KURUCZ | -1.429 | -0.434 | NIST/E |
| `gf_117815` | Si I | 6624.220 | 5.984 | -4.429 | K07 | KURUCZ | -4.020 | -0.409 | NIST/E |
| `gf_000077` | Al I | 8912.900 | 4.085 | -2.348 | K75 | KURUCZ | -1.960 | -0.388 | NIST/C+ |
| `gf_000079` | Al I | 8925.504 | 4.087 | -3.048 | K75 | KURUCZ | -2.660 | -0.388 | NIST/C |
| `gf_000078` | Al I | 8923.555 | 4.087 | -2.093 | K75 | KURUCZ | -1.710 | -0.383 | NIST/C+ |
| `gf_087247` | Mg I | 5183.604 | 2.717 | +0.180 | NIST ASD v5.11 grade A | NIST-C+ | -0.168 | +0.348 | Pehlivan17 |
| `gf_086153` | Hf II | 5071.200 | 1.497 | -2.665 | LNWLX | OTHER | -2.330 | -0.335 | LAWLER |
| `gf_117957` | Si I | 7083.949 | 5.984 | -1.943 | K07 | KURUCZ | -1.662 | -0.281 | NIST/E |
| `gf_117764` | Si I | 6518.733 | 5.954 | -1.982 | K07 | KURUCZ | -1.730 | -0.252 | NIST/D |
| `gf_117926` | Si I | 6976.513 | 5.954 | -1.170 | GARZ | OTHER | -0.924 | -0.246 | NIST/D |
| `gf_117783` | Si I | 6555.462 | 5.984 | -1.163 | VALD3 | VALD3 | -1.400 | +0.237 | NIST/D |
| `gf_068282` | Fe II | 4515.333 | 2.844 | -2.365 | NIST-C+ T7589 | NIST-C+ | -2.600 | +0.236 | MELBAR |
| `gf_118025` | Si I | 7272.034 | 6.099 | -1.377 | K07 | KURUCZ | -1.158 | -0.219 | NIST/D |
| `gf_117819` | Si I | 6631.048 | 5.984 | -2.668 | K07 | KURUCZ | -2.470 | -0.198 | NIST/E |
| `gf_117998` | Si I | 7210.315 | 6.083 | -2.096 | K07 | KURUCZ | -1.903 | -0.193 | NIST/E |
| `gf_000073` | Al I | 8772.865 | 4.022 | -0.170 | 1995JPhB.. | OTHER | -0.350 | +0.180 | NIST/B+ |
| `gf_087245` | Mg I | 5167.322 | 2.709 | -1.031 | NIST ASD v5.11 grade A | NIST-C+ | -0.854 | -0.177 | Pehlivan17 |
| `gf_000072` | Al I | 8076.289 | 4.087 | -2.636 | K75 | KURUCZ | -2.460 | -0.176 | NIST/C |
| `gf_000070` | Al I | 8065.968 | 4.085 | -1.936 | K75 | KURUCZ | -1.760 | -0.176 | NIST/C+ |
| `gf_069914` | Fe II | 5197.568 | 3.230 | -2.046 | NIST-C+ T7589 | NIST-C+ | -2.220 | +0.174 | MELBAR |
| `gf_000071` | Al I | 8075.353 | 4.087 | -1.681 | K75 | KURUCZ | -1.510 | -0.171 | NIST/C+ |
| `gf_117787` | Si I | 6560.566 | 5.964 | -1.409 | K07 | KURUCZ | -1.570 | +0.161 | NIST/D |

## What the cross-reference found in OUR pool — published values high by a whole isotope count

This was not what the audit was looking for, and it is the most consequential thing it found. **34 published `canonical_gf` values sit exactly log10(n_isotopes) above the correct total.**

| species | adjudicated lines | multiplier | offset | same signature, no linemake referee |
|---|---:|---:|---:|---:|
| Eu II | 4 | ×2 | +0.3005 dex | 1 |
| Nd II | 26 | ×7 | +0.8451 dex | 12 |
| Sm II | 4 | ×7 | +0.8451 dex | 0 |

The signature is self-evident and needed no isotope table to find: the published `log_gf` sits a **constant** offset above its own `gf_linelist_vald` sibling — identical to four decimals across every affected line of the species — and that constant is log10 of a small integer. What makes it a defect rather than a disagreement is the referee: **`gf_linelist_vald` is what matches linemake's primary-lab declared total, to 0.0000–0.0006 dex** (Lawler 2001 for Eu II, Den Hartog 2003 for Nd II, Lawler 2006 for Sm II). The published value is the one that is wrong, and it is wrong by exactly the number of stable isotopes the element has — 7 for Nd and Sm, 2 for Eu.

This is the same trap the audit guards against on linemake's own side: linemake lists HFS components *per isotope*, each isotope's set summing to the full gf, so a summation that ignores isotope identity inflates by log10(n_isotopes). The guard was built for the reference list and then found the defect in ours. Affected rows took the `gf_synth_ges` column rather than the VALD one; sibling rows of the same species that took the VALD column agree with linemake exactly, which rules out a species-wide scale error and localises it to a per-row column choice.

**Nothing here has been corrected** — RYA-161, validate-don't-tune. The affected rows are enumerated in `isotope_multiplier_suspects.csv` and logged in `data/audit/run_bug_ledger.csv` for adjudication.

## Where the EP key and the gf evidence point at different rows (1)

Matching on the physical transition is only as good as the excitation potential on *both* sides. These are matches where linemake's line landed on one of our rows by EP, while a **different** row of ours at the same wavelength carries linemake's log gf exactly — so the two lists disagree about the EP of the transition, and the row we matched is the wrong one. Their verdict is `EP_COLLISION`, not a promotion candidate: RYA-1034 refuses a match it cannot stand behind.

* **Fe I 5538.516 Å** — linemake gives EP 3.634 eV, log gf -1.540. Our `gf_058963` sits at that EP with log gf -5.097 (Δ -3.557), but `gf_058964` carries -1.540 — linemake's value to the digit — at EP 4.218 eV, **0.584 eV away**, sourced `PRIMARY LAB Ruffoni2014`. Wavelength and log gf say these are the same transition; the two lists disagree on its lower level by 0.58 eV. One of the two EPs is wrong and this audit cannot say which.

That there is exactly one such collision in 137,130 audited lines is itself the reassuring result — the EP key is essentially collision-free, which is what makes the wavelength+EP match trustworthy everywhere else.

## Orphaned published values (4)

`canonical_gf` keeps what each upstream delivery said for a line beside the value it publishes. These rows agree with **none** of their own provenance columns, and are not the RYA-945 laboratory upgrade (LAB-tier rows are excluded — those diverge from their stale siblings by design, and linemake reproduces all 396 matched Fe I ones to 0.000 dex, which is the positive control for this scan).

| line_id | species | λ_air (Å) | published | own siblings | linemake | source | tier |
|---|---|---:|---:|---|---:|---|---|
| `gf_087247` | Mg I | 5183.604 | +0.1800 | -0.239/-0.239/-0.239 | -0.168 (Pehlivan17) | NIST ASD v5.11 grade A | NIST-C+ |
| `gf_069752` | Fe II | 5169.028 | -0.8601 | -1.000/—/-1.250 | -1.000 (MELBAR) | NIST-C+ T7589 | NIST-C+ |
| `gf_055336` | Fe I | 4427.310 | -3.0443 | -2.924/—/-2.924 | -2.920 (OBR91) | NIST-C+ T3547 | NIST-C+ |
| `gf_070275` | Fe II | 5264.802 | -3.2306 | -3.130/-3.120/-3.120 | -3.130 (MELBAR) | NIST-C+ T7589 | NIST-C+ |

**`gf_087247`, Mg I 5183.604 (b3), is the one that does not survive scrutiny.** It publishes +0.180 tagged `NIST-C+` / "NIST ASD v5.11 grade A", while all three of its own provenance columns say −0.239 and linemake's primary-lab value (Pehlivan Rhodin et al. 2017) says −0.168 — a 0.35 dex gap. The line-to-line spacing across the Mg b triplet settles it without appeal to any external number: our two independent references agree with each other on the spacing (VALD/GES give 0.211 and 0.692 dex; linemake/Pehlivan17 give 0.195 and 0.686) and `data/linelists/nist_reference.csv` rows 29–31, the source of the published values, give **0.630 and 1.211**. Two references agreeing with each other and disagreeing with the third localises the defect to that file — which already carries a documented correction of the same class in its own header (RYA-592 fixed the grade and `aki_s-1` columns on two other Mg I rows). Both `gf_087245` (b1, 0.10 dex from its siblings) and `gf_087247` (b3) come from those rows.

The three Fe rows are 0.10–0.14 dex and are ordinary NIST-C+ versus compilation differences; they are listed for completeness, not flagged.

## Refusals, and why they are results

* **`AMBIGUOUS_MATCH` (18)** — two linemake candidates inside the tolerance, or one linemake line that is the only candidate for two distinct lines of ours. RYA-1034: a tolerance that cannot separate two candidates is a fact about the pool, not a tie to be broken by proximity. Refused, never argmin'd.
* **`HFS_AMBIGUOUS` (6)** — the HFS expansion could not be reconciled. linemake lists HFS components *per isotope*, each isotope's set summing to the full gf, so a naive component sum inflates by log10(n_isotopes) — **+0.301 dex for Eu II, +0.699 for Ba II, +0.845 for Nd II**. We compare against linemake's own declared total and verify it by per-isotope summation; only blocks that fail that verification are refused.
* **`LINEMAKE_NONPRIMARY` (51)** — linemake's own README flags the value as not a primary lab result, so a disagreement there is not evidence that we are wrong. Detected by parsing the README, not asserted:
  * `Cu I` — non-primary: README flags the source as to-be-treated-with-caution
  * `Nd II` — non-primary: solar-derived (astrophysical) log gf; scoped to 4314.50 A only
  * `CO` — non-primary: deliberate pragmatic gf offset
* **HFS asymmetry is kept, and measured.** Our `hfs_n_components > 1` rows are already collapsed physical lines (gf-weighted centroid, log10 Σgf), and linemake's comparable row is likewise always a total — so an asymmetric expansion is not by itself a refusal. It does cost agreement, and here is how much: pairs where neither side collapsed agree 79.3% of the time (median |Δ| 0.0000 dex), against 50.0% for `ours_only` (n=816, median |Δ| 0.0060), 40.6% for `linemake_only` (n=540, median |Δ| 0.0100), 37.9% for `both` (n=348, median |Δ| 0.0190). Those offsets are at the 0.006–0.02 dex level — the two lists cluster slightly different component sets. A summation error would have shown up at log10(n_isotopes), i.e. 0.30–0.85 dex, and does not.
* **`NO_MATCH` (126,563)** — linemake simply does not cover the line. That is the dominant outcome by count and it is expected: linemake is a curated few-thousand-line database, our pool is a 178,680-line survey list.

## How the numbers were derived

* **Match key** — species AND wavelength (±0.02 Å) AND excitation potential (±0.01 eV). Both tolerances were *read off the data* against a **measured** coincidence null: the same nearest-neighbour scan re-run with our wavelengths displaced by each of 20 offsets, which destroys every true pair while leaving both line densities intact. The tolerance is the outer edge of the last bin whose real count still exceeds that null by 5×. The answer does not hang on that multiple — the measured real/null ratio falls off a cliff at the selected edge (wavelength 42 → 8.7 → 2.3, EP 315 → 9.8 → 2.4), and the same rule at 3×/5×/10× gives 0.02/0.02/0.01 Å and 0.01/0.01/0.005 eV. Nothing was hardcoded. As an *independent* check the project's own ratified match window (`gf_grades.WAVE_TOL_A` = 0.02 Å, `EP_TOL_EV` = 0.02 eV) lands on the same wavelength window and a looser EP one.
* **Agreement threshold** — derived **per line**, with no free parameter: half the printed quantum of each side, summed. If we print a log gf to 3 decimals and linemake prints it to 2, the largest difference the two can show while representing one underlying number is 0.0005 + 0.005. Median over the 10,555 matched lines is 0.0100 dex (range 0.0005–0.055; lines we hold at one decimal place get a correspondingly wider threshold, because at that precision we are not entitled to call a small gap a disagreement). The matched Δ histogram corroborates it: a spike of exact zeros, a shoulder decaying to ~0.005 dex, then a flat continuum. The repo's own ratified `gf_grades.LOGGF_MATCH_TOL` (0.02 dex) is the same idea one step looser.
* **MOOG decoding** — fixed-width 4 × F10 (λ Å, species code, EP eV, log gf) then a free-text source tag from column 41. Every record is decoded twice, fixed-width and free-form, and a disagreement raises. Species code: integer part Z, first decimal digit the ionisation stage, further digits the isotope mass number.
* **Curated manifest** — the 59 + 34 file list is parsed out of linemake's own `mergenohfs` / `mergehfs` scripts, so it cannot drift from the repo it describes. The uncurated `moogatom*` bulk files are deliberately NOT used as corroboration: they carry no README source attribution, and an unsourced Kurucz value agreeing with ours is not evidence of anything.
* **Element symbols** — from `canonical_gf` (`key_z` + `ion`), cross-checked against the species stem of every linemake filename. A disagreement is fatal.

## Reproducibility of this artifact

Regenerated on Sirius (CI venv, py3.12) and compared against the Mac run: **every verdict count and every reported number is identical**. The files are nonetheless not byte-identical, and it is worth saying why rather than claiming they are. pandas' CSV float parser differs by up to one ULP between versions, so pass-through columns can print differently — `our_EP` as `4.8271999999999995` on one machine and `4.8272` on the other, `our_wavelength` as `19280.14200091231` against `19280.142000912318`.

One informational count moves with it: `bulk_moogatom.beyond_curated` is 57,163 on one machine and 57,164 on the other — a single line of 66,008 sitting exactly on the bulk-coverage tolerance boundary. No verdict, no match, and no scientific number is affected. That boundary case is deliberately **not** papered over by rounding the comparison inputs: a rounded number is not an identity, and manufacturing agreement at a boundary would be worse than recording that the boundary is there.

Row ordering *is* pinned. Every count- or magnitude-ordered table breaks ties on a stable key, because the first cross-machine comparison found the species tied at 3 and 6 matches emitted in a different order — same numbers, different rows. The clone path is deliberately not recorded in `provenance.json` for the same reason: it is a scratch temp directory that differs per machine and carries no information.

## Scope

The uncurated `moogatom*` bulk lists **are** parsed — 215,514 atomic lines over 132 species — and are then deliberately excluded from every verdict. That exclusion is a measured decision, not an oversight: 66,008 of our audited lines have a bulk counterpart, **57,163 of them lines the curated database does not reach at all**, so including them would have multiplied apparent coverage six-fold. But those files carry no source tag and no README attribution. An unsourced value agreeing with our unsourced value corroborates nothing, and one disagreeing is not a finding. Coverage is not the question this audit asks; provenance is.

Molecular species are out of scope for the numeric cross-reference: step 2 of the ticket scopes this to the atomic MOOG lists, and our pool holds no CO at all. The CO Δv = 2 caveat is recorded from the README above for completeness. `linemake` lines with no counterpart in our pool are **counted only**, per species, in `species_summary.csv` — they are not enumerated as import candidates, which is a separate decision.

