# `asplund2021_al` — the AGSS21-lineage solar Al reference line set

**RYA-1173.** The Al answer to RYA-946's *Mandatory Solar reference-line-set census — AGSS21
lineage, all bands*, built the way RYA-1109 built `asplund2021_fe` — except that the Fe route does
not exist here.

## 🔴 AGSS21 publishes no Al line list

RYA-1109 could transcribe AGSS21's **own** Table A.2 for Fe. There is no such table for Al. AGSS21's
*Aluminium (Z = 13)* section adopts someone else's number, so the reference set has to be
reconstructed from the papers AGSS21 points at:

```
  Asplund, Amarsi & Grevesse 2021, A&A 653, A141      adopts A(Al) = 6.43 ± 0.03. Publishes
    10.1051/0004-6361/202140445                       NOTHING per line.
      │
      ├─► Nordlander & Lind 2017, A&A 607, A75        the analysis AGSS21 adopts.
      │     10.1051/0004-6361/201730427               Table A.1 = per-line level identity, E_low,
      │                                               log gf, σ, van der Waals.
      │                                               §3.1.5 = the telluric exclusion.
      │                                               Fig. 8 = the six used lines, named.
      │
      └─────► Scott, Grevesse, Asplund et al. 2015b   the SELECTION: seven Al i lines.
                A&A 573, A25                          Table 2 = EWs, weights, five-model
                10.1051/0004-6361/201424109           abundances. Table 3 = level J identity.
```

`SOURCE_LINE_LIST_NOT_PUBLISHED` does **not** apply. The primaries publish the rows; they are here.

⚠️ **`Scott et al. (2015b)` is A&A 573, A25, not A26.** AGSS21's reference list orders the two 2015
Scott papers by author list, so `2015a` is the iron-group paper and `2015b` is *The intermediate mass
elements Na to Ca* — the one containing Al. Reading the letters as volume order sends this census to
the wrong paper.

## The set is six used + one excluded — and AGSS21 says five

| source | count | what it is |
|---|---|---|
| Scott et al. 2015b | **7** | "We retained seven quite weak Al i lines (Table 2)" |
| Nordlander & Lind 2017 | **6** | §3.1.5 drops 10891 Å for telluric contamination; Fig. 8 names the six |
| AGSS21 prose | **5** | "for the Sun, departures from LTE are in fact largely unimportant for these five Al i lines" |

7 − 1 = 6, and the primary's own figure axis names all six. **AGSS21's "five" reproduces from
neither source it cites.** This set carries six on the authority of the primary and records the
conflict in `published_line_count_conflict`; it is not resolved quietly. Two candidate explanations
are written down in `raw/lineage_quotations.md` and *neither is adopted*, because no paper states
either. A replication measuring this set measures **six** lines.

`10891.732 Å` ships **in this file**, flagged `EXCLUDED_BY_SOURCE_ANALYSIS` with its published
reason, along with Scott's EW, weight and per-model abundances for it. RYA-946: preserve explicit
negative selections. Dropping it would make this set indistinguishable from one that never had it.

## 🔴 Every gf in this set is theory

Scott §5.3: *"The data for our adopted lines come from theoretical calculations by the OP (Mendoza
et al. 1995), under the assumption of LS-coupling."* NL2017's per-line reference code agrees on all
seven (code 2 = TOPbase/Mendoza). So the AGSS21-lineage Al set is **not** a laboratory-gf set, and a
Codex line matching into it has gained no gf grade — RYA-946: *"an AGSS21 abundance value is not a
gf grade."*

## Files

| file | what |
|---|---|
| `asplund2021_al_lines.csv` | **the reference set**: 7 rows, 6 `USED_BY_SOURCE_ANALYSIS` + 1 `EXCLUDED_BY_SOURCE_ANALYSIS` |
| `nordlander_lind_2017_analysis_lines.csv` | all **55** rows of NL2017 Table A.1, with a `role` column |
| `asplund2021_al_lines.prov.json` | provenance + the five extraction controls, with their measured numbers |
| `raw/nordlander_lind_2017_tableA1.tsv` | verbatim Table A.1 |
| `raw/scott2015b_table2_al.tsv` | verbatim Table 2, Al block (EW, weight, per-model abundances) |
| `raw/scott2015b_table3_al_hfs.tsv` | verbatim Table 3, Al block (level J + HFS constants) |
| `raw/lineage_quotations.md` | the chain in the sources' own words, and the count conflict |

### The 55 analysis lines are not the abundance set

NL2017 analyses far more than six transitions — the 3944/3961 Å resonance lines, centre-to-limb
variation in 7835 Å, the 12.33 µm emission line, the HFS-sensitive 13123/16750 Å IR lines, HST/STIS
UV lines in metal-poor stars, J-band lines in HD 122563. Table A.1's own note says only that it
lists *"lines used in the spectrum analyses"* — a statement about the table, **not about any row**.
So `role` is `SOLAR_DIAGNOSTIC_NAMED_IN_TEXT` only where the running text names that line, and
`ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE` for the other 43. Inferring a role from a wavelength's
familiarity is what RYA-946 forbids.

## Extraction controls

Rebuild and re-run them with `python3 scripts/rya1173_build_asplund_al_lineset.py`; the artifacts are
**not written** if any control fails. `--check` reproduces them and compares to what is committed.

| | test | result |
|---|---|---|
| **C1** | Scott's abundance summary table (published Table 1 = preprint Table 5) — five model columns, two differences, the recommended value — must reproduce as the `Wt`-weighted mean of (model + ΔNLTE) | **8/8 exact.** The LTE misreading reproduces **0/5** |
| **C2** | `a_nlte_3d − a_lte_3d == delta_nlte_3d` per line, a relation Scott never writes down | max residual **0.0 dex** |
| **C3** | NL2017 Table A.1 and Scott Table 2, joined on the **level**, must agree exactly on log gf and E_low | **7/7 exact** |
| **C4** | E_up = E_low + hc/λ_vac single-valued per upper level, across all 55 rows — an identity neither paper tabulates | **0.00037 eV** vs a derived 0.001 eV tolerance |
| **C5** | The derived six are exactly NL2017's Fig. 8 axis labels, and those labels are unambiguous within Table A.1 | **PASS** |

Three of these carry a **measured negative**, because a control that cannot fail is not a control:

* **C1** — reading the summary table as LTE columns (its caption says *"A summary of the NLTE
  results"*) puts all five exactly 0.01 dex low **while both difference columns still reproduce**.
  That is a wrong reading that looks nearly right, and it is why the caption decides, not the
  arithmetic. Recorded as `lte_reading_measured_negative`.
* **C1** — the weight direction is measured, not assumed: only `larger = better` gives the published
  6.43 (unweighted 6.4187, inverted 6.4126).
* **C4** — the wavelength **medium** is measured, not assumed. NL2017 never says, and its table runs
  2103 Å – 12.33 µm across the convention boundary. Read as air the levels close to 0.00037 eV; read
  as vacuum, 0.00167 eV — **4.5× worse, past tolerance**. The column is air, including its UV rows.

**C4 also reports where it is blind.** Its sensitivity is hc/λ², so at 21208 Å a wavelength would
have to be wrong by ~23 Å before the identity noticed, against ~0.03 Å at 2204 Å. `detection_floor_A`
carries that per level so a global PASS is not read as uniform coverage.

## Not the same as `asplund2021_fe`

That directory is AGSS21's **own** published table. This one is a reconstruction from AGSS21's cited
primaries. Their `line_set` axis values are `asplund` and `asplund-al` and must never be merged: they
are different provenance chains, and RYA-1127 put `line_set` in the product identity key.
