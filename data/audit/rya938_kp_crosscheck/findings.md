# RYA-938 — Kitt Peak 1984 vs Kurucz 2005: what the cross-check found

Three findings, in order of consequence. Two are defects; the third is the
answer RYA-929 was after.

## 1. Kurucz 2005 is on a VACUUM grid. RYA-929 read it as air.

RYA-929 produced two results that could not both be right: a clean window-level
GO, and a line-level table where Kurucz showed almost no absorption at lines the
1984 atlas and the telluric-free IAG atlas both call deep.

One cause explains both. Reading a vacuum grid as air displaces every line by the
air–vacuum offset — **~1.7 Å at 6600 Å, about 200 sampled pixels**. That empties a
±0.35 Å line window while barely moving a 17–300 Å window statistic.

Measured by cross-correlating each product onto the 1984 air grid:

| window (Å) | IAG lag | Kurucz **as read** | Kurucz **vac→air** | predicted air–vac offset |
|---|---:|---:|---:|---:|
| 5160–5200 | −0.000 | +1.445 | −0.000 | +1.443 |
| 6100–6160 | −0.000 | +1.700 | −0.000 | +1.697 |
| 6600–6650 | −0.000 | +1.835 | −0.000 | +1.830 |

The measured lag matches the independently predicted offset to **2–5 mÅ**, and it
**scales with wavelength** exactly as the refractive index does. A constant
registration error, an RV shift, or a grid bug would not reproduce that slope.

**The IAG column is the positive control.** It is a genuinely-air product and
returns −0.000 Å, so a zero from this method means something.

The Kurucz header states only `nm` and never states the medium — which is why
RYA-929 could not settle it from documentation and correctly deferred rather than
guessing.

### The lines RYA-929 flagged, resolved

Depth = 1 − min(flux) in ±0.35 Å.

| line | 1984 | IAG | Kurucz as read (RYA-929) | Kurucz vac→air |
|---|---:|---:|---:|---:|
| Al I 6696.185 | 0.2544 | 0.2572 | 0.0081 | **0.2483** |
| K I 7698.964 | 0.8078 | 0.8074 | 0.0394 | **0.8012** |
| Fe I 6881.442 | 0.2137 | 0.2016 | 0.0542 | **0.2003** |
| N I 8216.336 | 0.0424 | 0.0351 | 0.1142 | **0.0362** |
| Al I 6631.218 | 0.0269 | 0.0311 | 0.0142 | 0.0248 |

All four flagged discrepancies are **artifacts of the medium error**, not defects
in the Kurucz product.

Two notes on the remaining entries:

* **K I 7664.899** — 1984 reads 1.0015 (the O₂ A band saturates on top of it),
  IAG 0.6425, Kurucz vac→air 0.8338. No longer absurd, but IAG and Kurucz still
  disagree by 0.19 on a strong resonance line inside a saturated telluric band.
  That is a real difference between two correction methods and stays **OPEN**.
* **N I 8216** — RYA-929 also used 8216.000 as the centre; the line is at
  8216.336. Its "discrepancy" was part medium error, part wrong centre.

## 2. `lm0840` was a saved HTTP 500 page, and the loader called it missing coverage

The 1984 atlas segment covering **8400–8441 Å** was 714 bytes of archived HTML
from `nispdata.nso.edu` — in **both** staged copies, so a failed download was
propagated, not corrupted once. The bundled `README` was a second saved 500 page.

`measure_band_ew.kp_segments()` swallowed the parse failure in a bare
`except Exception: continue`. Measured: **251 files on disk, 250 inventoried**, and
`load_kp_window(segs, 8420.0, 1.0)` raised

```
LookupError: no Kitt Peak segment covers 8420.000 A
```

**A corrupt file presented as missing coverage** — the RYA-833 shape, and exactly
the failure `_resolve_kp_dir` was written to prevent one level up.

**Blast radius: no committed result.** No Kitt Peak product covers 8400–8440 Å;
the reddest red-optical product, `AlI_6910_9199`, holds six lines
(7361.6, 7362.3, 7835.3, 7836.1, 8772.9, 8773.9) and none in the gap.

**Latent exposure was large.** The solar line list holds **33 Fe candidates
(30 Fe I + 3 Fe II)**, 61 Co I, 55 CN I, 14 Mn II and 12 Ti I in 8400–8440 Å. The
planned Fe/Fe II Kitt Peak runs would have reported all of them as uncovered.

Repaired: re-fetched from `https://nispdata.nso.edu/ftp/pub/atlas/fluxatl/lm0840`
(154,546 bytes, 4,067 rows, 8400.004–8440.999 Å, sha256 `7823350289012fa4…`),
installed in both copies, corrupt originals quarantined. The real README came with
it, and it is what settles §3.

The guard now lives in `pipeline/kp_atlas_integrity.py` and pins what a segment
**is** — three numeric columns, strictly increasing wavelength, residual flux near
unity, stem consistent with content — so it keeps working now that this particular
file is fixed. `kp_segments(allow_corrupt=True)` remains available for deliberate
degraded operation, because a caller who *chooses* to skip data is different from
one who never knew.

### Full structural pass, after repair

251 files, **251 good, 0 corrupt, 0 coverage gaps**; every segment monotonic and
stem-consistent; span 40.01–41.00 Å on a 4 nm stride (the README's "4 nm plus
0.1 nm of overlap"); 1,166,056 rows spanning **2960.00–13000.02 Å**.

## 3. Medium, continuum and units — settled from the source

The recovered NSO Atlas No. 1 README states the 1984 columns outright:

> "The first column contains the wavelength **in air**, the second column contains
> the **pseudo-residual flux**, and the third column contains the corresponding
> observed **irradiance in units of micro-watts per square centimeter per
> nanometer**."

| | 1984 Kitt Peak | Kurucz 2005 | IAG/Baker 2020 |
|---|---|---|---|
| medium | **air** (documented) | **vacuum** (measured) | vacuum wavenumber → air |
| flux | pseudo-residual, unity IS the continuum | absolute irradiance W/m²/nm | normalised, telluric-free |
| continuum | ships one — nothing re-fitted | **ships none** | ships one |
| coverage | 2960–13000 Å | 3000–10000 Å | 5000–10000 Å |

Only Kurucz needs a continuum imposed, and it gets a **diagnostic** 8 Å running
envelope, stated as such. The 1984 and IAG products are consumed as delivered —
no local re-fit, per RYA-911, where re-fitting put the continuum *below* the
observed flux and cost a median 23.8% of the EW.

## 4. Telluric state, on the RYA-805/905 instrument

`pct_below_0.5` / min flux, comparable to the existing HARPS/IAG/Elgueta figures.

| band (Å) | 1984 | Kurucz 2005 | IAG |
|---|---:|---:|---:|
| O₂ B 6867–6884 | 15.80 / −0.006 | 0.00 / 0.764 | 0.00 / 0.760 |
| H₂O 7160–7340 | 2.77 / 0.013 | 0.24 / 0.378 | 0.18 / 0.369 |
| O₂ A 7594–7685 | 32.83 / −0.007 | 0.28 / 0.166 | 0.18 / 0.357 |
| H₂O 8100–8400 | 3.58 / 0.001 | 0.29 / 0.193 | 0.29 / 0.192 |
| H₂O 9280–9600 | 29.12 / −0.002 | 0.02 / 0.469 | 0.00 / 0.547 |
| H₂O 11120–11560 | 33.86 / −0.004 | *not covered* | *not covered* |
| clean 6600–6650 | 0.14 / 0.419 | 0.14 / 0.425 | 0.15 / 0.415 |
| clean 6690–6702 | 0.00 / 0.746 | 0.00 / 0.752 | 0.00 / 0.743 |
| clean 5160–5200 | 12.03 / 0.068 | 11.62 / 0.071 | 12.16 / 0.065 |

RYA-929's window-level conclusion **stands**: 1984 retains substantial terrestrial
absorption; Kurucz 2005 behaves like the telluric-free reference.

The clean controls agree across all three to <1 percentage point — including
5160–5200 Å, where all three read ~12% below half continuum because the Mg b
lines are genuinely that deep. A control that moved would have invalidated the
band comparison.

**1984 min flux goes negative** (−0.002 to −0.007) in saturated telluric cores.
That is a real property of the FTS pseudo-residual flux, not a defect.

## 5. What this does and does not license

- **Does**: it establishes that Kurucz 2005 is usable at line level *once read in
  vacuum*, and that the two products can be compared honestly.
- **Does not**: substitute any Kurucz line into an abundance product. That is
  RYA-932/933/934, and the K I 7665 disagreement above is a reason to keep
  per-line dispositions rather than accepting the arm wholesale.
- **1984 remains the only arm beyond 10000 Å**, and it is uncorrected there
  (33.86% below half continuum in H₂O 11120–11560). A corrected Kitt Peak IR arm
  cannot come from Kurucz 2005.

Artifacts: `crosscheck.json`, `registration.csv`, `telluric_bands.csv`,
`diagnostic_lines.csv`, `atlas_integrity.json`, `rya938_registration.png`.
Reproduce with `scripts/rya938_kp_crosscheck.py` then `scripts/rya938_kp_plot.py`.
