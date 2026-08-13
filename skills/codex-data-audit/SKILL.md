---
name: codex-data-audit
description: The intake discipline for a held dataset — what must be determined from the FILES before any measurement runs on them, and what must be RECORDED so the next session does not re-derive it. Use this skill whenever data is acquired, staged, re-pulled, or first measured; whenever a ticket says DATA AUDIT; and whenever asked "can we measure X on this data". Owns the mandatory per-holding telluric_applied determination (RYA-806) that gates the RYA-424 telluric stage, plus the product-level, frame, coverage and SNR checks that decide GO / CAVEAT / NO-GO. A dataset that has not been through this is not measurable, and "it looks reduced" is not a determination.
---

# Codex data audit — intake skill

## Why this exists

RYA-370 asserted "no telluric correction" for the CRIRES+ Vesta set. RYA-373's spec
repeated it. RYA-794 and RYA-796 built on it. **Four tickets inherited a claim that nobody
had shown the keyword for.** It happened to be true (RYA-805 confirmed it), but that was
not knowable without looking, and the cost of looking was one afternoon against four
tickets of exposure.

The rule this skill enforces: **a property of the data is determined from the data, once,
and recorded** — never inherited from a sibling ticket, never inferred from the
instrument, never assumed from "it came from a pipeline so it must be reduced."

## The six steps

Run all six. Step 1 and Step 2 produce values that get WRITTEN, not just reported.

### Step 1 — Identity and product level
`OBJECT`, `INSTRUME`, `DATE-OBS`, `PRO CATG` / `PRODCATG`, `PROCSOFT` / `PIPE ID`.
Establish **what the product IS**, not what the folder is called. Vintage comes from the
pipeline id (`cr2res` is CRIRES+, not the 2012-13 CRIRES). Never trust the directory name
or the program id in a sibling ticket — RYA-794 found `PROG ID` was `60.A-9051(A)`, not
the `1102.22KH` a prior ticket recorded.

### Step 2 — Telluric state — MANDATORY for every IR dataset (RYA-806)
**Determine `telluric_applied` ∈ {`applied`, `not-applied`, `unknown`} from the headers,
and WRITE it to `data/catalog/holdings_manifest_registry.csv`.** Not a Step-6 glyph — a
recorded value, because it is a software switch the loader reads.

    from pipeline.telluric_intake import from_headers, from_many
    value, evidence = from_many(files)     # -> ('not-applied', [TelluricEvidence, ...])
    print(evidence[0].citation())          # paste into the holdings `notes`

It checks three places, because a pipeline records a telluric step in any of them:

1. **The recipe chain** — `ESO PRO REC*n* ID` walked to full depth. A molecfit / telluric
   / `corr_tell` / `calctrans` recipe is direct evidence of APPLIED. A chain that ends at
   extraction is evidence of NOT-APPLIED.
2. **A transmission extension** — `TRANS` / `RECON` / `TELL` / `MTRANS`. ⚠️ Presence is
   **not** proof: **a transmission array that is all 1.0 means the correction was NOT
   applied** — the model was computed and nothing was divided by it.
3. **A second flux column** — `FLUX_TELL*` beside `FLUX`.

⚠️ **`applied` can be true of the PRODUCT and false of its DEFAULT COLUMN.** NIRPS
`S1D_FINAL_A` carries `FLUX_TELL_EL` *and* an uncorrected `FLUX` / `FLUX_EL` / `FLUX_CAL`
in the same file. Record which column must be read; a consumer reading `FLUX` out of a
holding labelled `applied` gets uncorrected flux. `TelluricEvidence.required_column` says
which.

⚠️ **`unknown` is a real answer and it is never defaulted.** HARPS Phase-3 ADPs carry no
`PRO REC` chain at all, so their headers genuinely do not speak to telluric state.
Assuming `applied` fabricates a correction (forbidden, RYA-786); assuming `not-applied`
risks correcting an already-corrected product twice. Record `unknown` and let
`telluric_policy.gate_holding` refuse.

**Do not conflate the two axes.** `telluric_basis` (per-INSTRUMENT: does this band NEED
correction) and `telluric_applied` (per-HOLDING: has THIS product got one) are orthogonal.
Measured proof they are: alpha Cen CRIRES+ and alpha Cen NIRPS are both
`telluric_required=yes`, and one product is corrected while the other is not.

### Step 3 — Wavelength frame
`SPECSYS`, and the velocity keyword that backs it. ⚠️ A frame keyword is a **claim**, not
a measurement: `ESO TEL TARG RADVEL = 0.0` in all 18 Vesta IDPs is a placeholder, not a
zero velocity (RYA-796). ⚠️ For a **reflected** source the shift is the two-leg
Sun→body→observer rate, not a BERV (RYA-372). ⚠️ `astropy` strips `HIERARCH `, so a
lookup of `HIERARCH ESO ...` returns empty and **manufactures an absence** (RYA-791).

### Step 4 — Coverage
Test against **real pixels** — `QUAL == 0`, finite non-zero flux — never `WAVELMIN` /
`WAVELMAX`, which hide detector and inter-order gaps (RYA-377). Coverage on a tiled
instrument is a comb of settings, not a span.

### Step 5 — SNR and normalisation
Per-frame SNR against the floor; `FLUXCAL`; `CONTNORM`. ⚠️ Whether the product is
continuum-normalised is a property of the product, not a preference — treating
un-normalised adu as normalised is the RYA-713 defect (EWs low by a median 11.7 %, worst
71.4 %).

### Step 6 — The verdict table, and what gets written
Post as a Linear comment, each axis **GO / CAVEAT / NO-GO**:

| Axis | Finding | Verdict |
|---|---|---|
| Product type | `PRO CATG`, pipeline + version | |
| Telluric correction | the Step-2 value **and the keyword evidence** | |
| SNR | range, median, vs floor | |
| Wavelength | span, gaps, units | |
| BERV / rest frame | `SPECSYS` + the velocity keyword's credibility | |

**And WRITE, in the same PR:**
- `telluric_applied` + its evidence on the holding's row in
  `holdings_manifest_registry.csv` — this is the switch `gate_holding` reads.
- The holding itself, if new, on that registry (it is the anti-reinvent surface —
  **read it before proposing any download**).
- A register + `SEQUENCE.md` bump: the holdings registry is a **state surface**, so
  touching it without bumping fails CI's RYA-659 guard. `check_register_freshness.py
  --since-main origin/main` reads **committed** state — run it after committing.

## The rules (do not break)

- **An absence needs a positive control.** "I found no telluric keyword" is
  indistinguishable from a broken lookup. Pair it with something that WOULD carry the
  marker — our own corrected products carry `mtrans`; the Vesta IDPs carry none.
- **Enumerate the false positives.** A naive telluric regex over the Vesta headers returns
  **162 hits, all 162 spurious** (`ESO DET DEV1 BOARD*n* TRANS` shift registers, the
  `ESO OBS AMBI TRANS` sky-transparency *constraint*, FITS boilerplate). Filter them in
  code, with a comment, and test the filter.
- **Verify the test discriminates before trusting it.** A depth comparison in a window
  where both spectra are clean "confirms" whichever answer you expected. Run the control
  through the same test and check it SEPARATES (RYA-805: a 980–1080 nm comparison agreed
  to 0.05 pp and proved nothing; the O₂ 1.27 µm band separated 4–5×).
- **Prefer a falsifiable covariate.** Correlating absorption depth against the header's
  own water vapour (r = +0.996) forecloses "it might be stellar" instead of merely being
  consistent with "it's telluric".
- **Assert on CONTENT, not presence.** An unmarked directory link once saved a spectrum
  directory as an HTML page with every guard green (RYA-789). Check `SIMPLE` on a FITS,
  check the array, check the units.
- **Units are per-file and hostile.** CRIRES+ IDP `WAVE` is nm; Elgueta `sp/` is labelled
  `0.1nm` and is nm; `atomicy.dat` carries the same label and really is Å (RYA-794).
- **Never re-download before reading the registry.** RYA-805 established that no
  telluric-corrected CRIRES+ product exists anywhere in ESO's 4131-file collection —
  that answer is recorded so nobody pulls 37 MB again to re-learn it.

## Smoke test

    python3 -c "from pipeline.telluric_intake import from_headers; \
                print(from_headers('<one product file>').citation())"
    python3 scripts/rya806_backfill_telluric_applied.py     # report; --write to record

Then confirm the switch is live: a `load_window` call on an IR holding that is
`not-applied` raises `TelluricNotCorrected`, and one that is `unknown` raises
`TelluricStateUnknown`.
