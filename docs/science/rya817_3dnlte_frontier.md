# RYA-817 — the Amarsi 2022 3D-NLTE Fe engine, reactivated and domain-checked

**Ticket:** RYA-817 (frontier, post-calibration). **Star:** Sun (Teff 5772 K, log g 4.438,
vmic 1.00 km/s), Kitt Peak flux atlas, RYA-783 profile-fit EWs.
**Standing rule set here:** *"we use every tool we can"* — recorded in
`docs/SCIENCE_STANDARDS.md`.

Reproduce:

```
python3 scripts/rya817_recover_amarsi_training_set.py --show-discrimination
python3 scripts/rya817_run_3dnlte_bands.py --band-products-dir <RYA-783 band_products dir>
python3 scripts/rya817_check_3dnlte_availability.py
```

---

## The headline

**The engine is real, it works, and it does not reach the near-IR.**

Reactivated, the Amarsi+2022 3D-NLTE MLP reproduces its own paper's published solar row to
≤ 0.005 dex on both ions. Applied to the solar near-IR Fe band, **0 of 94 classifiable Fe I
lines and 0 of 4 Fe II lines are inside its training domain**, and every one of them is rejected on the
same axis: the transition energy.

The ticket's premise was that the MLP is wavelength-agnostic because it splits on
excitation potential rather than wavelength. That is true of its **routing** and false of
its **inputs**. The network takes seven features, two of which are E_low and E_up — and
E_up is defined as E_low + hc/λ_vac (verified to 4×10⁻⁶ eV against the vendored
`test_data.csv`). So **(E_up − E_low) IS the wavelength**, entering as a derived feature.
The training set spans 1.820–2.589 eV of transition energy; the measured IR band spans
1.355–1.792 eV. The two do not overlap at all.

A line can sit comfortably inside the trained E_low range, inside the trained E_up range,
inside the trained log gf range — and still be an extrapolation, because it is the PAIR
that encodes a photon energy the network never saw. The right-hand panel of
`data/results/rya817/rya817_domain_map.png` is that statement in one picture.

---

## What the training domain actually is, and how we know

`vendor/1L-3NErrors` ships three trained `MLPRegressor`s and nothing else — no training
data, no line list, no feature bounds. Its own driver checks only the stellar box
(Teff 5000–6500 K, log g 4.0–4.5, vmic 0–3 km/s, A(Fe;3N) 4.5–7.5) and will return a
confident float for any line parameters you hand it. `pipeline/nlte_corrections.py`
inherited exactly that blind spot.

The training line list was therefore **recovered, not assumed**, through the paper chain:

- Amarsi, Liljegren & Nissen 2022 (A&A 668, A68) Sect. 4 names it — the "golden" list of
  Jofré et al. 2014 (A&A 564, A133), Tables 4 and 5 — and quotes its wavelength range,
  478.783–681.026 nm.
- Those tables ship as `golden_Fe1.tex` / `golden_Fe2.tex` in the Jofré arXiv source
  (arXiv:1309.1099). VizieR J/A+A/564/A133 carries only table6 (which lines each star
  used), and there is **no VizieR catalogue for J/A+A/668/A68 at all**.

Four controls, all asserted by `scripts/rya817_recover_amarsi_training_set.py`, and the
script refuses to write the artifact if any fails:

| control | result |
|---|---|
| count | 171 Fe I + 12 Fe II — the counts Amarsi 2022 states, exactly |
| wavelength | 4787.83–6810.26 Å — the 478.783–681.026 nm the paper quotes, to the digit |
| E_up convention | E_low + hc/λ_vac reproduces the vendored `test_data.csv` E_up column on 13 rows to 4.1×10⁻⁶ eV |
| scaler moments | the vendored `StandardScaler`s carry `mean_`/`scale_` **of the training data**; the recovered subsets reproduce them (gt02 E_low 3.7135 vs 3.7155; fe2 log gf −3.0514 vs −3.0614) |

The scaler control is shown to **discriminate** (`--show-discrimination`, the RYA-805
rule): the same list with E_low shifted +1 eV, and one with randomised log gf, are both
rejected.

### The recovered domain

| network | n | λ_air (Å) | E_low (eV) | E_up (eV) | log gf | E_up − E_low (eV) |
|---|---|---|---|---|---|---|
| `lt02` (Fe I, E_low < 2) | 17 | 4994.1–6739.5 | 0.000–1.557 | 2.449–3.984 | −5.918…−1.988 | 1.8392–2.4819 |
| `gt02` (Fe I, E_low ≥ 2) | 154 | 4787.8–6810.3 | 2.176–5.099 | 4.103–7.293 | −4.038…+0.643 | 1.8200–2.5888 |
| `fe2` (Fe II) | 12 | 4923.9–6456.4 | 2.807–3.903 | 4.818–5.876 | −3.881…−1.260 | 1.9198–2.5173 |

Artifact: `data/reference/amarsi2022_training/amarsi2022_training_lines.csv` (+ `.prov.json`).

---

## The reactivation control — a published result, reproduced

Amarsi+2022 Table 6, solar row, reproduced from **Amarsi's own inputs**: the per-line 1D
LTE abundances of Allende Prieto et al. 2002 (ApJ 567, 544) Table 2, restricted to weak
lines with REW < −4.9 exactly as Sect. 6.1 specifies.

| | computed | published | n |
|---|---|---|---|
| A(Fe I) 1D LTE | 7.4671 | 7.47 | 41 |
| A(Fe I) 3D non-LTE | 7.4646 | 7.46 | 41 |
| A(Fe II) 1D LTE | 7.4054 | 7.41 | 13 |
| A(Fe II) 3D non-LTE | 7.4710 | 7.47 | 13 |
| Fe I d(Δ₁ₗ₃ₙ) over E_low 0→5 eV at [M/H]=0 | −0.0905 | "−0.05 to −0.1" | 171 |

**Which line list is the trap.** Run the same control over the *training* list and Fe I
misses by ≈ 0.04 dex — because the golden set is dominated by high-E_low lines while
Amarsi's solar set reaches down to 0.05 eV, and the correction tracks E_low with r = +0.94.
That miss looks exactly like a defect and is not one. It is the RYA-785 wrong-referee
failure, and it is why the control names its line list in the output.

### A wrong number, corrected at the source

`pipeline/nlte_corrections.py` asserted in a docstring that *"for the Sun: A(3D) = 7.58 +
(−0.127) ≈ 7.45"*. **The network says no such thing.** The solar Fe I correction is
**−0.002 dex**; the Fe II correction is **+0.066 dex**. The −0.127 is the gap between our
own 1D-LTE 7.58 and Asplund's 7.45 — a number reverse-engineered from a target rather than
computed from the model. The docstring is corrected in this change.

**The consequence is worth stating plainly:** the Amarsi 3D-NLTE correction is *incapable*
of carrying a 1D-LTE 7.58 to 7.466. It never was. Our VIS 1D-LTE excess is a gf zero-point
matter (RYA-161's ungraded-K07 0.17 dex floor), not a 3D-NLTE one, and the gold anchor's
7.466 comes from RYA-553's Magic-2013 1D→3D offset on a different base — not from this
engine. Anyone reading the old docstring would have concluded otherwise.

---

## Deliverable A — the products

`ENGINE-A-3DNLTE` is a new member of `pipeline.band_products.TREATMENTS`: the same EW
route as ENGINE-A, but the per-line departure comes from the Amarsi 3D-NLTE MLP instead of
the MPIA/Bergemann 1D-NLTE grid. Different **dimensionality**, so RYA-712 makes it its own
product; it is never merged with, or relabelled as, ENGINE-A.

| product | in-domain | measured | A(Fe) | σ (line scatter) |
|---|---|---|---|---|
| Fe I VIS `ENGINE-A-3DNLTE` | 114 | 159 | **7.604** | 0.342 |
| Fe II VIS `ENGINE-A-3DNLTE` | 7 | 11 | **7.642** | 0.205 |
| Fe I IR | **0** | 103 | — **NOT PRODUCED** | — |
| Fe II IR | **0** | 4 | — **NOT PRODUCED** | — |

Median per-line correction: Fe I VIS **+0.038 dex**, Fe II VIS **+0.054 dex**. Compare the
RYA-783 1D-LTE cells: Fe I VIS 7.586, Fe II VIS 7.568.

### Out-of-domain, by axis

Counts are of the lines that could be classified at all (both a 1D-LTE abundance and a
matched VALD transition); a line can fail more than one axis.

| band | classified | transition energy | feature box | level representation | stellar box |
|---|---|---|---|---|---|
| Fe I VIS | 146 | 14 | 6 | 17 | 0 |
| Fe II VIS | 11 | 0 | 4 | 2 | 0 |
| **Fe I IR** | **94** | **94** | 11 | 25 | 0 |
| **Fe II IR** | **4** | **4** | 0 | 0 | 0 |

The VIS rejections are the parts of the 3800–6910 Å band that fall outside the training
list's own 4788–6810 Å span, which is the check behaving correctly on data it should
partly accept — the positive half of the discrimination test. **The IR rejection is total
and it is one axis.**

### What the refused extrapolation would have been

Recorded as a diagnostic, in a column no aggregate reads: the network, forced past its
domain, returns a median **+0.026 dex** for the IR lines (range −0.058…+0.143). That
number is *not small enough to be harmless and not trustworthy enough to publish* — it is
the same order as the VIS correction we do trust, which is precisely why quoting it would
be indistinguishable from a real result.

### The A(Fe;3N) axis, and one trap worth recording

The grid's fourth axis is the **3D non-LTE** iron abundance, ceiling 7.5. Our solar
per-line 1D-LTE abundances run to ~7.8 (the RYA-161 gf systematic), so:

- Testing the **initial guess** against that ceiling rejects every strong line and keeps
  the weak ones. The first run of this ticket did that and the VIS control came out
  **0.07 dex low by selection alone**. The axis test belongs on the **converged** A(Fe;3N).
- The self-consistent solar axis lands at 7.604 (Fe I) and **rails** on the 7.5 ceiling.
  That is a property of our gf zero point, not of the Sun (whose A(Fe) ≈ 7.46 sits
  comfortably inside). The clamp is recorded, never silent, and its cost is measured: the
  product moves **0.007 dex across the entire 5.0–7.5 axis**, so it is immaterial.

The archived per-line-axis convention (RYA-207: each line starts from its own 1D-LTE
value) is carried as a diagnostic column. On the same in-domain lines it reads 7.386 for
Fe I VIS against the star-axis 7.604 — a **0.22 dex** difference that comes entirely from
which quantity is placed on a stellar axis. The vendor README describes the star-axis
procedure, so that is the product; the difference is on the record rather than buried.

---

## The comparison the ticket asked for

> *"IR-3D-NLTE vs the VIS-3D-NLTE gold anchor (7.466) vs Asplund 7.46. Asplund 2021
> predicts NO wavelength dependence in proper 3D-NLTE."*

**This test cannot be run with this engine, and that is the result.** One side of the
comparison does not exist: there is no in-domain IR line to build it from. Producing an IR
number here would require extrapolating a neural network into a region of transition energy
containing zero training points — inventing a number the training never saw, which is the
one thing the ticket forbids.

What *can* be said: the engine's solar correction is ≈ 0.00 dex for Fe I and +0.07 for
Fe II, so on the VIS band it moves 1D-LTE 7.586 → 7.604, nowhere near 7.466. The
band-dependence question in RYA-783 (IR − VIS = +0.053 in 1D-LTE, −0.008 on the flux-fit
route) is **untouched** by this ticket and stays where RYA-783 left it: an EW-route
artefact, chased through RYA-782/780, not through 3D-NLTE.

**Nothing here moves the gold anchor.** This adds a product; 7.466 is unchanged.

---

## Deliverable B — the 3D-NLTE availability matrix

`data/curation/threednlte_availability.csv`, verified by
`scripts/rya817_check_3dnlte_availability.py` (coverage + Crossref DOI resolution;
22/22 DOIs verified). The classification follows Asplund, Amarsi & Grevesse 2021 (A&A 653,
A141) Table 1, which is the current authority on what treatment each solar abundance rests
on.

Note the count: the canonical set is **28 species across 27 elements** — RYA-109 counts
Fe I and Fe II separately, RYA-757 adds Zn. The checker asserts this, because a matrix
listing "28 elements" would be quietly wrong.

### Solar treatment

| status | n | elements |
|---|---|---|
| **FULL_3D_NLTE** | 12 | Li, C, N, O, Na, Mg, Al, Si, K, Ca, Fe, Ba |
| 3D_NLTE_PUBLISHED_NOT_ADOPTED | 1 | **Mn** |
| MEAN3D_NLTE (3D LTE + 1D/⟨3D⟩ non-LTE) | 9 | S, Sc, Ti, Cr, Co, Cu, Zn, Sr, Eu |
| 3D_LTE_ONLY | 5 | P, V, Ni, Y, Zr |

### Off-solar 3D-NLTE grid — where 3D is actually *usable* on other stars

| status | n | elements |
|---|---|---|
| GRID_VENDORED | 3 | **C, Fe, O** |
| GRID_PUBLIC | 1 | **Li** (Wang et al. 2021 — across the full STAGGER grid) |
| GRID_PUBLIC (ion II only) | 1 | **Ca** (Lagae et al. 2025, Ca II, metal-poor FGK) |
| GRID_MEAN3D | 1 | **Al** (Nordlander & Lind 2017 — ⟨3D⟩, not 3D) |
| SOLAR_ONLY | 6 | Ba, K, Mg, N, Na, Si |
| NONE | 15 | Co, Cr, Cu, Eu, Mn, Ni, P, S, Sc, Sr, Ti, V, Y, Zn, Zr |

### Three rows worth reading

- **Mn is contested, not simply available.** Bergemann et al. 2019 published a full
  3D-NLTE solar Mn abundance (5.52); Asplund+2021 declined to adopt it, suspecting
  incomplete UV line blanketing overestimated the photoionisation rates, and kept 5.42.
  *A 3D-NLTE model existing is not the same as it being ratified.*
- **Ti is the strongest lead in the table.** Asplund+2021 state outright that the adopted
  1D non-LTE corrections are probably **not appropriate in 3D**, and that a consistent
  3D-NLTE analysis would **raise** Ti I — which is the direction that would close our own
  Ti I / Ti II discordance. No such calculation exists yet.
- **⟨3D⟩ is not 3D.** Al and K have off-solar grids computed on temporally and spatially
  averaged 3D models. That keeps the mean stratification and discards the inhomogeneity
  term — the term that carries most of the 3D effect. They are recorded under their own
  status for that reason.

### The map, in one line

Outside Fe, C and O, **3D-NLTE is a solar-abundance technology, not a stellar-survey one.**
Only Li has a full off-solar 3D-NLTE grid over the STAGGER parameter space. Any plan to
"go 3D" across the canonical set is, today, a plan to commission calculations that do not
exist.

---

## What this ticket did NOT do

- It did not move the gold anchor, re-freeze anything, or touch the live MPIA Fe leg.
- It did not build a grid (Deliverable B is a survey, as specified).
- It did not produce an IR 3D-NLTE number. There is no honest one to produce.

## Owed / follow-ups

1. **An IR 3D-NLTE Fe product needs new 3D-NLTE calculations on IR lines** — an Amarsi-side
   ask, not a run we can do. Worth raising: the near-IR Fe I lines are exactly the ones
   large surveys (APOGEE, CRIRES+, NIRPS) depend on, and the published 3D-NLTE grid does
   not cover them.
2. **Ti 3D-NLTE** is the single most valuable missing calculation for our board (see above).
3. `pipeline/amarsi3d.py` is solar-agnostic but has only been exercised on the Sun. Procyon
   sits at Teff 6554 K, **above** the 6500 K grid ceiling, so the stellar box will refuse
   it — correctly, and loudly.
